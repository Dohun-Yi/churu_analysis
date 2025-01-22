#!/usr/bin/env python

# TODO: write liscence
# TODO: check python2 and python3

"""
    TODO: Write LISCENCE Here
"""

version = "churu_who_am_i.py Version 0.0"

program_info=\
"""
identify.py: 
    Calculate posterior probability of cell match, and fraction of contamination

Version 0.0

=============================================================
"""

import os
os.environ["MKL_NUM_THREADS"] = "1" 
os.environ["NUMEXPR_NUM_THREADS"] = "1" 
os.environ["OMP_NUM_THREADS"] = "1" 
import sys, argparse, csv, math
import numpy as np
#  from scipy.stats import binom
from scipy.optimize import minimize
from scipy.special import betaln, comb



def print_header():
    """ Print the heaer line. """
    header = [
        "cell_id",      #1 - DepMap ID of cell line
        "cell_name",    #2 - Stripped name of cell line
        "patient_id",   #3 - Patient of origin
        "n_snp",        #4 - Detected SNPs from the cell
        "frac_overlap", #5 - n_snp / total detected SNPs
        "posterior",    #6 - Posterior probability of cell match
        "estimated_fraction", #7 - Estimated fraction of contaminating cell
        ]
    print('\t'.join(map(str,header)))


def parse_arguments():
    # Define parameters
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=program_info,
        #  version=version, # TODO: how to make this work on python3?
        epilog= \
"""
=============================================================
Example 1) - minimal output, only cell match
    $ identify.py [sample.vcf] [CCLE_variant_file.csv]

Example 2) - maximal output, include fractiontion
    $ identify.py [sample.vcf] [CCLE_variant_file.csv] \\
            --ccle-expression-file [CCLE_expression_file.csv] \\
            --ccle-model-file [CCLE_model_file.csv]

For support or inquiries, please contact: kutarballoon@gmail.com
""")

    ## Positional arguments
    # sample file in VCF format
    parser.add_argument('sample',
                        help='Sample variant file in VCF format')

    # CCLE reference files to compare against
    parser.add_argument('reference',metavar='reference',
                        help='ccle variant file to match against')


    ## Optional arguments
    parser.add_argument("-e", "--error-rate", type=float, metavar="E",
                        default=0.001,
                        help="sequencing error rate. DEPRECATED. default: %(default)s")

    parser.add_argument("-p", "--initial-prior", type=float, default=1e-10,
                        help="initial prior probability. default: %(default)s")

    parser.add_argument("--minimum-probability", type=float, default=1e-5,
                        help="minimum probability for each snp match, "\
                        "for log calculation. 1e-300=deactivation. default: %(default)s")

    parser.add_argument("--downsample", type=float, default=100,
                        help="maximum reads to calculate binomial. default: %(default)s")

    parser.add_argument("--minimum-popAF", type=float, default=1e-10,
                        help="minimum population allele frequency. default: %(default)s")

    parser.add_argument("--genotype-error", type=float, default=1e-5,
                        help="error rate for genotyping. default: %(default)s")

    parser.add_argument("-c", "--match-cutoff", type=float, default=0.5,
                        help="match cutoff of posterior probability. default: %(default)s")

    parser.add_argument("--popAF-mode", choices=['pop','ccle'], default='ccle',
                        help="popAF to use. pop: popAF from CCLE data, " \
                        "ccle: recalculated based on ccle variants. default: %(default)s")
    
    parser.add_argument("-n", "--nucl", choices=['rna','dna'], default='rna',
                        help=\
                        "Using RNA, cell fraction will be calculated considering gene expression (from --cell-expression-file). " \
                        "Using DNA, gene expression will be ignored. " \
                        "default: %(default)s")

    parser.add_argument("--ccle-model-file", metavar="FILE", type=str,
                        help="CCLE reference files for cell name parsing (Model.csv). " \
                            "If not provided, the output's field 2 (cell_name) and 3 (patient_id) will be empty.")

    parser.add_argument("--ccle-expression-file", metavar="FILE", type=str,
                        help="CCLE reference files for gene expression profile "\
                            "(OmicsExpressionProteinCodingGenesTPMLogp1.csv). "\
                            "if not provided, fraction will not be calculated.")

    parser.add_argument('--debug', action='store_true', help='debug mode')

    args = parser.parse_args()

    # validate values
    if args.error_rate < 0 or args.error_rate > 1:
        sys.exit("error: invalid error rate '%s'" \
                 % args.error_rate)

    if args.initial_prior < 0 or args.initial_prior > 1:
        sys.exit("error: invalid initial-prior '%s'" \
                 % args.initial_prior)

    if args.nucl == 'dna' and args.ccle_expression_file:
        sys.stderr.write("warning: nucl is set to DNA, --ccle-expression-file is ignored\n")

    if args.nucl == 'rna' and not args.ccle_expression_file:
        sys.stderr.write("error: option rna require --ccle-expression-file\n")

    if not os.path.isfile(args.sample):
        sys.exit("error: input file not found '%s'" \
                 % args.sample)

    if not os.path.isfile(args.reference):
        sys.exit("error: input file not found '%s'" \
                 % args.reference)

    if args.ccle_model_file and not os.path.isfile(args.ccle_model_file):
        sys.exit("error: input file not found '%s'" \
                 % args.ccle_model_file)

    if args.ccle_expression_file and not os.path.isfile(args.ccle_expression_file):
        sys.exit("error: input file not found '%s'" \
                 % args.ccle_expression_file)

    return args


def load_CCLE_model(path):
    """
    Loads a CCLE model file.
    
    required fields: 
        CellID (0) - For Unique Key generation
        PatientID (1) - To infer cell lines from same individuals
        CellName (2) - Name of the cell (e.g. THP-1) (Note that this field is often omitted)
        StrippedName (3) - Stripped name of the cell (e.g. THP1)

    return dictionary:
        nameD[CellID] = [PatientID, CellName, StrippedName]

    """
    nameD = {}
    try:
        f = open(path,'r')
        reader = csv.reader(f)
        #  reader.next() # header
        next(reader) # header
        for i, items in enumerate(reader):
            err = "input error in %s:%s: " % (path, i+1)

            if len(items)!=24:
                msg = "%s expecting 24 fields, got %d" % (err, len(items))
                raise Exception(msg)

            cell_id = items[0]
            patient_id = items[1]
            cell_name = items[2]
            stripped_name = items[3]

            nameD[cell_id] = [patient_id, cell_name, stripped_name]

        f.close()

    except IOError as e:
        msg = "I/O error while reading '%s': %s" % (path, str(e))
        raise Exception(msg)

    return nameD


def load_CCLE_snps(path, min_popAF=None, popAF_mode=None):
    """
    Loads a CCLE variant file.
    
    required fields: 
        Chrom (0), Pos (1), RefBase (2), AltBase (3) - For Unique Key generation
        VAF (4) - for H1 hypothesis
        GeneSymbol (13) - Gene Symbol, for expression match
        MAF (39) - -log10(MAF) value, for H2 hypothesis
        CellID (54) - Cell id, for table lookup

    return dictionary:
        refD[CellID] = {Key : [VAF, GeneSymbol]}
        mafD[key] = maf 
        # overwirte to maximum if multiple contrasting MAF exists

    """
    refD = {}
    mafD = {}
    try:
        f=open(path,'r')
        reader = csv.reader(f)
        #  reader.next() # header
        next(reader) # header
        for i, items in enumerate(reader):
            err = "input error in %s:%s: " % (path, i+1)

            if len(items)!=56:
                msg = "%s expecting 56 fields, got %d" % (err, len(items))
                raise Exception(msg)

            chrom, pos, ref_base, alt_base = items[0:4]
            gt = int(items[7][0]) + int(items[7][2])
            vaf = float(items[4])
            gene = items[13]
            maf = 1 / 10 ** float(items[39]) # maf is given in -log10 scale
            if min_popAF > maf:
                maf = min_popAF
            cell_id = items[54]
            key = '%s:%s:%s%s' % (chrom, pos, ref_base, alt_base)

            """ We only consider single nucleotide variant """
            if not (len(ref_base) == len(alt_base) == 1):
                continue

            refD.setdefault(cell_id, {})[key] = [vaf, gene, gt]
            maf_prior = mafD.get(key, -1)
            mafD[key] = max(maf_prior, maf)
        f.close()

        if popAF_mode == 'ccle':
            mafD = {} # renew
            for keys in refD.values():
                for k in keys:
                    mafD.setdefault(k, 0)
                    mafD[k] = 1.0 / len(refD)


    except IOError as e:
        msg = "I/O error while reading '%s': %s" % (path, str(e))
        raise Exception(msg)

    return refD, mafD


def parse_vcf(path):
    """
    Loads a sample VCF file.
    
    required fields: 
        Chrom (0), Pos (1), RefBase (3), AltBase (4) - For Unique Key generation
        Format (8) - for supporting reads (format: GT:GQ:SDP:DP:....)
        Sample (9) - for supporting reads (format: 0/1:66:81:81:....)

    return dictionary:
        varD[key] = [AltReads, TotalReads, VAF]
    """
    varD = {}
    try:
        f = open(path, 'r')
        for i, line in enumerate(f):
            err = "input error in %s:%s: " % (path, i+1)
            items = line.rstrip('\n').split('\t')

            if items[0].startswith('#'):
                continue

            if len(items) < 10:
                msg = "%s expecting at least 10 fields, got %d" % (err, len(items))
                raise Exception(msg)

            chrom, pos, _, ref_base, alt_base = items[0:5]
            info = {x:y for x, y in zip(items[8].split(':'), items[9].split(':'))}

            key = '%s:%s:%s%s' % (chrom, pos, ref_base, alt_base)
            alt_read = int(info['AD'])
            total_read = int(info['RD']) + alt_read
            vaf = float(info['FREQ'].rstrip('%')) / 100.0

            varD[key] = [alt_read, total_read, vaf]
        f.close()
        if not varD:
            msg = "%s file contain zero variants" % path
            raise Exception(msg)
    except IOError as e:
        msg = "I/O error while reading '%s': %s" % (path, str(e))
        raise Exception(msg)
    return varD


def get_ccle_expr(path, cell_id_list):
    """
    Loads a CCLE expression file.

    required fields:
        Every Gene Expression, from selected cell ID

    return dictionary:
        exprD[CellID] = {gene : exprssion}

    """
    cell_id_list = set(cell_id_list)
    exprD = {} # { cell_id : { gene : TPM } }
    try:
        f = open(path, 'r')
        reader = csv.reader(f)
        #  genes = reader.next()[1:] # header
        genes = next(reader)[1:] # header
        for items in reader:

            if items[0] not in cell_id_list:
                continue

            cell_id = items[0]
            expressions = items[1:]
            exprD[cell_id] = {x.split()[0]:y for x, y in zip(genes, expressions)}

        f.close()

        # Validate that all cells are there
        for cell_id in cell_id_list:
            if cell_id not in exprD:
                msg = "error while reading '%s', cell ID not found: %s" % (path, cell_id)
                raise Exception(msg)

    except IOError as e:
        msg = "I/O error while reading '%s': %s" % (path, str(e))
        raise Exception(msg)
    return exprD


def calc_posterior(lines, varD, refD, mafD, error_rate=None, initial_prior=None, min_prob=None, downsample=None, error_genotype=None):
    # TODO: update instruction
    """
    Calculate posterior


    Return:
        updated lines with number of SNP match, fraction of match, and posterior probability of match.

    Theory:
        H1: sample cell is the same with given reference cell.
        H2: sample cell is NOT the same with given reference cell.

        The chance that we observe given set of SNPs in H1 is... 
        p(D|H1) = PIE ( binomial (alt_supporting_read, total_read, prob) )
        
        Here, prob is defined as below. 
        if the obserrved SNP is present on CCLE variant of the cell, 
        (1) read is from altered base, with no error, 
        (2) read is from reference base, with sequencing error.
            prob = VAF * (1 - e) + e ( 1 - VAF ) 
        Or, if the SNP is abscent in corresponding cell, we assume this an error
            prob = e

        And, the chance what we observe given set of SNPs in H2 is...
        p(D|H2) = PIE ( p_Hetero * binom(alt_supporting_read, total_read, 0.5) +
                        p_Homo   * binom(alt_supporting_read, total_read, 1-e) +
                        p_None   * binom(alt_supporting_read, total_read,   e))

        probability of observing each allele in random data is defined as... (Hardy-Weinberg Equillibrium)
            p_Hetero = 2 * maf * ( 1 - maf )
            p_Homo   = maf ** 2
            p_None   = ( 1 - maf ) ** 2

        NOTE1: sequencing error can result in three different bases. But we don't consider that.
        NOTE2: in H1, we consider CNA by directly using VAF from CCLE rather than using diploid assumption. 
            However, H2 does not consider CNA, and simply assume diploid genome

    """

    """ define beta-binomial """
    betabinom_pmf = lambda k, n, a, b: \
        np.exp(np.log(comb(n, k)) + betaln(k + a, n - k + b) - betaln(a,b))
    a_hetero, b_hetero = 100, 100
    a_homo, b_homo = 100, 1
    a_none, b_none = 1, 100

    """ downsample reads """
    if downsample != None:
        for key in varD:
            x, n = varD[key][:2]
            if n > downsample:
                scale = 1.0 * downsample / n
                varD[key][0] = int(x * scale)
                varD[key][1] = downsample

    """ Pre-calculate probabilities that every cell share """
    e = error_rate
    e_gt = error_genotype
    p1_match_hetero_dict = {}
    p1_match_homo_dict = {}
    p1_mismatch_dict = {}
    p2_dict = {} # DEBUG
    p_h2 = 0
    n_total = len(varD)
    # remove snps that not present on CCLE variant file.
    # of note, this case contain some true dinucleotide polymolphisms, which we do not count.
    keylist = [key for key in varD.keys() if key in mafD]
    for key in keylist:
        x, n = varD[key][:2]
        maf = mafD[key]
        p_hetero = betabinom_pmf(x, n, a_hetero, b_hetero)
        p_homo   = betabinom_pmf(x, n, a_homo, b_homo)
        p_none   = betabinom_pmf(x, n, a_none, b_none)
        p1_match_hetero_dict[key] = (1-e_gt) * p_hetero + \
                e_gt/2 * p_homo + \
                e_gt/2 * p_none
        p1_match_homo_dict[key] = e_gt/2 * p_hetero + \
                (1-e_gt) * p_homo + \
                e_gt/2 * p_none
        p1_mismatch_dict[key] = e_gt/2 * p_hetero + \
                e_gt/2 * p_homo + \
                (1-e_gt) * p_none
        p2 = 2 * maf * (1 - maf) * p_hetero + \
                (maf ** 2) * p_homo + \
                ((1 - maf) ** 2) * p_none
        p_h2 += math.log10(max(p2, min_prob))
        p2_dict[key] = p2 # DEBUG

    """ Calculate posterior probability of match for each cell """
    log_lines = []
    for i, line in enumerate(lines):
        cell_id = line[0]
        refD_cell = refD[cell_id]

        n_snp_match = 0
        p_h1 = 0
        p_h2_tmp = 0 # DEBUG

        for key in keylist:
            x, n = varD[key][:2]
            gt = 0
            if key not in refD_cell:
                """ assume reference allele """
                p1 = p1_mismatch_dict[key]
            else:
                n_snp_match += 1
                vaf = refD_cell[key][0]
                gt = refD_cell[key][2]
                #  p1 = binom.pmf(x, n, vaf * (1 - e) + e * (1 - vaf) )
                if gt == 1:
                    p1 = p1_match_hetero_dict[key]
                else: # gt == 2
                    p1 = p1_match_homo_dict[key]
            p_h1 += math.log10(max(p1, min_prob))
            p_h2_tmp += math.log10(max(p2_dict[key], min_prob)) # DEBUG
            p_d = (10 ** p_h1 * initial_prior + 10 ** p_h2_tmp * (initial_prior-1)) # DEBUG
            posterior_tmp = 10 ** p_h1 * initial_prior / p_d if not p_d == 0 else 0  # DEBUG
            log_lines.append([cell_id, key, gt, x, n, str(key in refD_cell), mafD[key], p1, p2_dict[key], p_h1, p_h2_tmp, posterior_tmp]) # DEBUG

        # scale the values to avoid log(0)
        # because probability can drop under 1e-300
        # which python considers equal to zero
        p_max = max(p_h1, p_h2)
        p_h1_norm = max(p_h1 + (-1 * p_max), -300)
        p_h2_norm = max(p_h2 + (-1 * p_max), -300)

        p_d_h1 = (10 ** p_h1_norm) * initial_prior
        p_d_h2 = (10 ** p_h2_norm) * (1 - initial_prior)
        posterior = p_d_h1 / (p_d_h1 + p_d_h2) if not (p_d_h1 == p_d_h2 == 0) else 0
        lines[i][3] = n_snp_match
        lines[i][4] = 100.0 * n_snp_match / n_total
        lines[i][5] = posterior
    return lines, log_lines # DEBUG


def calc_fraction2(varD, refD, cell_id_list, mode=None, exprD=None):
    # TODO: update instruction
    best, second = cell_id_list
    target_snp = set(varD.keys()) & \
            (set(refD.get(best,{}).keys()) | \
            set(refD.get(second,{}).keys()))

    v1, v2 = [], []
    e1, e2 = [], []
    o = []
    for key in target_snp:
        vaf_observed = varD[key][2]

        vaf1, symbol1, _ = refD[best].get(key, (0.0, '-', 0))
        vaf2, symbol2, _ = refD[second].get(key, (0.0, '-', 0))

        if symbol1 == symbol2 == '-':
            continue
        elif symbol1 != '-' and symbol2 != '-' and \
                symbol1 != symbol2:
            continue
        else:
            symbol = symbol1 if symbol2 == '-' else symbol2

        expr1 = float(exprD[best].get(symbol, 0.0))
        expr2 = float(exprD[second].get(symbol, 0.0))

        v1.append(vaf1)
        v2.append(vaf2)
        e1.append(expr1)
        e2.append(expr2)
        o.append(vaf_observed)

    v1, v2, e1, e2, o = [np.array(x) for x in [v1, v2, e1, e2, o]]
    
    # Objective function to minimize
    #  def objective(f):
    #      prediction = (f * v1 * e1 + (1-f) * v2 * e2) / (f * e1 + (1-f) * e2)
    #      return np.mean((o - prediction) ** 2)
    #
    #  # Define the initial guess and bounds for f1 and f2
    #  initial_guess = [0.5]
    #  bounds = [(0, 1)]
    #
    #  # Perform constrained optimization
    #  result = minimize(objective, initial_guess, method='SLSQP', bounds=bounds)
    #  def objective(f1):
    #      f2 = 1 - f1
    #      prediction = (f1 * v1 * e1 + f2 * v2 * e2) / (f1 * e1 + f2 * e2)
    #      return np.mean((o - prediction) ** 2)

    def objective(f1):
        f2 = 1 - f1
        prediction = f1 * (e1 * o - e1 * v1) + f2 * (e2 * o - e2 * v2)
        return np.sum(prediction ** 2)

    # Define the initial guess and bounds for f1 and f2
    initial_guess = [0.5]
    bounds = [(0, 1)]

    # Perform constrained optimization
    result = minimize(objective, initial_guess, method='SLSQP', bounds=bounds)

    # Extract the optimized values of b1 and b2
    f1_pred = result.x[0]
    f2_pred = 1 - f1_pred

    return [f1_pred, f2_pred]


def calc_fraction(varD, refD, cell_id_list, nucl=None, exprD=None):
    # TODO: update instruction
    """
    Calculate fraction


    Return:
        updated lines with estimated fraction of contaminated cells

    Theory:
        if multiple cells have posterior probability of match over 0.5 (default match cutoff),
        we assume multiple cells are mixed in the sample.
        After removing cells from same patients (which have highly similar variants set),
        we estimate two best matches via deconvolution

        Deconvolution is done as followings.

        O = E*Vf
        O: vector - observed VAFs
        f: vector - cell fraction
        E: matrix - Expression profile of CCLE, row-normalized
        V: matrix - VAF from CCLE
        *: Hadamard product

        C: E*V is denoted as C in code.
        C: n*2 matrix, with cij = row-normalized expression * variant allele frequency

        Given O, E, and V, we estimate f by Sequential Least SQuare Programming (SLSQP)
        with constrain of sum(f) == 1

        in DNA mode, simply O = Vf

    """

    best, second = cell_id_list
    target_snp = set(varD.keys()) & \
            (set(refD.get(best,{}).keys()) | \
            set(refD.get(second,{}).keys()))

    rows_O = []
    rows_C = []
    for key in target_snp:
        vaf_observed = varD[key][2]

        vaf1, symbol1, _ = refD[best].get(key, (0.0, '-', 0))
        vaf2, symbol2, _ = refD[second].get(key, (0.0, '-', 0))

        if nucl == 'rna':
            if symbol1 == symbol2 == '-':
                continue
            elif symbol1 != '-' and symbol2 != '-' and \
                    symbol1 != symbol2:
                continue
            else:
                symbol = symbol1 if symbol2 == '-' else symbol2

            expr1 = float(exprD[best].get(symbol, 0.0))
            expr2 = float(exprD[second].get(symbol, 0.0))

            expr_frac1 = (expr1 / (expr1 + expr2)) if expr1 + expr2 != 0 else 0
            expr_frac2 = (expr2 / (expr1 + expr2)) if expr1 + expr2 != 0 else 0

            # Hadamard product
            c1 = expr_frac1 * vaf1
            c2 = expr_frac2 * vaf2
        elif nucl == 'dna':
            c1 = vaf1
            c2 = vaf2

        rows_O.append(vaf_observed)
        rows_C.append([c1, c2])

    O = np.array(rows_O)
    C = np.array(rows_C)
    
    # Objective function to minimize
    def objective(f):
        f1, f2 = f
        predictions = f1 * C[:, 0] + f2 * C[:, 1]
        return np.mean((O - predictions) ** 2)

    # Constraint function (b1 + b2 = 1)
    def constraint(f):
        return 1.0 - np.sum(f)

    # Define the initial guess and bounds for f1 and f2
    initial_guess = [0.5, 0.5]
    bounds = [(0, 1), (0, 1)]

    # Perform constrained optimization
    result = minimize(objective, initial_guess, method='SLSQP', bounds=bounds, constraints={'type': 'eq', 'fun': constraint})

    # Extract the optimized values of b1 and b2
    f1_pred, f2_pred = result.x

    return [f1_pred, f2_pred]


def add_name(lines, nameD):
    """ Add cell name and patient id """
    for i, line in enumerate(lines):
        cell_id = line[0]
        patient_id, cell_name, stripped = nameD.get(cell_id, ['-','-','-'])
        lines[i][1] = stripped
        lines[i][2] = patient_id
    return lines


def reformat(line):
    """ Reformat output to string """
    cell_id = line[0]
    cell_name = line[1]
    patient_id = line[2]
    n_match = '%i' % line[3]
    fraction_match = '%.2f%%' % line[4]
    posterior = ('%.2f' if line[5] > 1e-3 else '%.1e') % line[5]
    fraction = '%.1f%%' % (line[6] * 100) if line[6] != '-' else '-'
    return [cell_id, cell_name, patient_id, n_match, fraction_match, posterior, fraction]



if __name__ == "__main__":
    args = parse_arguments()

    """
        load reference files
        load (1) reference SNP, (2) CCLE model file
        But, (3) CCLE expression file is not opend until list of matched cells are acquired
    """
    try:
        refD, mafD = load_CCLE_snps(args.reference, args.minimum_popAF, args.popAF_mode)
        if args.ccle_model_file:
            nameD = load_CCLE_model(args.ccle_model_file)
    except Exception as e:
        sys.stderr.write("failed to load reference file : %s\n" % str(e))
        exit()

    """ Load sample VCF file """
    try:
        varD = parse_vcf(args.sample)
    except Exception as e:
        sys.stderr.write("failed to load sample file : %s\n" % str(e))
        exit(1)

    """ Initialize output """
    _ = '-'
    lines = [[cell_id, _, _, _, -1, -1, _] for cell_id in refD.keys()]

    """ Update posterior probability of cell match """
    lines, log_lines = calc_posterior(lines, varD, refD, mafD,
            args.error_rate,
            args.initial_prior,
            args.minimum_probability,
            args.downsample,
            args.genotype_error)
    #  debug = True
    if args.debug == True:
        for line in log_lines:
            line = '\t'.join(map(str,line))
            print(line)
        exit(1)


    """ Update cell name and patient name """
    if args.ccle_model_file:
        lines = add_name(lines, nameD)

    """ Sort by score """
    lines.sort(key=lambda x: x[4], reverse=True) # 1. sort by fraction
    lines.sort(key=lambda x: x[5], reverse=True) # 2. sort by posterior

    """ Update estimated fraction of contamination """
    c = args.match_cutoff # match cutoff for posterior
    best = lines[0][0] # cell_id of best match
    best_patient = lines[0][2]
    match_idx = list(filter(lambda i: lines[i][5] > c and \
            lines[i][2] != best_patient, range(len(lines))))

    try:
        if len(match_idx) > 0:
            second_idx = match_idx[0]
            second = lines[second_idx][0] # cell_id of second match
            if args.nucl=='rna' and args.ccle_expression_file:
                exprD = get_ccle_expr(args.ccle_expression_file, [best, second])
            else:
                exprD = None
            frac_list = calc_fraction2(varD, refD, [best, second], \
                    nucl=args.nucl, exprD=exprD)

            lines[0][-1] = frac_list[0] # manual
            lines[second_idx][-1] = frac_list[1]
    except Exception as e:
        sys.stderr.write("warning: failed to load expression, skipping fractionation: %s\n" % str(e))


    """ Print output """
    print_header()
    for line in lines:
        line = reformat(line)
        print('\t'.join(line))
