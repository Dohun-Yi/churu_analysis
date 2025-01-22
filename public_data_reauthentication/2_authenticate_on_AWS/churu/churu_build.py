#!/usr/bin/env python

# TODO: check python2 and python3
# TODO: make reference for DNA mode
version = '0.0.0'

program_info=\
"""
churu_build:
    Build a reference for cell line matching. This only needs to be done once.
    The process takes approximately 5 minutes.

Version %s

=============================================================
""" % version

import sys, os
import subprocess, time
import argparse, pyfaidx, csv, re

class transcript:
    def __init__(self, chrom, strand, tid, info):
        self.tid    = tid
        self.strand = strand
        self.chrom  = chrom
        self.info   = info
        self.exon   = []
    def add(self, typ, start, end):
        s, e = sorted([int(start), int(end)])
        if typ == 'exon':
            self.exon.append( (s, e) )
        else:
            pass
    def getlength(self):
        calc = lambda x: x[1] - x[0] + 1
        l_exon = sum([calc(x) for x in self.exon])
        return l_exon
    def gettpos(self, pos):
        pos = int(pos)
        tpos = 0
        getlen = lambda x, y: abs(x - y) + 1
        strand = {'+':1, '-':-1}[self.strand]
        for exon in sorted(self.exon)[::strand]:
            if strand * (exon[::strand][1] - pos) < 0:
                tpos += getlen(*exon)
            elif exon[0] <= pos <= exon[1]:
                tpos += abs(exon[::strand][0] - pos)
                return tpos
        return -1


class annotation:
    def __init__(self, gtf):
        suff = gtf.lower().split('.')[-1]
        if suff in set(['gtf', 'gff', 'gff3']):
            self.format = suff
        else:
            sys.stderr.write('Warning: suffix not ends with gtf or gff. assuming gtf...')
            self.format = suff
        self.tD = {}
        self._parse_gtf(gtf)
    def _parse_gtf(self, gtf):
        delim1 = 'transcript_id "' if self.format=='gtf' else 'transcript:'
        delim2 = '";' if self.format=='gtf' else ';'
        tD = {}
        with open(gtf) as f:
            lines = f.readlines()
        for i, line in enumerate(lines):
            if line.startswith('#'): continue
            items = line.rstrip('\n').split('\t')
            chrom, _, typ, start, end, _, strand, _, info = items
            if typ != 'exon': continue
            tid = items[8].split(delim1)[-1].split(delim2)[0]
            if tid.endswith('_PAR_Y'):
                continue
            tid = tid.split('.')[0] # remove version
            if tid not in tD:
                tD[tid] = transcript(chrom, strand, tid, info)
            tD[tid].add(typ, start, end)
        self.tD = tD

def get_parser(header): # this will replace identify_version
    # try to identify CCLE version, and return parser accordingly
    mapper = {name:i for i, name in enumerate(header)}
    # 22Q4: DepMap_ID, Later: ModelID
    if 'DepMap_ID' in mapper:
        mapper['ModelID'] = mapper['DepMap_ID']
    return mapper

def parse_ccle(row, mapper):
    fieldnames = ['Chrom', 'Pos', 'DNAChange', 'ModelID']
    chrom, gpos, tid_HGVS, cellid = [row[mapper[field]] for field in fieldnames]

    # DNAChange is given as HGVS in -23Q2, and TID:HGVS in 23Q4-
    if not tid_HGVS: # if none, skip this line
        tid, HGVS = ('-', 'c.e')
    elif ':' in tid_HGVS: # if format is TID:HGVS
        tid, HGVS = tid_HGVS.split(':')
    else: # if format is HGVS
        tid = row[mapper['Transcript']]
        HGVS = tid_HGVS
    return chrom, int(gpos), tid, HGVS, cellid

def parse_arguments():
    # Define parameters
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=program_info,
        #  version=version, # TODO: how to make this work on python3?
        epilog= \
"""
=============================================================

Example 1) RNA mode
    $ churu_build \\
            --transcript-fa [transcripts.fa] \\
            --gtf [annotation.gtf]  \\
            --ccle-variant-file [CCLE_variant_file.csv] \\
            --output [outputDirectory]

Example 2) DNA mode
    $ churu_build \\
            --genome-fa [genome.fa] \\
            --ccle-variant-file [CCLE_variant_file.csv] \\
            --output [outputDirectory]

Example 3) RNA and DNA
    $ churu_build \\
            --genome-fa [genome.fa] \\
            --transcript-fa [transcripts.fa] \\
            --gtf [annotation.gtf]  \\
            --ccle-variant-file [CCLE_variant_file.csv] \\
            --output [outputDirectory]

For support or inquiries, please contact: kutarballoon@gmail.com
""")
    
    parser.add_argument("--transcript-fa", metavar="FILE",
                        help=\
                        "FASTA file for the reference transcriptome. "\
                        "must match the GTF file")

    parser.add_argument("--genome-fa", metavar="FILE",
                        help=\
                        "FASTA file for the reference sequence. "\
                        "on RNA mode, FASTA must match the GTF file")

    parser.add_argument("--gtf", metavar="FILE",
                        help=\
                        "GTF file for the reference transcriptome. "\
                        "Must match the FASTA file")

    parser.add_argument("--ccle-variant-file", metavar="FILE", required=True,
                        help=\
                        "CCLE reference files for mutations. " \
                        "(OmicsSomaticMutations.csv)")

    parser.add_argument("-k", "--kmer", metavar="INT", type=int, default=31,
                        help="size of k-mer; must match the -k value in churu (default: %(default)s)")

    parser.add_argument("--output", metavar="DIR", default="./ChuruReference/",
                        help=\
                        "output directory (default: %(default)s)")

    parser.add_argument("-t", "--threads", metavar="INT", type=int, default=8,
                        help="number of threads to use for KMC (default: %(default)s)")

    parser.add_argument('-v', '--version', action='version',
                        version='%(prog)s ' + version)

    args = parser.parse_args()

    if not args.transcript_fa and not args.genome_fa:
        sys.exit("Error: none of --transcript-fa and --genome_fa is set")

    if args.transcript_fa and not args.gtf:
        sys.exit("Error: --transcript-fasta is set, but --gtf is not set")

    if args.genome_fa and not os.path.isfile(args.genome_fa):
        sys.exit("Error: input file not found '%s'" \
                % args.genome_fa)

    if args.transcript_fa and not os.path.isfile(args.transcript_fa):
        sys.exit("Error: input file not found '%s'" \
                % args.transcript_fa)

    if args.gtf and not os.path.isfile(args.gtf):
        sys.exit("Error: input file not found '%s'" \
                % args.gtf)

    if not os.path.isfile(args.ccle_variant_file):
        sys.exit("Error: input file not found '%s'" \
                % args.ccle_variant_file)

    return args

def revcomp(kmer):
    comp={'A':'T','T':'A','G':'C','C':'G', 'N':'N'}
    kmer_new = ''
    for nt in kmer[::-1]:
        kmer_new += comp[nt]
    return kmer_new

def get_snp_kmer(k, seq, pos, HGVS, strand, mode='rna'):
    """
    Generate mutated k-mer and paired reference k-mer.
    """
    # type of variation
    #   substitution c.nnX>X, c.nnX>XX or c.nnXX>X
    #   deletion c.nn_nndelXXXXXXX
    #   insertion c.nn_nninsXXXXXXX
    mode = {'rna':'rna', 'dna':'dna'}[mode]
    HGVS = HGVS.lstrip('c.')
    posL = re.findall(r"\d+", HGVS)
    var = re.sub("\d+|_","",HGVS)

    seq = str(seq[max(pos - 50, 0):(pos + 50)]).upper()
    pos = min(50, pos)

    pos_s, pos_e = 0, 0
    kmer_snp, kmer_ref = '.', '.'

    if '>' in var:
        bef, aft = var.split('>')
        varlen = 1 if len(posL) == 1 else int(posL[1]) - int(posL[0]) + 1

        if mode == 'dna':
            # TODO: Better to directly use ref/alt base in DNA mode
            pos_s, pos_e = pos - 1, pos + varlen - 2
            varseq_plus = seq[pos_s:(pos_e+1)]
            bef_rev, aft_rev = revcomp(bef), revcomp(aft)
            if varseq_plus == bef:
                strand='+'
            elif varseq_plus == bef_rev:
                strand='-'
                bef, aft = bef_rev, aft_rev
            else:
                sys.exit('Error: variant does not match %s' % HGVS)
        if mode == 'rna':
            pos_s = pos if strand == '+' else pos - (varlen-1)
            pos_e = pos + (varlen-1) if strand == '+' else pos
            if bef != str(seq[pos_s:(pos_e+1)]):
                msg = 'Error: sequence mismatches: %s\n' % HGVS
                msg += 'please check if fasta and variant are in same version'
                sys.exit(msg)

        seq_new = seq[:pos_s] + aft + seq[pos_e+1:]
        # total length should be equal to k
        mid = (pos_s + pos_e) // 2
        seq_s = max(mid - k // 2, 0)
        seq_e = seq_s + k

        kmer_snp = seq_new[seq_s:seq_e]
        kmer_ref = seq[seq_s:seq_e]
    elif var.startswith('del'):
        # deletions are not covered
        pass
    elif var.startswith('ins'):
        # insertions are not covered
        pass
    else:
        pass
        #  sys.exit('Error: unknown HGVS: %s' % HGVS)
    return pos_s, pos_e, kmer_snp, kmer_ref


def find_executable(executable):
    """
    Find executable file from directory of this script, or from PATH
    """
    # Check the directory of this Python script
    script_dir = os.path.dirname(os.path.abspath(__file__))
    local_path = os.path.join(script_dir, executable)
    if os.path.isfile(local_path) and os.access(local_path, os.X_OK):
        return local_path

    # If not found, check PATH
    paths = os.environ.get("PATH", "").split(os.pathsep)
    for directory in paths:
        executable_path = os.path.join(directory, executable)
        if os.path.isfile(executable_path) and os.access(executable_path, os.X_OK):
            return executable_path

    sys.exit("Error: required executable '%s' not found in PATH or install directory"\
            % executable)


def run_kmc(kmc, input_file, output_prefix, temp_dir, k_value=31, min_count=1, threads=8):
    """Wrapper function for the `kmc` command.

    Args:
        input_file (str): Path to the input file.
        output_prefix (str): Prefix for the output file.
        temp_dir (str): Directory for temporary files.
        k_value (int): Size of k-mer (default is 31).
        min_count (int): Minimum count for k-mers (default is 1).
        threads (int): Number of threads to use (default is 30).
    """
    command = [
        kmc,
        "-k{}".format(k_value),
        "-fa",
        "-ci{}".format(min_count),
        "-t{}".format(threads),
        input_file,
        output_prefix,
        temp_dir
    ]

    if not os.path.isdir(temp_dir):
        os.mkdir(temp_dir)
    try:
        result = subprocess.check_output(command, stderr=subprocess.STDOUT)
        print("kmc completed successfully:")
        print(result.decode("utf-8"))
    except subprocess.CalledProcessError as e:
        print("Error occurred during kmc execution:")
        print(e.output.decode("utf-8"))
    if os.path.isdir(temp_dir):
        os.rmdir(temp_dir)


if __name__=='__main__':
    tic = time.time()
    args = parse_arguments()

    path_csv = args.ccle_variant_file
    path_gtf = args.gtf
    path_tfa = args.transcript_fa
    path_gfa = args.genome_fa
    path_outdir = args.output
    cpu = args.threads
    k = args.kmer
    span = 0 # 0nt for span
    rna = True if args.transcript_fa else False
    dna = True if args.genome_fa else False

    kmc = find_executable('kmc')

    if not os.path.isdir(path_outdir):
        os.mkdir(path_outdir)

    """ Step 1. Read input fils """
    if rna:
        print('reading annotation...') 
        anno = annotation(path_gtf)
        tD = anno.tD

        print('reading transcript fasta...')
        tfa = pyfaidx.Fasta(path_tfa,
                key_function = lambda x: x.split('|')[0])

    if dna:
        print('reading genome fasta...')
        gfa = pyfaidx.Fasta(path_gfa,
                key_function = lambda x: x.split('|')[0])

    print('reading csv...')
    items = []
    skipped = set()
    with open(path_csv) as f:
        reader = csv.reader(f, delimiter=',',quotechar='"')
        mapper = get_parser(next(reader))
        for row in reader:
            chrom, gpos, tid, HGVS, cellid = parse_ccle(row, mapper)
            strand = '.' # missing on specific CCLE version
            if HGVS.startswith('c.e') or \
                    HGVS.startswith('n.') or \
                    '-' in HGVS or \
                    '+' in HGVS or \
                    '*' in HGVS:
                # intronic, non-coding, splice site all skipped
                # HGVS not present in 22Q2 is also skipped
                continue
            if rna:
                tid_clean = tid.split('.')[0]
                if tid_clean not in tD:
                    if tid_clean not in skipped:
                        sys.stderr.write('Warning: Transcript ID %s is not found on GTF, please check if annotation and CCLE version matches. skipping...\n' % tid_clean)
                        skipped.add(tid_clean)
                    continue
                strand = tD[tid_clean].strand
                tpos = tD[tid_clean].gettpos(gpos) # convert genomic position to transcript position
                s, e = max(tpos-int(span/2), 0), min(tpos+int(span/2)+1, tD[tid_clean].getlength())
            else:
                s, e = 0, 0
            items.append([tid, s, e, HGVS, chrom, gpos, strand, cellid])

    """ Step 2. Generate mutated k-mer and paired matched k-mer """
    print('processing kmers...')
    lines_tfa, lines_gfa = [], []
    lines_ttwin, lines_gtwin = [], []
    uniqpos = set()
    for item in items:
        tid, s, e, HGVS, chrom, gpos, strand, cellid = item
        uniqpos.add('%s:%i' % (chrom, gpos))
        if rna:
            pos_s, pos_e, kmer_snp, kmer_ref = get_snp_kmer(k, tfa[tid], s, HGVS, strand, mode='rna')
            if kmer_snp == '.': continue
            lines_ttwin.append('%s\t%i\t%s\t%s\t%s\t%s\n' % (chrom, gpos, strand, HGVS, kmer_snp, kmer_ref))
            lines_tfa.append('>%s:%i|%s:%i-%i|%s|%s|SNP\n%s\n' % (chrom, gpos, tid, pos_s, pos_e, HGVS, cellid, kmer_snp))
            lines_tfa.append('>%s:%i|%s:%i-%i|%s|%s|REF\n%s\n' % (chrom, gpos, tid, pos_s, pos_e, HGVS, cellid, kmer_ref))
        if dna:
            # TODO: this function will not work
            pos_s, pos_e, kmer_snp, kmer_ref = get_snp_kmer(k, gfa[chrom], gpos, HGVS, strand, mode='dna')
            if kmer_snp == '.': continue
            lines_gtwin.append('%s\t%i\t%s\t%s\t%s\t%s\n' % (chrom, gpos, strand, HGVS, kmer_snp, kmer_ref))
            lines_gfa.append('>%s:%i|%s:%i-%i|%s|%s|SNP\n%s\n' % (chrom, gpos, tid, pos_s, pos_e, HGVS, cellid, kmer_snp))
            lines_gfa.append('>%s:%i|%s:%i-%i|%s|%s|REF\n%s\n' % (chrom, gpos, tid, pos_s, pos_e, HGVS, cellid, kmer_ref))


    """ Step 3. Write outputs """
    print('writing output...')
    if rna:
        path_out_tfa = path_outdir + '/rna.kmerdb.fa'
        path_out_tbed = path_outdir + '/rna.kmerdb.bed'
        path_out_tkmc = path_outdir + '/rna.kmerdb.fa.kmc'
        path_ttemp = path_outdir + '/rna.kmerdb.fa.kmc.tmp'
        with open(path_out_tbed, 'w') as f:
            f.writelines(lines_ttwin)
        with open(path_out_tfa, 'w') as f:
            f.writelines(lines_tfa)
    if dna:
        path_out_gfa = path_outdir + '/dna.kmerdb.fa'
        path_out_gbed = path_outdir + '/dna.kmerdb.bed'
        path_out_gkmc = path_outdir + '/dna.kmerdb.fa.kmc'
        path_gtemp = path_outdir + '/dna.kmerdb.fa.kmc.tmp'
        with open(path_out_gbed, 'w') as f:
            f.writelines(lines_gtwin)
        with open(path_out_gfa, 'w') as f:
            f.writelines(lines_gfa)

    path_out_pos = path_outdir + '/common.kmerdb.pos.bed'
    with open(path_out_pos, 'w') as f:
        positions = list(map(lambda x: x.split(':'), list(uniqpos)))
        positions.sort(key=lambda x: int(x[1]))
        positions.sort(key=lambda x: x[0])
        positions = list(map(lambda x: "%s\t%s\n" % (x[0],x[1]), positions))
        f.writelines(positions)
    
    print('building KMC index..')
    if rna:
        run_kmc(kmc, path_out_tfa, path_out_tkmc, path_ttemp, k_value=k, min_count=1, threads=cpu)
    if dna:
        run_kmc(kmc, path_out_gfa, path_out_gkmc, path_gtemp, k_value=k, min_count=1, threads=cpu)
    # remove temp file
    if rna and os.path.isfile(path_out_tfa):
        os.remove(path_out_tfa)
    if dna and os.path.isfile(path_out_gfa):
        os.remove(path_out_gfa)

    print('output written') 
    print('')
    print('# run time: %i sec' % (time.time() - tic))
    print('# finished')
