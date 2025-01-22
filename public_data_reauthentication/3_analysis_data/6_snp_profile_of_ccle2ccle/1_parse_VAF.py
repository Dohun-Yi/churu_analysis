import sys, os
import csv
from churu_identify import load_CCLE_snps, parse_vcf, calc_fraction_dna, calc_fraction_rna

path_root = '../../'
path_workd = path_root + 'analysis/6_snp_profile_of_ccle2ccle/'
path_match = path_workd + '/1697001520.compare_result.final.txt' 
path_meta = path_workd + '/SRA_churu_output_compare.1697001520.meta'
path_ccle = path_workd + '/OmicsSomaticMutations.csv'
path_ccle_model = path_workd + '/Model.tsv'
path_ccle_expr = path_workd + '/OmicsExpressionProteinCodingGenesTPMLogp1.csv'
path_cello = path_workd + '/cello2ccle.table'
path_output = path_workd + '/parsed_variants.txt'
dir_churu = path_workd + '/1697001520/'

def read_match(path):
    sam2items = {}
    with open(path) as f:
        lines = f.readlines()
    for line in lines:
        items = line.rstrip('\n').split('\t')
        sam2items[items[0]] = items[1:]
    return sam2items


def read_meta(path):
    sam2srr = {}
    sam2source = {}
    with open(path) as f:
        lines = f.readlines()
    for line in lines:
        items = line.rstrip('\n').split('\t')
        srr, sam = items[0], items[17]
        source = items[21]
        sam2srr[sam] = srr
        sam2source[sam] = source
    return sam2srr, sam2source


def read_cello(path):
    name2achs = {}
    name2achs_lineage = {}
    with open(path) as f:
        lines = f.readlines()
    for line in lines[1:]: # header
        items = line.rstrip('\n').split('\t')
        child = items[0]
        achs = items[6].split('; ') # list
        achs_lineage = items[8].split('; ') # list
        name2achs[child] = achs
        name2achs_lineage[child] = achs_lineage
    return name2achs, name2achs_lineage


def parse_churu(path):
    ach2nSNP = {}
    with open(path) as f:
        lines = f.readlines()
    for line in lines[1:]:
        items = line.rstrip('\n').split('\t')
        ach, nSNP = items[0], items[3]
        ach2nSNP[ach] = int(nSNP)
    return ach2nSNP

def get_ccle_expr(path):
    exprD = {} # { cell_id : { gene : TPM } }
    try:
        f = open(path, 'r')
        reader = csv.reader(f)
        genes = next(reader)[1:] # header
        for items in reader:
            cell_id = items[0]
            expressions = items[1:]
            exprD[cell_id] = {x.split()[0]:y for x, y in zip(genes, expressions)}
        f.close()
    except IOError as e:
        msg = "I/O error while reading '%s': %s" % (path, str(e))
        raise Exception(msg)
    return exprD


if __name__=='__main__':
    sam2items = read_match(path_match)
    sam2srr, sam2source = read_meta(path_meta)
    name2achs, name2achs_lineage = read_cello(path_cello)
    refD, _ = load_CCLE_snps(path_ccle, 1e-10, 'ccle')
    exprD = get_ccle_expr(path_ccle_expr)


    sam_of_interest = [sam for sam, items in sam2items.items() if items[3] == 'contam-unknown-ccle']
    lines = []

    for sam in sam_of_interest:
        # read related files
        srr = sam2srr[sam]
        churu = f"{dir_churu}/output/{srr}/churu.out"
        vcf = f"{dir_churu}/output/{srr}/output.pileup.vcf"
        ach2nSNP = parse_churu(churu)
        varD = parse_vcf(vcf)

        # select the cell with best SNP match if multiple
        call, claim = sam2items[sam][:2]
        _calls = name2achs[call]
        _claims = name2achs_lineage[claim]
        
        _call = max(_calls, key=lambda x: ach2nSNP[x])
        _claim = max(_claims, key=lambda x: ach2nSNP[x])
        print(f"{sam}\t{call}\t{_call}\t{claim}\t{_claim}")

        var_call = refD[_call]
        var_claim = refD[_claim]

        key_var = set(varD.keys())
        key_call = set(var_call.keys())
        key_claim = set(var_claim.keys())
        key_intersection = key_call & key_claim

        key_call = key_var & key_call - key_intersection
        key_claim = key_var & key_claim - key_intersection

        source = sam2source[sam]
        if source == 'TRANSCRIPTOMIC' and _claim in exprD and _call in exprD:
            frac_list = calc_fraction_rna(varD, refD, [_call, _claim], \
                        exprD=exprD)
            mode = 'RNAmode'
        else:
            frac_list = calc_fraction_dna(varD, refD, [_call, _claim])
            mode = 'DNAmode'

        frac_call, frac_claim = frac_list
        # format: chr:pos:var:VAF_obs:VAF_ref:genotype
        for key in list(key_call):
            vaf_obs = varD[key][2]
            vaf_ref, _, genotype = var_call[key]
            lines.append(f"{sam}\tcall\t{call}\t{_call}\t{vaf_obs}\t{vaf_ref}\t{genotype}\t{frac_call}\t{source}\t{mode}\n")
        for key in list(key_claim):
            vaf_obs = varD[key][2]
            vaf_ref, _, genotype = var_claim[key]
            lines.append(f"{sam}\tclaim\t{claim}\t{_claim}\t{vaf_obs}\t{vaf_ref}\t{genotype}\t{frac_claim}\t{source}\t{mode}\n")

    with open(path_output, 'w') as f:
        f.writelines(lines)

    
