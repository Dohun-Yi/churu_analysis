import sys, os
from churu_identify import load_CCLE_snps

path_root = '../../'
path_workd = path_root + '/analysis/6_snp_profile_of_ccle2ccle/'
path_ccle = path_workd + '/OmicsSomaticMutations.csv'
path_out = path_workd + '/reference_snp_numbers.txt'

if __name__=='__main__':
    refD, _ = load_CCLE_snps(path_ccle, 1e-10, 'ccle')

    lines_out = []
    for ach in refD:
        lines_out.append(f'{ach}\t{len(refD[ach])}\n')

    with open(path_out, 'w') as f:
        f.writelines(lines_out)
