import sys, os
from churu_identify import parse_vcf

path_root = '../../'
path_workd = path_root + 'analysis/7_snp_profile_of_nonccle2ccle/'
path_target1 = path_workd + '/nonccle_records.full.txt'
path_target2 = path_workd + '/related_records.full.txt'
path_output = path_workd + '/parsed_variants.txt'
dir_churu = path_workd + '/1697001520/'

def read_target(path):
    sam2items = {}
    with open(path) as f:
        lines = f.readlines()
    for line in lines:
        items = line.rstrip('\n').split('\t')
        sam2items[items[0]] = items[1:]
    return sam2items

if __name__=='__main__':
    sam2items1 = read_target(path_target1)
    sam2items2 = read_target(path_target2)
    sam2items = {**sam2items1, **sam2items2} # python >= 3.5
    lines = []

    for i, sam in enumerate(sam2items):
        print(sam, i)
        srr = sam2items[sam][4]
        churu = f"{dir_churu}/output/{srr}/churu.out"
        vcf = f"{dir_churu}/output/{srr}/output.pileup.vcf"
        try:
            varD = parse_vcf(vcf)
        except:
            continue # zero variant

        key_var = set(varD.keys())

        for key in list(key_var):
            vaf_obs = varD[key][2]
            lines.append(f"{sam}\t{key}\t{vaf_obs}\n")

    with open(path_output, 'w') as f:
        f.writelines(lines)

    
