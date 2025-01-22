import sys, os

path_root = '../../'
path_workd = path_root + '/analysis/7_snp_profile_of_nonccle2ccle/'

path_institute = path_workd + '/institute_selected.parsed.corrected.txt'
path_meta = path_workd + '/SRA_churu_output_compare.1697001520.meta'
path_records1 = path_workd + '/nonccle_records'
path_records2 = path_workd + '/related_records'
path_out1 = path_workd + '/nonccle_records.full.txt'
path_out2 = path_workd + '/related_records.full.txt'

def read_meta(path):
    sam2srr = {}
    sam2proj = {}
    with open(path) as f:
        lines = f.readlines()
    for line in lines:
        items = line.rstrip('\n').split('\t')
        srr, sam, proj = items[0], items[17], items[18]
        sam2srr[sam] = srr
        sam2proj[sam] = proj
    return sam2srr, sam2proj

def read_institute(path):
    sam2institute = {}
    with open(path) as f:
        lines = f.readlines()
    for line in lines:
        items = line.rstrip('\n').split('\t')
        sam, ror, institute = items[0], items[2], items[3]
        sam2institute[sam] = [ror, institute]
    return sam2institute

def read_and_write(path_in, path_out, params):
    sam2srr, sam2proj, sam2institute = params
    lines_new = []
    for line in open(path_in):
        items = line.rstrip('\n').split('\t')
        sam = items[0]
        items += [sam2srr[sam], sam2proj[sam]] + sam2institute.get(sam, ['DISCARD']*2)
        lines_new.append('\t'.join(items) + '\n')

    with open(path_out, 'w') as f:
        f.writelines(lines_new)

if __name__=='__main__':
    sam2srr, sam2proj = read_meta(path_meta)
    sam2institute = read_institute(path_institute)
    params = (sam2srr, sam2proj, sam2institute)

    read_and_write(path_records1, path_out1, params)
    read_and_write(path_records2, path_out2, params)
