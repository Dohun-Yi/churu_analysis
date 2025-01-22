import sys, os

path_root = '../../'
path_ref = path_root + '/resources/ROR/v1.50-2024-07-29-ror-data.filt.tsv'

path_workd = path_root + '/analysis/3_country_institute/'
path_owner = path_workd + '/all_in_one.txt'
path_in = path_workd + '/institute_selected.txt'
path_out = path_workd + '/institute_selected.parsed.txt'

def read_ror(path):
    id2info = {}
    with open(path) as f:
        lines = f.readlines()

    for line in lines[1:]: # header
        items = line.rstrip('\n').split('\t')
        ror_id, name, country = items[0], items[1], items[19]
        relation = items[-1]
        id2info[ror_id] = [name, country, relation]
    return id2info

def read_owner(path):
    sam2owner = {}
    with open(path) as f:
        lines = f.readlines()
    for line in lines:
        items = line.rstrip('\n').split('\t')
        sam = items[0]
        owner = items[6]
        sam2owner[sam] = owner
    return sam2owner

if __name__=='__main__':
    id2info = read_ror(path_ref)
    sam2owner = read_owner(path_owner)

    lines_out = []
    for line in open(path_in):
        sam, ror_id = line.rstrip('\n').split('\t')
        owner = sam2owner[sam]
        info = id2info.get(ror_id, ['-'] * 3)
        info_line = '\t'.join(info)
        line_out = f"{sam}\t{owner}\t{ror_id}\t{info_line}\n"
        lines_out.append(line_out)

    with open(path_out, 'w') as f:
        f.writelines(lines_out)
