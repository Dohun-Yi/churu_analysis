import sys, os

path_root = '../../'
path_workd = path_root + '/analysis/3_country_institute/'
path_error1 = path_workd + '/institute_error1_unidentified.txt'
path_error2 = path_workd + '/institute_error2_multiple_candidates.txt'
path_raw = path_workd + '/all_in_one.txt'
path_gpt = path_workd + '/institute_selected.parsed.txt'
path_out = path_workd + '/institute_error_parsed_for_manual_correction.txt'

def get_error(path):
    owners = set()
    for line in open(path):
        owners.add(line.rstrip('\n').lstrip(' ').split(' ',1)[-1])
    return owners

def get_GEO(path):
    sam2geo = {}
    for line in open(path):
        items = line.rstrip('\n').split('\t')
        sam, info = items[0], items[7:]
        geo = '-' if len(info) == 1 else info[1]
        sam2geo[sam] = geo
    return sam2geo

def get_institute(path):
    sam2institute = {}
    sam2id = {}
    owner2sam = {}
    for line in open(path):
        items = line.rstrip('\n').split('\t')
        sam, owner, ror_id, institute = items[0], items[1], items[2], items[3]
        sam2institute[sam] = institute
        sam2id[sam] = ror_id
        owner2sam.setdefault(owner, []).append(sam)
    return sam2institute, sam2id, owner2sam

if __name__=='__main__':
    owners_error1 = get_error(path_error1)
    owners_error2 = get_error(path_error2)
    sam2geo = get_GEO(path_raw)
    sam2institute, sam2id, owner2sam = get_institute(path_gpt)

    owner_union = owners_error1 | owners_error2
    to_uniq_string = lambda x: ';'.join(sorted(list(set(x))))
    lines_out = []
    for owner in list(owner_union):
        sams = owner2sam[owner]
        is_unidentified = 'O' if owner in owners_error1 else 'X'
        is_conflicting = 'O' if owner in owners_error2 else 'X'
        geo_uniq = to_uniq_string([sam2geo[s] for s in sams])
        institute_uniq = to_uniq_string([sam2institute[s] for s in sams])
        id_uniq = to_uniq_string([sam2id[s] for s in sams])
        sams_str = ';'.join(sams)
        line = f"{owner}\t{len(sams)}\t{is_unidentified}\t{is_conflicting}\t{geo_uniq}\t{institute_uniq}\t{id_uniq}\t{sams_str}\n"
        lines_out.append(line)

    with open(path_out, 'w') as f:
        f.writelines(lines_out)
