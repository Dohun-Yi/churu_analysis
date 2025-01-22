import sys, os

path_root = '../../'
path_workd = path_root + '/analysis/4_check_misidentified/'
path_ror = path_root + '/resources/ROR/v1.50-2024-07-29-ror-data.filt.tsv'
path_match = path_workd + '/1697001520.compare_result.final.txt'
path_institute = path_workd + '/institute_selected.parsed.corrected.txt'
path_out = path_workd + '/misidentified_institute.txt'

def read_match(path):
    # sam: match
    skips = ['skip', 'contam-known-nonccle', 'contam-unknown-nonccle']
    with open(path) as f:
        lines = f.readlines()
        lines = map(lambda x: x.rstrip('\n').split('\t'), lines)
        sam2match = {l[0]:l[4] for l in lines if l[4] not in skips}
    return sam2match

def read_institute(path, conv=lambda x:x):
    # sam: [owner, ror_id, name, country, relation]
    with open(path) as f:
        lines = f.readlines()
        lines = list(map(lambda x: x.rstrip('\n').split('\t'), lines))
        sam2institute = {l[0]:l[1:4]+[conv(l[4]), l[5]] for l in lines}
    return sam2institute

def read_ror(path, cov=lambda x:x):
    # rorid: [name, country, relation]
    with open(path_ror) as f:
        lines = f.readlines()
        lines = list(map(lambda x: x.rstrip('\n').split('\t'), lines[1:]))
        ror2info = {l[0]:[l[1], conv(l[19]), l[28]] for l in lines}
    return ror2info

def count_contam(sam2match, sam2institute, ror2info):
    ror_counter = {}
    for sam in sam2institute:
        owner, ror_id, name, country, relation = sam2institute[sam]
        if ror_id not in ror2info:
            ror2info[ror_id] = [name, country, relation]

        if sam in sam2match:
            ror_counter.setdefault(ror_id, [0,0,0,0]) # match, mismatch, total, fraction
            
            idx = 0 if sam2match[sam] == 'match' else 1
            ror_counter[ror_id][idx] += 1 # match, mismatch
            ror_counter[ror_id][2] += 1 # total
            ror_counter[ror_id][3] = 1.0 * ror_counter[ror_id][1] / ror_counter[ror_id][2] # fraction
    return ror_counter

def get_parentD(ror2info, sam2institute):
    parentD = {}
    # all ROR informations
    for child_id in ror2info:
        items = ror2info[child_id]
        relation = items[-1]
        if "Parent" not in relation:
            parentD[child_id] = [items[1]] # Country
        else:
            parent = [x for x in relation.split('; ') if x.startswith('Parent')][0]
            parentD[child_id] = parent.lstrip('Parent: ').split(', ')
    # for CUSTOM records (manual imputation)
    for sam in sam2institute:
        _, child_id, _, country, relation = sam2institute[sam]
        if "Parent" not in relation:
            parentD[child_id] = [country] # Country
        else:
            parent = [x for x in relation.split('; ') if x.startswith('Parent')][0]
            parentD[child_id] = parent.lstrip('Parent: ').split(', ')
    return parentD

def crawler_list(child, parentD):
    if child in parentD and [child] != parentD[child]:
        parents = parentD[child] # is a list
        res = []
        for parent in parents:
            res += crawler_list(parent, parentD)
        return parents + res
    else:
        return []

def get_agg(ror_counter, ror_of_interest, parentD):
    ror2agg_count = {} # This one include country names now
    for ror_id in ror_of_interest:
        ror_list = list(set([ror_id] + crawler_list(ror_id, parentD)))
        for _ror_id in ror_list:
            ror2agg_count.setdefault(_ror_id, [0,0,0,0])
            ror2agg_count[_ror_id][0] += ror_counter[ror_id][0]
            ror2agg_count[_ror_id][1] += ror_counter[ror_id][1]
            ror2agg_count[_ror_id][2] += ror_counter[ror_id][2]
            ror2agg_count[_ror_id][3] = ror2agg_count[_ror_id][1] / ror2agg_count[_ror_id][2]
    return ror2agg_count

if __name__=='__main__':
    country_converter = {"The Netherlands": "Netherlands",
            "Hong Kong": "China",
            "Macao": "China",
            "Taiwan": "China"}
    conv = lambda x: country_converter.get(x,x)

    sam2match = read_match(path_match)
    sam2institute = read_institute(path_institute, conv)
    ror2info = read_ror(path_ror, conv)

    parentD = get_parentD(ror2info, sam2institute)
    ror_of_interest = list(set([sam2institute.get(sam,['.','DISCARD'])[1] for sam in sam2match]))

    ror_counter = count_contam(sam2match, sam2institute, ror2info)
    agg_counter = get_agg(ror_counter, ror_of_interest, parentD)

    fmt = lambda x: "%i\t%i\t%i\t%.4f" % tuple(x)
    lines_out = []
    for ror in agg_counter:
        institute, country = ror2info.get(ror, ['-', '-'])[:2]
        c_raw = fmt(ror_counter.get(ror,[0,0,0,0]))
        c_agg = fmt(agg_counter[ror])
        parents = '; '.join(parentD.get(ror,[ror]))
        lines_out.append(f"{ror}\t{institute}\t{country}\t{c_raw}\t{c_agg}\t{parents}\n")

    with open(path_out, 'w') as f:
        f.writelines(lines_out)
