import sys, os

path_root = '../../'
path_workd = path_root + '/analysis/4_check_misidentified/'
path_match = path_workd + '/1697001520.compare_result.final.txt'
#  path_match = path_workd + '/1697001520.compare_result.final.representative.txt'
# > representative should not be used. hierarchical information is used here.
path_proj = path_workd + '/SRA_churu_output_compare.1697001520.meta'
path_cello = path_workd + '/cello2ccle.table'
path_out = path_workd + '/inproject_swap.txt'

def read_match(path):
    sam2items = {}
    with open(path) as f:
        lines = f.readlines()
    for line in lines:
        items = line.rstrip('\n').split('\t')
        sam2items[items[0]] = items[1:]
    return sam2items


def read_proj(path):
    sam2proj, proj2sam = {}, {}
    with open(path) as f:
        lines = f.readlines()
    for line in lines:
        items = line.rstrip('\n').split('\t')
        sam, proj = items[17:19]
        sam2proj[sam] = proj
        proj2sam.setdefault(proj,[]).append(sam)
    return sam2proj, proj2sam


def read_cello(path):
    child2parents = {}
    name2cvcls = {}
    name2synonymus = {}
    ach2name = {}
    ach2name_lineage = {}
    with open(path) as f:
        lines = f.readlines()
    for line in lines[1:]: # header
        items = line.rstrip('\n').split('\t')
        child = items[0]
        parents = items[1].split('; ')
        synonymus_to = items[4].split('; ')
        cvcls = items[5].split('; ')
        achs = items[6].split('; ')
        achs_lineage = items[8].split('; ')
        child2parents[child] = parents
        name2cvcls[child] = cvcls
        name2synonymus[child] = synonymus_to
        for ach in achs:
            ach2name.setdefault(ach, set()).add(child)
        for ach in achs_lineage:
            ach2name_lineage.setdefault(ach, set()).add(child)
    return child2parents, name2synonymus, name2cvcls, ach2name, ach2name_lineage


if __name__=="__main__":
    sam2items = read_match(path_match)
    sam2proj, proj2sam = read_proj(path_proj)
    child2parents, name2synonymus, name2cvcls, ach2name, ach2name_lineage = read_cello(path_cello)

    lines_new = []
    counter = {'match': 0, 'lineage': 0, 'synonym': 0}
    for sam in sam2items:
        call, claim, _, curation = sam2items[sam]
        if not curation.startswith('contam-unknown'): # 409+169 = 578
        #  if not curation.startswith('contam-unknown-ccle'): # 298+131 = 429
        #  if not curation.startswith('contam-unknown-nonccle'): # 111+38 = 149
            continue

        parents_call = set(child2parents[call])
        synonymus_call = set([y for x in [call] + list(parents_call) for y in name2synonymus[x]])

        proj = sam2proj[sam]
        sams = proj2sam[proj]
        claims = set([_claim for _sam in sams for _claim in sam2items[_sam][1].split('; ')]) - set(['', '_unsure', '_primary'])
        # HUVEC- and HBMEC- is common chatGPT mistake for primary cell line
        claims = set(list(filter(lambda x: not (x.startswith('HUVEC-') or x.startswith('HBMEC-') ), claims)))

        parents_ref = set([y for x in list(claims) for y in child2parents[x]])

        claims_text = '; '.join(list(claims))

        match = False
        if call in claims:
            match = True
            counter['match'] += 1
        elif parents_call & parents_ref:
            match = True
            counter['lineage'] += 1
        elif synonymus_call & (parents_ref | claims):
            match = True
            counter['synonym'] += 1
        if match:
            lines_new.append(f"{sam}\t{call}\t{claim}\t{proj}\t{claims_text}\n")
    with open(path_out, 'w') as f:
        f.writelines(lines_new)
    print(counter)
