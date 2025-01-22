import sys, os
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
import numpy as np
from scipy.spatial import distance


if len(sys.argv) != 2:
    exit()
batch = int(sys.argv[1])

path_root = '../../'
path_workd = path_root + 'analysis/7_snp_profile_of_nonccle2ccle/'
path_sam1 = path_workd + '/nonccle_records.full.txt'
path_sam2 = path_workd + '/related_records.full.txt'
path_var = path_workd + '/parsed_variants.txt'
path_out = path_workd + f'/result_group.{batch}.txt'

def read_var(path):
    sam2var = {}
    for line in open(path):
        items = line.rstrip('\n').split('\t')
        sam2var.setdefault(items[0],{})[items[1]] = items[2]
    return sam2var

def read_sam(path):
    sam2items = {}
    for line in open(path):
        items = line.rstrip('\n').split('\t')
        sam2items[items[0]] = items[1:]
    return sam2items

def find_related(patteern2sam, sam2items):
    # pattern: [
    #   (1) cell A > B (contam-unknown-nonccle)
    #   (2) cell A (skip, other projects)
    #   (3) cell B (match)
    # ]
    claimproj = lambda sam: f"{sam2items[sam][1]} , {sam2items[sam][5]}"
    match_call = lambda sam, call: sam2items[sam][3] == 'match' and call in sam2items[sam][0:2]
    pattern2related = {}
    for pattern in pattern2sam:
        claim, call = pattern.split(' > ')
        sams_same_claim = [s for s in sam2items if sam2items[s][1] == claim]
        # (1) group 1: cell A > B (contam)
        sams1 = pattern2sam[pattern]
        # (2) group 2: cell A (skip)
        project_set = set([sam2items[s][5] for s in sams1]) # projects of sams1
        sams2 = filter(lambda s: sam2items[s][5] not in project_set, sams_same_claim)
        sams2 = filter(lambda s: sam2items[s][3] == 'skip', sams2)
        # (3) group: cell B (match)
        sams3 = filter(lambda s: sam2items[s][3] == 'match', sam2items.keys())
        sams3 = filter(lambda s: call in sam2items[s][0:2], sams3)
        pattern2related[pattern] = [sams1, list(sams2), list(sams3)]
    return pattern2related

def make_vaf_vector(sams14, sams, sam2var):
    variants = list(set([v for s in sams14 for v in sam2var[s].keys()]))
    sam2vec = {}
    for sam in sams:
        sam2vec[sam] = np.array([float(sam2var[sam].get(v,0)) for v in variants])
    return variants, sam2vec

def get_centroid(sams, sam2vec):
    centroid = np.array([0.0] * len(sam2vec[sams[0]]))
    for sam in sams:
        vec = sam2vec[sam]
        centroid += np.array(vec)
    return centroid / len(sams)


if __name__=='__main__':
    print('reading sam...', flush=True)
    sam2items1 = read_sam(path_sam1)
    sam2items2 = read_sam(path_sam2)
    sam2items = {**sam2items1, **sam2items2}
    print('reading var...', flush=True)
    sam2var = read_var(path_var)

    print('searching for patterns...', flush=True)
    pattern2sam = {}
    for sam in sam2items:
        items = sam2items[sam]
        if items[3].startswith('contam-unknown-nonccle'):
            pattern = f"{items[1]} > {items[0]}"
            pattern2sam.setdefault(pattern,[]).append(sam)

    print('searching for related samples...', flush=True)
    pattern2related = find_related(pattern2sam, sam2items)

    print('making groups per pattern...', flush=True)
    lines_group = []
    for i, pattern in enumerate(pattern2related):
        if i != batch: continue
        samlist = pattern2related[pattern]
        size = ','.join(map(lambda x: str(len(x)),samlist))
        print(f"{i}/{len(pattern2related)} - {pattern} size={size}")
        for j, sams in enumerate(samlist):
            group = str(j+1)
            for sam in sams:
                line_base = f'{pattern}\t{sam}\t{group}\t'
                lines_group += [line_base + f"{k}\t{v}\n" for k, v in sam2var.get(sam,{}).items()]

    with open(path_out, 'w') as f:
        f.writelines(lines_group)
