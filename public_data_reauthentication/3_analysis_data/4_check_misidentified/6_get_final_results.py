import sys,os

path_root = '../../'
path_workd = path_root + '/analysis/4_check_misidentified/'
path_match = path_workd + '1697001520.compare_result.before_verification.txt'
path_manual_1 = path_workd + 'manual_inspection_result.2_last_verification.tsv'
path_manual_2 = path_workd + 'manual_check_known_contam_from_match_lineage.txt'
path_out = path_workd + '1697001520.compare_result.final.txt'

def read_manual_1(path):
    sam2items = {}
    with open(path) as f:
        lines = f.readlines()
    for line in lines[1:]: # has header
        items = line.rstrip('\n').split('\t')
        sam = items[0]
        call, claim = items[3:5]
        claim_corrected = items[8]
        curation = items[6]
        if claim_corrected:
            claim = claim_corrected
        sam2items[sam] = [call, claim, curation]
    return sam2items

def read_manual_2(path):
    sam2items = {}
    converter = {
            'x': 'match_lineage',
            '.': 'contam-known-ccle',
            'p': 'contam-unknown-nonccle'
            }
    with open(path) as f:
        lines = f.readlines()
    for line in lines: # has header
        items = line.rstrip('\n').split('\t')
        sam = items[0]
        call, claim = items[1:3]
        curation = converter[items[3]]
        sam2items[sam] = [call, claim, curation]
    return sam2items

if __name__=='__main__':
    converter = {
            'match': 'match',
            'match_same': 'match',
            'match_lineage': 'match',
            'match_synonym': 'match',
            'match-like': 'match',
            'contam-known-ccle': 'contam-known-ccle',
            'contam-known-nonccle': 'contam-known-nonccle',
            'contam-unknown-ccle': 'contam-unknown-ccle',
            'contam-unknown-nonccle': 'contam-unknown-nonccle ',
            'skip1': 'skip',
            'skip2': 'skip',
            'skip3': 'skip',
            'skip4': 'skip',
            'unsure': 'skip' # this should be counted at the paper
            }
    sam2items_1 = read_manual_1(path_manual_1)
    sam2items_2 = read_manual_2(path_manual_2)

    lines_out = []
    compare = lambda x, y: set(x.split('; ')) == set(y.split('; '))
    for line in open(path_match):
        sam, call, claim, curation = line.rstrip('\n').split('\t')
        if sam in sam2items_1:
            _call, _claim, _curation = sam2items_1.get(sam)
            claim = _claim # use corrected version
            curation2 = converter[_curation]
        elif sam in sam2items_2:
            _call, _claim, _curation = sam2items_2.get(sam)
            curation2 = converter[_curation]
        else:
            curation2 = converter[curation]
        line_out = f"{sam}\t{call}\t{claim}\t{curation}\t{curation2}\n"
        lines_out.append(line_out)

    with open(path_out, 'w') as f:
        f.writelines(lines_out)
