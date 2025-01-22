import sys, os
import re

path_root = '../../'
path_workd = path_root + '/analysis/4_check_misidentified/'
path_name1 = path_workd + '/cellname_guesse.txt'
path_name2 = path_workd + '/cellname_by_access_number.txt'
path_cello = path_workd + '/cello2ccle.table'
path_churu = path_workd + '/1697001520.parse.table'
path_out = path_workd + '/1697001520.compare_result.before_verification.txt'
#  path_out = path_workd + '/1697001520.compare_result.before_verification.test.txt'

def read_churu(path):
    sam2calls = {}
    sam2support= {}
    with open(path) as f:
        lines = f.readlines()
    for line in lines:
        items = line.rstrip('\n').split('\t')
        sam = items[10]
        n_snp = int(items[8]) if items[8] != '-' else 0 # best cell matched SNP
        n_snp_tatal = int(items[9]) if items[9] != '-' else 0 # total SNP
        if items[14] != '-':
            calls = items[14].split(';')
            achs = [x.split('|')[0] for x in calls]
        else:
            achs = []
        sam2calls[sam] = achs # ACHs
        sam2support[sam] = n_snp # n snp support of top call
    return sam2calls, sam2support

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

def read_name(path):
    sam2name = {}
    with open(path) as f:
        lines = f.readlines()
    for line in lines:
        items = line.rstrip('\n').split('\t')
        sam = items[0]
        name = items[1]
        sam2name[sam] = name
    return sam2name

if __name__=='__main__':
    sam2calls, sam2support = read_churu(path_churu)
    child2parents, name2synonymus, name2cvcls, ach2name, ach2name_lineage = read_cello(path_cello)
    sam2name1 = read_name(path_name1)
    sam2name2 = read_name(path_name2)

    # Some correction for errornouse cases
    #  child2parent['DOV13'] = 'HEK293' # some how, DOV13 matches to HEK293
    ach2name['ACH-001134'] = set(['MYLA']) # just some imputation
    #  name2synonyms['MYLA'] = [] # just some imputation

    counter = {'raw':0,
            'skip1 - noguess':0,
            'skip2 - nSNP':0,
            'skip3 - Daudi':0,
            'skip4 - score':0,
            'analysis_target':0,
            'match_same':0,
            'match_lineage':0,
            'match_synonym':0,
            'mismatch':0}
    lines_out = []
    for sam in sam2calls:
        # cell name from document
        name_gpt = sam2name1[sam].replace('none of these' ,'')
        names_acc = sam2name2[sam].split(';;;')
        names_ref = set([name_gpt] + names_acc) - set([''])

        # cell identified by churu
        call_achs = sam2calls[sam]

        # Skip
        counter['raw'] += 1
        if len(names_ref) == 0: # no guess from document
            match = "skip1"
            counter['skip1 - noguess'] += 1
        elif sam2support[sam] < 10:
            match = "skip2"
            counter['skip2 - nSNP'] += 1
        elif 'ACH-000786' in call_achs: # "Daudi" is a common false positive
            match = "skip3"
            counter['skip3 - Daudi'] += 1
        elif not call_achs: # cell must match at least 10 SNPs
            match = "skip4"
            counter['skip4 - score'] += 1
        else:

            # cello2ccle match information
            match_equal = set([name for ach in call_achs for name in ach2name[ach]])
            match_lineage = set([name for ach in call_achs for name in ach2name_lineage[ach]])

            # synonymous match information
            parents_ref = set([p for n in list(names_ref) for p in child2parents[n]])
            parents_call = set([p for n in list(match_equal) for p in child2parents[n]])
            synonyms_call = set([s for p in list(match_equal | parents_call) for s in name2synonymus[p]])

            is_equal = bool(names_ref & match_equal)
            is_lineage = bool(names_ref & match_lineage)
            has_shared_synonym = bool((names_ref | parents_ref) & synonyms_call)

            counter['analysis_target'] += 1
            if is_equal:
                match = "match_same"
                counter['match_same'] += 1
            elif is_lineage:
                match = "match_lineage"
                counter['match_lineage'] += 1
            elif has_shared_synonym: # SS
                match = "match_synonym"
                counter['match_synonym'] += 1
            else:
                match = "mismatch"
                counter['mismatch'] += 1
                #  print(sam, parents_set_cal, parents_set_ref)

        call_best = list(ach2name[call_achs[0]])[0] if call_achs else '-'
        refnames = '; '.join(sorted(list(names_ref)))
        line_out = f"{sam}\t{call_best}\t{refnames}\t{match}\n"
        lines_out.append(line_out)

    with open(path_out, 'w') as f:
        f.writelines(lines_out)

    for k in counter:
        print(f"{k}:\t{counter[k]}")
