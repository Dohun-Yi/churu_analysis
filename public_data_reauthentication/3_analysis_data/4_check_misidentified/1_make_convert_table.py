import sys, os
import csv

path_root = '../../'
path_workd = path_root + '/analysis/4_check_misidentified/'
path_cello = path_workd + '/cellosaurus.table.reform.acc.cvcl'
path_ccle_model = path_workd + '/Model.tsv'
path_ccle_var = path_workd + '/OmicsSomaticMutations.csv'
path_out = path_workd + '/cello2ccle.table'

def read_cello(path):
    cvcl2name = {}
    name2cvcls = {}
    child2parent = {}
    parent2childs = {}
    name2synonyms = {}
    synonym2names = {}
    with open(path) as f:
        lines = f.readlines()
    for line in lines:
        items = line.rstrip('\n').split('\t')
        child = items[0]
        child_clean = items[1]
        parent = items[2]
        synonyms = [] if items[5] == '-' else items[5].split('; ')
        synonyms += [child_clean] # add itself
        items[11] = items[11].strip(' ')
        #  depmap_ids = [x.split(':')[1] for x in items[9].split(';') if 'DepMap:' in x]
        cvcls = [items[10]] + items[11].split('; ')
        for cvcl in cvcls:
            if cvcl == '': continue
            if cvcl in cvcl2name:
                exit(f'error {cvcl} is already in cvcl2name dictionary')
            cvcl2name[cvcl] = child
            name2cvcls.setdefault(child,[]).append(cvcl)
        child2parent[child] = parent
        name2synonyms[child] = synonyms
    for child, parent in child2parent.items():
        parent2childs.setdefault(parent, []).append(child)
    for name, synonyms in name2synonyms.items():
        for synonym in synonyms:
            synonym2names.setdefault(synonym, []).append(name)
    return name2cvcls, child2parent, parent2childs, name2synonyms, synonym2names

def read_ccle_model(path, targets=set()):
    # there are omitted records in CCLE - manual imputation
    ach2cvcl = {
            'ACH-001134': 'MYLA', # doesn't match to cello, special case
            'ACH-001172': 'CVCL_0021', # U-251 MG
            'ACH-002002': 'CVCL_C8YB', # A375-ER-2
            'ACH-002349': 'CVCL_DG57', # GBM001
            'ACH-002462': 'CVCL_C8DQ', # RPE1SS48
            'ACH-002463': 'CVCL_C8DT', # RPE1SS77
            'ACH-002464': 'CVCL_C8DS', # RPE1SS6
            'ACH-002465': 'CVCL_C8DN', # RPE1SS119
            'ACH-002466': 'CVCL_C8DM', # RPE1SS111
            'ACH-002467': 'CVCL_C8DR', # RPE1SS51
            'ACH-002471': 'CVCL_C8FR', # PSS008
            'ACH-002680': 'CVCL_LM81', # 170MGBA
            'ACH-002834': 'CVCL_C8FU', # PSS131R
            'ACH-002176': 'CVCL_1588' # NCIH748
            }
    ach2items = {
            'ACH-002176': ['NCI-H748', 'NCIH748', 'CVCL_1588']
            }
    cvcl2achs = {}
    with open(path) as f:
        lines = f.readlines()
    for line in lines[1:]:
        items = line.rstrip('\n').split('\t')
        ach  = items[0]
        name_raw = items[2]
        name_strip = items[3]
        cvcl = items[7]
        if ach not in targets: continue
        if ach not in ach2cvcl:
            ach2cvcl[ach] = cvcl
            if cvcl == '':
                exit(f"error on line\n> {line}")
        if ach not in ach2items:
            if cvcl == '':
                cvcl = ach2cvcl[ach]
            ach2items[ach] = [name_raw, name_strip, cvcl]

    for ach, cvcl in ach2cvcl.items():
        cvcl2achs.setdefault(cvcl, []).append(ach)
    return ach2items, ach2cvcl, cvcl2achs

def read_ccle_var(path):
    ach_with_var = set()
    with open(path) as f:
        reader = csv.reader(f, delimiter=',', quotechar='"')
        next(reader) # pass header
        for row in reader:
            ach = row[-2]
            ach_with_var.add(ach)
    return ach_with_var

def get_lineage(child, child2parent):
    # this dict contain information about "originate from same individual"
    parents = child2parent.get(child,'').split('; ')
    if [child] == parents:
        return parents
    else:
        lst = [pp for p in parents for pp in get_lineage(p, child2parent)]
        return parents + lst

def get_progenitors(child, child2parent):
    # child2parent contain information about "originate from same individual"
    # some cells has multiple parents
    parents = child2parent.get(child,'').split('; ')
    if [child] == parents:
        return parents
    else:
        lst = [pp for p in parents for pp in get_progenitors(p, child2parent)]
        return lst

def reformat_achs(achs, ach2items):
    text_achs = '; '.join(achs) # achs
    text_cvcls = '; '.join([ach2items[x][2] for x in achs]) # cvcl
    return f"{text_achs}\t{text_cvcls}"

def parse_progenitors(child2parent):
    child2progenitors = {}
    progenitor2childs = {}
    for child in child2parent:
        protenitors = get_progenitors(child, child2parent)
        child2progenitors[child] = protenitors
        for progenitor in protenitors:
            progenitor2childs.setdefault(progenitor, []).append(child)
    return child2progenitors, progenitor2childs


if __name__=='__main__':
    ach_with_var = read_ccle_var(path_ccle_var)
    ach2items, ach2cvcl, cvcl2achs = read_ccle_model(path_ccle_model, targets=ach_with_var)

    name2cvcls, child2parent, parent2childs, name2synonyms, synonym2names = read_cello(path_cello)
    child2parent['DOV13'] = 'HEK293' # some how, DOV13 matches to HEK293
    name2cvcls['MYLA'] = ['MYLA'] # just some imputation
    name2synonyms['MYLA'] = [] # just some imputation

    child2progenitors, progenitor2childs = parse_progenitors(child2parent)

    converter = lambda cvcl: [ach for ach in cvcl2achs.get(cvcl, [])]

    header = "cell\tparents\tprogenitors\tsynonyms\tsynonymous_to\tCVCLs\t"
    header += "CCLE\tCCLE(CVCL)\tCCLE_lineage\tCCLE_lineage(CVCL)\n"
    lines_new = [header]
    for line in open(path_cello):
        items = line.rstrip('\n').split('\t')
        child = items[0]

        parents = list(set(get_lineage(child, child2parent)))

        progenitors = child2progenitors[child]
        lineage_members = [c for p in progenitors for c in progenitor2childs[p]]

        synonyms = name2synonyms[child]
        synonymus_to = list(set([n for s in synonyms for n in synonym2names[s]]) - set([child]))

        # get cello > ccle match
        cvcls = [items[10]] + (items[11].split('; ') if items[11] != '' else [])
        achs = [ach for cvcl in cvcls for ach in converter(cvcl)]

        # get cello > ccle in lineage level
        cvcls_lineage = [cvcl for child in lineage_members for cvcl in name2cvcls[child]]
        achs_lineage = [ach for cvcl in cvcls_lineage for ach in converter(cvcl)]
        
        # Child, Parents, Progenitors, Synonyms, CCLEs (ACH, name, cvcl), CCLE_lineages (ACH, name, cvcl)
        items = [child, '; '.join(parents)]
        items += ['; '.join(progenitors), '; '.join(synonyms)]
        items += [ '; '.join(synonymus_to), '; '.join(cvcls)]
        items += [reformat_achs(achs, ach2items), reformat_achs(achs_lineage, ach2items)]
        line_new = '\t'.join(items) + '\n'
        lines_new.append(line_new)

    with open(path_out, 'w') as f:
        f.writelines(lines_new)

