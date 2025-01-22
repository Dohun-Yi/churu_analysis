import sys, os
import csv

path_root = '../../'
path_ror = path_root + '/resources/ROR/v1.50-2024-07-29-ror-data.filt.tsv'
path_workd = path_root + '/analysis/3_country_institute/'
path_cor = path_workd + '/institute_error_parsed_for_manual_correction.csv' # Manually corrected file
path_in = path_workd + '/institute_selected.parsed.txt'
path_out = path_workd + '/institute_selected.parsed.corrected.txt'
"""
 - "path_cor" contains the manually corrected institute names
 - further, correction for some of international institutes needed
 - lastely, manual imputation for few omitted insitutes needed
"""

def read_cor(path):
    owner2cor = {}
    with open(path) as f:
        reader = csv.reader(f, delimiter=',', quotechar='"')
        for row in reader:
            if 'OwnerName' in row[0]:
                continue # header
            owner = row[0]
            institute = row[9]
            country = row[10]
            ror_id = row[11]
            relation = ''
            if institute in ['-', '']:
                institute, country, ror_id, relation = ['DISCARD'] * 4
            owner2cor[owner] = [institute, country, ror_id, relation]
    return owner2cor

def read_ror(path):
    id2items = {}
    with open(path) as f:
        lines = f.readlines()
    for line in lines[1:]:
        items = line.rstrip('\n').split('\t')
        id2items[items[0]] = items[1:]
    return id2items

def imputation(owner2cor):
    imputations = {
            "AbLife. Inc.": ["China", "CUSTOM_1", ""],
            "Variant Bio": ["United States", "CUSTOM_2", ""],
            "Cancer Science Institute": ["Singapore", "CUSTOM_3", "Parent: https://ror.org/01tgyzw49"],
            "ENCODE": ["United States", "CUSTOM_4", "Parent: https://ror.org/00baak391"],
            "Harvard Medical School": ["United States", "CUSTOM_5", "Parent: https://ror.org/03vek6s52"],
            "Laboratory of Molecular Medicine and Genomics": ["Italy", "CUSTOM_6", "Parent: https://ror.org/02kqnpp86, https://ror.org/02kg21794, https://ror.org/0192m2k53"],
            "Sidra Medicine": ["Catarrh", "CUSTOM_7", "Parent: https://ror.org/05v5hg569"],
            }
    counter = 0
    for owner in owner2cor:
        items = owner2cor[owner]
        name = items[0]
        if name in imputations:
            owner2cor[owner] = [name] + imputations[name]
            print(f'record manually corrected: {name} (originally "{owner}")')
            counter += 1
    print(f'owner manual imputation: {counter}')
    return owner2cor

def get_main_headquarters(path):
    headquarters = {
            "Ludwig Cancer Research": "https://ror.org/05qdwtz81",
            "Pfizer": "https://ror.org/01xdqrp08",
            "European Molecular Biology Laboratory": "https://ror.org/03mstc592",
            "Bayer": "https://ror.org/04hmn8g73",
            "Roche": "https://ror.org/00by1q217",
            "GlaxoSmithKline": "https://ror.org/01xsqw823",
            "Novartis": "https://ror.org/02f9zrr09"
            }
    institute2head = {}
    with open(path) as f:
        lines = f.readlines()
    for line in lines[1:]:
        items = line.rstrip('\n').split('\t')
        ror_id, name = items[:2]
        for h in headquarters:
            # All pattern matches are visually checked
            if not (name.startswith(h + ' ') or name == h):
                continue
            print(f'"{name}" record is converted to {h} main headquarter ({headquarters[h]})')
            institute2head[ror_id] = headquarters[h]
    return institute2head

if __name__=="__main__":
    institute2head = get_main_headquarters(path_ror)
    id2items = read_ror(path_ror)
    owner2cor = read_cor(path_cor)
    owner2cor = imputation(owner2cor)

    lines_new = []
    counter = {'owner':0, 'international':0}
    for line in open(path_in):
        items = line.rstrip('\n').split('\t')
        sam, owner = items[:2] # not change
        ror_id, institute, country, relation = items[2:] # may change
        # Correction and imputation
        if owner in owner2cor:
            counter['owner'] += 1
            institute, country, ror_id, relation = owner2cor[owner]
            #  print(owner2cor[owner])
            if not ror_id.startswith('CUSTOM') and ror_id != "DISCARD":
                relation = id2items[ror_id][-1]
        # if international institute, chagne to headquarter
        if ror_id in institute2head:
            counter['international'] += 1
            ror_id = institute2head[ror_id]
            items = id2items[ror_id]
            #  print(items)
            institute, country, relation = items[0], items[18], items[27]
        if ror_id == '-':
            # if ROR-ID is still "-", then discard
            institute, country, ror_id, relation = ['DISCARD'] * 4

        line_new = f"{sam}\t{owner}\t{ror_id}\t"
        line_new += f"{institute}\t{country}\t{relation}\n"
        lines_new.append(line_new)

    print(f"n of samples: {len(lines_new)}")
    print(f"n of samples with corrected owner: {counter['owner']}")
    print(f"n of samples with international institute: {counter['international']}")
    with open(path_out, 'w') as f:
        f.writelines(lines_new)

