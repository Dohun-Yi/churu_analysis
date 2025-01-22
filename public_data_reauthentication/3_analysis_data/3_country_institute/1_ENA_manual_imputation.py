import sys, os
# This code take some times (few minutes)

root = '../../'
workd = root + '/analysis/3_country_institute'
path_interest = workd + '/SRA_Accessions_with_expr.uniq.tab.filter_f25_dedup'
path_S2Own = workd + '/biosample_set.clean.summary.interest.txt'
path_ENA = workd + '/ena_sample_human.table'
path_out = workd + '/pre_imputation_biosample'

"""
    from output of this code, i found that below samples has their
    center name omitted, and manually found each center name
    from https://www.ebi.ac.uk/
    some of them had proper biosample record, but omited in the middle.
    - it was omitted on xml

SAMEA9926956    Human Embryo and Stem Cell Laboratory
SAMEA9926954    Human Embryo and Stem Cell Laboratory
SAMEA9926953    Human Embryo and Stem Cell Laboratory
SAMEA9926949    Human Embryo and Stem Cell Laboratory
SAMEA6870306    Epithelial Carcinogenesis Group
SAMEA6870307    Epithelial Carcinogenesis Group
SAMEA6870308    Epithelial Carcinogenesis Group
SAMEA6870305    Epithelial Carcinogenesis Group
SAMEA6870310    Epithelial Carcinogenesis Group
SAMEA6870311    Epithelial Carcinogenesis Group
SAMEA6870309    Epithelial Carcinogenesis Group


"""

def makeD(path, idx1=0, idx2=1, empty='-', interest=False):
    print(f"making dictionary from {path}")
    d = {}
    with open(path) as f:
        lines = f.readlines()
    for line in lines:
        items = line.rstrip('\n').split('\t')
        k = items[idx1]
        if interest and k not in interest:
            continue
        v = items[idx2] if len(items) > idx2 else empty
        d[k] = v
    return d

def get_ena_institute(path, idx1=0, idx2=None, empty='-', interest=False):
    print(f"parsing ENA institute from {path}")
    d = {}
    with open(path) as f:
        lines = f.readlines()
    for line in lines:
        items = line.rstrip('\n').split('\t')
        k = items[idx1]
        if not k.startswith('SAMEA'):
            continue
        if interest and k not in interest:
            continue
        infoL = [x.split('==EQUAL==') for x in items[4].split(';;DELIMITER;;')]
        infoD = {x1:x2 for x1, x2 in infoL}
        v = infoD.get('INSDC center name', empty)
        d[k] = v
    return d

def get_whole_gsm(path):
    d = {}
    for line in open(path):
        items = line.rstrip('\n').split('\t')
        d[items[0]] = items
    return d

if __name__=='__main__':
    d_S2P = makeD(path_interest,17,18)
    interest = set(d_S2P.keys())
    d_S2ENA = get_ena_institute(path_S2Own, empty='-', interest=interest)
    d_ENA = makeD(path_ENA, empty='-', interest=interest)

    lines_out = []
    for s in d_S2P: # samples of interest
        if s.startswith('SAMEA'):
            owner1 = d_S2ENA.get(s,'-')
            owner2 = d_ENA.get(s,'-')
            if owner1 != '-':
                continue
            line_out = f"{s}\t{owner1}\t{owner2}\n"
            lines_out.append(line_out)

    with open(path_out,'w') as f:
        f.writelines(lines_out)
