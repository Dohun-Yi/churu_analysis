import sys, os
# This code take some times (few minutes)

root = '../../'
workd = root + '/analysis/3_country_institute'
path_interest = workd + '/SRA_Accessions_with_expr.uniq.tab.filter_f25_dedup'
path_P2GSE = workd + '/bioproject.gse.interest.txt'
path_P2PMID = workd + '/bioproject.journals.interest.txt'
path_S2GSM = workd + '/biosample_set.clean.geo.txt'
path_S2Own = workd + '/biosample_set.clean.summary.interest.txt'
path_ENA = workd + '/ena_sample_human.table'
path_gsm = workd + '/parsed_GSE.txt'
path_out = workd + '/all_in_one.txt'

"""
Desired format:
    biosample bioproject gsm gse PMID Owner_Country Onwer Country Institute State City PostalCode
"""

def manual_imputation(d_ENA):
    # Below records were ommited and manually filled
    d_ENA["SAMEA9926956"]="Human Embryo and Stem Cell Laboratory"
    d_ENA["SAMEA9926954"]="Human Embryo and Stem Cell Laboratory"
    d_ENA["SAMEA9926953"]="Human Embryo and Stem Cell Laboratory"
    d_ENA["SAMEA9926949"]="Human Embryo and Stem Cell Laboratory"
    d_ENA["SAMEA6870306"]="Epithelial Carcinogenesis Group"
    d_ENA["SAMEA6870307"]="Epithelial Carcinogenesis Group"
    d_ENA["SAMEA6870308"]="Epithelial Carcinogenesis Group"
    d_ENA["SAMEA6870305"]="Epithelial Carcinogenesis Group"
    d_ENA["SAMEA6870310"]="Epithelial Carcinogenesis Group"
    d_ENA["SAMEA6870311"]="Epithelial Carcinogenesis Group"
    d_ENA["SAMEA6870309"]="Epithelial Carcinogenesis Group"
    return d_ENA

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
    d_P2GSE = makeD(path_P2GSE)
    d_P2PMID = makeD(path_P2PMID)
    d_S2GSM = makeD(path_S2GSM, interest=interest)
    d_S2Own = makeD(path_S2Own, 0, 3, interest=interest)
    d_S2ENA = get_ena_institute(path_S2Own, empty='EBI', interest=interest)
    d_ENA = makeD(path_ENA)
    d_gsm = get_whole_gsm(path_gsm)

    d_ENA = manual_imputation(d_ENA) # manual

    lines_out = []
    for s in d_S2P: # samples of interest
        p = d_S2P[s]
        url = d_S2GSM.get(s,'-')
        gsm = url.split('?acc=')[1] if 'acc=' in url else '-'
        gse = d_P2GSE.get(p, '-')
        pmid = d_P2PMID.get(p, '-')
        owner = d_S2Own.get(s, '-')
        if s.startswith('SAMEA'):
            owner = d_S2ENA.get(s,'EBI')
        country = '-'
        others = '\t'.join(d_gsm[gsm][2:]) if gsm in d_gsm else '-'
        line_out = f"{s}\t{p}\t{gsm}\t{gse}\t{pmid}\t{country}\t{owner}\t{others}\n"
        lines_out.append(line_out)

    with open(path_out,'w') as f:
        f.writelines(lines_out)
