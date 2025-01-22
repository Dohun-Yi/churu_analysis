import sys,os
import csv
import re

# samn --> proj --> pubmed --> impact

root = '../../'
workd = root + '/analysis/5_journals/'
path_in = workd + '/1697001520.compare_result.final.txt'
path_out = workd + '/sam2journal.txt'

path_sra = root + '/resources/SRA_Accessions_with_expr.uniq.tab.filter_f25_dedup' # samn_proj
path_bioproject = root + '/resources/bioproject.journals.interest.txt' # proj_pubmed
path_if2year = root + '/resources/bioproject.journals.interest.pubmedids.name.if2year.txt' # if2year

def get_samn_proj(path):
    samn_proj = {}
    for line in open(path):
        items = line.rstrip('\n').split('\t')
        samn_proj[items[17]] = items[18]
    return samn_proj

def get_proj_pubmed(path):
    proj_pubmed = {}
    for line in open(path):
        items = line.rstrip('\n').split('\t')
        if items[1].startswith('10.1101/'): # bioRxiv
            items[1] = 'bioRxiv'
        proj_pubmed[items[0]] = items[1].lstrip('-').split(',,')
    return proj_pubmed

def get_if2year(path):
    if2yearD, titleD = {}, {}
    for line in open(path):
        items = line.rstrip('\n').split('\t')
        pmids = [items[0].split(':')[-1].split(',,')[0]]
        dois = items[1].split(',,')
        issns = items[2].split(',,')
        titles = items[3].split(',,')
        if2year = items[5]
        for k in pmids+dois+issns+titles:
            if2yearD[k] = if2year
            titleD[k] = titles[0]
    return if2yearD, titleD

def is_float(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

if __name__=='__main__':
    samn_proj = get_samn_proj(path_sra)
    proj_pubmed = get_proj_pubmed(path_bioproject)
    if2yearD, titleD = get_if2year(path_if2year)
    # Manually update exception cases
    if2yearD.update({
        '10.5483': if2yearD['BMB reports'], # PRJNA294025, BMB reports
        '10.1158/0008-5472': if2yearD['10.1158/0008-5472.CAN-18-3962'], # PRJNA272556, Cancer research
        '10.1186/s12864-016-3135-y':if2yearD['BMC genomics'],# PRJNA305679, BMC Genomics
        'bioRxiv':'bioRxiv'})
    titleD.update({
        '10.5483': 'BMB reports', # PRJNA294025, BMB reports
        '10.1158/0008-5472': 'Cancer research', # PRJNA272556, Cancer research
        '10.1186/s12864-016-3135-y':'BMC genomics',# PRJNA305679, BMC Genomics
        'bioRxiv':'bioRxiv'})

    lines_new = []
    for line in open(path_in):
        line = line.rstrip('\n')
        items = line.split('\t')
        samn = items[0]
        bioproject = samn_proj[samn]
        pubmeds = proj_pubmed.get(bioproject,'')
        try:
            if2years = [if2yearD[p] for p in pubmeds if p != '']
            titles   = [  titleD[p] for p in pubmeds if p != '']
            tup_max = [(float(s),t) for s,t in zip(if2years,titles) if is_float(s)]
            if len(tup_max) > 0:
                if2years_max, title_max = max(tup_max, key=lambda x:x[0])
            else:
                if2years_max, title_max = '', ''
            line_new = '\t'.join([samn,  bioproject, ';'.join(pubmeds), ';'.join(titles), ';'.join(if2years), title_max, str(if2years_max)])
            lines_new.append(line_new + '\n')
        except Exception as e:
            print ('CHECK THIS LINE >>>>', line, bioproject, pubmeds, e)
            #  exit()
    with open(path_out, 'w') as f:
        f.writelines(lines_new)
