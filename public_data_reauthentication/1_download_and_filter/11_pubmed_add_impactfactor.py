import sys, os
import csv
import re

root = '../'
workd = root + '/resources/'
path_if = workd + '/scimagojr_2022.csv'
path_in = workd + '/bioproject.journals.interest.pubmedids.name.txt'
path_out = workd + '/bioproject.journals.interest.pubmedids.name.if2year.txt'

def stratify_journal(text):
    text = text.lower().lstrip('the')
    text = text.replace('&amp;','').replace('&','').replace(' and ','')
    text = re.sub(r'[^a-zA-Z0-9\s]', '', text).replace(' ', '')
    return text

def read_if(path):
    issnD = {}
    titleD = {}
    with open(path) as csvfile:
        reader = csv.reader(csvfile, delimiter=';')
        next(reader) # header
        for row in reader:
            title = stratify_journal(row[2])
            issns = row[4].split(', ')
            totalcites3 = float(row[11])
            citabledoc3 = float(row[12])
            if2year = totalcites3 / citabledoc3 if citabledoc3 > 10 else 0 # actually it is 3 year
            titleD[title] = str(if2year)
            for issn in issns:
                if issn in issnD:
                    issnD[issn] = str(max(issnD[issn], float(if2year)))
                else:
                    issnD[issn] = str(if2year)
    return issnD, titleD

if __name__=='__main__':
    issnD, titleD = read_if(path_if)

    lines_new = []
    for line in open(path_in):
        items = line.rstrip('\n').split('\t')
        issn = items[2].replace('-','')
        title = stratify_journal(items[3])

        line_new = line.rstrip('\n') + '\t'
        if issn == '':
            line_new += 'preprint'
        elif issn in issnD:
            line_new += issnD[issn]
        elif title in titleD:
            line_new += titleD[title]
        else:
            line_new += 'N.A'

        lines_new.append(line_new + '\n')
    with open(path_out, 'w') as f:
        f.writelines(lines_new)
