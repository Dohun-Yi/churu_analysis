import sys, os
import re

root = '../../'
workd = root + '/analysis/1_check_output/'
batch = '1697001520'

path_table = workd + 'SRA_churu_output_compare.%s.categorize.txt' % batch
path_churu = workd + '/%s.churu.out' % batch
path_model = root + '/reference/CCLE/22Q4/Model.tsv'
path_sra = root + '/resources/SRA_Accessions_with_expr.uniq.tab.filter_f25_dedup'
path_sam = root + '/resources/biosample_set.clean.summary.interest.txt'
path_out = workd + '/1697001520.parse.table'

def stratify(name):
    if name == 'T-1':
        return 'TDASH1'
    elif name == 'H-4':
        return 'HDASH4'
    else:
        return re.sub(r'[^a-zA-Z0-9\s]', '', name).replace(' ', '').upper()

def read_churu(path):
    with open(path) as f:
        lines = f.readlines()
    srr_now = lines[0].rstrip('\n')

    srr_match = {}
    for line in lines[1:]:
        line = line.rstrip('\n')

        if line == '':
            srr_now = ''
        elif srr_now == '':
            srr_now = line
        elif line.startswith('cell_id'):
            continue
        else:
            items = line.split('\t')
            if float(items[5]) >= 0.5: # posterior
                srr_match.setdefault(srr_now,[]).append(items)
    return srr_match

def read_model(path):
    with open(path) as f:
        lines = f.readlines()
    cell_patient = {}
    for line in lines[1:]:
        items = line.split('\t')
        cell = items[3]
        patient = items[1]
        cell_patient[cell] = patient
    return cell_patient

def read_sra(path):
    srr_sam = {}
    with open(path) as f:
        lines = f.readlines()
    for line in lines:
        items = line.rstrip('\n').split('\t')
        srr = items[0]
        sam = items[17]
        srr_sam[srr] = sam
    return srr_sam
    
def read_sam(path):
    #  pattern = r"^cell[ _\-]?line$"
    pattern = r"cell[ _\-]?line"
    contains_pattern = lambda t: re.search(pattern, t[0])
    sam_info = {}
    with open(path) as f:
        lines = f.readlines()
    for line in lines:
        items = line.rstrip('\n').split('\t')
        sam = items[0]

        if len(items) == 6:
            items[3] = ' '.join([items[3], items[4]])
            items[4] = items[5]

        cellinfos = [x.split('==EQUAL==') for x in items[4].split(';;DELIMITER;;')]
        #  cellinfos = filter(contains_pattern, cellinfos)

        info = [items[2]] # date
        info += [items[3]] # organization
        info += [cellinfos]
        sam_info[sam] = info
    return sam_info

if __name__=='__main__':
    print 'reading churu output...'
    srr_match = read_churu(path_churu) # can have multiple match
    print 'reading reference model...'
    cell_patient = read_model(path_model) # only patient
    print 'reading sra list...'
    srr_sam = read_sra(path_sra) # only sam
    print 'reading sample list...'
    sam_info = read_sam(path_sam) # can have multiple cell info
    print 'reading done.'

    with open(path_out, 'w') as fo:
        with open(path_table) as fi:
            lines = fi.readlines()
        for line in lines:
            items = line.rstrip('\n').split('\t')[:10]
            srr = items[0]
            sam = srr_sam[srr]
            info = sam_info[sam]
            date = info[0]
            institute = info[1]
            cells_ori = info[2]
            matches = srr_match.get(srr,[])
            try:
                lst = items
                lst += [sam, date, institute]
                lst += [';'.join(['%s:%s' % (stratify(x), stratify(y)) for x,y in cells_ori])]
                lst += ['-'] if len(matches) == 0 else \
                        [';'.join(['|'.join(x) for x in matches])]
                line_new = '\t'.join(lst)
                #  print (line_new)
                fo.write(line_new + '\n')
            except Exception as e:
                print(line, sam, cells_ori)
                exit()



