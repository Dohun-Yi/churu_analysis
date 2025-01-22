import sys, os
import xml.etree.ElementTree as ET
import re

path_root = '../'
path_workd = path_root + '/resources/'
path_interest = path_workd + '/biosample_set.clean.summary.interest.txt'
path_sra = path_workd + '/SRA_Accessions.tab'
path_metadir = path_workd + '/NCBI_SRA_Metadata_Full_20230821/'
path_out = path_workd + '/SRA_Accessions_with_expr.tab'

def read_interest(path):
    # no header
    # need to remove attribute fields
    get_k = lambda x: x.split('==EQUAL==')[0].lower()
    get_v = lambda x: x.split('==EQUAL==')[1].upper()
    samnD = {}
    with open(path) as f:
        lines = f.readlines()
    for line in lines:
        items = line.rstrip('\n').split('\t')
        items[-1] = re.sub(' |-|_','', items[-1])
        attributes = {get_k(x):get_v(x) for x in items[-1].split(';;DELIMITER;;')}
        items[-1] = attributes.get('cellline','-')
        samn = items[0]
        samnD[samn] = items
    return samnD


def read_sra(path):
    samn_srx = {}
    srx_srr = {}
    srrD = {}
    with open(path) as f:
        lines = f.readlines()
    for line in lines[1:]: # has header
        items = line.rstrip('\n').split('\t')
        if items[6] != 'RUN': continue
        srr, srx, samn = items[0], items[10], items[17]
        if samn == '-':
            continue
        samn_srx.setdefault(samn, set()).add(srx)
        srx_srr.setdefault(srx, []).append(srr)
        srrD[srr] = items
    return samn_srx, srx_srr, srrD


def read_meta(srx_list, srx_srr, srrD, path_dir):
    metaD = {}
    n_total = len(srx_list)
    i = 0
    for srx in srx_list:
        #  if i % 100 == 0:
            #  print('\tprocessing %i / %i...' % (i, n_total))
        i += 1
        sra_list = list(set(map(lambda x: srrD[x][1], srx_srr[srx])))
        if len(sra_list) != 1:
            print('WARNING %s, %s' % (srx, sra_list))
        sra = sra_list[0]
        path_xml = path_dir + '/%s/%s.experiment.xml' % (sra, sra)
        if not os.path.isfile(path_xml):
            metaD[srx] = ['NULL','NULL','NULL','NULL','NULL']
            continue
        tree = ET.parse(path_xml)
        root = tree.getroot()
        for ee in root.findall(".//EXPERIMENT"):
            # Extract information
            experiment_accession = ee.get("accession")
            if experiment_accession != srx:
                continue
            library_strategy = ee.find(".//LIBRARY_STRATEGY").text
            library_source = ee.find(".//LIBRARY_SOURCE").text
            library_selection = ee.find(".//LIBRARY_SELECTION").text
            library_layout = ee.find(".//LIBRARY_LAYOUT/*").tag
            platform = ee.find(".//INSTRUMENT_MODEL").text
            metaD[srx] = [library_strategy, \
                    library_source, library_selection, library_layout, platform]
    return metaD

if __name__=='__main__':
    if len(sys.argv) != 1:
        print('usage: python ~.py')
        exit()
    
    print('reading interest...')
    samnD = read_interest(path_interest)
    print('reading sra...')
    samn_srx, srx_srr, srrD = read_sra(path_sra)
    
    srx_list = []
    for samn in samnD.keys():
        srx_list += list(samn_srx.get(samn, []))
    num_total_srx = len(srx_list)

    print('writing...')
    lines = []
    i = 0

    srx_done = set([])
    if os.path.isfile(path_out):
        with open(path_out) as f:
            lines = f.readlines()
            srx_done = set(map(lambda x: x.split('\t')[10], lines))

    i_srx = 0
    for samn in samnD.keys():
        for srx in list(samn_srx.get(samn, [])):
            i_srx += 1
            if i_srx % 1000 == 0:
                print('processing srx %i / %i...' % (i_srx, num_total_srx))
            if srx in srx_done:
                continue
            metaD = read_meta([srx], srx_srr, srrD, path_metadir)
            for srr in srx_srr.get(srx, []):
                i += 1
                value_srr = srrD.get(srr, ['NULL'] * 20)
                value_meta = metaD.get(srx, ['NULL','NULL','NULL','NULL'])
                lines.append('\t'.join(value_srr + value_meta) + '\n')
                if len(lines) > 10000:
                    print 'saving results on iteration %i' % i
                    with open(path_out, 'a') as f:
                        f.writelines(lines)
                    lines = []
