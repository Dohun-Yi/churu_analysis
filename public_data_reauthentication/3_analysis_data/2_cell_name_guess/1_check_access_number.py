import sys, os
import re
import json
import ahocorasick

# Input and output files
root = '../../'
workd = root + 'analysis/2_cell_name_guess/'
path_cello_reform_acc = workd + '/cellosaurus.table.reform.acc'
path_title = workd + '/biosample_set.clean.titles.filt.txt'
path_bioraw = workd + '/biosample_set.clean.summary.interest.filt.txt'
path_out = workd + '/cellname_by_access_number.txt'
path_log = workd + '/cellname_by_access_number.txt.log'

# Readers
def read_bioraw(path):
    rawD = {}
    for line in open(path):
        items = line.rstrip('\n').split('\t')
        biosam = items[0]
        info_raw = items[-1] # from 4 to -1 for error handling
        rawD[biosam] = [x.split('==EQUAL==') for x in info_raw.split(';;DELIMITER;;')]
    return rawD
def read_title(path):
    titleD = {}
    for line in open(path):
        items = line.rstrip('\n').split('\t')
        if len(items) == 2:
            items.append('-')
        elif len(items) > 3:
            print(items)
            exit()
        samn, title, alias = items
        titleD[samn] = (title, alias)
    return titleD

# Parser
def make_doc(titleD, rawD):
    docD = {}
    for key in rawD.keys():
        converted_doc = f"sample id: {key}" + '\n'
        converted_doc += '\n'.join([f"{field}: {value}" for field, value in zip(['title', 'sample name'], titleD[key])]) + '\n'
        converted_doc += '\n'.join([f"{field}: {value}" for field, value in rawD[key]])
        docD[key] = converted_doc
    return docD

# Define Target DB
def get_target_db():
    # CelloSaurus cross-reference list https://www.cellosaurus.org/description.html
    # Cell line catalogs/collections
    lst1 = "Abcam, Abeomics, ABM, AddexBio, ATCC, BCRC, BCRJ, BEI Resources, CancerTools, CBA, CCLV, CCTCC, Cell Biolabs, CLS, Coriell, DGRC, DiscoverX, DSHB, DSMZ, ECACC, FCDI, GeneCopoeia, HIV Reagent Program, Horizon Discovery, Hysigen, IBRC, ICLC, Imanis, Innoprot, IZSLER, JCRB, KCB, KCLB, Kerafast, KYinno, Millipore, MMRRC, NCBI_Iran, NCI-DTP, NHCDR, NIHhESC, NISES, NRFC, PerkinElmer, RCB (Riken), RIKEN_BRC_EPD, Rockland, RSCB, TKG, TNGB, Ubigene, WiCell, Ximbio"
    lst1 = lst1.split(', ')
    # Cell line databases/resources
    lst2 = "cancercelllines, CCRID, Cell_Model_Passport, CGH-DB, CLDB, ColonAtlas, Cosmic-CLP, dbMHC, DepMap, DSMZCellDive, ESTDAB, FCS-Free, FlyBase_Cell_line, GDSC, hPSCreg, ICLDB, IGRhCellID, IHW, IPD-IMGT/HLA, ISCR, LIMORE, LINCS_HMS, LINCS_LDP, Lonza, SKIP, SKY/M-FISH/CGH, SLKBase, TOKU-E"
    lst2 = lst2.split(', ')
    # below DB have duplicated names, either internallly or externally
    lst3 = ['CCTCC', 'hPSCreg', 'SKIP', 'Ubigene', 'NCBI_Iran', 'Innoprot', 'CLDB', 'ColonAtlas', 'Imanis']
    # these DB has IDs consist of short pure number
    lst4 = ['CGH-DB','BCRJ', 'ISCR','Lonza', 'TOKU-E', 'IPD-IMGT/HLA', 'BCRC', 'KCLB', 'SLKBase', 'SKY/M-FISH/CGH', 'dbMHC']
    target_db = (set(lst1) | set(lst2)) - set(lst3) - set(lst4)
    return target_db

# test for access number
def read_cello(path, target_db):
    accD = {}
    sources = {}
    for line in open(path):
        items = line.rstrip('\n').split('\t')
        child = items[0]
        acc_list = items[9].split(';')
        for acc in acc_list:
            if acc=='': continue
            try:
                source, cell = acc.split(':',1)
                cell_clean = re.sub("[_,-]"," ",cell).upper()
                if source not in target_db: continue
                if cell_clean in accD and child != accD[cell_clean]:
                    exit(f'ERROR: {source} {cd} is {child} and in accDc from source {sources[cell_clean]} ({accD[cell_clean]})')
                sources.setdefault(cell_clean,[]).append(source)
                accD[cell_clean] = child
            except Exception as e:
                exit(f"ERROR:{acc}, {e}")
    return accD, sources

if __name__=='__main__':

    titleD = read_title(path_title)
    rawD = read_bioraw(path_bioraw)

    target_db = get_target_db()
    docD = make_doc(titleD, rawD)
    accD, sources = read_cello(path_cello_reform_acc, target_db)

    # Create an Aho-Corasick automaton
    keywords = accD.keys()
    A = ahocorasick.Automaton()
    for idx, keyword in enumerate(keywords):
        A.add_word(keyword, (idx, keyword))
    A.make_automaton()

    # Search for keywords in the text
    matchcounter = {True:0, False:0}
    lines_out = []
    lines_log = []
    for biosample, text in docD.items():
        hasmatch = False
        text_clean = re.sub("[_,-]"," ",text).upper()
        unique_matches = set()
        unique_cell_name = set()
        # Iterate all search results
        for end_index, (idx, acc) in A.iter(text_clean):
            start_index = end_index - len(acc) + 1
            state='pass'
            if bool(re.match("\S", text[start_index-1])):
                state='skip-start'
            elif end_index+1 != len(text) and bool(re.match("\S", text[end_index+1])):
                state='skip-end'
            else:
                hasmatch = True
                unique_cell_name.add(accD[acc])
                unique_matches.add(f"{accD[acc]}:::{acc}:::{sources[acc]}")
                subcontext = text[start_index-15:end_index+15].replace('\n',' ')
                line_log = f">{biosample}: [{state}] Found an ID '{acc}' of cell '{accD[acc]}' (source:{sources[acc]} from index {start_index} to {end_index}, context: '{subcontext}'\n"
                lines_log.append(line_log)
        line = f"{biosample}\t" + ';;;'.join(list(unique_cell_name)) + "\t" + ';;;'.join(list(unique_matches)) + '\n'
        lines_out.append(line)
        matchcounter[hasmatch] += 1

    with open(path_out, 'w') as f:
        f.writelines(lines_out)
    with open(path_log, 'w') as f:
        f.writelines(lines_log)

    print(f"Successfully finished")
    print(f"number of all documents: {len(docD)}")
    print(f"at least one access number found: {matchcounter[True]}")
    print(f"no access number found: {matchcounter[False]}")
