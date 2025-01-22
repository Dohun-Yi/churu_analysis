#!/usr/bin/env python
# coding: utf-8

import sys, os
import re
import random
import time
import json
from thefuzz import fuzz, process
from langchain_openai import ChatOpenAI

if len(sys.argv) != 3:
    exit("python ~.py n_idx n_batch")
n_idx = int(sys.argv[1])
n_batch = int(sys.argv[2])

root = '../../'
workd = root + 'analysis/1_check_output'
path_cello_reform = workd + '/cellosaurus.table.reform'
path_title = workd + '/../../resources/biosample_set.clean.titles.filt.txt'
path_bioraw = workd + '/../../resources/biosample_set.clean.summary.interest.filt.txt'

path_outdir = root + '/analysis/2_cell_name_guess/'
path_out_full = path_outdir + '/guess/cellname_guesse_full.txt'
path_out = path_outdir + '/guess/cellname_guesse.txt'

# Setup openAI API
#  os.environ["OPENAI_API_KEY"] = ""
model = ChatOpenAI(model="gpt-4o-mini", temperature=0.0)


# # Parse Cellosaurus Name
def read_cello_test(path):
    parentD, problemD, synD = {}, {}, {}
    allnames = set()
    syn2regular = {}
    regular2syn = {}
    for line in open(path):
        items = line.rstrip('\n').split('\t')
        child = items[0]
        parent = items[2]
        problem, synlist = items[4], items[5]
        
        allnames.add(child)
        syn2regular.setdefault(child,[]).append(child)
        for syn in synlist.split('; '):
            if syn not in set(['-','']):
                allnames.add(syn)
                regular2syn.setdefault(child,[]).append(syn)
                syn2regular.setdefault(syn,[]).append(child)
    # add itself
    allnames_clean = {}
    for name in list(allnames):
        name_clean = re.sub('\W','',name).upper()
        allnames_clean.setdefault(name_clean,set()) 
        allnames_clean[name_clean] |= set(syn2regular[name])
    return allnames_clean, allnames, syn2regular, regular2syn

def get_candidates(name_clean, allnames_clean_dict, s2r, r2s):
    names = list(allnames_clean_dict[name_clean])
    names_regular = set()
    for n in names:
        names_regular |= set(s2r[n])
    return [[ns, r2s.get(ns,'-')] for ns in list(names_regular)]

allnames_clean_dict, _, s2r, r2s = read_cello_test(path_cello_reform)
allnames_clean_set = set(allnames_clean_dict.keys())


# # Parse Biosample Documents
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

titleD = read_title(path_title)
rawD = read_bioraw(path_bioraw)

docD = {}
for key in rawD.keys():
    converted_doc = f"sample id: {key}" + '\n'
    converted_doc += '\n'.join([f"{field}: {value}" for field, value in zip(['title', 'sample name'], titleD[key])]) + '\n'
    converted_doc += '\n'.join([f"{field}: {value}" for field, value in rawD[key]])
    docD[key] = converted_doc


# # Run OpenAI API
def guess_name(myid):
    prompt1 = """
    Question: "From below context, try to find unique identifier of cell line name and return. remove any descriptive words other than name. don't explain."
    Given context to you: "%s"
    Cell line name:
    """
    
    prompt2 = """
    Question: "From the list below, select the option that is most likely to match the given cell line. Use very strict criteria for similarity, considering exact matches or very close names only. give only an alphabet, not the contents. always check full document for accuracy."
    Given cell line to you: "%s"
    Given options: "A. not in this list\n%s"
    Full document: "%s"
    Answer:
    """

    cleaner = lambda text: re.sub('\W','',text)
    
    # (1). extract cell line name from BioSample document
    guess_gpt_1 = model.invoke(prompt1 % docD[myid]).content
    # (2-1). find similar cell line name from Cellosaurus DB
    similar5 = process.extract(cleaner(guess_gpt_1), allnames_clean_set, scorer=fuzz.ratio, limit=5)
    # (2-2). remove redundant ID and select four
    candidates = []
    for sim, score in similar5: # remove redundant keys
        _tmps = get_candidates(sim, allnames_clean_dict, s2r, r2s)
        for _tmp in _tmps:
            if _tmp[0] not in [c[0] for c in candidates]:
                candidates.append(_tmp + [score])
    candidates = candidates[:4]
    # (3). select the cell line name most liekly to match
    choices = {'BCDE'[i] : x[0] for i, x in enumerate(candidates)}
    choices['A'] = 'none of these'
    context = '\n'.join(['BCDE'[i] + '.' + x[0] + ' (synonyms: %s)' % ', '.join(x[1]) for i, x in enumerate(candidates)])
    guess_gpt_2 = model.invoke(prompt2 % (guess_gpt_1, context, docD[myid])).content
    final_choice = choices.get(guess_gpt_2, 'ERROR')
    
    log_data = {
        "id": myid,
        "prompt1": prompt1 % docD[myid],
        "gpt-4o-mini-response-1": guess_gpt_1,
        "cellosaurus-candidates": str(candidates),
        "prompt2": prompt2 % (guess_gpt_1, context, docD[myid]),
        "gpt-4o-mini-response-2": guess_gpt_2,
        "final-choice": final_choice
    }
    return (final_choice, log_data)


# Function test

# myid='SAMN14686168' # where cell line is specified correctly
# myid='SAMN22236237' # in title
# myid='SAMN11893676' # testing
# myid='SAMN11893676' # problematic case > CVCL may solve it?
#myid = random.choice(list(docD.keys())) # random
#guess_name(myid)


# Run
buff_out = ''
buff_json = {}
total_records = len(docD)

max_retry=10
retry_delay=5

for i, myid in enumerate(docD):
    if i % n_batch != n_idx: continue
    for attempt in range(max_retry):
        try:
            final_choice, log_data = guess_name(myid)
            buff_out += f"{myid}\t{final_choice}\n"
            buff_json[myid] = log_data
            break
        except Exception as e:
            time.sleep(retry_delay)

# last batch
with open(f"{path_out_full}_{n_idx}_{n_batch}", 'w') as f:
    json.dump(buff_json, f, indent=4)
with open(f"{path_out}_{n_idx}_{n_batch}", 'w') as f:
    f.writelines(buff_out)
