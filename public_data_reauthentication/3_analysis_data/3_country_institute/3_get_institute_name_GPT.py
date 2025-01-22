# import packages and specify input files
import pysqlite3, sys
sys.modules["sqlite3"] = sys.modules.pop("pysqlite3")
import sys, os
import numpy as np
from thefuzz import fuzz, process
import random
import json

if len(sys.argv) != 3:
    exit("python ~.py n_idx n_batch")
n_idx = int(sys.argv[1])
n_batch = int(sys.argv[2])

root = '../../'
workd = root + '/analysis/3_country_institute'
path_in = workd + '/all_in_one.txt'
path_ror = root + '/resources/ROR/v1.50-2024-07-29-ror-data.filt.tsv' # the version that "withdrawn" records are removed
path_out_full = workd + '/select_gpt/institute_selected_full.txt'
path_out = workd + '/select_gpt/institute_cellname_guesse.txt'

# Process input file
# unique owner information is required for candidate selection
# and whole record is required for contextual information about the sample (like GEO contact information)

owner2items = {}
owner2queries = {}
query2owners = {}
query2has_perfect_match = {}
owners_target = set()
sam2candidates = {}
sam2candidates_1 = {}
sam2candidates_2 = {}
sam2candidates_3 = {}
sam2items = {}
for line in open(path_in):
    items = line.rstrip('\n').split('\t')
    sam, owner = items[0], items[6]
    queries = [owner] + list(map(lambda x: x.strip(' '), owner.split(',')))
    sam2items[sam] = items
    sam2candidates[sam] = [] # Empty for now
    sam2candidates_1[sam] = []
    sam2candidates_2[sam] = []
    sam2candidates_3[sam] = []
    owner2items.setdefault(owner, []).append(items)
    owner2queries[owner] = queries
    owners_target.add(owner)
    for query in queries:
        query2owners.setdefault(query,set()).add(owner)
        query2has_perfect_match[query] = False

# Process reference file (ROR - filtered to remove withdrawn records)
# these records will be used for...
#   1. uncased perfect match
#   2. match based on levenshtein distance(theffuzz)
#   3. match based on semantic similarity (embedding vector)

# name can't be used as key, because there are duplications
any_name2ids = {} # str:lst
id2doc = {}      # str:str
id2items = {}    # str:lst

# aliases: seperated by "; "
# labels: seperated by "; ". and starts with language like "en: ", "fr: "
# acronyms: seperated by "; "
for line in open(path_ror):
    if line.startswith('id\t'): continue # header
    items = line.rstrip('\n').split('\t')
    ror_id, name, aliases, labels, acronyms = items[0], items[1], items[5], items[6], items[7]
    city, region, country = items[12], items[16], items[19]
    
    doc = f"name:{name}\naliases:{aliases}\nlabels:{labels}\nacronyms:{acronyms}\ncountry:{country}\nregion:{region}\ncity:{city}"
    id2doc[ror_id] = doc
    id2items[ror_id] = items

    all_names = [name]
    all_names += aliases.split('; ')
    all_names += list(map(lambda x: x.split(': ')[-1], labels.split('; ')))
    all_names += acronyms.split('; ')
    for n in all_names:
        if n == '': continue
        any_name2ids.setdefault(n,[]).append(ror_id)

# Setup OpenAI Embedding, vectorDB
from langchain_community.vectorstores import Chroma
from langchain_openai import OpenAIEmbeddings
#  os.environ["OPENAI_API_KEY"] = ""
embedding = OpenAIEmbeddings(model='text-embedding-3-large')

# DO NOT MAKE DB AGAIN: this one takes a lot of times!!, it is already final (removed withdrawn records)
db_path = "./chroma_db_openAIembed3large_ror_v150_latest"
if not os.path.exists(db_path):
    print('this code should not be executed...')
    db = Chroma.from_texts(list(id2doc.values()), embedding, persist_directory=db_path) # text-embedding-3-large
else:
    print(f'Found chroma DB {db_path}.')
    db = Chroma(persist_directory=db_path, embedding_function=embedding)


### Candidate 1 - uncased perfect match
uncase2ids = {} # str(upper): set
for name, ror_id_list in any_name2ids.items():
    name_upper = name.upper()
    uncase2ids.setdefault(name_upper, set())
    for ror_id in ror_id_list:
        uncase2ids[name_upper].add(ror_id)

print(f"n before uncase : {len(any_name2ids)}")
print(f"n after uncase  : {len(uncase2ids)}")

# counter = {'owner':0, 'sample':0}
queries_all = [q for q, m in query2has_perfect_match.items()]
print(f"n quries not matched yet - before: {len(queries_all)} / {len(query2has_perfect_match)}")

for query in queries_all:
    if query.upper() not in uncase2ids:
        continue
    query2has_perfect_match[query] = True
    ror_id_list = uncase2ids[query.upper()]
    owners = query2owners[query]
    for owner in owners:
        sams = [items[0] for items in owner2items[owner]]
        for sam in sams:
            sam2candidates[sam] += ror_id_list
            sam2candidates_1[sam] += ror_id_list

counter = {}
counter['query'] = len([v for v in query2has_perfect_match.values() if v])
counter['owner'] = len([o for o, ql in owner2queries.items() if any([query2has_perfect_match[q] for q in ql])])
counter['sample'] = len([v for v in sam2candidates.values() if len(v) > 0])

print(f"n queries with uncased perfect match: {counter['query']} / {len(query2has_perfect_match)}")
print(f"n owners with uncased perfect match: {counter['owner']} / {len(owners_target)}")
print(f"n samples with uncased perfect match: {counter['sample']} / {len(sam2candidates)}")


### Candidate 2 - match based on levenshtein distance (will make a lot of poor matches)

import multiprocessing, time
from tqdm import tqdm

def fuzzy_search(query, docs, limit=10):
    return (query, process.extract(query, docs, scorer=fuzz.token_set_ratio, limit=10))

def fuzzy_runner(queries, docs, cpus=1):
    PROCESSES = cpus
    with multiprocessing.Pool(PROCESSES) as pool:
        params = [(query, docs, ) for query in queries]
        results = [pool.apply_async(fuzzy_search, p) for p in params]
        output = [r.get() for r in tqdm(results)]
    return output

# This takes like 10 minutes
import pickle
path_save_fuzz = './path_save_fuzz.pickle'

# Make query
docs = list(id2doc.values())
# queries_not_m = [q for q, m in query2hasmatch.items() if not m] # IS IT OK TO REMOVE ALL PERFECT MATCHING QUERIES??? NOT SURE
# print(f"n quries not matched yet - before: {len(queries_not_m)} / {len(query2hasmatch)}")

# Make result, or load premade files
if os.path.exists(path_save_fuzz):
    print(f"pickle for fuzz results found {path_save_fuzz}")
    with open(path_save_fuzz, 'rb') as f:
        results_fuzz = pickle.load(f)
else:
    print(f"pickle for fuzz results not found, making new one...")
    results_fuzz = fuzzy_runner(queries_all, docs, cpus=48)
    with open(path_save_fuzz, 'wb') as f:
        pickle.dump(results_fuzz, f)

doc2id = {v:k for k,v in id2doc.items()}

for query, result in results_fuzz:
    owners = query2owners[query]
    ror_id_list = [doc2id[doc] for doc, _ in result]
    ror_id_and_score_list = [[doc2id[doc], score] for doc, score in result] # debug
    for owner in owners:
        sams = [items[0] for items in owner2items[owner]]
        for sam in sams:
            sam2candidates[sam] += ror_id_list
            sam2candidates_2[sam] += ror_id_and_score_list # debug

### Candidate 3 - match based on embedding
from concurrent.futures import ThreadPoolExecutor # Using this because vectorDB can not be pickled (so, can not run by multiprocessing.Pool)
import time, warnings
from tqdm import tqdm
warnings.filterwarnings('ignore')

def embedding_search(query, limit=10):
    return (query, db.similarity_search_with_relevance_scores(query, k=10))

def embedding_runner(queries, db, PROCESSES=1):
    with ThreadPoolExecutor(max_workers=PROCESSES) as executor:
        futures = {executor.submit(embedding_search, q): q for q in queries}
        output = [future.result() for future in tqdm(futures, total=len(queries))]
    return output

# This takes like 20 minutes
import pickle
path_save_embedding = './path_save_embedding.pickle'

# Make query
# ignore perfectly matches keywords
# queries_not_m = [q for q, m in query2hasmatch.items() if not m]
# print(f"n quries not matched yet - before: {len(queries_not_m)} / {len(query2hasmatch)}")

# Make result, or load premade files
if os.path.exists(path_save_embedding):
    print(f"pickle for embedding results found {path_save_embedding}")
    with open(path_save_embedding, 'rb') as f:
        results_embedding = pickle.load(f)
else:
    print(f"pickle for embedding results not found, making new one...")
    results_embedding = embedding_runner(queries_all, db, PROCESSES=48)
    with open(path_save_embedding, 'wb') as f:
        pickle.dump(results_embedding, f)

doc2id = {v:k for k,v in id2doc.items()}

for query, result in results_embedding:
    owners = query2owners[query]
    ror_id_list = [doc2id[doc.page_content] for doc, _ in result]
    ror_id_and_score_list = [[doc2id[doc.page_content], score] for doc, score in result] # debug
    for owner in owners:
        sams = [items[0] for items in owner2items[owner]]
        for sam in sams:
            sam2candidates[sam] += ror_id_list
            sam2candidates_3[sam] += ror_id_and_score_list # debug

### Utils
def crawler(child2parents, child):
    if child in child2parents:
        parents = child2parents[child]
        hierarchy = []
        for parent in parents:
            hierarchy += [parent] + crawler(child2parents, parent)
        return hierarchy
    else:
        return []
        
def show_candidates(mysam, candidates, has_score=False, topN=None):
    owner = sam2items[mysam][6]
    conv = lambda x: id2items[x][1]
    scores = ['-'] * len(candidates)
    if has_score:
        candidates.sort(key=lambda x: x[1], reverse=True)
        scores = [c[1] for c in candidates]
        candidates = [c[0] for c in candidates]
        
    for i, candidate in enumerate(candidates[:topN]):
        parents_of_candidate = crawler(child2parents, candidate)
        candidate_name = conv(candidate)
        parents_of_candidate = list(map(lambda x: conv(x), parents_of_candidate))
        print(f"{scores[i]}\t{owner}:\t{candidate}\t{candidate_name}\t{parents_of_candidate}")

child2parents = {}
for child_id, items in id2items.items():
    if "Parent" not in items[-1]:
        continue
    parent_ids = [x for x in items[-1].split('; ') if x.startswith('Parent')][0]
    child2parents[child_id] = parent_ids.lstrip('Parent: ').split(', ')


### Find matching insititue from candidates
from langchain_openai import ChatOpenAI
model = ChatOpenAI(model="gpt-4o-mini", temperature=0.0)

def make_options(sam, sam2candidates_1, sam2candidates_2, sam2candidates_3, id2items, topN=None):
    conv = lambda x: id2items[x][1]
    c1 = sam2candidates_1[sam]
    c2 = [ror_id for ror_id, score in sorted(sam2candidates_2[sam], key=lambda x:x[1], reverse=True)][:topN]
    c3 = [ror_id for ror_id, score in sorted(sam2candidates_3[sam], key=lambda x:x[1], reverse=True)][:topN]
    ror_id_list = list(set(c1 + c2 + c3))
    options = [("-","none of these")]
    ror_id_list = sorted(ror_id_list, key=lambda x: conv(x)) # debug, just for visual inspection
    for ror_id in ror_id_list:
        items = id2items[ror_id]
        name, aliases, labels, acronyms = items[1], items[5], items[6], items[7]
        city, region, country = items[12], items[16], items[19]
        parents_of_candidate = crawler(child2parents, ror_id)
        
        aka = aliases.split('; ') 
        aka += list(map(lambda x: x.split(': ')[-1], labels.split('; ')))
        aka += acronyms.split('; ')
        aka = [x for x in aka if x != '']
        aka = ', '.join(aka) if aka else '-'
        address = ', '.join([city, region, country])
        parents_of_candidate = ', '.join(list(map(lambda x: conv(x), parents_of_candidate)))
        options.append((ror_id, f"{name} (also known as: {aka}) (address: {address}) (parent organizations: {parents_of_candidate})"))
    return options

def select_institute(sam):
    prompt = """
Question: "Given organization name consists of a single organization and its hierarchy.
From the list below, select the number that corresponds to the option that is exactly equal to the given organization.
If the exact lower level of organization does not exist in the given options, then find the higher level organization in the hierarchy that matches.
Do not make any guess about their roles based on their name, just compare only organization names.
Always strictly double-check all informations including full hierarchy, synonyms, parental organization, additional context, and address to avoid confusions for homonym organization.
Give only number and do not explain."
Given organization name to you: "%s"
Additional context: "%s"
Options: "%s"
Answer:
"""
    options = make_options(mysam, sam2candidates_1, sam2candidates_2, sam2candidates_3, id2items, topN=20)
    
    institute_name = sam2items[mysam][6]
    additional_context = "Most recent affiliation:" + ', '.join(sam2items[mysam][7:])
    options_text = '\n'.join([f"{i+1}. {o[1]}" for i, o in enumerate(options)])
    
    prompt_sample = prompt % (institute_name, additional_context, options_text)
    response = model.invoke(prompt_sample).content
    ror_id, final_choice = options[int(response)-1]
    log_data = {
        "id": sam,
        "prompt": prompt_sample,
        "gpt-4o-mini-response": response,
        "final-choice": final_choice,
        "ror-id": ror_id
    }
    return (ror_id, log_data)

### Run
buff_out = ''
buff_json = {}
total_records = len(sam2items)

max_retry=10
retry_delay=5

for i, mysam in enumerate(sam2items):
    if i % n_batch != n_idx: continue
    for attempt in range(max_retry):
        try:
            final_choice, log_data = select_institute(mysam)
            buff_out += f"{mysam}\t{final_choice}\n"
            buff_json[mysam] = log_data
            break
        except Exception as e:
            time.sleep(retry_delay)

# last batch
with open(f"{path_out_full}_{n_idx}_{n_batch}", 'w') as f:
    json.dump(buff_json, f, indent=4)
with open(f"{path_out}_{n_idx}_{n_batch}", 'w') as f:
    f.writelines(buff_out)
