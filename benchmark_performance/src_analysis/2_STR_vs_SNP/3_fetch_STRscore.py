# require python3
import requests
import json
import re


root = '../../'
workd = root + '/analysis/2_STR_vs_SNP/'
file_path = workd + '/cellosaurus.txt'
file_cell = workd + '/cell_list.txt'
path_out = workd + '/STR_score.txt'
url = "https://www.cellosaurus.org/str-search/api/batch"

def stratify(name):
    return re.sub(r'[^a-zA-Z0-9\s]', '', name).replace(' ', '').upper()

def get_cell_list(path):
    cells = set()
    for line in open(path):
        cell_name = line.rstrip('\n')
        cells.add(cell_name)
    return cells

def read_and_split_file(file_path):
    with open(file_path, 'r') as file:
        content = file.read()
        content = content.split('_____\n')[-1]
        blocks = content.split('//\n')
        return blocks

def parse_block(block):
    data = {}
    source='NULL'
    for line in block.split('\n'):
        if line.startswith('ID'):
            cell = line.replace('ID   ','')
        elif line.startswith('AC'):
            access = line.replace('AC   ','')
        elif line.startswith("ST"):
            parts = line.split(": ")
            if len(parts) == 2:
                key, value = parts[0].strip().replace("ST   ", ""), parts[1].strip()
                if key.startswith('Source'):
                    source = value.split(':')[-1].split('; ')[0]
                if '(' in value:
                    if source in value.split('(')[-1]:
                        data[key] = value.split('(')[0].rstrip(' ')
                else:
                    data[key] = value

    return cell, access, source, data

def format_data_for_api(cell, data):
    formatted_data = {
        "description": cell,
        "algorithm": 1,
        "scoringMode": 1,
        "scoreFilter": 20,
        "includeAmelogenin": True
    }
    formatted_data.update(data)
    return formatted_data

def make_api_request(data, url):
    response = requests.post(url, json=[data])  # Wrap data in a list
    return response

def parse_and_format_response(response_text):
    try:
        response_data = json.loads(response_text)
        formatted_results = []

        # Iterate through each item in the response array
        for item in response_data:
            results = item.get("results", [])
            
            # Iterate through each result in the results array
            for result in results:
                formatted_result = {
                    "accession": result.get("accession", "N/A"),
                    "name": result.get("name", "N/A"),
                    "species": result.get("species", "N/A"),
                    "bestScore": result.get("bestScore", "N/A"),
                    "problematic": result.get("problematic", "N/A")
                }
                scores = []
                for profile in result.get("profiles", []):
                    scores.append(str(profile.get("score", "N/A")))
                scores = ';'.join(scores)
                formatted_result["scores"] = scores
                formatted_results.append(formatted_result)

        return formatted_results
    except json.JSONDecodeError:
        return "Error parsing JSON"


if __name__=='__main__':
    cells = get_cell_list(file_cell)
    blocks = read_and_split_file(file_path)

    i = 0
    with open(path_out, 'w') as f:
        pass
    for block in blocks:
        lines = []
        if block.strip():  # Check if the block is not just empty spaces
            cell, access, source, parsed_data = parse_block(block)
            cell_st = stratify(cell)
            if cell_st not in cells:
                continue
            if not parsed_data:
                continue
            i += 1
            print(f'processing {cell_st} {i}')
            formatted_data = format_data_for_api(cell_st, parsed_data)
            response = make_api_request(formatted_data, url)
            formatted_response = parse_and_format_response(response.text)
            for item in formatted_response:
                #  print(f"{cell_st}\t{stratify(item['name'])}\t{item['accession']}\t{item['name']}\t{item['species']}\t{item['bestScore']}\t{item['problematic']}")
                line = f"{cell_st}\t{access}\t{stratify(item['name'])}\t{cell}\t{source}\t{item['accession']}\t{item['name']}\t{item['species']}\t{item['bestScore']}\t{item['problematic']}\t{item['scores']}"
                lines.append(line + '\n')
        with open(path_out, 'a') as f:
            f.writelines(lines)

