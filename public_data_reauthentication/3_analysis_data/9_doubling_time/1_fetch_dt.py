# require python3
import re

root = '../../'
workd = root + '/analysis/9_doubling_time/'
file_path = workd + '/cellosaurus.txt'
path_out = workd + '/doubling_time.txt'

def read_and_split_file(file_path):
    with open(file_path, 'r') as file:
        content = file.read()
        content = content.split('_____\n')[-1]
        blocks = content.split('//\n')
        return blocks

def remove_text_within_parentheses(text):
    pattern = r'\([^()]*\)'
    while re.search(pattern, text):
        text = re.sub(pattern, '', text)
    return text.rstrip('.').strip()

def parse_block(block):
    dtlist = []
    for line in block.split('\n'):
        if line.startswith('ID'):
            cell = line.replace('ID   ','')
        elif line.startswith("CC"):
            line = remove_text_within_parentheses(line.replace('CC   ',''))
            if 'Doubling time:' not in line:
                continue
            line = line.replace('Doubling time: ','').strip()
            for a in line.split(','):
                for b in a.split(';'):
                    b=b.strip().rstrip('s')
                    try:
                        hour = b.lstrip('~><=').split(' ')[0]
                        if '-' in hour:
                            hour = sum(map(float,hour.split('-')))/2
                        else:
                            hour = float(hour)
                    except Exception as e:
                        print(e, cell, line)
                        continue
                    if b.endswith('hour'):
                        hour = hour
                    elif b.endswith('day'):
                        hour = hour*24
                    elif b.endswith('week'):
                        hour = hour*24*7
                    else:
                        print ('error: time not specified', b, cell)
                        continue
                    dtlist.append(hour)
    return cell, dtlist

if __name__=='__main__':
    blocks = read_and_split_file(file_path)

    i = 0
    lines = []
    for block in blocks:
        if block.strip():  # Check if the block is not just empty spaces
            i += 1
            cell, dtlist = parse_block(block)
            dt_all = ';'.join(map(str,dtlist))
            if dtlist:
                dt_mean = str(sum(dtlist)/len(dtlist))
                line = f'{cell}\t{dt_mean}\t{dt_all}'
                lines.append(line + '\n')
    with open(path_out, 'w') as f:
        f.writelines(lines)

