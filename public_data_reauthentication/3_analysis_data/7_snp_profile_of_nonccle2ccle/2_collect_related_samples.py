import sys, os

path_root = '../../'
path_workd = path_root + '/analysis/7_snp_profile_of_nonccle2ccle/'
path_match = path_workd + '/1697001520.compare_result.final.representative.txt'
path_nonccle = path_workd + '/nonccle_records'
path_out = path_workd + '/related_records'

def read_it(path):
    table = []
    for line in open(path):
        items = line.rstrip('\n').split('\t')
        table.append(items)
    return table

if __name__=='__main__':
    table_all = read_it(path_match)
    table_nonccle = read_it(path_nonccle)

    calls = set()
    for items in table_nonccle:
        calls.add(items[1])

    lines_out = []
    for items in table_all:
        if items[4] != 'match': continue
        if items[1] in calls or items[2] in calls:
            lines_out.append('\t'.join(items)+'\n')

    with open(path_out, 'w') as f:
        f.writelines(lines_out)


    
