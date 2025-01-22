import sys,os

path_root = '../../'
path_workd = path_root + '/analysis/4_check_misidentified/'
path_in = path_workd + '1697001520.compare_result.final.txt'
path_out = path_workd + '1697001520.compare_result.final.representative.txt'
path_ref = path_workd + 'cello2ccle.table'

def read_ref(path):
    child2parent = {}
    for line in open(path):
        items = line.rstrip('\n').split('\t')
        child2parent[items[0]] = items[2]
    return child2parent

if __name__=='__main__':
    child2parent = read_ref(path_ref)

    lines_out = []
    for line in open(path_in):
        sam, call, claim, cur1, cur2 = line.rstrip('\n').split('\t')
        if call not in ['-','MYLA']:
            call = child2parent[call]
        if call == "__Parent_cell_line_of_DLD-1/HCT 8/HCT 15/HRT-18":
            call = 'DLD-1'
        line_out = f"{sam}\t{call}\t{claim}\t{cur1}\t{cur2}\n"
        lines_out.append(line_out)

    with open(path_out, 'w') as f:
        f.writelines(lines_out)
