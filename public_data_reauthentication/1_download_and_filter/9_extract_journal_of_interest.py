import sys,os

root = '../'
workd = root + '/resources/'
input_interest = workd + '/SRA_Accessions_with_expr.uniq.tab.filter_f25_dedup'
# journal
input_bioproject = workd + '/bioproject.journals.txt'
output = workd + '/bioproject.journals.interest.txt'
# gse
input_bioproject2 = workd + '/bioproject.gse.txt'
output2 = workd + '/bioproject.gse.interest.txt'

def read_interest(path):
    proj_set = set()
    for line in open(path):
        proj = line.rstrip('\n').split('\t')[18]
        proj_set.add(proj)
    return proj_set

def write_interest(proj_set, path_in, path_out):
    lines_new = []
    for line in open(path_in):
        if line.rstrip('\n').split('\t')[0] in proj_set:
            lines_new.append(line)
    with open(path_out,'w') as f:
        f.writelines(lines_new)

if __name__=='__main__':
    proj_set = read_interest(input_interest)
    write_interest(proj_set, input_bioproject, output)
    write_interest(proj_set, input_bioproject2, output2)

