import sys, os

root = '../../'
workd = root + '/analysis/1_check_output/'
path_in = workd + '/SRA_churu_output_compare.1697001520.txt'
path_out = workd + '/SRA_churu_output_compare.1697001520.categorize.txt'

if __name__=='__main__':
    lines_new = []
    for line in open(path_in):
        items = line.rstrip('\n').split('\t')
        log = items[-1]
        if items[3] == 'true':
            lines_new.append(line)
            continue
        if 'Kill' in log or 'FATAL' in log:
            items[3] = 'sysfault'
        elif log == '' or 'Segmentation' in log:
            items[3] = 'malform'
        elif 'samtools sort' in log:
            items[3] = 'unknown'
        else:
            items[3] = 'novar'
        lines_new.append('\t'.join(items) + '\n')

    with open(path_out, 'w') as f:
        f.writelines(lines_new)

