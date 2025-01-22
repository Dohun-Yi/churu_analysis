import sys,os

root = '../'
workd = root + '/resources/'

path_finalset = workd + '/SRA_Accessions_with_expr.uniq.tab.filter_f25_dedup'
path_biosample1 = workd + '/biosample_set.clean.summary.interest.txt'
path_biosample1_filt = workd + '/biosample_set.clean.summary.interest.filt.txt'
path_biosample2 = workd + '/biosample_set.clean.titles.txt'
path_biosample2_filt = workd + '/biosample_set.clean.titles.filt.txt'

def get_final(path):
    samples = set()
    for line in open(path):
        items = line.rstrip('\n').split('\t')
        samples.add(items[17])
    return samples


def filt_by_final(final, path_in, path_out):
    lines_new = []
    for line in open(path_in):
        samn = line.split('\t')[0]
        if samn in final:
            lines_new.append(line)
    with open(path_out, 'w') as f:
        f.writelines(lines_new)


if __name__=='__main__':
    final = get_final(path_finalset)
    filt_by_final(final, path_biosample1, path_biosample1_filt)
    filt_by_final(final, path_biosample2, path_biosample2_filt)
