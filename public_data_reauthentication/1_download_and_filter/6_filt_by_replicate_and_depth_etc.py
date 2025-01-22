import sys, os

root = '../'
workd = root + '/resources'
path_in = workd + '/SRA_Accessions_with_expr.uniq.tab.filter_f25'
path_out = workd + '/SRA_Accessions_with_expr.uniq.tab.filter_f25_dedup'

# filtration criteria
#
# remove read length > 500bp or < 40bp
# remove depth > 10x or < 0.1x
# if left sample has multiple sequencing run,
#   select TRANSCRIPTOMIC data if any
#   if TRANSCRIPTOMIC, select data whose depth is close to 0.5x
#   if GENOMIC, select data whose depth is close to 5x

def filter_by_depth_and_read_length(lines):
    lines_new = []
    for line in lines:
        items = line.rstrip('\n').split('\t')
        spots = int(items[14] if items[14] != '-' else 1e100)
        bases = int(items[15] if items[15] != '-' else 1e-100)
        cov = 1.0 * bases / (3.2 * 10 ** 9)
        readlen = 1.0 * bases / spots
        if 0.1 < cov < 10 and 40 < readlen < 500:
            lines_new.append(line)
    return lines_new

def pick_replicate(lines):
    sampleD = {}
    lines_new = []
    for line in lines:
        items = line.rstrip('\n').split('\t')
        sample = items[11]
        source = items[21]
        bases = int(items[15])
        cov = 1.0 * bases / (3.2 * 10 ** 9)
        infos = (source, cov)
        sampleD.setdefault(sample, []).append((infos, line))
    for sample in sampleD:
        runs = sampleD[sample]
        if len(runs) == 1:
            lines_new.append(runs[0][1])
            continue
        # select TRANSCRIPTOMIC if any
        if any([run[0][0] == 'TRANSCRIPTOMIC' for run in runs]):
            runs = filter(lambda x: x[0][0] == 'TRANSCRIPTOMIC', runs)
            optimal_depth = 0.5 # rna
        else:
            optimal_depth = 5 # dna
        runs.sort(key=lambda x: abs(optimal_depth - float(x[0][1])))
        lines_new.append(runs[0][1])
    return lines_new

if __name__=='__main__':
    if len(sys.argv) != 1:
        print('usage: python ~.py')
        exit()

    with open(path_in) as f:
        lines = f.readlines()

    lines_filt1 = filter_by_depth_and_read_length(lines)
    lines_filt2 = pick_replicate(lines_filt1)

    print ('raw: %i' % len(lines))
    print ('filt1: %i' % len(lines_filt1))
    print ('filt2: %i' % len(lines_filt2))
    with open(path_out, 'w') as f:
        f.writelines(lines_filt2)
