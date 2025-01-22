import sys, os
# This code take some times (few minutes)

root = '../../'
workd = root + '/analysis/3_country_institute'
path_S2date = workd + '/biosample_set.clean.summary.interest.txt'
path_in = workd + '/institute_selected.parsed.corrected.txt'
path_out = workd + '/institute_selected.parsed.corrected.add_time.txt'

def makeD(path, idx1=0, idx2=1, empty='-', interest=False):
    print(f"making dictionary from {path}")
    d = {}
    with open(path) as f:
        lines = f.readlines()
    for line in lines:
        items = line.rstrip('\n').split('\t')
        k = items[idx1]
        if interest and k not in interest:
            continue
        v = items[idx2] if len(items) > idx2 else empty
        d[k] = v
    return d

if __name__=='__main__':
    interest = set(makeD(path_in, 0, 1).keys())
    d_S2date = makeD(path_S2date, 0, 2, interest=interest)

    lines_out = []
    for line in open(path_in):
        items = line.rstrip('\n').split('\t')
        sam = items[0]
        date = d_S2date[sam]
        line_out = line.rstrip('\n') + '\t' + date + '\n'
        lines_out.append(line_out)

    with open(path_out,'w') as f:
        f.writelines(lines_out)
