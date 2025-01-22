import sys, os
import math

root = '../../'
workd = root + '/analysis/9_doubling_time/'
path_dt = workd + '/doubling_time.txt'
path_pub = workd + '/1697001520.compare_result.final.txt'
path_out = workd + '/doubling_time_fold_change.txt'

def get_dt(path):
    dtD = {}
    for line in open(path):
        items = line.rstrip('\n').split('\t')
        name = items[0]
        dt = float(items[1]) # mean hours
        dtD[name] = dt
    return dtD

if __name__=="__main__":
    dtD = get_dt(path_dt)
    lines = []
    for line in open(path_pub):
        items = line.rstrip('\n').split('\t')
        sam = items[0]
        call = items[1]
        claim = items[2]
        curation = items[4]

        call_dt = dtD.get(call, 0)
        claim_dt = dtD.get(claim, 0)
        fc = 'N/A'
        if call_dt and claim_dt:
            fc = '%.4f' % math.log(call_dt / claim_dt, 2)
        line = '\t'.join(map(str,[sam, call, claim, curation, call_dt, claim_dt, fc]))
        lines.append(line + '\n')
    with open(path_out, 'w') as f:
        f.writelines(lines)
