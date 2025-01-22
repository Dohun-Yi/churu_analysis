import sys,os
from commands import getoutput
import glob
# warning: this code is bad

batch = '1697001520' # real run

root = '../../'
workd = root + '/analysis/1_check_output/'
dir_churu = workd + '/' + batch + '/output/'
metadata = workd + 'SRA_Accessions_with_expr.uniq.tab'
path_out1 = workd + 'SRA_Accessions_with_expr.uniq.depth_length.txt'
path_out2 = workd + 'SRA_churu_output_compare.' + batch + '.txt'

def parse_meta(id_list, path):
    metaD = {}
    for line in open(path):
        items = line.rstrip('\n').split('\t')
        sra = items[0]
        spots = int(items[14]) if items[14] != '-' else -1
        bases = int(items[15]) if items[15] != '-' else -1
        readlen = -1
        depth = -1
        if spots != -1 and bases != -1:
            readlen = 1.0 * bases / spots
        if bases != -1:
            depth = 1.0 * bases / (3.2 * 10**9)
        metaD[sra] = [spots, bases, '%.1f' % readlen, '%.4f' % depth]
    return metaD

def write_meta(metaD, path_out):
    lines = []
    for key in metaD:
        line = '\t'.join(map(str,[key] + metaD[key])) + '\n'
        lines.append(line)
    with open(path_out, 'w') as f:
        f.writelines(lines)

def scan_log(dir_churu):
    logD = {}
    dir_samples = glob.glob(dir_churu + '/*/')
    for sample in dir_samples:
        sra = sample.rstrip('/').split('/')[-1]
        log_err = glob.glob(sample + '/*err')
        log_log = glob.glob(sample + '/*log')
        log_vcf = glob.glob(sample + '/*vcf')
        log_time = glob.glob(sample + '/*time')
        log_out= glob.glob(sample + '/churu.out')
        if len(log_time) != 1:
            continue # the upload is not complete

        time = getoutput('cut -d" " -f8 %s' % log_time[0]) if len(log_time) == 1 else "??"
        #  time = getoutput('cut -d" " -f9 %s' % log_time[0]) if len(log_time) == 1 else "??"
        time = time.replace('s','').split()
        nucl = getoutput('head -n 1 %s' % log_log[0]) if len(log_log) == 1 else "??"
        dying_message = getoutput('tail -n 1 %s' % log_err[0]) if len(log_err) == 1 else "??"
        n_variants = getoutput('grep -v "^#" -c %s' % log_vcf[0]) if len(log_vcf) == 1 else 0
        n_topmatch = getoutput('head -n 2 %s | tail -n 1 | cut -f 4' % log_out[0]) if len(log_out) == 1 else '-'
        if not n_topmatch: n_topmatch = '-'
        if n_topmatch != '-': success = 'true'
        else: success = 'false'
        logD[sra] = [success, time[0], time[1], time[2], nucl, n_topmatch, n_variants, dying_message]
    return logD

if __name__=='__main__':
    if len(sys.argv) != 1:
        print('usage: python ~.py')
        exit()

    metaD = parse_meta(1, metadata)
    logD = scan_log(dir_churu)

    write_meta(metaD, path_out1)

    lines = []
    for sra in logD:
        log = map(str, logD[sra])
        meta = map(str,metaD[sra])
        line = '\t'.join([sra] + meta[2:] + log) + '\n'
        lines.append(line)
    with open(path_out2, 'w') as  f:
        f.writelines(lines)
        


