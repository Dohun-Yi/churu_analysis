import sys, os

root = '../../'
workd = root + '/analysis/1_benchmark/'
path_in = workd + 'output.txt'

def parse_line(path):
    n = 0
    tp_churu = 0
    tp_celid = 0
    tp_uniqu = 0
    fp_churu = 0
    fp_celid = 0
    fp_uniqu = 0
    ptp_churu = 0
    ptp_celid = 0
    ptp_uniqu = 0
    pfp_churu = 0
    pfp_celid = 0
    pfp_uniqu = 0
    for line in open(path):
        items = line.rstrip('\n').split('\t')
        if items[-3:] != ['O','O','O']:
            continue
        n += 1
        tp_cell = items[1]
        call_churu = filter(lambda x: x != '-', items[2].split('|'))
        call_celid = filter(lambda x: x != '-', items[3].split('|'))
        call_uniqu = filter(lambda x: x != '-', items[4].split('|'))

        if tp_cell in call_churu:
            tp_churu += 1
            fp_churu += len(call_churu) - 1
        else:
            fp_churu += len(call_churu)

        if tp_cell in call_celid:
            tp_celid += 1
            fp_celid += len(call_celid) - 1
        else:
            fp_celid += len(call_celid)

        if tp_cell in call_uniqu:
            tp_uniqu += 1
            fp_uniqu += len(call_uniqu) - 1
        else:
            fp_uniqu += len(call_uniqu)

        tp_cell = items[5]
        call_churu = set(filter(lambda x: x != '-', items[6].split('|')))
        call_celid = set(filter(lambda x: x != '-', items[7].split('|')))
        call_uniqu = set(filter(lambda x: x != '-', items[8].split('|')))

        if tp_cell in call_churu:
            ptp_churu += 1
            pfp_churu += len(call_churu) - 1
        else:
            pfp_churu += len(call_churu)

        if tp_cell in call_celid:
            ptp_celid += 1
            pfp_celid += len(call_celid) - 1
        else:
            pfp_celid += len(call_celid)

        if tp_cell in call_uniqu:
            ptp_uniqu += 1
            pfp_uniqu += len(call_uniqu) - 1
        else:
            pfp_uniqu += len(call_uniqu)

    lst = [tp_churu, tp_celid, tp_uniqu, \
            fp_churu, fp_celid, fp_uniqu]
    lst += [ptp_churu, ptp_celid, ptp_uniqu, \
            pfp_churu, pfp_celid, pfp_uniqu]
    lst += [n]
    return lst

if __name__=='__main__':
    lst = parse_line(path_in)
    print 'tp, fp, fn'
    print '\tchuru:',lst[0],lst[3], lst[12]-lst[0]
    print '\tcelid:',lst[1],lst[4], lst[12]-lst[1]
    print '\tuniqu:',lst[2],lst[5], lst[12]-lst[2]

    print 'tp, fp, fn (patient level)'
    print '\tchuru:',lst[6],lst[9], lst[12]-lst[6]
    print '\tcelid:',lst[7],lst[10], lst[12]-lst[7]
    print '\tuniqu:',lst[8],lst[11], lst[12]-lst[8]

    print 'n:', lst[12]
