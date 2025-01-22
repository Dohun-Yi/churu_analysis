import sys, os
import re

root = '../../'
workd = root + '/analysis/1_check_output/'
path_in = root + '/reference/cellosaurus/cellosaurus.txt'
path_out = workd + '/cellosaurus.table'
path_out_reform = workd + '/cellosaurus.table.reform'
path_out_reform_acc = workd + '/cellosaurus.table.reform.acc'
path_out_reform_acc_cvcl = workd + '/cellosaurus.table.reform.acc.cvcl'

def stratify(name):
    if name == 'T-1':
        return 'TDASH1'
    elif name == 'H-4':
        return 'HDASH4'
    else:
        return re.sub(r'[^a-zA-Z0-9\s]', '', name).replace(' ', '').upper()

def read_cello(path):
    with open(path) as f:
        lines = f.readlines()
    cell_id = ''
    synD = {} # synonyms
    hiD = {} # hierarchy
    oiD = {} # originate from same individual
    oxD = {} # species of origin
    prD = {} # problematic
    acD = {} # access numbers - external
    cvD = {} # access numbers - CVCL
    csD = {} # access numbers - CVCL, Alternative
    for line in lines:
        line = line.rstrip('\n')
        if line.startswith('//'):
            cell_id = ''
            continue
        elif line.startswith('ID'):
            cell_id = line.split('   ')[-1]
            hiD[cell_id] = []
            prD[cell_id] = False
        elif cell_id == '':
            continue
        elif line.startswith('SY'):
            k, v = line.split('   ')
            synD[cell_id] = v
        elif line.startswith('HI'):
            k, v = line.split('   ')
            parent = v.split(' ! ')[-1]
            hiD.setdefault(cell_id,[]).append(parent)
        elif line.startswith('OI'):
            k, v = line.split('   ')
            neighbor = v.split(' ! ')[-1]
            oiD.setdefault(cell_id,[]).append(neighbor)
        elif line.startswith('OX'):
            k, v = line.split('   ')
            species = v.split(' ! ')[-1]
            oxD.setdefault(cell_id,[]).append(species)
        elif line.startswith('CC'):
            if 'problematic' in line.lower():
                prD[cell_id] = True
        elif line.startswith('DR'):
            k, v = line.split('   ')
            source, access = v.split('; ')
            text_access = '%s:%s' % (source, access)
            acD.setdefault(cell_id,[]).append(text_access)
        elif line.startswith('AC'):
            k, v = line.split('   ')
            cvcl = v.split(' ! ')[-1]
            cvD[cell_id] = cvcl
        elif line.startswith('AS'):
            k, v = line.split('   ')
            cvcl_as = v.split(' ! ')[-1]
            csD.setdefault(cell_id,[]).append(cvcl_as)
    return synD, hiD, oiD, oxD, prD, acD, cvD, csD

def pick_one(oiD):
    oiD_uniq = {}
    for k, v in oiD.items():
        k_uniq = sorted(sorted([k] + v), key=lambda x: len(x))[0]
        oiD_uniq[k] = k_uniq
        for kk in v:
            oiD_uniq[kk] = k_uniq
    return oiD_uniq

if __name__=='__main__':
    synD, hiD, oiD, oxD, prD, acD, cvD, csD = read_cello(path_in)
    oiD_uniq = pick_one(oiD)
    lines = []
    lines_r = []
    lines_ac = []
    lines_ac_cvcl = []
    for c in hiD.keys():
        h = hiD[c]
        cv = cvD[c]
        o = oxD[c] # species
        s = '; '.join(set(map(stratify, synD.get(c, '').split('; '))))
        ac = ';'.join(acD.get(c,[]))
        cs = ';'.join(csD.get(c,[]))
        if s == '':
            s = '-'
        if len(h) == 1:
            u = oiD_uniq.get(h[0], h[0])
        elif len(h) == 0:
            h = '-'
            u = oiD_uniq.get(c, c)
            #  o = '-'
        else:
            u = '; '.join([oiD_uniq.get(hh, hh) for hh in h])
            h = '; '.join(h)
        if len(o) > 1 or ('Human' not in o[0]):
            continue
        items = map(str,[stratify(c), stratify(u), prD[c], s, c, h[0], u])
        items_r = map(str,[c, stratify(c), u, stratify(u), prD[c], s, c, h[0], u])
        items_ac = map(str,[c, stratify(c), u, stratify(u), prD[c], s, c, h[0], u, ac])
        items_ac_cvcl = map(str,[c, stratify(c), u, stratify(u), prD[c], s, c, h[0], u, ac, cv, cs])
        lines.append('\t'.join(items) + '\n')
        lines_r.append('\t'.join(items_r) + '\n')
        lines_ac.append('\t'.join(items_ac) + '\n')
        lines_ac_cvcl.append('\t'.join(items_ac_cvcl) + '\n')
    #  with open(path_out, 'w') as f:
    #      f.writelines(lines)
    #  with open(path_out_reform, 'w') as f:
    #      f.writelines(lines_r)
    #  with open(path_out_reform_acc, 'w') as f:
    #      f.writelines(lines_ac)
    with open(path_out_reform_acc_cvcl, 'w') as f:
        f.writelines(lines_ac_cvcl)
