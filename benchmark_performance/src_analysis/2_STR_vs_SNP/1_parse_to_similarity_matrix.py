import sys, os

root = '../../'
workd = root + '/analysis/2_STR_vs_SNP/'
path_datalist = workd + '/datalist.txt'
path_identity = workd + '/identity/'
path_out = workd + '/matrix.txt'
path_out_lineage = workd + '/matrix_lineage.txt'

def get_datalist(path):
    sra_name = {}
    for line in open(path):
        items = line.rstrip('\n').split('\t')
        sra, name = items[0], items[1]
        sra_name[sra] = name
    return sra_name

def fetch_identity(path):
    name_posterior = {}
    name_patient = {}
    for line in open(path):
        if line.startswith('cell_id'): continue
        items = line.rstrip('\n').split('\t')
        name, patient, posterior = items[1], items[2], items[5]
        name_posterior[name] = float(posterior)
        name_patient[name] = patient # this is static
    return name_posterior, name_patient

if __name__=='__main__':
    sra_name = get_datalist(path_datalist)
    name_name_posterior = {}
    for sra, name in sra_name.items():
        print('reading %s...' % sra)
        path = path_identity + '/%s.churu.out' % sra
        name_posterior, name_patient = fetch_identity(path)
        name_name_posterior[name] = name_posterior
    
    cell_intersection = list(set(sra_name.values()) & set(name_name_posterior.values()[0].keys()))
    cell_intersection.sort(key=lambda x: name_patient[x])
    cell_concat = ['%s:%s' % (c, name_patient[c]) for c in cell_intersection]

    # all matrix line
    lines_c = ['cell\t' + '\t'.join(cell_concat) + '\n']
    for c1, cc1 in zip(cell_intersection, cell_concat):
        print ('calculating %s...' % c1)
        line = cc1
        for c2, cc2 in zip(cell_intersection, cell_concat):
            line += '\t' + str(name_name_posterior[c1][c2])
        lines_c.append(line + '\n')

    with open(path_out, 'w') as f:
        f.writelines(lines_c)

    # lineage
    patients_count = {}
    for c in cell_intersection:
        p = name_patient[c]
        patients_count.setdefault(p,0)
        patients_count[p] += 1
    patients_with_lineage = [p for p,c in patients_count.items() if c > 1]
    cell_intersection_p = [c for c in cell_intersection \
            if name_patient[c] in patients_with_lineage]
    cell_concat_p = ['%s:%s' % (c, name_patient[c]) for c in cell_intersection_p]
    lines_p = ['cell\t' + '\t'.join(cell_concat_p) + '\n']
    for c1, cc1 in zip(cell_intersection_p, cell_concat_p):
        print ('calculating %s...' % c1)
        line = cc1
        for c2, cc2 in zip(cell_intersection_p, cell_concat_p):
            line += '\t' + str(name_name_posterior[c1][c2])
        lines_p.append(line + '\n')

    with open(path_out_lineage, 'w') as f:
        f.writelines(lines_p)


