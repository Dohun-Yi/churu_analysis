import sys, os
import re

root = '../../'
workd = root + '/analysis/1_benchmark/'
path_datalist = root + '/datalist.txt'
path_model = workd + '/Model.tsv'
path_cell_list_celid = workd + '/cell_list_celid.txt'
path_cell_list_uniqu = workd + '/cell_list_uniqu.txt'
path_cell_list_churu = workd + '/cell_list_churu.txt'
path_out = workd + '/output.txt'

def stratify(name):
    return re.sub(r'[^a-zA-Z0-9\s]', '', name).replace(' ', '').upper()


def parse_cell_list(path, db='celid'):
    cells = set()
    with open(path) as f:
        lines = f.readlines()
    for dirty in lines:
        dirty = dirty.rstrip('\n')
        if db == 'celid':
            cells.add(stratify('.'.join(dirty.split('.')[1:-1])))
        elif dirty in modelD:
            cells.add(modelD[dirty])
    return cells


def parse_model(path):
    patientD = {}
    modelD = {}
    for line in open(path):
        if line.startswith("ModelID"):
            continue
        items = line.rstrip('\n').split('\t')
        model = items[0]
        patient_id = items[1]
        stripped_name = items[3]
        patientD[stripped_name] = patient_id
        modelD[model] = stripped_name
    return patientD, modelD

def parse_datalist(path):
    cellD = {}
    for line in open(path):
        sra, cell, _, _s = line.rstrip('\n').split('\t')
        cellD[sra] = cell
    return cellD

def parse_uniqu(path):
    try:
        with open(path) as f:
            lines = f.readlines(10)
        if len(lines) <= 1:
            return '-'
        cells_matched = [l.split('\t')[0] for l in lines[1:] if l.split('\t')[5] == 'TRUE']
        if len(cells_matched) < 1:
            return '-'
        cells_matched = '|'.join(map(stratify,cells_matched))
        return cells_matched
    except Exception as e:
        return '-'

def parse_celid(path):
    try:
        with open(path) as f:
            lines = f.readlines(10)
        if len(lines) <= 1:
            return '-'
        cell_dirty = lines[1].split('with')[0].split('is')[1]
        cell_dirty = cell_dirty.split('.',1)[-1][::-1].split('.',1)[1][::-1]
        cell = stratify(cell_dirty)
        return cell
    except Exception as e:
        return '-'

def parse_churu(path):
    try:
        with open(path) as f:
            lines = f.readlines(10)
        if len(lines) <= 1:
            return '-'
        cells_matched = [l.split('\t')[1] for l in lines[1:] if float(l.split('\t')[5]) >= 0.5]
        if len(cells_matched) < 1:
            return '-'
        cells_matched = '|'.join(map(stratify,cells_matched))
        return cells_matched
    except Exception as e:
        return '-'

if __name__=='__main__':
    cellD = parse_datalist(path_datalist)
    patientD, modelD = parse_model(path_model)
    cells_celid = parse_cell_list(path_cell_list_celid)
    cells_uniqu = parse_cell_list(path_cell_list_uniqu, db=modelD)
    cells_churu = parse_cell_list(path_cell_list_churu, db=modelD)
    n_match = 0 # only uniqu
    n_mismatch = 0 # only uniqu
    lines = []
    parse_patient = lambda x: patientD.get(x,'-')
    for sra in cellD:
        path_churu = workd + '/identity/' + sra + '.churu.out'
        path_celid = workd + '/identity/' + sra + '.celid.out'
        path_uniqu = workd + '/identity/' + sra + '.uniqu.out'

        call_churu = parse_churu(path_churu)
        call_celid = parse_celid(path_celid)
        call_uniqu = parse_uniqu(path_uniqu)

        p_true = patientD.get(cellD[sra], '-')
        p_churu = '|'.join([ parse_patient(c) for c in call_churu.split('|')])
        p_celid = '|'.join([ parse_patient(c) for c in call_celid.split('|')])
        p_uniqu = '|'.join([ parse_patient(c) for c in call_uniqu.split('|')])

        output = [sra]
        output += [cellD[sra], call_churu, call_celid, call_uniqu]
        output += [ p_true, p_churu, p_celid, p_uniqu]
        output += ['O' if cellD[sra] in cells_churu else 'X']
        output += ['O' if cellD[sra] in cells_celid else 'X']
        output += ['O' if cellD[sra] in cells_uniqu else 'X']
        print('\t'.join(output))
        lines.append('\t'.join(output) + '\n')

    with open(path_out ,'w') as f:
        f.writelines(lines)
