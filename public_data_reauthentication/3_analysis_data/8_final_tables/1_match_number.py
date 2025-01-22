import sys, os

path_root = '../../'
path_workd = path_root + '/analysis/8_final_tables/'
path_base = path_workd + '/raw/1697001520.compare_result.final.txt'

def get_sam_of_interest(path, mode='all'):
    mode = {'all':'all', 'nonskip':'nonskip'}[mode]
    sams = set()
    for line in open(path):
        items = line.rstrip('\n').split('\t')
        curation = items[4]
        if mode == 'all':
            sams.add(items[0])
        elif mode == 'nonskip' and curation != 'skip':
            sams.add(items[0])
    return sams

def filter_file(path_in, path_out, sams, column=0, header=False):
    lines_new = []
    with open(path_in) as f:
        lines = f.readlines()
    for line in lines[header:]:
        items = line.rstrip('\n').split('\t')
        if items[column] in sams:
            lines_new.append(line)
    with open(path_out, 'w') as f:
        f.writelines(lines_new)

if __name__=='__main__':
    sams1 = get_sam_of_interest(path_base, mode='all')
    sams2 = get_sam_of_interest(path_base, mode='nonskip')
    path_from = path_workd + '/raw/'
    path_to1 = path_workd + '/filt1/'
    path_to2 = path_workd + '/filt2/'

    if not os.path.isdir(path_to1):
        os.mkdir(path_to1)
    if not os.path.isdir(path_to2):
        os.mkdir(path_to2)

    params = [
            ('cellname_by_access_number.txt', 0, False),
            ('cellname_guesse.txt', 0, False),
            ('institute_selected.parsed.corrected.add_time.txt', 0, False),
            ('sam2journal.txt', 0, False),
            ('SRA_Accessions_with_expr.uniq.tab.filter_f25_dedup', 17, False)
            ]
    for f, c, h in params:
        filter_file(path_from + f, path_to1 + f, sams1, c, h)
        filter_file(path_from + f, path_to2 + f, sams2, c, h)
