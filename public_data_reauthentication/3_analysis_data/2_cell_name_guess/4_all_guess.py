import sys, os

path_root = '../../'
path_workd = path_root + '/analysis/2_cell_name_guess'
path_guess1 = path_workd + '/cellname_guesse.txt'
path_guess2 = path_workd + '/cellname_by_access_number.txt'
path_output = path_workd + '/cellname_guess_merge.txt'
path_output2 = path_workd + '/cellname_guess_all.txt'

def readit(path):
    sam2guess = {}
    for line in open(path):
        items = line.rstrip('\n').split('\t')
        sam = items[0]
        guess = items[1].split(';;;')
        sam2guess[sam] = set(guess) - set(['', 'none of these'])
    return sam2guess

if __name__ == '__main__':
    sam2guess1 = readit(path_guess1)
    sam2guess2 = readit(path_guess2)

    sams = sorted(sam2guess1.keys())
    lines_out = []
    lines_out2 = []
    reform = lambda x: ';;;'.join(list(x))
    for sam in sams:
        guess1 = sam2guess1[sam]
        guess2 = sam2guess2[sam]
        guess = reform(guess1 | guess2)
        line_out = f"{sam}\t{guess}\n"
        lines_out.append(line_out)
        line_out2 = f"{sam}\t{reform(guess1)}\t{reform(guess2)}\t{guess}\n"
        lines_out2.append(line_out2)

    with open(path_output, 'w') as f:
        f.writelines(lines_out)

    with open(path_output2, 'w') as f:
        f.writelines(lines_out2)

