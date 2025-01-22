import sys, os

path_root = '../../'
path_workd = path_root + '/analysis/8_final_tables/'
path_srr = path_workd + '/filt1/SRA_Accessions_with_expr.uniq.tab.filter_f25_dedup'
path_journal = path_workd + '/filt1/sam2journal.txt'
path_institute = path_workd + '/filt1/institute_selected.parsed.corrected.add_time.txt'
path_guess = path_workd + '/raw/cellname_guess_all.txt'
path_curation = path_workd + '/raw/1697001520.compare_result.final.txt'
path_out = path_workd + '/final_table.txt'

def parse_journal(path):
    sam2journal = {}
    for line in open(path):
        items = line.rstrip('\n').split('\t')
        sam = items[0]
        pmids = items[2]
        journals = items[3]
        journal1 = items[5]
        if3 = items[6]
        sam2journal[sam] = [pmids, journals, journal1, if3]
    return sam2journal

def parse_institute(path):
    sam2institute = {}
    for line in open(path):
        items = line.rstrip('\n').split('\t')
        sam = items[0]
        owner = items[1]
        ror_id = items[2]
        institute_name = items[3]
        country = items[4]
        submission = items[6]
        sam2institute[sam] = [owner, ror_id, institute_name, country, submission]
    return sam2institute

def parse_guess(path):
    sam2guess = {}
    for line in open(path):
        items = line.rstrip('\n').split('\t')
        sam = items[0]
        sam2guess[sam] = items[1:]
    return sam2guess

def parse_curation(path):
    sam2curation = {}
    for line in open(path):
        items = line.rstrip('\n').split('\t')
        sam = items[0]
        sam2curation[sam] = items[1:]
    return sam2curation

if __name__=='__main__':
    sam2journal = parse_journal(path_journal)
    sam2institute = parse_institute(path_institute)
    sam2guess = parse_guess(path_guess)
    sam2curation = parse_curation(path_curation)

    header = 'biosample\tsrr\tbioproject\tguess1\tguess2\tguess\tcall\tclaim\tcuration1\tcuration2\towner\trorid\tinstitute\tcountry\tdate\tpmids\tjournals\tjournal1\tif3\n'
    lines_out = [header]
    for line in open(path_srr):
        items = line.rstrip('\n').split('\t')
        srr = items[0]
        sam = items[17]
        bioproject = items[18]
        journal_info = sam2journal.get(sam,['-']*4)
        institute_info = sam2institute.get(sam,['-']*5)
        guess_info = sam2guess[sam]
        curation_info = sam2curation[sam]
        items_new = [sam, srr, bioproject]
        items_new += guess_info + curation_info + institute_info + journal_info
        lines_out.append('\t'.join(items_new) + '\n')

    with open(path_out, 'w') as f:
        f.writelines(lines_out)
