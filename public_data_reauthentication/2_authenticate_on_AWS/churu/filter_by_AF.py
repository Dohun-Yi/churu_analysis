#!/usr/bin/env python
import sys, os

comp={'A':'T','T':'A','G':'C','C':'G'}

def revcomp(kmer):
    kmer_new = ''
    for nt in kmer[::-1]:
        kmer_new += comp[nt]
    return kmer_new

def get_canonical_kmer(kmer):
    return sorted([kmer, revcomp(kmer)]) # canonical, non-canonical


if __name__=='__main__':
    if len(sys.argv) != 5:
        print('usage: python ~.py input.kmer ref.kmer out.kmer out.bed')
        exit()
    path_input = sys.argv[1]
    path_ref = sys.argv[2]
    out_kmer = sys.argv[3]
    out_bed = sys.argv[4]

    with open(path_input) as f:
        counts = map(lambda x: x.rstrip('\n').split('\t'), f.readlines())
        kmerD = {k:int(n) for k, n in counts}

    kmers_new = set()
    gpos_set = set()
    for line in open(path_ref):
        chrom, pos, strand, HGVS, k1, k2 = line.rstrip('\n').split('\t')
        #  if chrom == 'chrM':
        #      continue # I HATE mitochondrial genes. they are the bottleneck
        k1c, k1cn = get_canonical_kmer(k1) # snp kmer
        gpos = '%s\t%s' % (chrom, pos)
        if k1c in kmerD:
            k2c, k2cn = get_canonical_kmer(k2) # ref kmer
            n1, n2 = kmerD[k1c], kmerD.get(k2c, 0)
            varfrac = 1.0 * n1 / (n1 + n2)
            if chrom == 'chrM':
                logline = line.rstrip('\n') + '\tSKIP_chrM\t'
            elif varfrac >= 0.01:
                kmers_new.add(k1c) # canonical kmer
                kmers_new.add(k2c) # canonical kmer
                kmers_new.add(k1cn) # non-canonical kmer for cookie
                kmers_new.add(k2cn) # non-canonical kmer for cookie
                gpos_set.add(gpos)
                # log
                logline = line.rstrip('\n') + '\t%s\t%s\t%i\t%i\t%.4f' % (k1c, k2c, n1, n2, varfrac)
            else:
                logline = line.rstrip('\n') + '\tSKIP_VARFRAC\t' + str(varfrac)
            print(logline)

    with open(out_kmer, 'w') as f:
        for kmer in kmers_new:
            f.write('%s\t1\n' % kmer)
    with open(out_bed, 'w') as f:
        gpos_sort = list(gpos_set)
        gpos_sort.sort(key=lambda x: int(x.split('\t')[1]))
        gpos_sort.sort(key=lambda x: x.split('\t')[0])
        for gpos in gpos_sort:
            f.write('%s\n' % gpos)

