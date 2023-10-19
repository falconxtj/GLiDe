import pandas as pd
import regex as re
from numpy import random
import os
import numpy as np

def get_reverse_comp(seq):
    code = {'A':'T','T':'A','C':'G','G':'C','N':'N','a':'t','t':'a','c':'g','g':'c','n':'n'}
    reverse_comp_seq = ''.join([code[s] for s in seq[::-1]])
    return reverse_comp_seq

def quality_trim(seq, GCcontent_min, GCcontent_max):
    i = 0
    for s in seq:
        if s in ['G', 'C']:
            i += 1
    ratio = i / 20 * 100
    if (ratio < GCcontent_min) | (ratio > GCcontent_max):
        return False
    if 'AAAA' in seq:
        return False
    if 'TTTT' in seq:
        return False
    if 'CCCC' in seq:
        return False
    if 'GGGG' in seq:
        return False
    return True

def specificity_penalty(seq1, seq2, spacer_length, NAG=False):
    if NAG:
        penalty_vector = np.array([10]*7+[7]*5+[3]*(spacer_length-12))
    else:
        penalty_vector = np.array([8]*7+[4.5]*5+[2.5]*(spacer_length-12))
    mismatch_bool = np.array([seq1[i]!=seq2[i] for i in range(spacer_length)][::-1])
    penalty = penalty_vector.dot(mismatch_bool)
    return penalty

def ncRNA_design(genome_seq, nc_num, GCcontent_min, GCcontent_max, prefix):
    genome_seq_r = get_reverse_comp(genome_seq)

    target1 = '[ATCG]{21}GG'
    target2 = '[ATCG]{21}AG'
    NGG_temp = re.findall(re.compile(target1, flags=re.IGNORECASE), genome_seq, overlapped=True)
    NGG_nontemp = re.findall(re.compile(target1, flags=re.IGNORECASE), genome_seq_r, overlapped=True)
    NAG_temp = re.findall(re.compile(target2, flags=re.IGNORECASE), genome_seq, overlapped=True)
    NAG_nontemp = re.findall(re.compile(target2, flags=re.IGNORECASE), genome_seq_r, overlapped=True)
    NGG = NGG_temp + NGG_nontemp
    NAG = NAG_temp + NAG_nontemp
    NGG = [s[:-3] for s in NGG]
    NAG = [s[:-3] for s in NAG]

    with open('%s/NGG.fa'%prefix, 'w') as f:
        for i in range(len(NGG)):
            f.write('>%d\n%s\n' % (i, NGG[i][-12:]))
    with open('%s/NAG.fa'%prefix, 'w') as f:
        for i in range(len(NAG)):
            f.write('>%d\n%s\n' % (i, NAG[i][-12:]))

    bases = ['A', 'T', 'C', 'G']
    seqs = {}

    f = open('%s/NC.fa'%prefix, 'w')
    k = 0
    for i in range(50000):
        seq = ''.join(random.choice(bases, 20))
        if quality_trim(seq, GCcontent_min, GCcontent_max):
            seqs['NC%d' % k] = seq
            f.write('>NC%d\n%s\n' % (k, seq[-12:]))
            k += 1
    f.close()

    os.system('seqmap 2 %s/NC.fa %s/NGG.fa %s/a.txt /output_all_matches /forward_strand /available_memory:4096'%(prefix, prefix, prefix))
    os.system('seqmap 2 %s/NC.fa %s/NAG.fa %s/b.txt /output_all_matches /forward_strand /available_memory:4096'%(prefix, prefix, prefix))

    ids = []
    a = open('%s/a.txt'%prefix, 'r')
    lines = a.readlines()
    for line in lines[1:]:
        id1, id2 = int(line.split('\t')[0]), line.split('\t')[3]
        penalty = specificity_penalty(NGG[id1], seqs[id2], 20)
        if penalty < 25:
            ids.append(id2)
    a.close()
    b = open('%s/b.txt'%prefix, 'r')
    lines = b.readlines()
    for line in lines[1:]:
        id1, id2 = int(line.split('\t')[0]), line.split('\t')[3]
        penalty = specificity_penalty(NAG[id1], seqs[id2], 20)
        if penalty < 25:
            ids.append(id2)
    b.close()
    ids = list(set(ids))
    #len(ids)
    index = set(seqs.keys()) - set(ids)
    #len(index)
    NC_chosed = {i: seqs[i] for i in list(index)[:nc_num]}
    NC_chosed = pd.Series(NC_chosed)
    #NC_chosed.to_excel('NC.xlsx', header=False)
    return NC_chosed