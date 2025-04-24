
'''
Author: Huibao Feng & Tongjun Xiang
Description: this code can be used to design the sgRNA library for pooled CRISPRi screening for bacteria.
Date: 2025-04-24
'''
import configparser
import shutil

import numpy as np
import pandas as pd
import regex as re
import os
import sys
import datetime

import NC_design

import smtplib
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart

#import D3GB

ISOTIMEFORMAT = '%Y%m%d'

'''
Function: read out configure parameters which contain:
1. Gene annotation file (format:ptt/rnt/)
2. Genome sequence file (format:fna)
3. Coding gene file (format:ffn/frn/)
4. Offtarget threshold
5. Lower limit of GC-content
6. Upper limit of GC-content
Input: text (configure file name)
Return: var_dic (configure variables, dictionary)
'''


def configure_process(text):
    config = configparser.ConfigParser()
    config.read(text)
    section = ''.join(config.sections())
    var_dic = {option: config.get(section, option) for option in config.options(section)}
    return var_dic


'''
Function: read out the input genome sequence file
Input: genome_file (genome sequence file, fasta format)
Return: genome_seq (genome sequence, string)
'''


def get_genome_seq(prefix, genome_file, format):
    genome_seq = ''
    sequence_dic = {}
    with open('./' + genome_file, 'r') as file:
        #with open('%s/map_view/result.fna' % prefix, 'a') as fna:
        lines = file.readlines()
        name = lines[0][1:].split(' ')[0]
        for line in lines:
            if line[0] == '>':
                sequence_name = line[1:].split(' ')[0]
                sequence_dic[sequence_name] = len(genome_seq)
                #fna.write('>' + sequence_name + '\n')
            elif line[0] != '>':
                genome_seq += line.rstrip()
                #fna.write(line)
    return genome_seq, sequence_dic, name


'''
Function: read out the input gene annotation file
Input: annotation_file (gene annotation file), format , CDS (Boolean)
1. For the old version of NCBI refseq database (FTP site: ftp://ftp.ncbi.nlm.nih.gov/genomes/archive/old_refseq/Bacteria/),
   the format can be either ptt (for protein-coding genes) or rnt (for RNA-coding genes)
2. For the new version of NCBI refseq database (FTP site: ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/), the format
   is
Return: a dataframe containing the information of coding genes
'''


def fasta_process(annotation_file, format, name, CDS):
    if format in ['ptt', 'rnt']:
        ptt_rnt_info = pd.read_csv('./' + annotation_file, skiprows=2, sep='\t')
        fasta_info = pd.DataFrame(columns=['Seqid', 'Source', 'Type', 'Start', 'End', 'Score', 'Strand', 'Phase', 'Attributes'])
        fasta_info['Strand'] = ptt_rnt_info['Strand']
        fasta_info['Start'] = [item.split('..')[0] for item in ptt_rnt_info['Location']]
        fasta_info['End'] = [item.split('..')[1] for item in ptt_rnt_info['Location']]
        for index, row in ptt_rnt_info.iterrows():
            fasta_info.loc[index]['Attributes'] = 'PID=' + str(row['PID']) + ';Parent=gene-' + row[
                'Synonym'] + ';COG=' + row['COG'] + ';Name=' + row['Synonym'] + ';Gene=' + row['Gene'] + ';Product=' + \
                                                  row[
                                                      'Product']
        if CDS:
            fasta_info['Type'] = 'CDS'
        else:
            fasta_info['Type'] = 'RNA'
        fasta_info['Source'] = '.'
        fasta_info['Score'] = '.'
        fasta_info['Phase'] = '.'
        fasta_info['Seqid'] = name
    else:
        fasta_info = pd.read_csv(annotation_file, sep='\t', header=None, comment='#').dropna(axis=0)
        fasta_info.columns = ['Seqid', 'Source', 'Type', 'Start', 'End', 'Score', 'Strand', 'Phase', 'Attributes']
        if CDS:
            fasta_info = fasta_info[fasta_info['Type'] == 'CDS']
        else:
            fasta_info = fasta_info[['RNA' in fasta_info['Type'][i] for i in fasta_info.index]]
    fasta_info = fasta_info.reset_index(drop=True)
    return fasta_info


def annotation_process(annotation_file, format, name, CDS=True):
    if format in ['ptt', 'rnt']:
        ptt_rnt_info = pd.read_csv('./' + annotation_file, skiprows=2, sep='\t')
        gene_info = pd.DataFrame(columns=['Seqid', 'Source', 'Type', 'Start', 'End', 'Score', 'Strand', 'Phase', 'Attributes'])
        gene_info['Strand'] = ptt_rnt_info['Strand']
        gene_info['Start'] = [item.split('..')[0] for item in ptt_rnt_info['Location']]
        gene_info['End'] = [item.split('..')[1] for item in ptt_rnt_info['Location']]
        for index, row in ptt_rnt_info.iterrows():
            gene_info.loc[index]['Attributes'] = 'PID=' + str(row['PID']) + ';Parent=gene-' + row[
                'Synonym'] + ';COG=' + row['COG'] + ';Name=' + row['Gene'] + ';Gene=' + row['Gene'] + ';Product=' + row[
                                                     'Product']
        if CDS:
            gene_info['Type'] = 'CDS'
        else:
            gene_info['Type'] = 'RNA'
        gene_info['Source'] = '.'
        gene_info['Score'] = '.'
        gene_info['Phase'] = '.'
        gene_info['Length'] = ptt_rnt_info['Length']
        gene_info['Seqid'] = name
    else:
        gene_info = pd.read_csv(annotation_file, sep='\t', header=None, comment='#').dropna(axis=0)
        gene_info.columns = ['Seqid', 'Source', 'Type', 'Start', 'End', 'Score', 'Strand', 'Phase', 'Attributes']
        if CDS:
            gene_info = gene_info[gene_info['Type'] == 'CDS']
        else:
            gene_info = gene_info[['RNA' in gene_info['Type'][i] for i in gene_info.index]]
        gene_info['Length'] = gene_info[["Start", "End"]].apply(
            lambda x: (max(int(x["Start"]), int(x["End"])) - min(int(x["Start"]), int(x["End"])) + 1) / 3 - 1, axis=1)
    return gene_info


'''
Function: generate a dictionary to label each nucleotide in the genome
Input: gene_info (information of coding gene, dataframe), genome_seq (genome sequence, string)
Return: position_dic
1. For nucleotides positioned in coding regions (which are not overlapped with each other), it would be {postition of
   nucleotide: index of coding gene}
2. For nucleotides positioned in overlapping regions, it would be {postition of nucleotide: None}
3. Otherwise, it would be {postition of nucleotide: None}
'''


def get_position_dic(gene_info, genome_seq):
    position_dic = {i + 1: None for i in range(len(genome_seq))}
    for i in gene_info.index:
        start, end = gene_info.loc[i, 'Start'], gene_info.loc[i, 'End']
        overlap = []
        for k in range(int(start), int(end) + 1):
            if position_dic[k]:
                position_dic[k] = None
                overlap.append(k)
            elif k in overlap:
                continue
            else:
                position_dic[k] = i
    return position_dic


'''
Function: generate a fasta file which contains all coding gene candidate
Input: coding_file (coding gene file), position_dic (label of each nucleotide), format
1. For the old version of NCBI refseq database (FTP site: ftp://ftp.ncbi.nlm.nih.gov/genomes/archive/old_refseq/Bacteria/),
   the format can be either ffn (for protein-coding genes) or frn (for RNA-coding genes)
2. For the new version of NCBI refseq database (FTP site: ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/), the format
   is
Output: coding gene file, fasta format
'''


def coding_file_process(prefix, gene_info, genome_seq, position_dic):
    with open('%s/coding_file' % prefix, 'w') as outfile:
        for i in gene_info.index:
            pos = [gene_info.loc[i, 'Start'], gene_info.loc[i, 'End']]
            start, end = min(int(pos[0]), int(pos[1])), max(int(pos[0]), int(pos[1]))
            for j in range(start, end + 1):
                if position_dic[j] != None:
                    outfile.write('>' + str(position_dic[j]) + '\n')
                    seq = genome_seq[start - 1: end]
                    if gene_info.loc[i, 'Strand'] != '+':
                        seq = get_reverse_comp(seq)
                    outfile.write(seq + '\n')
                    break


'''
Function: check if two genes are in the same cluster
Input: gene_info (information of coding genes, dataframe), ID1, ID2 (index of the two genes), percIdentity (similarity of
       two genes), alnLength (length of consensus sequences)
Return: True if in a cluster, otherwise False
'''


def is_cluster(gene_info, ID1, ID2, percIdentity, alnLength):
    if ID1 >= ID2:
        return False
    elif percIdentity / 100 < 0.95:
        return False
    elif alnLength / (gene_info.loc[ID1, 'Length']) < 0.95:
        return False
    elif alnLength / (gene_info.loc[ID2, 'Length']) < 0.95:
        return False
    else:
        return True


'''
Function: align all the coding genes with each other, cluster genes with high similarity
Requirement: blast software (ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
Input: gene_info (information of coding genes, dataframe)
Return: cluster list (genes in the same cluster would be in the same list)
'''


def blastn_process(prefix, gene_info):
    os.system(
        'makeblastdb -in %s/coding_file -dbtype nucl -title %s/blast_db -out %s/blast_db' % (prefix, prefix, prefix))
    os.system('blastn -query %s/coding_file -db %s/blast_db -evalue 0.001 -out %s/blast_result -outfmt 6' % (
        prefix, prefix, prefix))
    file = open('%s/blast_result' % prefix, 'r')
    lines = file.readlines()
    cluster_info = open('%s/cluster_info' % prefix, 'w')
    cluster_lst = []
    for line in lines:
        ID1, ID2, percIdentity, alnLength = line.split('\t')[:4]
        ID1, ID2 = int(ID1), int(ID2)
        if is_cluster(gene_info, ID1, ID2, float(percIdentity), float(alnLength)):
            cluster_info.write(line)
            flag = True
            for ID_lst in cluster_lst:
                if (ID1 in ID_lst) or (ID2 in ID_lst):
                    flag = False
                    ID_lst.extend([ID1, ID2])
            if flag:
                cluster_lst.extend([[ID1, ID2]])
    cluster_lst = [list(set(cluster)) for cluster in cluster_lst]
    return cluster_lst


'''
Function: get the reverse complementary sequences
Input: sequence of forward strand
Return: sequence of reverse strand
'''


def get_reverse_comp(seq):
    code = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N', 'K': 'M', 'M': 'K', 'R': 'Y', 'Y': 'R', 'S': 'S',
            'W': 'W', 'B': 'V', 'V': 'B', 'D': 'H', 'H': 'D', 'a': 't', 't': 'a', 'c': 'g', 'g': 'c', 'n': 'n',
            'k': 'm', 'm': 'k', 'r': 'y', 'y': 'r', 's': 's', 'w': 'w', 'b': 'v', 'v': 'b', 'd': 'h', 'h': 'd'}
    reverse_comp_seq = ''.join([code[s] for s in seq[::-1]])
    return reverse_comp_seq


'''
Function: search for possible targets
Input: seq (genome sequence), PAM (protospacer adjacent motif), spacer_length (length of spacer), strand (True if forward,
       otherwise False)
Return: target_dic (a dictionary which contains all possible targets, {position: spacer})
'''


def find_target(seq, PAM, spacer_length, strand):
    PAM = PAM.replace('N', '[ATCG]').replace('n', '[ATCG]').replace('R', '[AG]').replace('r', '[AG]').replace('Y',
                                                                                                              '[TC]').replace(
        'y', '[TC]').replace('V', '[ACG]').replace('v', '[ACG]')
    target = '[ATCG]{%d}' % spacer_length + PAM
    if not strand:
        seq = get_reverse_comp(seq)
    target_dic = {}
    for i in re.finditer(re.compile(target, flags=re.IGNORECASE), seq, overlapped=True):
        if strand:
            target_dic[i.start() + 1] = i.group()[:spacer_length]
        else:
            seq_length = len(seq)
            target_dic[seq_length - i.start()] = i.group()[:spacer_length]
    return target_dic


'''
Function: get target candidates
Input: gene_info (information of coding genes, dataframe), position_dic (label of each nucleotide), temp_candidate (targets
       which positioned in template strand), nontemp_candidate (targets which positioned in non-template strand), strand (True
       if targets positioned in template strand, otherwise False)
Return: target_candidate_dic (target candidates), nontarget_candidate_dic (other targets)
'''


def get_target_candidate_dic(prefix, gene_info, position_dic, temp_candidate, nontemp_candidate, strand):
    target_candidate_dic, nontarget_candidate_dic = {}, {}
    pseudo_pos = max(position_dic) + 1
    for i in temp_candidate:
        if position_dic[i] != None:
            flag = True if gene_info.loc[position_dic[i], 'Strand'] == '+' else False
            if flag == strand:
                target_candidate_dic[i] = temp_candidate[i]
            else:
                nontarget_candidate_dic[pseudo_pos] = temp_candidate[i]
                pseudo_pos += 1
        else:
            nontarget_candidate_dic[pseudo_pos] = temp_candidate[i]
            pseudo_pos += 1
    for i in nontemp_candidate:
        if position_dic[i] != None:
            flag = True if gene_info.loc[position_dic[i], 'Strand'] == '-' else False
            if flag == strand:
                target_candidate_dic[i] = nontemp_candidate[i]
            else:
                nontarget_candidate_dic[pseudo_pos] = nontemp_candidate[i]
                pseudo_pos += 1
        else:
            nontarget_candidate_dic[pseudo_pos] = nontemp_candidate[i]
            pseudo_pos += 1
    with open('%s/target_candidate' % prefix, 'w') as target_candidate:
        for i in target_candidate_dic:
            target_candidate.write('>' + str(i) + '\n')
            target_candidate.write(target_candidate_dic[i][-12:] + '\n')
    with open('%s/nontarget_candidate' % prefix, 'w') as nontarget_candidate:
        for i in nontarget_candidate_dic:
            nontarget_candidate.write('>' + str(i) + '\n')
            nontarget_candidate.write(nontarget_candidate_dic[i][-12:] + '\n')
    return target_candidate_dic, nontarget_candidate_dic


'''
Function: calculate the penalty score to evaluate the off-target probability
Input: seq1 (sequence of spacer1), seq2 (sequence of spacer2), spacer_length(length of spacer), NAG (True if PAM == 'NAG',
       otherwise False)
Output: penalty score (intager)
'''


def specificity_penalty(seq1, seq2, spacer_length, NAG=False):
    if NAG:
        penalty_vector = np.array([10] * 7 + [7] * 5 + [3] * (spacer_length - 12))
    else:
        penalty_vector = np.array([8] * 7 + [4.5] * 5 + [2.5] * (spacer_length - 12))
    mismatch_bool = np.array([seq1[i] != seq2[i] for i in range(spacer_length)][::-1])
    penalty = penalty_vector.dot(mismatch_bool)
    return penalty


'''
Function: check if two coding gene are in the same cluster
Input: ID1 (coding gene ID1), ID2 (coding gene ID2), cluster_lst (list of all clusters)
Output: flag (True if in the same cluster, otherwise False)
'''


def in_a_cluster(ID1, ID2, cluster_lst):
    flag = False
    if ID1 == ID2:
        return True
    for cluster in cluster_lst:
        if (ID1 in cluster) & (ID2 in cluster):
            flag = True
    return flag


'''
Function: filter the potential offtarget candidates
Requirement: seqmap (http://www-personal.umich.edu/~jianghui/seqmap/#download)
Input: off_threshold (threshold of offtarget penalty), spacer_length (length of spacer), target_candidate_dic (target candidates),
       nontarget_candidate_dic (other candidates), position_dic (label of each nucleotide), cluster_lst (list of all clusters),
       target_for_NAG (True if PAM == 'NAG', otherwise False)
Output: offtarget_dic
'''


def offtarget_filter(prefix, off_threshold, spacer_length, target_candidate_dic, nontarget_candidate_dic, position_dic,
                     cluster_lst, target_for_NAG=None):
    os.system(
        'seqmap 2 %s/target_candidate %s/target_candidate %s/target_mismatch.txt /output_all_matches /forward_strand /available_memory:4096' % (
            prefix, prefix, prefix))
    os.system(
        'seqmap 2 %s/target_candidate %s/nontarget_candidate %s/nontarget_mismatch.txt /output_all_matches /forward_strand /available_memory:4096' % (
            prefix, prefix, prefix))
    offtarget_log = pd.DataFrame(index=target_candidate_dic.keys(),
                                 columns=['ID', 'sequence', 'offtarget_sequence', 'penalty', 'off_num'])
    offtarget_log['ID'] = target_candidate_dic.keys()
    offtarget_log['sequence'] = target_candidate_dic.values()
    offtarget_sequence = {i: [] for i in target_candidate_dic.keys()}
    offtarget_penalty = {i: [] for i in target_candidate_dic.keys()}
    with open('%s/target_mismatch.txt' % prefix, 'r') as file:
        lines = file.readlines()
        for line in lines[1:]:
            ID1, ID2 = int(line.split('\t')[0]), int(line.split('\t')[3])
            if not in_a_cluster(position_dic[ID1], position_dic[ID2], cluster_lst):
                penalty = specificity_penalty(target_candidate_dic[ID1], target_candidate_dic[ID2], spacer_length)
                if penalty < off_threshold:
                    offtarget_sequence[ID2].append(target_candidate_dic[ID1])
                    offtarget_penalty[ID2].append(penalty)
    with open('%s/nontarget_mismatch.txt' % prefix, 'r') as file:
        lines = file.readlines()
        for line in lines[1:]:
            ID1, ID2 = int(line.split('\t')[0]), int(line.split('\t')[3])
            penalty = specificity_penalty(target_candidate_dic[ID2], nontarget_candidate_dic[ID1], spacer_length)
            if penalty < off_threshold:
                offtarget_sequence[ID2].append(nontarget_candidate_dic[ID1])
                offtarget_penalty[ID2].append(penalty)
    if target_for_NAG:
        os.system(
            'seqmap 2 %s/target_candidate %s/NAG_file %s/NAG_mismatch.txt /output_all_matches /forward_strand /available_memory:4096' % (
                prefix, prefix, prefix))
        with open('%s/NAG_mismatch.txt' % prefix, 'r') as file:
            lines = file.readlines()
            for line in lines[1:]:
                ID1, ID2 = int(line.split('\t')[0]), int(line.split('\t')[3])
                penalty = specificity_penalty(target_candidate_dic[ID2], target_for_NAG[ID1], spacer_length, NAG=True)
                if penalty < off_threshold:
                    offtarget_sequence[ID2].append(target_for_NAG[ID1])
                    offtarget_penalty[ID2].append(penalty)
    offtarget_log['offtarget_sequence'] = offtarget_sequence.values()
    offtarget_log['penalty'] = offtarget_penalty.values()
    off_num = [len(off_penalty) for off_penalty in offtarget_penalty.values()]
    offtarget_log['off_num'] = off_num
    offtarget_log = offtarget_log[offtarget_log['off_num'] > 0]
    offtarget_log.to_excel('%s/result/offtarget_log.xlsx' % prefix, index=0)
    for i in offtarget_log['ID']:
        target_candidate_dic.pop(i)


def GC_content_filter(target_candidate_dic, GCcontent_min, GCcontent_max, spacer_length):
    unqualified_seq = []
    for i in target_candidate_dic:
        GC_content = 0
        seq = target_candidate_dic[i]
        for j in seq:
            if j in ['C', 'G']:
                GC_content += 1
        GC_content = GC_content / spacer_length * 100
        if (GC_content < GCcontent_min) | (GC_content > GCcontent_max):
            unqualified_seq.append(i)
    for i in unqualified_seq:
        target_candidate_dic.pop(i)


def main(prefix, configureFile):
    os.makedirs('%s/result' % prefix, exist_ok = True)

    # configure process
    starttime = datetime.datetime.now()
    var_dic = configure_process(configureFile)
    reference_file = var_dic['reference_file'].split(",")
    off_threshold = float(var_dic['off_threshold'])
    GCcontent_min = float(var_dic['gccontent_min'])
    GCcontent_max = float(var_dic['gccontent_max'])
    spacer_length = int(var_dic['spacer_length'])
    target_strand = True if var_dic['strand'] == 'template' else False
    CDS = True if var_dic['target'] == 'cds' else False
    NC_number = int(var_dic['nc_number'])
    genome_seq = ''
    gene_info = pd.DataFrame(columns=['Seqid', 'Source', 'Type', 'Start', 'End', 'Score', 'Strand', 'Phase', 'Attributes', 'Length'])
    fasta_info = pd.DataFrame(columns=['Seqid', 'Source', 'Type', 'Start', 'End', 'Score', 'Strand', 'Phase', 'Attributes'])
    sequence_dic = {}
    for i in range(0, int(len(reference_file) / 2)):
        genome_file = reference_file[2 * i + 1]
        genome_format = genome_file.split('.')[-1]
        genome_seq_temp, sequence_dic_temp, name = get_genome_seq(prefix, genome_file, genome_format)
        sequence_dic.update(sequence_dic_temp)
        annotation_file = reference_file[2 * i]
        if annotation_file != 'None':
            annotation_format = annotation_file.split('.')[-1]
            gene_info = pd.concat([gene_info, annotation_process(annotation_file, annotation_format, name, CDS)])
            fasta_info = pd.concat([fasta_info, fasta_process(annotation_file, annotation_format, name, CDS)])
        genome_seq += genome_seq_temp

    gene_info = gene_info.reset_index(drop=True)
    fasta_info = fasta_info.reset_index(drop=True)
    for i in gene_info.index:
        gene_info.loc[i, 'Start'] = int(gene_info.loc[i, 'Start']) + sequence_dic[gene_info.loc[i, 'Seqid']]
        gene_info.loc[i, 'End'] = int(gene_info.loc[i, 'End']) + sequence_dic[gene_info.loc[i, 'Seqid']]

    # coding process and cluster identification
    gene_info['Guide_seq'] = [[] for i in range(len(gene_info))]
    gene_info['Guide_pos'] = [[] for i in range(len(gene_info))]
    position_dic = get_position_dic(gene_info, genome_seq)
    coding_file_process(prefix, gene_info, genome_seq, position_dic)
    cluster_lst = blastn_process(prefix, gene_info)

    temp_candidate = find_target(genome_seq, 'NGG', spacer_length, True)
    nontemp_candidate = find_target(genome_seq, 'NGG', spacer_length, False)
    target_for_NAG1 = find_target(genome_seq, 'NAG', spacer_length, True)
    target_for_NAG2 = find_target(genome_seq, 'NAG', spacer_length, False)
    target_for_NAG = list(set(list(target_for_NAG1.values()) + list(target_for_NAG2.values())))
    target_for_NAG = {i: target_for_NAG[i] for i in range(len(target_for_NAG))}
    with open('%s/NAG_file' % prefix, 'w') as NAG_file:
        for i in target_for_NAG:
            NAG_file.write('>' + str(i) + '\n')
            NAG_file.write(target_for_NAG[i][-12:] + '\n')
    target_candidate_dic, nontarget_candidate_dic = get_target_candidate_dic(prefix, gene_info, position_dic, temp_candidate, nontemp_candidate, target_strand)

    offtarget_filter(prefix, off_threshold, spacer_length, target_candidate_dic, nontarget_candidate_dic, position_dic, cluster_lst, target_for_NAG)

    GC_content_filter(target_candidate_dic, GCcontent_min, GCcontent_max, spacer_length)

    sgRNA_list = []
    for i in target_candidate_dic:
        Seqid = gene_info.loc[position_dic[i], 'Seqid']
        startpoint = sequence_dic[Seqid]
        gene_info.loc[position_dic[i], 'Guide_seq'].insert(0, target_candidate_dic[i])
        gene_info.loc[position_dic[i], 'Guide_pos'].insert(0, i - startpoint)
        if gene_info.loc[position_dic[i], 'Strand'] == '+':
            end = i - startpoint
            start = end - spacer_length - 2
            strand = '+' if target_strand else '-'
        else:
            start = i - startpoint
            end = start + spacer_length + 2
            strand = '-' if target_strand else '+'
        attributes = 'Sequence=' + str(target_candidate_dic[i]) + ';Name=sgRNA-' + Seqid + '-' + str(i - startpoint)
        sgRNA_list.append([Seqid, '.', 'sgRNA', start, end, '.', strand, '.', attributes])
    sgRNA_info = pd.DataFrame(sgRNA_list, columns=['Seqid', 'Source', 'Type', 'Start', 'End', 'Score', 'Strand', 'Phase', 'Attributes'])
    for i in gene_info.index:
        gene_info.loc[i, 'guide_num'] = len(gene_info.loc[i, 'Guide_pos'])
    gene_info['Start'] = fasta_info['Start']
    gene_info['End'] = fasta_info['End']
    gene_info = gene_info.reset_index(drop=True)

    # NC RNA design
    if NC_number > 0:
        NC_list = NC_design.ncRNA_design(genome_seq, NC_number, GCcontent_min, GCcontent_max, prefix)
        gene_info.loc[len(gene_info.index)] = ['.', '.', 'NC RNA', '.', '.', '.', '.', '.', '.', '.', NC_list.tolist(),
                                               '.', NC_number]
    gene_info.to_excel('%s/result/gene_info.xlsx' % prefix, index=False)
    endtime = datetime.datetime.now()
    print((endtime - starttime).seconds)


if __name__ == '__main__':
    main('.', sys.argv[1])