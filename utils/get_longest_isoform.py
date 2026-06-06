#! /usr/bin/python3

import pandas as pd
import numpy as np
import re
from Bio import SeqIO
import argparse

parse = argparse.ArgumentParser()
parse.add_argument('-g', '--gff', dest = 'gff', type = str, help = 'The gff file')
parse.add_argument('-p', '--protein', dest='protein', type = str, help = 'The protein fasta file, the ID of fasta should corresponded to the ID of feature CDS in gff file')
parse.add_argument('-c', '--cds', dest = 'cds', type = str, help = 'The CDS sequence fasta')
parse.add_argument('-o', '--output', dest = 'output', type = str, help = 'The output directory, there will two output files, one is longest protein file named "*.faa", and another is corresponding cds file named "*.fna"')
parse.add_argument('-s', '--species', dest = 'species', type = str, help = 'The species suffix, which will be texted beheind the fasta id, for example, \n>A1BG|NycAvi\nMTIEDLPDFPLEGNSLIGRYSFLFSDTPVTFSISAAPMPSDCEFSF')
args = parse.parse_args()

gene_lst = {'NCBI_ID':[],'gene':[], 'Parent':[]}
mrna_lst = {'Parent':[], 'mRNA':[]}
cds_lst = {'Parent':[], 'CDS':[]}
gff_content = open(args.gff, 'r')
for line in gff_content:
    line = line.strip()
    if not line.startswith('#'):
        array = line.split('\t')
        feature = array[2]
        attributes = array[8]
        if feature == 'gene' and ('gene_biotype=protein_coding' in attributes):
             #parent = re.search(r'ID=gene-[a-zA-z0-9-]*', attributes).group().replace('ID=', '')
             parent = attributes.split(';')[0].replace('ID=', '')
             gene = re.search(r'gene=[a-zA-Z0-9-]*', attributes).group().replace('gene=','')
             ncbi_id = re.search(r'Dbxref=GeneID:[0-9]*', attributes).group().replace('Dbxref=GeneID:', '')
             gene_lst['gene'].append(gene)
             gene_lst['Parent'].append(parent)
             gene_lst['NCBI_ID'].append(ncbi_id)
        if feature == 'mRNA':
            #mrna = re.search(r'ID=rna-[a-zA-Z0-9_]*.[0-9]*', attributes).group().replace('ID=', '')
            mrna = attributes.split(';')[0].replace('ID=', '')
            #parent = re.search(r'Parent=gene-[a-zA-Z0-9-]*', attributes).group().replace('Parent=', '')
            parent = attributes.split(';')[1].replace('Parent=', '')
            mrna_lst['mRNA'].append(mrna)
            mrna_lst['Parent'].append(parent)
        if feature == 'CDS':
            #cds = re.search(r'ID=cds-[a-zA-Z0-9_]*.[0-9]*', attributes).group().replace('ID=', '')
            cds = attributes.split(';')[0].replace('ID=', '')
            #parent = re.search(r'Parent=[a-zA-Z0-9_-]*.[0-9]*', attributes).group().replace('Parent=', '')
            parent = attributes.split(';')[1].replace('Parent=', '')
            cds_lst['CDS'].append(cds)
            cds_lst['Parent'].append(parent)
gff_content.close()

records = SeqIO.parse(args.protein, 'fasta')
seq_lst = {'CDS':[], 'SEQ':[], 'length':[]}
for record in records:
    seq_lst['CDS'].append(record.id)
    seq_lst['SEQ'].append(str(record.seq))
    seq_lst['length'].append(len(str(record.seq)))

rnaseq_records = SeqIO.parse(args.cds, 'fasta')
rnaseq_lst = {'mRNA':[], 'SEQ':[]}
for record in rnaseq_records:
    rnaseq_lst['mRNA'].append(record.id)
    rnaseq_lst['SEQ'].append(str(record.seq))
    
gene_tab = pd.DataFrame(gene_lst)
mrna_tab = pd.DataFrame(mrna_lst)
cds_tab = pd.DataFrame(cds_lst) 
seq_tab = pd.DataFrame(seq_lst)
rnaseq_tab = pd.DataFrame(rnaseq_lst)

final_tab = pd.merge(gene_tab, mrna_tab, how = 'left', on = 'Parent')
final_tab = pd.merge(final_tab, cds_tab, how='left', left_on='mRNA', right_on='Parent')
final_tab = final_tab.drop(columns=['Parent_x', 'Parent_y']).drop_duplicates().reset_index()
final_tab.mRNA = final_tab.mRNA.str.replace('rna-', '')
final_tab.CDS = final_tab.CDS.str.replace('cds-', '')
final_tab = pd.merge(final_tab, seq_tab, on = 'CDS', how='left')
final_tab = pd.merge(final_tab, rnaseq_tab, on = 'mRNA', how = 'left', suffixes=['_protein', '_mrna'])
maxindex = final_tab.groupby(by = 'gene').length.idxmax()
final_tab = final_tab.loc[maxindex]
print(final_tab)

outfile_protein = open(f'{args.output}/{args.species}.longest.faa', 'w')
outfile_cds = open(f'{args.output}/{args.species}.longest.fna', 'w')
for index, data in final_tab.iterrows():
    outfile_protein.write(f'>{data.gene}|{args.species}\n{data.SEQ_protein}\n')
    outfile_cds.write(f'>{data.gene}|{args.species}\n{data.SEQ_mrna}\n')
outfile_protein.close()
outfile_cds.close()

        