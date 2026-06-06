#! /usr/bin/python3

import pandas as pd
import numpy as np
from Bio import SeqIO
import argparse

parse = argparse.ArgumentParser()
parse.add_argument('-g', '--gff', dest = 'gff', type = str, help = 'The gff file')
parse.add_argument('-p', '--protein', dest='protein', type = str, help = 'The protein fasta file, the ID of fasta should corresponded to the ID of feature CDS in gff file')
parse.add_argument('-c', '--cds', dest = 'cds', type = str, help = 'The CDS/cDNA sequence fasta')
parse.add_argument('-o', '--output', dest = 'output', type = str, help = 'The output directory, there will two output files, one is longest protein file named "*.faa", and another is corresponding cds file named "*.fna"')
parse.add_argument('-s', '--species', dest = 'species', type = str, default = 'NycAvi',help = 'The species suffix, which will be texted beheind the fasta id, for example, \n>A1BG|NycAvi\nMTIEDLPDFPLEGNSLIGRYSFLFSDTPVTFSISAAPMPSDCEFSF')
args = parse.parse_args()

gene_lst = {'ID':[],'Parent':[]}
mrna_lst = {'Parent':[], 'mRNA':[]}
gff_content = open(args.gff, 'r')
for line in gff_content:
    line = line.strip()
    if (not line.startswith('#')) and len(line) > 0:
        array = line.split('\t')
        feature = array[2]
        attributes = array[8]
        if feature == 'gene':
             parent = attributes.split(';')[0].replace('ID=', '')
             geneid = attributes.split(';')[0].replace('ID=', '')
             gene_lst['Parent'].append(parent)
             gene_lst['ID'].append(geneid)
        if feature == 'mRNA':
            mrna = attributes.split(';')[0].replace('ID=', '')
            parent = attributes.split(';')[1].replace('Parent=', '')
            mrna_lst['mRNA'].append(mrna)
            mrna_lst['Parent'].append(parent)
gff_content.close()

records = SeqIO.parse(args.protein, 'fasta')
seq_lst = {'mRNA':[], 'SEQ':[], 'length':[]}
for record in records:
    seq_lst['mRNA'].append(record.id)
    seq_lst['SEQ'].append(str(record.seq))
    seq_lst['length'].append(len(str(record.seq)))

rnaseq_records = SeqIO.parse(args.cds, 'fasta')
rnaseq_lst = {'mRNA':[], 'SEQ':[]}
for record in rnaseq_records:
    rnaseq_lst['mRNA'].append(record.id)
    rnaseq_lst['SEQ'].append(str(record.seq))
    
gene_tab = pd.DataFrame(gene_lst)
mrna_tab = pd.DataFrame(mrna_lst)
seq_tab = pd.DataFrame(seq_lst)
rnaseq_tab = pd.DataFrame(rnaseq_lst)

final_tab = pd.merge(gene_tab, mrna_tab, how = 'left', on = 'Parent')
final_tab = pd.merge(final_tab, seq_tab, how = 'left', on = 'mRNA')
final_tab = pd.merge(final_tab, rnaseq_tab, how = 'left', on = 'mRNA', suffixes=['_protein', '_mrna'])
final_tab = final_tab.drop(columns=['Parent']).drop_duplicates().reset_index()
maxindex = final_tab.groupby(by = 'ID').length.idxmax()
final_tab = final_tab.loc[maxindex]
print(final_tab)

outfile_protein = open(f'{args.output}/{args.species}.longest.faa', 'w')
outfile_cds = open(f'{args.output}/{args.species}.longest.fna', 'w')
for index, data in final_tab.iterrows():
    outfile_protein.write(f'>{data.ID}|{args.species}\n{data.SEQ_protein}\n')
    outfile_cds.write(f'>{data.ID}|{args.species}\n{data.SEQ_mrna}\n')
outfile_protein.close()
outfile_cds.close()

        