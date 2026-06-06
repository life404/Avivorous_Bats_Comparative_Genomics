#! /usr/bin/python3

from pathlib import Path
import numpy as np
import pandas as pd
import argparse
from Bio import SeqIO

def make_parse():
    parse = argparse.ArgumentParser()
    parse.add_argument('-i', '--input', help='The gene count file of orthofinder')
    parse.add_argument('-o', '--output', help='The output directory')
    args = parse.parse_args()
    return args

def concatenated(path, output):
    input_path = Path(path)
    output_path = Path(output)
    input_files = list(input_path.glob('*.pal2nal'))
    seq_dir = {}
    for input_file in input_files:
        records = SeqIO.parse(input_file, 'fasta')
        for record in records:
            species = record.id.split('|')[1]
            if species in seq_dir.keys():
                seq_dir[species].append(str(record.seq))
            else:
                seq_dir[species] = []
                seq_dir[species].append(str(record.seq))
    
    output_file = output_path / 'concatenated.fa'
    output_open = open(output_file, 'w')
    for key, value in seq_dir.items():
        output_open.write(f">{key}\n{''.join(value)}\n")
    output_open.close()
    return output_file

def ultramatic_tree(output):
    root = Path(output)
    
    output = root / 'MCMC'
    output.mkdir(exist_ok=True, parents=True)
    
    input = output / 'input'
    input.mkdir(exist_ok=True, parents=True)
    
    # using MCMC tree 
    cds_path = '/home/panda2bat/Avivorous_bat/output/12_evolution-species_tree_orthofinder/pal2nal'
    concatenated_fa = concatenated(cds_path, input)
    third_file = open((input / '3rd_sites.fa'), 'w')
    for record in SeqIO.parse(concatenated_fa, 'fasta'):
        third_file.write(f">{record.id}\n{str(record.seq)[2::3]}\n")
    third_file.close()
    

def main():
    args = make_parse()
    output = args.output
    
    ultramatic_tree(output)

if __name__ == '__main__':
    main()
    