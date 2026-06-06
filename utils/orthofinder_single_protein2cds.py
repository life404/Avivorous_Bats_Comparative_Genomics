#! /usr/bin/python3

from pathlib import Path
from Bio import SeqIO
import argparse


def make_parse():
    parse = argparse.ArgumentParser()
    parse.add_argument(
        "-p",
        "--protein",
        help="The directory contains single copy files",
        default="/home/panda2bat/Avivorous_bat/output/11_evolution-single_copy_gene/orthofinder/output/Results_Jun19/Single_Copy_Orthologue_Sequences",
    )
    parse.add_argument(
        "-c",
        "--cds",
        help="The directory contains all cds files",
        default="/home/panda2bat/Avivorous_bat/output/11_evolution-single_copy_gene/inparanoid/cds_filtered",
    )
    parse.add_argument(
        "-o",
        "--output",
        help="The output directory",
        default="/home/panda2bat/Avivorous_bat/output/11_evolution-single_copy_gene/orthofinder/output/Results_Jun19/Single_Copy_Orthologue_CDS",
    )
    args = parse.parse_args()
    return args

def obtain_cds(protein_path, cds_path, output):
    root = Path(output)
    root.mkdir(exist_ok=True, parents = True)
    
    protein_files = list(protein_path.glob('*.faa'))
    cds_files = list(cds_path.glob('*.fna'))
    cds_records = {}
    for fa in cds_files:
        for record in SeqIO.parse(fa, 'fasta'):
            cds_records[record.id] = record
    
    for fa in protein_files:
        tmp_cds = []
        for record in SeqIO.parse(fa, 'fasta'):
            tmp_cds.append(cds_records[record.id])
        SeqIO.write(tmp_cds, (root / f"{fa.stem}.fna"), 'fasta')

def main():
    args = make_parse()
    protein_path = Path(args.protein)
    cds_path = Path(args.cds)
    output = args.output
    
    obtain_cds(protein_path, cds_path, output)
    
if __name__ == "__main__":
    main()