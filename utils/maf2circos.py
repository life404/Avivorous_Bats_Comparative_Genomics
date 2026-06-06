#! /usr/bin/env python3

import subprocess
from pathlib import Path
import re
from utils import cmd_run_multiple

def sort(maf, ou_path):
    maf = Path(maf)
    maf_sort = f"conda run -n maf-sort -n 2 {maf} > {maf.stem}.sort.maf"
    cmd_run_multiple(maf_sort, ou_path, check = 'maf-sort')


def convert(maf, ou_path):
    maf = Path(maf)
    maf_convert = f"conda run -n maf-convert gff {maf.stem}.sort.maf > {maf.stem}.sort.gff"
    cmd_run_multiple(maf_convert, ou_path, check = 'maf-convert')
    return(f'{maf.stem}.sort.gff')

def parse_gff(gff, ou_path):
    gff = Path(gff)
    link_file = open(ou_path / f'{gff.stem.stem}.link', 'w')
    link_lst = []
    with open(gff, 'r') as file:
        for line in file:
            if not line.startswith('#'):
                line = line.strip()
                array = line.split('\t')
                target = f'{array[0]}\t{array[3]}\t{array[4]}'
                array = array[8].split(';')[0].replace('Target=', '').split(' ')
                query = f'{array[0]}\t{array[1]}\t{array[2]}'
                link_lst.append(f"{target}\t{query}\n")
    link_file.writelines(link_lst)
    
def  