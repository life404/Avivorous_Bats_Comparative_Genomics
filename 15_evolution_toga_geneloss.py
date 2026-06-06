#! /usr/bin/python3

import argparse
from pathlib import Path
from utils import cmd_run_multiple, cmd_run
import concurrent.futures as cf
from Bio import SeqIO
import pandas as pd
import numpy as np


#def target_split(target, bed, selected, results):
#    root = Path(results)
#    tParts_path = root / "tParts"
#    tParts_path.mkdir(exist_ok=True, parents=True)
#
#    target_info = twoBitInfo(target, tParts_path)
#    target_fa = twoBitToFa(target, tParts_path)
#    bed_lst = {}
#    with open(bed, "r") as file:
#        for line in file:
#            line = line.strip()
#            array = line.split("\t")
#            bed_lst[array[3]] = line
#    info_lst = {
#        i.strip().split("\t")[0]: i.strip().split("\t")[1] for i in open(target_info, "r")
#    }
#
#    tParts_lst = []
#    selected_genes = [i.strip() for i in open(selected, "r")]
#    for gene in selected_genes:
#        out_dir = tParts_path / f"{gene}"
#        out_dir.mkdir(exist_ok=True, parents=True)
#        raw_bed = bed_lst[gene]
#        plank_bed_path, toga_bed_path = modified_bed(raw_bed, 5000, info_lst, out_dir)
#        target_plank_fa_path = getFasta(
#            plank_bed_path, target_fa, out_dir, f"target_plank.fa"
#        )
#        record = SeqIO.read(target_plank_fa_path, 'fasta')
#        record.id = record.id.split(':')[0]
#        SeqIO.write([record], target_plank_fa_path, 'fasta')
#        target_plank_twoBit_path = faToTwoBit(target_plank_fa_path, out_dir)
#        tParts_lst.append(target_plank_twoBit_path)
#    return tParts_lst
#
#def target_split(target, target_bed, selected, results):
#    root = Path(results)
#    
#    target_bed_df = pd.read_csv(
#        target_bed, header = None, 
#        names = [
#            'chrom',
#            'chromStart',
#            'chromEnd',
#            'name',
#            'score',
#            'strand',
#            'thickStart',
#            'thickEnd',
#            'itemRgb',
#            'blockCount',
#            'blockSizes',
#            'blockStarts'
#        ],
#        sep = '\t'
#        )
#    
#    selected_genes = [i.strip() for i in open(selected, 'w')]
#    selected_target_bed_df = target_bed_df[target_bed_df.name.isin(selected_genes)]
#    
#    for index, value in selected_target_bed_df.iterrows():
#        out_dir = root / 'togaParts' / value.
#        out_dir.mkdir(exist_ok=True, parents=True)
#

def modified_bed(raw_bed, plank, info_lst, directory):
    directory = Path(directory)
    directory.mkdir(exist_ok=True, parents=True)

    plank_left = int(plank)
    plank_right = int(plank)
    array = raw_bed.split("\t")
    chr_start = 0
    chr_end = int(info_lst[array[0]])

    if (int(array[1]) - plank_left) <= 0:
        plank_left = int(array[1]) - chr_start
    if (int(array[2]) + plank_right) >= chr_end:
        plank_right = chr_end - int(array[2])

    # The plank bed file used to extrac sequences, which will be used as target sequence in LASTZ
    plank_bed_path = directory / "plank.bed"
    tmp_bed_buff = open(plank_bed_path, "w")
    tmp_bed_buff.write(
        "\t".join([str(x) for x in 
                [
                    array[0],
                    (int(array[1]) - plank_left),
                    (int(array[2]) + plank_right),
                    array[3],
                    array[4],
                    array[5],
                ]
            ]
        )
    )
    tmp_bed_buff.close()

    # The second bed file will be modified to used in toga, the position infomation of gene in bed file will be modified based on the sequence files extracted based on the previous bed file
    toga_bed_path = directory / "toga.bed"
    toga_bed_buff = open(toga_bed_path, "w")
    toga_bed_buff.write(
        "\t".join([str(x) for x in 
                [
                    array[0],
                    plank_left,
                    (int(array[2]) - int(array[1]) + plank_left),
                    array[3],
                    array[4],
                    array[5],
                    (int(array[6]) - int(array[1]) + plank_left),
                    (int(array[7]) - int(array[1]) + plank_left),
                    array[8],
                    array[9],
                    array[10],
                    array[11],
                ]
            ]
        )
    )
    toga_bed_buff.close()

    return (plank_bed_path, toga_bed_path)

def getFasta(bed, genome, directory, out):
    out = Path(directory) / f"{out}"
    cmds = f"~/TOOLS/bedtools getfasta -fi {genome} -fo {out} -bed {bed}"
    check = cmd_run_multiple(cmds, ou_path=directory, check=".getFasta")
    return out

def twoBitInfo(twoBit, dir):
    size = dir / twoBit.name.replace(".2bit", ".info")
    cmds = f"twoBitInfo {twoBit} {size}"
    check = cmd_run_multiple(cmds, ou_path=dir, check="null")
    return size

def twoBitToFa(twoBit, dir):
    twoBit = Path(twoBit).absolute()
    fa = dir / twoBit.name.replace(".2bit", ".fa")
    check = f".twoBitToFa.{twoBit.stem}"
    cmds = f"twoBitToFa {twoBit} {fa}"
    check = cmd_run_multiple(cmds, check=check, ou_path=dir)
    return fa

def faToTwoBit(fa, dir):
    if isinstance(fa, list):
        twoBit = [dir / Path(i).name.replace(".fa", ".2bit") for i in fa]
        cmds = [f"faToTwoBit {i} {dir}/{Path(i).name.replace('.fa', '.2bit')}" for i in fa]
        check = cmd_run_multiple(cmds, ou_path=dir, check=".faToTwoBit")
        return twoBit
    else:
        twoBit = dir / Path(fa).name.replace(".fa", ".2bit")
        cmds = f"faToTwoBit {fa} {twoBit}"
        check = cmd_run_multiple(cmds, ou_path=dir, check=".faToTwoBit")
        return twoBit
    
def miniprot_index(fa, dir):
    dir = Path(dir)
    fa = Path(fa)
    miniprot_path = "/home/panda2bat/TOOLS/miniprot/miniprot-0.11_x64-linux/miniprot"
    cmds = f"{miniprot_path} -t 16 -d query.index {fa}"
    check = cmd_run_multiple(cmds, ou_path=dir, check='.miniprot_index')
    return dir / 'query.index'

def miniprot_map(fa, index, dir):
    dir = Path(dir)
    fa = Path(fa)
    miniprot_path = "/home/panda2bat/TOOLS/miniprot/miniprot-0.11_x64-linux/miniprot"
    cmds = f"{miniprot_path} --gff {index} {fa} > {dir}/query.gff"
    cmd_run(cmds, dir)
    return dir / 'query.gff'

def query_split_sub(gene, results, query_fa_records, qplank, isoform_df, index_path):
    out_path = Path(results) / 'togaParts' / gene
    out_path.mkdir(exist_ok=True, parents=True)
    
    #generate the splited isoform file
    isoform_tmp = out_path / 'isoform.tsv'
    selected_isoform_df = isoform_df[isoform_df.GeneID == gene]
    selected_isoform_df.to_csv(isoform_tmp, index=False, sep = '\t', columns = ['GeneID', 'TransID'], header = True)
    #generate the splited protein file
    pep_tmp_path = out_path / 'tmp.pep.fa'
    pep_tmp_buff = open(pep_tmp_path, 'w')
    for index, value in selected_isoform_df.iterrows():
        pep_tmp_buff.write(
            f">{value.TransID}\n{value.seq}\n"
        )
    pep_tmp_buff.close()
    ## using miniprot to identified the potential region, then modified the region, and extrac coressponding sequence
    gff_path = miniprot_map(pep_tmp_path, index_path, out_path)
    
    
    gff_info = {'chr':[], 'start':[], 'end':[]}
    with open (gff_path, 'r') as file:
        for line in file:
            if not line.startswith('#'):
                array = line.strip().split('\t')
                if array[2] == 'mRNA':
                    gff_info['chr'].append(array[0])
                    gff_info['start'].append(int(array[3]))
                    gff_info['end'].append(int(array[4]))
    
    gff_info_df = pd.DataFrame(gff_info).groupby('chr', as_index = False).agg({'start':np.min,'end': np.max})
    query_split_path = out_path / 'query.fa'
    query_split_buff = open(query_split_path, 'w')
    for index, value in gff_info_df.iterrows():
        start = value.start - qplank
        end = value.end + qplank
        chromEnd = len(str(query_fa_records[value.chr].seq))
        start = 0 if start <= 0 else start
        end = chromEnd if end >= chromEnd else end
        seq = str(query_fa_records[value.chr][start:end].seq)
        query_split_buff.write(f'>{value.chr}\n{seq}\n')
    query_split_buff.close() 
    query_split_twoBit = faToTwoBit(query_split_path, out_path)
    return out_path


def query_split(query, target_pep, selected, results, qplank, isoform):
    root = Path(results)
    
    pep_records = {}
    pep_records['id'] = [record.id for record in SeqIO.parse(target_pep, 'fasta')]
    pep_records['seq'] = [record.seq for record in SeqIO.parse(target_pep, 'fasta')]
    pep_df = pd.DataFrame(pep_records)
    isoform_df = pd.read_csv(isoform, sep = '\t', header=0)
    isoform_pep_df = pd.merge(isoform_df, pep_df, left_on='TransID', right_on='id', how='inner')
    selected_genes = [i.strip() for i in open(selected, 'r')]
    
    # generate miniprot index
    query_fa_path = twoBitToFa(query, results)
    index_path = miniprot_index(query_fa_path, results)
    query_fa_records = {record.id:record for record in SeqIO.parse(query_fa_path, 'fasta')}
    
    with cf.ProcessPoolExecutor(max_workers=64) as e:
        process_lst = [
                e.submit(query_split_sub, gene, results, query_fa_records, qplank, isoform_pep_df, index_path) for gene in selected_genes
            ]
    output_lst = list()
    for process in cf.as_completed(process_lst):
        output = process.result()
        output_lst.append(output)
    
    return output_lst

  
def target_split(target, target_bed, output_lst, tplank, results):
    root = Path(results)
    
    target_fa = twoBitToFa(target, root)
    target_records = {record.id:record for record in SeqIO.parse(target_fa, 'fasta')}
    
    target_bed_df = pd.read_csv(
        target_bed, header = None, 
        names = [
            'chrom',
            'mRNAStart',
            'mRNAEnd',
            'name',
            'score',
            'strand',
            'exonStart',
            'exonEnd',
            'itemRgb',
            'cdsCount',
            'cdsSizes',
            'cdsStarts'
        ],
        sep = '\t'
        )
    
    for out_path in output_lst:
        isoform_tmp = out_path / 'isoform.tsv'
        isoform_tmp_df = pd.read_csv(isoform_tmp, header = 0, sep = '\t')
        target_bed_df_tmp = target_bed_df[target_bed_df.name.isin(isoform_tmp_df.TransID)]
        region = target_bed_df_tmp.groupby('chrom', as_index = False).agg({'mRNAStart': np.min, 'mRNAEnd': np.max})
        chrom = region.chrom[0]
        chromStart = 0
        chromEnd = len(str(target_records[chrom].seq))
        tplank_left = int(tplank)
        tplank_right = int(tplank)
        mRNAStart = int(region.mRNAStart[0])
        mRNAEnd = int(region.mRNAEnd[0])
        if mRNAStart - tplank_left <= chromStart:
            tplank_left = mRNAStart - chromStart
        else:
            chromStart = mRNAStart - tplank_left
        if mRNAEnd + tplank_right >= chromEnd:
            tplank_right = chromEnd - mRNAEnd
        else:
            chromEnd = mRNAEnd + tplank_right
        
        target_fa = out_path / 'target.fa'
        target_fa_buff = open(target_fa, 'w')
        target_fa_buff.write(
            f">{chrom}\n{str(target_records[chrom][chromStart:chromEnd].seq)}\n"
        )
        faToTwoBit(target_fa, out_path)
        target_fa_buff.close()
        target_bed_df_tmp.loc[:, 'mRNAEnd'] = target_bed_df_tmp.loc[:, 'mRNAEnd'] - chromStart
        target_bed_df_tmp.loc[:, 'exonStart'] = target_bed_df_tmp.loc[:, 'exonStart'] - chromStart
        target_bed_df_tmp.loc[:, 'exonEnd'] = target_bed_df_tmp.loc[:, 'exonEnd'] - chromStart
        target_bed_df_tmp.loc[:, 'mRNAStart'] = target_bed_df_tmp.loc[:, 'mRNAStart'] - chromStart
        modified_bed_path = out_path / 'toga.bed'
        target_bed_df_tmp.to_csv(modified_bed_path, sep = '\t', index = False, header=False)
    
    return output_lst

def toga_run(query, t_2bit):
    run_dir = Path(t_2bit).absolute().parent
    bash_path = Path(__file__).absolute().parent / 'TOGA_run' / 'toga_run.sh'
    cmd = f"cd {run_dir} && {bash_path} -t {t_2bit} -q {query}\n"
    return cmd

def toga_run(out_path, results):
    bash_path = Path(__file__).absolute().parent / 'bin' / 'toga_run.sh'
    cmds = [
        f"cd {path} && {bash_path} -t {path / 'target.2bit'} -q {path / 'query.2bit'}\n"
        for path in out_path
    ]
    
    toga_jobs_buff = open(results / 'toga.jobs', 'w')
    toga_jobs_buff.writelines(cmds)
    toga_jobs_buff.close()
    
def make_parse():
    parse = argparse.ArgumentParser()
    parse.add_argument("--target", "-t", dest="target", help="target 2bit file")
    parse.add_argument(
        "--bed", "-b", dest="bed", help="the bed annotation of target species"
    )
    parse.add_argument(
        "--selected",
        "-s",
        dest="selected",
        help='the selected single copy genes of target file. Each row is a gene. If this parameter is "all", all genes will be used in analysis.',
    )
    parse.add_argument('--isoform', dest='isoform', help='the isoforms file. The first column is the GeneID, tht second forms is the corresponding isoform id')
    parse.add_argument("--target_plank", dest='tplank', type = int, default = 5000, help='the plank of split target')
    parse.add_argument('--target_pep', dest = 'target_pep', help = 'The pep file contains the protein sequences, this file will be used to extrac protein sequences to use in miniprot, therefore, the name of sequences should corresponding the name of transcripts in isoform file and bed file')
    parse.add_argument("--query", "-q", dest="query", help="query 2bit file")
    parse.add_argument("--query_plank", type = int, dest='qplank', default = 10 * 1000 * 1000, help='the plank of split query')
    parse.add_argument("--results", "-r", dest="results", help="results directory")
    args = parse.parse_args()
    return args

def main():
    #args = make_parse()
    #target = Path(args.target).absolute()
    #query = Path(args.query).absolute()
    #results = Path(args.results).absolute()
    #bed = Path(args.bed).absolute()
    #selected = Path(args.selected).absolute()
    #isoform = Path(args.isoform).absolute()
    #tplank = int(args.tplank)
    #qplank = int(args.qplank)
    #target_pep = Path(args.target_pep).absolute()
    
    target = Path('/home/panda2bat/Avivorous_bat/script/test/hg38.2bit')
    query = Path('/home/panda2bat/Avivorous_bat/script/test/NycAvi.2bit')
    results = Path('/home/panda2bat/Avivorous_bat/script/test/output')
    bed = Path('/home/panda2bat/Avivorous_bat/script/test/toga.transcripts.bed')
    selected = Path('/home/panda2bat/Avivorous_bat/script/test/selected_genes')
    isoform = Path('/home/panda2bat/Avivorous_bat/script/test/toga.isoforms.tsv')
    tplank = int(5000)
    qplank = int(10 * 1000 * 1000)
    target_pep = Path('/home/panda2bat/Avivorous_bat/script/test/hg38.pep.fa')
    
    
    #q_split_out_lst = query_split(query, target_pep, selected, results, qplank, isoform)
    #tmp_buff = open('/home/panda2bat/Avivorous_bat/script/test/temp_lst', 'w')
    #for i in q_split_out_lst:
    #    tmp_buff.write(f"{i}\n")
    #tmp_buff.close()
    q_split_out_lst = [Path(i.strip()) for i in open('/home/panda2bat/Avivorous_bat/script/test/temp_lst', 'r')]
    out_path = target_split(target, bed, q_split_out_lst, tplank, results)
    toga_run(out_path, results)

    #tParts_lst = target_split(target, bed, selected, results)
    #toga_parallel_cmds = [toga_run(query, t_2bit) for t_2bit in tParts_lst]
    #toga_parallel_cmds_buff = open((results) / 'toga.jobs', 'w')
    #toga_parallel_cmds_buff.writelines(toga_parallel_cmds) 
    #toga_parallel_cmds_buff.close()
   
if __name__ == "__main__":
    main()
