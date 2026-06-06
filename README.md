# Avivorous Bat Genome Analysis Pipeline

**Associated manuscript:** *Comparative Genomics Reveals Convergent Evolution Between Avivorous Bats (*Ia io* and *Nyctalus aviator*)*

This repository contains the full workflow for genome assembly, annotation, and evolutionary analysis of the avivorous bats *Ia io* and *Nyctalus aviator*, supporting the comparative genomics study above. 

## Overview

This pipeline supports *Comparative Genomics Reveals Convergent Evolution Between Avivorous Bats (*Ia io* and *Nyctalus aviator*)* and spans raw read quality control through genome assembly, Hi-C–assisted chromosome scaffolding, repeat and gene annotation, synteny analysis, endogenous retrovirus (ERV) detection, and downstream evolutionary analyses including phylogenetics, positive selection, and gene loss.

## Experimental Design

**Short description (for Dryad metadata):**

> We collected tissue samples from two distantly related avivorous bats, *Ia io* and *Nyctalus aviator*, to investigate convergent genomic evolution associated with bird-eating ecology. For each species, we generated multi-platform sequencing data comprising short-read whole-genome sequencing, long-read sequencing, and Hi-C proximity ligation sequencing to assemble chromosome-level reference genomes. RNA sequencing was additionally performed to support transcript-based gene annotation. These genomic resources were used in comparative analyses with other bat species to examine synteny, gene family evolution, positive selection, and gene loss, with the goal of identifying molecular signatures of convergent adaptation to an avivorous diet.

**Data types generated:**

| Data type | Purpose |
|-----------|---------|
| Short-read WGS (paired-end) | Genome survey, hybrid assembly, and polishing |
| Long-read sequencing | Contig assembly and gap filling |
| Hi-C sequencing | Chromosome-level scaffolding |
| RNA-seq | Gene structure annotation and expression support |
| Comparative genomic datasets | Ortholog inference, phylogenomics, and convergent evolution analyses across bats |

## Directory Structure

```
.
├── README.md                          # This documentation
├── utils.py                           # Shared utilities (parallel execution, checkpointing)
├── 01.quality.py                      # Read quality control
├── 02.genome-survey.py                # Genome survey (k-mer analysis)
├── 03.genome-assembly.py              # Genome assembly
├── 04.genome-polish.py                # Assembly polishing
├── 05.genome-scaffold.py              # Hi-C chromosome scaffolding
├── 06.genome-completes.py             # Genome/annotation completeness (BUSCO)
├── 07.genome-repeat.py                # Repeat annotation
├── 08.genome-annotation.py            # Gene structure annotation
├── 09.genome-systeny.py               # Genome synteny analysis
├── 10.genome-ERVdetect.py             # Endogenous retrovirus detection
├── 11.evolution-single_copy_gene.py   # Single-copy ortholog identification
├── 12.evolution-species_tree.py       # Species tree inference
├── 13.genome-family_analysis.py       # Gene family expansion/contraction
├── 14_evolution_batch_codeml.py       # Batch CodeML positive selection
├── 15_evolution_toga_geneloss.py      # TOGA gene loss analysis
├── toga.def                           # TOGA Singularity container definition
├── panther.def                        # PANTHER Singularity container definition
└── utils/                             # Helper scripts and Jupyter notebooks
    ├── assign_id_to_gff3.py
    ├── get_longest_isoform.py
    ├── get_longest_isoform_for_pasa.py
    ├── orthofinder_single_protein2cds.py
    ├── orthoMam_divided_by_species.py
    ├── pal2nal.check.py
    ├── maf2circos.py
    └── *.ipynb                        # Result parsing and visualization notebooks
```

## Analysis Workflow

| Step | Script | Main Tools | Description |
|------|--------|------------|-------------|
| 01 | `01.quality.py` | fastp | Filter low-quality reads (discard sequences with >20% bases Q≤20 or >10% N) |
| 02 | `02.genome-survey.py` | Jellyfish | k-mer frequency profiling and genome size estimation |
| 03 | `03.genome-assembly.py` | nextDenovo, wtdbg2 | Hybrid assembly from long and short reads |
| 04 | `04.genome-polish.py` | nextPolish | Polish assembly with long and short reads |
| 05 | `05.genome-scaffold.py` | chromap, YaHS, Juicer, 3D-DNA | Hi-C–assisted chromosome-level scaffolding |
| 06 | `06.genome-completes.py` | BUSCO | Assess genome and annotation completeness |
| 07 | `07.genome-repeat.py` | EDTA | Repeat element annotation and soft masking |
| 08 | `08.genome-annotation.py` | BRAKER3, miniprot, PASA, EvidenceModeler | Gene prediction and integrated annotation |
| 09 | `09.genome-systeny.py` | LASTZ, maf2circos | Inter-species synteny alignment and visualization |
| 10 | `10.genome-ERVdetect.py` | BLAST, ERVin | Endogenous retroviral element detection |
| 11 | `11.evolution-single_copy_gene.py` | InParanoid, OrthoFinder | Single-copy ortholog identification |
| 12 | `12.evolution-species_tree.py` | PRANK, pal2nal, MCMCtree | Multi-gene concatenated phylogeny |
| 13 | `13.genome-family_analysis.py` | CAFE5 | Gene family expansion and contraction |
| 14 | `14_evolution_batch_codeml.py` | CodeML | Batch branch/site model positive selection tests |
| 15 | `15_evolution_toga_geneloss.py` | TOGA | Cross-species gene loss analysis |

## Requirements

### System

- Linux server (recommended: ≥128 cores, ≥256 GB RAM)
- Python ≥ 3.8
- Singularity / Apptainer (for containerized tools such as TOGA and InParanoid)

### Python Dependencies

```bash
pip install numpy pandas biopython tqdm rpy2
```

### Key Bioinformatics Tools

| Category | Tools |
|----------|-------|
| QC & Assembly | fastp, Jellyfish, nextDenovo, wtdbg2, nextPolish |
| Hi-C Scaffolding | chromap, YaHS, Juicer, 3D-DNA |
| Annotation | BUSCO, EDTA, BRAKER3, miniprot, PASA, EvidenceModeler |
| Evolution | InParanoid, OrthoFinder, PRANK, pal2nal, MCMCtree, CodeML, CAFE5, TOGA |
| Other | LASTZ, BLAST+, PANTHER |

> **Note:** Some scripts contain hard-coded paths (e.g. `/home/panda2bat/TOOLS/`). Update these to match your local installation, or bind-mount directories when running Singularity containers.

## Usage

Each script is invoked via command-line arguments and either prints shell commands for execution or runs the analysis directly. Typical examples:

### Example: Read Quality Control

```bash
python3 01.quality.py \
  --flist sample_list.tsv \
  --in /path/to/raw_reads/ \
  --ou /path/to/clean_reads/ \
  --na sample_prefix
```

`sample_list.tsv` format (tab-separated; column 1 = R1 filename, column 2 = R2 filename):

```
sample_R1.fq.gz    sample_R2.fq.gz
```

### Example: k-mer Survey

```bash
python3 02.genome-survey.py \
  --in /path/to/clean_reads/ \
  --ou /path/to/survey/ \
  --na bat \
  --range 21:141:10
```

### Example: Stepwise Execution (Hi-C Scaffolding)

```bash
# Step 1: Build chromap index
python3 05.genome-scaffold.py --genome genome.fasta --ou_path output/ --step 1

# Step 2: Hi-C alignment
python3 05.genome-scaffold.py --genome genome.fasta --hic_path hic_reads/ --ou_path output/ --step 2

# Steps 3–5: Scaffolding, Juicer conversion, and manual curation
```

### Checkpointing

The `cmd_run_multiple` function in `utils.py` supports checkpointing: if a `{check}.ok` file already exists in the output directory, that step is skipped automatically.

## Depositing on Dryad

### Recommended Upload Contents

| Type | Contents | Notes |
|------|----------|-------|
| **Code** | All `.py` scripts and the `utils/` directory | Ensures reproducibility |
| **Container definitions** | `toga.def`, `panther.def` | Singularity image build files |
| **Config files** | Step-specific config templates (e.g. `.nextPolish.config`) | Include if present in the repository |
| **Raw/processed data** | Genome assembly, annotation GFF3/GTF, protein/CDS sequences | Release per journal requirements |
| **Intermediate results** | OrthoFinder output, phylogenetic trees, CodeML results, etc. | Optional; useful for reviewer verification |
| **Documentation** | This README | Usage guide for data and code |

