

# Overview

This repository contains the scripts used to analyze data and create figures for the manuscript "Evolutionary Genomics Guides Scalable Coral Probiotics for Climate Resilience"

This repository contains scripts and pipelines for genomic data processing, assembly, annotation, pseudogene analysis, amplicon analysis, and large-scale comparative genomics using PopCOGenT.
The codebase appears to be organized by workflow modules, each corresponding to a major step in the analysis.


## Installation

For some third-party scripts, please run the following code in your terminal to add `PYTHONPATH`

`git clone https://github.com/444thLiao/evol_tk target_dir/`

`export PYTHONPATH=$target_dir/evol_tk;$PYTHONPATH`

For the other python library,

`pip install pandas==1.3.5 biopython==1.79 tqdm==4.65.0 ete3==3.1.2 plotly==4.14.3 scipy==1.7.3`

The installtion time should be less than one hour if network is fine.

Installtion of `flye canu unicycler pseudofinder Wtp iqtree FastTreeMP` and all individual anntoation software, please refer to corresponding websites.


## Metadata

Metdata folder contains multiple files that used in the scripts and not sequencing or genomic data.



# Directory and Script Descriptions

## amplicon_analysis

`genetree_method_bin.py`

This script performs **gene tree-based population (MC) assignment** for amplicon sequence variants (ASVs) using reference genome markers and phylogenetic placement.

**Main Functions:**

- Uses **reference gene alignments** and phylogenetic trees (e.g., `parC`, `ATP5B`, `nirS2`) to guide ASV placement.
- Runs **Papara** and **EPA-ng** to align ASVs to reference sequences and place them on the reference tree.
- Generates **subset trees** with both reference genomes and ASVs using **IQ-TREE** or **FastTree**.
- Identifies **monophyletic clades** of ASVs associated with a specific main cluster (MC) based on reference genome assignments.
- Produces an `assigned_asv.txt` file mapping each ASV to its MC.

## annotations
    
`complicated.py`
This script integrates multiple genome annotation and analysis tasks into a single workflow.  
Key functions include:
- **Genomic island detection** via IslandViewer API (submit, monitor status, download results).
- **Duplicated region detection** in genomes.
- **MacSyFinder** runs for protein function modules (CasFinder, CONJScan, TFFscan, TXSScan).
- **Phage detection** using `What_the_Phage` and prophage mapping via `Prophage_Tracer`.
- Automatic calls to `single_anno.py` for running annotation modules (e.g., interproscan, COG, pseudofinder, antiSMASH, RGI, ISEScan, KOFAMSCAN).

`single_anno.py`
This is a **general-purpose genome annotation driver** supporting ~15 different bioinformatics tools.  
Given genome files (`.faa`, `.gbk`, `.ffn`, `.fna`), it can:
- Run functional annotation tools such as **KOFAMSCAN**, **InterProScan**, **antiSMASH**, **COG**, **CARD**, **VFDB**.
- Identify mobile genetic elements (**ISEScan**, **digIS**), CRISPR-Cas systems, antimicrobial resistance genes (**ResFinder**, **RGI**, **BacARscan**), plasmids (**Plascad**), phages (**What_the_Phage**), alien DNA (**AlienHunter**), and pseudogenes (**pseudofinder**).
- Manage command generation, execution, and logging, with optional SLURM job submission.
It supports `--auto_imp` to auto-detect related genome files from the protein file path according to prokka outputing format.

## Genomeassembly

`run_canu.py`
Automates **genome assembly** from Nanopore sequencing data using **Canu**, followed by polishing with **Pilon** using Illumina paired-end reads.  
It handles:
- Merging raw FASTQ files from sequencing runs.
- Running Canu with genome size estimates.
- Performing three rounds of Pilon polishing.
- Renaming and formatting final FASTA output.

`run_flye.py`
Runs **Flye** assembler on selected Nanopore sequencing datasets.  
Key features:
- Selects target strains/genomes based on metadata.
- Prepares and merges raw FASTQ files if needed.
- Calculates genome size from reference FASTA.
- Submits Flye jobs for assembly.

`run_unicycler.py`
Performs genome assembly using **Unicycler**, supporting both:
- Long-read-only mode (Nanopore data).
- Hybrid mode (Nanopore + Illumina paired-end data).
Also includes a preprocessing step with **fastp** to clean Illumina reads before hybrid assembly.

## largescale_popcogeneT

`sA.autoBatchpopForNew.py`
This script automates the process of running **PopCOGenT** for newly sequenced genomes.  
It identifies new genomes not yet assigned to main clusters (MCs), groups them based on phylogenetic proximity, and runs PopCOGenT in batch mode (with Mugsy alignment) for each group.  
The results are then used to assign or create MC IDs for these new genomes.

`# sB.merge_popcogeneT.py`
This script merges PopCOGenT clustering results from multiple groups into a unified main cluster (MC) table.  
It cross-references previous MC assignments, updates MC IDs for new genomes, and writes an updated MC mapping file for downstream analyses.

## pseudogenes/classify_bin.py

`classify_bin.py` is a Python script designed to **classify and annotate pseudogenes** in bacterial genomes, with a focus on identifying their potential inactivation mechanisms.  It integrates pseudogene predictions, insertion sequence (IS) element data, and genome annotations to provide a detailed classification for each pseudogene.

**Classification categories:**
- **IS-mediated** – Inactivated by insertion sequence elements.
- **Frameshift mutation** – Inactivated due to reading frame shifts.
- **Nonsense mutation** – Inactivated due to premature stop codons.
- **Partial remains** – Only part of the original gene remains.
- **Close to sequence edge** – Located near contig ends.
- **Unclassified** – Mechanism not determined.

## pseudogenes/merged_pseudogenes_bin.py

This script automates identification, merging, and annotation of pseudogenes in microbial genome datasets, especially for coral-associated bacteria. The core logic prevents fragmented annotation by merging pseudogene pieces that are close together physically and in homology, aiming for accurate gene loss characterization.

The script supports batch processing of many genome files and annotation outputs for high-throughput comparative genomics studies.

This tool is designed for robust pseudogene analysis, enabling researchers to track gene decay and genome evolution in microbial populations.

Based on [Pseudofinder](https://github.com/filip-husnik/pseudofinder) but improvded 

## search_MC_distinguishable_genes.py

This script identifies and analyzes marker genes distinguishable across microbial populations(MCs). Identified ATP5B, parC. 


# Notes
Within the scripts, if you find any import like `from bin.format_newick import renamed_tree`. Please referred to the other repo `https://github.com/444thLiao/evol_tk`.

# Data availability
Raw sequencing data are available on NCBI under several BioProjects with the private access link for reviewers. 

The raw reads of 409 and assembled genomes of 419 Ruegeria isolates (10 genomes missing raw data) sequenced on the secondgeneration platform (DNBseq) and 34 Ruegeria isolates sequenced on the third-generation platform (Nanopore) are available under NCBI BioProject accessions [PRJNA1264799](https://dataview.ncbi.nlm.nih.gov/object/PRJNA1264799?reviewer=mp7rulknapinfb95c9733nf72r) and [PRJNA1275854](https://dataview.ncbi.nlm.nih.gov/object/PRJNA1275854?reviewer=g4mffnb8b965tq6b2662cpt6i8), respectively. 

The parC and ATP5B amplicon sequencing data for natural coral samples are available under [PRJNA1275610](https://dataview.ncbi.nlm.nih.gov/object/PRJNA1275610?reviewer=nst2vifa08hqr8otcjvv39ge4l) and [PRJNA1275576](https://dataview.ncbi.nlm.nih.gov/object/PRJNA1275576?reviewer=g865emc830spbkko9b86qt). 

The 16S rRNA and parC amplicon sequencing data for outplanted corals are 436 available under [PRJNA1275585](https://dataview.ncbi.nlm.nih.gov/object/PRJNA1275585?reviewer=9ne88342f7at8ndur1klq05mo8) and [PRJNA1275617](https://dataview.ncbi.nlm.nih.gov/object/PRJNA1275617?reviewer=e6g3aatk7pg8r6m6ao0riovjnu), respectively.

# Publication

Under review and submission.


# Contact Us
If you have any questions or suggestions on these scripts, you are welcome to contact us via email: l0404th@gmail.com.
If you have any questions on the paper and experimental parts, you are welcome to contact us via email: hluo2006@gmail.com.