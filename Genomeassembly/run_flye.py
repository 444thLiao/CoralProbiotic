import pandas as pd
import os
from os.path import *
from glob import glob
from Bio import SeqIO
from bin.multiple_sbatch import sbatch_all # see evol_tk repo.

meta = []
for m in glob('/home-user/thliao/data/sequencing_data/nanopore_*/metadata'):
    meta.append(pd.read_csv(m,sep='\t'))
meta = pd.concat(meta,axis=0)
strain2g = dict(zip(meta['strain'],meta['NGS genome'],))

f = '/mnt/ivy/thliao/project/coral_ruegeria/ms_figures/table SX2.xlsx'
ms_used = pd.read_excel(f,index_col=0)

used_ids = list(ms_used.index)

used_ids = ['H9','AP11','BG7','AB11']

cmds = []
odir = "/mnt/ivy/thliao/project/coral_ruegeria/nanopore_processing/assembly_flye"
for strain,sid in strain2g.items():
    if sid not in used_ids:
        continue
    _odir = f"{odir}/{sid}/"
    if not exists(_odir):
        os.makedirs(_odir)
    in_fastq = f"/mnt/ivy/thliao/project/coral_ruegeria/nanopore_processing/rawdata/{strain}/tmp.fastq"
    if not exists(in_fastq):
        fqs = glob(join(dirname(in_fastq),'*.fastq'))
        if not len(fqs)==1:
            os.system(f"cat " + ' '.join(fqs) +f' > {in_fastq}')
        else:
            in_fastq = glob(join(dirname(in_fastq),'*.fastq'))[0]
    fna = realpath(f"/mnt/ivy/thliao/project/coral_ruegeria/data_processing/all_faa/{sid}.faa").replace('.faa','.fna')
    r = list(SeqIO.parse(fna,'fasta'))
    size = round(sum([len(_.seq) for _ in r])/10**6,1)
    cmd = f"/home-user/thliao/anaconda3/bin/flye -g {size}m -o {_odir} -t 10 --nano-raw {in_fastq}"
    # v2.9
    if not exists(f'{_odir}/assembly.fasta'):
        cmds.append(cmd)
sbatch_all(cmds,thread_per_tasks=10,prefix_name='flye')

    