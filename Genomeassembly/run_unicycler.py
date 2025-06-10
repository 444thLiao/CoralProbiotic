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
# ! 3. assembly the nanopore sequencing using unicycler
cmds = []
odir = "/mnt/ivy/thliao/project/coral_ruegeria/nanopore_processing/unicycler_o"
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
    cmd = f"/home-user/thliao/anaconda3/bin/unicycler -l {in_fastq} -o {_odir}"
    if not exists(f'{_odir}/assembly.gfa'):
        cmds.append(cmd)
sbatch_all(cmds,thread_per_tasks=8,prefix_name='unicycler')

cmds = []
odir = "/mnt/ivy/thliao/project/coral_ruegeria/nanopore_processing/unicycler_o/mixed_mode"
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
    r1 = f'/mnt/ivy/thliao/project/coral_ruegeria/data_processing/clean_NGS/{sid}_R1.clean.fq.gz'
    r2 = f'/mnt/ivy/thliao/project/coral_ruegeria/data_processing/clean_NGS/{sid}_R2.clean.fq.gz'
    if not exists(r1):
        print(sid,r1)
    cmd = f"/home-user/thliao/anaconda3/bin/unicycler -1 {r1} -2 {r2} -l {in_fastq} -o {_odir}"
    if not exists(f'{_odir}/assembly.gfa'):
        cmds.append(cmd)
sbatch_all(cmds,thread_per_tasks=8,prefix_name='mixed')



# if not exists clean NGS
d1 = pd.read_csv('/mnt/ivy/thliao/project/coral_ruegeria/data_processing/20220325/input_data.tab',sep='\t',index_col=0)
d2 = pd.read_csv('/mnt/ivy/thliao/project/coral_ruegeria/data_processing/20211224/input_data.tab',sep='\t',index_col=0)
d = pd.concat([d1,d2],axis=0)

d = d.loc[['AP11', 'BG7', 'AB11'],:]
cmds = []
for _,row in d.iterrows():
    if not exists(row['R1']):
        print(_)
    cmd = f"conda run -n wgs fastp  -i {row['R1']} -I {row['R2']} -o /mnt/ivy/thliao/project/coral_ruegeria/data_processing/clean_NGS/{_}_R1.clean.fq.gz -O /mnt/ivy/thliao/project/coral_ruegeria/data_processing/clean_NGS/{_}_R2.clean.fq.gz -j /mnt/ivy/thliao/project/coral_ruegeria/data_processing/clean_NGS/{_}.json -h /mnt/ivy/thliao/project/coral_ruegeria/data_processing/clean_NGS/{_}.html -w 10"
    cmds.append(cmd)
sbatch_all(cmds,thread_per_tasks=10,prefix_name='fastp')
