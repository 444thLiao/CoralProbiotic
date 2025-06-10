from tqdm import tqdm
from glob import glob
import os
from os.path import *
import pandas as pd
from Bio import SeqIO
from bin.multiple_sbatch import sbatch_all
from collections import defaultdict
# bin.multiple_sbatch may require you run the following command in the terminal first
# export PYTHOHPATH=/home-user/thliao/script/evol_tk
os.chdir('/mnt/ivy/thliao/project/coral_ruegeria/nanopore_processing')


######################## parameter
error_ids = []  # will pass them
basedir = '/home-user/thliao/project/coral_ruegeria/nanopore_processing/'
indir = '/mnt/ivy/thliao/sequencing_data/nanopore_20230215'

metadata = f'{indir}/metadata'

######################## functions
def get_fq(fq):
    cmd = f"seqtk fqchk  -q0 {fq} |head -n3"
    return int(str(os.popen(cmd).read()).split('\n')[2].split('\t')[1])


pilon_jar = '/home-user/thliao/software/pilon/pilon-1.24.jar'
def run_pilon(raxon_fa,r1,r2,_odir,i=1):
    cmd = f"""bwa index {raxon_fa}
    bwa mem -t 8 {raxon_fa} {r1} {r2} > {_odir}/pilon{i}.sam
    samtools faidx {raxon_fa}
    samtools view -@ 8 -b -S {_odir}/pilon{i}.sam | samtools sort -@ 8 -Obam > {_odir}/pilon{i}.bam
    samtools index {_odir}/pilon{i}.bam
    java -Xmx16G -jar {pilon_jar} --genome {raxon_fa} --frags {_odir}/pilon{i}.bam --output {_odir}/pilon{i} --threads 8"""
    return cmd

######


TIME = indir.strip('/').split('_')[-1]
#       barcode  strain NGS genome size
# 0   barcode01   M2-49         B1 123456
# 1   barcode02   M2-57       AA12 123456
# 2   barcode03   M2-64        AB1 123456
mdf = pd.read_csv(metadata,sep='\t')
for _,row in mdf.iterrows():
    barcode = row['barcode']
    sid = row['strain']
    _odir = f"{basedir}/rawdata/{sid}/"
    if not exists(_odir):os.makedirs(_odir)
    if not exists(f"{_odir}/{sid}.fastq") or sid in error_ids:
        fastqgz = glob(f'{indir}/{barcode}/*.gz')
        if not fastqgz:continue
        all_files = ' '.join(fastqgz)
        cmd = f"zcat {all_files} > '{_odir}/{sid}_{TIME}.fastq' "
        os.system(cmd)   


gid2raw_inputs = {}
for _,row in tqdm(mdf.iterrows()):
    strain = row['strain']
    genome = row['NGS genome']
    size_raw = sum([get_fq(fq)
                    for fq in glob(f'{basedir}/rawdata/{strain}/*.fastq')                   
                    ])
    gid2raw_inputs[genome] = round(size_raw/1024/1024,2)
    
gid2NGSraw_inputs = {}
for sid in tqdm(gid2raw_inputs.keys()):
    if sid in gid2NGSraw_inputs:continue
    r1 = f'/mnt/ivy/thliao/project/coral_ruegeria/data_processing/clean_NGS/{sid}_R1.clean.fq.gz'
    r2 = f'/mnt/ivy/thliao/project/coral_ruegeria/data_processing/clean_NGS/{sid}_R2.clean.fq.gz'
    size_raw = sum([get_fq(fq) for fq in [r1 ,r2]])
    gid2NGSraw_inputs[sid] = round(size_raw/1024/1024,2)
    
    

for _,row in mdf.iterrows():
    faa = realpath(f"/mnt/ivy/thliao/project/coral_ruegeria/data_processing/all_faa/{row['NGS genome']}.faa")
    fna = faa.replace('.faa','.fna')
    size = sum([len(_.seq) for _ in SeqIO.parse(fna,'fasta')])
    mdf.loc[_,'size'] = size

ofiles = []
cmds = []
odir = f"{basedir}/canu_o"
for sid,row in mdf.iterrows():
    sid = row['NGS genome']
    strain = row['strain']
    each_odir = f"{odir}/{sid}/"
    if not exists(each_odir):os.makedirs(each_odir)
    in_fastq = f"{odir}/../rawdata/{strain}/{strain}_{TIME}.fastq"
    s = round(row['size']/10**6,1)
    if not exists(in_fastq):
        print(strain,in_fastq)
    cmd = f"canu -p {sid} -d {each_odir} genomeSize={s}m -nanopore '{realpath(in_fastq)}' "
    if not exists(f'{each_odir}/{sid}.contigs.fasta'):
        cmds.append(cmd)
        ofiles.append(f'{each_odir}/{sid}.contigs.fasta')
for _ in cmds:
    os.system(_)

import time
while 1:
    if all([exists(of) for of in ofiles]):
        break
    time.sleep(30)


ofiles = []
cmds = []
odir = f"{basedir}/canu_o"
for _,row in mdf.iterrows():
    sid = row['NGS genome']
    size = int(row['size']/10**6)
    r1 = f'/mnt/ivy/thliao/project/coral_ruegeria/data_processing/clean_NGS/{sid}_R1.clean.fq.gz'
    r2 = f'/mnt/ivy/thliao/project/coral_ruegeria/data_processing/clean_NGS/{sid}_R2.clean.fq.gz'
    if not exists(r1):
        print(sid,r1)
        continue
    input_fa = f"{odir}/{sid}/{sid}.contigs.fasta"
    each_odir = f"{odir}/{sid}/08_pilon"
    if not exists(each_odir):os.makedirs(each_odir)
    if not exists(input_fa) or exists(f"{each_odir}/pilon3.fasta"):
        continue
    cmd1 = run_pilon(input_fa,r1,r2,_odir=each_odir,i=1)
    cmd2 = run_pilon(f"{each_odir}/pilon1.fasta",r1,r2,_odir=each_odir,i=2)
    cmd3 = run_pilon(f"{each_odir}/pilon2.fasta",r1,r2,_odir=each_odir,i=3)
    cmds.append(';'.join([cmd1,cmd2,cmd3]))
    ofiles.append(f"{each_odir}/pilon3.fasta")
sbatch_all(cmds,thread_per_tasks=8,prefix_name='pilon')


import time
while 1:
    if all([exists(of) for of in ofiles]):
        break
    time.sleep(30)
    
for fna in ofiles:
    sid = fna.split('/')[-3]
    ofna = fna.replace('pilon3.fasta',f"{sid}.fna")
    if exists(ofna):
        continue
    info = f'{basedir}/canu_o/{sid}/{sid}.contigs.layout.tigInfo'
    _df = pd.read_csv(info,sep='\t',index_col=0)
    _df = _df.loc[_df['tigClass']=='contig']
    i2des = {}
    for i,row in _df.iterrows():
        des = f"{row['coverage']};{row['tigLen']};{row['sugCirc']}"
        i2des[i]=des
    seqs = [_ for _ in SeqIO.parse(fna,'fasta')]
    for i,s in enumerate(seqs):
        i = i+1
        ori_name = int(s.id.split('_')[0].replace('tig',''))
        des = i2des[ori_name]
        s.id = f"{sid}_{i}"
        s.name =  ''
        s.description = des
    with open(ofna,'w') as f1:
        SeqIO.write(seqs,f1,format='fasta')




    
    