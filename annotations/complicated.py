"""
Including several taks that can be put together with the other annotating tasks.

such as
islandviewer that needs submit and download

"""
import requests, sys
import os
from requests_toolbelt.multipart.encoder import MultipartEncoder

token = ''
def islandviewer(mygenome):
    #mygenome="/mnt/ivy/thliao/project/coral_ruegeria/nanopore_processing/canu_o/AL3/09_prokka/AL3.gbk"
    server = "https://www.pathogenomics.sfu.ca/islandviewer"
    ext = "/rest/submit/"
    genome2jobid = {}
    gid = mygenome.split("/")[-3]
    multipart_data = MultipartEncoder(fields={"format_type": "GENBANK",
                                              'email_addr': 'XXXX',
    #  For incomplete genomes include a reference accession
                'ref_accnum': 'NZ_CP031946.1',
                'genome_name':gid,
                'genome_file': ('filename', open(mygenome, 'rb'), 'text/plain')})
    headers={'Content-Type': multipart_data.content_type,
            'x-authtoken': token,
            # token will be expires on 2023-07-10 02:51
            # go to https://www.pathogenomics.sfu.ca/islandviewer/user/token/ to get a new token
            }
    r = requests.post(server+ext, headers=headers, data=multipart_data)
    if not r.ok:
        r.raise_for_status()
        sys.exit("failed to submit job")
    decoded = r.json()
    genome2jobid[gid] = decoded['token']
    return genome2jobid

# Check status
def check_status(jobid):
    return os.popen(f"curl https://www.pathogenomics.sfu.ca/islandviewer/rest/job/{jobid}/ -H 'X-authtoken:{token}'").read()
# Check status
def download_result(jobid,ofile,format='tab'):
    cmd = f"curl https://www.pathogenomics.sfu.ca/islandviewer/rest/job/{jobid}/download/{format}/ -H 'x-authtoken:{token}' > {ofile}"
    os.system(cmd)

genome2jobid = {row.split('\t')[0]:row.split('\t')[1].strip()
                for row in open('/mnt/ivy/thliao/project/coral_ruegeria/nanopore_processing/annotations/islandviewer/gid2jobid.tsv').read().split('\n')
                if row }
g2success = {row.split('\t')[0]:row.split('\t')[1].strip()
                for row in open('/mnt/ivy/thliao/project/coral_ruegeria/nanopore_processing/annotations/islandviewer/gid2jobid.tsv').read().split('\n')
                if row and 'SUCCESS' in row}
from glob import glob
from os.path import exists
from tqdm import tqdm
import time
gid2fna = {
    faa.split("/")[-3]: faa
    for faa in glob(f"/mnt/ivy/thliao/project/coral_ruegeria/nanopore_processing/canu_o/*/09_prokka/*.fna")}
print(len(gid2fna),len(genome2jobid),len(set(gid2fna).difference(g2success)))

process_ids = set(gid2fna).difference(set(g2success))
for genome in tqdm(process_ids):
    ofile = f'/mnt/ivy/thliao/project/coral_ruegeria/nanopore_processing/annotations/islandviewer/{genome}.tsv'
    _genome2jobid = islandviewer(mygenome=f'/mnt/ivy/thliao/project/coral_ruegeria/nanopore_processing/canu_o/{genome}/09_prokka/{genome}.gbk')
    
    jobid = _genome2jobid[genome]
    while 1:
        if 'status": "Complete' in check_status(jobid):
            break
        time.sleep(180)
    download_result(jobid, f'/mnt/ivy/thliao/project/coral_ruegeria/nanopore_processing/annotations/islandviewer/{genome}.tsv')
    download_result(jobid, f'/mnt/ivy/thliao/project/coral_ruegeria/nanopore_processing/annotations/islandviewer/{genome}.gbk','genbank')
    genome2jobid.update(_genome2jobid)

with open('/mnt/ivy/thliao/project/coral_ruegeria/nanopore_processing/annotations/islandviewer/gid2jobid.tsv','w') as f1:
    for gid,jobid in genome2jobid.items():
        ofile = f'/mnt/ivy/thliao/project/coral_ruegeria/nanopore_processing/annotations/islandviewer/{gid}.tsv'
        if exists(ofile):
            if not ('Server Error (500)' in open(ofile).read() or getsize(ofile)==0):
                status = 'SUCCESS'
            else:
                status = 'FAILED'
        f1.write(f"{gid}\t{jobid}\t{status}\n")

### duplicated region
from .calculate_duplication import get_G2dup_regions
for ml,mi in [(1000,99),(500,99),
              (2000,99),(3000,99),
              (500,98),(500,95),(500,90),(1000,98),(1000,95),(1000,90)]:
    genome2dup_ratio,genome2duplication_regions,sub_dupdf = get_G2dup_regions(ml,mi)
    sub_dupdf.to_csv(f'/mnt/ivy/thliao/project/coral_ruegeria/nanopore_processing/annotations/manual_duplication/dup_{ml}.{mi}.tsv',sep='\t',index=0)
    with open(f'/mnt/ivy/thliao/project/coral_ruegeria/nanopore_processing/annotations/manual_duplication/dup_{ml}.{mi}.list','w') as f1:
        for contig,region_list in genome2duplication_regions.items():
            for r in region_list:
                f1.write(f"{contig}\t{r[0]}\t{r[1]}\n")


###
gid2nano_fna = {
    fna.split("/")[-3]: fna
    for fna in glob(
        f"/mnt/ivy/thliao/project/coral_ruegeria/nanopore_processing/canu_o/*/09_prokka/*.fna"
    )
    if fna.split('/')[-1].replace('.fna','') == fna.split('/')[-3]
}

from Bio import SeqIO
import pandas as pd
from collections import defaultdict
basedir = '/home-user/thliao/project/coral_ruegeria/nanopore_processing/'
odir = '/home-user/thliao/project/coral_ruegeria/nanopore_processing/annotations/macsyfinder_output'
circular_contigs = defaultdict(list)
for row in open('/home-user/thliao/project/coral_ruegeria/nanopore_processing/canu_o/chromosome.txt').read().strip().split('\n'):
    circular_contigs[row.split('\t')[0]].append(row.split('\t')[1] )

cmds = []
for gid,fna in tqdm(gid2nano_fna.items()):
    records = {_.id:_ for _ in SeqIO.parse(fna,'fasta')}
    with open(f'{odir}/{gid}.t','w') as f1:
        for contig_name,seq in records.items():
            if contig_name in circular_contigs[gid]:
                v = 'circular'
            else:
                v = 'linear'
            f1.write(f"{contig_name.replace('_','x')}:{v}\n")
    faa = fna.replace('.fna','.faa')
    records = {_.id:_ for _ in SeqIO.parse(faa,'fasta')}
    gff = pd.read_csv(fna.replace('.fna','.gff'),sep='\t',index_col=0,comment='#',header=None,low_memory = False)
    l2c = {}
    for _,row in gff.iterrows():
        if _.startswith(f"{gid}_"):
            l2c[row[8].split(';')[0].split('=')[-1]]=_
    new_records = []
    for cid,seq in records.items():
        oid = seq.id
        seq.id = seq.name = seq.description = ''
        seq.id = f"{l2c[oid].replace('_','x')}_{oid.replace('_','x')} {oid}"
        new_records.append(seq)
    with open(f'{odir}/{gid}_new.faa','w') as f1:
        SeqIO.write(new_records,f1,'fasta-2line')
    for d in ['CasFinder','CONJScan','TFFscan','TXSScan']:
        cmd = f"/home-user/thliao/anaconda3/envs/macsyfinder/bin/macsyfinder --models-dir /home-user/thliao/db/protein_db/ARG_related/macsyfinder --models {d} all  --sequence-db {odir}/{gid}_new.faa --db-type gembase  --topology-file {odir}/{gid}.t -o {odir}/{gid}_{d}"
        if not exists(f"{odir}/{gid}_{d}"):
            cmds.append(cmd)




################################### For public database ones
from Bio import SeqIO
import pandas as pd
def prepare_macsyfinder(fna,faa,odir):
    cmds = []
    records = {_.id:_ for _ in SeqIO.parse(fna,'fasta')}
    with open(f'{odir}/{genome}.t','w') as f1:
        for contig_name,seq in records.items():
            v = 'circular'
            f1.write(f"{contig_name.replace('_','x')}:{v}\n")
    records = {_.id:_ for _ in SeqIO.parse(faa,'fasta')}
    gff = pd.read_csv(realpath(fna).replace('.fna','.gff'),
                      sep='\t',index_col=0,comment='#',header=None,low_memory = False)
    l2c = {}
    for _,row in gff.iterrows():
        if _.startswith(f"{genome}_"):
            l2c[row[8].split(';')[0].split('=')[-1]]=_
    new_records = []
    for cid,seq in records.items():
        oid = seq.id
        seq.id = seq.name = seq.description = ''
        seq.id = f"{l2c[oid].replace('_','x')}_{oid.replace('_','x')} {oid}"
        new_records.append(seq)
    with open(f'{odir}/{genome}_new.faa','w') as f1:
        SeqIO.write(new_records,f1,'fasta-2line')
    for d in ['CasFinder','CONJScan','TFFscan','TXSScan']:
        cmd = f"/home-user/thliao/anaconda3/envs/macsyfinder/bin/macsyfinder --models-dir /home-user/thliao/db/protein_db/ARG_related/macsyfinder --models {d} all --sequence-db {odir}/{genome}_new.faa --db-type gembase --topology-file {odir}/{genome}.t -o {odir}/{genome}_{d}"
        if not exists(f"{odir}/{genome}_{d}"):
            cmds.append(cmd)
    return cmds



################################### For phage
import os
from os.path import *
import pandas as pd
from glob import glob
from bin.multiple_sbatch import sbatch_all,batch_iter
from subprocess import check_call
import os

def is_completed(path):
    return bool(glob(join(path, '*.html')))
os.chdir('/mnt/ivy/thliao/project/ruegeria_prophage/singularity_images')
aid2fna = {fna.split('/')[-1].replace('.fna',''):fna
           for fna in glob('/mnt/ivy/thliao/project/coral_ruegeria/nanopore_processing/canu_o/*/09_prokka/*.fna')
           if '_' not in fna.split('/')[-1]}
for genome in {'GNM004006175' ,
'GNM003443535',
'GNM000011965',
'GNM000014065'}:
    fna = f'/mnt/ivy/thliao/project/coral_ruegeria/data_processing/pub_dataset/prokka_o/{genome}/{genome}.fna'
    aid2fna[genome] = fna

odir = '/mnt/ivy/thliao/project/coral_ruegeria/nanopore_processing/annotations/wtp_output'
batch_ids = batch_iter([aid
                        for aid in list(aid2fna)
                        if not is_completed(join(odir,aid))],5)

print(len(batch_ids))
# export SINGULARITY_LOCALCACHEDIR=/mnt/ivy/thliao/project/ruegeria_prophage/sin_tmp , prior to run the nextflow.
m = 'docker'
tmp_dir = "/mnt/ivy/thliao/project/coral_ruegeria/nanopore_processing/.tmp_data"
# cmds = []
for aids in batch_ids:
    if exists(f"{tmp_dir}"):
        os.system(f"rm -r {tmp_dir}")
    for _ in aids:
        if exists(f"{odir}/{_}"):
            os.system(f"rm -r {odir}/{_}")
    os.system(f'mkdir {tmp_dir}')
    for aid in aids:
        cmd = f"ln -s {aid2fna[aid]} {tmp_dir}/"
        check_call(cmd,shell=1)
    try:
        #if using singularity, you should add two more parameters
        cmd1 = f"nextflow run 444thLiao/What_the_Phage -r v1.0.6 --cores 10 --max-cores 10 -profile local,singularity -w /mnt/ivy/thliao/tmp/ --fasta '{tmp_dir}/*.fna'  --output {odir} --dv --databases /mnt/ivy/thliao/project/ruegeria_prophage/nextflow-autodownload-databases --cachedir /mnt/ivy/thliao/project/ruegeria_prophage/singularity_images "
        cmd2 = f"nextflow run 444thLiao/What_the_Phage -r v1.0.6 --cores 10 --max-cores 10 -profile local,docker -w /mnt/ivy/thliao/tmp/ --fasta '{tmp_dir}/*.fna'  --output {odir} --dv --databases /mnt/ivy/thliao/project/ruegeria_prophage/nextflow-autodownload-databases"
        if m == 'docker':
            check_call(cmd2,shell=1)
        elif m == 'singularity':
            check_call(cmd1,shell=1)
    except KeyboardInterrupt:
        break
    except:
        pass


# ! prophagetracer
""" 
bwa mem B4 /home-user/thliao/project/coral_ruegeria/data_processing/clean_NGS/B4_R1.clean.fq.gz /home-user/thliao/project/coral_ruegeria/data_processing/clean_NGS/B4_R2.clean.fq.gz > ./B4.sam

#/mnt/ivy/thliao/project/coral_ruegeria/nanopore_processing/bowtie2_map


#export INDIR=/mnt/ivy/thliao/project/coral_ruegeria/nanopore_processing/bowtie2_map
samtools view -S -b ./B4.sam -o ./B4.bam
samtools sort -Obam -o ./B4.sort.bam ./B4.bam

samtools markdup -r ./B4.sort.bam ./B4.rmdup.bam
samtools view ./B4.rmdup.bam -o ./B4.rmdup.sam
"""

from os.path import *
import os
from glob import glob
from subprocess import check_call
gid2fna = {}
for fna in glob('/mnt/ivy/thliao/project/coral_ruegeria/nanopore_processing/canu_o/*/09_prokka/*.fna'):
    gid = fna.split('/')[-1].replace('.fna','')
    if '_' not in fna.split('/')[-1]:
        gid2fna[gid] = fna
cmds = []        
_odir = '/mnt/ivy/thliao/project/coral_ruegeria/nanopore_processing/annotations/prophage_tracer/'
for gid,fna in gid2fna.items():
    odir = join(_odir,gid)
    if not exists(join(odir,f'{gid}.prophage.out')):
        r1 = f"/home-user/thliao/project/coral_ruegeria/data_processing/clean_NGS/{gid}_R1.clean.fq.gz"
        r2 = f"/home-user/thliao/project/coral_ruegeria/data_processing/clean_NGS/{gid}_R2.clean.fq.gz"
        if not exists(r1):
            print(gid)
            continue
        #os.system(f"bwa index {fna} -p {dirname(realpath(fna))}/{gid}")
        c1 = f"""mkdir -p {odir} && cd {odir} && /home-user/thliao/anaconda3/envs/wgs/bin/bwa mem -t 20 {dirname(realpath(fna))}/{gid} {r1} {r2} > {odir}/{gid}.sam && bash /home-user/thliao/software/Prophage_Tracer/prophage_tracer.sh -m {odir}/{gid}.sam -r {fna} -p {gid} -t 20"""
        cmds.append(c1)
    
    

gid2faa = {
    fna.split("/")[-3]: fna
    for fna in glob(
        f"/mnt/ivy/thliao/project/coral_ruegeria/nanopore_processing/canu_o/*/09_prokka/*.faa"
    )
    if fna.split('/')[-1].replace('.faa','') == fna.split('/')[-3]
}
for genome in {'GNM004006175' ,
'GNM003443535',
'GNM000011965',
'GNM000014065'}:
    fna = f'/mnt/ivy/thliao/project/coral_ruegeria/data_processing/pub_dataset/prokka_o/{genome}/{genome}.faa'
    gid2faa[genome] = fna
        
cmds = []
for sid,faa in gid2faa.items():
    os.system(f"python /mnt/ivy/thliao/project/coral_ruegeria/single_anno.py -p {faa} -o /mnt/ivy/thliao/project/coral_ruegeria/nanopore_processing/annotations/ -auto_imp -sl digIS #interproscan,cog,digIS,antismash,rgi,ISEScan,KOFAMSCAN")
    cmds.extend(open('/mnt/ivy/thliao/project/coral_ruegeria/nanopore_processing/annotations/cmd.log').read().strip().split('\n'))
cmds = [_ for _ in cmds if _]
cmds = []
for i in sids:
    os.system(f"python /mnt/ivy/thliao/project/coral_ruegeria/single_anno.py -p /mnt/ivy/thliao/project/coral_ruegeria/nanopore_processing/canu_o/{i}/09_prokka/{i}.faa -o /mnt/ivy/thliao/project/coral_ruegeria/nanopore_processing/annotations/ -auto_imp -sl pseudofinder")
    cmds.extend(open('/mnt/ivy/thliao/project/coral_ruegeria/nanopore_processing/annotations/cmd.log').read().strip().split('\n'))


cmds = []
for i in sids:
    os.system(f"python ../single_anno.py -p /mnt/ivy/thliao/project/coral_ruegeria/Merged_popcogeneT/inproteins/{i}.faa -o /mnt/ivy/thliao/project/coral_ruegeria/data_processing/annotations -auto_imp -sl digIS,antismash,rgi,isescan,kofamscan")
    cmds.extend(open('/mnt/ivy/thliao/project/coral_ruegeria/data_processing/annotations/cmd.log').read().strip().split('\n'))
