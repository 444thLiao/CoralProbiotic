#!/home-user/thliao/anaconda3/bin/python3.7

"""
binary file for annotating genomeic sequences with ~15 softwares 

"""
import sys
sys.path.insert(0,'/home-user/thliao/script/evol_tk')
from bin.multiple_sbatch import sbatch_all

import os
from os.path import *
from glob import glob

################################# setting ####################################
KOFAMSCAN_exe = '/home-user/thliao/software/kofamscan/exec_annotation'
KOFAMSCAN_profiles = '/mnt/home-db/pub/protein_db/kegg/v20230301/profiles'
KOFAMSCAN_ko_list = '/mnt/home-db/pub/protein_db/kegg/v20230301/ko_list'
remove_intermediates = True
ISESCAN_exe = "python /home-user/thliao/software/ISEScan-1.7.2.3/isescan.py"
ISESCAN_thread = "20"
digIS_exe = "python /home-user/thliao/script/digIS/digIS_search.py"

antishmash_exe = "source activate /home-user/thliao/anaconda3/envs/antismash; antismash "
antishmash_db = "/home-user/thliao/db/antimash_db"

resfinder_exe = "/home-user/thliao/anaconda3/envs/resfinder4/bin/run_resfinder.py "
resfinder_db = "/home-user/thliao/db/protein_db/ARG_related/resfinder/db_resfinder"

rgi_exe = "/home-user/thliao/anaconda3/envs/rgi/bin/rgi"

bacarscan_db = '/mnt/maple/thliao/data/protein_db/ARG_related/BacARscan/BacARscan_standalone/pARGhmm/254_pARG.hmm'

Plascad_exe = '/home-user/thliao/anaconda3/bin/Plascad'
abricate_exe = '/home-user/thliao/anaconda3/bin/abricate'
abricate_db_list = ['virctors','megares','card','vfdb','argannot','resfinder','plasmidfinder','ncbi','megaresv3']

interpro_exe = '/home-user/thliao/software/interproscan-5.63-95.0/interproscan.sh'
interpro_application_list = "CDD,COILS,Gene3D,HAMAP,MobiDBLite,PANTHER,Pfam,PIRSF,PRINTS,SFLD,SMART,SUPERFAMILY,TIGRFAM"
interpro_cpu = '5'

CRISPRCasFinder_exe = 'source activate /home-user/thliao/anaconda3/envs/macsyfinder && /home-user/thliao/anaconda3/envs/macsyfinder/bin/perl /home-user/thliao/software/CRISPRCasFinder/CRISPRCasFinder.pl'
CRISPRCasFinder_set1 = '/home-user/thliao/software/CRISPRCasFinder/sel392v2.so'

CARD_db = '/mnt/maple/thliao/data/protein_db/ARG_related/CARD_3.2.6/protein_fasta_protein_homolog_model.fasta'

VF_db = '/mnt/maple/thliao/data/protein_db/VFDB/20230303/VFDB_setB_pro.fas'

alienhunter_exe = '/home-user/thliao/download/alien_hunter-1.7/alien_hunter'

COG_diamond_db = '/mnt/maple/thliao/data/protein_db/COG/2020/cog-20_renamed.dmnd'
##############################################################################



global executed_cmd
global LOGGER
LOGGER = []
executed_cmd = []

def check(ofile,cmd,name,dry_run=True):
    if exists(ofile):
        LOGGER.append(f'{name} existed')
    else:
        if not dry_run and not exists(dirname(ofile)):
            os.makedirs(dirname(ofile))
        LOGGER.append(f'{name} run {cmd}')
    if (not dry_run) and (not exists(ofile)):
        #print(name)
        executed_cmd.append(cmd)
    
def run_cmd(faa,fna,gbk,ffn,
            genome,odir,
            sections = [],
            dry_run=True):
    # !KOFAMSCAN
    if 'KOFAMSCAN'.lower() in sections:
        oname = join(odir,'KOFAMSCAN',genome+'.kofamout')
        cmd = f"{KOFAMSCAN_exe} -p {KOFAMSCAN_profiles} -k {KOFAMSCAN_ko_list} --tmp-dir .{genome}_tmp -o {oname} -f mapper-one-line --no-report-unannotated {faa} "
        if remove_intermediates:
            cmd += f' && rm -rf .{genome}_tmp'
        check(oname,cmd,'KOFAMSCAN',dry_run=dry_run)

    # !ISEScan
    if 'ISEscan'.lower() in sections:
        oname = join(odir,'ISEscan',genome)
        cmd = f"{ISESCAN_exe} --seqfile {fna} --output {oname} --removeShortIS --nthread {ISESCAN_thread}"
        check(oname,cmd,'ISEscan',dry_run=dry_run)
        
    # !digIS
    if 'digIS'.lower() in sections:
        oname = join(odir,'digIS',genome)
        cmd = f"{digIS_exe} -i {fna} -g {gbk} -o {oname}"
        check(oname,cmd,'digIS',dry_run=dry_run)

    # !CrisprCasFinder
    if 'CrisprCasFinder'.lower() in sections:
        oname = join(odir,'CrisprCasFinder',genome)
        cmd = f"mkdir -p {odir} && cp {fna} {odir}/ && cd {odir} && sudo singularity exec -B $PWD /home-user/thliao/software/CRISPRCasFinder/CrisprCasFinder.simg perl /usr/local/CRISPRCasFinder/CRISPRCasFinder.pl -so /usr/local/CRISPRCasFinder/sel392v2.so -cf /usr/local/CRISPRCasFinder/CasFinder-2.0.3 -drpt /usr/local/CRISPRCasFinder/supplementary_files/repeatDirection.tsv -rpts /usr/local/CRISPRCasFinder/supplementary_files/Repeat_List.csv -cas -def G -out Result -in {genome}.fna && rm {genome}.fna"
        check(oname,cmd,'CrisprCasFinder',dry_run=dry_run)
    # !pseudofinder
    if 'pseudofinder'.lower() in sections:
        finalname = join(odir,'pseudofinder',genome,f'{genome}_nr_pseudos.gff')
        oname = join(odir,'pseudofinder',genome)
        #_odir = f"/mnt/ivy/thliao/project/coral_ruegeria/nanopore_processing/annotations/pseudofinder/"
        #odir = '/mnt/ivy/thliao/project/coral_ruegeria/data_processing/pub_dataset/pseudogenes/'
        cmd = f"mkdir -p {oname} ; python3 /home-user/thliao/software/pseudofinder/pseudofinder.py annotate -db /home-db/pub/protein_db/nr/v20230523/diamond_index/nr.dmnd -g {gbk} -t 4 -skpdb -op {oname}/{genome}_nr -di -l 0.8 --compliant -dip /home-db/pub/protein_db/nr/v20230523/diamond_index/diamond"
        check(finalname,cmd,'pseudofinder',dry_run=dry_run)
    # !antismash
    if 'antismash'.lower() in sections:
        oname = join(odir,'antismash',genome+'_out')
        cmd = f"{antishmash_exe} --cb-general --cb-knownclusters --cb-subclusters --asf --pfam2go --smcog-trees {gbk} --output-dir {oname} --databases {antishmash_db} -c 5 --genefinding-tool none"
        check(oname+f'/{genome}.gbk',cmd,'antismash',dry_run=dry_run)

    # !COG
    if 'COG'.lower() in sections:
        oname = join(odir,'COG_diamond',genome+'.out')
        cmd = f"diamond blastp -o {oname} -p 20 --db {COG_diamond_db} --query {faa}"
        check(oname,cmd,'COG_diamond',dry_run=dry_run)

    # !interproscan
    if 'interproscan'.lower() in sections:
        oname = join(odir,'ipr',f'{genome}.anno')
        cmd = f"""export LD_LIBRARY_PATH='' && {interpro_exe} --disable-precalc -i {faa} -d {dirname(oname)} -cpu {interpro_cpu} -pa -goterms -iprlookup -appl {interpro_application_list} """
        check(oname,cmd,'interpro',dry_run=dry_run)

    # !CRISPRCasFinder
    if 'CRISPRCas'.lower() in sections:
        oname = join(odir,'CRISPRCas',genome)
        cmd = f"{CRISPRCasFinder_exe} -in {fna} -cas -log -quiet -force -html -keep -so {CRISPRCasFinder_set1} -ms 10 -outdir {oname}"
        check(f"{oname}",cmd,'CRISPRCas',dry_run=dry_run)

    # !card
    if 'card'.lower() in sections:
        oname = join(odir,'blast','CARD',genome)
        ofile = f"{oname}/{genome}.blastout"
        cmd = f"blastp -num_threads 3 -query {faa} -db {CARD_db} -outfmt 6 > {ofile}"
        check(ofile,cmd,'card',dry_run=dry_run)

    #! VF
    if 'VF'.lower() in sections:
        oname = join(odir,'blast','VFDB',genome)
        ofile = f"{oname}/{genome}.blastout"
        cmd = f"blastp -query {faa} -db {VF_db} -num_threads 2 -max_target_seqs 100000 -outfmt 6 -evalue 1e-3 -out {ofile}"
        check(ofile,cmd,'VFDB',dry_run=dry_run)


    #! GI annotations
    if 'alien_hunter'.lower() in sections:
        oname = join(odir,'alienhunter',genome)
        ofile = f"{oname}/{genome}.out"
        cmd = f"{alienhunter_exe} {fna} {ofile}"
        check(ofile,cmd,'alienhunter',dry_run=dry_run)


    #! resfinder
    if 'resfinder'.lower() in sections:
        oname = join(odir,'resfinder',genome)
        cmd = f"{resfinder_exe} -o {oname} -s 'Other' --acquired -ifa {fna} -db_res {resfinder_db}"
        check(f'{odir}/{genome}',cmd,'resfinder',dry_run=dry_run)


    #!rgi
    if 'rgi'.lower() in sections:
        oname = join(odir,'rgi',genome)
        cmd = f"{rgi_exe} main --input_sequence {faa} --output_file {oname}/{genome}  -t protein --clean "
        check(f'{oname}/{genome}.json',cmd,'rgi',dry_run=dry_run)

    #! BacARscan
    if 'BacARscan'.lower() in sections:
        oname = join(odir,'BacARscan',genome)
        cmd = f"hmmsearch --tblout {oname}/{genome}.tab --acc --noali --notextw --cpu 5 {bacarscan_db} {faa}"
        check(f"{oname}/{genome}.tab",cmd,'BacARscan',dry_run=dry_run)

    # !plasmid
    if 'Plascad'.lower() in sections:
        oname = join(odir,'plascad',genome)
        cmd = f"/usr/bin/cp -f {fna} {oname}/ ; cd {oname} && Plascad -i ./{basename(fna)} -n > /dev/null 2>&1 ; rm {oname}/{genome}.fna"
        check(f"{oname}/{genome}_mob_unconj_plasmids_loc_sum.txt",cmd,'plascad',dry_run=dry_run)

    # ! abricate
    if 'abricate'.lower() in sections:
        oname = join(odir,'abricate',genome)
        for db in abricate_db_list:
            cmd = f"{abricate_exe} --db {db} {ffn} > {oname}/{genome}_{db}.out"
            check(f"{oname}/{genome}_{db}.out",cmd,'abricate',dry_run=dry_run)

    # ! WTP
    if 'wtp'.lower() in sections:
        oname = join(odir,'wtp_output')
        cmd = f"source activate wtp; export SINGULARITY_LOCALCACHEDIR=/mnt/ivy/thliao/project/ruegeria_prophage/sin_tmp; nextflow run 444thLiao/What_the_Phage -r v1.0.6 --cores 20 -profile local,singularity -w /mnt/ivy/thliao/tmp/ --fasta {fna} --output {oname} --dv --cachedir $SINGULARITY_LOCALCACHEDIR --databases /mnt/ivy/thliao/project/ruegeria_prophage/nextflow-autodownload-databases --cachedir /mnt/ivy/thliao/project/ruegeria_prophage/singularity_images"
        check(f"{oname}/{genome}/sample_overview_small.html",cmd,'WTP (phage detect)',dry_run=dry_run)
        
                
                
                
# case insensitive
software_list = ['KOFAMSCAN',
                 'interproscan',
                 'antismash',
                 'macsyfinder',
                 'abricate','card','VF','resfinder','rgi','BacARscan'
                 'ISESCAN','digIS',
                 'pseudofinder',
                 'WTP',
                 'alien_hunter','plascad']
software_list = [_.lower() for _ in software_list]

# file requirement
software_requirements = {'KOFAMSCAN':['protein'],'interproscan':['protein'],
                         'antismash':['genome'],
                         'bacarscanisescan':['protein'],
                 'macsyfinder':['genome'],
                 'abricate':["cds"],'card':['protein'],'VF':['protein'],'resfinder':['genome'],'rgi':['protein'],'BacARscan':['protein'],
                 'alien_hunter':["genome"],
                 'ISESCAN':['genome'],'digIS':['genome'],
                 'pseudofinder':["genbank"],'plascad':["genome"],
                 "COG":["protein"],
                 'WTP':['genome'],}
software_requirements = {k.lower():v for k,v in software_requirements.items()
                         }

def check_requirement(protein,genbank,cds,genome,
                      software_list):
    found = []
    found += ['protein'] if protein is not None else []
    found += ['genbank'] if genbank is not None else []
    found += ['genome'] if genome is not None else []
    found += ['cds'] if cds is not None else []

    run_s = []
    for f in software_list.split(','):
        f = f.lower()
        if f not in software_requirements:
            continue
        r = software_requirements[f]
        if r[0] not in found:
            print(f"{f} require {r[0]}. {f} is skipped.")
        run_s.append(f)
    return run_s

import click
@click.command()
@click.option("-p", '--protein', "protein", help='path of protein file')
@click.option("-o", '--output-dir', "output_dir", help='The directory you want to output to')
@click.option("-sl", '--software_list', "software_list",
              help=f' comma separated software list. default run all software including [' + ','.join(software_list) + ' ]',
              default=','.join(software_list))
@click.option("-gbk", '--genbank', "genbank", type=str,default=None,
              help='path of protein file')
@click.option("-cds", '--CDS', "cds", type=str,default=None,
              help='path of protein-encoding gene file')
@click.option("-fna", '--genome', "genome", type=str,default=None,
              help='path of genome files')
@click.option("-dry_run", '--dry_run', "dry_run", default=False,is_flag=True, show_default=True,
              help='default is No')
@click.option("-slurm",  "slurm",default=False,is_flag=True, show_default=True,
              help='use slurm to submit jobs.')
@click.option("-auto_imp",  "auto_implement",default=False,is_flag=True, show_default=True,
              help='give a protein path and automatically implement the other files. such as genbank,cds,fna')
def cli(protein,genbank,cds,genome,
        output_dir,
        software_list,
        dry_run,slurm,auto_implement):
    genome_name = basename(protein).rsplit('.',1)[0]
    if auto_implement:
        idir = dirname(realpath(protein))
        genbank = join(idir,genome_name+'.gbk')
        cds = join(idir,genome_name+'.ffn')
        genome = join(idir,genome_name+'.fna')
    run_s = check_requirement(protein,genbank,cds,genome,
                              software_list)
    run_cmd(protein,genome,genbank,cds,
            genome_name,odir=output_dir,
            sections = run_s,dry_run=dry_run)
    #print(executed_cmd)
    with open(join(output_dir,'cmd.log'),'w') as f1:
        f1.write('\n'.join(executed_cmd))
    if dry_run:
        return
    if slurm:
        sbatch_all(executed_cmd,thread_per_tasks=3,batch_size=5)
    
if __name__ == '__main__':
    cli()

## python single_anno.py -p /mnt/ivy/thliao/project/coral_ruegeria/nanopore_processing/canu_o/FI1/09_prokka/FI1.faa -o /mnt/ivy/thliao/project/coral_ruegeria/nanopore_processing/annotations/ -auto_imp -sl KOFAMSCAN

## python single_anno.py -p /mnt/ivy/thliao/project/coral_ruegeria/nanopore_processing/canu_o/BG7/09_prokka/BG7.faa -o /mnt/ivy/thliao/project/coral_ruegeria/nanopore_processing/annotations/ -auto_imp -sl interproscan






