
#####
# popcogeneT key
# automatic select MC and its genomes to perform another popcogeneT
from ete3 import Tree
from tqdm import tqdm
import pandas as pd
from collections import defaultdict
import re
from glob import glob
import os
os.chdir('/mnt/ivy/thliao/project/coral_ruegeria/Merged_popcogeneT')
df = pd.read_csv('./2305Ruegeria_MCs.tsv',sep='\t',index_col=0)
odir = '/mnt/ivy/thliao/project/coral_ruegeria/Merged_popcogeneT/updated_popcogenT20241022_b'
tre = Tree('/mnt/ivy/thliao/project/coral_ruegeria/Merged_popcogeneT/phylogeny/bac120/bac120_MVrooted_2305Ruegeria.newick')
final_mcdf = f'/mnt/ivy/thliao/project/coral_ruegeria/Merged_popcogeneT/2305Ruegeria_MCs.tsv'
######## change above params


g2MC = df['MC'].to_dict()
n2node = {_.name:_ for _ in tre.get_leaves()}
# g = open('/mnt/ivy/thliao/project/coral_ruegeria/data_processing/20231227/ruegeria_filtered.list').read().strip().split('\n')
new_gids = [_ for _ in tre.get_leaf_names() if _ not in g2MC ]

mc2g = defaultdict(list)
for g,mc in g2MC.items():
    mc2g[mc].append(g)
c = 1
gid2group = defaultdict(list)
group2mc = defaultdict(list)
for g in tqdm(new_gids):
    n = n2node[g]
    if n.name not in gid2group:
        gid2group[n.name] = c
        c+=1
    other_nodes = sorted(tre.get_leaves(),key=lambda x:x.get_distance(n,topology_only=True))[1:]
    top5MC = []
    for _n in other_nodes:
        if _n.name not in g2MC:
            if _n.name in gid2group:
                gid2group[n.name] = gid2group[_n.name]
            else:
                gid2group[_n.name] = gid2group[n.name]
            continue
        mc = g2MC[_n.name]
        if mc in top5MC:
            continue
        else:
            top5MC.append(mc)
            if len(top5MC)==5:
                break
    group2mc[gid2group[n.name]].extend(top5MC)
group2mc = {k:set(v) for k,v in group2mc.items()}

import random
for g,mc_l in group2mc.items():
    genomes = []
    for mc in mc_l:
        gs = mc2g[mc]
        if len(gs)>=50:
            gs = random.sample(gs,k=50)
        genomes.extend(gs)
    g2 = [k for k,v in gid2group.items() if v == g]
    genomes = set(genomes+g2)
    print(g,len(genomes))

#### ln some protein files
from os.path import exists
for i in new_gids:
    if not exists(f"/mnt/ivy/thliao/project/coral_ruegeria/Merged_popcogeneT/inproteins/{i}.faa") and exists(f"/mnt/ivy/thliao/project/coral_ruegeria/Merged_popcogeneT/tmp_inproteins/{i}.faa"):
        os.system(f" ln -sf `realpath /mnt/ivy/thliao/project/coral_ruegeria/Merged_popcogeneT/tmp_inproteins/{i}.faa` /mnt/ivy/thliao/project/coral_ruegeria/Merged_popcogeneT/inproteins/{i}.faa")
assert all([exists(f"/mnt/ivy/thliao/project/coral_ruegeria/Merged_popcogeneT/inproteins/{i}.faa") for i in new_gids])


#### for each group, perform a small popcogeneT
from os.path import *
import os
cmds = []
idir = '/mnt/ivy/thliao/project/coral_ruegeria/Merged_popcogeneT/inproteins'
for g,mc_l in group2mc.items():
    genomes = [_ for mc in mc_l for _ in mc2g[mc]]
    g2 = [k for k,v in gid2group.items() if v == g]
    genomes = set(genomes+g2)
    _odir = f"{odir}/subg{g}"
    if not exists(_odir):
        os.makedirs(_odir)
    genome_dir = f"{_odir}/ingenomes"
    if not exists(genome_dir):os.makedirs(genome_dir)
    for gid in genomes:
        fna = realpath(join(idir,f"{gid}.faa")).replace('.faa','.fna')
        if not exists(fna):
            i = '/mnt/ivy/thliao/project/coral_ruegeria/Merged_popcogeneT/ingenomes'
            fna = realpath(join(i,f"{gid}.fna"))
        os.system(f"ln -sf {fna} {genome_dir}/{gid}.fna")
    with open(f"{_odir}/known_mc.tsv",'w') as f1:
        for g in genomes:
            if g in g2MC:
                f1.write(f"{g}\t{g2MC[g]}\n")
            else:
                f1.write(f"{g}\tnew\n")
    cmd = f"python /home-user/thliao/software/PopCOGenT/src/PopCOGenT/get_alignment_and_length_bias.py --genome_dir {genome_dir} --genome_ext .fna --alignment_dir {_odir}/proc/ --mugsy_path /home-user/software/mugsy_v1r2.3/mugsy --mugsy_env /home-user/software/mugsy_v1r2.3/mugsyenv.sh --base_name output --final_output_dir {_odir}/output --num_threads 20 --keep_alignments --slurm "
    cmds.append(cmd)

#### for all new sequencing genomes , perform a small popcogeneT
genomes = list(gid2group)
_odir = f"{odir}/all"
if not exists(_odir):
    os.makedirs(_odir)
genome_dir = f"{_odir}/ingenomes"
if not exists(genome_dir):
    os.makedirs(genome_dir)
for gid in genomes:
    fna = realpath(join(idir,f"{gid}.faa")).replace('.faa','.fna')
    if not exists(fna):
        i = '/mnt/ivy/thliao/project/coral_ruegeria/Merged_popcogeneT/ingenomes'
        fna = realpath(join(i,f"{gid}.fna"))
    os.system(f"ln -sf {fna} {genome_dir}/{gid}.fna")
with open(f"{_odir}/known_mc.tsv",'w') as f1:
    for g in genomes:
        if g in g2MC:
            f1.write(f"{g}\t{g2MC[g]}\n")
        else:
            f1.write(f"{g}\tnew\n")
cmd = f"python /home-user/thliao/software/PopCOGenT/src/PopCOGenT/get_alignment_and_length_bias.py --genome_dir {genome_dir} --genome_ext .fna --alignment_dir {_odir}/proc/ --mugsy_path /home-user/software/mugsy_v1r2.3/mugsy --mugsy_env /home-user/software/mugsy_v1r2.3/mugsyenv.sh --base_name output --final_output_dir {_odir}/output --num_threads 20 --keep_alignments --slurm "
cmds.append(cmd)

all_cmds = []
for c in cmds:
    os.system(c)
    all_cmds += open('./cmds').read().strip().split('\n')
all_cmds = [_ for _ in all_cmds if _]

from bin.multiple_sbatch import sbatch_all
sbatch_all(all_cmds,thread_per_tasks=4,prefix_name='mugsy',batch_size=500)


for a in cmds:
    c = a.split(' ')
    d = c[c.index('--final_output_dir')+1]
    if not exists(d):
        os.makedirs(d)
    try:
        if not exists(f"{d}/output.length_bias.txt"):
            os.system(a.replace(' --slurm',' --merged'))
        if not exists(f"{d}/output_0.000355362.txt.cluster.tab.txt"):
            os.system(f"/home-user/thliao/anaconda3/envs/PopCOGenT/bin/python3.6 /home-user/thliao/software/PopCOGenT/src/PopCOGenT/cluster.py --base_name output --length_bias_file {d}/output.length_bias.txt --output_directory {d}/ --infomap_path /home-user/software/PopCOGenT/Infomap/Infomap > {d}/cluster.log 2>&1")
    except:
        print('not completed')



## merge all output results
nf = f"{odir}/all/output/output_0.000355362.txt.cluster.tab.txt"
newdf = pd.read_csv(nf,sep='\t',index_col=0)
this_time_g2MC = newdf['Main_cluster'].to_dict()
this_time_MC2g = defaultdict(list)
for g,mc in this_time_g2MC.items():
    this_time_MC2g[mc].append(g)

_df = pd.read_csv('/mnt/ivy/thliao/project/coral_ruegeria/Merged_popcogeneT/2307Ruegeria_MCs.tsv',sep='\t',index_col=0)
largest_mc = sorted([int(re.findall('\d+',_)[0]) for _ in _df['MC'].unique()])[-1]
new_g2mc = {}
for f in sorted(glob(f"{odir}/sub*/known_mc.tsv")):
    knownmc = pd.read_csv(f,header=None,sep='\t')
    nf = join(dirname(f),'output','output_0.000355362.txt.cluster.tab.txt')
    #if not exists(nf):continue
    newdf = pd.read_csv(nf,sep='\t',index_col=0)
    g2mc = dict(zip(knownmc[0],knownmc[1]))
    newgids = [g for g,mc in g2mc.items() if mc =='new']
    for g in newgids:
        if g in new_g2mc:
            continue
        cid = newdf.loc[g,'Cluster_ID']
        snewdf = newdf.loc[newdf['Cluster_ID']==cid,:]
        existedgids = set(snewdf.index).difference(newgids)
        if not existedgids:
            new_mc = f"MC{largest_mc+1}"
            if len(this_time_MC2g[this_time_g2MC[g]])==1:
                new_g2mc[g] = new_mc
            else:
                for _ in this_time_MC2g[this_time_g2MC[g]]:
                    new_g2mc[_] = new_mc
            largest_mc+=1
        else:
            # assign previous MC to new sequencing genomes
            new_g2mc[g] = g2mc[list(existedgids)[0]]

mdf = pd.read_csv(final_mcdf,sep='\t',index_col=0)
for k,v in new_g2mc.items():
    mdf.loc[k,'MC'] = v
mdf.to_csv(final_mcdf,sep='\t',index=1)


################ no need to run the sB.



