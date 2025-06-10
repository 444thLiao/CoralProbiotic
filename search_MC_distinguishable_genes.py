import pandas as pd
import os
os.chdir('/mnt/ivy/thliao/project/coral_ruegeria/data_processing_rpoB/CRF_figures')

g2pop = {_.split(',')[0]:int(_.split(',')[-1] )
         for _ in open('./final_pop.txt').read().strip().split('\n') if not _.startswith('#') and len(_.split(','))==3}
pop2g = defaultdict(list)
for k,v in g2pop.items():
    pop2g[int(v)].append(k)
pop2gs = {
    10: 'G1',
    23: 'G2',
    4: 'G3',
    0: 'G4',
    3: 'S6',
    5: 'S1',
    7: 'S5',
    13: 'S2',
    19: 'S3',
    44: 'S4'}
    
from ete3 import Tree
tre = Tree('/mnt/ivy/thliao/project/coral_ruegeria/data_processing/of_out20230527/OG_tree/iqtree/OG_concat_656G.newick',3)
outgroups = tre.get_common_ancestor('GNM000156135,GNM900108275'.split(',')).get_leaf_names()
og_df = pd.read_csv('/mnt/ivy/thliao/project/coral_ruegeria/data_processing/of_out20230527/of_full/Results_Jun04/Orthogroups/Orthogroups.tsv',sep='\t',index_col=0)
gids = open('/mnt/ivy/thliao/project/coral_ruegeria/ms_figures/pruned_ids_fig1.txt').read().strip().split('\n')[2:]
gids = [_ for _ in gids if _ not in outgroups]
sub_df = og_df.loc[:,gids]
num_df = sub_df.applymap(lambda x: 0 if pd.isna(x)
                         or not x else len(x.split(',')))
candidate_genes = num_df.index[(num_df==1).all(1)]
print(len(candidate_genes),len(gids))


locus2nuc_seq = {}
for genome in gids:
    ffn = realpath(f"/mnt/ivy/thliao/project/coral_ruegeria/data_processing/all_faa/{genome}.faa").replace('.faa','.ffn')
    if not exists(ffn):print(genome)
    locus2nuc_seq.update({_.id:_ for _ in SeqIO.parse(ffn,'fasta')})

#! prepare large number of candidate genes/KO
anno_file = '/mnt/ivy/thliao/project/coral_ruegeria/data_processing/annotations/KEGG_latest/kegg_anno.tsv'
kegg_annodf = pd.read_csv(anno_file, sep='\t', index_col=0)
sub_annodf = kegg_annodf.loc[gids,:]
locus2ko = {locus: ko
            for ko, _d in sub_annodf.to_dict().items()
            for genome, locus_list in _d.items()
            for locus in str(locus_list).split(',')}

from os.path import *
import os
from bin.multiple_sbatch import sbatch_all
from Bio import SeqIO 
from collections import defaultdict,Counter
from tqdm import tqdm
import copy

odir = './primer_desgins_aims/'
if not exists(odir):
    os.system(f"mkdir -p {odir}")
og2locus = {}
for OG in list(candidate_genes):
    faa = f"/mnt/ivy/thliao/project/coral_ruegeria/data_processing/of_out20230527/of_full/Results_Jun04/Orthogroup_Sequences/{OG}.fa"
    records = [_ for _ in SeqIO.parse(faa,'fasta') 
               if _.id.split('_')[0] in gids]
    og2locus[OG] = [_.id for _ in records]
    aln = join(odir,f'{OG}.aln')
    
og2ko = defaultdict(list)
for og,l_list in og2locus.items():
    og2ko[og].extend([locus2ko[_] for _ in l_list if _ in locus2ko]  )
og2ko = {k:Counter(v) for k,v in og2ko.items()}
og2ko = {k:list(v)[0] for k,v in og2ko.items() if len(v)==1}
functional_og = [og for og in og2ko.items()]
print(len(functional_og))  

cmds = []
for OG in tqdm(functional_og):
    new_faa = join(odir,f"{OG}_425isolates.ffn")
    seqs = [copy.deepcopy(locus2nuc_seq[l])
            for l in og2locus[OG]]
    for r in seqs:
        r.name = r.description = ''
        r.id = r.id.split('_')[0]
    with open(new_faa,'w') as f1:
        SeqIO.write(seqs,f1,'fasta-2line')
    aln = new_faa.replace('.ffn','.aln')
    cmd = f"mafft --auto {new_faa} > {aln}"
    cmds.append(cmd)
sbatch_all(cmds,thread_per_tasks=10,batch_size=50,prefix_name='aln')


#! pairwise calculcate the similarity
from glob import glob
import numpy as np
import pandas as pd
import itertools
def get_pwdf(df):
    s2s_iden = defaultdict(dict)
    for s1,s2 in tqdm(itertools.combinations(df.index,2),
                      total=90100):
        seq1 = df.loc[s1,:]
        seq2 = df.loc[s2,:]
        iden = (seq1==seq2).sum()/df.shape[1]*100
        s2s_iden[s1][s2] = iden
        s2s_iden[s2][s1] = iden
    pw_df = pd.DataFrame.from_dict(s2s_iden)
    pw_df = 100-pw_df.fillna(100)
    pw_df = pw_df.reindex(pw_df.columns)
    return pw_df



gid2pop_i = {_.split(',')[0]:_.split(',')[2] for _ in open('./final_pop.txt').read().split('\n') 
         if len(_.split(','))==3 and _.split(',')[1].startswith('#') }
pop2gids = defaultdict(list)
for g,p in gid2pop_i.items():
    pop2gids[p].append(g)


for aln in tqdm(glob("primer_desgins_aims/*.aln")):
    records = list(SeqIO.parse(aln,'fasta'))
    seqs = np.array([list(_.seq) for _ in records])
    alndf = pd.DataFrame(seqs,index =[_.id for _ in records] )
    final_c = alndf.columns[(alndf!='-').all(0)]
    subdf = alndf.loc[:,final_c]
    #print(subdf.shape[1])
    pw_df = get_pwdf(subdf)
    pw_df.to_csv(aln.replace('.aln','.dist'),sep='\t',index=1)
    clus.fit(pw_df.values)
    
clus = AgglomerativeClustering(n_clusters =10,
                               affinity='precomputed',linkage='single')
clus.fit(pw_df.values)    
clus.labels 

from sklearn.cluster import DBSCAN
from sklearn.neighbors import NearestNeighbors

og2pop2ratio = {}
for d in tqdm(glob("primer_desgins_aims/*.dist")):
    og = d.split('/')[-1].split('_')[0]
    pw_df = pd.read_csv(d,sep='\t',index_col=0)
    nn = NearestNeighbors(metric='precomputed')
    nn.fit(pw_df.values)    
    pop2ratio = defaultdict(list)
    for g,p in gid2pop_i.items():
        other_g = [_ for _ in pop2gids[p][::] if _ !=g]
        num = len(other_g)
        dist,idx = nn.kneighbors(pw_df.loc[[g],:],n_neighbors=num)
        found_gids = pw_df.columns[idx.flatten()]
        overlapped_gids = set(other_g).intersection(set(found_gids))
        num_overlapped = len(overlapped_gids)
        ratio_overlapped = num_overlapped/num
        pop2ratio[p].append(ratio_overlapped)
    og2pop2ratio[og] = {k:np.mean(v) for k,v in pop2ratio.items()}
tmp_df = pd.DataFrame.from_dict(og2pop2ratio,orient='index')


og2len = {}
for d in tqdm(glob("primer_desgins_aims/*.aln")):
    og = d.split('/')[-1].split('_')[0]
    records = [_ for _ in SeqIO.parse(d,'fasta')]
    length = [len([_ for _ in r.seq if _!='-']) for r in records]
    og2len[og] = np.mean(length)
fdf = tmp_df.sort_values(list(tmp_df.iloc[:,:-1].min(0).sort_values().index),ascending=False)
fdf.loc[:,'length'] = [int(og2len[_]) for _ in fdf.index]

sorted_df = fdf.loc[fdf['length']>=1000,:].sort_values('length',ascending=False)
sorted_df = sorted_df.loc[(sorted_df.iloc[:,:10]>=0.80).all(1),:]
sorted_df.loc[:,'ko ID'] = [og2ko[_] for _ in sorted_df.index]
sorted_df.to_excel('./single_copy_OG.xlsx')

# sorted_df = pd.read_excel('./single_copy_OG.xlsx',index_col=0)
final_kos = list(sorted_df['ko ID'])
## left 43 kos

cmds = []
for og in sorted_df.index:
    aln = f"./primer_desgins_aims/{og}_425isolates.aln"
    cmd = f"FastTree {aln} > {aln.replace('.aln','.newick')} "
    if not exists(aln.replace('.aln','.newick')):
        cmds.append(cmd)
        
for og in ['OG0001239','OG0000833']:
    t = f"./primer_desgins_aims/{og}_425isolates.newick"
    cmd = f"FastRoot.py -i {t} -m MV -o {t.replace('.newick','.rooted.newick')}"
    cmd = f"format_newick.py rename -i {t.replace('.newick','.rooted.newick')} -o {t.replace('.newick','.renamed.newick')}"
    os.system(cmd)

from ete3 import Tree
from api_tools import *
for og in ['OG0001239','OG0000833']:
    t = f"./primer_desgins_aims/{og}_425isolates.renamed.newick"
    tre = Tree(t,3)
    name2g = {}
    for n in tre.traverse():
        for pop,g in pop2g.items():
            if set(n.get_leaf_names()) == set(g):
                name2g[n.name] = pop2gs[pop]
    text = to_color_range(name2g,{v: '#b71c1c' if v[0]=='G' else '#f4cccc' for k,v in name2g.items()})
    with open(f"./primer_desgins_aims/{og}_colorrange.txt",'w') as f1:
        f1.write(text)



sbatch_all(cmds,thread_per_tasks=10,batch_size=1,prefix_name='tree')
from api_tools import get_itoltree
for og in sorted_df.index:
    treefile = f"./primer_desgins_aims/{og}_425isolates.newick"
    get_itoltree(treefile,
                 outfile=f"./primer_desgins_aims/{og}_425isolates.png",
                 anno_files=['./final_pop.txt',
                             './genome2location_colorstrip.txt',
                             './genome2compartment_colorstrip.txt'])
sub_keggdf = kegg_annodf.loc[:,final_kos]
num_df = sub_keggdf.applymap(lambda x:0 if pd.isna(x) else str(x).count(',')+1 )
single_ko_annodf = sub_keggdf.loc[:,(num_df==1).sum(0)>=700] 
final_annodf = sorted_df.loc[sorted_df['ko ID'].isin(single_ko_annodf.columns)]
# final_annodf.to_excel('./single_ko_OG.xlsx')
# final_annodf = pd.read_excel('./single_ko_OG.xlsx',index_col=0)
## left 29 kos
used_locus2ko = {locus: ko
            for ko, _d in single_ko_annodf.to_dict().items()
            for genome, locus_list in _d.items()
            for locus in str(locus_list).split(',')}

locus2seq = {}
for g in tqdm(single_ko_annodf.index):
    ffn = realpath(f'/mnt/ivy/thliao/project/coral_ruegeria/data_processing/all_faa/{g}.faa').replace('.faa', '.ffn')
    n2seq = {n.id: n for n in SeqIO.parse(ffn, 'fasta') if n.id in used_locus2ko}
    locus2seq.update(n2seq)
locus2nucl = {}
for g in tqdm(single_ko_annodf.index):
    ffn = realpath(f'/mnt/ivy/thliao/project/coral_ruegeria/data_processing/all_faa/{g}.faa').replace('.faa', '.ffn')
    if not exists(ffn):
        ffn = realpath(f'/mnt/ivy/thliao/project/coral_ruegeria/data_processing/pub_dataset/prokka_o/{g}/{g}.ffn')
    n2seq = {n.id: n for n in SeqIO.parse(ffn, 'fasta') if n.id in used_locus2ko}
    locus2seq.update(n2seq)
import copy
for ko,col in tqdm(single_ko_annodf.iteritems()):
    locus_list = [str(_).split(',')[0] for _ in col.values if not pd.isna(_)]
    seqs = [copy.deepcopy(locus2seq[locus])
            for locus in locus_list]
    for s in seqs:
        s.name = s.description = ''
        s.id = s.id.split('_')[0]
    with open(f'/mnt/ivy/thliao/project/coral_ruegeria/data_processing_rpoB/CRF_figures/primer_kos/{ko}.fna', 'w') as f1:
        SeqIO.write(seqs, f1, 'fasta-2line')

cmds = []
for fna in glob("primer_kos/*.fna"):
    aln = fna.replace('.fna','.aln')
    cmd = f"mafft {fna} > {aln}; FastTree {aln} > {aln.replace('.aln','.newick')} "
    if not exists(aln.replace('.aln','.newick')):
        cmds.append(cmd)
        

for ko in tqdm(single_ko_annodf.columns):
    treefile = f"./primer_kos/{ko}.newick"
    get_itoltree(treefile,
                 outfile=f"./primer_kos/{ko}.png",
                 anno_files=['./final_pop.txt',
                             './genome2location_colorstrip.txt',
                             './genome2compartment_colorstrip.txt'])
    
    
    
## draw the positions
from ..nanopore_data.IS_results import parse_gbk

## show the primer gene
targeted_genomes = ['H9','O6','G7','AF3','AL3','H4','AO1','B4']

used_kos = 'K02621 K02314 K01965 K02519 K02133'.split(' ')
genome2locus2info = {}
for s in targeted_genomes:
    gbk = f'/mnt/ivy/thliao/project/coral_ruegeria/nanopore_processing/canu_o/{s}/09_prokka/{s}.gbk'
    genome2locus2info[s] = parse_gbk(gbk)


from api_tools.tk import *

nano_locus2ko = {}
for s in targeted_genomes:
    tab = f'/mnt/ivy/thliao/project/coral_ruegeria/nanopore_processing/canu_o/{s}/09_prokka/{s}.tab'
    locus2kegg_list,kegg2locus = parse_hmmscan(tab)
    for ko in used_kos:
        locus = kegg2locus[ko][0][0]
        nano_locus2ko[locus] = ko

import plotly.graph_objects as go
fig = go.Figure()
for idx,genome in enumerate(targeted_genomes):
    ref_fna = f"/mnt/ivy/thliao/project/coral_ruegeria/nanopore_processing/canu_o/{genome}/09_prokka/{genome}.fna"
    ref_records = {_.id:_ for _ in SeqIO.parse(ref_fna,'fasta')}
    seqs = sorted(ref_records.values(),key=lambda x:len(x.seq),reverse=True)
    base_start = [0,]
    for _ in list(seqs)[:-1]:
        base_start.append(base_start[-1]+len(_.seq))
    base_start = dict(zip([_.id for _ in seqs],base_start))
    
    for contig in ref_records:
        _s = base_start[contig]
        fig.add_scatter(x=[_s,_s+len(ref_records[contig].seq),None],
                        y=[idx,idx,None],
                        line=dict(width=10),
                        mode='lines',
                        showlegend=False,)
        
    for l,ko in nano_locus2ko.items():
        if ko not in used_kos:continue
        if not l.startswith(genome): continue
        info = genome2locus2info[genome][l]
        _s = base_start[info[0]]
        pos = (info[1] + info[2])/2 + _s
        if ko not in used_kos:
            fig.add_scatter(x=[pos],y=[idx+0.5],text=ko,mode='text',legendgroup=ko,name=ko,
                            textposition='top left',)
            used_kos.append(ko)
        else:
            fig.add_scatter(x=[pos],y=[idx+0.5],text=ko,mode='text',legendgroup=ko,name=ko,showlegend=False,
                            textposition='top left',)
        fig.add_annotation(x=pos, 
                           y=idx,
                           text='',
                           showarrow=True,arrowhead=5)
fig.layout.yaxis.ticktext = list([f"{_} (MC{g2pop[_]})" if _ in g2pop else _ for _ in targeted_genomes])     
fig.layout.yaxis.tickvals = list(range(len(targeted_genomes)))
fig.write_html(f'/mnt/ivy/thliao/project/coral_ruegeria/data_processing_rpoB/CRF_figures/show_primer_kos/5target_kos.html')



    # fig.write_html(f'./show_primer_kos/{genome}.html')
    # fig.write_image(f'./show_primer_kos/{genome}.png')