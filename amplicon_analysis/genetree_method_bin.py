"""
Define a gene tree-based MC assignment

"""

from glob import glob
from os.path import *
import os
from Bio import SeqIO
from ete3 import Tree
import numpy as np
from subprocess import check_call
import pandas as pd
import json
import itertools
from collections import Counter,defaultdict
from glob import glob
from bin.format_newick import renamed_tree
import warnings
warnings.filterwarnings('ignore')
from scipy.spatial.distance import cdist
from tqdm import tqdm

def cal_jaccard(tre,g2pop,threshold=1):
    """
    Given a tree and a dictinoary which assign each tip to a particular group.
    0. transform the dictionary into a dataframe(row: groups, cols:tips)
    1. prune tree according to the dictinoary
    2. transform subtree into a node2tips dataframe (row: internal node, cols:tips)
    3. calculate the distance between "groups" and "internal node"
        distance metric (0-1) 
            K: Number of tips in this mc & in this node
            V: Number of tips in this node (might contain tips that assigned to the other MC)
            calculate the ratio (K/V)
        0: no descending tips of this nodes belong to this mc. 
        1: all descending tips of this nodes belong to this mc. 
    4.  
    # note: it can not differentiate the node containing ASV at the basal or in the inside.

    """
    def dfun(a,b):
        ns = ((a==b) & (a==1)).sum()
        return ns/b.sum()
        
    ### convert g2pop and mc2genomes
    mc2genomes = defaultdict(dict)
    for g,mc in g2pop.items():
        mc2genomes[mc][g] = 1
    mc2genomes['clade1'] = {k:1 
                            for k in 'AL10,AJ11,AM3,AL9,GNM6454,GNM4293,AD4,AJ8,AJ4,AJ10,GNM3532,BA3,AC2,AJ7,AD3,AC5,AC12,M5,AD7,CI4,AK9,BF1,L5,L7,L1,AD8,M4,L4,AK12,N6,AE8,BJ8,AE2,BJ9,AG8,M7,L2'.split(',')}
    mc2genome_df = pd.DataFrame.from_dict(mc2genomes).T.fillna(0)
    ### prune 
    # 1.benefit: also removed the node containing ASV as the basal and the isolates
    # 2.not raise error
    subtre = tre.copy()
    n = set(subtre.get_leaf_names())
    subtre.prune([k for k in g2pop if k in n])
    name2node = {n.name:n for n in subtre.traverse()}
    node2tips = {n.name:{k:1 
                         for k in n.get_leaf_names()} 
                 for n in subtre.traverse()}
    node2tip_df = pd.DataFrame.from_dict(node2tips).T.fillna(0)
    node2tip_df = node2tip_df.reindex(columns=mc2genome_df.columns).fillna(0)
    ### 

    d = pd.DataFrame(cdist(mc2genome_df,node2tip_df,metric=dfun),
                     index=mc2genome_df.index,
                     columns=node2tip_df.index)
    # row is MC, columns is each internal node (defined by its descending nodes)
    # values are the similarities  
    mc2best_lca = {}
    for mc,row in d.iterrows():
        s = pd.DataFrame(row.sort_values())
        ss = s.loc[s.iloc[:,0]>=threshold]
        ss.loc[:,'dist'] = [subtre.get_distance(name2node[_],topology_only=True) for _ in ss.index]
        ss.loc[:,'num'] = [len(node2tips[_]) for _ in ss.index]
        ss = ss.sort_values('num',ascending=False)
        if len(mc2genomes[mc])!=1:
            # if this mc contain multiple genomes, it should remove the tip nodes.
            ss = ss.loc[ss['num']!=1,:]
        if ss.shape[0]==0:
            continue
            #mc2best_lca[mc] = s.index[0]
        else:
            mc2best_lca[mc] = [name2node[k] for k in ss.index]
    return mc2best_lca,name2node

# redundant
def refine_nodes(nodes,fullname2nodes,g2pop):
    """
    remove nodes that are descendent to the other nodes
    """
    sorted_nodes = sorted(nodes,
               key=lambda n: len(n.get_leaf_names()),reverse=True)
    s = []
    for n in sorted_nodes:
        if n.is_leaf():
            s.append(n)
            continue
        a,b = fullname2nodes[n.name].children
        _n = [_ for _ in [a,b] if _.is_leaf()]
        if _n:
            # if one of the children is leaf node and it has MC, keep it.
            _n = _n[0]
            if _n in g2pop:
                s.append(n)
        else:
            # if both children are internal nodes
            if any([all([_ not in g2pop 
                         for _ in subn.get_leaf_names()])
                    for subn in [a,b]]):
                # if any of the internal node contain all ASVs, exclude it 
                continue
            else:
                s.append(n)
    
    final_nodes = []
    while len(s)!=0:
        f = s.pop(0)
        final_nodes.append(f)
        desc = set([n for n in f.traverse()]+[f])
        s = set(s).difference(desc)
        if len(s)==0:
            break
        s = sorted(list(s),
               key=lambda x: len(x.get_leaf_names()),reverse=True)
    return final_nodes


def tight_cluster(core_node,flexible_node,prefix='ASV',ratio_threshold=90):
    tips = core_node.get_leaves()
    tip_ASV_names = [_ for _ in core_node.get_leaf_names() if _.startswith(prefix)]
    dist_tips = [a.get_distance(b) 
                 for a,b in itertools.combinations(tips,2)]
    outasv = [_ for _ in flexible_node.get_leaves() if _.name.startswith(prefix) and _.name not in tip_ASV_names]
    outasv2dist = {asv.name:np.mean([asv.get_distance(t) 
                                     for t in tips]) 
                   for asv in outasv}
    threshold = np.max(dist_tips) - (np.max(dist_tips) - np.min(dist_tips) )/100 * ratio_threshold
    finalasv = {asv:(dist,threshold) 
                for asv,dist in outasv2dist.items() 
                if dist <=threshold}
    for asv in tip_ASV_names:
        finalasv[asv] = None
    return finalasv

def get_dist_matrix(tpath,g2pop):
    tre = Tree(tpath)
    tre.resolve_polytomy()
    named = renamed_tree(tre)
    for n in named.traverse():
        if n.is_root():
            n.name = 'ROOT'
        elif not n.is_leaf():
            n.name = n.name+'_INODE'
    n2nodes = {n.name:n for n in named.traverse()}
    #print(f"Get distance matrix between groups and all internal nodes")
    group2IN,sub_name2node = cal_jaccard(named,g2pop)
    # n2nodes: from tree containing ASV and ref
    # sub_name2node: from tree containing ref only
    # group2IN: Internal nodes which all descending nodes belong to a specific group of reference genomes.
    
    # group2IN = {MC:refine_nodes(nodes,fullname2nodes=n2nodes,g2pop=g2pop) 
    #             for MC,nodes in group2IN.items()}
    return group2IN,n2nodes,named

def get_asv_group(group2IN,n2nodes,named,pop2g,
                  prefix='ASV',ratio_threshold=None):
    """_summary_
    get candidate 
    Args:
        tpath (_type_): _description_
        g2pop (_type_): _description_
        prefix (str, optional): _description_. Defaults to 'ASV'.
        percentile_threshold (int, optional): _description_. Defaults to 100.

    Returns:
        _type_: _description_
    """
    # nodes in here are from the tree pruning ASVs
    asv_l = set([_.name 
                    for _ in named.get_leaves() 
                    if _.name.startswith(prefix)])
    ### 
    processed_nodes = set()
    asv2groups = {}
    asv2reasons = defaultdict(list)
    for g,nodes in tqdm(group2IN.items()):
        # g: MC
        # nodes: internal nodes (without ASVs)
        for node_noASV in nodes:
            node_withASV = n2nodes[node_noASV.name] # nodes: internal nodes (with ASVs)
            processed_nodes.add(node_withASV)
            static_node_withASV = node_withASV
            tips = set(node_noASV.get_leaf_names())
            if ratio_threshold is not None and len(tips)>=2:
                # for retrieve sister lineage
                ## keep get its parent node until it meet a parent node that will add new MC. Stop
                while 1:
                    nn = node_withASV._up
                    if nn in processed_nodes:
                        node_withASV = nn
                        continue
                    if all([_.startswith(prefix) 
                            for _ in set(nn.get_leaf_names()).difference(pop2g[g])]):
                        # all tips descending this node either from this pop or ASV
                        node_withASV = nn
                    else:
                        # unclassfied genome or genome from the other pop found
                        break
            if static_node_withASV is not node_withASV:
                # it is another node now and percentils_threshold != 100
                asv2selected = tight_cluster(static_node_withASV,
                                             node_withASV,prefix=prefix,
                                             ratio_threshold=ratio_threshold)
                for asv,r in asv2selected.items():
                    asv2groups[asv] = g
                    if r is None:
                        asv2reasons[asv].append(('100% Internal node',g))
                    else:
                        asv2reasons[asv].append((f'sister. dist: {r[0]} < {r[1]}',g))
            else:
                tips = set(static_node_withASV.get_leaf_names())
                final_asv = tips & asv_l
                for asv in final_asv:
                    asv2reasons[asv].append(('100% Internal node',g))
                    asv2groups[asv] = g
    remaining_asv = set(asv_l).difference(set(asv2groups.keys()))
    for a in remaining_asv:
        asv2groups[a] = 'Unassigned'
    return asv2groups,asv2reasons
    



from glob import glob
from os.path import *
import os
from Bio import SeqIO
from ete3 import Tree
import numpy as np
from subprocess import check_call
import pandas as pd
import json

idir = '/mnt/ivy/thliao/project/coral_ruegeria/data_otheramplicons/markers_trees_add1302xyfeng/fasttree/'
gene2besttre = {'nirS2':f"{idir}/nirS2_trimmed.nucl.formatted_newick",
                "ATP5B":f"{idir}/ATP5B_trimmed.nucl.formatted_newick",
                "parC":f"{idir}/parC_trimmed.nucl.formatted_newick",
                #"pccA":f"{idir}/pccA_trimmed.nucl.fasttree.newick",
                }

def batch_iter(iter, batch_size):
    # generating batch according batch_size
    iter = list(iter)
    n_iter = []
    batch_d = 0
    for batch_u in range(0, len(iter), batch_size):
        if batch_u != 0:
            n_iter.append(iter[batch_d:batch_u])
        batch_d = batch_u
    n_iter.append(iter[batch_d : len(iter) + 1])
    return n_iter

import subprocess
def run(args):
    if isinstance(args, str):
        cmd = args
        log = '/dev/null'
    else:
        cmd, log = args
    try:
        # subprocess.check_output(cmd, shell=1)
        check_call(cmd,
                   shell=True,
                   stdout=open(log, 'w'))

    except subprocess.CalledProcessError as e:
        pass
        # print('error', e.output)
    if log != '/dev/null':
        t = open(log, 'r', newline='\n').read().replace('\r', '\n')
        with open(log, 'w') as f1:
            f1.write(t)

if __name__ == '__main__':
    ############## parameters
    # ref tree should include all reference tips in the tree (target_tree) generated using ASV and reference
    #ref2group = ''
    #target_tree = ""
    # target_tree should be unrooted one. The process will root it with MV method.
    ungroupped_seq_prefix = 'ASV'
    rep_otu = '/mnt/ivy/thliao/project/coral_ruegeria/data_otheramplicons/outplant_amplicon/merged_parC_20241220/results/Pdada2_output/Positive_rep.fasta'
    gene = 'parC'
    mode = 'multiprocessing' # or multiprocessing
    ############### example parameters

    ## 1.0 placement
    ref_aln = f'{idir}/../{gene}_trimmed.nucl.aln'
    ref_newick = gene2besttre[gene]
    final_odir = join(dirname(rep_otu),'placement')
    if not exists(final_odir):
        os.makedirs(final_odir)
    pre_sequences = list(SeqIO.parse(ref_aln,'fasta'))
    ref_phy = f"{final_odir}/ref_trimmed.phy"
    records = list(SeqIO.parse(ref_aln,'fasta'))
    rs = []
    names = {}
    for r in records:
        if r.id in ['GNM000473225_04317','GNM900116455_04551']:
            continue
        if r.id.startswith('GCA') or r.id.startswith('GCF'):
            name = r.id
        else:
            name = r.id.split('_')[0]
        if name in set(names.values()):
            continue
        names[r.id] = name
        rs.append(f"{name}{' '*10}{r.seq.upper()}")
    n1 = len(rs)
    n2 = len(records[0].seq)
    rs = [f"{n1} {n2}"] + rs
    with open(ref_phy,'w') as f1:
        f1.write('\n'.join(rs))
    os.system(f"cp {ref_newick} {join(final_odir,'placement.newick')}")
    fixed_ref_newick = join(final_odir,'placement.newick')
    tre = Tree(fixed_ref_newick,3)  
    # if not exists(f"{final_odir}/epa_result.jplace"):
    if exists(join(final_odir,'papara_log.aln')):
        os.system(f"rm -r {join(final_odir,'papara_log.aln')}") 
    cmd1 = f"cd {final_odir} && papara -t {fixed_ref_newick} -s {ref_phy} -q {rep_otu} -r -n aln -j 30"
    check_call(cmd1,shell=1)
    cmd2 = f"cd {final_odir} && /home-user/thliao/anaconda3/bin/epa-ng --split {ref_phy} papara_alignment.aln --redo"
    check_call(cmd2,shell=1)

    ### tree
    cmds = []
    ofile_name = join(dirname(rep_otu),'placement')
    if gene not in gene2besttre:
        raise IOError('weird') #continue
    if 'formatte' in gene2besttre[gene]:
        tre = Tree(gene2besttre[gene],3)
    else:
        tre = Tree(gene2besttre[gene])
    ref_ids = tre.get_leaf_names()
    all_aln = {_.id:_ for _ in SeqIO.parse(f"{ofile_name}/query.fasta",'fasta')}
    for k,v in all_aln.items():
        if ';' in v.id:
            v.id = v.id.split(';')[0]

    used_asv_ids = list(all_aln)
    all_aln.update({_.id:_ for _ in SeqIO.parse(f"{ofile_name}/reference.fasta",'fasta')})

    tre.prune([_ for _ in tre.get_leaf_names() if _ in all_aln])
    tre.unroot()
    tre.write(outfile=f"{dirname(rep_otu)}/guide_ref.newick",format=3)

    n_iter = batch_iter(used_asv_ids,len(used_asv_ids)//10)
    print(len(used_asv_ids)//10)

    for idx,n in tqdm(enumerate(n_iter),total=len(n_iter)):
        aln_subset = [all_aln[asv] for asv in n] + [all_aln[_] for _ in ref_ids]

        tre_dir = realpath(f"{ofile_name}/../asv_fasttree")
        if not exists(tre_dir):
            os.makedirs(tre_dir)
        with open(f'{tre_dir}/ASV_{idx}.aln','w') as f:
            SeqIO.write(aln_subset,f,'fasta-2line')
        cmd2 = f"/home-user/thliao/anaconda3/bin/FastTreeMP -gtr {tre_dir}/ASV_{idx}.aln > {tre_dir}/ASV_{idx}.tree"
        #cmds.append(cmd2)

        tre_dir = realpath(f"{ofile_name}/../asv_IQtree")
        if not exists(tre_dir):
            os.makedirs(tre_dir)
        with open(f'{tre_dir}/ASV_{idx}.aln','w') as f:
            SeqIO.write(aln_subset,f,'fasta-2line')
        cmd = f"iqtree -nt 3 -m MFP -mset GTR,HKY85,K80 -mrate E,I,G,I+G -mfreq FU -g {dirname(rep_otu)}/guide_ref.newick -pre {tre_dir}/ASV_{idx} -s {tre_dir}/ASV_{idx}.aln -fast "
        cmds.append(cmd)
    
    if mode=='slurm':
        from bin.multiple_sbatch import sbatch_all
        sbatch_all(cmds,thread_per_tasks=3,prefix_name='tree')
    elif mode=='multiprocessing':
        import multiprocessing as mp
        with mp.Pool(processes=5) as tp:
            r = list(tqdm(tp.imap(run, cmds), total=len(cmds)))
    ## 3.0 assign MC


    g2pop = pd.read_csv('/mnt/ivy/thliao/project/coral_ruegeria/Merged_popcogeneT/2307Ruegeria_MCs.tsv',sep='\t',index_col=0)
    g2pop = g2pop['MC'].to_dict()
    pop2g = defaultdict(list)
    for k,v in g2pop.items():
        pop2g[v].append(k)
            
    
    tre_dir = join(dirname(rep_otu),'asv_IQtree')
    _asv2groups = {}
    for tpath in tqdm(glob(f"{tre_dir}/ASV_*.treefile")):
        ntpath = tpath.replace('.treefile','.MVroot.iqtree_newick')
        if not exists(ntpath):
            os.system(f"python /home-user/thliao/anaconda3/envs/soft/bin/FastRoot.py -i {tpath} -m MV -o {ntpath} > /dev/null")
        group2IN,n2nodes,named = get_dist_matrix(ntpath,g2pop)
        asv2groups,asv2reasons = get_asv_group(group2IN,n2nodes,named,pop2g,'ASV',ratio_threshold=80)
        _asv2groups.update(asv2groups)
    with open(f'{tre_dir}/assigned_asv.txt','w') as f1:
        for asv,group in _asv2groups.items():
            f1.write(f'{asv}\t{group}\n')



