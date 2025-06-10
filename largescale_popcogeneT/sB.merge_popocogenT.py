##### merged MCs and add previously MC
import pandas as pd
from collections import defaultdict
pre_MC = '/mnt/ivy/thliao/project/coral_ruegeria/Merged_popcogeneT/2307Ruegeria_MCs.tsv'
d = pd.read_csv(pre_MC,sep='\t')

genome2mc = dict(zip(d.iloc[:,0],d['MC']))
genome2mc = {k:v.split('.')[0] for k,v in genome2mc.items()}
mc2genomes = defaultdict(list)
for genome,mc in genome2mc.items():
    mc2genomes[mc.split('.')[0]].append(genome)


def get_number(pre_mc=None,existed_mcs = []):
    if pre_mc is None:
        nid = sorted(existed_mcs,key=lambda x:float(x.replace('MC','').replace('_','.').split('s')[0]))[-1]
        nid = int(nid.replace('MC','').replace('_','.').split('s')[0])
        new_mc =f"MC{int(nid+1)}"
        return new_mc
    else:
        #if pre_mc =='MC8_52':return 
        nid = [m for m in existed_mcs if m.startswith(pre_mc) and 's' in m]
        if not nid:
            return pre_mc+'s1'
        else:
            nid = sorted(nid,key=lambda x:int(x.split('s')[-1]))[-1]
            return pre_mc + f"s{int(nid.split('s')[-1])+1}"
        
syn_genome2mc = {}
syn_mc2genomes = defaultdict(set)
new_mc_otab = ['/mnt/ivy/thliao/project/coral_ruegeria/Merged_popcogeneT/output/group1_output_0.000355362.txt.cluster.tab.txt',
               '/mnt/ivy/thliao/project/coral_ruegeria/Merged_popcogeneT/output/group2_output_0.000355362.txt.cluster.tab.txt',
               '/mnt/ivy/thliao/project/coral_ruegeria/Merged_popcogeneT/output/group3_output_0.000355362.txt.cluster.tab.txt',]
for otab in new_mc_otab:
    d = pd.read_csv(otab,sep='\t')
    s2m = dict(zip(d['Strain'],d['Main_cluster']))
    m2s = defaultdict(list)
    mapping = {}
    for genome,mc in s2m.items():
        m2s[mc].append(genome)
    
    for s,now_mc in s2m.items():
        if now_mc in mapping: 
            continue
        ### in what situation, we need to 
        existed_mcs = list(mc2genomes)+list(set(mapping.values())) + list(syn_mc2genomes)
        if s in genome2mc: 
            pre_mc = genome2mc[s]
            pre_genomes_sameMC = [_ for _ in mc2genomes[pre_mc] if _ !=s]
            nowMCs_genomessameMC = set([s2m[_] for _ in pre_genomes_sameMC if _ in s2m])
            if len(pre_genomes_sameMC)==0:
                mapping[now_mc]= pre_mc
            elif len(nowMCs_genomessameMC)==1:
                ## the previous mc only contain a single genome
                ## genomes assigned to the same MC in the last time are all assigned to the same MC again
                mapping[now_mc]= pre_mc
                #mapping[now_mc]= get_number(existed_mcs=existed_mcs)
            elif len(nowMCs_genomessameMC) !=1:
                for mc in nowMCs_genomessameMC:
                    existed_mcs = list(mc2genomes)+list(set(mapping.values())) + list(syn_mc2genomes)
                    mapping[mc]= get_number(pre_mc,existed_mcs)
        else:
            others = [_ for _ in m2s[now_mc] if _ in genome2mc]
            if others:continue 
            mapping[now_mc]= get_number(existed_mcs=existed_mcs)
            
    for now_mc,pre_mc in mapping.items():
        syn_mc2genomes[pre_mc] = m2s[now_mc]
        for _ in m2s[now_mc]:
            syn_genome2mc[_] = pre_mc
with open('/mnt/ivy/thliao/project/coral_ruegeria/Merged_popcogeneT/1561Ruegeria_MCs.tsv','w') as f1:
    f1.write(f"Genome ID\tMC\tprevious MC used in Coral project\tEcosystem\n")
    for k,v in syn_genome2mc.items():
        f1.write(f"{k}\t{v}\t{genome2mc.get(k,'NA')}\t{g2eco.get(k,'NA')}\n")
        