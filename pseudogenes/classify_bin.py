
import os
from os.path import exists,basename,getsize
import pandas as pd
import portion as P
from tqdm import tqdm
from Bio import SeqIO
from collections import defaultdict
from geneblocks import DiffBlocks
import time

############################# requirement ###########################################
pseudo_table_file = "/mnt/ivy/thliao/project/coral_ruegeria/nanopore_processing/annotations/pseudofinder/pseudofinderALL_w_annot_MERGEDcontinuous.tsv"
IS_file = "/mnt/ivy/thliao/project/coral_ruegeria/nanopore_processing/canu_o/IS_info.tsv"
genomic_component_file = '/mnt/ivy/thliao/project/coral_ruegeria/nanopore_processing/annotations/DFD/s2out/100Genomes_pos.tsv'

all_genomic_component_df = pd.read_csv(genomic_component_file,sep='\t',index_col=0)
all_IS_df = pd.read_csv(IS_file,sep='\t',index_col=0)
all_pseudo_df = pd.read_csv(pseudo_table_file,sep='\t',index_col=0)

############################# params ###########################################
# query_fna = "/mnt/ivy/thliao/project/coral_ruegeria/nanopore_processing/canu_o/AG8/09_prokka/AG8.fna"
# pseudo_fna = '/mnt/ivy/thliao/project/coral_ruegeria/nanopore_processing/annotations/pseudogenes/AG8_pseudos_extracted.fasta'

# subject_fna = '/mnt/ivy/thliao/project/coral_ruegeria/nanopore_processing/canu_o/BJ9/09_prokka/BJ9.fna'

# odir = "/mnt/ivy/thliao/project/coral_ruegeria/nanopore_processing/annotations/pseudogenes/pairwise/"



############################# functions ########################################
def load_contig2seq(fna_list):
    if type(fna_list) == str:
        fna_list = [fna_list]
    contig2seq = {}
    for fna in fna_list:
        for record in SeqIO.parse(fna,'fasta'):
            contig2seq[record.id] = record
    return contig2seq

def get_fna(s,e,c,f,contig2seq):
    full_seq = contig2seq[c]
    seq = full_seq[s:e]
    with open(f,'w') as f1:
        SeqIO.write([seq],f1,'fasta-2line')
    return seq

def get_pseudo(pseudo_df,c,s,e):
    s,e = min([s,e]),max([s,e])
    i = P.closed(s,e)
    sub_df = pseudo_df.loc[pseudo_df['contig']==c,:]
    idx = []
    for _idx, row in sub_df.iterrows():
        if P.closed(row['start'],row['end']).overlaps(i):
            idx.append(_idx)
    return idx

def anno_pseudo(blastdf,sub_pseudo_df):
    tqdm.write('annotating blast results...')
    xs = []
    for _,row in tqdm(blastdf.iterrows()):
        contig,start,end = row[[1,8,9]].values.tolist()
        n = f"{contig}:{start-1}-{end}"
        if n == row[0]:
            xs.append('same')
            continue
        idx = get_pseudo(sub_pseudo_df,contig,start,end)
        if len(idx) !=0:
            xs.append('also pseudo')
        else:
            xs.append('functional')
    blastdf.loc[:,'type'] = xs
    return blastdf

def exec_blast(pseudo_fna,subject_fna,odir,force=False):
    """
    the query must be the extracted pseudogene sequences
    the subject must be the whole genome sequences
    """
    name1 = pseudo_fna.split('/')[-1].split('_')[0].rsplit('.',1)[0]
    name2 = subject_fna.split('/')[-1].split('_')[0].rsplit('.',1)[0]
    ofile = f"{odir}/{name1}_{name2}_blastN.tsv"
    cmd = f"blastn -query {pseudo_fna} -subject {subject_fna} -outfmt 6 > {ofile}"
    #print(cmd)
    if not exists(ofile) or force:
        os.system(cmd)
    return ofile

def read_blast(blast_ofile,):
    blast_df = pd.read_csv(blast_ofile,sep='\t',header=None)
    anno_blastdf = anno_pseudo(blast_df,all_pseudo_df)
    prepared_blastdf = anno_blastdf.loc[anno_blastdf['type']=='functional',:]
    return prepared_blastdf

def overlapped_with_IS(contig,start,end):
    start,end = min([start,end]),max([start,end])
    _df = all_IS_df.loc[all_IS_df['contig']==contig,:]
    sub_df = _df.loc[((_df['start'] >= start) & (_df['stop'] <= end)) | \
         ((_df['start'] <= start) & (_df['stop'] >= start)) | \
         ((_df['start'] <= end) & (_df['stop'] >= end))  ,:]
    return sub_df

def get_genomic_parts(contig,start,end):
    
    start,end = min([start,end]),max([start,end])
    _df = all_genomic_component_df.loc[all_genomic_component_df['contig']==contig,:]
    sub_df = _df.loc[((_df['start'] >= start) & (_df['end'] <= end)) | \
         ((_df['start'] <= start) & (_df['end'] >= start)) | \
         ((_df['start'] <= end) & (_df['end'] >= end))  ,:]
    i1 = P.closed(start,end)
    idx = []
    for _idx, row in sub_df.iterrows():
        i2 = P.closed(row['start'],row['end'])
        if i2.overlaps(i1):
            i = i1.intersection(i2)
            cov = (i.upper - i.lower)/(i2.upper-i2.lower)
            idx.append((_idx,cov))
    def get_v(name_cov):
        name,cov = name_cov
        if 'ign' in name:
            return (0,cov)
        elif 'ign' not in name:
            return (1,cov)
    idx = sorted(idx,key=lambda x:get_v(x))
    if len(idx)==0:
        return 'ign'
    return idx[-1][0]


def identified_diff(c1,s1,e1,
                    c2,s2,e2,
                    strand2=1,
                    contig2seq={},
                    match_length=0,
                    extend_bp=0):
    # s must be smaller than e. the strand2 is the match strand between query and subject/
    if type(extend_bp)==int:
        s1,e1 = int(s1)-extend_bp,int(e1)+extend_bp
        s2,e2 = int(s2)-extend_bp,int(e2)+extend_bp
    else:
        cds_extendbp = extend_bp[2],extend_bp[3]
        pseudo_extendbp = extend_bp[0],extend_bp[1]
        s1,e1 = int(s1)-pseudo_extendbp[0],int(e1)+pseudo_extendbp[1]
        # if strand2==1:
        s2,e2 = int(s2)-cds_extendbp[0],int(e2)+cds_extendbp[1]
        # else:
        #     s2,e2 = int(s2)-extend_bp[1],int(e2)+extend_bp[0]+1
    
    seq1 = get_fna(s1,e1,c1,'t1.fna',contig2seq)
    seq2 = get_fna(s2,e2,c2,'t2.fna',contig2seq)
    if strand2 == -1:
        seq2 = seq2.reverse_complement()
    try:
        diff_blocks = DiffBlocks.from_sequences(seq1, seq2)
    except:
        diff_blocks = None

    if diff_blocks is None:
        return 'unknown error'
    feas = diff_blocks.diffs_as_features()
    feas = [_ for _ in feas
            if _.type!='diff_equal' and _.location.start.real!=0 and _.location.end.real!=match_length]
    # equal, starting diff and ending diff are not considered
    frameshifts_labels = [len(_.qualifiers['label'][1:]) for _ in feas if _.type!='diff_replace']
    mutation_labels = [_ for _ in feas if _.type=='diff_replace']
    if len(frameshifts_labels)==0 and len(mutation_labels)!=0:
        return 'mutation drive'
    elif len(frameshifts_labels)!=0 :
        #  any([_%3!=0 for _ in frameshifts_labels]) not necessary....
        return 'frameshift drive'
    return 'others'


def classify_each_pseudo(pseudo_idx,row,contig2seq):
    query_contig,(query_start,query_end) = pseudo_idx.split(':')[0], pseudo_idx.split(':')[1].split('-')
    query_start,query_end = int(query_start),int(query_end) 
    
    pseudo_idx,subject_contig,identity,match_length,_,_,pseudo_start,pseudo_end,subject_start,subject_end = row[:10].values.tolist()
            
    strand = 1 if subject_start - 1<=subject_end else -1
    if strand ==1:
        subject_start,subject_end = subject_start-1,subject_end
    else:
        subject_start,subject_end = subject_end-1,subject_start
    # -1 is transform the 1-based to 0-based
    
    ###### check if the subject contig is overlapped with IS
    subject_IS = overlapped_with_IS(subject_contig,subject_start,subject_end)
    query_IS = overlapped_with_IS(query_contig,query_start,query_end)
    if subject_IS.shape[0]!=0 or query_IS.shape[0]!=0:
        return ('IS-mediated',f"{subject_contig}:{subject_start}-{subject_end}",0,strand)
        # continue
    ########## match and no IS found ###############
    query_part = get_genomic_parts(query_contig,query_start,query_end)
    subject_part = get_genomic_parts(subject_contig,subject_start,subject_end)

    if 'ign' not in query_part and 'ign' in subject_part:
        return ('IS-mediated', f"{subject_contig}:{subject_start}-{subject_end}",0,strand)
        #continue
    elif 'ign' not in subject_part:
        cds_start,cds_end,cds_strand = all_genomic_component_df.loc[subject_part,['start','end','strand',]]
        left_length = subject_start - cds_start
        right_length = cds_end - subject_end
        right_length = 0 if right_length <=0 else right_length

        cds_extendbp = (left_length,right_length)
        if strand == 1:
            pseudo_extendbp = (left_length,right_length)
            # ex_bp = (0,right_length)
        else:
            pseudo_extendbp = (right_length,left_length)
            # ex_bp = (left_length,0)
            # if pseudo_start - left_length<=0 and right_length>=200:
            #     # exception
            #     ex_bp = (0,right_length)
        ex_bp = tuple(list(pseudo_extendbp)+list(cds_extendbp))
        diff_type = identified_diff(query_contig,query_start,query_end,
                                    subject_contig,subject_start,subject_end,
                                    strand2=strand,contig2seq=contig2seq,match_length=match_length,
                                    extend_bp=ex_bp)
        return (diff_type,f"{subject_contig}:{subject_start}-{subject_end}",ex_bp,strand)
                  
    else:
        diff_type = identified_diff(query_contig,query_start,query_end,
                                    subject_contig,subject_start,subject_end,
                                    strand2=strand,contig2seq=contig2seq,match_length=match_length,
                                    extend_bp=0)
        return (diff_type,f"{subject_contig}:{subject_start}-{subject_end}",0,strand)
    

def classify_pseudos(prepared_blastdf,contig2seq,top_n = 2,verbose=False): 
    pseudo2type_list = defaultdict(list)
    if verbose:
        iter_obj = tqdm(prepared_blastdf.groupby(0))
    else:
        iter_obj = prepared_blastdf.groupby(0)
    for pseudo_idx,sub_blast_df in iter_obj:
        for _,row in list(sub_blast_df.iterrows())[:top_n]:
            pseudo2type_list[pseudo_idx].append(classify_each_pseudo(pseudo_idx,row,contig2seq))  
    return pseudo2type_list


def main(pseudo_fna,query_fna,subject_fna,odir,overwrite=False):
    contig2seq = load_contig2seq([query_fna,subject_fna])
    blast_ofile = exec_blast(pseudo_fna,subject_fna,odir,force=overwrite)
    prepared_blast_ofile = blast_ofile.replace('_blastN.tsv',"_preparedblastN.tsv")
    if not exists(prepared_blast_ofile) or overwrite:
        while getsize(blast_ofile)==0:
            blast_ofile = exec_blast(pseudo_fna,subject_fna,odir,force=True)
            tqdm.write(blast_ofile)
            time.sleep(1)
        prepared_blastdf = read_blast(blast_ofile)
        prepared_blastdf.to_csv(prepared_blast_ofile,sep='\t',index=False)
    else:
        prepared_blastdf = pd.read_csv(prepared_blast_ofile,sep='\t',index_col=False)
        prepared_blastdf.columns = [int(_) for _ in prepared_blastdf.columns[:-1]] + [prepared_blastdf.columns[-1]]
    pseudo2type_list = classify_pseudos(prepared_blastdf,contig2seq,top_n=2,verbose=True)
    
    ofile = f"{odir}/{basename(blast_ofile).split('_blastN')[0]}_rawtypes.txt"
    with open(ofile,'w') as f1:
        for pid,v in pseudo2type_list.items():
            for (ft,subject_pos,exbp,strand) in v:
                f1.write(f"{pid}\t{ft}\t{subject_pos}\t{exbp}\t{strand}\n")
    
    pseudo2final_type = {}
    for k,v in pseudo2type_list.items():
        type_set = ';'.join(sorted(set([_[0] for _ in v])))
        if 'IS-mediated' in type_set:        
            pseudo2final_type[k] = 'IS-mediated'
        elif 'frameshift drive' in type_set:
            pseudo2final_type[k] = 'frameshift mutation'
        elif 'mutation drive' in type_set:
            pseudo2final_type[k] = 'nonsense mutation'
        elif 'Close to sequence edge' in type_set:
            pseudo2final_type[k] = 'Close to sequence edge'
        elif 'Remains due to Truncation' in type_set:
            pseudo2final_type[k] = 'Partial remains'
        else:
            pseudo2final_type[k] = 'Unclassified'
    ofile = f"{odir}/{basename(blast_ofile).split('_blastN')[0]}.txt"
    with open(ofile,'w') as f1:
        for pid,final_type in pseudo2final_type.items():
            f1.write(f"{pid}\t{final_type}\n")


############################# main ############################################

import multiprocessing as mp
def run(args):
    pseudo_fna,query_fna,subject_fna = args
    main(pseudo_fna,query_fna,subject_fna,odir,overwrite=True)

if __name__ == "__main__":
    # g2pop = pd.read_csv('/mnt/ivy/thliao/project/coral_ruegeria/Merged_popcogeneT/1837Ruegeria_MCs.tsv',sep='\t',index_col=0)
    # g2pop = g2pop['MC'].to_dict()

    comp2name = {'GNM004006175':"Parasedimentitalea marina", 
             'GNM003443535':"R. sp. AD91A", 
             'GNM000011965':"R. pomeroyi DSS-3", 
             'GNM000014065':"R. sp. TM1040"}

    from glob import glob
    gid2nano_fna = {fna.split("/")[-3]: fna
                    for fna in glob(f"/mnt/ivy/thliao/project/coral_ruegeria/nanopore_processing/canu_o/*/09_prokka/*.fna")
                    if '_' not in fna.split('/')[-1]}
    for c in comp2name:
        gid2nano_fna[c] = f"/mnt/ivy/thliao/project/coral_ruegeria/data_processing/pub_dataset/prokka_o/{c}/{c}.fna"
        
    gid2pseudo = {f.split('/')[-1].split('_')[0]:f 
                  for f in glob('/mnt/ivy/thliao/project/coral_ruegeria/nanopore_processing/annotations/pseudofinder/pairwise/*_pseudos_extracted.fasta')}
    
    odir = "/mnt/ivy/thliao/project/coral_ruegeria/nanopore_processing/annotations/pseudofinder/pairwise/"

    # sids = ['GNM526588003432', 'GNM526588003444', 'GNM526588001514', 'GNM526588001567']
    import itertools
    cmds = []
    for gid,gid2 in itertools.product(gid2nano_fna.keys(),gid2nano_fna.keys()):
        pseudo_fna = gid2pseudo[gid]
        subject_fna = gid2nano_fna[gid2]
        # if gid in sids or gid2 in sids:
        #     cmds.append((pseudo_fna,gid2nano_fna[gid],subject_fna))
        #if g2pop[gid].startswith('MC59') and g2pop[gid2].startswith('MC59'):
        if not exists(f"{odir}/{gid}_{gid2}_preparedblastN.tsv"):
            cmds.append((pseudo_fna,gid2nano_fna[gid],subject_fna))
    tqdm.write(f"{len(cmds)}")
    with mp.Pool(processes=20) as tp:
        r = list(tqdm(tp.imap(run, cmds), total=len(cmds)))




# from tqdm import tqdm
# import multiprocessing as mp
# from subprocess import check_call
# def run(cmd):
#     check_call(cmd,shell=1)        
# with mp.Pool(processes=5) as tp:
#     r = list(tqdm(tp.imap(run, cmds), total=len(cmds)))      