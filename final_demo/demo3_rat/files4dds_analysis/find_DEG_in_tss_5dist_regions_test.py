#this script is used to find DEG genes in TSS and 5distnace regions, then export all of them with a bed format such
#as chr, start, end, mr_ID, Pvalues et al.
import pandas as pd
import numpy as np
import os

#find MRs mapped to TSS and 5Distance region files 
in_folder='../../../final_demo_data/rat_data/out_data/DMR_CpG_context/out_map2genome/'
tss_file='5_chroms_all_mr_data_range_dmrRanking_TSS_Up5000_Down1000_removedShort_overlap1e-09.bed'
in_tss_file=os.path.join( in_folder , tss_file)
print('Read: ')
print(in_tss_file)

out_result_folder='../../../final_demo_data/rat_data/out_data/DMR_CpG_context'

dist_file='5_chroms_all_mr_data_range_dmrRanking_noGenes_5dist_Up1000000_Up5000removedShort_overlap1e-09.bed'
in_5dist_file=os.path.join( in_folder , dist_file)
print(in_5dist_file)

#DEG file
in_deg_file='../../../final_demo_data/rat_data/in_data/DEG/Adrenal1vsAdrenal2_DEG_genes_zscores.tsv'
print(in_deg_file)

#read all data
tss_df=pd.read_csv(in_tss_file,sep='\t',header=None)
dist_df=pd.read_csv(in_5dist_file,sep='\t',header=None)
deg_df=pd.read_csv(in_deg_file,sep='\t')

#combine all data to datagrame
combined_df=pd.concat([tss_df,dist_df]).copy()
combined_df.reset_index(inplace=True,drop=True)
combined_df.columns=['chrom','start_pos','end_pos','gene_type','mr_chrom','mr_start_pos','mr_end_pos','mr_id','pval']

#find DMR and DEG regions are overlapping
dmr_combined_df=combined_df[combined_df.mr_id.apply(lambda x: ':D' in x)].copy()
dmr_combined_df.reset_index(inplace=True, drop=True)

deg_genes=set(deg_df.gene_name.to_list())
dmr_combined_df['gene_name']=dmr_combined_df.gene_type.apply(lambda x: x.split('||')[2].split(':')[0] )

deg_dmr_combined_df=dmr_combined_df[dmr_combined_df.gene_name.apply(lambda x : x in deg_genes)].copy()
deg_dmr_combined_df.reset_index(inplace=True,drop=True)

#export DEG, DMR TSS and 5dist genes
out_df=deg_dmr_combined_df[['mr_chrom','mr_start_pos','mr_end_pos','mr_id','pval','gene_type','gene_name']]
out_file=os.path.join(out_result_folder, 'dmr_regions_in_deg_tss_5dist_rat.bed')
out_df.to_csv(out_file,sep='\t',header=None,index=False)
print('Export:')
print('\t ', out_file)


