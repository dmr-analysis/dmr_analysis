import pandas as pd
import os
#exec(open('map_dmr2genome.py').read())

def count_dmrs_not_mapped2genome(in_sorted_dmr_file,record_out_files,dmr_min_cutoff, is_greater=True):
  #count how many DMRs are not mapped to annotated geneomic regions
  #coumt MR or DMR in genomic files
  all_data_df=[]
  for fil in record_out_files:
     tmp_data_df=pd.read_csv(fil,header=None, sep='\t')
     all_data_df.append(tmp_data_df.copy())

  lines=0
  for i in all_data_df:
     lines+= len(i)

  all_indata_df=pd.concat(all_data_df)
  #if all_indata_df.shape[1]==9:
  all_indata_df.columns=['chrs','start_pos','end_pos','genome_info','mr_chrs','mr_start_pos','mr_end_pos','mr_sites','mr_logReg_proba']
  #elif all_indata_df.shape[1]=8:
  #   all_indata_df.columns=['chrs','start_pos','end_pos','genome_info','mr_chrs','mr_start_pos','mr_end_pos','mr_sites']

  uq_indata_df=all_indata_df.drop_duplicates(subset=['mr_sites'])
  total_uq_mrs=uq_indata_df
  min_cutoff=dmr_min_cutoff
  
  if is_greater:
    print('DMR selection based on >= '+ str(min_cutoff))
    total_uq_dmrs=uq_indata_df[uq_indata_df['mr_logReg_proba']>=min_cutoff]
  else:
    print('DMR selection based on <= ' + str(min_cutoff))
    total_uq_dmrs=uq_indata_df[uq_indata_df['mr_logReg_proba']<=min_cutoff]

  in_dmrs_df=pd.read_csv(in_sorted_dmr_file,header=None,sep='\t')
  in_dmrs_df.columns=['mr_chrs','mr_start_pos','mr_end_pos','mr_sites','mr_logReg_proba']
  total_in_mrs=in_dmrs_df
  if is_greater:
    total_in_dmrs=in_dmrs_df[in_dmrs_df['mr_logReg_proba']>=min_cutoff]
  else:
    total_in_dmrs=in_dmrs_df[in_dmrs_df['mr_logReg_proba']<=min_cutoff]

  print('Number of MR or DMR do not find mapped genome information')
  print(total_uq_mrs.shape[0]-total_in_mrs.shape[0])
  print(total_uq_dmrs.shape[0]-total_in_dmrs.shape[0])

  print('Perentage of MR or DMR mapped to genome info')
  print(total_uq_mrs.shape[0]/total_in_mrs.shape[0])
  print(total_uq_dmrs.shape[0]/total_in_dmrs.shape[0])

  diff_in_mr=set(total_in_mrs.mr_sites.to_list())- set(total_uq_mrs.mr_sites.to_list()) 
  diff_in_dmr=set(total_in_dmrs.mr_sites.to_list())- set(total_uq_dmrs.mr_sites.to_list())

  #check the distribution of unmapped MRs in chromes
  chrs=[]
  for i in range(1,25):
    if i<23:
      chrs.append('chr'+str(i))
    elif i==23:
      chrs.append('chrX')
    elif i==24:
      chrs.append('chrY')

  diff_in_dmr_df=pd.DataFrame(data=list(diff_in_dmr),columns=['dmr_sites'])
  dict_chr={}
  total=0
  for ii in chrs:
    dict_chr[ii]=diff_in_dmr_df[diff_in_dmr_df.dmr_sites.str.contains(ii+':')]
    total += dict_chr[ii].shape[0]
  return total, dict_chr

if __name__=='__main__':
  dmr_min_cutoff=0.8
  top_percent=0.91
  isSmooth=0
  in_sorted_dmr_file='out_raw_dmr/24_chroms_high_miniPercentChange_gt_0.0001_Pcutoff_0.05_isSmooth_'+ str(isSmooth) +'_isModTest_0__all_dmrRanking_top_' +str(top_percent) +'_minLogReg_proba_'+ str(dmr_min_cutoff)+'.bed'
  in_region_files='out/DMR/list_region_files.txt'

  region_files=pd.read_csv(in_region_files,header=None)
  region_files=region_files.loc[:,0].to_list()

  methylation_file=in_sorted_dmr_file
  reference_file='out/DMR/data/hg19.refFlat_clean_sorted.bed'
  min_overlap=1E-9
  out_folder='out_raw_dmr/'

  #for each region to find its DMR or MRs
  record_out_files=[]
  for fil in region_files:
    region_file=fil
    region_name=os.path.basename(fil).split('_')[0].lower()
    print(region_name)
    out_methylation_name = out_folder + '/' + os.path.basename(methylation_file)[:-4]
    out_region_name = os.path.basename(region_file)[:-4]

    #print(methylation_file)
    #print(region_file)
    dist5_methylation_file = out_methylation_name + '_' + 'noGenes.bed'
    print(dist5_methylation_file)

    #print(out_methylation_name,out_region_name)

    if region_name == '5dist':
        # For 5distance we first remove genes(TSS, geneBody and TES) from the two methylation files
        os.system('bedtools intersect -a ' + methylation_file + ' -b ' + reference_file + \
                  ' -v > ' + dist5_methylation_file)

        out = dist5_methylation_file[:-4] + '_' + out_region_name + '_overlap' + str(min_overlap) + '.bed'

        os.system('bedtools intersect -a ' + region_file + ' -b ' + \
                  dist5_methylation_file + ' -wa -wb -f ' + str(min_overlap) + ' > ' + out)

        os.system('rm ' + dist5_methylation_file)  # removes temporary file
    else:
        out = out_methylation_name + '_' + out_region_name + '_overlap' + str(min_overlap) + '.bed'
        os.system('bedtools intersect -a ' + region_file + ' -b ' + methylation_file + \
                  ' -wa -wb -f ' + str(min_overlap) + ' > ' + out)
    print(out)
    record_out_files.append(out)


  #count how many DMRs are not mapped to annotated geneomic regions
  #coumt MR or DMR in genomic files
  #dmr_min_cutoff=0.8
  total, dict_chr=count_dmrs_not_mapped2genome(in_sorted_dmr_file,record_out_files,dmr_min_cutoff)

















