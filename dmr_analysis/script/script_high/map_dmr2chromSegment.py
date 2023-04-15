#map DMR to predicted chromatin segment from 6 cell lines with seven states
#https://genome.ucsc.edu/cgi-bin/hgFileUi?db=hg19&g=wgEncodeAwgSegmentation
#TSS     Bright Red      Predicted promoter region including TSS
#PF      Light Red       Predicted promoter flanking region
#E       Orange  Predicted enhancer
#WE      Yellow  Predicted weak enhancer or open chromatin cis regulatory element
#CTCF    Blue    CTCF enriched element
#T       Dark Green      Predicted transcribed region
#R       Gray    Predicted Repressed or Low Activity region
#exec(open('map_dmr2chromSegment.py').read())
import os
import glob
import pandas as pd
#import combine_multChrs_dmr4ranking as combin_dmr_rank
from .map_dmr2genome import count_dmrs_not_mapped2genome
from .dmr_utility import chrStr2numeric



def read_bed_file_and_export2type_bed_file(in_files, column_name, out_folder,min_length):
  #input list of bed file in_files , export type to different bed file in out_folder with min_length
  record_out_file_names=[]
  for fil in in_files:
    print(fil)
    fil_name=fil.split('.')[0]
    tmp_data_df=pd.read_csv(fil,sep='\t',header=None,compression='gzip')
    tmp_data_df.columns=column_name
    #record_data_df[fil_name]=tmp_data_df.copy()
    print(tmp_data_df.shape)

    #export each type of region to separate files
    uq_type_in_df=tmp_data_df['type'].drop_duplicates().to_list()
    for ti in uq_type_in_df:
       sub_tmp_data_df=tmp_data_df.copy()
       sub_tmp_data_df=sub_tmp_data_df[ (sub_tmp_data_df['type']==ti) ]
       sub_tmp_data_df=sub_tmp_data_df[['chrs','start_pos','end_pos','type']]
       out_file=os.path.join(out_folder,os.path.basename(fil_name)+ '_min'+str(min_length) + '_'+ti+'.bed.gz')
       print(out_file)

       #filter regions< min_length
       sub_tmp_data_df['len_of_region']=sub_tmp_data_df['end_pos']-sub_tmp_data_df['start_pos']
       print(sub_tmp_data_df.shape)
       sub_tmp_data_df=sub_tmp_data_df[sub_tmp_data_df['len_of_region']>=min_length]
       print(sub_tmp_data_df.shape)
       sub_tmp_data_df.to_csv(out_file,sep='\t',header=False,index=False,compression='gzip' )
       record_out_file_names.append(out_file)
  return record_out_file_names, uq_type_in_df

def sort_bed_file_df(combined_data_df,column_names,ishuman):
  #sort bed file by chromsome and start position and return sorted bed file df
  #here assume the first and second column names are chrs, start_pos
  #sort comnbined bed file
  combined_data_df.columns=column_names
  combined_data_df['order_by_chr']=combined_data_df['chrs'].apply(chrStr2numeric, human=ishuman)
  combined_data_df['order_by_chr']=combined_data_df['order_by_chr'].astype(int)
  sorted_combined_data_df=combined_data_df.copy()
  sorted_combined_data_df=sorted_combined_data_df.sort_values(['order_by_chr','start_pos'],ascending=True)
  sorted_combined_data_df=sorted_combined_data_df.reset_index(drop=True)
  sorted_combined_data_df=sorted_combined_data_df.drop(columns=['order_by_chr'])
  return sorted_combined_data_df

if __name__=='__main__':
  #generate chromatin segmant files:   
  in_folder='in_data/chromatin_segment/in/'
  in_file_name='*.bed.gz'
  in_files=glob.glob(os.path.join(in_folder,in_file_name))
  column_name=['chrs','start_pos','end_pos','type','len','strand','starts','ends','code']

  out_folder='out/DMR/data/chromatin_segment/'
  #read all files in a dictionary as dataframe
  #set minimum lenght of each region
  min_length=10
  out_file_names,type_names =read_bed_file_and_export2type_bed_file(in_files,column_name, out_folder,min_length)

  #combine the same type of bed file from different cell lines then sort the combined file to export in bed format
  out_file_name_df=pd.DataFrame(columns=['name'],data=out_file_names)
  sorted_out_file_name=[]
  for ty in type_names:
    out_file='combined_six_cells_chromatin_segment_min'+str(min_length)+'_'+ty + '.bed'

    #combine bed files
    combined_data_df=[]
    extracted_files=out_file_name_df['name'][out_file_name_df['name'].str.contains('_'+ty+'.bed')].to_list() 
    tmp_df=[]
    for fil in extracted_files:
       tmp_data_df=pd.read_csv(fil, sep='\t',header=None, compression='gzip')
       tmp_df.append(tmp_data_df)
    combined_data_df=pd.concat(tmp_df)
    columns_name=['chrs','start_pos','end_pos','type','len_of_region'] 
    ishuman=True
    sorted_combined_data_df=sort_bed_file_df(combined_data_df,columns_name,ishuman)

    out_file=os.path.join(out_folder,out_file)
    print(out_file)
    sorted_combined_data_df.to_csv(out_file,sep='\t',index=False, header=False)
    sorted_out_file_name.append(out_file)

  #remove exported type bed files after combining them
  for fil in out_file_names:
    os.remove(fil)

  #use bed tools to merger regions in each type
  merged_out_files=[]
  for fil in sorted_out_file_name:
     out_fil=fil.replace('.bed','_merged.bed')
     command='bedtools merge -i ' + fil +  ' -c 4 -o distinct >' + out_fil
     print(out_fil)
     os.system(command)
     os.remove(fil)
     merged_out_files.append(out_fil)


  #map DMR to chromatin segmant files
  in_dmr_folder='out_raw_dmr/'
  dmr_min_cutoff=0.8
  top_percent=0.91
  isSmooth=0
  in_sorted_dmr_file=in_dmr_folder +'24_chroms_high_miniPercentChange_gt_0.0001_Pcutoff_0.05_isSmooth_' + str(isSmooth) \
                       +'_isModTest_0__all_dmrRanking_top_'+ str(top_percent)+'_minLogReg_proba_'+str(dmr_min_cutoff) +'.bed'
  methylation_file=in_sorted_dmr_file
  out_folder=in_dmr_folder #'out/DMR/data/chromatin_segment/'
  min_overlap=1E-9

  record_out_files=[]
  for fil in merged_out_files:
    region_file=fil
    out_methylation_name = out_folder + '/' + os.path.basename(methylation_file)[:-4]
    out_region_name = os.path.basename(fil)[:-4]
    out = out_methylation_name + '_' + out_region_name + '_overlap' + str(min_overlap) + '.bed'
    os.system('bedtools intersect -a ' + region_file + ' -b ' + methylation_file + \
                  ' -wa -wb -f ' + str(min_overlap) + ' > ' + out)
    print(out)
    record_out_files.append(out)

  #count how many DMRs are not mapped to annotated geneomic regions
  #coumt MR or DMR in genomic files
  #coumt MR or DMR in genomic files
  #dmr_min_cutoff=0.8
  total, dict_chr=map_dmr2genome.count_dmrs_not_mapped2genome(in_sorted_dmr_file,record_out_files,dmr_min_cutoff)








