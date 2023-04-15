#here we assign both genome annotatoin and  chromatinSegment annotation to the DMR or MRs.
import os
import pandas as pd
import glob
import multiprocessing as mp
import numpy as np

#exec(open("combine_annota2dmrs.py").read())
def merge_multiple_files2df(record_out_files):
  all_data_df=[]
  for fil in record_out_files:
     tmp_data_df=pd.read_csv(fil,header=None, sep='\t')
     all_data_df.append(tmp_data_df.copy())

  lines=0
  for i in all_data_df:
     lines+= len(i)
     #print(i)
  all_indata_df=pd.concat(all_data_df)
  all_indata_df=all_indata_df.reset_index(drop=True)
  return lines, all_indata_df


#def find_genome_chromSegment4dmr(selected_dmr_sorted_df0, block_chunks_index,\
#                   all_in_genome_df,all_in_chromSegment_df,num_of_processes):
def find_genome_chromSegment4dmr(args):
#def find_genome_chromSegment4dmr(args):
#  selected_dmr_sorted_df0,\
#  block_chunks_index,all_in_genome_df,all_in_chromSegment_df,\
#  num_of_processes=args

  global gb_selected_dmr_sorted_df0 , gb_all_in_genome_df, gb_all_in_chromSegment_df 
  block_chunks_index,num_of_processes=args

  #find matached genome and crhomSegmant annotationfor each MR or DMR in df
  selected_dmr_sorted_df=gb_selected_dmr_sorted_df0.copy()
  selected_dmr_sorted_df=selected_dmr_sorted_df.loc[block_chunks_index,:]
  record_index=[None]*selected_dmr_sorted_df.shape[0]
  record_genome_info=[None]*selected_dmr_sorted_df.shape[0]
  record_chromSegment_info=[None]*selected_dmr_sorted_df.shape[0]
  loop=0
  for index, row in selected_dmr_sorted_df.iterrows():
    #print(index)
    tmp_row_mr_sites=row.mr_info
    tmp_rows=gb_all_in_genome_df[gb_all_in_genome_df.mr_sites==tmp_row_mr_sites].copy()
    if tmp_rows.genome_info.shape[0]==0:
      tmp_genome_row='NA'
    elif tmp_rows.genome_info.shape[0]>1 :
      tmp_genome_row='~'.join(tmp_rows.genome_info.to_list())
    elif tmp_rows.genome_info.shape[0]==1:
      tmp_genome_row=tmp_rows.genome_info.values[0]
    tmp_crows=gb_all_in_chromSegment_df[gb_all_in_chromSegment_df.mr_sites==tmp_row_mr_sites].copy()
  #  tmp_chromSegment_row='$'.join(tmp_crows[['chrs','start_pos','end_pos','genome_info']].apply(lambda x: ':'.join(x.astype(str)),axis=1))
    if tmp_crows.shape[0]==0:
      tmp_chromSegment_row='NA'
    elif tmp_crows.genome_info.unique().shape[0]>1:
      tmp_chromSegment_row='~'.join(list(tmp_crows.genome_info.unique()))
    elif tmp_crows.genome_info.unique().shape[0]==1:
      tmp_chromSegment_row=tmp_crows.genome_info.unique()[0]
    record_index[loop]=index
    record_genome_info[loop]=tmp_genome_row
    record_chromSegment_info[loop]=tmp_chromSegment_row
    loop +=1
    #print(loop)
  #insert to df
  selected_dmr_sorted_df2=selected_dmr_sorted_df.copy()
  selected_dmr_sorted_df2.loc[record_index,'genome_info']=record_genome_info
  selected_dmr_sorted_df2.loc[record_index,'chromSegment_info']=record_chromSegment_info
  if num_of_processes is not None:
     out_file='selected_dmr_sorted_df2_' +str(num_of_processes) +'.csv'
     print(out_file)
     selected_dmr_sorted_df2.to_csv(out_file,index=False, sep='\t')
     return out_file
  else:
     return selected_dmr_sorted_df2

def find_parallel_process_and_index(num_of_processes,selected_dmr_sorted_df0):
  #find index for ecah chunk block of the data and the matched number of processe
  all_index=selected_dmr_sorted_df0.index
  num_in_chunks=int(all_index.stop/num_of_processes) 
  all_blocks=np.linspace(all_index.start,all_index.stop-1,num=all_index.stop)
  all_blocks=np.array(all_blocks,int)
  block_chunks = [all_blocks[x:x+num_in_chunks] for x in range(0, len(all_blocks), num_in_chunks)]
  if len(block_chunks) > num_of_processes:
      num_of_processes=len(block_chunks)
  return num_of_processes, block_chunks

def initpool(t_selected_dmr_sorted_df0, t_all_in_genome_df,t_all_in_chromSegment_df):
    global gb_selected_dmr_sorted_df0, gb_all_in_genome_df, gb_all_in_chromSegment_df
    gb_selected_dmr_sorted_df0=t_selected_dmr_sorted_df0
    gb_all_in_genome_df=t_all_in_genome_df
    gb_all_in_chromSegment_df=t_all_in_chromSegment_df

def call_parallel_genome_function(num_of_processes,selected_dmr_sorted_df0,block_chunks,all_in_genome_df,all_in_chromSegment_df,in_dmr_sorted_files,logReg_proba_cutoff):
  #useing parallel computation for calling genome info function
  pool=mp.Pool(processes=num_of_processes,initializer=initpool, initargs=(selected_dmr_sorted_df0,all_in_genome_df,all_in_chromSegment_df) )
  
  #there is a memory issue for very large data set, the pool will crash if there is a block_chunk with too large size!!
  #only solution is to use global variable for df for reducing memory usage in each process 
  #files= pool.map(find_genome_chromSegment4dmr,[(selected_dmr_sorted_df0, block_chunks[loop], all_in_genome_df, all_in_chromSegment_df,loop) for loop in range(0,num_of_processes)],1)
  #files= pool.imap(find_genome_chromSegment4dmr,[(selected_dmr_sorted_df0, block_chunks[loop], all_in_genome_df, all_in_chromSegment_df,loop) for loop in range(0,num_of_processes)],2)
  #files= pool.starmap(find_genome_chromSegment4dmr,[(selected_dmr_sorted_df0, block_chunks[loop], all_in_genome_df, all_in_chromSegment_df,loop) for loop in range(0,num_of_processes)])

  files= pool.map(find_genome_chromSegment4dmr,[(block_chunks[loop],loop) for loop in range(0,num_of_processes)],1)
  pool.close()

  files=list(filter(None,files))
  record_pd=[]
  for fi in files:
     record_pd.append(pd.read_csv(fi,sep='\t'))
 
  all_pd=pd.concat(record_pd)
  
  out_file=in_dmr_sorted_files.replace('.bed','_genome_chromSegment_info_ProbaGt_'+ str(logReg_proba_cutoff)+'.bed')
  print("Export all results at : " + out_file)
  all_pd.to_csv(out_file,sep='\t',index=False, float_format='%g')
  all_pd=all_pd.reset_index(drop=True)
  #remove the rest of files
  for fi in files:
      os.remove(fi) 
  return all_pd



if __name__=='__main__':

  dmr_mini_cutoff=0.7
  top_percent=0.93
  isSmooth=0
  in_dmr_sorted_files='24_chroms_high_miniPercentChange_gt_0.0001_Pcutoff_0.05_isSmooth_' + str(isSmooth)+ \
                       '_isModTest_0__all_dmrRanking_top_'+str(top_percent) +'_minLogReg_proba_' +str(dmr_mini_cutoff) +'.bed'
  in_out_genome_annot_folder='out_raw_dmr/'
  in_dmr_sorted_df=pd.read_csv(os.path.join(in_out_genome_annot_folder,in_dmr_sorted_files),sep='\t',header=None)
  in_dmr_sorted_df.columns= ['mr_chrs','mr_start_pos','mr_end_pos','mr_info','mr_logReg_proba']
  

  logReg_proba_min_cutoff=dmr_mini_cutoff
  intergenetic_min_len=10
  in_out_genome_file='24*proba_'+str(logReg_proba_min_cutoff) +'*.bed'
  tmp_in_genome_files=glob.glob(os.path.join(in_out_genome_annot_folder+'genomes/',in_out_genome_file))

  #load all genome annotated files
  for fi in tmp_in_genome_files:
     if 'intergenic' in fi:
       if not ('minLen'+str(intergenetic_min_len) in fi):
          tmp_in_genome_files.remove(fi)

  print(tmp_in_genome_files)
  #read genome info files
  in_genome_lines,all_in_genome_df=merge_multiple_files2df(tmp_in_genome_files)
  all_in_genome_df.columns=['chrs','start_pos','end_pos','genome_info','mr_chrs','mr_start_pos','mr_end_pos','mr_sites','mr_logReg_proba']


  #load all chromSegment annotated files
  min_merge_len=10
  in_out_chromSegment_folder=in_out_genome_annot_folder+'chromSegment/'
  in_out_chromSegment_file='24*proba_'+str(logReg_proba_min_cutoff)+ '*_min' +str(min_merge_len) + '*.bed'
  tmp_in_chromSegment_files=glob.glob(os.path.join(in_out_chromSegment_folder,in_out_chromSegment_file))
  #read chromSegment info files
  all_in_chromSegment_lines, all_in_chromSegment_df=merge_multiple_files2df(tmp_in_chromSegment_files) 
  all_in_chromSegment_df.columns=['chrs','start_pos','end_pos','genome_info','mr_chrs','mr_start_pos','mr_end_pos','mr_sites','mr_logReg_proba']

  #do parallel in genome annotation
  num_of_processes=15
  #selected_dmr_sorted_df0= in_dmr_sorted_df[in_dmr_sorted_df.mr_logReg_proba>=logReg_proba_min_cutoff].copy()
  selected_dmr_sorted_df0=in_dmr_sorted_df.copy()
  selected_dmr_sorted_df0=selected_dmr_sorted_df0.reset_index(drop=True)
  #selected_dmr_sorted_df0=selected_dmr_sorted_df0.loc[0:1000,:]
  num_of_processes,block_chunks=find_parallel_process_and_index(num_of_processes,selected_dmr_sorted_df0)
  print(num_of_processes)
  
  #input('Click to start run ')
  all_out_df=call_parallel_genome_function(num_of_processes,selected_dmr_sorted_df0,block_chunks,all_in_genome_df,all_in_chromSegment_df,in_dmr_sorted_files,logReg_proba_min_cutoff)  

