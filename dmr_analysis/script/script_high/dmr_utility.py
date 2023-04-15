#Utility functions in DMR pipeline
import numpy as np
import glob
import os
import pandas as pd
import logging
import json
import argparse
import subprocess
import math
import time

#load in-house library
from .others.convert2hmst_seq_input_bed import read_bed_gz_files 
#import plot_test
from .others.find_regions_test import find_methylation_blocks

#np.seterr(all='raise')
np.seterr(all='ignore')

def divide_data2chunk(num_of_processes, list_of_ids):
 '''divide all id to blocks based on num of processes
   num_of_processes=15
   all_blocks= gene_df.gene_id.unique()
 '''
 all_blocks=list_of_ids
 num_in_chunks=int(math.ceil(len(all_blocks)/num_of_processes))
 block_chunks= [all_blocks[x:x+num_in_chunks] for x in range(0, len(all_blocks), num_in_chunks)]

 #do find_gene_bins_position in each chunk of blocks
 if len(block_chunks) > num_of_processes:
      num_of_processes=len(block_chunks)
 elif len(block_chunks) < num_of_processes:
      print('Number of blocks smaller than the number of available processes ')
      num_of_processes=len(block_chunks)
 elif len(block_chunks)==1:
      num_of_processes= 1
 return block_chunks, num_of_processes

def submit_job_in_subprocess(record_cmds):
   '''submit command line job in multiple tasks'''
#  processes=[]
#  for cmd in record_cmds:
#        processes.append(subprocess.Popen(cmd,shell=True))
#  for p in processes:
#        p.communicate()
   processes=[]
   num_processes=50
   cmd_processed=0
   print('Running ', len(record_cmds), ' jobs in ', min(num_processes, len(record_cmds)), ' parallel processes!' )
   while len(record_cmds) >0:
       while len(processes) >= num_processes :
            for i in range(len(processes)):
               if processes[i].poll() is not None:
                  del processes[i]
                  break

            if len(processes) >= num_processes:
               time.sleep(0.01)
       
       processes.append(subprocess.Popen(record_cmds[0],shell=True))
       cmd_processed +=1 
       del record_cmds[0]

   for p in processes:
       p.communicate()
 
def prepare_result_folder(res_folder):
    if os.path.exists(res_folder):
       print("The result folder, ", res_folder, " , exits ")
    else:
       print("Create result folder, ", res_folder)
       os.mkdir(res_folder)

    #if os.path.exists(res_folder): 
    #    print("The result folder,", os.path.abspath(res_folder), ", exists and will be erased")
    #
    #clear_folder(res_folder)

def clear_folder(folder):
    if os.path.exists(folder):
        shutil.rmtree(folder)

    os.mkdir(folder)


def chrStr2numeric(string, human):
   #convert string chr to numeric value for sorting
   string=string.lower()
   string=string.replace('chr','')
   if human :
     if string=='x':
        string='23'
     elif string=='y':
        string='24'
     elif string=='m':
        string='25'
   else:
     if string =='x':
        string='20'
     elif string=='y':
        string='21'
     elif string=='m':
        string='22'
   return string


def combine_multiple_pd(in_files, record_raw_df, needed_columns):
   ''' 
     Input: in_files is a list of input file names
            record_raw_df is a dictionary of dataframe where keys are input file names
            needed_columns is the column name of dataframes that need be combined from multiple file/dataframes
     Output: combined datafram
     in_files is input column names/key of dictionary record_raw_df,
     record_raw_df is a dictionary of dataframe, 
     needed_columns is the exporting column key words from dataframe
     here combine multipe dataframe based on chrs, starts and ends columns, and only export the selected columns in needed_columns
     Export a data frame, where the first three columns are Chrs, Starts,Ends, then the rest of columns are selected columns from multiple input dataframes
   '''
   left_pd=[]
   right_pd=[]
   tmp_pd=[]
   for i in range(len(in_files)-1):
     if i ==0:
       left_pd=record_raw_df[in_files[i]]
       left_name=os.path.basename(in_files[i]).replace('.bed.gz','')
       right_pd=record_raw_df[in_files[i+1]]
       right_name=os.path.basename(in_files[i+1]).replace('.bed.gz','')
     else:
       left_pd=tmp_pd
       left_name=os.path.basename(in_files[i]).replace('.bed.gz','')
       right_pd=record_raw_df[in_files[i+1]]
       right_name=os.path.basename(in_files[i+1]).replace('.bed.gz','')
     #print left_name, right_name
     #raw_input(' ')
     tmp_pd=pd.merge(left_pd,right_pd,how='outer',on=['Chrs','Starts','Ends'],suffixes=('_'+left_name,'_'+right_name))
     needed_names=[]
     #needed_columns='Methylation'
     for cl in tmp_pd.columns:
       if needed_columns in cl:
          needed_names.append(cl)
     tmp_pd2=tmp_pd[['Chrs','Starts','Ends']+needed_names]
     #fill nan as 0
     tmp_pd2=tmp_pd2.fillna(0)
   return tmp_pd2

def sort_and_export_pd(in_chrm,input_df,data_start_pos):
  '''
    Input: input_df is a merged dataframe
           data_start_pos is an index of column that shall be sorted before exporting, this index start from 1
                          here we assume that before this index is chromosome position but after the index are all data columns that will be exported latter.
    Output: export sorted dataframe in a temp file 'tmp_out_pd.csv', and return sorted dataframe, numpy matrix of position and all data  
    Sort position and export the merged data frame
  '''
  df_columns=input_df.columns
  start_pos=data_start_pos
  num_of_rows_in_df=input_df.shape[0]

  #BUG in input_df 
  #print(input_df.Ends.astype('string')) 
  #print(start_pos,df_columns[start_pos-1], input_df[df_columns[start_pos-1]].astype(str) )

  tmp_pd2= input_df.sort_values(by=[df_columns[start_pos-1]],ascending=True)
  positions= tmp_pd2[df_columns[start_pos-1]].values
  methylations=tmp_pd2.iloc[:, start_pos: ].values

  out_tmp_pd2_file=in_chrm+ '_tmp_out_pd.csv'
  tmp_pd2.to_csv(out_tmp_pd2_file,index=False, sep='\t',float_format='%g')
  print('Export ',num_of_rows_in_df, ' lines of data at ' + out_tmp_pd2_file)
  return tmp_pd2,positions, methylations, out_tmp_pd2_file


def make_new_bin_range(bin_range):
   #here assume input bin_range is already sorted and the last one is the true maximum
   #the new_bin_range will remove any bin range smaller than the true maxiumum 

   max_range=bin_range[-1]
   #remove duplicated ranges 
   bin_range=list(set(bin_range))
   sorted_bin_range=np.sort(bin_range)
   index2range=int(np.where(sorted_bin_range==max_range)[0])+1
   new_bin_range=sorted_bin_range[:index2range]
   return new_bin_range


def combine_multiple_files_to_one_dataframe(in_folder,in_chrm, group_key, needed_columns,logger_num) :
   #first combine multiple files of the same chromosome data into one dataframe
   #in_folder='in_data/WGBS-data/'
   #in_chrm='chr22'
   #group_key='chr22'
   #needed_columns= 'Methylation'
   ''' combine multiple files into one dataframe
      all input files assue in the same bed format
      Export dataframe sorted by position
   '''

   in_files=[]
   in_files=glob.glob(os.path.join(in_folder,in_chrm, '*'+group_key+'*'))
   #maxium number of rows to read
   num_of_data_points=None

   record_raw_df={}
   for ff in in_files:
     logger_num.info( ff)
     tmp_df= read_bed_gz_files(ff)
     record_raw_df[ff]=tmp_df[:num_of_data_points]

   #combine multiple dataframe based on needed_columns
   input_df = combine_multiple_pd(in_files, record_raw_df, needed_columns)

   #chr_position=1
   data_start_pos=3
   data_end_pos=15
   #num_of_samples=12
   #sort and export dataframe
   sorted_tmp_pd2, positions, methylations, out_tmp_pd2_file = sort_and_export_pd(in_chrm,input_df,data_start_pos)
   return sorted_tmp_pd2, positions, methylations, out_tmp_pd2_file


def find_blocks_in_chrom(in_chrm,num_of_data_points, data_start_pos, data_end_pos, max_read_length, mini_block_size, out_mr_file_name,out_mr_file_path,logger_num) :
   #Then try to find all possible blocks in the chromome based on available data points inthe chromoeome
   #num_of_data_points=None
   #data_start_pos=3
   #data_end_pos=15
   ''' find all possinle MR blocks in a chromosomes
      Input: num_of_data_points is the number of rows in the dataframe will be considered, if it is None then all rows are considered.
             data_start_pos-1 is the chrom position, data_start_pos and data_end_pos are the start and end position of data columns that will be used latter
             max_read_length is the maximum distance between two adjacent MRs/blocks, mini_block_size is the minimum number of data points within a MR/block.
      Output: all block, block size, block_methylation (data matrix), block_lenght, block_methylation_columns (names)
   '''
   tmp_out_pd_file=in_chrm+ '_tmp_out_pd.csv'
   tmp_pd=pd.read_csv(tmp_out_pd_file,sep='\t')
   tmp_pd2=tmp_pd[0:num_of_data_points]
   tmp_position=tmp_pd2.iloc[:,data_start_pos-1].values
   tmp_methylation=tmp_pd2.iloc[:,data_start_pos:data_end_pos].values
   tmp_column_names=tmp_pd2.columns[data_start_pos:].values

   #1000 sees have the best MR distribution
   #max_read_length=200
   #mini_block_size=5
   logger_num.info( 'Blocks with distance greater than %g ', max_read_length)
   logger_num.info( ' and minimum data points in block %g', mini_block_size)
   block, block_size, block_methylation, block_length = find_methylation_blocks(tmp_position,tmp_methylation, tmp_column_names, max_read_length,mini_block_size)
   logger_num.info('block size %g ',  len(block_size))

   block_methylation_columns=tmp_column_names

   #export data
   outf1=out_mr_file_name.replace('.csv','_data.txt')
   logger_num.info(  'Export data in  %s', outf1)
   os.rename(tmp_out_pd_file,outf1)
   os.system('gzip -f '+ outf1)
   if not os.path.exists(out_mr_file_path):
     os.makedirs(out_mr_file_path)
   os.system('mv ' + outf1+'.gz' + ' ' + out_mr_file_path)

   return block, block_size, block_methylation, block_length, block_methylation_columns

def convert_position2range(all_data_df_position):
   ''' this function is to convert long CpG position to a range with the number of data pointes'''
   all_positions=all_data_df_position.split(',')
   new_position_range=all_positions[0]+','+all_positions[1]+'-'+ all_positions[-1]+','+ str(len(all_positions)-1)
   return new_position_range

def export_data_to_file(logger_n,max_loop,  record_passed_mr_hyper, record_passed_mr_hypo, record_passed_mr_mix,out_file_name,is_modT ) :
   #dump all data in json file
   #may not need to do it
   if False:
     with open(out_file_name +'_'+ str(max_loop) + '_hyper.json', 'w') as outfile:
       json.dump(record_passed_mr_hyper, outfile)
  
     with open(out_file_name +'_' + str(max_loop) +  '_hypo.json', 'w') as outfile:
       json.dump(record_passed_mr_hypo, outfile)

     with open(out_file_name + '_' + str(max_loop) + '_mix.json', 'w') as outfile:
       json.dump(record_passed_mr_mix, outfile)
     #print "Export hyper, hypo, mix, MRs in file ", outfile

   if is_modT==0:
       testStr='T-test'
   elif is_modT==1:
       testStr='modT-test'
   elif is_modT==2:
       testStr='KS-test'

   #export DMR in csv file
   column_names= ['mr_id',testStr+'_pval_smoothed_data',testStr+'_pval_interpolated_data','percent_data_passed_ttest',
                   'gcb_vs_grpsDist_pval','tumor_vs_grpsDist_pval','gcb_vs_grpsDist_tval','tumor_vs_grpsDist_tval','cluster_accuracy',
                  'low_negative_tumor_vs_gcb_percent','median_negative_tumor_vs_gcb_percent','high_negative_tumor_vs_gcb_percent',
                  'low_positive_tumor_vs_gcb_percent','median_positive_tumor_vs_gcb_percent','high_positive_tumor_vs_gcb_percent',
                   'is_DMR', 'position','DMR_type']

   all_mr_id=[]
   all_mr_type=[]
   all_mr_data=[]
   #here hyper ,hypo is based on the high perentage cutoff changes
   for k in record_passed_mr_hyper.keys():
      all_mr_id.append(k)
      all_mr_type.append('hyper')
      all_mr_data.append([k]+record_passed_mr_hyper[k][1]+ [','.join(record_passed_mr_hyper[k][2])]  + ['hyper'])

   for k in record_passed_mr_hypo.keys():
      all_mr_id.append(k)
      all_mr_type.append('hypo')
      all_mr_data.append([k] + record_passed_mr_hypo[k][1] +  [','.join(record_passed_mr_hypo[k][2])]  +  ['hypo'])


   for k in record_passed_mr_mix.keys():
      all_mr_id.append(k)
      all_mr_type.append('mix')
      all_mr_data.append([k] + record_passed_mr_mix[k][1] +  [','.join(record_passed_mr_mix[k][2])]  +  ['mix'])

   all_data_df=pd.DataFrame(data=all_mr_data,columns=column_names)
   #reduce positions to range of positions
   #all_data_df.position=all_data_df.position.apply(lambda x: convert_position2range(x) )
   all_data_df.to_csv(out_file_name,sep='\t', index=False,float_format='%g')
   logger_n.info( 'Export results at : %s ', out_file_name)

   return all_data_df


def read_json_data(in_hyper_file, in_hypo_file, in_mix_file):
   with open(in_mix_file) as f:
     record_passed_mr_mix = json.load(f)

   for kk in record_passed_mr_mix.keys() :
     tmp_data= np.array( record_passed_mr_mix[kk])
     tmp_df= pd.DataFrame(data=tmp_data, columns=['count','percentage','bin_start','bin_end'])
     print(tmp_df)

   with open(in_hyper_file) as f:
     record_passed_mr_hyper = json.load(f)
   for kk in record_passed_mr_hyper.keys() :
     tmp_data= np.array( record_passed_mr_hyper[kk])
     tmp_df= pd.DataFrame(data=tmp_data, columns=['count','percentage','bin_start','bin_end'])
     print(tmp_df)
     #raw_input('')

   with open(in_hypo_file) as f:
     record_passed_mr_hypo = json.load(f)
   for kk in record_passed_mr_hypo.keys() :
     tmp_data= np.array( record_passed_mr_hypo[kk])
     tmp_df= pd.DataFrame(data=tmp_data, columns=['count','percentage','bin_start','bin_end'])
     print(tmp_df)
     #raw_input('')
   print("Read hyper, hypo, mix MR data from ", in_hyper_file, in_hypo_file , in_mix_file)


def setup_logger(out_folder, in_chrm,max_read_length,mini_block_size,log_number):
   ''' set up logger files for different processes
   '''
   log_file= os.path.join(out_folder , in_chrm, in_chrm + '_maxlength_'+str(max_read_length)+ '_miniSize_'+ str(mini_block_size)+'_myapp.log_'+ str(log_number))
   if not os.path.exists(os.path.join(out_folder ,in_chrm)):
     try:
       os.makedirs(os.path.join(out_folder , in_chrm))
     except Exception:
       pass
  
   logger= logging.getLogger(str(log_number))
   file_handler=logging.FileHandler(log_file,mode='w')
   formatter= logging.Formatter('%(asctime)s %(levelname)-8s %(message)s', datefmt='%a, %d %b %Y %H:%M:%S')
   file_handler.setFormatter(formatter)  
   logger.setLevel(logging.INFO)
   logger.propagate=False
   logger.addHandler(file_handler)
   
   #logging.addHandler(logging.FileHandler(log_file))
   #logging.basicConfig(level=logging.INFO,
   #                 format='%(asctime)s %(levelname)-8s %(message)s',
   #                 datefmt='%a, %d %b %Y %H:%M:%S',
   #                 filename=log_file,
   #                 filemode='w')
   print("Log file:", log_file)
   #logging print to screen
   #logging.getLogger().addHandler(logging.StreamHandler(sys.stdout))
   #loggin print to a file
   #logger1= logging.getLogger(str(log_number))
   #logger1.addHandler(logging.FileHandler(log_file))
   return logger


