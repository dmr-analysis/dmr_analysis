import logging
import numpy as np
import pandas as pd
import matplotlib as mlt
mlt.use('Agg')
import matplotlib.pyplot as plt 
import argparse
import os
import glob
import logging
import multiprocessing as mp
import math
import warnings
#from multiprocessing.pool import Pool
warnings.filterwarnings("ignore")

from .script_high.dmr_utility import setup_logger, combine_multiple_files_to_one_dataframe,find_blocks_in_chrom, export_data_to_file, convert_position2range
from .script_high.dmr_data_analysis import do_DMR_analysis_in_a_block, cacl_neg_pos_percent_changes, replace_nan_by_median, calculate_mean_and_std, differential_analaysis_in_two_groups, difference_of_eucdist, do_test_between_two_groups, do_PCA_analysis,do_tSNE_analysis,do_kmeans_on_tSNE
#import dmr_utility
#import dmr_data_analysis
#import dmr_plots
from .script_high.dmr_plots import plot_raw_tumor_data, plot_raw_gcb_data, plot_3D_PCA, plot_2D_tSNE, plot_smoothed_data, plot_hist_of_difference_in_groups, plot_hist_for_block_length_and_data_size

EPS=np.finfo(np.float32).eps
#np.seterr(all='raise')
np.seterr(all='ignore')

def my_parser(parser):
   required=parser.add_argument_group("Required")
   required.add_argument('-in','--in_file_folder',
                         required=True,
                         help='input file folder for DNA methylation data such as WGBS. '
                              'In this folder, all samples are provided in each chromosome folder with BED format where sample name is indicated by file name. ')
   required.add_argument('-chr','--chromosome',
                          required=True,
                          help='select a chromosome for running dmr_analysis_block.'
                               ' For example, a chromosome (e.g., chrY) file folder under --in_file_folder that contains methylation data of samples in a chromosome (e.g., chrY)')
   required.add_argument('-gkey','--group_key',
                          required=True,
                          help='group key name, all bed files name contain this group_key will be combined together for dmr_analysis_block.'
                               'In other words, only bed file name with this group_key will be selected in analysis.'
                              ' Usually, it is the same as the file folder (or chromosome name) in --chromosome. ')

   optional = parser.add_argument_group('Optional, has default values')
   optional.add_argument('-out','--out_file_folder',
                         default='out',
                         metavar='',
                         help='output file folder that store all results computed by dmr_analysis_block, default is out/')
   optional.add_argument('-ncol','--need_column',
                         default='Methylation',
                         type=str,
                         metavar='',
                         help='select a columan name from bed file that will be used by dmr_analysis_block. '
                              ' For example, there are only six columns allowed in a bed file and the column labels will be added automatically '
                              ' such as (Chrs, Starts, Ends, Methylation, Total_counts, Strand) after loading the data.'
                              ' Here, if we assumes the fourth column of the input bed file is the methylation level, then --need_column = Methylation, default = Methylation ')

   optional.add_argument('-wtStr','--wildType_fileString',
                         default='gcb',
                         type=str,
                         metavar='',
                         help='Use the first few character string of a file name to indicate it is a normal/wide type sample, or to labele this file/sample as wide type/normal condition.'
                              ' For example, if a file name under a chromosome folder of --in_file_folder starts with gcb_meht1_* is a wild type/control/normal sample'
                              ', then --wildType_fileString is gcb. Default is gcb in the program '
                        )

   optional.add_argument('-dstart','--data_start_position',
                         default=3,
                         type=int,
                         metavar='',
                         help=' data start position, a start column position for input dataframe after methylation levels of multiple samples are combined into one file , default is 3 for chrY demo.'
                              ' Please note, it starts from 0 index at here.'
                              ' for example, in file chrY_MR_data4maxBlockDistanc*.gz, column labels are started from "chrs", "pos_start", "pos_end", "data1"'
                              ' "data2", ..., "datam", where the first 3 columns are chromosome positions, then the --data_start_position is 3. '
                              ' From the data_start_position to the data_end_position are all availble data columns of combined samples. '
                              ' For example, there are 12 samples in chrY demo, then th data_start_position=3 and data_end_position=15, which is the default setting in the program.'
                              ' Usually, these two parameters do not need to be changed if data size of input file is the same format as chrY_MR_data4maxBlockDistance*.gz ')
   optional.add_argument('-dend','--data_end_position',
                         default=15,
                         type=int,
                         metavar='',
                         help=' data end position, an end column position for input dataframe after methylation levels of multiple samples are combined into one file, default is 15 for chrY demo. '
                              ' More information please refer to "--data_start_position". Please note that data_end_position - data_start_position = total number of input samples' )
   optional.add_argument('-mxL','--maximum_adjacency_length',
                         default=250,
                         type=int,
                         metavar='',
                         help='maximum length of adjancey CpG sites is allowed in a methylation region (MR) , default = 250 ')
   optional.add_argument('-minS','--minimum_block_size',
                         default=5,
                         type=int,
                         metavar='',
                         help='minimum number of CpG sites is requested in a block/methylatoin region (MR), default = 5')
   optional.add_argument('-nDP','--number_of_data_points',
                          default=0,
                          type=int,
                          metavar='',
                          help='the number of data points (or rows from a dataframe) will be considered in analysis, ' 
                               ' default=0 that means all data points (or the combined dataframe from all samples and CpG sites) will be used. '
                               ' if it sets =1000 then only the first 1000 rows of dataframe will be used in the analysis. This option is for debuging purforse when using a small sample for testing ')
   optional.add_argument('-pC','--P_cutoff',
                          default=0.05,
                          type=float,
                          metavar='',
                          help='P cutoff value for T-test or other statistic significance test used in the analysis , default is 0.05 ')
   optional.add_argument('-aC','--accuracy_cutoff',
                         default=0.5,
                         type=float,
                         metavar='',
                         help='Clustering accuracy cutoff value for binary classification (e.g., DMR or not DMR; tumor or normal samples) of 2D t-SNE map, '
                              ' default is > 0.5. This option will be removed latter because it does not affect the overall results!')
   optional.add_argument('-mPer','--minimum_percentage_changes',
                         default=0.0001,
                         type=float,
                         metavar='',
                         help=' minimum percengate of data points in a MR that passed a predifined filtering condition such as the methylation changes greater than a predefined cutoff value (e.g. -perC=0.1), '
                               ' default is 0.0001 = 0.1 percentage ')
   optional.add_argument('-perC','--percentage_cutoff',
                         default='0.07,0.15,0.2',
                         type=str,
                         metavar='',
                         help=' a comma separated string list for predefined cutoff values (e.g., 0.07,0.15,0.2) of percentage of methylation changes between the two groups is low, median, high, respectively. '
                              ' For example, in smoothed data points, a difference of methylation leveles between the tumor and the normal samples can be set >0.07, 0.15, >0.2 for low, median, high changes, respectively. '
                              '  where methylation levels are between 0 and 1. Default parameter is 0.07,0.15,0.2 ')
   optional.add_argument('-lmhC','--low_median_high_cutoff',
                         default=2,
                         type=int,
                         metavar='',
                         help=' Use low, median, or high percentage changes (e.g., one of the inputs from --percentage_cutoff) to predict high confidence DMRs.  '
                              ' Here, use 0, 1, or 2 to represent low, median, or high percentage cutoff for the prediction. '
                              ' default is 2 (high percentage changes in --percentage_cutoff)')
   optional.add_argument('-numPro', '--number_of_processes',
                          default=1,
                          type=int,
                          metavar='',
                          help='Number of parallel rocesses will be used in calculation, default is 1')
   optional.add_argument('-isSmd','--is_smoothed_data',
                          default=2,
                          type=int,
                          metavar='',
                          help =' Here, use 0, 1, or 2 to represent raw data, interpolated data, or smoothed data are used, respectively, in predicting DMRs when compared the two groups, '
                                ' default is 2 - use smoothed data points to predict DMRs' )
   optional.add_argument('-isMT','--is_moderate_ttest',
                          default=0,
                          type=int,
                          metavar='',
                          help =' 0 is standard T-test, 1 is moderate T-test, and 2 is KS-test for evaluting the significance of differential methylation between the two groups'
                                ' across all CpG sites in a MR. Default is 0, uses a T-test to test the significance. The moderate T-test is slower than T-test' )
   optional.add_argument('-isExt','--is_export_data',
                          default=0,
                          type=int,
                          metavar='',
                          help ='0 is not exporting data, 1 is exporting data during the analysis, where both raw and smoothed methylation levels of all DMRs will be exported, default is 0' )
   optional.add_argument('-dotOrUnderscore','--column_splitBy_dotOrUnderscore', default=0,
                         type=int, metavar='', help= '0 for using dot . to split column labels, 1 for using underscore _ to split column labels, default=0 using dot to split column labels')

   return parser


def generate_output_file_names(tmp_out_file_name, tmp_out_path,log_nr_str):
   ''' Generate output file names and a dataframe for parameters 
	based on old format of file names
   '''
   tmp_out_file0=tmp_out_file_name
   tmp_out_data_string=tmp_out_file0.split('_')
   tmp_out_file=tmp_out_data_string[0]+'_'+log_nr_str + '_mr_data.tsv'
   tmp_range_out_file=tmp_out_data_string[0]+ '_' + log_nr_str + '_mr_data_range.tsv'
   tmp_parameter_out_file=tmp_out_data_string[0]+ '_'+ log_nr_str + '_mr_parameters.tsv'

   #make dataframe for parameters
   parameter_df=pd.DataFrame(data=tmp_out_data_string+[tmp_out_file0])

   #export data and parameters
   out_file=os.path.join(tmp_out_path, tmp_out_file)
   range_out_file=os.path.join(tmp_out_path, tmp_range_out_file)
   parameter_out_file=os.path.join(tmp_out_path, tmp_parameter_out_file)
   return out_file, range_out_file, parameter_out_file, parameter_df.copy()


def do_block_finding(log_number,in_folder, out_folder, in_chrm, group_key, needed_columns, num_of_data_points, data_start_pos, data_end_pos, max_read_length,mini_block_size,wildType_fileString):
   '''
      do MR lock finding from input data, this part does not need to perform parallel computation. 
   '''

   logger0=setup_logger(out_folder, in_chrm,max_read_length,mini_block_size,log_number)

   #first combine multiple files of the same chromosome data into one dataframe
   #data file is the same for all type of max_read_length and mini_block_size
   tmp_data_file=os.path.join(out_folder , in_chrm,in_chrm+ '_MR_data4maxBlockDistance_*')
   tmp_in_file=glob.glob(tmp_data_file)
   if len(tmp_in_file)==0:
      sorted_tmp_pd2, positions, methylations , out_tmp_pd2_file = combine_multiple_files_to_one_dataframe(in_folder,in_chrm, group_key, needed_columns, logger0)
   else:
      logger0.info('File load %s ', tmp_in_file)
      os.system('gunzip '+ tmp_in_file[0])
      os.rename(tmp_in_file[0].replace('.gz',''),in_chrm+ '_tmp_out_pd.csv')

   #Then try to find all possible blocks in the chromome based on available data points in the chromoeome
   out_mr_file_name= os.path.join(out_folder , in_chrm + '_MR_data4maxBlockDistance_' + str(max_read_length)+ '_minBlockSize_' +str(mini_block_size) + '.csv')
   out_mr_file_path= os.path.join(out_folder , in_chrm )
   block, block_size, block_methylation, block_length, block_methylation_columns = find_blocks_in_chrom(in_chrm, num_of_data_points, data_start_pos, data_end_pos, 
                                             max_read_length, mini_block_size,out_mr_file_name,out_mr_file_path,logger0)

   #plot historgram of block length, and block data points size
   out_fig_name=os.path.join(out_mr_file_path,in_chrm+'hist_mr_length_'+ str(max_read_length) + '_and_size_'+str(mini_block_size) +'.pdf')
   plot_hist_for_block_length_and_data_size(block_length,block_size,max_read_length,out_fig_name,logger0)

   #passed_mr=0
   #find gcb in columns
   tmp_cols=block_methylation_columns.astype(str)
   #jbw
   #print(tmp_cols, wildType_fileString)
   gcb_col_idx=np.where(np.char.find(tmp_cols,wildType_fileString)>=0)
   tumor_col_idx=np.where(np.char.find(tmp_cols,wildType_fileString)<0)
   #print(gcb_col_idx, tumor_col_idx)
   logger0.info("Wild type /control sample file name is %s ", wildType_fileString)
   logger0.info("Wild/control sample %s , ", str(len(gcb_col_idx[0])))
   logger0.info("Tumor/KO sample %s , ", str(len(tumor_col_idx[0])))

   #for each MR perfrom data interpretation
   max_loop= len(block)
   #loop =0
   export_fig_path=os.path.join(out_folder , in_chrm , "plots")
   if not os.path.exists(export_fig_path):
      os.makedirs(export_fig_path)
   logger0.info("DMR export path %s", export_fig_path )

   export_data_path=os.path.join(out_folder, in_chrm, "data")
   if not os.path.exists(export_data_path):
      os.makedirs(export_data_path)
   logger0.info("DMR export MR data path %s", export_data_path)

   return logger0, max_loop,tmp_cols, export_fig_path, block_methylation, block, block_size, block_length, tumor_col_idx, gcb_col_idx, block_methylation_columns

def main(in_folder, out_folder, in_chrm, group_key, needed_columns,wildType_fileString,
         num_of_data_points, data_start_pos, data_end_pos, max_read_length,mini_block_size,
         P_cutoff, accuracy_cutoff, mini_percent, percent_cutoff, isSmooth, isExport,dotOrUnderscore, is_modT , low_median_high_cutoff):

   log_number=0
   logger_n,max_loop,tmp_cols, export_fig_path, block_methylation, block, block_size, block_length, tumor_col_idx, \
       gcb_col_idx, block_methylation_columns = do_block_finding(log_number,in_folder, out_folder, in_chrm,\
                                                  group_key, needed_columns, num_of_data_points, data_start_pos, data_end_pos, max_read_length,mini_block_size,wildType_fileString)

   is_parallel=False
   call_do_DMR_analysis_in_block((logger_n, block, max_loop, tmp_cols, export_fig_path,
                       block_methylation,block,block_size,block_length,
                       tumor_col_idx,gcb_col_idx,percent_cutoff,
                       block_methylation_columns,isSmooth,isExport,dotOrUnderscore,P_cutoff,is_modT,low_median_high_cutoff,
		out_folder,in_chrm,max_read_length,mini_block_size,mini_percent,accuracy_cutoff, wildType_fileString, is_parallel) ) 

def call_do_DMR_analysis_in_block(args):
   '''
     do DMR for each chunk of blocks
   '''
   logger_n, chunk_block, max_loop, tmp_cols, export_fig_path, \
              block_methylation,block,block_size,block_length,\
              tumor_col_idx,gcb_col_idx,percent_cutoff,\
              block_methylation_columns,isSmooth,isExport,dotOrUnderscore,P_cutoff, \
              is_modT,low_median_high_cutoff, out_folder,in_chrm,max_read_length,mini_block_size,mini_percent,accuracy_cutoff, wildType_fileString, is_parallel = args
   record_passed_mr_hyper={}
   record_passed_mr_hypo={}
   record_passed_mr_mix={}

   if isinstance(logger_n, int):
      log_number=logger_n
      logger_n=setup_logger(out_folder, in_chrm,max_read_length,mini_block_size,logger_n)   
   else:
      log_number=0

   #loop in each block 
   passed_mr=0
   loop=0
   for kk in chunk_block:
     if loop< max_loop:
        loop +=1
     else:
        break
     tmp_id=kk
     tmp_data=block_methylation[kk]
     tmp_pos=block[kk]
     logger_n.info('\n')
     logger_n.info(kk)
     #print(kk)

     #print tmp_pos
     tmp_num_of_data=block_size[kk]
     tmp_len=block_length[kk]

     #plot raw data points
     #here prepare for sorting the data based on position if it is needed
     #this may not needed if interp1d function will sort x before fitting

     #plot tumor
     selected_data= tmp_data[:,tumor_col_idx[0]]
     selected_data2tumor= replace_nan_by_median(selected_data,tmp_id)
     selected_cols= tmp_cols[tumor_col_idx[0]]
     selected_pos2= tmp_pos

     #plot gcb
     selected_data= tmp_data[:,gcb_col_idx[0]]
     #replace nan by zero or median of other availale data
     selected_data2gcb= replace_nan_by_median(selected_data,tmp_id)
     selected_cols= tmp_cols[gcb_col_idx[0]]
     #selected_pos2gcb= tmp_pos

     #DMR analysis
     fig=plt.figure(figsize=(18,7))
     statistic_test_p, cluster_accuracy, negative_percent2, positive_percent2,results2 , total_percent2 = do_DMR_analysis_in_a_block(logger_n, percent_cutoff, selected_data2tumor,
                selected_pos2, selected_data2gcb,tmp_len, tmp_id, tmp_num_of_data, gcb_col_idx, tumor_col_idx, block_methylation_columns,
                fig,isSmooth,isExport,dotOrUnderscore,P_cutoff,is_modT,low_median_high_cutoff, wildType_fileString,out_folder,in_chrm)

     #in statistic_test_p, the first is smoothed P-value and the second is interpolated P-value, the third is percenage of data points passed T-test filtering, the fourth and fifth P-values
     #are the signifcance of euclain-distance between within group and across groups for gcb and tumor group, respectively
     logger_n.info('P values of KS-test on interpolated mean group data: %g', statistic_test_p[1])
     logger_n.info('P values of KS-test on smoothed mean group data: %g', statistic_test_p[0])
     logger_n.info('The significance of euclidean distance (P value) gcb vs the between groups: %g', statistic_test_p[3])
     logger_n.info('The significance of euclidean distance (P value) tumor vs the between groups: %g',statistic_test_p[4])
     str_cutoff=['low','median','high']
     logger_n.info('Minimum ' +str_cutoff[low_median_high_cutoff]+ ' '+ str(percent_cutoff[low_median_high_cutoff]) + ' percentage of the mean methylation changes between two groups across all positions >: %g', mini_percent)

     #always use P-value from smoothed data for selection
     #first filtering by using low, median, or high cutoff value for percenage changes bewteen the groups 
     if (total_percent2[low_median_high_cutoff] > mini_percent):
        if  ~(statistic_test_p[3]>=P_cutoff and statistic_test_p[4]>=P_cutoff) or (statistic_test_p[2]>0) :
            #second filtering by using the significant difference of euclidean distance withitn group vs. between groups
                  isDMR='D'
                  passed_mr +=1
        else:
            isDMR='U'
            logger_n.info("Not pass: %s, %g , %g, %g , %g", tmp_id , P_cutoff , cluster_accuracy, statistic_test_p[0], statistic_test_p[1])
     else:
        isDMR='U'
        logger_n.info("Not pass: %s, %g , %g, %g , %g", tmp_id , P_cutoff , cluster_accuracy, statistic_test_p[0], statistic_test_p[1])
     plt.close()
     #input(" ")

     #record all MR, DMR marked as D, not DMR marked as U, this is a primary result which will be further filtered by a next function!
     #here only use high cutoff percenage as the cutoff for hyper, hypo
     if negative_percent2[low_median_high_cutoff] ==0 and positive_percent2[low_median_high_cutoff]>0 :
          record_passed_mr_hyper[tmp_id]= (results2.tolist(), statistic_test_p + [ cluster_accuracy]+ negative_percent2.tolist() + positive_percent2.tolist() +[isDMR], [in_chrm] + selected_pos2.astype(str).tolist() )
     elif negative_percent2[low_median_high_cutoff] >0 and positive_percent2[low_median_high_cutoff]==0 :
          record_passed_mr_hypo[tmp_id]= (results2.tolist(),  statistic_test_p + [ cluster_accuracy] +negative_percent2.tolist()+ positive_percent2.tolist()+[isDMR], [in_chrm] +selected_pos2.astype(str).tolist()  )
     else :
          record_passed_mr_mix[tmp_id]= (results2.tolist(),  statistic_test_p + [ cluster_accuracy] +negative_percent2.tolist()+ positive_percent2.tolist()+[ isDMR], [in_chrm] + selected_pos2.astype(str).tolist()  )

   #export
   if isSmooth==2:
       logger_n.info("Export percentage changes based on smoothed data .")
   elif isSmooth==1:
       logger_n.info("Export percentage changes based on interpolated data .")
   else:
       logger_n.info("Export percentage changs based on raw data .")

   #here out file name shall be changed and the parameters in original file name shall be moved to a new file
   str_percent_cutoff=[str(ii) for ii in percent_cutoff]
   str_percent_cutoff="_".join(str_percent_cutoff + [ str_cutoff[low_median_high_cutoff]])
   out_file_name=in_chrm +'_maxDist_' + str(max_read_length)+ '_minSize_' +str(mini_block_size) + '_DMR_clusterAccuracy_gt_' \
                + str(accuracy_cutoff) + '_miniMethyChange_gt_'+ str_percent_cutoff+ '_miniPercentChange_gt_'  + str(mini_percent) + '_Pcutoff_' \
                + str(P_cutoff)+ '_isSmooth_'+str(isSmooth) + '_isModTest_'+ str(is_modT)+ '_'+str(max_loop) +'_'+str(log_number)
   
   #add jbw 2023
   #generate output filename and path for data and parameters
   #here out file name shall be changed to a new name and the parameters in original file name shall be moved to a new file
   #out_file_name=os.path.join(export_fig_path, out_file_name)
 
   tmp_out_path=export_fig_path
   tmp_out_file0=out_file_name
   out_file_name, range_out_file, parameter_out_file, parameter_df= generate_output_file_names(tmp_out_file0, tmp_out_path, str(log_number))

   #export all data add jbw 2023
   all_dmr_data_df=export_data_to_file(logger_n,max_loop, record_passed_mr_hyper, record_passed_mr_hypo, record_passed_mr_mix, out_file_name,is_modT)
  
   #only export range data for single process calculation 
   if not is_parallel:
     range_pd=all_dmr_data_df.copy()
     range_pd.position= range_pd.position.apply(lambda x: convert_position2range(x)) 
     range_pd.to_csv(range_out_file,sep='\t',index=False, float_format='%g')

   #only export parameter_df in single process calculation or the first process in parallel calculation
   if not is_parallel or log_number==1 :
     parameter_df.to_csv(parameter_out_file, sep='\t',index=False, header=None, float_format='%g')
     
   return out_file_name

def parallel_main(in_folder, out_folder, in_chrm, group_key, needed_columns,wildType_fileString, num_of_data_points, data_start_pos, data_end_pos, max_read_length,mini_block_size, 
               P_cutoff, accuracy_cutoff, mini_percent, percent_cutoff, isSmooth, isExport,dotOrUnderscore, is_modT , low_median_high_cutoff, num_of_processes):
   '''parallel version of predicting DMRs

   '''
   log_number=0
   logger_n,max_loop,tmp_cols, export_fig_path, block_methylation, block, block_size, block_length, tumor_col_idx, \
       gcb_col_idx, block_methylation_columns = do_block_finding(log_number,in_folder, out_folder, in_chrm,\
                                                  group_key, needed_columns, num_of_data_points, data_start_pos, data_end_pos, max_read_length,mini_block_size,wildType_fileString)

   #equally distribut block keys to num_of_processes
   all_blocks=list(block.keys())
   num_in_chunks=int(math.ceil(len(all_blocks)/num_of_processes))
   block_chunks = [all_blocks[x:x+num_in_chunks] for x in range(0, len(all_blocks), num_in_chunks)]

   #do parallel DMR in each chunk of blocks
   if len(block_chunks) > num_of_processes:
      num_of_processes=len(block_chunks)
   elif len(block_chunks) < num_of_processes:
      logger_n.info('Number of blocks smaller than the number of available processes ')
      num_of_processes=len(block_chunks)
   elif len(block_chunks)==1:
      num_of_processes= 1
   logger_n.info("Do parallel calculation by using %g processes ", num_of_processes) 

   pool = mp.Pool(processes=num_of_processes)
   log_number=range(1,num_of_processes+1)
   #jbw for debug
   #print(log_number, num_of_processes, len(block_chunks))
   #for loop in range(0, num_of_processes):
   #    print(log_number[loop], block_chunks[loop])
   #end debug
   is_parallel=True
   files= pool.map(call_do_DMR_analysis_in_block, [ (log_number[loop], block_chunks[loop], max_loop, tmp_cols, export_fig_path,
                       block_methylation,block,block_size,block_length,
                       tumor_col_idx,gcb_col_idx,percent_cutoff,
                       block_methylation_columns,isSmooth,isExport,dotOrUnderscore, P_cutoff,
                       is_modT,low_median_high_cutoff,out_folder,in_chrm,max_read_length,
                       mini_block_size,mini_percent,accuracy_cutoff,wildType_fileString,is_parallel) for loop in range(0, num_of_processes) ],1)

   #start to merge files
   files= list(filter(None, files))
   record_pd=[]
   for fi in files:
      record_pd.append(pd.read_csv(fi,sep='\t'))
   all_pd=pd.concat(record_pd)
   range_pd=all_pd.copy()
   range_pd.position= range_pd.position.apply(lambda x: convert_position2range(x))

   #below two file names shall be changed and the parameters in original file name shall be moved into a separted file
   out_file_str=files[0].split('_')
   #out file at 'position' column contains all positions of methylation sites in a MR region, which are separated by ,
   out_file='_'.join( out_file_str[0:-1]+ ['all'])
   #range_out_file at 'position' column contains a range of MR position such as chr1,start-end,number of methylation sites
   range_out_file='_'.join(out_file_str[0:-1]+ ['range'])

   #added jbw 2023
   #find file path and name
   out_paths =os.path.split(out_file)
   tmp_out_path=out_paths[0]
   tmp_out_file0=out_paths[1]

   #in parallel calculation does not export the same parameter_df repeatedly, only once in process 1 
   out_file, range_out_file, parameter_out_file, parameter_df= generate_output_file_names(tmp_out_file0, tmp_out_path,'all')
   logger_n.info("Export all position results at : " + out_file)
   logger_n.info("Export range position results at : " + range_out_file)

   all_pd.to_csv(out_file,sep='\t',index=False, float_format='%g')
   range_pd.to_csv(range_out_file,sep='\t',index=False, float_format='%g')

   #remove the rest of files
   for fi in files:
      os.remove(fi) 
   #end testing
   return out_file

def run(args):
  in_folder=args.in_file_folder
  out_folder=args.out_file_folder
  in_chrm=args.chromosome
  group_key=args.group_key
  wildType_fileString=args.wildType_fileString
  dotOrUnderscore=args.column_splitBy_dotOrUnderscore


  needed_columns=args.need_column
  if args.number_of_data_points ==0:
     num_of_data_points=None
  else:
     num_of_data_points=args.number_of_data_points
  data_start_pos=args.data_start_position
  data_end_pos=args.data_end_position
  max_read_length=args.maximum_adjacency_length
  mini_block_size=args.minimum_block_size

  P_cutoff=args.P_cutoff
  accuracy_cutoff=args.accuracy_cutoff
  mini_percent=args.minimum_percentage_changes
  percent_cutoff=args.percentage_cutoff
  percent_cutoff=np.asfarray(percent_cutoff.split(','))
  #here contains Low, median, High percentage of changes cutoff values
  #print(percent_cutoff)
  low_median_high_cutoff=args.low_median_high_cutoff
  num_of_processes=args.number_of_processes

  #if args.is_smoothed_data.lower() in ('yes','y'): 
  #   isSmooth=True
  #else:
  #   isSmooth=False
  isSmooth= args.is_smoothed_data
  is_modT=args.is_moderate_ttest
  if args.is_export_data ==0:
     isExport=False
  else:
     isExport=True
  #print(in_folder, out_folder, in_chrm, group_key, needed_columns,
  #       num_of_data_points, data_start_pos, data_end_pos, max_read_length,mini_block_size,
  #       P_cutoff, accuracy_cutoff, mini_percent, percent_cutoff, isSmooth)

  if num_of_processes==1:
     main(in_folder, out_folder, in_chrm, group_key, needed_columns,wildType_fileString,
         num_of_data_points, data_start_pos, data_end_pos, max_read_length,mini_block_size,
         P_cutoff, accuracy_cutoff, mini_percent, percent_cutoff, isSmooth, isExport,dotOrUnderscore, is_modT, low_median_high_cutoff )
  else:
     parallel_main(in_folder, out_folder, in_chrm, group_key, needed_columns,wildType_fileString,
         num_of_data_points, data_start_pos, data_end_pos, max_read_length,mini_block_size,
         P_cutoff, accuracy_cutoff, mini_percent, percent_cutoff, isSmooth, isExport,dotOrUnderscore, is_modT , low_median_high_cutoff, num_of_processes)

if __name__== '__main__':
  #here for test run
  #in_folder='in_data/WGBS-data/'
  #out_folder='out/DMR/'
  #in_chrm='chrY'
  #group_key='chrY'
  #needed_columns='Methylation'
  #num_of_data_points=None
  #data_start_pos=3
  #data_end_pos=15
  #max_read_length=250
  #mini_block_size=5
  #P_cutoff=1.0
  #accuracy_cutoff=0.5
  #mini_percent=0.05
  #percent_cutoff=0.05
  #isSmooth=False

  #here for real run
  #np.warnings.filterwarnings('ignore', category=np.VisibleDeprecationWarning)

  args=my_parser(argparse.ArgumentParser('python dmr_analysis_block.py')).parse_args()
  run(args)




