#combine results from multiple chromosomes DMR for ranking by using dmr_rank.py
import glob
import os
import pandas as pd
import numpy as np

#in-house functions
#import rank_dmr
#import dmr_data_analysis
#import plot_dmr
from .script_high.rank_dmr import logRegress_feature_score
from .script_high.plot_dmr import select_DMRs_for_analysis
from .script_high.dmr_data_analysis import accuracy
from .script_high.dmr_utility import chrStr2numeric


#exec(open('combine_multChrs_dmr4ranking.py').read())
def my_parser(parser):
  required=parser.add_argument_group("Required")
  required.add_argument('-inChrs','--in_chroms_number', required=True, help="a string of chromosomes that need be combined such as chr1,chr2,chr3")
  required.add_argument('-inFold','--in_file_folder', required=True, help="a file folder of MR data in each chromosome that is exported by dmr_analysis_block ")

  optional=parser.add_argument_group("Optional, has default values")
  optional.add_argument('-inSFold','--in_File_subFolder', default="plots",metavar='',help="a subfolder of exported MR data in each chromosome based on dmr_analysis_block, default = plots that is under each chromosome folder ")
  optional.add_argument('-inMPer','--in_minimum_percentage_changes',default=0.0001,metavar='',type=float,help='parameter of minimum percentage of data points in a MR with a predefined methylation changes greater than a cutoff, default=0.0001 ')
  optional.add_argument('-inPC','--in_Pvalue_cutoff', default=0.05, metavar='',type=float, help="P value cutoff for significance test, default = 0.05")
  optional.add_argument('-inIsSmd','--in_is_smoothed_data', default=0, metavar='',type=int, help="is result based on raw =0, interpolated =1 , or smoothed =2 data, default= 0")
  optional.add_argument('-inAC','--in_accuracy_cutoff_range', default='0.0,1.1',metavar='',type=str, help="range of clustering accurancy, default include all (e.g. from 0.0  to 1.0) is 0.0,1.1")
  optional.add_argument('-inDMRst','--in_DMR_string_type',default='D',metavar='', type=str,help='a string used to represent predicted DMR in the file, default is D for prediced DMR')
  optional.add_argument('-inLRpb','--in_LogReg_proba', default=0.8, metavar='',type=float, help="a probability value used by logistic Regression to select DMRs, default =0.8")
  optional.add_argument('-inMT','--in_moderate_ttest', default=0,metavar='', type=int, help="0 for standard T-test, 1 for moderate T-test, and 2 for KS-test for evaluating the test of significance, default=0")
  optional.add_argument('-inLMH','--in_low_median_high_cutoff',default='high', type=str, help="use high, median, or low minimum percentage change for DMR ranking, default = high ")
  optional.add_argument('-inFEstr','--in_file_ending_string',default='_all',type=str,help="file name ending string that will be used to search for results file in result folder (e.g., chrY/plots/), default = _all ")
  return parser

def read_chr_dmr_files(in_chrs,in_folder,in_sub_folder, file_string ):
  #in_chrs: selected chromosome name
  #in_folder: input data folder
  #file_string: file string that will be used in search file names for combining
  #read input files
  #in_chrs=('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12',\
  #'chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY')
  #in_folder='out/DMR/'
  
  files2combine={}
  loop=0
  for i in in_chrs:
     #tmp_folder=os.path.join(in_folder,i,'plots','*_high_miniPercentChange_gt_0.0001_Pcutoff_0.05_isSmooth_2_isModTest_0_*_all')
     tmp_folder=os.path.join(in_folder,i,in_sub_folder,file_string)
     #print(tmp_folder)
     tmp_file=glob.glob(tmp_folder)
     files2combine[i]=tmp_file
     loop += len(tmp_file)

  print(str(len(in_chrs))+ ' chromes '+ str(loop) + ' '+ file_string + ' files find in '+ in_folder )

  #print(files2combine)
  return files2combine

def combine_dmr_files2df(files2combine):
  #combine files to one dataframe and add a new column for chroms
  all_chrs_df=[]
  for key in files2combine.keys():
    tmp_file=files2combine[key]
    #print(tmp_file)
    tmp_df=[]
    if len(tmp_file)==1:
       tmp_df=pd.read_csv(tmp_file[0],sep='\t')
       tmp_df['chroms']=key 
    else:
       print(tmp_file)
       print('Multiple files are found please check ' + tmp_file)
    all_chrs_df.append(tmp_df)
    #print(tmp_df.shape)

  all_dmr_df=pd.concat(all_chrs_df)
  return all_dmr_df


def export_logReg_score_in_dmr(tmp_dmr_df,in_chrs,files2combine, out_percent, logReg_proba):
  '''export results with the final rank score to each chroms file folder'''
  all_logReg_score=[]
  for i in in_chrs:
     tmp_dmr_df2=tmp_dmr_df.copy()
     tmp_dmr_df2=tmp_dmr_df2[tmp_dmr_df2.chroms==i]
     #here the out file name shall be changed and the parameters in original file names shall be moved to a new file
     #tmp_out_file=files2combine[i][0]+'_dmrRanking_top_' + '{:8.2f}'.format(out_percent).strip() + '_minLogReg_proba_' + '{:8.1f}'.format(logReg_proba).strip()
  
     #add jbw 2023 generate out file name
     tmp_out_file0= files2combine[i][0].replace('.tsv','')
     tmp_out_file= tmp_out_file0 + '_dmrRanking.tsv'
     tmp_parameter=['dmrRanking_top' , '{:8.2f}'.format(out_percent).strip() , 'minLogReg_proba' , '{:8.1f}'.format(logReg_proba).strip()]
     tmp_parameter_file=tmp_out_file0 +'_dmrRanking_parameters.tsv'
     tmp_parameter_df=pd.DataFrame(data=tmp_parameter)
     
     print(tmp_dmr_df2.shape)
     #sort reaults
     tmp_dmr_df2=tmp_dmr_df2.sort_values(by='logReg_score')
     tmp_dmr_df2.reset_index(drop=True)
     print('Export data: ' + tmp_out_file)
     tmp_dmr_df2.to_csv(tmp_out_file,index=False,sep='\t')

     print('Export parameter: ' + tmp_parameter_file)
     tmp_parameter_df.to_csv(tmp_parameter_file, index=False,sep='\t',header=None, float_format='%g')
     #store data
     all_logReg_score.append(tmp_dmr_df2)
  all_logReg_df=pd.concat(all_logReg_score)
  return all_logReg_df
 
def find_start_end_position(string):
   #find chr, start, end position from a string
   #add 10 bp flank sequence to two sides of MR
   flank_post=1
   tmp_post=string.split(',')
   tmp_chr=tmp_post[0]
   if '-' in tmp_post[1]:
     #mr range
     tmp_start,tmp_end=np.array(tmp_post[1].split('-')).astype(int)
     tmp_start -= flank_post
     tmp_end += flank_post
   else:
     #mr cg positions
     tmp_posts=np.asarray(tmp_post[1:]).astype(int)
     tmp_start=tmp_posts.min() -flank_post
     tmp_end=tmp_posts.max() + flank_post

   #out_string=':'.join([str(tmp_chr),str(tmp_start),str(tmp_end)])
   return tmp_chr, tmp_start,tmp_end

#def chrStr2numeric(string, human):
#   #convert string chr to numeric value for sorting
#   string=string.lower()
#   string=string.replace('chr','')
#   if human :
#     if string=='x':
#        string='23'
#     elif string=='y':
#        string='24'
#     elif string=='m':
#        string='25'
#   else:
#     if string =='x':
#        string='20'
#     elif string=='y':
#        string='21'
#     elif string=='m':
#        string='22'
#   return string

def sort_selected_columns(all_logReg_df):
  #export selected dataframe columns to a file
  #here * for unpack the elements in return variable
  selected_cols_logReg_df=all_logReg_df.copy()
  selected_cols_logReg_df=selected_cols_logReg_df[['mr_id','chroms','combined_id','is_DMR','position','DMR_type','logReg_predicted_dmr']]
  selected_cols_logReg_df['new_chr'],selected_cols_logReg_df['new_start'],selected_cols_logReg_df['new_end']=zip(*selected_cols_logReg_df['position'].apply(find_start_end_position))
  selected_cols_logReg_df=selected_cols_logReg_df[['new_chr','new_start','new_end','combined_id','position','DMR_type','is_DMR','logReg_predicted_dmr']]

  selected_cols_logReg_df['name_combined']=selected_cols_logReg_df[selected_cols_logReg_df.columns[3:-1]].apply(lambda x: ':'.join(x.astype(str)), axis=1)
  selected_cols_logReg_df['order_by_chr']=selected_cols_logReg_df['new_chr'].apply(chrStr2numeric,human=True)
  selected_cols_logReg_df['order_by_chr']=selected_cols_logReg_df['order_by_chr'].astype(int)
  #change the format of column position by removing all cg sites positions
  #jbw
  #print(selected_cols_logReg_df['position'] )
  selected_cols_logReg_df['name_combined']=selected_cols_logReg_df['name_combined'].apply(lambda x: ':'.join(np.array(x.split(':'))[[0,1,-2,-1]]) )

  selected_cols_logReg_df_sorted=selected_cols_logReg_df.copy()
  selected_cols_logReg_df_sorted=selected_cols_logReg_df.sort_values(['order_by_chr','new_start'],ascending=True)

  selected_cols_logReg_df_sorted=selected_cols_logReg_df_sorted.drop(['combined_id','position','DMR_type','is_DMR','order_by_chr'],axis=1)
  selected_cols_logReg_df_sorted=selected_cols_logReg_df_sorted[['new_chr','new_start','new_end','name_combined','logReg_predicted_dmr']]
  return selected_cols_logReg_df_sorted

def main(in_chrs, in_folder, in_sub_folder, in_Mpert,in_PVCut, in_Smoot, in_ModTt, file_string,ac_cutoff1,ac_cutoff2,str2dmr,logReg_proba,other_condition,parameter_string_list):
  files2combine= read_chr_dmr_files(in_chrs, in_folder, in_sub_folder, file_string)

  #combine files to one dataframe and add a new column for chroms
  #print(files2combine)
  all_dmr_df= combine_dmr_files2df(files2combine)

  #rank all dmr in the same dataframe
  tmp_dmr_df, conf_matrix =logRegress_feature_score(all_dmr_df)
  prediction_accuracy=accuracy(conf_matrix)
  print('Prediction accuracy:')
  print(prediction_accuracy)
  print('Confusion matrix:')
  print(conf_matrix)

  #check combined results
  in_not_passed_df2, needs_check,  (len_dmr,len_d,percent_recovered,num_of_top1,percent_of_top) = select_DMRs_for_analysis(ac_cutoff1, ac_cutoff2, str2dmr, other_condition, tmp_dmr_df)

  out_percent=percent_recovered #percent_of_top
  out_min_logReg_score=in_not_passed_df2.logReg_score.min()
  print('Selected mini-logReg-score and logReg probability')
  #print(out_min_logReg_score, in_not_passed_df2.loc[in_not_passed_df2.logReg_score==out_min_logReg_score,'logReg_predicted_dmr'].tolist())
  print(out_min_logReg_score, in_not_passed_df2.loc[in_not_passed_df2.logReg_predicted_dmr>=logReg_proba,['logReg_score','logReg_predicted_dmr']].min().tolist())
  #input('Please click any kew on keyboard to continue ?? ' + str(out_percent))
  print(out_percent)

  #export results with the final rank score to each chroms file folder
  all_logReg_df= export_logReg_score_in_dmr(tmp_dmr_df,in_chrs,files2combine,out_percent,logReg_proba)
  all_logReg_df['combined_id']=all_logReg_df['chroms'].str.cat(all_logReg_df['mr_id'].astype(str),sep=':')

  #here out file name shall be changed and the parameters in original file name shall be moved into  a new file
  #out_file1=str(len(in_chrs))+'_chroms'+ file_string + '_dmrRanking_top_' + '{:8.2f}'.format(out_percent).strip() + '_minLogReg_proba_' + '{:8.1f}'.format(logReg_proba).strip()
  #add jbw 2023
  out_file1= str(len(in_chrs))+'_chroms' + file_string.replace('.tsv','') +  '_dmrRanking.tsv'
  out_file1=out_file1.replace('*','')
  out_file1=out_file1.replace('.tsv','.bed')
  out_parameter_file1= out_file1.replace('.bed','_parameter.bed')
  out_parameter=[ str(len(in_chrs)) , 'chroms'] + parameter_string_list +  [ 'dmrRanking_top' , '{:8.2f}'.format(out_percent).strip() , 'minLogReg_proba' ,  '{:8.1f}'.format(logReg_proba).strip()]
  out_parameter=out_parameter + ['_'.join(out_parameter)]
  out_parameter_df=pd.DataFrame(data=out_parameter)

  out_file=os.path.join(in_folder, out_file1)
  print('Export selected column data of combined chromosomes at :', out_file)

  out_parameter_file=os.path.join(in_folder, out_parameter_file1)
  print('Export parameters of combined chromosomes at :', out_parameter_file)
  out_parameter_df.to_csv(out_parameter_file, index=False, header=None, sep='\t')

  selected_cols_logReg_df_sorted= sort_selected_columns(all_logReg_df)
  print(selected_cols_logReg_df_sorted.shape)
  selected_cols_logReg_df_sorted.to_csv(out_file,index=False,header=None, sep='\t')
  return out_file

def run(args):
  in_chrs=args.in_chroms_number
  in_chrs=tuple(in_chrs.split(',') )
  print('input Chromosomes: ',in_chrs)

  in_folder=args.in_file_folder
  in_sub_folder=args.in_File_subFolder
  print('input File folders: ',in_folder,in_sub_folder)

  in_Mpert=args.in_minimum_percentage_changes
  in_PVCut=args.in_Pvalue_cutoff
  in_Smoot=args.in_is_smoothed_data
  in_ModTt=args.in_moderate_ttest
  in_MPchange=args.in_low_median_high_cutoff
  in_endfile_str=args.in_file_ending_string

  #here file name shall be changed
  #file_string='*_'+ str(in_MPchange) + '_miniPercentChange_gt_'+ str(in_Mpert) +'_Pcutoff_'+str(in_PVCut)+'_isSmooth_' + str(in_Smoot) +'_isModTest_'+ str(in_ModTt) +'_*'+in_endfile_str 
  #'_*_all'
  parameter_string_list=[ str(in_MPchange) , 'miniPercentChange_gt', str(in_Mpert) ,'Pcutoff', str(in_PVCut), 'isSmooth'  , str(in_Smoot) , 'isModTest',  str(in_ModTt)]
  file_string= '*_all_mr_data' +  in_endfile_str
  print('search File strings: ',file_string)

  ac_cutoff=np.array(args.in_accuracy_cutoff_range.split(','),float)
  ac_cutoff1=ac_cutoff[0]
  ac_cutoff2=ac_cutoff[1]

  str2dmr=args.in_DMR_string_type
  logReg_proba=args.in_LogReg_proba
  other_condition= 'logReg_proba_'+ str(logReg_proba)
  main(in_chrs, in_folder, in_sub_folder, in_Mpert,in_PVCut, in_Smoot, in_ModTt, file_string,ac_cutoff1,ac_cutoff2,str2dmr,logReg_proba,other_condition,parameter_string_list)

if __name__== '__main__':

  args=my_parser(argparse.ArgumentParser('python dmr_combine_multChrs4rank.py')).parse_args()
  run(args)



 


