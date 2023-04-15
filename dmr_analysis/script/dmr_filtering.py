#this script is used to filter DMR only
#exec(open("filter_dmr.py").read())

#public package
import os
import pandas as pd
import numpy as np
import matplotlib as mlt
#mlt.use('TkAgg')
mlt.use('Agg')
import matplotlib.pyplot as plt
import argparse

#in-house package
#import dmr_data_analysis
#import dmr_analysis_block
#import dmr_utility
#import plot_dmr
from .script_high.plot_dmr import select_DMRs_for_analysis, get_parameters_for_analysis


def my_parser(parser):
   required =parser.add_argument_group('Required')
   required.add_argument('-in_DMRfile','--in_dmr_file',required=True,help='input DMR file from a chromosome, chrY_*_all')
   required.add_argument('-in_DTfile','--in_data_file',required=True, help='input DMR data file for a chromosome, chrY_*_data.txt.gz')
   #required.add_argument('-in_folder','--in_file_folder',required=True, help='input DMR and data file folder')
   required.add_argument('-in_DMRlevel','--in_selected_DMR_level',required=True,help='filtering level for DMRs, low, median, or high ')

   optional = parser.add_argument_group('Optional, has default values')
   optional.add_argument('-wtStr','--wildType_fileString',
                         default='gcb',
                         type=str,
                         metavar='',
                         help='First few character string that labeled in file name as Wide Type condition, for example, if file name start with gcb_meht1_* as wild type control sample'
                              ', then --wildType_fileString is gcb, which is the default setting in the program '
                        )
   return parser

def filtering_DMRs_by_dmr_parameters(args):
    '''
      filtering a dataframe based on a tuple of input file parameters
      return filtred dataframe that passed these parameters
    '''
    #clustering accuracy
    selected_not_passed_df0, str2level, mini_percentage_cutoff, ac_cutoff1, ac_cutoff2, \
    P_cutoff,percent_cutoff4Ttest=args

    #here the three steps filtering conditions are the same as that in function  call_do_DMR_analysis_in_block()
    #filtering by clustering accuracy
    selected_not_passed_df=selected_not_passed_df0[((selected_not_passed_df0.cluster_accuracy>ac_cutoff1) &  (selected_not_passed_df0.cluster_accuracy<=ac_cutoff2))]

    #filtering by mini percentage of changes in delt(tumor-gcb) 
    tmp_passed_df=eval('selected_not_passed_df[selected_not_passed_df.'+str2level+ '_positive_tumor_vs_gcb_percent + selected_not_passed_df.'+ str2level+ '_negative_tumor_vs_gcb_percent > mini_percentage_cutoff]')

    tmp_passed_df2=tmp_passed_df.copy()
    #filtering by euclian distance Ttest pvalues and percentage of data points passed T-test
    #tmp_passed_df2=tmp_passed_df2[ ~( (tmp_passed_df2.gcb_vs_grpsDist_pval >= P_cutoff) & (tmp_passed_df2.tumor_vs_grpsDist_pval>=P_cutoff)) | (tmp_passed_df2.percent_data_passed_ttest>percent_cutoff4Ttest)]
    #tmp_passed_df2.is_DMR='D'
    tmp_passed_df2.loc[~( (tmp_passed_df2.gcb_vs_grpsDist_pval >= P_cutoff) & (tmp_passed_df2.tumor_vs_grpsDist_pval>=P_cutoff)) | (tmp_passed_df2.percent_data_passed_ttest>percent_cutoff4Ttest), 'is_DMR']='D'

    return tmp_passed_df2

def update_DMR_type_based_on_newLevel(median_passed_df0,str2level):
    median_passed_df=median_passed_df0.copy()
    exec('median_passed_df.loc[(median_passed_df.' + str2level + '_negative_tumor_vs_gcb_percent==0) & (median_passed_df.'+ str2level + '_positive_tumor_vs_gcb_percent>0),"'"DMR_type"'"]= "'"hyper"'"')
    exec('median_passed_df.loc[(median_passed_df.' + str2level + '_negative_tumor_vs_gcb_percent>0) & (median_passed_df.'+ str2level + '_positive_tumor_vs_gcb_percent==0),"'"DMR_type"'"]= "'"hypo"'"')
    exec('median_passed_df.loc[(median_passed_df.' + str2level + '_negative_tumor_vs_gcb_percent>0) & (median_passed_df.'+ str2level + '_positive_tumor_vs_gcb_percent>0),"'"DMR_type"'"]= "'"mix"'"')
    return median_passed_df

def filtering_DMR_by_miniPercentageChange_level(in_data_file,in_dmr_file, str2miniPercentageChange_level, wildType_fileString):
   '''
     filtering DMR export file based on miniPercentageChange_level low, median, high
     which is based on file name
     all other filterings are the same as function call_do_DMR_analysis_in_block when update DMR_type (hyper, hypo, mix)
   '''
   #read data files
   print('Read data file ' + in_data_file)
   print('Read DMR file ' +in_dmr_file)
   in_data_df=pd.read_csv(in_data_file,sep='\t',compression='gzip')
   in_dmr_df=pd.read_csv(in_dmr_file,sep='\t')

   #use cluster_accuracy to select DMRs  
   ac_cutoff_low=0.0
   ac_cutoff_high=1.9
   other_condition=None
   #passed DMR in input file
   str2dmr='D'
   selected_passed_df, passed_needs_check, (len_dmr,len_d,percent_recovered,num_of_top1,percent_of_top) = select_DMRs_for_analysis(ac_cutoff_low, ac_cutoff_high, str2dmr,other_condition, in_dmr_df)
   #not passesd DMR in input file
   str2dmr='U'
   selected_not_passed_df, not_passed_need_checks,(len_dmr,len_d,percent_recovered,num_of_top1,percent_of_top) =  select_DMRs_for_analysis(ac_cutoff_low, ac_cutoff_high, str2dmr,other_condition, in_dmr_df)
   print(selected_not_passed_df.shape)

   #find all dmr paramseters from file names, which can be used in the new filtering of DMRs
   data_start_col=3
   gcb_col_idx, tumor_col_idx, percent_cutoff, low_median_high_cutoff, in_chrm, \
       max_read_length, mini_block_size, P_cutoff, isSmooth, is_modT , mini_percentage_cutoff = get_parameters_for_analysis(data_start_col, in_data_df, in_dmr_file, wildType_fileString)
   #do filtering for a new level (e.g., median) of percentage of changes between tumor and gcb groups
   str2level=str2miniPercentageChange_level
   percent_cutoff4Ttest=0
   #median_passed_df=filtering_DMRs_by_dmr_parameters((selected_not_passed_df, str2level, mini_percentage_cutoff,ac_cutoff_low,ac_cutoff_high,P_cutoff,percent_cutoff4Ttest))
   median_not_passed_df=filtering_DMRs_by_dmr_parameters((selected_not_passed_df, str2level, mini_percentage_cutoff,ac_cutoff_low,ac_cutoff_high,P_cutoff,percent_cutoff4Ttest))
   median_passed_df=filtering_DMRs_by_dmr_parameters((selected_passed_df,str2level,mini_percentage_cutoff,ac_cutoff_low,ac_cutoff_high,P_cutoff,percent_cutoff4Ttest))
   print(median_passed_df.shape)
   print(median_not_passed_df.shape)
 
   #update hyper, hypo, mix type based on new level of cutoff values
   median_passed_df=update_DMR_type_based_on_newLevel(median_passed_df,str2level)
   #selected_passed_df=update_DMR_type_based_on_newLevel(selected_passed_df,str2level)
   median_not_passed_df=update_DMR_type_based_on_newLevel(median_not_passed_df,str2level)

   #combine newly filtered and old filtered results
   #new_selected_passed_df=pd.concat([selected_passed_df,median_passed_df])
   new_selected_passed_df=pd.concat([median_passed_df,median_not_passed_df])

   #find MR IDs are not included in new_selected_passed_df
   #new_selected_not_passed_df= selected_not_passed_df.copy()
   new_selected_not_passed_df=in_dmr_df.copy()
   new_selected_not_passed_df=new_selected_not_passed_df[~new_selected_not_passed_df.mr_id.isin(new_selected_passed_df.mr_id.to_list())]
   #update parameter for new not passed df
   new_selected_not_passed_df.is_DMR='U'
   new_selected_not_passed_df=update_DMR_type_based_on_newLevel(new_selected_not_passed_df,str2level)

   new_all_df=pd.concat([new_selected_passed_df,new_selected_not_passed_df])
   new_all_df=new_all_df.sort_index()

   out_file=in_dmr_file + '_' + str2miniPercentageChange_level
   print('Export at : ' + out_file)
   new_all_df.to_csv(out_file,index=False, sep='\t', float_format='%g')
   return new_all_df

def run(args):
  in_dmr_file=args.in_dmr_file
  in_data_file=args.in_data_file
  #in_data_folder=args.in_file_folder
  str2miniPercentageChange_level=args.in_selected_DMR_level
  wildType_fileString=args.wildType_fileString

  #in_dmr_folder=os.path.join(in_data_folder,'plots')
  #in_dmr_file=os.path.join(in_dmr_folder,in_dmr_file)
  #in_data_file=os.path.join(in_data_folder,in_data_file)
  print(in_dmr_file)
  print(in_data_file)

  #str2miniPercentageChange_level='high'
  new_df=filtering_DMR_by_miniPercentageChange_level(in_data_file,in_dmr_file, str2miniPercentageChange_level,wildType_fileString)
  print(new_df[new_df.is_DMR=='D'].shape)


if __name__ == '__main__':
  args=my_parser(argparse.ArgumentParser('python dmr_filtering.py')).parse_arg()
  run(args)

  #for python3
  #np.warnings.filterwarnings('ignore', category=np.VisibleDeprecationWarning)

  #input file names
  #in_dmr_file='chrY_maxDist_250_minSize_5_DMR_clusterAccuracy_gt_0.5_miniMethyChange_gt_0.07_0.15_0.2_high_miniPercentChange_gt_0.0001_Pcutoff_0.05_isSmooth_2_isModTest_0_1756_all'
  #in_data_file='chrY_MR_data4maxBlockDistance_250_minBlockSize_5_data.txt.gz'
  #in_data_folder='out/DMR/chrY/'


  #in_dmr_file='chr1_maxDist_250_minSize_5_DMR_clusterAccuracy_gt_0.5_miniMethyChange_gt_0.07_0.15_0.2_high_miniPercentChange_gt_0.0001_Pcutoff_0.05_isSmooth_2_isModTest_0_93816_all'
  #in_data_file='chr1_MR_data4maxBlockDistance_250_minBlockSize_5_data.txt.gz'
  #in_data_folder='out/DMR/chr1/'

  #in_dmr_folder=os.path.join(in_data_folder,'plots')
  #in_dmr_file=os.path.join(in_dmr_folder,in_dmr_file)
  #in_data_file=os.path.join(in_data_folder,in_data_file)
  #print(in_dmr_file)
  #print(in_data_file)

  #str2miniPercentageChange_level='high'
  #new_df=filtering_DMR_by_miniPercentageChange_level(in_data_file,in_dmr_file, str2miniPercentageChange_level)
  #print(new_df[new_df.is_DMR=='D'].shape)



#  if False:
#    #chrY passed
#    needs_check=('mr3', 'mr12','mr36','mr38','mr45')
#    loop=0
#    len_needs_check=len(needs_check)
#
#    #plot DMR only
#    for index, row in in_dmr_df.iterrows():
#       tmp_pos=row.position.split(',')
#       tmp_chr=tmp_pos.pop(0)
#       tmp_pos=list(np.array(tmp_pos,int))
#       tmp_id=row.mr_id
#
#       #print(tmp_id)
#       loop += 1
#
#       if tmp_id in needs_check :
#          show_selected_dmr(tmp_chr, tmp_pos, in_data_df, data_start_col, gcb_col_idx, tumor_col_idx ,isSmooth,mini_percentage)
#          list_check=list(needs_check)
#          list_check.remove(tmp_id)
#          needs_check=tuple(list_check)
#          
#    print('checked total DMRs ', loop ) 
#    print('needs check DMRs ', len_needs_check)
#    print('Remaining ', list(needs_check))




