#this script is used to plot DMR only
#exec(open("plot_dmr.py").read())

#import analyze_DMR_in_multiple_samples as analyze_mlt_files
import os
import pandas as pd
import numpy as np
import matplotlib as mlt
#mlt.use('TkAgg')
mlt.use('Agg')
import matplotlib.pyplot as plt

#in-house 
#import dmr_data_analysis
#import dmr_analysis_block
#import dmr_utility
from .dmr_data_analysis import replace_nan_by_median 
from .dmr_utility import setup_logger
from .dmr_data_analysis import do_DMR_analysis_in_a_block


def show_selected_dmr(logger_n,in_chrm, tmp_chr, tmp_pos, in_data_df, \
          data_start_col, gcb_col_idx, tumor_col_idx ,isSmooth,isExport, dotOrUnderscore, isPlot, \
          percent_cutoff,P_cutoff, is_modT, is_pause, low_median_high_cutoff,tmp_id,\
          wildType_fileString,in_resolution_dpi,in_figure_format,out_folder ) :
    #look for data
    selected_data_df=[]
    if len(tmp_pos)>2:
       selected_data_df=in_data_df[(in_data_df.Chrs==tmp_chr) & (in_data_df['Starts'].isin(tmp_pos)) ]
       tmp_num_of_data=len(tmp_pos)
    else:
       selected_data_df=in_data_df[(in_data_df.Chrs==tmp_chr) & (in_data_df['Starts'].between(tmp_pos[0],tmp_pos[1])) ]
       tmp_num_of_data=selected_data_df.shape[0]

    selected_all_data_df=selected_data_df.iloc[:,data_start_col:]
    block_methylation_columns=selected_all_data_df.columns.values

    selected_pos2,selected_all_data,selected_data2tumor, selected_data2gcb=[],[],[],[]

    selected_all_data=selected_all_data_df.to_numpy()

    selected_data2tumor=selected_all_data[:,tumor_col_idx[0]]
    selected_data2tumor= replace_nan_by_median(selected_data2tumor,tmp_id)

    selected_pos2=selected_data_df.iloc[:,1]

    selected_data2gcb=selected_all_data[:,gcb_col_idx[0]]
    selected_data2gcb= replace_nan_by_median(selected_data2gcb,tmp_id)

    #DMR length and number of data points in a DMR
    tmp_len=max(tmp_pos)-min(tmp_pos)
    #tmp_num_of_data=len(tmp_pos)
 

    fig=plt.figure(figsize=(18,8))
    #is use smoothed data for analysis
    #isSmooth=False
    if isSmooth==2:
       logger_n.info('analysis of smoothed data, P cutoff= %g ', P_cutoff)
    elif isSmooth==1:
       logger_n.info('analysis of interpolated data, P cutoff= %g ', P_cutoff)
    else:
       logger_n.info('analysis of raw data, P cutoff= %g ', P_cutoff)
    
    #np.warnings.filterwarnings('ignore', category=np.VisibleDeprecationWarning)
    #there is a warning bug about narray at this function??
    statistic_test_p, cluster_accuracy, negative_percent2, positive_percent2,results2, total_percent2 = do_DMR_analysis_in_a_block(logger_n,percent_cutoff, selected_data2tumor,
                selected_pos2, selected_data2gcb,tmp_len, tmp_id, tmp_num_of_data, gcb_col_idx, \
                tumor_col_idx, block_methylation_columns,fig,isSmooth,isExport, dotOrUnderscore, P_cutoff,is_modT,low_median_high_cutoff, wildType_fileString,out_folder,tmp_chr) 

    logger_n.info(statistic_test_p)
    logger_n.info('Start to plot ')
    out_fig_name=in_chrm+'_DMR_'+ tmp_id + '.' + in_figure_format
    if not os.path.exists(out_folder):
       os.mkdir(out_folder)
       print('Create, ', out_folder)
    out_fig_folder=os.path.join(out_folder, tmp_chr)
    if not os.path.exists(out_fig_folder):
       os.mkdir(out_fig_folder)
       print('Create, ', out_fig_folder)
    out_fig_name=os.path.join(out_fig_folder,out_fig_name)
    logger_n.info(out_fig_name)
    if is_pause:
       plt.show()
    elif isPlot :
       #plt.show(block=False)
       #input('')
       logger_n.info('Export figure resolution ' + str(in_resolution_dpi) + ' dpi')
       plt.savefig(out_fig_name,dpi=in_resolution_dpi)
       #plt.pause(2)
    plt.close()

def compute_D_and_U_percent_in_selection(other_condition,is_not_passed_df, str2type, str2dmr,ac_cutoff1, ac_cutoff2):
   '''select rows based on input condtion
   '''
   #if 'top_to_down' in other_condition:
   percent_of_top=int(other_condition.split('_')[-1])/100    
   num_of_top1=int(np.ceil(is_not_passed_df.shape[0]*percent_of_top))
   
   #sort data by logReg_score
   is_not_passed_df=is_not_passed_df.sort_values(by='logReg_score')
   is_not_passed_df=is_not_passed_df.reset_index(drop=True)

   print('Overall prediction')
   len_dmr=np.where(is_not_passed_df.iloc[num_of_top1:].is_DMR=='D')[0].shape 
   len_d=np.where(is_not_passed_df.is_DMR=='D')[0].shape
   percent_recovered= len_dmr[0]/len_d[0]
   print(len_dmr,len_d,percent_recovered)

   in_not_passed_df2=is_not_passed_df.copy()
   if 'top_to_down' in str2type or 'top_rank' in str2type:
       num_of_top1=-1*num_of_top1
       in_not_passed_df2=in_not_passed_df2[num_of_top1:] 
   elif 'down_rank' in str2type:
       in_not_passed_df2=in_not_passed_df2[:num_of_top1]

   print('Selected prediction ')
   print(str2type,num_of_top1,percent_of_top)
   sl_len_dmr=np.where(in_not_passed_df2.is_DMR=='D')[0].shape 
   sl_len_d=np.where(in_not_passed_df2.is_DMR=='U')[0].shape
   sl_percent_recovered= [sl_len_dmr[0]/abs(num_of_top1), sl_len_d[0]/abs(num_of_top1)]
   print('D','U','D_percent','U_percent', ' in selection') 
   print(sl_len_dmr,sl_len_d,sl_percent_recovered)

   if str2dmr==None:
       in_not_passed_df2=in_not_passed_df2[ ((in_not_passed_df2.cluster_accuracy>ac_cutoff1) & (in_not_passed_df2.cluster_accuracy<=ac_cutoff2)) ]  
   else:
       in_not_passed_df2=in_not_passed_df2[ (in_not_passed_df2.is_DMR==str2dmr)  & ((in_not_passed_df2.cluster_accuracy>ac_cutoff1) & (in_not_passed_df2.cluster_accuracy<=ac_cutoff2)) ]
    
   return in_not_passed_df2, (len_dmr,len_d,percent_recovered, num_of_top1, percent_of_top)

def select_DMRs_for_analysis(ac_cutoff1, ac_cutoff2, str2dmr, other_condition, in_dmr_df):
  '''selet DMRs from in_dmr_df dataframe based on range of clustering accurancy [ac_cutoff1, ac_cutoff2] and DMR type str2dmr[D/U]
     return a tuple of DMRs are selected for further analysis -> needs_check
  '''
  #use cluster_accuracy to select DMRs for plot 
  #ac_cutoff1=0.9
  #ac_cutoff2=1.9
  #str2dmr='D'
  is_not_passed_df=in_dmr_df.copy()
  if str2dmr != None:
     in_not_passed_df=is_not_passed_df[is_not_passed_df.is_DMR==str2dmr]
  else:
     in_not_passed_df=is_not_passed_df
  
  #if other_condition==None:
  in_not_passed_df=in_not_passed_df[((in_not_passed_df.cluster_accuracy>ac_cutoff1) &  (in_not_passed_df.cluster_accuracy<=ac_cutoff2))]
  print('Cluster accuracy range :',ac_cutoff1, ac_cutoff2)
  len_dmr,len_d, percent_recovered, num_of_top1, percent_of_top =[0],[0],0,0,0
  print(str2dmr)
  print(other_condition)

  if other_condition !=None and 'logReg_proba' in other_condition:
     #selection based on logReg probability greater than a defined value and cluster accuracy range
     proba=float(other_condition.split('_')[-1])
     in_not_passed_df2=in_not_passed_df[((in_not_passed_df.cluster_accuracy>ac_cutoff1) &  (in_not_passed_df.cluster_accuracy<=ac_cutoff2)) & (in_not_passed_df.logReg_predicted_dmr>=proba) ]
     if str2dmr != None:
        len_dmr=np.where(in_not_passed_df2.is_DMR==str2dmr)[0].shape
     else:
        len_dmr=np.where(in_not_passed_df2.is_DMR=='D')[0].shape
     all_len_d=np.where(is_not_passed_df.is_DMR=='D')[0].shape
     len_d=np.where(in_not_passed_df2.is_DMR=='D')[0].shape
     percent_recovered=len_dmr[0]/all_len_d[0]
     print('Number of selected DMR, ','Recovered number of DMR, ','Total number of DMR, ','Percentage of recovered DMR')
     print(len_dmr, len_d, all_len_d, percent_recovered)

  elif other_condition=='min_dmr':
     #selection based on minimum logReg_score with D class plus cluster accuracy
     is_not_passed_df=is_not_passed_df.sort_values(by='logReg_score')
     is_not_passed_df=is_not_passed_df.reset_index()
     mini_logReg_score= is_not_passed_df.loc[is_not_passed_df[is_not_passed_df.is_DMR=='D' ].index.min(),].logReg_score
     mini_logReg_id= is_not_passed_df.loc[is_not_passed_df[is_not_passed_df.is_DMR=='D' ].index.min(),].mr_id
     in_not_passed_df2=in_not_passed_df[((in_not_passed_df.cluster_accuracy>ac_cutoff1) &  (in_not_passed_df.cluster_accuracy<=ac_cutoff2)) & (in_not_passed_df.logReg_score<mini_logReg_score) ]
  elif other_condition !=None and 'max_dmr' in other_condition:
     #selection based on cluster accuracy and logReg_score greater than a prefined value such as a minimumm logReg_score 
     mini_logReg_score=float(other_condition.split('_')[-1]) 
     in_not_passed_df2=in_not_passed_df[((in_not_passed_df.cluster_accuracy>ac_cutoff1) &  (in_not_passed_df.cluster_accuracy<=ac_cutoff2)) & (in_not_passed_df.logReg_score>=mini_logReg_score) ]
  elif other_condition != None and 'top_to_down' in other_condition :
     str2type='top_to_down'
     in_not_passed_df2,(len_dmr,len_d,percent_recovered,num_of_top1, percent_of_top)= compute_D_and_U_percent_in_selection(other_condition,is_not_passed_df, str2type, str2dmr,ac_cutoff1,ac_cutoff2)
  elif other_condition !=None and 'top_rank' in other_condition :
     str2type='top_rank'
     in_not_passed_df2,(len_dmr,len_d,percent_recovered,num_of_top1,percent_of_top)= compute_D_and_U_percent_in_selection(other_condition,is_not_passed_df, str2type, str2dmr, ac_cutoff1,ac_cutoff2)
  elif other_condition !=None and 'down_rank' in other_condition :
     str2type='down_rank'
     in_not_passed_df2,(len_dmr,len_d,percent_recovered, num_of_top1,percent_of_top)= compute_D_and_U_percent_in_selection(other_condition,is_not_passed_df, str2type, str2dmr,ac_cutoff1,ac_cutoff2)
  elif other_condition==None:
     #only filtering by the clustering accuracy range
     in_not_passed_df2=in_not_passed_df.copy()

  print(in_not_passed_df2.shape)

  #all_names=in_not_passed_df2.name.str.split(':')
  #needs_check=[ i[-1].replace(' ','') for i in all_names ]
  #needs_check=tuple(needs_check)
  all_names=in_not_passed_df2.mr_id
  needs_check=tuple(all_names.tolist())

  return in_not_passed_df2, needs_check, (len_dmr[0],len_d[0],percent_recovered, num_of_top1,percent_of_top)

def get_parameters_for_analysis(data_start_col, in_data_df, in_dmr_file,wildType_fileString):
  '''
    Obtain parameters from in_data_df dataframe and in_dmr_file name for the further analysis
    return:  gcb_col_idx, tumor_col_idx, percent_cutoff, low_median_high_cutoff, in_chrm, max_read_length, mini_block_size, P_cutoff, isSmooth, is_modT
  '''
  if in_data_df.empty :
     tmp_data_cols=[]
     gcb_col_idx=[]
     tumor_col_idx=[]
  else:
     tmp_data_cols=in_data_df.columns.to_list()
     gcb_col_idx=np.where(np.char.find(tmp_data_cols[data_start_col:],wildType_fileString)>=0)
     tumor_col_idx=np.where(np.char.find(tmp_data_cols[data_start_col:],wildType_fileString)<0)

  #get input parameters from the file name
  #here file name has to be changed jbw 2023 because in the new version  , 
  #these parameters are stored in a file such as chr11_1_mr_parameters.tsv instead of at the file name!!
  t_in_dmr_file=os.path.basename(in_dmr_file)
  all_paramt=t_in_dmr_file.split('_')
  percent_cutoff=[float(i) for i in all_paramt[11:14]]
  percent_str={'low':0,'median':1,'high':2}
  low_median_high_cutoff=percent_str[all_paramt[14]]
  in_chrm=os.path.basename(t_in_dmr_file).split('_')[0]
  max_read_length=int(all_paramt[2])
  mini_block_size=int(all_paramt[4])
  mini_percentage_cutoff=float(all_paramt[17])
  P_cutoff=float(all_paramt[19])
  isSmooth=int(all_paramt[21])
  is_modT=int(all_paramt[23])

  return gcb_col_idx, tumor_col_idx, percent_cutoff, low_median_high_cutoff, in_chrm, max_read_length, mini_block_size, P_cutoff, isSmooth, is_modT, mini_percentage_cutoff

if __name__ == '__main__':

  #for python3
  #np.warnings.filterwarnings('ignore', category=np.VisibleDeprecationWarning)
  
  #in_dmr_file='chrY_maxDist_250_minSize_5_DMR_clusterAccuracy_gt_0.5_miniMethyChange_gt_0.07_0.15_0.2_high_miniPercentChange_gt_0.0001_Pcutoff_0.05_isSmooth_2_isModTest_0_1756_all_test'
  #in_data_file='chrY_MR_data4maxBlockDistance_250_minBlockSize_5_data.txt.gz'
  #in_data_folder='out/DMR/chrY/'

  #in_dmr_file='chr1_maxDist_250_minSize_5_DMR_clusterAccuracy_gt_0.5_miniMethyChange_gt_0.07_0.15_0.2_high_miniPercentChange_gt_0.0001_Pcutoff_0.05_isSmooth_2_isModTest_0_93816_all_dmrRanking_top_0.3_minLogReg_score_30.8'
  #in_data_file='chr1_MR_data4maxBlockDistance_250_minBlockSize_5_data.txt.gz'
  #in_data_folder='out/DMR/chr1/'

  in_dmr_file='chr3_maxDist_250_minSize_5_DMR_clusterAccuracy_gt_0.5_miniMethyChange_gt_0.07_0.15_0.2_high_miniPercentChange_gt_0.0001_Pcutoff_0.05_isSmooth_2_isModTest_0_65052_all_dmrRanking_top_0.94_minLogReg_proba_0.8'
  in_data_file='chr3_MR_data4maxBlockDistance_250_minBlockSize_5_data.txt.gz'
  in_data_folder='out/DMR/chr3/'
 
  #in_dmr_file='chr18_maxDist_250_minSize_5_DMR_clusterAccuracy_gt_0.5_miniMethyChange_gt_0.07_0.15_0.2_high_miniPercentChange_gt_0.0001_Pcutoff_0.05_isSmooth_2_isModTest_0_28038_all_dmrRanking_top_0.94_minLogReg_proba_0.8'
  #in_data_file='chr18_MR_data4maxBlockDistance_250_minBlockSize_5_data.txt.gz'
  #in_data_folder='out/DMR/chr18/'


  #in_data_file='chr2_MR_data4maxBlockDistance_250_minBlockSize_5_data.txt.gz'
  #in_data_folder='out/DMR/chr2/'

  in_dmr_folder=os.path.join(in_data_folder,'plots')

  in_dmr_file=os.path.join(in_dmr_folder,in_dmr_file)
  in_data_file=os.path.join(in_data_folder,in_data_file)
  print(in_dmr_file)
  print(in_data_file)

  in_data_df=pd.read_csv(in_data_file,sep='\t',compression='gzip')
  in_dmr_df=pd.read_csv(in_dmr_file,sep='\t')
   
  #use cluster_accuracy to select DMRs for plot 
  ac_cutoff1=0.6
  ac_cutoff2=1.55
  str2dmr='D'
  other_condition='down_rank_10' # 'top_rank_10' #'max_dmr_34.4'    #'top_rank_20' #'down_rank_20', 'top_to_down_10' #'max_dmr' #'mini_dmr' #0
  in_not_passed_df2, needs_check, (len_dmr,len_d,percent_recovered,num_of_top1,percent_of_top) = select_DMRs_for_analysis(ac_cutoff1, ac_cutoff2, str2dmr, other_condition, in_dmr_df)

  #find tumor and gcb data column index
  data_start_col=3
  gcb_col_idx, tumor_col_idx, percent_cutoff, low_median_high_cutoff, in_chrm, \
      max_read_length, mini_block_size, P_cutoff, isSmooth, is_modT, mini_percentage_cutoff = get_parameters_for_analysis(data_start_col, in_data_df, in_dmr_file)

  #set parameter for show DMRs
  isExport=0
  is_pause=False
  log_number=0
  out_folder='out'
  logger_n=setup_logger(out_folder, in_chrm,max_read_length,mini_block_size,log_number)

  logger_n.info(in_dmr_file)
  logger_n.info(in_data_file)
  logger_n.info(in_not_passed_df2.shape)
  str2dmr_type={'D': 'DMRs','U': 'not DMRs'}
  logger_n.info('Type of MRs is : ' +  str2dmr_type[str2dmr])
  logger_n.info('Cluster accuracy from %g to %g', ac_cutoff1, ac_cutoff2)

  #loop in each row of dataframe to plot figure
  loop=0
  num_of_find=0
  len_needs_check=len(needs_check)

  #here is for debugging
  #chr1
  #needs_check=('mr74704','mr85400')
  #chry
  #needs_check=('mr800','mr1372')

  #chr3
  #bcl6
  needs_check=('mr60384','mr60385','mr60386','mr60387')

  #chr18
  #bcl2
  #needs_check=('mr22255','mr22253','mr22254')

  #chr2
  #needs_check=('mr17880','mr17881')
  #ńeeds_check=('mr183')
  #logger_n.info(needs_check)
  #logger_n.info(low_median_high_cutoff)
  ##logger_n.info( mini_percentage_cutoff)
  #logger_n.info(percent_cutoff)
  #low_median_high_cutoff=0

  #plot DMR only
  logger_n.info(needs_check)
  input('Is continue ? ')
  for index, row in in_dmr_df.iterrows():
     tmp_pos=row.position.split(',')
     tmp_chr=tmp_pos.pop(0)
     tmp_pos=list(np.array(tmp_pos,int))
     tmp_id=row.mr_id

     #print(tmp_id)
     loop += 1

     if tmp_id in needs_check :
        logger_n.info(tmp_id)
        #start debug
        #tmp_chr='chr2'
        #in_chrm='chr2'
        #print(tmp_chr, in_chrm, tmp_pos, in_data_df, data_start_col, gcb_col_idx, tumor_col_idx ,isSmooth, isExport, percent_cutoff,P_cutoff, is_modT, is_pause, low_median_high_cutoff)
        #end of debug

        show_selected_dmr(logger_n,tmp_chr, in_chrm, tmp_pos, in_data_df, data_start_col,\
                gcb_col_idx, tumor_col_idx ,isSmooth, isExport,dotOrUnderscore, percent_cutoff,P_cutoff,\
                 is_modT, is_pause, low_median_high_cutoff,tmp_id) 
        list_check=list(needs_check)
        list_check.remove(tmp_id)
        needs_check=tuple(list_check)
        num_of_find +=1
        print(num_of_find)
 
  logger_n.info('checked total DMRs %g ', loop ) 
  logger_n.info('needs check DMRs %g ', len_needs_check)
  logger_n.info('Remaining %s ', list(needs_check))




