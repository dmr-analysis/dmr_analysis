#this function is to plot selected DMR/MR 
import os
import pandas as pd
import numpy as np
import matplotlib as mlt
import glob
#mlt.use('TkAgg')
mlt.use('Agg')
import matplotlib.pyplot as plt

#in-house function
from .script_high.plot_dmr import show_selected_dmr, compute_D_and_U_percent_in_selection, select_DMRs_for_analysis, get_parameters_for_analysis
from .script_high.dmr_utility import setup_logger

def my_parser(parser):
   required =parser.add_argument_group('Required')
   required.add_argument('-inDMR','--in_DMR_file',required=True,help='A file generated after performing "dmr_analysis_block" and "dmr_combine_multChrs4rank" ' 
                       ' which contains DMR information and logReg probability values for each MR in addtion to the main feature scores exported from dmr_analysis_block. '
                       ' usually, this file is exported for each chromosome under folder "plots" at the output folder of dmr_analysis_block. '
                       ' For example, a file "chrY_maxDist_250_minSize_5_DMR*_top_0.95_minLogReg_proba_0.7" exported in "plots" of output folder "chrY" after running "dmr_combine_multChrs4rank".   ')

   required.add_argument('-inData','--in_data_file', required=True, help='A file contains raw DNA methylation data of all MRs in a chromosome after running "dmr_analysis_block" '
                         ' this file is exported for each chromosome at the output folder after running  dmr_analysis_block.'
                         ' For example, a file (chrY_MR_data4maxBlockDistance_250_minBlockSize_5_data.txt.gz) is stored at chrY output folder after running  "dmr_analysis_block".')

   required.add_argument('-inFolder','--in_data_folder',required=True, help='Path of a file folder that contains all files mentioned in -inDMR and -inData. '
                         ' For example, "out/chrY" is the path that store DMR/MR raw data for "chrY" after running dmr_analysis_block and dmr_combine_multChrs4rank. ')
  
   optional= parser.add_argument_group('Optional , has default values')
   optional.add_argument('-inDMRFolder','--in_DMR_folder', default='plots',metavar='',
                           help='A subfolder for storing exported MR/DMRs in each chromosome after running dmr_analysis_block, '
                                                          ' default is plots/ under --in_data_folder' )
   optional.add_argument('-inAcCut','--in_Accuracy_cutoff_range', default='0.0,1.1', metavar='',type=str, 
                           help='A range of clustering accuracy will be considered in the program, the default is "0.0,1.0" which means between 0.0 and 1.0')

   optional.add_argument('-inDmrSt','--in_DMR_string_type',default='D', metavar='',type=str, 
                            help='A string used in --in_DMR_file that represents the predicted DMR, default is "D" representing prediced DMR')

   optional.add_argument('-inOtCod','--in_Other_condition', default='down_rank_10',type=str,
                           help= ' A string type for selecting DMR from the whole data such as, "logReg_proba_", "min_dmr_", "max_dmr_", '
                                 ' "top_to_down_", "top_rank_", "down_rank_" , "None". default is "down_rank_10", 10 means 10 percentage. Usually, '
                                 ' it shall use the default parameter to select DMRs for plotting.')

   optional.add_argument('-inDstart','--in_Data_start_column', default=3, metavar='', type=int,
                        help=' The postion of data start column in a file inputted from --in_data_file, default is 3 because the first three columns are chromosome positions '
                             ' such as chrom, start_post, end_post, data1, data2, ....')

   optional.add_argument('-inIsExt','--is_export', default=0, metavar='', type=int, help='whether to export data 1 or not export data 0 for selected DMR/MRs, default is 0 for not exporting data')

   optional.add_argument('-dotOrUnderscore','--column_splitBy_dotOrUnderscore', default=0,
                          type=int, metavar='', help= '0 for dot . split column label, 1 for underscore _ split column label in file --in_data_file , default=0 dot split column labels')

   optional.add_argument('-inIsPlt','--is_plot', default=0, metavar='', type=int, help='whetehr to plot figure "1" or not plot figure "0" for selected DMR/MRs, default is 0 not plot figure')

   optional.add_argument('-inIsPas','--is_pause', default=0, metavar='', type=int, help='whether to pause 1 or continue  0 in loop, default is 0 for continue in the loop of all input DMR/MRs')

   optional.add_argument('-inLogNm','--in_Logger_number', default=0, metavar='', type=int, help='logger number for recording screen displays , default is 0')

   optional.add_argument('-inOutFd','--out_folder', default='out', metavar='', type=str, help='Path for output file folder, default is out/')

   optional.add_argument('-inNdChk','--needs_check_mr', default='not', metavar='', type=str,
                          help='String ID of selected DMR/MRs that will be checked/plotted/data exported which is comma separated DMR/MR IDs such as "mr10,mr111" default is not.'
                               'Currently, this function only supports for plotting/exporting data for DMR/MRs  within the same chromosome in each run. ')

   optional.add_argument('-dpi','--figure_dpi', default=30, metavar='', type=int, help='Export figure resolution in DPI, default dpi=30')

   optional.add_argument('-format','--figure_format',default='pdf',metavar='',type=str, 
                           help='File format of exported figure such as jpg, svg, pdf, or png, default = pdf')

   optional.add_argument('-wtStr','--wildType_fileString',
                         default='gcb',
                         type=str,
                         metavar='',
                         help='A file name of the wild type condition file which shall start with these characters. For example, if a file name starts with "gcb_meht1_*"  is a wild type/control sample'
                              ', then --wildType_fileString is gcb, which is the default setting in the program '
                        )
   return parser

def main(logger_n, needs_check, in_dmr_df, in_chrm,in_data_df,\
        data_start_col,gcb_col_idx, tumor_col_idx, isSmooth, isExport, dotOrUnderscore, isPlot, percent_cutoff, \
        P_cutoff, is_modT, is_pause, low_median_high_cutoff, wildType_fileString ,\
        in_resolution_dpi,in_figure_format,out_folder):
   #loop in each row of dataframe to plot figure
   loop=0
   num_of_find=0
   len_needs_check=len(needs_check)

   #plot DMR only
   logger_n.info(needs_check)
   #input('Click any key to continue ? ')
   for index, row in in_dmr_df.iterrows():
     tmp_pos=row.position.split(',')
     tmp_chr=tmp_pos.pop(0)
     if '-' in tmp_pos[0]:
       tmp_pos=list(np.array(tmp_pos[0].split('-'),int))
     else:
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
                gcb_col_idx, tumor_col_idx ,isSmooth, isExport,dotOrUnderscore,isPlot, percent_cutoff,P_cutoff,\
                 is_modT, is_pause, low_median_high_cutoff,tmp_id, wildType_fileString, \
                 in_resolution_dpi,in_figure_format,out_folder)
        list_check=list(needs_check)
        list_check.remove(tmp_id)
        needs_check=tuple(list_check)
        num_of_find +=1
        print(num_of_find)

   logger_n.info('checked total DMRs %g ', loop )
   logger_n.info('needs check DMRs %g ', len_needs_check)
   logger_n.info('Remaining %s ', list(needs_check))


def run(args):
   in_resolution_dpi=args.figure_dpi
   in_figure_format=args.figure_format
   in_dmr_file=args.in_DMR_file
   in_data_file=args.in_data_file 
   in_data_folder=args.in_data_folder
   in_dmr_folder=args.in_DMR_folder
   in_dmr_folder=os.path.join(in_data_folder,in_dmr_folder)
   wildType_fileString=args.wildType_fileString
   dotOrUnderscore= args.column_splitBy_dotOrUnderscore

   in_dmr_file=os.path.join(in_dmr_folder,in_dmr_file)
   in_data_file=os.path.join(in_data_folder,in_data_file)
   print(in_dmr_file)
   print(in_data_file)

   in_data_df=pd.read_csv(in_data_file,sep='\t',compression='gzip')
   in_dmr_df=pd.read_csv(in_dmr_file,sep='\t')
  
   #here has to find paramteter file for the in_dmr_file in each chromosome jbw 2023
   in_parameter_file=glob.glob(os.path.join(in_dmr_folder,'*_mr_parameters.tsv'))[0]
   print(in_parameter_file)
   in_parameter_df=pd.read_csv(in_parameter_file,sep='\t',header=None)
   in_parameter_file_str=in_parameter_df[0].to_list()[-1]
   print(in_parameter_file_str)

   #use cluster_accuracy to select DMRs for plot 
   ac_cutoff=np.array(args.in_Accuracy_cutoff_range.split(','),float)
   ac_cutoff1=ac_cutoff[0]
   ac_cutoff2=ac_cutoff[1]
   str2dmr=args.in_DMR_string_type
   other_condition=args.in_Other_condition

   data_start_col=args.in_Data_start_column
   #isExport=args.is_export
   if args.is_export==0:
      isExport=False
   else:
      isExport=True

   if args.is_plot ==0:
      isPlot=False
   else:
      isPlot=True

   if args.is_pause ==0:
    is_pause=False
   else:
    is_pause=True
   
   log_number=args.in_Logger_number
   out_folder=args.out_folder

 
   #find tumor and gcb data column index ?? here file name need be changed jbw
   #in the old version in_dmr_file name is the parameter, but in the new version, the parameter has to be obtained from a parameter file  
   #gcb_col_idx, tumor_col_idx, percent_cutoff, low_median_high_cutoff, in_chrm, \
   #   max_read_length, mini_block_size, P_cutoff, isSmooth, is_modT, mini_percentage_cutoff = get_parameters_for_analysis(data_start_col, in_data_df, in_dmr_file, wildType_fileString)
   gcb_col_idx, tumor_col_idx, percent_cutoff, low_median_high_cutoff, in_chrm, \
      max_read_length, mini_block_size, P_cutoff, isSmooth, is_modT, mini_percentage_cutoff = get_parameters_for_analysis(data_start_col, in_data_df, in_parameter_file_str, wildType_fileString)


   #set parameter for show DMRs
   logger_n=setup_logger(out_folder, in_chrm,max_read_length,mini_block_size,log_number)

   logger_n.info(in_dmr_file)
   logger_n.info(in_data_file)

   if args.needs_check_mr != 'not':
     needs_check=tuple(args.needs_check_mr.split(','))
   else:
     in_not_passed_df2, needs_check, (len_dmr,len_d,percent_recovered,num_of_top1,percent_of_top) = select_DMRs_for_analysis(ac_cutoff1, ac_cutoff2, str2dmr, other_condition, in_dmr_df)
     logger_n.info(in_not_passed_df2.shape)

   str2dmr_type={'D': 'DMRs','U': 'not DMRs'}
   logger_n.info('Type of MRs is : ' +  str2dmr_type[str2dmr])
   logger_n.info('Cluster accuracy from %g to %g', ac_cutoff1, ac_cutoff2)

   #loop in each row of dataframe to plot figure
   main(logger_n, needs_check, in_dmr_df, in_chrm,in_data_df, \
          data_start_col,gcb_col_idx, tumor_col_idx, isSmooth, isExport, dotOrUnderscore,isPlot, \
          percent_cutoff, P_cutoff, is_modT, is_pause, low_median_high_cutoff,wildType_fileString,in_resolution_dpi,in_figure_format,out_folder )

if __name__== '__main__':

   args=my_parser(argparse.ArgumentParser('python dmr_selected4plot.py')).parse_args()
   run(args)








