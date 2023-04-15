#this function is used to calculate percentage of DMR in selected genomic regions.
#calculate percentage of DMR in 5 genomic regions
#first count the total number of MR in each type of genomic regions
#Then count the DMR of hyper, hypo, and mix type in the region, respectively
#finally calculate the percentage of DMR in each genomic reigons
#such as       
#  Total_MR        Hyper_DMR       Hypo_DMR
#GM12878.hg38.chr17_vs_H1.hg38.chr17_5mC_geneBody_imputedWith_zeros_DMRs_all     2041    221     761
#GM12878.hg38.chr17_vs_H1.hg38.chr17_5mC_5dist_imputedWith_zeros_DMRs_all        2329    141     1142
#GM12878.hg38.chr17_vs_H1.hg38.chr17_5mC_TSS_imputedWith_zeros_DMRs_all  1820    115     384
#GM12878.hg38.chr17_vs_H1.hg38.chr17_5mC_TES_imputedWith_zeros_DMRs_all  1542    186     531
#GM12878.hg38.chr17_vs_H1.hg38.chr17_5mC_intergenic_imputedWith_zeros_DMRs_all   125     9       46


import glob
import os
import pandas as pd


def my_parser(parser):
  required = parser.add_argument_group('Required') 
  required.add_argument('-inOFolder','--in_outFile_folder',required=True, help='Path of file folder that stores the exported files from dmr_map2genome ' 
                        ', for example, exported MR/DMRs intersected to predefined genomic regions from "dmr_map2genome"')
  required.add_argument('-inOFile','--in_outFile_name', required=True, help='output file name for storing the percentage of DMRs in predifined genmoic regions, respectively. ')
  required.add_argument('-inFileStr','--in_fileName_string', required=True, help='Input file name string that will be searched by program to obtain mapped MR/DMR information '
                        'For example, an output bed format file from "dmr_combine_multChrs4rank" is ' 
                        ' 24_chroms_high_miniPercentChange_gt_0.0001_Pcutoff_0.05_isSmooth_2_isModTest_0__all_dmrRanking_top_0.95_minLogReg_proba_0.7*, '
                        ' which is usually stored in folder --in_outFile_folder. In the new version, it becomes a file such as 5_chroms_all_mr_data_range_dmrRanking_*') 
 
  optional= parser.add_argument_group('Optional, has default values')
  optional.add_argument('-LsOrGt','--in_Ls_or_Gt_pval', default=1, metavar='', 
                         type=int, help='Use less (0) or greater (1) than Probability value cutoff to search for DMRs, default=1 use greater than Probability value as cutoff to select DMRs ' ) 
  optional.add_argument('-inLRpb', '--in_LogReg_proba', default=0.8, metavar='',
			type=float, help='a probability value used by logistic Regression to select DMRs, default =0.8')
  return parser


def main(out_file, out_dmr_files,count_data_df, is_less_or_great_pval,logProb_cutoff):
  #in the new version , we have to find paramter file to get cutoff values
  #but not the file name as old version added jbw 2023
  #print(out_dmr_files)
  #tmp_file=os.path.basename(out_dmr_files[0])
  #tmp_file=tmp_file.split('minLogReg_proba_')[-1]
  #tmp_file=tmp_file.split('_')[0]
  logReg_proba_cutoff=float(logProb_cutoff)

  for fil in out_dmr_files:
    #fil=out_dmr_files[0]
    if is_less_or_great_pval==1:
      print('LogReg_proba>=', str(logReg_proba_cutoff))
    else:
      print('LogReg_proba<', str(logReg_proba_cutoff))

    print(fil)
    in_file_df=pd.read_csv(fil,header=None, sep='\t')
    in_file_df.columns=['chrs','start_pos','end_pos','genome_info','mr_chrs','mr_start_pos','mr_end_pos','mr_info','mr_logReg_proba']

    #count mr and dmr in the region
    #genome_name=fil.split('_')[-4].lower()
    print(fil)
    if '_gene_u' in fil.lower():
      genome_name='gene'
    elif '_tss_u' in fil.lower():
      genome_name='tss'
    elif '_tes_u' in fil.lower():
      genome_name='tes'
    elif '_5dist_u' in fil.lower():
      genome_name='5distUp'
    elif '_intergenic_u' in fil.lower():
      genome_name='intergenic'
    elif '_enhancer_' in fil.lower() or '_enhancers_' in fil.lower():
      genome_name='enhancer'
    elif '_5dist_d' in fil.lower():
      genome_name='5distDown'

    #input('Click ')
    uq_genome_info=in_file_df.genome_info.unique().shape
    uq_mr_info=in_file_df.mr_info.unique().shape

    uq_genome_info_df, uq_mr_info_df, uq_dmr_info_df=[],[],[]
    uq_genome_info_df=in_file_df.drop_duplicates(subset=['genome_info'])
    uq_mr_info_df=in_file_df.drop_duplicates(subset=['mr_info'])
    uq_dmr_info_df=uq_mr_info_df.copy()
    if is_less_or_great_pval==1:
       uq_dmr_info_df=uq_dmr_info_df[uq_dmr_info_df.mr_logReg_proba>=logReg_proba_cutoff]
    else:
       uq_dmr_info_df=uq_dmr_info_df[uq_dmr_info_df.mr_logReg_proba<logReg_proba_cutoff]

    total_mr=uq_mr_info_df.shape[0]
    total_dmr=uq_dmr_info_df.shape[0]
    total_hyper=uq_dmr_info_df[uq_dmr_info_df.mr_info.str.contains('hyper')].shape[0]
    total_hypo=uq_dmr_info_df[uq_dmr_info_df.mr_info.str.contains('hypo')].shape[0]
    total_mix=uq_dmr_info_df[uq_dmr_info_df.mr_info.str.contains('mix')].shape[0]

    count_data_df.at[genome_name,'Total_MR']=total_mr
    count_data_df.at[genome_name,'Hyper_DMR']=total_hyper
    count_data_df.at[genome_name,'Hypo_DMR']=total_hypo
    count_data_df.at[genome_name,'Mix_DMR']=total_mix

  print('Export at: ', out_file)
  count_data_df.to_csv(out_file,sep='\t')

def run(args):
  out_folder=args.in_outFile_folder
  file_string=args.in_fileName_string
  out_file_name=args.in_outFile_name
  out_file=os.path.join(out_folder,out_file_name)

  is_less_or_great_pval=args.in_Ls_or_Gt_pval
  logProb_cutoff=args.in_LogReg_proba

  out_dmr_file_path=os.path.join(out_folder,file_string+'_*.bed')
  out_dmr_files=glob.glob(out_dmr_file_path)
  print(out_dmr_file_path)

  count_data_df=pd.DataFrame(columns=['Total_MR','Hyper_DMR','Hypo_DMR','Mix_DMR'],index=['gene','tss','tes','5dist','intergenic','enhancer'])
  main(out_file, out_dmr_files,count_data_df, is_less_or_great_pval,logProb_cutoff)

if __name__=='__main__':

  args=my_parser(argparse.ArgumentParser('python dmr_cal2genome_percent.py')).parse_args()
  run(args)

 

 


