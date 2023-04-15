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

#exec(open("calc_dmr_percent.py").read())
import glob
import os
import pandas as pd

def my_parser(parser):
  required = parser.add_argument_group('Required') 
  required.add_argument('-inOFolder','--in_outFile_folder',required=True, help='input file folder for chromatin segment out file, for example, export files from "dmr_map2chromSegment"')
  required.add_argument('-inOFile','--in_outFile_name', required=True, help='output file name of percentage DMR that will be calculated in chromatin segment regions ')
  required.add_argument('-inFileStr','--in_fileName_string', required=True, help='input file name string that will be searched '
                        '(e.g., 24_chroms_high_miniPercentChange_gt_0.0001_Pcutoff_0.05_isSmooth_2_isModTest_0__all_dmrRanking_top_0.95_minLogReg_proba_0.7*) under folder --in_outFile_folder ')

  optional= parser.add_argument_group('Optional, has default values')
  optional.add_argument('-inLRpb', '--in_LogReg_proba', default=0.8, metavar='',
                        type=float, help='a probability value used by logistic Regression to select DMRs, default =0.8')

  return parser


def main(out_file, out_dmr_files, count_data_df,logProb_cutoff):
  #tmp_file=os.path.basename(out_dmr_files[0])
  #tmp_file=tmp_file.split('minLogReg_proba_')[-1]
  #tmp_file=tmp_file.split('_')[0]
  #logReg_proba_cutoff=float(tmp_file)
  logReg_proba_cutoff=float(logProb_cutoff)

  for fil in out_dmr_files:
    #fil=out_dmr_files[0]
    print('LogReg_proba>=', str(logReg_proba_cutoff))
    print(fil)
    in_file_df=pd.read_csv(fil,header=None, sep='\t')
    in_file_df.columns=['chrs','start_pos','end_pos','genome_info','mr_chrs','mr_start_pos','mr_end_pos','mr_info','mr_logReg_proba']
  
    #count mr and dmr in the region
    #genome_name=fil.split('_')[-4].lower()
    if '_T_' in fil:
      genome_name='T'
    elif '_R_' in fil:
      genome_name='R'
    elif '_TSS_' in fil:
      genome_name='TSS'
    elif 'PF' in fil:
      genome_name='PF'
    elif '_E_' in fil:
      genome_name='E'
    elif '_WE_' in fil:
      genome_name='WE'
    elif '_CTCF_' in fil:
      genome_name='CTCF'

    uq_genome_info=in_file_df.genome_info.unique().shape
    uq_mr_info=in_file_df.mr_info.unique().shape
    print(genome_name,uq_genome_info,uq_mr_info)

    uq_genome_info_df, uq_mr_info_df, uq_dmr_info_df=[],[],[]
    uq_genome_info_df=in_file_df.drop_duplicates(subset=['genome_info'])
    uq_mr_info_df=in_file_df.drop_duplicates(subset=['mr_info'])
    uq_dmr_info_df=uq_mr_info_df.copy()
    uq_dmr_info_df=uq_dmr_info_df[uq_dmr_info_df.mr_logReg_proba>=logReg_proba_cutoff]

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

  logProb_cutoff=args.in_LogReg_proba

  out_dmr_file_path=os.path.join(out_folder,file_string+'_*.bed')
  out_dmr_files=glob.glob(out_dmr_file_path)
  print(out_dmr_file_path)
  if len(out_dmr_files)<1:
     print('No DMR file is found please check input file name and path, I stop !', out_dmr_files)
     exit(1)

  count_data_df=pd.DataFrame(columns=['Total_MR','Hyper_DMR','Hypo_DMR','Mix_DMR'],index=['CTCF','T','R','TSS','PF','E','WE'])
  main(out_file, out_dmr_files,count_data_df,logProb_cutoff)
 


if __name__=='__main__':

  args=my_parser(argparse.ArgumentParser('python dmr_cal2chromSegment_percent.py')).parse_args()
  run(args)


