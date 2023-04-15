#bar plot of percentage of DMR in different genomic regions
import pandas as pd
import os
import numpy as np
import matplotlib as mlt
import matplotlib.pyplot as plt
mlt.use('Agg')
#mlt.use('TkAgg')

def my_parser(parser):
  required = parser.add_argument_group('Required')
  required.add_argument('-inCFolder','--in_countFile_folder',required=True, help='Input path of a file folder that contains a count table of DMRs in predefined genomic regions '
                            ' that exported by dmr_cal2genome_percent ')
  required.add_argument('-inCFname','--in_countFile_name', required=True, help='Input file name of the count table for DMR/MRs in predefined genomic regions, '
                       'for example, an exported file from dmr_cal2chromSegment_percent or dmr_cal2genome_percent ')
  return parser

def main(count_file_folder, count_file_name):
  count_file=os.path.join(count_file_folder,count_file_name)
  if not os.path.exists(count_file):
     print('Input count table file is not found please check input file name or path, I stop !', count_file)
     exit(1)

  count_data_df=pd.read_csv(count_file,index_col=0,sep='\t')
  count_data_matrix=count_data_df.to_numpy()
  percent_hyper=count_data_matrix[:,1]/count_data_matrix[:,0]*100
  percent_hypo=count_data_matrix[:,2]/count_data_matrix[:,0] *100
  percent_mix=count_data_matrix[:,3]/count_data_matrix[:,0] *100
  data_size=count_data_matrix[:,1].shape[0]

  percent_data_matrix=np.concatenate((percent_hyper.reshape(1,data_size),percent_hypo.reshape(1,data_size),percent_mix.reshape(1,data_size)),axis=0)
  percent_data_df=pd.DataFrame(data=percent_data_matrix.T,columns=['Hyper_percent','Hypo_percent','Mix_percent'],index=count_data_df.index.to_list())
  percent_data_df['genome']=percent_data_df.index.to_list()
  ax=percent_data_df.plot.bar(x='genome',y=['Hyper_percent','Hypo_percent','Mix_percent'], title=' tumor vs normal',figsize=(6,8))
  ax.set_ylabel('Percentage')
  out_fig_file=count_file.replace('.csv','.pdf')
  #out_fig_file=out_fig_file.replace('counts','percent')
  print(out_fig_file)
  plt.savefig(out_fig_file, format='pdf')
  #plt.show()

def run(args):
  count_file_folder=args.in_countFile_folder
  count_file_name= args.in_countFile_name
  main(count_file_folder, count_file_name)


if __name__=='__main__':

  #count_file_folder='out/DMR/'
  #logReg_proba_cutoff=0.8
  #count_file_name='fl_tumor_vs_gcb_counts_DMR_hyper_hypo_mix_'+str(logReg_proba_cutoff)+ '.csv'


  #chromatin map
  #count_file_folder='out_raw_dmr/genomes/'
#  count_file_folder='out_raw_dmr/chromSegment/'
#  logReg_proba_cutoff=0.8
  #count_file_name='fl_tumor_vs_gcb_counts_DMR_hyper_hypo_mix_'+str(logReg_proba_cutoff)+ '.csv'
#  count_file_name='fl_tumor_vs_gcb_counts_DMR_hyper_hypo_mix_in_chromSeg_'+str(logReg_proba_cutoff)+ '.csv'

#count_file=os.path.join(count_file_folder,count_file_name)
#count_data_df=pd.read_csv(count_file,index_col=0,sep='\t')
##count_data_matrix=count_data_df.to_numpy()
#percent_hyper=count_data_matrix[:,1]/count_data_matrix[:,0]*100
#percent_hypo=count_data_matrix[:,2]/count_data_matrix[:,0] *100
#percent_mix=count_data_matrix[:,3]/count_data_matrix[:,0] *100
#data_size=count_data_matrix[:,1].shape[0]

#percent_data_matrix=np.concatenate((percent_hyper.reshape(1,data_size),percent_hypo.reshape(1,data_size),percent_mix.reshape(1,data_size)),axis=0)
#percent_data_df=pd.DataFrame(data=percent_data_matrix.T,columns=['Hyper_percent','Hypo_percent','Mix_percent'],index=count_data_df.index.to_list())
#percent_data_df['genome']=percent_data_df.index.to_list()
#percent_data_df.plot.bar(x='genome',y=['Hyper_percent','Hypo_percent','Mix_percent'], title='FL tumor vs GCB normal',figsize=(6,8))

#out_fig_file=count_file.replace('.csv','.pdf')
#out_fig_file=out_fig_file.replace('counts','percent')
#print(out_fig_file)
#plt.savefig(out_fig_file, format='pdf')
#plt.show()
  args=my_parser(argparse.ArgumentParser('python dmr_percen2plott.py')).parse_args()
  run(args)
 





