import glob
import os
import pandas as pd
import time

from .script_high.combine_annota2dmrs import merge_multiple_files2df, find_genome_chromSegment4dmr, find_parallel_process_and_index,  call_parallel_genome_function

def my_parser(parser):
  required= parser.add_argument_group('Required')
  required.add_argument('-inSDFile','--sortedDMR_file',required=True, 
                         help='An exported bed format file from dmr_combine_multChrs4rank, which contains DMR/MR positions, IDs, and probability from logistic regression model. '
                              ' This file shall be sorted by chromosome position.' )

  required.add_argument('-inOFolder','--dmr_outFile_folder', required=True, 
                         help='Path of file folder that contains all results from dmr_analysis such as "DMR", "genome", "chromSegment" output files.') 
 
  optional= parser.add_argument_group('Optional, has default values')
  optional.add_argument('-inGenFolder','--dmr_genomeFile_folder', default='genome/',metavar='',type=str,
                         help='Name of a File folder that contains all output results from "dmr_map2genome" which assume under --dmr_outFile_folder , default= genome/ ')
  
  optional.add_argument('-inChSgFolder','--dmr_chromSegmentFile_folder',default='chromSegment/',metavar='',type=str,
                         help='Name of a file folder that contains all output results from "dmr_map2chromSegment" which assume under --dmr_outFile_folder, default= chromSegment/')

  optional.add_argument('-inGFString','--dmr_genomeFile_string', default='24*proba_',metavar='',type=str,
                      help='File name starts with this string in folder --dmr_genomeFile_folder will be considered in analysis, default=24*proba_ ')

  optional.add_argument('-inCFString','--dmr_chromSegmentFile_string',default='24*proba_',metavar='',type=str,
                      help='File name start with this string in folder --dmr_chromSegmentFile_folder will be considered in analysis , default=24*proba_')

  optional.add_argument('-numPro','--number_of_processes', default=10,metavar='',type=int, help='Number of parallel processes will be used to combine results , default=10')

  optional.add_argument('-inMinLen','--miniLength_of_intergenetic', default=10, metavar='',type=int,help='Minimum length of intergenetic regions will be considered, default=10')

  optional.add_argument('-inMinMeg','--miniLength_of_merge', default=10, metavar='',type=int, help='Minimum length will be merged by "bedtools merger", default=10')

  optional.add_argument('-inMinLogRegP','--miniLogReg_proba_cutoff', default=0.8,metavar='',type=float,help='minimym cutoff value of logReg_proba that will be used to select DMR for gene annotation, default = 0.8')

  return parser

def main(num_of_processes, logReg_proba_min_cutoff,in_dmr_sorted_df, all_in_genome_df, all_in_chromSegment_df, in_dmr_sorted_files):
  #do parallel in genome annotation
  selected_dmr_sorted_df0= in_dmr_sorted_df[in_dmr_sorted_df.mr_logReg_proba>=logReg_proba_min_cutoff].copy()
  print('\nData size:')
  print(selected_dmr_sorted_df0.shape, 'logReg_proba>=' +str(logReg_proba_min_cutoff) )

  #selected_dmr_sorted_df0=in_dmr_sorted_df.copy()
  selected_dmr_sorted_df0=selected_dmr_sorted_df0.reset_index(drop=True)
  #here is only for test
  #selected_dmr_sorted_df0=selected_dmr_sorted_df0.loc[0:1000,:]
  print(selected_dmr_sorted_df0.shape)

  num_of_processes,block_chunks=find_parallel_process_and_index(num_of_processes,selected_dmr_sorted_df0)
  print('\nNumber of parallel processes:')
  print(num_of_processes)

  #input('Click to start run ')
  time.sleep(10)
  all_out_df=call_parallel_genome_function(num_of_processes,selected_dmr_sorted_df0,block_chunks,all_in_genome_df,all_in_chromSegment_df, in_dmr_sorted_files,logReg_proba_min_cutoff)

def run(args):
  in_dmr_sorted_files=args.sortedDMR_file
  #jbw here file name has to change the parameters no long shown in the file name in new version
  #itmp_string=in_dmr_sorted_files.split('LogReg_proba_')[-1]
  #tmp_string=tmp_string.replace('.bed','')
  #parameter_file=os.path.join(inin_dmr_sorted_files.replace('.bed','_parameter.bed')
  #print(parameter_file)
  #parameter_df=pd.read_csv(parameter_file,sep='\t',header=None)
  #tmp_str=parameter_df[0].to_list()[-1]
  #print(tmp_str)
  #tmp_string=tmp_str.split('_')[-1]
  #print('Input files minimum LogReg probability is ',tmp_string)
  #dmr_mini_cutoff=float(tmp_string)
  
  in_out_genome_annot_folder=args.dmr_outFile_folder
  in_dmr_sorted_file= os.path.join(in_out_genome_annot_folder,in_dmr_sorted_files)
  in_dmr_sorted_df=pd.read_csv(in_dmr_sorted_file,sep='\t',header=None)
  in_dmr_sorted_df.columns= ['mr_chrs','mr_start_pos','mr_end_pos','mr_info','mr_logReg_proba']

  #jbw here file name has to change the parameters no long shown in the file name in new version
  parameter_file=in_dmr_sorted_file.replace('.bed','_parameter.bed')
  print('\nLoad parameter file:')
  print(parameter_file)
  parameter_df=pd.read_csv(parameter_file,sep='\t',header=None)
  tmp_str=parameter_df[0].to_list()[-1]
  print(tmp_str)
  tmp_string=tmp_str.split('_')[-1]
  print('Input files minimum LogReg probability is ',tmp_string)
  dmr_mini_cutoff=float(tmp_string)


  intergenetic_min_len=args.miniLength_of_intergenetic
  in_out_genome_file=args.dmr_genomeFile_string + '*_overlap'+'*.bed'
  #jbw 2024
  tmp_in_genome_folders=os.path.join(in_out_genome_annot_folder,args.dmr_genomeFile_folder,in_out_genome_file)
  tmp_in_genome_files=glob.glob(tmp_in_genome_folders)
  print(tmp_in_genome_folders)
  #tmp_in_genome_files=glob.glob(os.path.join(in_out_genome_annot_folder,args.dmr_genomeFile_folder,in_out_genome_file))

  #load all genome annotated files but remove intergenitc file
  for fi in tmp_in_genome_files:
     if 'intergenic' in fi:
       if not ('minLen'+str(intergenetic_min_len) in fi):
          tmp_in_genome_files.remove(fi)

  print('\nLoad genome files: ')
  [  print(os.path.basename(fi)) for fi in tmp_in_genome_files] 
  #read genome info files
  in_genome_lines,all_in_genome_df=merge_multiple_files2df(tmp_in_genome_files)
  all_in_genome_df.columns=['chrs','start_pos','end_pos','genome_info','mr_chrs','mr_start_pos','mr_end_pos','mr_sites','mr_logReg_proba']

  #load all chromSegment annotated files
  min_merge_len=args.miniLength_of_merge
  in_out_chromSegment_folder=os.path.join(in_out_genome_annot_folder , args.dmr_chromSegmentFile_folder)
  in_out_chromSegment_file=args.dmr_chromSegmentFile_string +'*_six_cells_'+ '*.bed'
  tmp_in_chromSegment_files=glob.glob(os.path.join(in_out_chromSegment_folder,in_out_chromSegment_file))

  print('\nLoad chromSegment files:')
  [ print(os.path.basename(fi)) for fi in tmp_in_chromSegment_files]

  #read chromSegment info files
  all_in_chromSegment_lines, all_in_chromSegment_df=merge_multiple_files2df(tmp_in_chromSegment_files)
  all_in_chromSegment_df.columns=['chrs','start_pos','end_pos','genome_info','mr_chrs','mr_start_pos','mr_end_pos','mr_sites','mr_logReg_proba']

  #do parallel in genome annotation
  num_of_processes=args.number_of_processes
  logReg_proba_min_cutoff=args.miniLogReg_proba_cutoff
  main(num_of_processes, logReg_proba_min_cutoff,in_dmr_sorted_df, all_in_genome_df, all_in_chromSegment_df, in_dmr_sorted_files)


if __name__=='__main__' :
  args=my_parser(argparse.ArgumentParser('python dmr_combine2geneAnnot.py')).parse_args()
  run(args)

