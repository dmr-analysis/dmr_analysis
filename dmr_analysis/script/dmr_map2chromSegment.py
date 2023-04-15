#map DMR to predicted chromatin segment from 6 cell lines
import glob
import os
import subprocess
import pandas as pd
from .script_high.map_dmr2genome import count_dmrs_not_mapped2genome
from .script_high.map_dmr2chromSegment import read_bed_file_and_export2type_bed_file, sort_bed_file_df

def my_parser(parser):
   required = parser.add_argument_group('Required')
   required.add_argument('-inFolder','--in_chromatinSegment_file_folder', required=True, help='Path of input file folder for combined chromation segment files from six cells in bed format,'
                         ' combined six cells *.bed.gz files download from https://genome.ucsc.edu/cgi-bin/hgFileUi?db=hg19&g=wgEncodeAwgSegmentation ') 
   required.add_argument('-outFolder','--in_outFile_folder', required=True, help='Path of output file folder for DMR/MRs mapped to chromSegments ' )
   required.add_argument('-inDFile','--in_DMR_file', required=True, help='input DMR file name, '
                         ' which is sorted export file from DMR_analysis with BED format such as export file from "dmr_combine_multChrs4rank')
   optional = parser.add_argument_group('Optional, has default values')
   optional.add_argument('-inMinLen','--in_minimum_length4region', default=10, type=int, help='minimum length for each selected chromatin segment ,default= 10 bp')
   optional.add_argument('-inFileStr','--in_fileName_string', default='*.bed.gz', type=str, help='String name for input chromatin segment files that are stored in folder, default is *.bed.gz')
   optional.add_argument('-OFileStr','--in_outFileName_string',default='combined_six_cells_chromatin_segment_min', type=str,help='Output file type string '
                         'that will be used for output combined chromatin segment files, default is "combined_six_cells_chromatin_segment_min"')
   optional.add_argument('-inMinOp','--in_minimum_overlap4bedtools', default=1e-9, type=float, help='minimmum overlap rate in bedtools intersection, default is 1e-9 or 1bp')
   optional.add_argument('-inExist','--in_combined_chromatinSegment_exist',default=0, type=int, help='Whether the combined six cell chromatin segment files exists or not, default=0 for not exist '
                         ' which means the -in_chromatinSegment_file_folder is the path for raw files from downloaded. '
                         ' if inExist=1 then it indicates the combined six cell chromatin segment files are ready and sorted by chromosome position,'
                         ' where the path of -inFolder is the location of thes combined files of six cells with file name started with combined_six_cells_*.gz ')
   optional.add_argument('-dmrCutoff','--dmr_min_cutoff', default= 0 , metavar='', type=float, help='minimum cutoff value for selecting DMR, default is 0 which means the cutoff value will be' \
         ' selected from accompanied parameter file of in_sortedDMR_file with >= cutoff, otherwise, it will use the input dmr_min_cutoff to select DMR by <=cutoff  ')

   return parser

def main2(in_files,out_folder,column_name,in_sorted_dmr_file, dmr_min_cutoff, min_overlap,is_greater=True):
  #combined six cells files already exist in in_files
  #there is a bug in recent bedtools but bedtools2.2 version works, now use _range instead of _all files for overlapping works for both old and new bedtools 
  methylation_file=in_sorted_dmr_file
  print(out_folder)
  merged_out_files=in_files
  record_gzip_files=[]
  processes=[]
  if len(merged_out_files)<1:
     print('No input bed file is found please check input file name or path, I stop! ', merged_out_files)
     exit(1)

  for fi in merged_out_files:
      if '.gz' in fi:
         processes.append(subprocess.Popen('gunzip -f '+ fi, shell=True)) 
         record_gzip_files.append(fi.replace('.gz',''))
      else:
         record_gzip_files.append(fi)

  for p in processes:
     p.communicate()
  if len(record_gzip_files)>0:
      merged_out_files=record_gzip_files

  record_out_files=[]
  for fil in merged_out_files:
    region_file=fil
    out_methylation_name = out_folder + '/' + os.path.basename(methylation_file)[:-4]
    out_region_name = os.path.basename(fil)[:-4]
    out = out_methylation_name + '_' + out_region_name + '_overlap' + str(min_overlap) + '.bed'
    os.system('bedtools intersect -a ' + region_file + ' -b ' + methylation_file + \
                  ' -wa -wb -f ' + str(min_overlap) + ' > ' + out)
    print(out)
    record_out_files.append(out)

  #count how many DMRs are not mapped to annotated geneomic regions
  #coumt MR or DMR in genomic files
  #coumt MR or DMR in genomic files
  #dmr_min_cutoff=0.8
  total, dict_chr=count_dmrs_not_mapped2genome(in_sorted_dmr_file,record_out_files,dmr_min_cutoff,is_greater)
  if len(record_gzip_files)>0:
     for fi in record_gzip_files:
         os.system('gzip -f ' + fi)

  return total, dict_chr


def main(type_names, out_file_string,out_file_names,min_length,out_folder,column_name,in_sorted_dmr_file, dmr_min_cutoff, min_overlap,is_greater=True):
  #combine the same type of bed file from different cell lines then sort the combined file to export in bed format
  out_file_name_df=pd.DataFrame(columns=['name'],data=out_file_names)
  sorted_out_file_name=[]
  for ty in type_names:
    out_file=out_file_string+str(min_length)+'_'+ty + '.bed'

    #combine bed files
    combined_data_df=[]
    extracted_files=out_file_name_df['name'][out_file_name_df['name'].str.contains('_'+ty+'.bed')].to_list()
    tmp_df=[]
    for fil in extracted_files:
       tmp_data_df=pd.read_csv(fil, sep='\t',header=None, compression='gzip')
       tmp_df.append(tmp_data_df)
    combined_data_df=pd.concat(tmp_df)
    columns_name=['chrs','start_pos','end_pos','type','len_of_region']
    ishuman=True
    sorted_combined_data_df=sort_bed_file_df(combined_data_df,columns_name,ishuman)

    out_file=os.path.join(out_folder,out_file)
    print(out_file)
    sorted_combined_data_df.to_csv(out_file,sep='\t',index=False, header=False)
    sorted_out_file_name.append(out_file)

  #remove exported type bed files after combining them
  for fil in out_file_names:
    os.remove(fil)

  #use bed tools to merger regions in each type
  merged_out_files=[]
  for fil in sorted_out_file_name:
     out_fil=fil.replace('.bed','_merged.bed')
     command='bedtools merge -i ' + fil +  ' -c 4 -o distinct >' + out_fil
     print(out_fil)
     os.system(command)
     os.remove(fil)
     merged_out_files.append(out_fil)

  #map DMR to chromatin segmant files
  #in_dmr_folder='out_raw_dmr/'
  #dmr_min_cutoff=0.8
  #top_percent=0.91
  #isSmooth=0
  #in_sorted_dmr_file=in_dmr_folder +'24_chroms_high_miniPercentChange_gt_0.0001_Pcutoff_0.05_isSmooth_' + str(isSmooth) \
  #                     +'_isModTest_0__all_dmrRanking_top_'+ str(top_percent)+'_minLogReg_proba_'+str(dmr_min_cutoff) +'.bed'
  methylation_file=in_sorted_dmr_file
  print(out_folder)
  #out_folder=in_dmr_folder #'out/DMR/data/chromatin_segment/'
  #min_overlap=1E-9
  #print(in_dmr_folder)
  

  record_out_files=[]
  for fil in merged_out_files:
    region_file=fil
    out_methylation_name = out_folder + '/' + os.path.basename(methylation_file)[:-4]
    out_region_name = os.path.basename(fil)[:-4]
    out = out_methylation_name + '_' + out_region_name + '_overlap' + str(min_overlap) + '.bed'
    os.system('bedtools intersect -a ' + region_file + ' -b ' + methylation_file + \
                  ' -wa -wb -f ' + str(min_overlap) + ' > ' + out)
    print(out)
    record_out_files.append(out)

  #count how many DMRs are not mapped to annotated geneomic regions
  #coumt MR or DMR in genomic files
  #coumt MR or DMR in genomic files
  #dmr_min_cutoff=0.8
  total, dict_chr=count_dmrs_not_mapped2genome(in_sorted_dmr_file,record_out_files,dmr_min_cutoff,is_greater)
  return total, dict_chr


def run(args):

  is_combined_chromatinSegment_exist=args.in_combined_chromatinSegment_exist
  in_folder=args.in_chromatinSegment_file_folder
  in_file_name=args.in_fileName_string
  out_folder=args.in_outFile_folder
  if not os.path.exists(out_folder):
     print("Create output folder " + out_folder)
     os.system("mkdir " + out_folder)

  min_length=args.in_minimum_length4region
  out_file_string=args.in_outFileName_string
   
  #in_dmr_folder=args.in_DMR_file_folder
  in_sorted_dmr_file=args.in_DMR_file
  methylation_file=in_sorted_dmr_file
  min_overlap=args.in_minimum_overlap4bedtools

  in_files=glob.glob(os.path.join(in_folder,in_file_name))
  column_name=['chrs','start_pos','end_pos','type','len','strand','starts','ends','code']


  if args.dmr_min_cutoff==0:
    #jbw here file name has to change the parameters no long shown in the file name in new version
    #tmp_str=os.path.basename(in_sorted_dmr_file )
    #tmp_str=tmp_str.split('_')
    #tmp_str=tmp_str[-1]
    #tmp_str=tmp_str.replace('.bed','')
    #tmp_str=tmp_str.replace('.gz','')
    #use default cutoff value from accompanied parameter file from dmr_analysis
    # >= cutoff value
    parameter_file=in_sorted_dmr_file.replace('.bed','_parameter.bed')
    print(parameter_file)
    parameter_df=pd.read_csv(parameter_file,sep='\t',header=None)
    tmp_str=parameter_df[0].to_list()[-1]
    print(tmp_str)
    tmp_str=tmp_str.split('_')[-1]
    dmr_min_cutoff=float(tmp_str)
    is_greater=True
  else:
    #use <= cutoff value for selecting DMRs
    dmr_min_cutoff=args.dmr_min_cutoff
    is_greater=False


  if is_combined_chromatinSegment_exist==0:
    print('Processing raw download files before map to DMR ')
    out_file_names,type_names =read_bed_file_and_export2type_bed_file(in_files,column_name, out_folder,min_length)
    main(type_names, out_file_string,out_file_names,min_length,out_folder,column_name,in_sorted_dmr_file, dmr_min_cutoff, min_overlap, is_greater)
  else:
    print('Combined six cells chromatin segment files already exist, map it to DMR now ')
    print('Skip preprocess raw chromatin segment files, ', in_files, ' but export files to ', out_folder )
    main2(in_files,out_folder,column_name,in_sorted_dmr_file, dmr_min_cutoff, min_overlap, is_greater)

if __name__== '__main__':
  args=my_parser(argparse.ArgumentParser('python dmr_map2chromSegment.py ')).parse_args()
  run(args)
