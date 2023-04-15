#map DMR to genome reference file
import os
import pandas as pd

from .script_high.map_dmr2genome import count_dmrs_not_mapped2genome
from .script_high.dmr_utility import prepare_result_folder


def my_parser(parser):
   required = parser.add_argument_group('Required')
   required.add_argument('-inSDFile','--in_sortedDMR_file',required=True, help='An exported file from DMR_analysis with BED format, such as exported file from "dmr_combine_multChrs4rank" ' 
                          ' where columns are chrom, start_position, end_position, DMR_ID, probability values. Pleaes note this input file shall be sorted by chromosome position ' 
                          ' before it is being intersected with predefined genomic regions such as TSS, gene, or TES et al.')
   required.add_argument('-inGRFile','--in_geneRegion_file', required=True, help='A list of gene region files that is generated by hmst_seq_analyszer such as '
                        ' "list_region_files.txt" where file name and path of genomic region files are listed. These predefined genomic region files will be intersected with sorted bed file of MR/DMR, '
                        ' please refer to demo of "run_hmsq.sh" for how to generate these predefined genomic region/annotation files. ')
   required.add_argument('-inOFolder','--in_outFile_folder', required=True, help='Path of output file folder for exporting intersected MR/DMR files with each genomic regions. ')
   required.add_argument('-inRefFile','--in_refFlat_file', required=True, help='Path of input refFlat reference file, such as an export file "hg19.refFlat_clean_sorted.bed" from  "hmst_seq_Analyzer" ')

   optional= parser.add_argument_group('Optional , has default values')
   optional.add_argument('-minOlp','--in_minimum_overlap4bedtools', default=1e-9, metavar='', type=float, help='minimum overlap rate in bedtools intersection, default is 1e-9 or 1bp overlap, which is float number between 0 and 1.')
   optional.add_argument('-dmrCutoff','--dmr_min_cutoff', default= 0 , metavar='', type=float, help='minimum cutoff value for selecting DMR, default is 0 which means the >= dmr_min_cutoff value will be' \
	 ' selected from accompanied parameter file of in_sortedDMR_file, otherwise, it will use the input <= dmr_min_cutoff   ')

   return parser


def main(region_files, methylation_file, reference_file, in_sorted_dmr_file, dmr_min_cutoff,out_folder,min_overlap,is_greater=True):
  #for each region to find its DMR or MRs
  #there is a bug in recent bedtools but bedtools2.2 version works, now use *_range instead of *_all file for regions overlapping, which works in both old and new version of bedtools
  record_out_files=[]
  for fil in region_files:
    region_file=fil
    region_name=os.path.basename(fil).split('_')[0].lower()
    print(region_name)
    out_methylation_name = out_folder + '/' + os.path.basename(methylation_file)[:-4]
    out_region_name = os.path.basename(region_file)[:-4]

    #print(methylation_file)
    #print(region_file)
    dist5_methylation_file = out_methylation_name + '_' + 'noGenes.bed'
    print(dist5_methylation_file)

    #print(out_methylation_name,out_region_name)
    if region_name == '5dist':
        # For 5distance we first remove genes(TSS, geneBody and TES) from the two methylation files
        cmd='bedtools intersect -a ' + methylation_file + ' -b ' + reference_file + \
                      ' -v > ' + dist5_methylation_file

        results=os.system(cmd)

        if results !=0:
           print('Error in 5dist bedtools, ', cmd )
           exit(1)

        out = dist5_methylation_file[:-4] + '_' + out_region_name + '_overlap' + str(min_overlap) + '.bed'

        cmd='bedtools intersect -a ' + region_file + ' -b ' + \
                  dist5_methylation_file + ' -wa -wb -f ' + str(min_overlap) + ' > ' + out
        results=os.system(cmd)
        if results !=0:
           print('Error in bedtools intersect, ', cmd)
           exit(1)

        os.system('rm ' + dist5_methylation_file)  # removes temporary file
    else:
        out = out_methylation_name + '_' + out_region_name + '_overlap' + str(min_overlap) + '.bed' 
        cmd= 'bedtools intersect -a ' + region_file + ' -b ' + methylation_file + \
                  ' -wa -wb -f ' + str(min_overlap) + ' > ' + out
        results=os.system(cmd)
        if results != 0:
           print('Error in bedtools intersect, ', cmd )
           exit(1)


    print(out)
    record_out_files.append(out)

  #count how many DMRs are not mapped to annotated geneomic regions
  #coumt MR or DMR in genomic files
  #dmr_min_cutoff=0.8
  if is_greater:
      #default set use >= min_cutoff for selecting DMRs
      total, dict_chr=count_dmrs_not_mapped2genome(in_sorted_dmr_file,record_out_files,dmr_min_cutoff)
  else:
      #optional one use <= min_cutoff for selecting DMRs
      total, dict_chr=count_dmrs_not_mapped2genome(in_sorted_dmr_file,record_out_files,dmr_min_cutoff,is_greater)
  return total, dict_chr
 
def run(args):
  in_sorted_dmr_file=args.in_sortedDMR_file
  in_region_files=args.in_geneRegion_file
  reference_file=args.in_refFlat_file
  min_overlap=args.in_minimum_overlap4bedtools
  out_folder=args.in_outFile_folder

  prepare_result_folder(out_folder)

  region_files=pd.read_csv(in_region_files,header=None)
  region_files=region_files.loc[:,0].to_list()
  methylation_file=in_sorted_dmr_file

  if args.dmr_min_cutoff ==0 :
    #jbw here file name has to change the parameters no long shown in the file name in new version
    #tmp_str=os.path.basename(in_sorted_dmr_file )
    #tmp_str=tmp_str.split('_')
    #tmp_str=tmp_str[-1]
    #tmp_str=tmp_str.replace('.bed','')
    #use default cutoff value from accompanied parameter file from dmr_analysis
    # >= cutoff value
    parameter_file=in_sorted_dmr_file.replace('.bed','_parameter.bed')
    print(parameter_file)
    parameter_df=pd.read_csv(parameter_file,sep='\t',header=None)
    tmp_str=parameter_df[0].to_list()[-1]
    print(tmp_str)
    tmp_str=tmp_str.split('_')[-1]
    dmr_min_cutoff=float(tmp_str)
    main(region_files, methylation_file, reference_file, in_sorted_dmr_file, dmr_min_cutoff,out_folder,min_overlap)
  else:
    #use <= cutoff value for selecting DMRs
    dmr_min_cutoff=args.dmr_min_cutoff
    is_greater=False
    main(region_files, methylation_file, reference_file, in_sorted_dmr_file, dmr_min_cutoff,out_folder,min_overlap,is_greater)

if __name__ == '__main__':
  args=my_parser(argparse.ArgumentParser('python dmr_map2genome.py ')).parse_args()
  run(args)





