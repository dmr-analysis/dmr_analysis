#
import argparse
import sys

class Main(object):
   def __init__(self):
       parser = argparse.ArgumentParser(
             description= 'DMR-Analysis: A Differentially Methylated Region Analysis Tool',
             usage= ''' dmr_analysis <task> [<args>]

Tasks available for using:
    dmr_analysis_block 	Predict Differentially Methylated Region (DMR) from genome-wide methylation regions such as by using WGBS or other similar techniques        
    dmr_combine_multChrs4rank 	Combine predicted DMRs/MRs from multiple chromosomes then rank the DMRs by using logistic regresssion model
    dmr_selected4plot 	Plot figure and export raw/smoothed methylation data for selected DMR or MR
    dmr_map2genome 	Map all DMR/MRs to reference genome
    dmr_map2chromSegment 	Map all DMR/MRs to chromation segments generated from 6 human cell lines
    dmr_cal2genome_percent 	Calculate percentage of DMRs intersected with predefined genomic regions such as TSS, TES, 5dist et al.
    dmr_cal2chromSegment_percent 	Calculate percentage of DMRs intersected with chromatin segments generated from 6 human celles 
    dmr_percent2plot 	Plot percentage of DMRs in predefined genomic or chromatin segment regions
    dmr_combine2geneAnnot 	Combine annotations from both predefined genomic regions and chromatin segments (This function is slow and requests both genome and chromatin segment results available)
    dmr_exportData	Plot and export data for DMRs/MRs located in specific regions (e.g., DMRs/MRs intersected with mutation block or enhancer regions)
''')

#dmr_filtering      Filtering DMR by using low, median , or high minimum percentage change of methylation level

       parser.add_argument('task',help='Pipeline task to run')
       args= parser.parse_args(sys.argv[1:2])
       if not hasattr(self, args.task):
          print('****Error: Unrecognized task ****')
          parser.print_help()
          exit(1)
       getattr(self,args.task)()

   def dmr_analysis_block(self):
       from .script.dmr_analysis_block import my_parser, run
       parser = my_parser(argparse.ArgumentParser(prog='dmr_analysis dmr_analysis_block',
            description='Genome-wide prediction of DMR/MRs based on bed format input file'))
       run(parser.parse_args(sys.argv[2:]))
   
   def dmr_combine_multChrs4rank(self):
       from .script.dmr_combine_multChrs4rank import my_parser, run
       parser = my_parser(argparse.ArgumentParser(prog='dmr_analysis dmr_combine_multChrs4rank',
             description='Combine DMR/MRs from multiple chromomoses'))
       run(parser.parse_args(sys.argv[2:]))

   def dmr_selected4plot(self):
       from .script.dmr_selected4plot import my_parser, run
       parser = my_parser(argparse.ArgumentParser(prog='dmr_analysis dmr_selected4plot',
             description='Plot figure or export methylation data for selected DMR/MRs')) 
       run(parser.parse_args(sys.argv[2:]))

   def dmr_map2genome(self):
       from .script.dmr_map2genome import my_parser, run
       parser= my_parser(argparse.ArgumentParser(prog='dmr_analysis dmr_map2genome',
               description = 'Map DMR/MRs to predefined genomic regions such as TSS, TES, 5dist et al. '))
       run(parser.parse_args(sys.argv[2:]))
  
   def dmr_map2chromSegment(self):
       from .script.dmr_map2chromSegment import my_parser, run
       parser= my_parser(argparse.ArgumentParser(prog='dmr_analysis dmr_map2chromSegment',
               description= 'Map DMR/MRs to predicted chromation segments from 6 human cell lines'))
       run(parser.parse_args(sys.argv[2:]))
      
   def dmr_cal2genome_percent(self):
       from .script.dmr_cal2genome_percent import my_parser, run
       parser= my_parser(argparse.ArgumentParser(prog='dmr_analysis dmr_cal2genome_percent',
                description='Calculate percentage of DMRs in predfined genomic regions (e.g., TSS, TES, Gene, 5Dist, et al)'))
       run(parser.parse_args(sys.argv[2:]))

   def dmr_cal2chromSegment_percent(self):
       from .script.dmr_cal2chromSegment_percent import my_parser, run
       parser = my_parser(argparse.ArgumentParser(prog='dmr_analysis dmr_cal2chromSegment_percent',
                 description='Calculate percentage of DMRs in chromatin segmemnt (e.g., generated from ENCODE predictions of 6 human cell lines)')) 
       run(parser.parse_args(sys.argv[2:]))
   
   def dmr_percent2plot(self):
       from .script.dmr_percent2plot import my_parser,run
       parser= my_parser(argparse.ArgumentParser(prog='dmr_analysis dmr_percent2plot',
               description='Plot percentage of DMRs in predefined genomic regions or predicted chromatin segment regions '))
       run(parser.parse_args(sys.argv[2:]))
  
   def dmr_combine2geneAnnot(self):
       from .script.dmr_combine2geneAnnot import my_parser, run
       parser= my_parser(argparse.ArgumentParser(prog='dmr_analysis dmr_combine2geneAnnot',
               description=' Combine results of both genome and chromSegment (e.g.,requests exported files from both dmr_map2genome and dmr_map2chromSegment)'
                           ' to a single gene annotation file for DMRs '))
       run(parser.parse_args(sys.argv[2:]))
   
#   def dmr_filtering(self):
#       from .script.dmr_filtering import my_parser, run
#       parser= my_parser(argparse.ArgumentParser(prog='dmr_analysis dmr_filtering',
#                description='Filter DMR by using low, median, or high minimum percehange change of methylation level'))  
#       run(parser.parse_args(sys.argv[2:])) 

   def dmr_exportData(self):
       from .script.dmr_exportData import my_parser, run
       parser= my_parser(argparse.ArgumentParser(prog='dmr_analysis dmr_exportData',
               description='Plot figure or export data for DMRs/MRs in specific regions such as DMRs overlapping to mutation blocks or enhancer regions etc.'))
       run(parser.parse_args(sys.argv[2:]))

def main():
   Main()

if __name__ == '__main__':
   Main()
