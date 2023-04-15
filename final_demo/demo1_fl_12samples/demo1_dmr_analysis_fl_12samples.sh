#!/bin/bash

#use bash script to call dmr_analysis
##Please note variable under "#-- " have to be manually adjusted in different run or in different input data! 


#path of input data folder 
#here, we assume all WGBS methylation profiles are already prepared in bed format at chromosome named folders under in_wgbs_folder,
#-- where file names indicate the sample name and conditions
in_wgbs_folder='../../final_demo_data/fl_12samples/in_data/WGBS-data/' 

#path of input reference genome folder 
#-- where genome size file and refFlat files will be used in making predefined genomic regions (e.g., TSS, TES, gene et al.) by using hmst-seq-analyzer
in_genome_folder='../../final_demo_data/genome/hg19/'
in_genome_refFlat='hg19.refFlat.txt'
in_genome_size='hg19.chrom.sizes.clear.sorted'
replace='_clean_sorted.bed'
finds='.txt'
in_sorted_refFlat=${in_genome_refFlat//$finds/$replace}

#path of input chromatin segment folder 
#such as combined chromatin segment from six human cell lines generated from Encode project https://genome.ucsc.edu/cgi-bin/hgFileUi?db=hg19&g=wgEncodeAwgSegmentation
#-- where each predicted chromatin type (e.g., TSS, Enhancers) is listed in a bed formated file and the file name indicate the predicted chromatin type.
in_chromSegment_folder='../../final_demo_data/chromSegment/hg19/'

#output data folder 
#path to output folders
#-- where predicted DMR/MRs will be exported
out_result_folder='../../final_demo_data/fl_12samples/out_data/DMR_CpG_context'

#-- Name of output folders and files that will be created in out_result_folder
out_folder4genome_map='out_map2genome'
logProb_cutoff=0.7
out_file4genome_map=gcb_vs_tumor_DMR_hyper_hypo_mix_${logProb_cutoff}.csv
out_folder4chromSegment_map='out_map2chromSegment'
out_file4chromSeg_map=gcb_vs_tumor_DMR_hyper_hypo_mix_in_chromSeg_${logProb_cutoff}.csv
#a file name that contains all ranked DMRs by combining results from all chromosomes
mr_IN_FILE='2_chroms_all_mr_data_range_dmrRanking'


#STEP 1. run dmr_analysis to predict DMRs
#a) do dmr_analysis in blocks
for in_chrom in chrY chr18
do 
dmr_analysis dmr_analysis_block --in_file_folder $in_wgbs_folder \
        --chromosome $in_chrom --group_key $in_chrom \
        --out_file_folder $out_result_folder \
        --wildType_fileString gcb \
        --data_start_position 3 --data_end_position 15 \
        --maximum_adjacency_length 1000 --minimum_block_size 5 \
        --P_cutoff 0.05 --minimum_percentage_changes 0.0001 \
        --percentage_cutoff 0.05,0.1,0.2 --low_median_high_cutoff 2 \
        --number_of_processes 15 \
        --is_smoothed_data 2 --is_moderate_ttest 0 --is_export_data 1 \
        --column_splitBy_dotOrUnderscore 1 
done
echo "dmr_analysis_block - Done"

#b) combine results from multiple chromosomes and rank the DMRs
dmr_analysis dmr_combine_multChrs4rank \
	--in_chroms_number chr18,chrY \
	--in_file_fold $out_result_folder \
	--in_is_smoothed_data 2 \
	--in_LogReg_proba 0.7 \
	--in_low_median_high_cutoff high \
	--in_file_ending_string _range.tsv 
echo dmr_combine_multChrs4rank - Done

#STEP 2. Plot and export data
chrom='chr18'
#-- please note the name of in_DMR_file may be changed in different run because of the parameters, the total number of input and the top percentage et al
  #in_DMR_file=${chrom}'_maxDist_1000_minSize_5_DMR_clusterAccuracy_gt_0.5_miniMethyChange_gt_0.05_0.1_0.2_high_miniPercentChange_gt_0.0001_Pcutoff_0.05_isSmooth_2_isModTest_0_1285_range_dmrRanking_top_0.78_minLogReg_proba_0.7'
in_DMR_file=${chrom}'_all_mr_data_range_dmrRanking.tsv'
in_data_file=${chrom}'_MR_data4maxBlockDistance_1000_minBlockSize_5_data.txt.gz'
in_wildType_string='gcb'

#some additional features for plotting and exporting data
# select DMR for ploting such as mr5,mr6,mr8 from selected chromosome
#here --in_DMR_file is exported by dmr_combine_multChrs4rank at folder "out_result_folder"/chrY/plots
##--in_data_file is exported by dmr_analysis_block at folder "out_result_folder"/chrY
dmr_analysis dmr_selected4plot --in_DMR_file ${in_DMR_file} \
        --in_data_file  ${in_data_file} \
        --in_data_folder ${out_result_folder}/${chrom}/ \
        --column_splitBy_dotOrUnderscore 1 --is_plot 1 --is_export 1 \
	--needs_check_mr mr5,mr6,mr8 --wildType_fileString ${in_wildType_string} \
	--out_folder ${out_result_folder}/out_selected4plot

echo "dmr_selected4plot -- Done"

#export selected DMR based on bed format file 0
##--input_file contains all MRs in bed foramt that need to extract their raw and smoothed methylation data
dmr_analysis dmr_exportData  \
                       --input_mr_data_folder ${out_result_folder} \
                       --output_file_folder ${out_result_folder}/out_exportData \
                       --input_file_format 0 \
                       --wildType_fileString ${in_wildType_string} --column_splitBy_dotOrUnderscore 1 --input_file test_mr.bed

echo "dmr_ExportData -- Done"



#STEP 3. mapp predicted DMR/MRs to predefined genomic regions (e.g., TSS, TES, 5dist etl al) or predicted chromatin segments for further analysis
#below is a result file generated from dmr_combine_multChrs4rank, where DMR/MRs from multiple chromosomes are combined and ranked them by logisitic regression model 
#-- Please note this file name needs to be input manually because it is generated after running "dmr_combine_multChrs4rank" and expored at folder "out_result_folder"
#mr_IN_FILE='2_chroms_high_miniPercentChange_gt_0.0001_Pcutoff_0.05_isSmooth_2_isModTest_0__range_dmrRanking_top_0.78_minLogReg_proba_0.7'

#a) generate predefined genomic regions (e.g., TSS, TES, gene et al.) by using hmst-seq-analyzer
#https://hmst-seq.github.io/hmst/
#Here, to edit exported "list_region_files.txt" for adding/removing predefined genomic regions.
#For example, to add file path for enhancer reginos in "list_region_files.txt" if user want to include enhancer in the analysis
hmst_seq_analyzer gene_annotation -F ${out_result_folder} -i no -l 10 \
        -xL 50000000 -X 5000 -Y 1000 -M 5000 -N 1000000 -hu yes -n no \
        -r ${in_genome_folder}/${in_genome_refFlat} \
        -g ${in_genome_folder}/${in_genome_size}
echo export genome annotation files at: ${out_result_folder}/data
echo gene_annotation-DON

#b) map DMR to predefined genomic regions such as TSS, TES, gene et al.
dmr_analysis dmr_map2genome --in_sortedDMR_file ${out_result_folder}/${mr_IN_FILE}.bed \
        --in_geneRegion_file ${out_result_folder}/list_region_files.txt \
        --in_outFile_folder ${out_result_folder}/${out_folder4genome_map} \
        --in_refFlat_file ${out_result_folder}/data/${in_sorted_refFlat}
echo dmr_map2genome - Done


#c) calculate percentage of DMR in annotated genomic regions
dmr_analysis dmr_cal2genome_percent --in_outFile_folder ${out_result_folder}/${out_folder4genome_map} \
        --in_outFile_name ${out_file4genome_map} --in_LogReg_proba ${logProb_cutoff} \
        --in_fileName_string $mr_IN_FILE  
echo dmr_cal2genome_percent - Done

#d) plot percentage of DMR in annotated genomic regions
dmr_analysis dmr_percent2plot --in_countFile_folder ${out_result_folder}/${out_folder4genome_map} \
        --in_countFile_name ${out_file4genome_map}
echo dmr_percent2plot - Done


#e) map DMR to predicated chromatin states such as predicated chromatin segment from 6 human cell lines.
dmr_analysis dmr_map2chromSegment --in_chromatinSegment_file_folder ${in_chromSegment_folder} \
        --in_fileName_string 'combined_six*bed.gz' --in_combined_chromatinSegment_exist 1 \
        --in_outFile_folder ${out_result_folder}/${out_folder4chromSegment_map} \
        --in_DMR_file ${out_result_folder}/${mr_IN_FILE}.bed
echo dmr_map2chromSegment - Done
  
#f) calculate percentage of DMRs in predicted chromatin states.
dmr_analysis dmr_cal2chromSegment_percent --in_outFile_folder ${out_result_folder}/${out_folder4chromSegment_map} \
        --in_outFile_name ${out_file4chromSeg_map} \
        --in_fileName_string ${mr_IN_FILE}_combined_six_cells_chromatin_segment_min10
echo dmr_cal2chromSegment_percent - Done

#g) plot percentage DMR in chromSegment
dmr_analysis dmr_percent2plot --in_countFile_folder ${out_result_folder}/${out_folder4chromSegment_map} \
        --in_countFile_name ${out_file4chromSeg_map} 
echo dmr_percent2plot - Done

#h) Combine annotated results from both genome and chromatin segment 
#please note both genome and chromatin segment have to be available before running this function.
#This function is slow and not recommend to use for large data but use dds_analysis instead.
dmr_analysis dmr_combine2geneAnnot --number_of_processes 10 --miniLogReg_proba_cutoff 0.7 \
	--sortedDMR_file ${mr_IN_FILE}.bed \
	--dmr_outFile_folder ${out_result_folder}/ \
	--dmr_genomeFile_folder ${out_folder4genome_map} \
	--dmr_chromSegmentFile_folder ${out_folder4chromSegment_map} \
	--dmr_genomeFile_string '2_chroms_*' \
	--dmr_chromSegmentFile_string '2_chroms_*'
echo dmr_combine2geneAnnot - Done




























