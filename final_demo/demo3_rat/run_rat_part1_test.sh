#!/bin/bash
#use bash script to call dmr_analysis
#this demo does not include chromatin segment data and no combination of restuls from genomic regions with chromatin segmentations.
##Please note variable under "#-- " have to be manually adjusted in different run or in different input data! 


#--IN parameters --
#path of input data folder 
#here, we assume all WGBS methylation profiles are already prepared in bed format at chromosome named folders under in_wgbs_folder,
#-- where file names indicate the sample name and conditions
in_wgbs_folder='../../final_demo_data/rat_data/in_data/WGBS-data/' 

#path of input reference genome folder 
#-- where genome size file and refFlat files will be used in making predefined genomic regions (e.g., TSS, TES, gene et al.) by using hmst-seq-analyzer
in_genome_folder='../../final_demo_data/genome/rn6/'
in_genome_refFlat='rn6.refFlat.txt'
in_genome_size='rn6.chrom.sizes.clear.sorted'
replace='_clean_sorted.bed'
finds='.txt'
in_sorted_refFlat=${in_genome_refFlat//$finds/$replace}

#--OUT parameters --
#output data folder 
#path to output folders
#-- where predicted DMR/MRs will be exported
out_result_folder='../../final_demo_data/rat_data/out_data/DMR_CpG_context'

#-- Name of output folders and files that will be created in out_result_folder
out_folder4genome_map='out_map2genome'
#cutoff value for selected DMRs in plot of map2genome
logProb_cutoff=0.7
out_file4genome_map=control_vs_test_DMR_hyper_hypo_mix_${logProb_cutoff}.csv
#a file name that contains all ranked DMRs by combining results from all chromosomes
mr_IN_FILE='*_chroms_all_mr_data_range_dmrRanking'

#STEP 1. run dmr_analysis to predict DMRs
#a) do dmr_analysis in blocks
for in_chrom in chr1 chr2 chr3 chrX chrY
do 
dmr_analysis dmr_analysis_block --in_file_folder $in_wgbs_folder \
        --chromosome $in_chrom --group_key $in_chrom \
        --out_file_folder $out_result_folder \
        --wildType_fileString _Ctrl \
        --data_start_position 3 --data_end_position 13 \
        --maximum_adjacency_length 1000 --minimum_block_size 5 \
        --P_cutoff 0.05 --minimum_percentage_changes 0.0001 \
        --percentage_cutoff 0.05,0.1,0.2 --low_median_high_cutoff 2 \
        --number_of_processes 10 \
        --is_smoothed_data 2 --is_moderate_ttest 0 --is_export_data 0 \
        --column_splitBy_dotOrUnderscore 0 
done
echo "dmr_analysis_block - Done"

#b) combine results from multiple chromosomes and rank the DMRs
dmr_analysis dmr_combine_multChrs4rank \
	--in_chroms_number chr1,chr2,chr3,chrX,chrY \
	--in_file_fold $out_result_folder \
	--in_is_smoothed_data 2 \
	--in_LogReg_proba 0.6 \
	--in_low_median_high_cutoff high \
	--in_file_ending_string _range.tsv 
echo dmr_combine_multChrs4rank - Done

#STEP 2. Plot and export data for selected DMRs
#-- please note the name of in_DMR_file may be changed in different run because of the parameters, the total number of input and the top percentage et al
chrom='chr3'
in_DMR_file=${chrom}'_all_mr_data_range_dmrRanking.tsv'
in_data_file=${chrom}'_MR_data4maxBlockDistance_1000_minBlockSize_5_data.txt.gz'
in_wildType_string='_Ctrl'

#a) some additional features for plotting and exporting data
# select DMR for ploting such as mr5,mr9,mr11
#here --in_DMR_file is exported by dmr_combine_multChrs4rank in "out_result_folder"/chrY/plots
##--in_data_file is exported by dmr_analysis_block in "out_result_folder"/chrY
dmr_analysis dmr_selected4plot --in_DMR_file ${in_DMR_file} \
	--in_data_file ${in_data_file} \
	--in_data_folder ${out_result_folder}/${chrom}/ \
	--column_splitBy_dotOrUnderscore 0 --is_plot 1 --is_export 1 \
	--needs_check_mr mr2,mr9,mr12 --wildType_fileString  ${in_wildType_string} \
	--out_folder ${out_result_folder}/out_selected4plot

echo plot selected MR - Done

#b) export selected DMR based on bed format file 0
##--input_file_name contains all MRs in bed foramt that need to extract their raw and smoothed methylation data
dmr_analysis dmr_exportData  \
                       --input_mr_data_folder ${out_result_folder} \
                       --output_file_folder ${out_result_folder}/out_exportData \
                       --input_file_format 0 \
                       --wildType_fileString ${in_wildType_string} --input_file test_mr.bed
echo export selected MR - Done

#STEP 3. mapp predicted DMR/MRs to predefined genomic regions (e.g., TSS, TES, 5dist etl al) or predicted chromatin segments for further analysis
#below is a result file generated from dmr_combine_multChrs4rank, where DMR/MRs from multiple chromosomes are combined and ranked them by logisitic regression model 
#-- Please note this file name needs to be input manually because it is generated after running "dmr_combine_multChrs4rank" and expored at "out_result_folder"
#mr_IN_FILE='5_chroms_high_miniPercentChange_gt_0.0001_Pcutoff_0.05_isSmooth_2_isModTest_0__range_dmrRanking_top_0.73_minLogReg_proba_0.6'
#mr_IN_FILE='*_chroms_all_mr_data_range_dmrRanking'

#a) generate predefined genomic regions (e.g., TSS, TES, gene et al.) by using hmst-seq-analyzer
#https://hmst-seq.github.io/hmst/
#Here, to edit exported "list_region_files.txt" for adding/removing predefined genomic regions
#For example, to add file path for enhancer reginos in "list_region_files.txt" if user want to include enhancer in the analysis
hmst_seq_analyzer gene_annotation -F ${out_result_folder} -i no -l 10 \
        -xL 50000000 -X 5000 -Y 1000 -M 5000 -N 1000000 -hu rat -n no \
        -r ${in_genome_folder}/${in_genome_refFlat} \
        -g ${in_genome_folder}/${in_genome_size}

echo export genome annotation files at: ${out_result_folder}/data
echo gene_annotation-Done
  
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
















