#!/bin/bash
#this script is used to prepare files for dds_analysis

OUT_PATH='../../../final_demo_data/rat_data/out_data/DMR_CpG_context/'
FILE_FOLD=${OUT_PATH}/out4mr_not_in_tss_enhancer
BACK_FILE=${OUT_PATH}/background_samples_list.tsv


a=10
b=10
if [ $a == $b ];
then
  ##15. export dmr regions for test
  dmr_analysis dmr_exportData \
        --input_mr_data_folder '../../../final_demo_data/rat_data/out_data/DMR_CpG_context/' \
        --output_file_folder ${OUT_PATH}/out4dmr_in_deg_tss_5dist \
        --input_file_format 0 \
        --number_of_processes 15 --input_file ${OUT_PATH}/uqdmr_regions_in_deg_tss_5dist_rat.bed -wtStr '_Ctrl'

  #export mr not in tss and enhancer for test
  dmr_analysis dmr_exportData  \
        --input_mr_data_folder '../../../final_demo_data/rat_data/out_data/DMR_CpG_context/' \
        --output_file_folder ${OUT_PATH}/out4mr_not_in_tss_enhancer \
        --input_file_format 0 \
        --number_of_processes 15 --input_file ${OUT_PATH}/mr_regions_not_in_enhancers_tss.bed -wtStr '_Ctrl'
fi

#create background file
if ! [ -f $BACK_FILE ];
then
  echo $BACK_FILE " not exists and create one ! "
  if [ -e $FILE_FOLD ];
  then
    	ls  ./${FILE_FOLD}/chr*/data/*raw*.* > $BACK_FILE
        echo "Create " $BACK_FILE 
  else
    	echo "Cannot create background file because no data folder find! " $FILE_FOLD
  fi
fi

#run dTarget in dds_analysis
gene_mr_file=${OUT_PATH}/uqGeneDmr_regions_in_deg_tss_rat.bed
gene_exp_file='../../../final_demo_data/rat_data/in_data/DEG/Adrenal1vsAdrenal2_DEG_genes_zscores_tab.tsv'
in_mr_data_folder=${OUT_PATH}/out4dmr_in_deg_tss_5dist
in_background_mr_file=$BACK_FILE

a=10
b=10
number_of_samples=100
if [ $a == $b ];
then
 #for TSS
 dds_analysis dTarget_methy_vs_express -inGeneMRfile $gene_mr_file  -mrTAB  \
	-inGeneEXPfile $gene_exp_file -expTAB \
	-inMRfolder $in_mr_data_folder -outName 'tss_region_' \
	-pathDepth 1 -inBackgroundList $in_background_mr_file -cutoff 0.05 -totalSamples $number_of_samples -numOfprocesses 20

 echo "Done with TSS target gene prediction"

 #for 5dist
 gene_mr_file=${OUT_PATH}/uqGeneDmr_regions_in_deg_5dist_rat_overlap_enhancer.bed
 dds_analysis dTarget_methy_vs_express -inGeneMRfile $gene_mr_file -mrTAB  \
	-inGeneEXPfile $gene_exp_file -expTAB \
	-inMRfolder $in_mr_data_folder -outName 'distance_region_'  \
	-pathDepth 1 -inBackgroundList $in_background_mr_file -cutoff 0.01 -totalSamples $number_of_samples -numOfprocesses 20

 echo "Done with 5distance target gene prediction"
fi



