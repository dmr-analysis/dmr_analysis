#!/bin/bash
#these scripts are used to prepare for files in dds_analysis
python find_DEG_in_tss_5dist_regions_test.py

echo "Find DMR regions are overlapping with DEG of TSS and 5distance regions - Done"
 
python preprocess_data4dds_test.py
echo "Prepare DMR files for gene-DMR target prediction - Done " 
