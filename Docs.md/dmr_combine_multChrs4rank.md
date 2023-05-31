# Differential Methylated Region Analysis Tool Documentation

dmr_analysis is a software tool for differentially Methylated Regions analysis to rank significant DMRs.



[Home](index.md) | [DMR Block Analysis](dmr_analysis_block.md) | [Combine MultiChr4Rank](dmr_combine_multChrs4rank.md) | [Selected4plot](dmr_selected4plot.md) | [map2genome](dmr_map2genome.md) | [map2chromSegment](dmr_map2chromSegment.md) | [cal2genome_percent](dmr_cal2genome_percent.md) | [cal2chromSegment_percent](dmr_cal2chromSegment_percent.md) | [percent2plot](dmr_percent2plot.md) | [combine2geneAnnot](dmr_combine2geneAnnot.md) | [exportData](dmr_exportData.md) | [gene annotation](dmr_gene_annotation.md)  

## Combine MultiChr4Rank
<p>This program is used to Combine DMR/MRs from multiple chromomoses.</p>
<strong>Required:</strong>

<ul>
  <li><code>-inChrs IN_CHROMS_NUMBER, --in_chroms_number IN_CHROMS_NUMBER</code>: a string of chromosomes that need be combined such as chr1,chr2,chr3</li>
  <li><code>-inFold IN_FILE_FOLDER, --in_file_folder IN_FILE_FOLDER</code>: a file folder of MR data in each chromosome that is exported by dmr_analysis_block</li>
</ul>
<strong>Optional, has default values:</strong>

<ul>
  <li><code>-inSFold , --in_File_subFolder</code>: a subfolder of exported MR data in each chromosome based on dmr_analysis_block, default = plots that is under each chromosome folder</li>
  <li><code>-inMPer , --in_minimum_percentage_changes</code>: parameter of minimum percentage of data points in a MR with a predefined methylation changes greater than a cutoff, default=0.0001</li>
  <li><code>-inPC , --in_Pvalue_cutoff</code>: P value cutoff for significance test, default = 0.05</li>
  <li><code>-inIsSmd , --in_is_smoothed_data</code>: is result based on raw =0, interpolated =1 , or smoothed =2 data, default= 0</li>
  <li><code>-inAC , --in_accuracy_cutoff_range</code>: range of clustering accurancy, default include all (e.g. from 0.0 to 1.0) is 0.0,1.1</li>
  <li><code>-inDMRst , --in_DMR_string_type</code>: a string used to represent predicted DMR in the file, default is D for prediced DMR</li>
  <li><code>-inLRpb , --in_LogReg_proba</code>: a probability value used by logistic Regression to select DMRs, default =0.8</li>
  <li><code>-inMT , --in_moderate_ttest</code>: 0 for standard T-test, 1 for moderate T-test, and 2 for KS-test for evaluating the test of significance, default=0</li>
  <li><code>-inLMH IN_LOW_MEDIAN_HIGH_CUTOFF, --in_low_median_high_cutoff IN_LOW_MEDIAN_HIGH_CUTOFF</code>: use high, median, or low minimum percentage change for DMR ranking, default = high</li>
  <li><code>-inFEstr IN_FILE_ENDING_STRING, --in_file_ending_string IN_FILE_ENDING_STRING</code>: file name ending string that will be used to search for results file in result folder (e.g., chrY/plots/), default = _all</li>
</ul>
