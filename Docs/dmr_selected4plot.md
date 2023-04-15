# Differential Methylated Region Analysis Tool Documentation

dmr_analysis is a software tool for differentially Methylated Regions analysis to rank significant DMRs.



[Home](index.md) | [DMR Block Analysis](dmr_analysis_block.md) | [Combine MultiChr4Rank](dmr_combine_multChrs4rank.md) | [Selected4plot](dmr_selected4plot.md) | [map2genome](dmr_map2genome.md) | [map2chromSegment](dmr_map2chromSegment.md) | [cal2genome_percent](dmr_cal2genome_percent.md) | [cal2chromSegment_percent](dmr_cal2chromSegment_percent.md) | [percent2plot](dmr_percent2plot.md) | [combine2geneAnnot](dmr_combine2geneAnnot.md) | [exportData](dmr_exportData.md)

## Selected4plot
<p>Plot figure or export methylation data for selected DMR/MRs
</p>

<ul>
  <li><code>-inDMR IN_DMR_FILE</code>: A file generated after performing "dmr_analysis_block" and "dmr_combine_multChrs4rank" which contains DMR information and logReg probability values for each MR in addtion to the main feature scores exported from dmr_analysis_block. usually, this file is exported for each chromosome under folder "plots" at the output folder of dmr_analysis_block. For example, a file "chrY_maxDist_250_minSize_5_DMR*_top_0.95_minLogReg_proba_0.7" exported in "plots" of output folder "chrY" after running "dmr_combine_multChrs4rank".</li>
  <li><code>-inData IN_DATA_FILE</code>: A file contains raw DNA methylation data of all MRs in a chromosome after running "dmr_analysis_block" this file is exported for each chromosome at the output folder after running dmr_analysis_block. For example, a file (chrY_MR_data4maxBlockDistance_250_minBlockSize_5_data.txt.gz) is stored at chrY output folder after running "dmr_analysis_block".</li>
  <li><code>-inFolder IN_DATA_FOLDER</code>: Path of a file folder that contains all files mentioned in <code>-inDMR</code> and <code>-inData</code>. For example, "out/chrY" is the path that store DMR/MR raw data for "chrY" after running dmr_analysis_block and dmr_combine_multChrs4rank.</li>
  
  <h2>Optional Paramters </h2>
  
  <li><strong>-inDMRFolder</code></strong>: A subfolder for storing exported MR/DMRs in each chromosome after running dmr_analysis_block, default is plots/ under --in_data_folder.</li>
  <li><strong>-inAcCut</strong>: A range of clustering accuracy will be considered in the program, the default is "0.0,1.0" which means between 0.0 and 1.0.</li>
  <li><strong>-inDmrSt</strong>: A string used in <code>--inDMR_file</code> that represents the predicted DMR, default is "D" representing prediced DMR.</li>
  <li><strong>-inOtCod IN_OTHER_CONDITION</strong>: A string type for selecting DMR from the whole data such as, "logReg_proba_", "min_dmr_", "max_dmr_", "top_to_down_", "top_rank_", "down_rank_" , "None". default is "down_rank_10", 10 means 10 percentage. Usually, it shall use the default parameter to select DMRs for plotting.</li>
  <li><strong>-inDstart</strong>: The postion of data start column in a file inputted from <code>--in_data_file</code>, default is 3 because the first three columns are chromosome positions such as chrom, start_post, end_post, data1, data2, ....</li>
  <li><strong>-inIsExt</strong>: Whether to export data 1 or not export data 0 for selected DMR/MRs, default is 0 for not exporting data.</li>
  <li><strong>-dotOrUnderscore</strong>: 0 for dot . split column label, 1 for underscore _ split column label in file <code>--in_data_file</code>, default=0 dot split column labels.</li>
  
  
  <li><strong>-inIsPlt, --is_plot</strong> - whetehr to plot figure "1" or not plot figure "0" for selected DMR/MRs, default is 0 not plot figure</li>
  <li><strong>-inIsPas, --is_pause</strong> - whether to pause 1 or continue 0 in loop, default is 0 for continue in the loop of all input DMR/MRs</li>
  <li><strong>-inLogNm, --in_Logger_number</strong> - logger number for recording screen displays, default is 0</li>
  <li><strong>-inOutFd, --out_folder</strong> - Path for output file folder, default is out/</li>
  <li><strong>-inNdChk, --needs_check_mr</strong> - String ID of selected DMR/MRs that will be checked/plotted/data exported which is comma separated DMR/MR IDs such as "mr10,mr111" default is not. Currently, this function only supports for plotting/exporting data for DMR/MRs within the same chromosome in each run.</li>
  <li><strong>-dpi, --figure_dpi</strong> - Export figure resolution in DPI, default dpi=30</li>
  <li><strong>-format, --figure_format</strong> - File format of exported figure such as jpg, svg, pdf, or png, default = pdf</li>
  <li><strong>-wtStr, --wildType_fileString</strong> - A file name of the wild type condition file which shall start with these characters. For example, if a file name starts with "gcb_meht1_*" is a wild type/control sample, then --wildType_fileString is gcb, which is the default setting in the program</li>
</ul>
