
# Differential Methylated Region Analysis Tool 
## dmr_analysis Documentation

dmr_analysis is a software tool for differentially Methylated Regions analysis to rank significant DMRs.



[Home](index.md) | [DMR Block Analysis](dmr_analysis_block.md) | [Combine MultiChr4Rank](dmr_combine_multChrs4rank.md) | [Selected4plot](dmr_selected4plot.md) | [map2genome](dmr_map2genome.md) | [map2chromSegment](dmr_map2chromSegment.md) | [cal2genome_percent](dmr_cal2genome_percent.md) | [cal2chromSegment_percent](dmr_cal2chromSegment_percent.md) | [percent2plot](dmr_percent2plot.md) | [combine2geneAnnot](dmr_combine2geneAnnot.md) | [exportData](dmr_exportData.md) | [gene annotation](dmr_gene_annotation.md) | [gene annotation](dmr_gene_annotation.md)


## cal2chromSegment_percent

Calculate percentage of DMRs in chromatin segmemnt (e.g., generated from ENCODE predictions of 6 human cell lines)
<p><strong>Required:</strong></p>
<ul>
  <li><code>-inOFolder IN_OUTFILE_FOLDER, --in_outFile_folder IN_OUTFILE_FOLDER</code>: input file folder for chromatin segment out file, for example, export files from "dmr_map2chromSegment"</li>
  <li><code>-inOFile IN_OUTFILE_NAME, --in_outFile_name IN_OUTFILE_NAME</code>: output file name of percentage DMR that will be calculated in chromatin segment regions</li>
  <li><code>-inFileStr IN_FILENAME_STRING, --in_fileName_string IN_FILENAME_STRING</code>: input file name string that will be searched (e.g., 24_chroms_high_miniPercentChange_gt_0.0001_Pcutoff_0.05_isSmooth_2_isModTest_0__all_dmrRanking_top_0.95_minLogReg_proba_0.7*) under folder --in_outFile_folder</li>
</ul>

<p><strong>Optional, has default values:</strong></p>
<ul>
  <li><code>-inLRpb , --in_LogReg_proba</code>: a probability value used by logistic Regression to select DMRs, default =0.8</li>
</ul>

