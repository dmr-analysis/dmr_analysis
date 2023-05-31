# Differential Methylated Region Analysis Tool Documentation

dmr_analysis is a software tool for differentially Methylated Regions analysis to rank significant DMRs.



[Home](index.md) | [DMR Block Analysis](dmr_analysis_block.md) | [Combine MultiChr4Rank](dmr_combine_multChrs4rank.md) | [Selected4plot](dmr_selected4plot.md) | [map2genome](dmr_map2genome.md) | [map2chromSegment](dmr_map2chromSegment.md) | [cal2genome_percent](dmr_cal2genome_percent.md) | [cal2chromSegment_percent](dmr_cal2chromSegment_percent.md) | [percent2plot](dmr_percent2plot.md) | [combine2geneAnnot](dmr_combine2geneAnnot.md) | [exportData](dmr_exportData.md) | [gene annotation](dmr_gene_annotation.md)

## percent2plot
<p>Plot percentage of DMRs in predefined genomic regions or predicted chromatin segment regions.</p>




<h2>Required:</h2>
<ul>
  <li><code>-inCFolder IN_COUNTFILE_FOLDER, --in_countFile_folder IN_COUNTFILE_FOLDER</code> - Input path of a file folder that contains a count table of DMRs in predefined genomic regions that exported by dmr_cal2genome_percent</li>
  <li><code>-inCFname IN_COUNTFILE_NAME, --in_countFile_name IN_COUNTFILE_NAME</code> - Input file name of the count table for DMR/MRs in predefined genomic regions, for example, an exported file from dmr_cal2chromSegment_percent or dmr_cal2genome_percent</li>
</ul>

<h2>Optional, has default values:</h2>
<ul>
  <li><code>-LsOrGt , --in_Ls_or_Gt_pval</code> - Use less (0) or greater (1) than Probability value cutoff to search for DMRs, default=1 use greater than Probability value as cutoff to select DMRs</li>
  <li><code>-inLRpb , --in_LogReg_proba</code> - a probability value used by logistic Regression to select DMRs, default =0.8</li>
</ul>

