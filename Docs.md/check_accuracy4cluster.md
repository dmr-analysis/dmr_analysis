
# Differential Methylated Region Analysis Tool 
## dmr_analysis Documentation

dmr_analysis is a software tool for differentially Methylated Regions analysis to rank significant DMRs.



[Home](index.md) | [DMR Block Analysis](dmr_analysis_block.md) | [Combine MultiChr4Rank](dmr_combine_multChrs4rank.md) | [Selected4plot](dmr_selected4plot.md) | [map2genome](dmr_map2genome.md) | [map2chromSegment](dmr_map2chromSegment.md) | [cal2genome_percent](dmr_cal2genome_percent.md) | [cal2chromSegment_percent](dmr_cal2chromSegment_percent.md) | [percent2plot](dmr_percent2plot.md) | [combine2geneAnnot](dmr_combine2geneAnnot.md) | [exportData](dmr_exportData.md)  | [gene annotation](dmr_gene_annotation.md) 


## check_accuracy4cluster

Significant test of TF binding affinity changes between the foreground and the background calculations

## Required Paramters
<ul>
  <li><code>background_folder: </code> Folder containing the computed background model "
                                                                "produced by background_affinity_changes.py",</li>
<li><code>foreground_folder: </code> Folder containing the computed delta-dba values for "
                                                                "patient data produced by bayespi_bar.py"</li>
 
  
</ul>

## Optional Parameters

<ul>
  <li><code>output_file: </code> Output file name, default is - </li>
<li><code>p_value: </code> P-value cutoff for the Wilcoxon test, default=0.05</li>
  <li><code>max_rank: </code> Maximum rank of PWM to consider in the patient results, default=15 </li>
<li><code>exact_test: </code> Use exact Wilcoxon rank-sum test from R. R needs to be installed, default is False"</li>
  <li><code>pval_correction: </code> Whether adjust P-values by bonferroni correction, default is False"</li>
<li><code>: </code></li>
  
</ul>
