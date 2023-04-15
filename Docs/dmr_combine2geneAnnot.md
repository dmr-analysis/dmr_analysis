# Differential Methylated Region Analysis Tool Documentation

dmr_analysis is a software tool for differentially Methylated Regions analysis to rank significant DMRs.



[Home](index.md) | [DMR Block Analysis](dmr_analysis_block.md) | [Combine MultiChr4Rank](dmr_combine_multChrs4rank.md) | [Selected4plot](dmr_selected4plot.md) | [map2genome](dmr_map2genome.md) | [map2chromSegment](dmr_map2chromSegment.md) | [cal2genome_percent](dmr_cal2genome_percent.md) | [cal2chromSegment_percent](dmr_cal2chromSegment_percent.md) | [percent2plot](dmr_percent2plot.md) | [combine2geneAnnot](dmr_combine2geneAnnot.md) | [exportData](dmr_exportData.md)
## combine2geneAnnot

Combine results of both genome and chromSegment (e.g.,requests exported files from both dmr_map2genome and dmr_map2chromSegment) to a single gene annotation file for DMRs.

### Required Parameters
<ul>
  <li><code>in_group1_str : </code>Group1 column label </li>
<li><code>in_group2_str : </code> Group2 column label</li>
  <li><code> in_folder: </code>  Input file folder where contains bpb3 exported
                        differential expression gene file</li>
<li><code> in_file: </code>Input file name of bpb3 exported differential
                        expression gene file</li>
    
</ul>

### Optional Parameters

<ul>
  <li><code> min_median_RPKM: </code>Minimum median of RPKM value in each group,
                        default=1.0 </li>
<li><code> rr_cutoff: </code> Cutpff value for rratios to filter differential
                        expression genes beteween two groups default=0.18,
                        rr>0.14, 0.18, 0.22, 0.4, 0.67 represent 1.15, 1.2,
                        1.25, 1.5, 2.0 folder changes, respectively.</li>
  <li><code> sort_by_column_str: </code> Column name of dataframe column that will be sorted
                        after filtering ,
                        default=differential_expressoin_T_test_p_value</li>
<li><code> is_median: </code> Use mean or median of group values to calculate folder
                        changes between groups, default=False use mean value
                        of the group</li>

</ul>  

[Home](index.md)
