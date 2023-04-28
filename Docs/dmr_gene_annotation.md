
# Differential Methylated Region Analysis Tool 
## dmr_analysis Documentation

dmr_analysis is a software tool for differentially Methylated Regions analysis to rank significant DMRs.



[Home](index.md) | [DMR Block Analysis](dmr_analysis_block.md) | [Combine MultiChr4Rank](dmr_combine_multChrs4rank.md) | [Selected4plot](dmr_selected4plot.md) | [map2genome](dmr_map2genome.md) | [map2chromSegment](dmr_map2chromSegment.md) | [cal2genome_percent](dmr_cal2genome_percent.md) | [cal2chromSegment_percent](dmr_cal2chromSegment_percent.md) | [percent2plot](dmr_percent2plot.md) | [combine2geneAnnot](dmr_combine2geneAnnot.md) | [exportData](dmr_exportData.md) | [gene annotation](dmr_gene_annotation.md) 


## dmr_gene_annotation	

Cleans reference file and creates genomic region files (TSS, geneBody, TES, 5dist and intergenic) from the reference

<p><strong>Required:</strong></p>
<ul>
  <li><code>-r file, --referenceFile file</code>: reference file</li>
  <li><code>-g file, --genomeFile file</code>: the genome file, one column with chromosome names and one column with the length of each chromosome</li>
  <li><code>-hu y/m/r, --human y/m/r</code>: is sample from human? [yes/y] or [mouse/m] or [rat/r]. if no, it must be from either mouse/m or rat/r genome</li>
  <li><code>-n y/n, --numericalChr y/n</code>: will you operate with numerical chromosome names or not, [yes/y] or [no/n]. the same as chosen here must be as in the genomeFile, default=no</li>
  
</ul>


<p><strong>Optional, has default values:</strong></p>
<ul>
  <li><code>-rm , --removeMir</code>: remove genes with names starting with 'Mir'? [yes/y] or [no/n]. default=yes</li>
  <li><code>-X</code>: the number of upstream bp, TSS, TES, gene. default=1000</li>
  <li><code>-Y</code>: the number of downstream bp, TSS, TES, gene. default=1000</li>
  <li><code>-M</code>: the number of bp from gene start site, 5dist. default=10000</li>
  <li><code>-N</code>: the number of bp from gene start site, 5dist. default=1000000</li>
  <li><code>-i , --intergenicBTGenes</code>: intergenic regions is between gene body regions [yes/y], or between TSS and TES [no/n]?. default=yes</li>
  <li><code>-l , --minIntergenicLen</code>: minimum intergenic region distance. default=2000</li>
  <li><code>-xL , --maxIntergenicLen</code>: maximum intergenic region distance, default=10000000</li>
  <li><code>-rem , --removeShort</code>: should regions be removed from TSS, TES and 5distance regions if gene body is removed for its size being less than 0? [yes/y] or [no/n]. default=yes</li>
  <li><code>-F , --folderOut</code>: what should be the name of the out-folder. default=out</li>
</ul>

