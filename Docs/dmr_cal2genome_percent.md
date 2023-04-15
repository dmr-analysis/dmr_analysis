# Differential Methylated Region Analysis Tool Documentation

dmr_analysis is a software tool for differentially Methylated Regions analysis to rank significant DMRs.



[Home](index.md) | [DMR Block Analysis](dmr_analysis_block.md) | [Combine MultiChr4Rank](dmr_combine_multChrs4rank.md) | [Selected4plot](dmr_selected4plot.md) | [map2genome](dmr_map2genome.md) | [map2chromSegment](dmr_map2chromSegment.md) | [cal2genome_percent](dmr_cal2genome_percent.md) | [cal2chromSegment_percent](dmr_cal2chromSegment_percent.md) | [percent2plot](dmr_percent2plot.md) | [combine2geneAnnot](dmr_combine2geneAnnot.md) | [exportData](dmr_exportData.md)

## cal2genome_percent

Calculate percentage of DMRs in predfined genomic regions (e.g., TSS, TES, Gene, 5Dist, et al)
<ul>
  <li>
    <strong>Required:</strong>
    <ul>
      <li>-inOFolder <code>IN_OUTFILE_FOLDER</code>, --in_outFile_folder <code>IN_OUTFILE_FOLDER</code><br>
          Path of file folder that stores the exported files from dmr_map2genome , for example, exported MR/DMRs intersected to predefined genomic regions from "dmr_map2genome"
      </li>
      <li>-inOFile <code>IN_OUTFILE_NAME</code>, --in_outFile_name <code>IN_OUTFILE_NAME</code><br>
          output file name for storing the percentage of DMRs in predifined genmoic regions, respectively.
      </li>
      <li>-inFileStr <code>IN_FILENAME_STRING</code>, --in_fileName_string <code>IN_FILENAME_STRING</code><br>
          Input file name string that will be searched by program to obtain mapped MR/DMR information For example, an output bed format file from "dmr_combine_multChrs4rank" is 24_chroms_high_miniPerc entChange_gt_0.0001_Pcutoff_0.05_isSmooth_2_isModTest_ 0__all_dmrRanking_top_0.95_minLogReg_proba_0.7*, which is usually stored in folder --in_outFile_folder. In the new version, it becomes a file such as 5_chroms_all_mr_data_range_dmrRanking_*
      </li>
    </ul>
  </li>
  <li>
    <strong>Optional, has default values:</strong>
    <ul>
      <li>-LsOrGt <strong>, --in_Ls_or_Gt_pval </strong><br>
          Use less (0) or greater (1) than Probability value cutoff to search for DMRs, default=1 use greater than Probability value as cutoff to select DMRs
      </li>
      <li>-inLRpb <strong>, --in_LogReg_proba </strong><br>
          a probability value used by logistic Regression to select DMRs, default =0.8
      </li>
    </ul>
  </li>
</ul>


[Home](index.md)
