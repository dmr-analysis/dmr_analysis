# Differential Methylated Region Analysis Tool Documentation

dmr_analysis is a software tool for differentially Methylated Regions analysis to rank significant DMRs.



[Home](index.md) | [DMR Block Analysis](dmr_analysis_block.md) | [Combine MultiChr4Rank](dmr_combine_multChrs4rank.md) | [Selected4plot](dmr_selected4plot.md) | [map2genome](dmr_map2genome.md) | [map2chromSegment](dmr_map2chromSegment.md) | [cal2genome_percent](dmr_cal2genome_percent.md) | [cal2chromSegment_percent](dmr_cal2chromSegment_percent.md) | [percent2plot](dmr_percent2plot.md) | [combine2geneAnnot](dmr_combine2geneAnnot.md) | [exportData](dmr_exportData.md)| [gene annotation](dmr_gene_annotation.md)

## map2chromSegment
<p>Map DMR/MRs to predicted chromation segments from 6 human cell lines. </p>

### Required Parameters:
<ul>
  <li><code>-inFolder IN_CHROMATINSEGMENT_FILE_FOLDER, --in_chromatinSegment_file_folder IN_CHROMATINSEGMENT_FILE_FOLDER</code> - Path of input file folder for combined chromation segment files from six cells in bed format, combined six cells *.bed.gz files download from <a href="https://genome.ucsc.edu/cgi-bin/hgFileUi?db=hg19&amp;g=wgEncodeAwgSegmentation">https://genome.ucsc.edu/cgi-bin/hgFileUi?db=hg19&amp;g=wgEncodeAwgSegmentation</a></li>
  <li><code>-outFolder IN_OUTFILE_FOLDER, --in_outFile_folder IN_OUTFILE_FOLDER</code> - Path of output file folder for DMR/MRs mapped to chromSegments</li>
  <li><code>-inDFile IN_DMR_FILE, --in_DMR_file IN_DMR_FILE</code> - input DMR file name, which is sorted export file from DMR_analysis with BED format such as export file from "dmr_combine_multChrs4rank"</li>
  <h2> Optional Parameters </h2>
  
  <li><strong>-inMinLen IN_MINIMUM_LENGTH4REGION, --in_minimum_length4region IN_MINIMUM_LENGTH4REGION</strong> - minimum length for each selected chromatin segment, default= 10 bp</li>
  <li><strong>-inFileStr IN_FILENAME_STRING, --in_fileName_string IN_FILENAME_STRING</strong> - String name for input chromatin segment files that are stored in folder, default is *.bed.gz</li>
  <li><strong>-OFileStr IN_OUTFILENAME_STRING, --in_outFileName_string IN_OUTFILENAME_STRING</strong> - Output file type string that will be used for output combined chromatin segment files, default is "combined_six_cells_chromatin_segment_min"</li>
  <li><strong>-inMinOp IN_MINIMUM_OVERLAP4BEDTOOLS, --in_minimum_overlap4bedtools IN_MINIMUM_OVERLAP4BEDTOOLS</strong> - minimmum overlap rate in bedtools intersection, default is 1e-9 or 1bp</li>
  <li><strong>-inExist IN_COMBINED_CHROMATINSEGMENT_EXIST, --in_combined_chromatinSegment_exist IN_COMBINED_CHROMATINSEGMENT_EXIST</strong> - Whether the combined six cell chromatin segment files exists or not, default=0 for not exist which means the -in_chromatinSegment_file_folder is the path for raw files from downloaded. if inExist=1 then it indicates the combined six cell chromatin segment files are ready and sorted by chromosome position, where the path of -inFolder is the location of thes combined files of six cells with file name started with combined_six_cells_*.gz</li>
</ul>

