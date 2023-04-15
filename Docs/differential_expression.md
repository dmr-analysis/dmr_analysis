# BayesPI-BAR in Python3 - bpb3 Documentation

bpb3 is a software tool for Bayesian method for protein-DNA interaction with binding affinity Ranking in Python3.


 
[Home](index.md) | [Differential Expression](differential_expression.md) | [Gene Regions](gene_regions.md) | [Mussd](mussd.md) | [Highly Mutated Blocks](highly_mutated_blocks.md) | [BayesPiBar](bayespi_bar.md) | [Choose Background Parameters](choose_background_parameters.md) | [Background Affinity Changes](background_affinity_changes.md) | [Affinity Change Significance](affinity_change_significance_test.md) | [Parallel](parallel.md) | [make_cluster4pwm](make_cluster4pwm.md) | [bpb3 SelectedPWM](bpb3selectedPWM.md)  | [Run Pipeline](run_pipeline.md) | [clean_tmp](clean_tmp.md)  


## differential_expression
<p>This program determines which genes are differentially expressed based on RNA-seq data for two groups of samples. RPKM values arecomputed for each sample, optionally normalized, and Kolmogorov-Smirnov test/T-test is then applied to them to determine significant difference between distributions of values of the two groups. Order of optional normalizations: 
    
  </p>
<ol> 
  <li>quantile normalization</li> 
  <li>log transform</li> 
  <li>z-score transform. </li> 
</ol>

### Required Parameters:
<ul>
  <li><code>group_1_count_files: </code>list of files with read counts for group 1 (e.g., tumor group). Each file must have at least two columns: gene name and count. </li>
  <li><code>group_2_count_files: </code>List of files with read counts for group 2 (e.g., normal group).Each file must have at least two columns: gene name and count.</li>
  <li><code>gene_lengths: </code>File with two columns: gene name and its length. Only genes listed in this file will be considered in the computation. Note that RPKM computation is affected by the set of considered genes".</li>
  
</ul>  
[Home](index.md)


## Optional Parameters

<ul>    
    
  <li><code>output_group_1_rpkm:  </code> Output file name for gene expression in group 1. Two columns will be written: gene name and the corresponding median RPKM (quantile-normalized if --quantile_normalization is specified across group 1 counts. default is None.</li>
  <li><code>output_group_2_rpkm: </code> Output file name for gene expression in group 2. Two columns will be written: gene name and the corresponding median RPKM (quantile-normalized if --quantile_normalization is specified ) across group 2 counts. default is None</li>
  <li><code>p_value: </code> (float) P-value cutoff for the KS test or T-test, default= 0.05.</li>
  <li><code>quantile_normalization: </code>Do not rr. This mode can handle partially computed results from an interrupted computation. Only the missing or corrupted output will be recomputed. default is False.</li>
  <li><code>integration: </code> (bool) Apply quantile normalization to RPKM values of all experiments, default is False. </li>
  <li><code>log_transform:  </code>(bool) Apply log-transformation (ln(1 + x)) to all values, default is False</li>
  <li><code>z_score:  </code> Apply z-score transformation values of each experiment separately, default is False</li>
  <li><code>min_fold_change:  </code>In addition to the KS-test or T-test, check that the minimum fold change in median RPKM between the two groups is above the specified number. If quantile normalization is activated, the quantile normalized values are compared, default is None</li>
  <li><code>min_medianRPKM: </code>(float) Check the median of RPKM in each group, and only keep genes with RPKM greater than the minimum value in bith groups.If quantile normalization is activated, then the quantile normalized values are checked, default is None" </li>
  <li><code>output_all_values: (float) </code>Put values (RPKM or their z-scores) for all input datasets as additional columns in the output file, default is False."</li>
  <li><code>test_method:  </code>(int) Differential expression test methods: 0 for KS-test, 1 for T-test, default is 0 for KS-test</li>
  <li><code>output_file:  </code>Output file name. Two columns will be written: gene name and the corresponding KS test P-value. Only genes with P-values smaller than the threshold given by --p_value parameter will be written, default is -.</li>

</ul>




    

   


    

   


