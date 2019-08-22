## Summary

This is a collection of scripts for the transcriptomic analysis of different human brain regions in type 2 diabetes (T2D). In this study, we observed prominent difference among transcriptomic profiles of diverse brain regions in T2D. We identified caudate as the brain region most involved in T2D and discovered a caudate-specific gene module  by network methods. Gene expressions of this module were disturbed in T2D and it was enriched with T2D SNPs. Hub genes of these module have also been reported to disregulated in both T2D and neurodegenerative disease. We believe this work might provide help for further studies. 

## Citation

A manuscript detailing our work has been accepted by Aging and will be online soon. 

## Code comments

| Name                        | Comments                                                     |
| --------------------------- | ------------------------------------------------------------ |
| 1_pre_sam_info.R            | Filter the samples                                           |
| 2_tissue_data.sh            | Extract counts for each brain regions                        |
| 3_matchit_optimal_outlier.R | Match the samples using optimal methods                      |
| 4_DE                        | 1) Identified differentially expressed genes using DESeq2,; <br />2) transformed raw counts to r-log values<br />3) plot the distributions of DAGs |
| 5_Samplesize                | plot the sample size and sample distribution                 |
| 6_WGCNA                     | 1) linear regression to regress out uninterested factors<br />2) heatmap<br />3) WGCNA network construction<br />4) WGCNA plot<br />5) module hub genes |
| 7_module_analysis           | 1) align T2D and height SNPs to nearest genes<br />2) regional specific marker enrichment<br />3) prepare SNP enrichment<br />4) GWAS plot<br />5) GO and KEGG pathway analysis |





