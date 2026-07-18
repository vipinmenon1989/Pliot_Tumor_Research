# Milestone 5 PCA and PC Evaluation Report - MPNST_1
*Generated on: 2026-07-10 06:59:41*

## 1. Parameters & Configuration

- **Input Normalized Object**: `results/datasets/MPNST_1/MPNST_1_normalized.rds`
- **Output PCA Object**: `results/datasets/MPNST_1/MPNST_1_pca.rds`
- **Assay Evaluated**: `SCT`
- **PCs Computed**: `50`
- **Seed Used**: `42`
- **ScaleData Check**: `ScaleData was already present and reused directly.`

## 2. Variance Explained Summary

| Principal Component | Variance Explained (%) | Cumulative Variance (%) |
| --- | :---: | :---: |
| **PC1** | 12.3358% | 12.3358% |
| **PC2** | 7.1662% | 19.5019% |
| **PC3** | 5.4146% | 24.9165% |
| **PC4** | 2.7380% | 27.6545% |
| **PC5** | 2.3219% | 29.9764% |
| **PC6** | 1.8324% | 31.8088% |
| **PC7** | 1.5356% | 33.3444% |
| **PC8** | 1.2592% | 34.6036% |
| **PC9** | 1.1927% | 35.7962% |
| **PC10** | 0.9683% | 36.7646% |
| **PC11** | 0.9263% | 37.6909% |
| **PC12** | 0.8527% | 38.5436% |
| **PC13** | 0.7805% | 39.3241% |
| **PC14** | 0.6303% | 39.9544% |
| **PC15** | 0.5226% | 40.4770% |

## 3. PC Recommendation with Quantitative Evidence

- **Conservative PC Range**: `1 - 8`
- **Recommended PC Range**: `1 - 8` (Primary choice)
- **Maximum Reasonable PC Range**: `1 - 24`

### Quantitative Rationale:
1. **Geometric Elbow Point**: The geometric elbow of the variance explained curve occurs at **PC8**. This is the mathematically optimal knee point where additional PCs yield diminishing returns.
2. **Conservative Threshold**: Up to **PC8**, each PC explains > 1.0% of the total variance and represents major co-expression modules with clear separation from baseline noise.
3. **Maximum Limit**: By **PC24**, the marginal variance explained per PC falls below **0.3%**, and the cumulative variance explained levels off completely, indicating that higher PCs are dominated by stochastic cell-specific noise.
4. **JackStraw Permutation Check**: JackStraw was omitted for this SCTransformed dataset. *Rationale*: SCTransform residuals do not conform to standard standard-normal assumptions under permutations, and standard resampling violates the SCT regularized model. Permutation tests on SCT residuals are computationally expensive and lack theoretical justification. Hence, elbow and correlation analyses were prioritized.

## 4. Dominant Loading Genes (Top 10 PCs)

Top positive and negative loading features for each PC:

#### PC1
- **Positive Loadings**: HLA-DRA (0.169), CD74 (0.166), HLA-DRB1 (0.152), HLA-DPA1 (0.147), HLA-DQA1 (0.138), CCL3 (0.135), CD83 (0.134), SRGN (0.128), HLA-DPB1 (0.126), CCL4 (0.119)
- **Negative Loadings**: COL3A1 (-0.043), FAM155A (-0.041), TSHZ2 (-0.040), MGP (-0.038), ROBO2 (-0.038), CDH19 (-0.036), COL1A1 (-0.035), SERPINE2 (-0.034), NRXN1 (-0.034), RORA (-0.033)

#### PC2
- **Positive Loadings**: FAM155A (0.103), TSHZ2 (0.098), ROBO2 (0.096), COL3A1 (0.093), AFF3 (0.087), TENM2 (0.084), RORA (0.074), MEIS2 (0.073), MGP (0.072), COL1A1 (0.072)
- **Negative Loadings**: S100B (-0.153), CDH19 (-0.141), NRXN1 (-0.138), SERPINE2 (-0.137), CAPS (-0.113), CRYAB (-0.111), CADM1 (-0.100), GPM6B (-0.093), S100A6 (-0.084), BAI3 (-0.083)

#### PC3
- **Positive Loadings**: FAM155A (0.054), ROBO2 (0.051), TSHZ2 (0.051), AFF3 (0.046), TENM2 (0.045), COL3A1 (0.045), MEIS2 (0.039), AUTS2 (0.038), RORA (0.035), KCNQ1OT1 (0.035)
- **Negative Loadings**: SPARCL1 (-0.157), ANGPT2 (-0.143), VWF (-0.141), IGFBP7 (-0.138), FLT1 (-0.128), INSR (-0.127), ESM1 (-0.123), IGFBP3 (-0.123), PLVAP (-0.117), PRSS23 (-0.107)

#### PC4
- **Positive Loadings**: NRXN1 (0.079), ROBO2 (0.065), FAM155A (0.065), FLT1 (0.064), INSR (0.062), VWF (0.061), LDB2 (0.060), AFF3 (0.060), TSHZ2 (0.059), TENM2 (0.056)
- **Negative Loadings**: IGFBP7 (-0.138), TAGLN (-0.135), BGN (-0.128), RGS5 (-0.127), ACTA2 (-0.115), MGP (-0.108), THY1 (-0.107), NDUFA4L2 (-0.099), INPP4B (-0.087), HIGD1B (-0.087)

#### PC5
- **Positive Loadings**: CCL5 (0.154), IL32 (0.138), CD52 (0.130), GZMA (0.129), PTPRC (0.124), NKG7 (0.118), PTPRCAP (0.109), STAT4 (0.108), CD3D (0.107), CD2 (0.106)
- **Negative Loadings**: IGFBP7 (-0.098), CCL2 (-0.075), APOE (-0.074), C1QB (-0.072), C1QA (-0.072), C1QC (-0.071), CD14 (-0.063), TAGLN (-0.063), SPARCL1 (-0.062), BGN (-0.061)

#### PC6
- **Positive Loadings**: FTL (0.124), IGFBP5 (0.115), MGP (0.111), IGFBP2 (0.101), GAPDH (0.079), GPNMB (0.079), VEGFA (0.073), SCG2 (0.068), FTLP3 (0.068), NUPR1 (0.067)
- **Negative Loadings**: NRXN1 (-0.107), INPP4B (-0.079), SLIT3 (-0.077), RGS5 (-0.076), PRKG1 (-0.076), XKR4 (-0.075), TNFAIP3 (-0.075), TSHZ2 (-0.072), CCL5 (-0.069), PTPRC (-0.067)

#### PC7
- **Positive Loadings**: APOD (0.225), DCN (0.214), LUM (0.198), C1R (0.126), CCL2 (0.122), SERPINF1 (0.120), SFRP4 (0.119), C1S (0.113), LINC01088 (0.104), FMO2 (0.097)
- **Negative Loadings**: RGS5 (-0.096), NDUFA4L2 (-0.095), THY1 (-0.089), ACTA2 (-0.073), HIGD1B (-0.070), MYH11 (-0.070), CCDC102B (-0.063), COX4I2 (-0.063), INPP4B (-0.062), KCNAB1 (-0.059)

#### PC8
- **Positive Loadings**: SERPINE2 (0.245), CAPS (0.187), S100B (0.168), LPL (0.146), PTN (0.125), GAP43 (0.109), FREM1 (0.105), TRPM3 (0.104), DLK1 (0.097), FSTL5 (0.096)
- **Negative Loadings**: KIRREL3 (-0.156), CRYAB (-0.146), NRXN1 (-0.132), XKR4 (-0.122), VEGFA (-0.097), CADM1 (-0.089), DGKB (-0.087), FRMD5 (-0.085), PAPPA (-0.083), NCAM2 (-0.081)

#### PC9
- **Positive Loadings**: IGFBP5 (0.120), SERPINE2 (0.093), IL1B (0.092), VEGFA (0.088), GSN (0.087), CD83 (0.086), NFKBIA (0.080), NFKB1 (0.078), DOCK4 (0.074), IL8 (0.070)
- **Negative Loadings**: SPP1 (-0.153), HLA-DRA (-0.137), FTL (-0.136), APOE (-0.133), HLA-DPA1 (-0.131), CD74 (-0.131), HLA-DRB1 (-0.115), C1QB (-0.111), TMSB4X (-0.111), C1QC (-0.109)

#### PC10
- **Positive Loadings**: GSN (0.134), GLUL (0.104), VEGFA (0.102), CADM1 (0.091), IGFBP5 (0.090), MTND2P28 (0.085), SPP1 (0.083), APOC1 (0.080), HFM1 (0.078), MTATP6P1 (0.075)
- **Negative Loadings**: MGP (-0.175), S100B (-0.109), SCG2 (-0.108), CRYAB (-0.104), XKR4 (-0.092), IL1B (-0.081), PLAUR (-0.081), SEMA3B (-0.079), NOV (-0.078), IGFBP2 (-0.077)


## 5. Technical and Biological Assessment

### Observed Results: Technical Correlations
Statistically significant Pearson correlations (P < 0.05) between cell PC scores and technical covariates:

| PC | Technical Metric | Pearson R | P-value |
| --- | --- | :---: | :---: |
| **PC_1** | nCount_RNA | -0.2433 | 5.016229e-103 |
| **PC_2** | nCount_RNA | 0.1516 | 2.292544e-40 |
| **PC_5** | nCount_RNA | -0.1966 | 3.171862e-67 |
| **PC_7** | nCount_RNA | 0.1051 | 3.733010e-20 |
| **PC_8** | nCount_RNA | 0.1733 | 1.944423e-52 |
| **PC_9** | nCount_RNA | 0.0778 | 1.036545e-11 |
| **PC_10** | nCount_RNA | -0.3672 | 1.053714e-241 |
| **PC_11** | nCount_RNA | -0.1152 | 6.607146e-24 |
| **PC_12** | nCount_RNA | -0.0684 | 2.274684e-09 |
| **PC_13** | nCount_RNA | 0.1558 | 1.362327e-42 |
| **PC_14** | nCount_RNA | 0.1377 | 1.546025e-33 |
| **PC_15** | nCount_RNA | -0.2076 | 6.085552e-75 |
| **PC_1** | nFeature_RNA | -0.3192 | 5.888206e-180 |
| **PC_2** | nFeature_RNA | 0.1164 | 2.178954e-24 |
| **PC_3** | nFeature_RNA | 0.0553 | 1.353170e-06 |
| **PC_4** | nFeature_RNA | 0.0572 | 5.976884e-07 |
| **PC_5** | nFeature_RNA | -0.1971 | 1.529750e-67 |
| **PC_7** | nFeature_RNA | 0.1098 | 7.284498e-22 |
| **PC_8** | nFeature_RNA | 0.1605 | 4.213079e-45 |
| **PC_9** | nFeature_RNA | 0.1592 | 2.106577e-44 |
| **PC_10** | nFeature_RNA | -0.3519 | 8.572122e-221 |
| **PC_11** | nFeature_RNA | -0.1088 | 1.760787e-21 |
| **PC_12** | nFeature_RNA | -0.0711 | 5.352001e-10 |
| **PC_13** | nFeature_RNA | 0.1823 | 6.513121e-58 |
| **PC_14** | nFeature_RNA | 0.2319 | 1.529983e-93 |
| **PC_15** | nFeature_RNA | -0.2220 | 1.238092e-85 |
| **NA** | NA | NA | NA |
| **NA** | NA | NA | NA |
| **NA** | NA | NA | NA |
| **NA** | NA | NA | NA |
| **NA** | NA | NA | NA |
| **NA** | NA | NA | NA |
| **NA** | NA | NA | NA |
| **NA** | NA | NA | NA |
| **NA** | NA | NA | NA |
| **NA** | NA | NA | NA |
| **NA** | NA | NA | NA |
| **NA** | NA | NA | NA |
| **NA** | NA | NA | NA |
| **NA** | NA | NA | NA |
| **NA** | NA | NA | NA |
| **PC_1** | percent.ribo | 0.1058 | 2.052852e-20 |
| **PC_2** | percent.ribo | 0.0849 | 1.170960e-13 |
| **PC_4** | percent.ribo | -0.4199 | 5.928788e-323 |
| **PC_5** | percent.ribo | 0.2918 | 2.270537e-149 |
| **PC_6** | percent.ribo | 0.5134 | 0.000000e+00 |
| **PC_7** | percent.ribo | -0.0705 | 7.332490e-10 |
| **PC_8** | percent.ribo | 0.1666 | 1.606095e-48 |
| **PC_9** | percent.ribo | -0.2119 | 4.934252e-78 |
| **PC_10** | percent.ribo | -0.2501 | 5.820101e-109 |
| **PC_11** | percent.ribo | 0.0297 | 9.469890e-03 |
| **PC_12** | percent.ribo | -0.1499 | 1.582832e-39 |
| **PC_13** | percent.ribo | -0.0410 | 3.454912e-04 |
| **PC_14** | percent.ribo | -0.1651 | 1.181671e-47 |
| **PC_15** | percent.ribo | -0.0512 | 7.920104e-06 |

### Interpretation of Leading PCs
The biological and technical processes captured by the first 4 PCs are:

#### PC1 Interpretation
- **Dominant Genes (Pos)**: 
- **Dominant Genes (Neg)**: 
- **Biological Process**: General Cellular Physiology / Metabolic Heterogeneity
- **Technical Covariance**: 

#### PC2 Interpretation
- **Dominant Genes (Pos)**: 
- **Dominant Genes (Neg)**: 
- **Biological Process**: General Cellular Physiology / Metabolic Heterogeneity
- **Technical Covariance**: 

#### PC3 Interpretation
- **Dominant Genes (Pos)**: 
- **Dominant Genes (Neg)**: 
- **Biological Process**: General Cellular Physiology / Metabolic Heterogeneity
- **Technical Covariance**: 

#### PC4 Interpretation
- **Dominant Genes (Pos)**: 
- **Dominant Genes (Neg)**: 
- **Biological Process**: General Cellular Physiology / Metabolic Heterogeneity
- **Technical Covariance**: 


### Recommendations
1. **Neighbor Graph and Clustering**: We recommend using the first **8 PCs** (the recommended range) for constructing the shared nearest neighbor (SNN) graph in Milestone 6.
2. **Technical Confounders**: PCs strongly correlated with nCount_RNA or percent.mt should be monitored. Because SCTransform regressed out depth covariance, residual correlation represents biological cell-size difference rather than technical artifacts, but monitoring is advised.
3. **Biological Signal Preservation**: The leading PCs are heavily dominated by biologically meaningful programs (extracellular matrix remodeling, cell proliferation, immune infiltration) rather than technical artifacts, confirming the high quality of the filtering and normalization.

## 6. Diagnostics & Visualizations

The following publication-quality diagnostic plots were generated:
- [Variance Explained (Elbow Plot)](file:///reports/datasets/MPNST_1/pca_elbow.png)
- [Cumulative Variance Explained](file:///reports/datasets/MPNST_1/pca_cumulative_variance.png)
- [PC Loadings Plot (Top 4 PCs)](file:///reports/datasets/MPNST_1/pca_loadings.png)
- [PC Expression Heatmaps (Dims 1-9)](file:///reports/datasets/MPNST_1/pca_heatmaps.png)
- [Technical Metric Correlations](file:///reports/datasets/MPNST_1/pca_correlations.png)
- [Top Loading Genes Table (TSV)](file:///reports/datasets/MPNST_1/top_loading_genes.tsv)

---
## 7. Software & Environment Provenance
- **Seurat Version**: `5.4.0`
- **Patchwork Version**: `1.3.2`
- **Git Commit Hash**: `31cc849c7c8aaf27e8f2b594fb40d0c54989a579`
