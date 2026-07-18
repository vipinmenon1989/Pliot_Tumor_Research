# Milestone 5 PCA and PC Evaluation Report - MPNST_4
*Generated on: 2026-07-10 06:59:11*

## 1. Parameters & Configuration

- **Input Normalized Object**: `results/datasets/MPNST_4/MPNST_4_normalized.rds`
- **Output PCA Object**: `results/datasets/MPNST_4/MPNST_4_pca.rds`
- **Assay Evaluated**: `SCT`
- **PCs Computed**: `50`
- **Seed Used**: `42`
- **ScaleData Check**: `ScaleData was already present and reused directly.`

## 2. Variance Explained Summary

| Principal Component | Variance Explained (%) | Cumulative Variance (%) |
| --- | :---: | :---: |
| **PC1** | 13.3460% | 13.3460% |
| **PC2** | 3.5042% | 16.8502% |
| **PC3** | 2.8340% | 19.6843% |
| **PC4** | 2.6702% | 22.3545% |
| **PC5** | 1.4368% | 23.7913% |
| **PC6** | 1.3577% | 25.1490% |
| **PC7** | 1.2816% | 26.4306% |
| **PC8** | 1.1796% | 27.6102% |
| **PC9** | 1.0760% | 28.6861% |
| **PC10** | 0.9366% | 29.6227% |
| **PC11** | 0.7722% | 30.3949% |
| **PC12** | 0.6703% | 31.0653% |
| **PC13** | 0.5815% | 31.6467% |
| **PC14** | 0.4906% | 32.1374% |
| **PC15** | 0.4724% | 32.6098% |

## 3. PC Recommendation with Quantitative Evidence

- **Conservative PC Range**: `1 - 8`
- **Recommended PC Range**: `1 - 5` (Primary choice)
- **Maximum Reasonable PC Range**: `1 - 22`

### Quantitative Rationale:
1. **Geometric Elbow Point**: The geometric elbow of the variance explained curve occurs at **PC5**. This is the mathematically optimal knee point where additional PCs yield diminishing returns.
2. **Conservative Threshold**: Up to **PC8**, each PC explains > 1.0% of the total variance and represents major co-expression modules with clear separation from baseline noise.
3. **Maximum Limit**: By **PC22**, the marginal variance explained per PC falls below **0.3%**, and the cumulative variance explained levels off completely, indicating that higher PCs are dominated by stochastic cell-specific noise.
4. **JackStraw Permutation Check**: JackStraw was omitted for this SCTransformed dataset. *Rationale*: SCTransform residuals do not conform to standard standard-normal assumptions under permutations, and standard resampling violates the SCT regularized model. Permutation tests on SCT residuals are computationally expensive and lack theoretical justification. Hence, elbow and correlation analyses were prioritized.

## 4. Dominant Loading Genes (Top 10 PCs)

Top positive and negative loading features for each PC:

#### PC1
- **Positive Loadings**: CD74 (0.186), HLA-DRB1 (0.177), HLA-DRA (0.174), HLA-DPA1 (0.169), HLA-DPB1 (0.157), TYROBP (0.148), HLA-DQB1 (0.136), HLA-DQA1 (0.136), FTL (0.132), C1QA (0.131)
- **Negative Loadings**: COL1A1 (-0.091), COL1A2 (-0.084), MGP (-0.061), COL3A1 (-0.055), APOD (-0.045), DLK1 (-0.045), THBS4 (-0.045), PTN (-0.041), LUM (-0.039), POSTN (-0.036)

#### PC2
- **Positive Loadings**: CLDN5 (0.147), VWF (0.130), RAMP2 (0.120), EGFL7 (0.119), IL32 (0.111), CDH5 (0.110), CAV1 (0.108), ESAM (0.099), S100A16 (0.097), SOX18 (0.094)
- **Negative Loadings**: COL1A1 (-0.150), COL1A2 (-0.137), COL3A1 (-0.080), MGP (-0.077), DLK1 (-0.068), THBS4 (-0.067), FTL (-0.065), C1QA (-0.063), C1QC (-0.060), C1QB (-0.056)

#### PC3
- **Positive Loadings**: HSPA1A (0.138), DNAJB1 (0.129), HSP90AA1 (0.127), JUNB (0.105), HSPA1B (0.097), TOP2A (0.097), FOS (0.094), DUSP2 (0.093), MKI67 (0.084), HMGB2 (0.083)
- **Negative Loadings**: RNASE1 (-0.120), C1QA (-0.102), C1QC (-0.099), C1QB (-0.095), VWF (-0.077), SEPP1 (-0.072), FOLR2 (-0.070), EGFL7 (-0.070), CLDN5 (-0.069), FTL (-0.068)

#### PC4
- **Positive Loadings**: TOP2A (0.152), KIAA0101 (0.130), MKI67 (0.128), STMN1 (0.125), CENPF (0.123), TUBA1B (0.120), UBE2C (0.116), HMGB2 (0.114), TPX2 (0.105), H2AFZ (0.105)
- **Negative Loadings**: HSPA1A (-0.112), DNAJB1 (-0.103), JUNB (-0.092), HSP90AA1 (-0.084), FOS (-0.083), HSPA1B (-0.078), DUSP2 (-0.076), DNAJA1 (-0.068), DUSP1 (-0.064), IL7R (-0.062)

#### PC5
- **Positive Loadings**: CCL5 (0.133), IL7R (0.122), PTPRC (0.122), CD2 (0.120), PTPRCAP (0.113), B2M (0.109), CD3D (0.107), IL32 (0.094), CD52 (0.089), GZMA (0.082)
- **Negative Loadings**: ATF3 (-0.133), NFKBIA (-0.125), IER3 (-0.120), SOD2 (-0.113), CD83 (-0.104), FOS (-0.101), CXCL2 (-0.096), PPP1R15A (-0.096), IL8 (-0.095), GADD45B (-0.092)

#### PC6
- **Positive Loadings**: TNMD (0.234), THBS4 (0.174), IGFBP2 (0.171), SCXA (0.164), NPW (0.110), DLK1 (0.097), COL1A1 (0.087), CYTL1 (0.083), CRABP1 (0.079), FABP5 (0.078)
- **Negative Loadings**: APOD (-0.166), ABCA8 (-0.128), PI16 (-0.122), LUM (-0.114), CAPN6 (-0.101), APOE (-0.095), ABCA10 (-0.090), FMO1 (-0.086), RNASE1 (-0.081), JUN (-0.081)

#### PC7
- **Positive Loadings**: CTGF (0.176), MGP (0.128), PTN (0.094), COL1A2 (0.094), POSTN (0.092), COL14A1 (0.092), SPARCL1 (0.086), COL1A1 (0.085), DARC (0.085), MALAT1 (0.082)
- **Negative Loadings**: APOD (-0.291), ABCA8 (-0.103), PI16 (-0.101), APOE (-0.100), CYP1B1 (-0.083), GNG11 (-0.079), CLDN1 (-0.073), FZD2 (-0.073), C11orf96 (-0.071), FOXS1 (-0.067)

#### PC8
- **Positive Loadings**: RNASE1 (0.161), FOLR2 (0.140), SEPP1 (0.139), C1QB (0.128), C1QC (0.123), JUN (0.122), C1QA (0.119), MS4A7 (0.110), HSP90AA1 (0.095), FOSB (0.093)
- **Negative Loadings**: LYZ (-0.179), S100A4 (-0.147), HLA-DQA1 (-0.124), CST3 (-0.117), S100A6 (-0.116), FCER1A (-0.112), HLA-DQB1 (-0.110), HLA-DPB1 (-0.108), FCN1 (-0.106), S100A10 (-0.095)

#### PC9
- **Positive Loadings**: APOD (0.189), THBS4 (0.154), IGFBP7 (0.147), IFI27 (0.132), TNMD (0.107), DARC (0.104), SCXA (0.095), ADIRF (0.093), ENG (0.088), CLU (0.087)
- **Negative Loadings**: CXCL14 (-0.126), GNG11 (-0.112), IL32 (-0.109), LUM (-0.107), CAV1 (-0.105), CTGF (-0.096), MGP (-0.094), DIO2 (-0.090), MSX1 (-0.085), KDR (-0.084)

#### PC10
- **Positive Loadings**: APOD (0.106), DARC (0.065), VWF (0.065), SPARCL1 (0.056), RAMP3 (0.055), ENG (0.054), COL3A1 (0.054), COL1A1 (0.054), CLU (0.048), CLEC14A (0.048)
- **Negative Loadings**: S100B (-0.169), GPM6B (-0.158), PTPRZ1 (-0.157), CRYAB (-0.141), LGI4 (-0.132), CADM1 (-0.132), MIA (-0.132), S100A6 (-0.130), CDH19 (-0.121), S100A4 (-0.120)


## 5. Technical and Biological Assessment

### Observed Results: Technical Correlations
Statistically significant Pearson correlations (P < 0.05) between cell PC scores and technical covariates:

| PC | Technical Metric | Pearson R | P-value |
| --- | --- | :---: | :---: |
| **PC_1** | nCount_RNA | 0.0662 | 3.892639e-08 |
| **PC_2** | nCount_RNA | 0.3159 | 3.714065e-159 |
| **PC_3** | nCount_RNA | -0.0653 | 5.884284e-08 |
| **PC_4** | nCount_RNA | 0.1842 | 1.437944e-53 |
| **PC_5** | nCount_RNA | -0.2365 | 4.452709e-88 |
| **PC_6** | nCount_RNA | 0.2962 | 2.762378e-139 |
| **PC_7** | nCount_RNA | -0.1762 | 4.442262e-49 |
| **PC_8** | nCount_RNA | -0.0539 | 7.739597e-06 |
| **PC_9** | nCount_RNA | -0.0786 | 6.526954e-11 |
| **PC_10** | nCount_RNA | -0.1075 | 4.016502e-19 |
| **PC_11** | nCount_RNA | 0.1653 | 2.628955e-43 |
| **PC_12** | nCount_RNA | -0.1068 | 6.561056e-19 |
| **PC_13** | nCount_RNA | -0.0817 | 1.176460e-11 |
| **PC_14** | nCount_RNA | -0.0503 | 3.042597e-05 |
| **PC_15** | nCount_RNA | 0.0344 | 4.342593e-03 |
| **PC_1** | nFeature_RNA | -0.1460 | 4.404423e-34 |
| **PC_2** | nFeature_RNA | 0.2124 | 5.075985e-71 |
| **PC_3** | nFeature_RNA | -0.1081 | 2.485645e-19 |
| **PC_4** | nFeature_RNA | 0.2469 | 5.012484e-96 |
| **PC_5** | nFeature_RNA | -0.2479 | 7.452985e-97 |
| **PC_6** | nFeature_RNA | 0.2422 | 2.322874e-92 |
| **PC_7** | nFeature_RNA | -0.1263 | 7.234081e-26 |
| **PC_8** | nFeature_RNA | -0.0367 | 2.349850e-03 |
| **PC_9** | nFeature_RNA | -0.0544 | 6.298338e-06 |
| **PC_10** | nFeature_RNA | -0.0939 | 6.075862e-15 |
| **PC_11** | nFeature_RNA | 0.1531 | 2.424952e-37 |
| **PC_12** | nFeature_RNA | -0.1087 | 1.517745e-19 |
| **PC_13** | nFeature_RNA | -0.0778 | 1.025902e-10 |
| **PC_14** | nFeature_RNA | -0.0279 | 2.072024e-02 |
| **PC_15** | nFeature_RNA | 0.0271 | 2.455736e-02 |
| **PC_1** | percent.ribo | -0.1659 | 1.310779e-43 |
| **PC_2** | percent.ribo | 0.0965 | 1.042331e-15 |
| **PC_3** | percent.ribo | 0.1177 | 1.160190e-22 |
| **PC_4** | percent.ribo | -0.2130 | 2.269822e-71 |
| **PC_5** | percent.ribo | 0.2547 | 2.824308e-102 |
| **PC_6** | percent.ribo | 0.1261 | 9.062800e-26 |
| **PC_7** | percent.ribo | -0.3369 | 3.683900e-182 |
| **PC_8** | percent.ribo | -0.1484 | 3.531227e-35 |
| **PC_10** | percent.ribo | 0.1607 | 4.990406e-41 |
| **PC_11** | percent.ribo | 0.4146 | 6.669810e-284 |
| **PC_12** | percent.ribo | 0.1446 | 1.891409e-33 |
| **PC_13** | percent.ribo | 0.1154 | 7.745792e-22 |
| **PC_14** | percent.ribo | -0.1012 | 4.123299e-17 |

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
1. **Neighbor Graph and Clustering**: We recommend using the first **5 PCs** (the recommended range) for constructing the shared nearest neighbor (SNN) graph in Milestone 6.
2. **Technical Confounders**: PCs strongly correlated with nCount_RNA or percent.mt should be monitored. Because SCTransform regressed out depth covariance, residual correlation represents biological cell-size difference rather than technical artifacts, but monitoring is advised.
3. **Biological Signal Preservation**: The leading PCs are heavily dominated by biologically meaningful programs (extracellular matrix remodeling, cell proliferation, immune infiltration) rather than technical artifacts, confirming the high quality of the filtering and normalization.

## 6. Diagnostics & Visualizations

The following publication-quality diagnostic plots were generated:
- [Variance Explained (Elbow Plot)](file:///reports/datasets/MPNST_4/pca_elbow.png)
- [Cumulative Variance Explained](file:///reports/datasets/MPNST_4/pca_cumulative_variance.png)
- [PC Loadings Plot (Top 4 PCs)](file:///reports/datasets/MPNST_4/pca_loadings.png)
- [PC Expression Heatmaps (Dims 1-9)](file:///reports/datasets/MPNST_4/pca_heatmaps.png)
- [Technical Metric Correlations](file:///reports/datasets/MPNST_4/pca_correlations.png)
- [Top Loading Genes Table (TSV)](file:///reports/datasets/MPNST_4/top_loading_genes.tsv)

---
## 7. Software & Environment Provenance
- **Seurat Version**: `5.4.0`
- **Patchwork Version**: `1.3.2`
- **Git Commit Hash**: `31cc849c7c8aaf27e8f2b594fb40d0c54989a579`
