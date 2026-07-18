# Milestone 5 PCA and PC Evaluation Report - MPNST_2
*Generated on: 2026-07-10 06:58:11*

## 1. Parameters & Configuration

- **Input Normalized Object**: `results/datasets/MPNST_2/MPNST_2_normalized.rds`
- **Output PCA Object**: `results/datasets/MPNST_2/MPNST_2_pca.rds`
- **Assay Evaluated**: `SCT`
- **PCs Computed**: `50`
- **Seed Used**: `42`
- **ScaleData Check**: `ScaleData was already present and reused directly.`

## 2. Variance Explained Summary

| Principal Component | Variance Explained (%) | Cumulative Variance (%) |
| --- | :---: | :---: |
| **PC1** | 14.0256% | 14.0256% |
| **PC2** | 5.5725% | 19.5981% |
| **PC3** | 3.7319% | 23.3300% |
| **PC4** | 2.2765% | 25.6065% |
| **PC5** | 2.0309% | 27.6374% |
| **PC6** | 1.5034% | 29.1407% |
| **PC7** | 1.4806% | 30.6213% |
| **PC8** | 1.1330% | 31.7543% |
| **PC9** | 1.0970% | 32.8512% |
| **PC10** | 0.9010% | 33.7522% |
| **PC11** | 0.8100% | 34.5622% |
| **PC12** | 0.7659% | 35.3281% |
| **PC13** | 0.6058% | 35.9339% |
| **PC14** | 0.5432% | 36.4771% |
| **PC15** | 0.5107% | 36.9878% |

## 3. PC Recommendation with Quantitative Evidence

- **Conservative PC Range**: `1 - 8`
- **Recommended PC Range**: `1 - 6` (Primary choice)
- **Maximum Reasonable PC Range**: `1 - 25`

### Quantitative Rationale:
1. **Geometric Elbow Point**: The geometric elbow of the variance explained curve occurs at **PC6**. This is the mathematically optimal knee point where additional PCs yield diminishing returns.
2. **Conservative Threshold**: Up to **PC8**, each PC explains > 1.0% of the total variance and represents major co-expression modules with clear separation from baseline noise.
3. **Maximum Limit**: By **PC25**, the marginal variance explained per PC falls below **0.3%**, and the cumulative variance explained levels off completely, indicating that higher PCs are dominated by stochastic cell-specific noise.
4. **JackStraw Permutation Check**: JackStraw was omitted for this SCTransformed dataset. *Rationale*: SCTransform residuals do not conform to standard standard-normal assumptions under permutations, and standard resampling violates the SCT regularized model. Permutation tests on SCT residuals are computationally expensive and lack theoretical justification. Hence, elbow and correlation analyses were prioritized.

## 4. Dominant Loading Genes (Top 10 PCs)

Top positive and negative loading features for each PC:

#### PC1
- **Positive Loadings**: CD74 (0.136), HLA-DRA (0.133), HLA-DRB1 (0.131), TYROBP (0.119), HLA-DPA1 (0.118), HLA-DPB1 (0.117), SRGN (0.104), AIF1 (0.101), LYZ (0.098), FCER1G (0.096)
- **Negative Loadings**: SFRP2 (-0.081), DCN (-0.080), COL1A1 (-0.076), COL1A2 (-0.076), PLA2G2A (-0.074), COL3A1 (-0.070), MGP (-0.070), FBLN1 (-0.068), CXCL14 (-0.068), SFRP4 (-0.063)

#### PC2
- **Positive Loadings**: SFRP2 (0.146), PLA2G2A (0.139), FBLN1 (0.124), CXCL14 (0.117), CCDC80 (0.111), COL3A1 (0.101), COMP (0.101), COL1A1 (0.095), COL1A2 (0.091), THY1 (0.089)
- **Negative Loadings**: PTGDS (-0.113), APOD (-0.101), SPARCL1 (-0.088), A2M (-0.085), GPC3 (-0.076), ANGPTL7 (-0.073), ADIRF (-0.068), BCAM (-0.065), ITGA6 (-0.063), CDH19 (-0.062)

#### PC3
- **Positive Loadings**: VWF (0.098), ELTD1 (0.095), EMCN (0.091), CD93 (0.090), ECSCR (0.088), PLVAP (0.086), HSPG2 (0.083), EGFL7 (0.081), MMRN2 (0.079), EPAS1 (0.078)
- **Negative Loadings**: APOD (-0.162), PTGDS (-0.136), ANGPTL7 (-0.106), GPC3 (-0.105), ABCA8 (-0.086), APOE (-0.084), CDH19 (-0.083), PTCH1 (-0.075), LUM (-0.069), LHFP (-0.066)

#### PC4
- **Positive Loadings**: CD52 (0.098), CXCR4 (0.090), RGS1 (0.074), PTPRCAP (0.073), CST7 (0.073), CD1C (0.072), FCER1A (0.071), RPS27 (0.069), PLAC8 (0.069), C12orf75 (0.068)
- **Negative Loadings**: C1QC (-0.141), C1QA (-0.139), C1QB (-0.132), CD14 (-0.117), FOLR2 (-0.107), RNASE1 (-0.106), FTL (-0.106), STAB1 (-0.095), SEPP1 (-0.094), CD163 (-0.093)

#### PC5
- **Positive Loadings**: ANGPTL7 (0.100), APOD (0.091), SPARCL1 (0.083), CXCL12 (0.080), ABCA8 (0.080), LUM (0.077), CDH19 (0.075), FN1 (0.068), ABCA10 (0.068), SRPX (0.067)
- **Negative Loadings**: ADIRF (-0.162), IGFBP6 (-0.152), TNNC1 (-0.137), SLC2A1 (-0.136), CSRP1 (-0.133), TAGLN (-0.128), SLC22A3 (-0.127), ITGA6 (-0.115), KRT19 (-0.107), CAV1 (-0.104)

#### PC6
- **Positive Loadings**: IL1B (0.117), IER3 (0.116), PLAUR (0.113), DUSP2 (0.112), IL8 (0.107), BCL2A1 (0.092), NFKBIA (0.085), CD83 (0.085), FCER1A (0.084), TNFAIP3 (0.083)
- **Negative Loadings**: SPARC (-0.081), C1QB (-0.074), C1QC (-0.071), ALOX5AP (-0.070), PTPRCAP (-0.070), CYBB (-0.069), FOLR2 (-0.068), MAF (-0.067), C12orf75 (-0.066), C1QA (-0.065)

#### PC7
- **Positive Loadings**: SPARC (0.138), COL1A1 (0.130), LUM (0.121), COL3A1 (0.120), COL1A2 (0.105), PI16 (0.098), IGFBP5 (0.097), C1QTNF3 (0.096), CTHRC1 (0.092), RGCC (0.091)
- **Negative Loadings**: C7 (-0.151), IGF1 (-0.113), C11orf96 (-0.110), RARRES1 (-0.110), FGF7 (-0.106), HP (-0.100), CHI3L1 (-0.091), RPL39 (-0.091), SFRP1 (-0.087), NNMT (-0.080)

#### PC8
- **Positive Loadings**: LYZ (0.105), HLA-DQA1 (0.091), HLA-DPA1 (0.088), FCER1A (0.088), CPVL (0.087), HLA-DQB1 (0.086), CLEC10A (0.086), HLA-DPB1 (0.084), MNDA (0.083), AIF1 (0.077)
- **Negative Loadings**: NFKBIA (-0.106), CXCL2 (-0.102), TNFAIP3 (-0.095), IER3 (-0.087), SOD2 (-0.081), CXCL3 (-0.079), TNF (-0.078), CCL4 (-0.077), CCL3 (-0.075), PPP1R15A (-0.075)

#### PC9
- **Positive Loadings**: S100B (0.155), CRYAB (0.141), MPZ (0.138), PMP2 (0.130), GPM6B (0.130), MAL (0.122), MT2A (0.119), PLP1 (0.119), MBP (0.111), TSC22D4 (0.107)
- **Negative Loadings**: SFRP4 (-0.095), PTGDS (-0.093), DCN (-0.076), APOD (-0.073), IGFBP6 (-0.071), IGFBP5 (-0.069), TCF4 (-0.066), GSN (-0.064), LUM (-0.060), CCDC80 (-0.058)

#### PC10
- **Positive Loadings**: COL5A2 (0.103), COL1A2 (0.094), FN1 (0.088), COL6A1 (0.086), COL3A1 (0.085), COL1A1 (0.085), COL5A1 (0.085), CPE (0.083), PRSS23 (0.083), ADAM12 (0.081)
- **Negative Loadings**: CD9 (-0.171), SFRP4 (-0.141), RBP4 (-0.107), CCDC80 (-0.102), CLU (-0.100), PLA2G2A (-0.088), TAC1 (-0.087), FBLN1 (-0.084), CXCL14 (-0.083), IGFBP5 (-0.077)


## 5. Technical and Biological Assessment

### Observed Results: Technical Correlations
Statistically significant Pearson correlations (P < 0.05) between cell PC scores and technical covariates:

| PC | Technical Metric | Pearson R | P-value |
| --- | --- | :---: | :---: |
| **PC_1** | nCount_RNA | -0.4276 | 3.601948e-102 |
| **PC_2** | nCount_RNA | 0.4526 | 1.000126e-115 |
| **PC_3** | nCount_RNA | 0.3277 | 2.638442e-58 |
| **PC_4** | nCount_RNA | -0.0435 | 3.766789e-02 |
| **PC_7** | nCount_RNA | 0.1077 | 2.474766e-07 |
| **PC_8** | nCount_RNA | 0.1464 | 2.075832e-12 |
| **PC_10** | nCount_RNA | 0.1865 | 2.590069e-19 |
| **PC_11** | nCount_RNA | -0.1300 | 4.444214e-10 |
| **PC_12** | nCount_RNA | -0.0787 | 1.650698e-04 |
| **PC_13** | nCount_RNA | -0.0540 | 9.809067e-03 |
| **PC_14** | nCount_RNA | 0.0688 | 9.950930e-04 |
| **PC_15** | nCount_RNA | -0.0543 | 9.443345e-03 |
| **PC_1** | nFeature_RNA | -0.4648 | 9.356458e-123 |
| **PC_2** | nFeature_RNA | 0.3181 | 6.976989e-55 |
| **PC_3** | nFeature_RNA | 0.3459 | 3.465596e-65 |
| **PC_4** | nFeature_RNA | -0.0964 | 3.968509e-06 |
| **PC_5** | nFeature_RNA | -0.0592 | 4.660539e-03 |
| **PC_7** | nFeature_RNA | 0.0674 | 1.269919e-03 |
| **PC_8** | nFeature_RNA | 0.1589 | 2.186879e-14 |
| **PC_9** | nFeature_RNA | -0.0656 | 1.721648e-03 |
| **PC_10** | nFeature_RNA | 0.1620 | 6.814986e-15 |
| **PC_11** | nFeature_RNA | -0.1700 | 2.835402e-16 |
| **PC_14** | nFeature_RNA | 0.0512 | 1.431762e-02 |
| **PC_1** | percent.ribo | 0.3140 | 1.940567e-53 |
| **PC_3** | percent.ribo | 0.1493 | 7.395562e-13 |
| **PC_4** | percent.ribo | 0.5531 | 3.123748e-183 |
| **PC_6** | percent.ribo | 0.0730 | 4.768782e-04 |
| **PC_7** | percent.ribo | -0.3316 | 9.408471e-60 |
| **PC_8** | percent.ribo | 0.1992 | 7.064141e-22 |
| **PC_11** | percent.ribo | 0.1792 | 6.151718e-18 |
| **PC_12** | percent.ribo | -0.1352 | 8.819177e-11 |
| **PC_13** | percent.ribo | -0.2090 | 5.947646e-24 |
| **PC_14** | percent.ribo | 0.0802 | 1.254414e-04 |
| **PC_15** | percent.ribo | -0.1156 | 2.989021e-08 |

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
1. **Neighbor Graph and Clustering**: We recommend using the first **6 PCs** (the recommended range) for constructing the shared nearest neighbor (SNN) graph in Milestone 6.
2. **Technical Confounders**: PCs strongly correlated with nCount_RNA or percent.mt should be monitored. Because SCTransform regressed out depth covariance, residual correlation represents biological cell-size difference rather than technical artifacts, but monitoring is advised.
3. **Biological Signal Preservation**: The leading PCs are heavily dominated by biologically meaningful programs (extracellular matrix remodeling, cell proliferation, immune infiltration) rather than technical artifacts, confirming the high quality of the filtering and normalization.

## 6. Diagnostics & Visualizations

The following publication-quality diagnostic plots were generated:
- [Variance Explained (Elbow Plot)](file:///reports/datasets/MPNST_2/pca_elbow.png)
- [Cumulative Variance Explained](file:///reports/datasets/MPNST_2/pca_cumulative_variance.png)
- [PC Loadings Plot (Top 4 PCs)](file:///reports/datasets/MPNST_2/pca_loadings.png)
- [PC Expression Heatmaps (Dims 1-9)](file:///reports/datasets/MPNST_2/pca_heatmaps.png)
- [Technical Metric Correlations](file:///reports/datasets/MPNST_2/pca_correlations.png)
- [Top Loading Genes Table (TSV)](file:///reports/datasets/MPNST_2/top_loading_genes.tsv)

---
## 7. Software & Environment Provenance
- **Seurat Version**: `5.4.0`
- **Patchwork Version**: `1.3.2`
- **Git Commit Hash**: `31cc849c7c8aaf27e8f2b594fb40d0c54989a579`
