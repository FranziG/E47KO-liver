# The role of the transcription factor E47 in modulating glucocorticoid receptor
## Bioinformatic analysis for RNA-Seq and ChIP-Seq performed in Hemmer et al. 2018

### Study design

The aim of this study was to analyze the interplay of the transcription factors E47, a splice variant of E2A, and glucocorticoid receptor (GR) in the metabolic fuction of glucocorticoid signalling. To address the research aim, we performen Chromatin-Immunoprecipitation for the transcription factors E2A and GR in untreated liver tissue (1). For functional studies on E47 loss-of-function mice, we performed RNA-Seq of metabolic tissues relevent for GR signalling: epididymal white adipose tissue (WAT), skeletal muscle (SM, quadriceps) and liver. For WAT (2) and SM (3) the samples from wildtype and total E47 knock-out (E47KO) where treated with dexamethason for 3 weeks via the drinking water and killed after 16hrs fasting and GTT. For liver, we compared wildtype and total E47KO male mice after long-term (3weeks) dexamethason (4) or corticosterone treatment (5) via the drinking water. These animals were sacrificed after 16hrs fasting and GTT as well. A small study on short-term (1h) dexamethason effect on wildtype and liver-specific E47KO male mice (Albcre/+;E47fl/fl) was included to identify direct GR target genes (7). Without treatment and with ad libidum food supply, expression in wildtype and total E47KO liver tissue was analysed to exclude any basal and developmental differences in E47KO livers (8).

# data accessibility via GEO
    (1) ChIP-Seq: GSE111526
    (2)-(8) RNA-Seq: 

### Files

(1) bash scripts used for mapping, peak calling and motif enrichment of ChIP-Seq data
    1. 2018_3_Hemmer_ChIP-Seq_E2A_bash.txt
    2. 2018_3_Hemmer_ChIP-Seq_GR_bash.txt

(1) R scripts used to analyze the ChIP-Seq data with the corresponding R session information
    1. reprodPeaks.R and RsessionInfo_reprodPeaks.txt
    
- bash scripts used for mapping and counting of RNA-Seq data
    (8) 2017-3-Hemmer_RNA-Seq_untreated_bash.txt 

- R scripts used to analyze the ChIP-Seq data with the corresponding R session information
