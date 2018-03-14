# The role of the transcription factor E47 in modulating glucocorticoid receptor
## Bioinformatic analysis for RNA-Seq and ChIP-Seq performed in Hemmer et al. 2018

### Study design

The aim of this study was to analyze the interplay of the transcription factors E47, a splice variant of E2A, and glucocorticoid receptor (GR) in the metabolic fuction of glucocorticoid signalling. To address the research aim, we performen Chromatin-Immunoprecipitation for the transcription factors E2A and GR in untreated liver tissue (1). For functional studies on E47 loss-of-function mice, we performed RNA-Seq of metabolic tissues relevent for GR signalling: epididymal white adipose tissue (WAT), skeletal muscle (SM, quadriceps) and liver. For WAT (2) and SM (3) the samples from wildtype and total E47 knock-out (E47KO) where treated with dexamethason for 3 weeks via the drinking water and killed after 16hrs fasting and GTT. For liver, we compared wildtype and total E47KO male mice after long-term (3weeks) dexamethason (4) or corticosterone treatment (5) via the drinking water. These animals were sacrificed after 16hrs fasting and GTT as well. A small study on short-term (1h) dexamethason effect on wildtype and liver-specific E47KO male mice (Albcre/+;E47fl/fl) was included to identify direct GR target genes (6). Without treatment and with ad libidum food supply, expression in wildtype and total E47KO liver tissue was analysed to exclude any basal and developmental differences in E47KO livers (7).

### data accessibility via GEO
    (1) ChIP-Seq: GSE111526
    (2)-(7) RNA-Seq: 

### Files

(1) Commands used for mapping, peak calling and motif enrichment of ChIP-Seq data

    1. 2018_3_Hemmer_ChIP-Seq_E2A_bash.txt
    2. 2018_3_Hemmer_ChIP-Seq_GR_bash.txt

(1) R scripts used to analyze the ChIP-Seq data with the corresponding R session information

    1. reprodPeaks.R and RsessionInfo_reprodPeaks.txt
    
(2)-(7) Commands used for mapping and counting of RNA-Seq data

    (2) 2018-3-Hemmer_RNA-Seq_WAT_bash.txt
    (3) 2018-3-Hemmer_RNA-Seq_SM_bash.txt
    (4) 2018-3-Hemmer_RNA-Seq_longtermDex_commands.txt
    (5) 2018-3-Hemmer_RNA-Seq_longtermCort_commands.txt
    (6) 2018-3-Hemmer_RNA-Seq_shorttermDex_commands.txt
    (7) 2017-3-Hemmer_RNA-Seq_untreated_bash.txt 

(2)-(7) R scripts used to analyze the RNA-Seq data with the corresponding R session information
    
    (2) Differential_WAT.R and Differential_WAT_RsessionInfo.txt
    (3) Differential_SM.R and Differential_SM_RsessionInfo.txt
    (4) Differential_liver_longtermDex.R and Differential_liver_longtermDex_RsessionInfo.txt
    (5) Differential_liver_longtermCort.R and Differential_liver_longtermCort_RsessionInfo.txt
    (6) Differential_liver_shorttermDex.R and Differential_liver_shorttermDex_RsessionInfo.txt
    (7) Differential_untreated.R and Differential_untreated_RsessionInfo.txt
