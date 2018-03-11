
library("GenomicRanges")
library ("ChIPpeakAnno")
library("ChIPseeker")
library("GenomicFeatures")

#function to convert MACS peaks into GenomiRanges object
macs2GRanges <-function(peaks)
{
  myrange <- GRanges(seqnames=peaks$chr,range=IRanges(start=peaks$start, end=peaks$end, names=paste(peaks$chr,peaks$start,sep=":")),
                     strand="*",
                     count=peaks$pileup,
                     score=peaks$X.log10.pvalue.,
                     FE=peaks$fold_enrichment,
                     fdr=peaks$X.log10.qvalue.#,
                     #maxpos=peaks$abs_summit
  )
  return(myrange)
}
#function to convert BED file into GenomicRanges object
bed2GRanges <-function(peaks)
{
  myrange <- GRanges(seqnames=peaks$V1,range=IRanges(start=peaks$V2, end=peaks$V3, names=paste(peaks$V1,peaks$V2,sep=":")),
                     strand="*"
  )
  return(myrange)
}

#mouse chromosomes
chr<-c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chrX","chrY")

#get blacklisted regions: blacklisted regions downloaded from http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/mm10-mouse/mm10.blacklist.bed.gz
blacklist <- read.table("./blacklist_mm10.bed", header=FALSE, sep='\t', stringsAsFactors=FALSE)
blacklist<-bed2GRanges(blacklist)

#1.load peaks

#E2A ChIP-Seq in untreated livers from male, wildtype C57BL/6 mice

E47wt_E2ASC_rep1<-read.table("./liver_wt_untreated_No58_aE2A_SCFreiburg_rep1_GAR1607_subsampled_mm10_MACS2_PE_FDR005_peaks.xls", header=TRUE, sep='\t', stringsAsFactors=FALSE)
E47wt_E2ASC_rep1<-E47wt_E2ASC_rep1[E47wt_E2ASC_rep1$chr%in%chr,]
#use 10th percentile as threshold for FE
E47wt_E2ASC_rep1<-E47wt_E2ASC_rep1[E47wt_E2ASC_rep1$fold_enrichment>=quantile(E47wt_E2ASC_rep1$fold_enrichment,.1),]

E47wt_E2ASC_rep2<-read.table("./liver_wt_untreated_No670_aE2A_SCFreiburg_rep2_GAR1608_subsampled_mm10_MACS2_PE_FDR005_peaks.xls", header=TRUE, sep='\t', stringsAsFactors=FALSE)
E47wt_E2ASC_rep2<-E47wt_E2ASC_rep2[E47wt_E2ASC_rep2$chr%in%chr,]
#use 10th percentile as threshold for FE
E47wt_E2ASC_rep2<-E47wt_E2ASC_rep2[E47wt_E2ASC_rep2$fold_enrichment>=quantile(E47wt_E2ASC_rep2$fold_enrichment,.1),]

#GR ChIP-Seq in untreated livers from male, wildtype C57BL/6 mice

E47KO_Gr_rep1<-read.table("./liver_E47KO_untreated_No673_aGRPT_rep1_GAR1594_subsampled_mm10_MACS2_PE_FDR005_peaks.xls", header=TRUE, sep='\t', stringsAsFactors=FALSE)
E47KO_Gr_rep1<-E47KO_Gr_rep1[E47KO_Gr_rep1$chr%in%chr,]
#use 10th percentile as threshold for FE
E47KO_Gr_rep1<-E47KO_Gr_rep1[E47KO_Gr_rep1$fold_enrichment>=quantile(E47KO_Gr_rep1$fold_enrichment,.1),]

E47KO_Gr_rep2<-read.table("./liver_E47KO_untreated_No55_aGRPT_rep2_GAR1593_subsampled_mm10_MACS2_PE_FDR005_peaks.xls", header=TRUE, sep='\t', stringsAsFactors=FALSE)
E47KO_Gr_rep2<-E47KO_Gr_rep2[E47KO_Gr_rep2$chr%in%chr,]
#use 10th percentile as threshold for FE
E47KO_Gr_rep2<-E47KO_Gr_rep2[E47KO_Gr_rep2$fold_enrichment>=quantile(E47KO_Gr_rep2$fold_enrichment,.1),] 

E47wt_Gr_rep1<-read.table("./liver_wt_untreated_No670_aGRPT_rep1_GAR1589_subsampled_mm10_MACS2_PE_FDR005_peaks.xls", header=TRUE, sep='\t', stringsAsFactors=FALSE)
E47wt_Gr_rep1<-E47wt_Gr_rep1[E47wt_Gr_rep1$chr%in%chr,]
#use 10th percentile as threshold for FE
E47wt_Gr_rep1<-E47wt_Gr_rep1[E47wt_Gr_rep1$fold_enrichment>=quantile(E47wt_Gr_rep1$fold_enrichment,.1),]

E47wt_Gr_rep2<-read.table("./liver_wt_untreated_NoM1_aGrPT_rep2_GAR1588_subsampled_mm10_MACS2_PE_FDR005_peaks.xls", header=TRUE, sep='\t', stringsAsFactors=FALSE)
E47wt_Gr_rep2<-E47wt_Gr_rep2[E47wt_Gr_rep2$chr%in%chr,]
#use 10th percentile as threshold for FE
E47wt_Gr_rep2<-E47wt_Gr_rep2[E47wt_Gr_rep2$fold_enrichment>=quantile(E47wt_Gr_rep2$fold_enrichment,.1),]

#2.convert to GenomicRanges

E47wt_E2ASC_rep1<-macs2GRanges(E47wt_E2ASC_rep1)
E47wt_E2ASC_rep2<-macs2GRanges(E47wt_E2ASC_rep2)
E47KO_Gr_rep1<-macs2GRanges(E47KO_Gr_rep1)
E47KO_Gr_rep2<-macs2GRanges(E47KO_Gr_rep2)
E47wt_Gr_rep1<-macs2GRanges(E47wt_Gr_rep1)
E47wt_Gr_rep2<-macs2GRanges(E47wt_Gr_rep2)

#3.remove blacklisted regions

E47wt_E2ASC_rep1_bl<-subsetByOverlaps(E47wt_E2ASC_rep1,blacklist)
E47wt_E2ASC_rep1<-GenomicRanges::setdiff(E47wt_E2ASC_rep1,E47wt_E2ASC_rep1_bl,ignore.strand=TRUE)
E47wt_E2ASC_rep2_bl<-subsetByOverlaps(E47wt_E2ASC_rep2,blacklist)
E47wt_E2ASC_rep2<-GenomicRanges::setdiff(E47wt_E2ASC_rep2,E47wt_E2ASC_rep2_bl,ignore.strand=TRUE)
E47KO_Gr_rep1_bl<-subsetByOverlaps(E47KO_Gr_rep1,blacklist)
E47KO_Gr_rep1<-GenomicRanges::setdiff(E47KO_Gr_rep1,E47KO_Gr_rep1_bl,ignore.strand=TRUE)
E47KO_Gr_rep2_bl<-subsetByOverlaps(E47KO_Gr_rep2,blacklist)
E47KO_Gr_rep2<-GenomicRanges::setdiff(E47KO_Gr_rep2,E47KO_Gr_rep2_bl,ignore.strand=TRUE)
E47wt_Gr_rep1_bl<-subsetByOverlaps(E47wt_Gr_rep1,blacklist)
E47wt_Gr_rep1<-GenomicRanges::setdiff(E47wt_Gr_rep1,E47wt_Gr_rep1_bl,ignore.strand=TRUE)
E47wt_Gr_rep2_bl<-subsetByOverlaps(E47wt_Gr_rep2,blacklist)
E47wt_Gr_rep2<-GenomicRanges::setdiff(E47wt_Gr_rep2,E47wt_Gr_rep2_bl,ignore.strand=TRUE)

#4. determine overlapping peak ranges = reproducible peaks

E47wt_E2ASC<-subsetByOverlaps(E47wt_E2ASC_rep1,E47wt_E2ASC_rep2,type = "any",minoverlap = .2, maxgap = 50)
E47KO_Gr<-subsetByOverlaps(E47KO_Gr_rep1,E47KO_Gr_rep2,type = "any",minoverlap = .2, maxgap = 50)
E47wt_Gr<-subsetByOverlaps(E47wt_Gr_rep1,E47wt_Gr_rep2,type = "any",minoverlap = .2, maxgap=50)

#5. annotate peak to closest TSS

#load anntoation information for mm10
mm10KG<-makeTxDbFromUCSC(genome="mm10",tablename="knownGene",transcript_ids=NULL,circ_seqs=DEFAULT_CIRC_SEQS,url="http://genome-euro.ucsc.edu/cgi-bin/",goldenPath_url="http://hgdownload.cse.ucsc.edu/goldenPath",miRBaseBuild="GRCm38")

#annotate and write bed files for reproducible peaks (used for HOMER motif enrichment and GREAT functional enrichment analysis); promoter region is defined as -1000 to +500 around TSS

peakAnno<-as.data.frame(annotatePeak(E47wt_E2ASC,tssRegion=c(-1000,500),TxDb=mm10KG, level=c("transcript","gene"), genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron", "Intergenic"),annoDb="org.Mm.eg.db"))
write.table(peakAnno, file="./liver_wt_untreated_aE2A_SCFreiburg_OverlapReplicates_mm10_MACS2_PE_FDR005_peaks_annotated.tab", sep="\t",row.names = FALSE)
E47wt_E2ASC_rep1_bed<-data.frame(cbind(as.character(peakAnno$seqnames),peakAnno$start,peakAnno$end,paste(peakAnno$SYMBOL,peakAnno$start,sep="_"),"",rep(0,nrow(peakAnno))),stringsAsFactors=FALSE)
write.table(E47wt_E2ASC_rep1_bed, "./liver_wt_untreated_No58_aE2A_SCFreiburg_OverlapReplicates_mm10_MACS2_PE_FDR005_peaks.bed",row.names = FALSE,col.names=FALSE,quote = FALSE,sep="\t")

peakAnno<-as.data.frame(annotatePeak(E47wt_Gr,tssRegion=c(-1000,500),TxDb=mm10KG, level=c("transcript","gene"), genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron", "Intergenic"),annoDb="org.Mm.eg.db"))
write.table(peakAnno, file="./liver_wt_untreated_aGr_overlap_PE_FDR005_peaks_annotated.tab", sep="\t",row.names = FALSE)
E47wt_Gr_bed<-data.frame(cbind(as.character(peakAnno$seqnames),peakAnno$start,peakAnno$end,paste(peakAnno$SYMBOL,peakAnno$start,sep="_"),"",rep(0,nrow(peakAnno))),stringsAsFactors=FALSE)
write.table(E47wt_Gr_bed, "./liver_wt_untreated_aGr_overlap_PE_FDR005_peaks.bed",row.names = FALSE,col.names=FALSE,quote = FALSE,sep="\t")

peakAnno<-as.data.frame(annotatePeak(E47KO_Gr,tssRegion=c(-1000,500),TxDb=mm10KG, level=c("transcript","gene"), genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron", "Intergenic"),annoDb="org.Mm.eg.db"))
write.table(peakAnno, file="./liver_E47KO_untreated_aGr_overlap_PE_FDR005_peaks_annotated.tab", sep="\t",row.names = FALSE)
E47KO_Gr_bed<-data.frame(cbind(as.character(peakAnno$seqnames),peakAnno$start,peakAnno$end,paste(peakAnno$SYMBOL,peakAnno$start,sep="_"),"",rep(0,nrow(peakAnno))),stringsAsFactors=FALSE)
write.table(E47KO_Gr_bed, "./liver_E47KO_untreated_aGr_overlap_PE_FDR005_peaks.bed",row.names = FALSE,col.names=FALSE,quote = FALSE,sep="\t")


