library("GenomicRanges")
library(VennDiagram)

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
bed2GRanges2 <-function(peaks)
{
  myrange <- GRanges(seqnames=peaks$seqnames,range=IRanges(start=peaks$start, end=peaks$end, names=paste(peaks$seqnames,peaks$start,sep=":")),
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
E47wt_Gr_rep1<-macs2GRanges(E47wt_Gr_rep1)
E47wt_Gr_rep2<-macs2GRanges(E47wt_Gr_rep2)

#3.remove blacklisted regions

E47wt_E2ASC_rep1_bl<-subsetByOverlaps(E47wt_E2ASC_rep1,blacklist)
E47wt_E2ASC_rep1<-GenomicRanges::setdiff(E47wt_E2ASC_rep1,E47wt_E2ASC_rep1_bl,ignore.strand=TRUE)
E47wt_E2ASC_rep2_bl<-subsetByOverlaps(E47wt_E2ASC_rep2,blacklist)
E47wt_E2ASC_rep2<-GenomicRanges::setdiff(E47wt_E2ASC_rep2,E47wt_E2ASC_rep2_bl,ignore.strand=TRUE)
E47wt_Gr_rep1_bl<-subsetByOverlaps(E47wt_Gr_rep1,blacklist)
E47wt_Gr_rep1<-GenomicRanges::setdiff(E47wt_Gr_rep1,E47wt_Gr_rep1_bl,ignore.strand=TRUE)
E47wt_Gr_rep2_bl<-subsetByOverlaps(E47wt_Gr_rep2,blacklist)
E47wt_Gr_rep2<-GenomicRanges::setdiff(E47wt_Gr_rep2,E47wt_Gr_rep2_bl,ignore.strand=TRUE)

#4. Merge E2A replicates to generate a E2A "peak universe" in wt liver tissue

E47_wtUniverse<-rbind(data.frame(E47wt_E2ASC_rep1),data.frame(E47wt_E2ASC_rep2))
E47_wtUniverse<-bed2GRanges2(E47_wtUniverse)
E47_wtUniverse<-reduce(E47_wtUniverse)

#5. Merge GR replicates to generate a GR "peak universe" in wt liver tissue

E47wt_Gruniverse<-rbind(data.frame(E47wt_Gr_rep1,stringsAsFactors=FALSE),data.frame(E47wt_Gr_rep2,stringsAsFactors=FALSE))
E47wt_Gruniverse<-bed2GRanges2(E47wt_Gruniverse)
E47wt_Gruniverse<-reduce(E47wt_Gruniverse)

#E2A and Gr overlap in untreated wt liver by Venn diagram

draw.pairwise.venn(area1 = length(E47_wtUniverse),area2 = length(E47wt_Gruniverse),cross.area = length(subsetByOverlaps(E47_wtUniverse,E47wt_Gruniverse)),fill = c("Red","DodgerBlue2"), col=c("Red","DodgerBlue2"),cex=2.5, cat.cex = 1.5,alpha = 0.8, category = c("E2A","GR"))
