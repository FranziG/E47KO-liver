library("DESeq2")
library("AnnotationDbi")
library("org.Mm.eg.db")
library("biomaRt")
library("ggrepel")
library("ggplot2")
library("dplyr")

# data preprocessing: make count table. The raw count data is diposited at GEO under the accession GSE.

# 1. read STAR count files: 1st column: gene ID (UCSC), 2nd column: counts for both strands, 3rd: 5´-->3´counts; 4th: 3`-->5`counts
E47KO_Dex_GAR0179<-read.table("./liver_Dex_E47KO_215_GAR0179_mm10_STAR_count.txt", header=FALSE, sep="\t",stringsAsFactors = FALSE)
E47KO_Dex_GAR0179<-E47KO_Dex_GAR0179[c(1:(nrow(E47KO_Dex_GAR0179)-5)),]
E47KO_Dex_GAR0181<-read.table("./liver_Dex_E47KO_230_GAR0181_mm10_STAR_count.txt", header=FALSE, sep="\t",stringsAsFactors = FALSE)
E47KO_Dex_GAR0181<-E47KO_Dex_GAR0181[c(1:(nrow(E47KO_Dex_GAR0181)-5)),]
E47KO_Dex_GAR0183<-read.table("./liver_Dex_E47KO_382_GAR0183_mm10_STAR_count.txt", header=FALSE, sep="\t",stringsAsFactors = FALSE)
E47KO_Dex_GAR0183<-E47KO_Dex_GAR0183[c(1:(nrow(E47KO_Dex_GAR0183)-5)),]
E47wt_Dex_GAR0182<-read.table("./liver_Dex_wt_231_GAR0182_mm10_STAR_count.txt", header=FALSE, sep="\t",stringsAsFactors = FALSE)
E47wt_Dex_GAR0182<-E47wt_Dex_GAR0182[c(1:(nrow(E47wt_Dex_GAR0182)-5)),]
E47wt_Dex_GAR0184<-read.table("./liver_Dex_wt_384_GAR0184_mm10_STAR_count.txt", header=FALSE, sep="\t",stringsAsFactors = FALSE)
E47wt_Dex_GAR0184<-E47wt_Dex_GAR0184[c(1:(nrow(E47wt_Dex_GAR0184)-5)),]

#2. generate one table for raw countdata; We used the counts ignoring strandedness as not all RNA-Seq samples of this study were stranded.

counts<-cbind(E47KO_Dex_GAR0179,E47KO_Dex_GAR0181$V2,E47KO_Dex_GAR0183$V2,E47wt_Dex_GAR0182$V2,E47wt_Dex_GAR0184$V2)
rownames(counts)<-counts[,1]
counts<-counts[,-1]
colnames(counts)<-c("E47KO_Dex_GAR0179","E47KO_Dex_GAR0181","E47KO_Dex_GAR0183","E47wt_Dex_GAR0182","E47wt_Dex_GAR0184")
counts[is.na(counts)]<-0
#counts_all<-counts_all[rowMeans(counts_all)>5,]
#filter lowest 25% counts (rowmean below )
#quantile(counts_all$E47wt_Dex_GAR0184,.25)
#179:22; 181:25; 183: 17; 182:22; 184:15 => take lowest as threshold: 3x15
counts<-counts[rowMeans(counts)>45,]

#3. prepare metadata

genotype<-c("KO","KO","KO","wt","wt")
samples<-colnames(counts)
pdata<-as.data.frame(cbind(samples,genotype))

#4. Differential expression analysis using DESeq2

cds<-DESeqDataSetFromMatrix(countData=counts,colData=pdata,design=~genotype)
#set wildtype as basal level
cds$genotype <- relevel(cds$genotype, "wt")
cds<-DESeq(cds)

# get differential expression wt versus E47KO
res_Dex<-data.frame(results(cds,contrast=c("genotype","KO","wt")),stringsAsFactors = FALSE)

#5. annotated gene IDs

#load annotation DB
mart <- useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")
#annotate results object
res_Dex$ensembl<-rownames(res_Dex)
symbol<-getBM(c("mgi_symbol","ensembl_gene_id","entrezgene"),"ensembl_gene_id",res_Dex$ensembl, mart)
res_Dex<-merge(res_Dex,symbol,by.x="ensembl",by.y="ensembl_gene_id")
#add normalized counts to results table
counts_Dex2<-data.frame(DESeq2::counts(cds,normalized=TRUE),stringsAsFactors = FALSE)
counts_Dex2<-cbind(counts_Dex2,ensembl=rownames(counts_Dex2))
res_Dex<-merge(counts_Dex2,res_Dex,by="ensembl")
#remove genes that have NA for log2FC and adjP
res_Dex_f<-res_Dex[!is.na(res_Dex$log2FoldChange)&!is.na(res_Dex$padj),]
#sort differential expressed genes by log2FC in descending order
res_Dex_f<-arrange(res_Dex_f,desc(abs(log2FoldChange)))
#remove duplicated entries (By sorting data first, we remove the once with lower fold change)
res_Dex_f<-res_Dex_f[!duplicated(res_Dex_f$mgi_symbol),]

# write results to table. These data is published as processed data along with the raw counts on the NCBI Gene Omnibus with the accession number GSE.

write.table(res_Dex_f,file="./liver_totalE47KOvsWt_Dex_unfiltered.txt",row.names = FALSE,sep="\t")

#6. Vulcano plot to display differential expression

#add color flag for genes significantly (P<0.05) changing more the 1.3x
res_Dex_f$color_flag <- factor(ifelse(res_Dex_f$log2FoldChange>0.38&res_Dex_f$pvalue<0.05, 1, ifelse(res_Dex_f$log2FoldChange<(-0.38)&res_Dex_f$pvalue<0.05,2,0)))
#plot
ggplot(data=res_Dex_f, aes(x=log2FoldChange, y=-log10(pvalue),col=color_flag)) +
  geom_point(alpha=0.7, size=4) +
  xlim(c(-5, 3)) + ylim(c(0, 12)) +
  xlab("log2 fold change") + ylab("-log10 p-value")+
  geom_text_repel(aes(label=ifelse(abs(log2FoldChange)>0.56&padj<0.1, mgi_symbol,"")), size=5)+
  theme_bw()+
  scale_color_manual(values=c("grey","red","blue"))+
  theme(legend.position = "none",panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
