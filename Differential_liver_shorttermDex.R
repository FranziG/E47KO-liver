library("DESeq2")
library("AnnotationDbi")
library("org.Mm.eg.db")
library("biomaRt")
library("ggrepel")
library("ggplot2")
library("dplyr")

# data preprocessing: make count table. The raw count data is diposited at GEO under the accession GSE.

# 1. read STAR count files: 1st column: gene ID (UCSC), 2nd column: counts for both strands, 3rd: 5´-->3´counts; 4th: 3`-->5`counts

GAR1221<-read.table("./liver_Dex_E47KO_554_GAR1221_mm10_STARReadsPerGene.out.tab", header=FALSE, sep="\t")
GAR1221<-GAR1221[c(5:(nrow(GAR1221))),]
GAR1220<-read.table("./liver_Dex_E47wt_553_GAR1220_mm10_STARReadsPerGene.out.tab", header=FALSE, sep="\t")
GAR1220<-GAR1220[c(5:(nrow(GAR1220))),]
GAR1219<-read.table("./liver_Dex_E47KO_551_GAR1219_mm10_STARReadsPerGene.out.tab", header=FALSE, sep="\t")
GAR1219<-GAR1219[c(5:(nrow(GAR1219))),]
GAR1217<-read.table("./liver_Dex_E47wt_544_GAR1217_mm10_STARReadsPerGene.out.tab", header=FALSE, sep="\t")
GAR1217<-GAR1217[c(5:(nrow(GAR1217))),]
GAR1216<-read.table("./liver_Dex_E47wt_543_GAR1216_mm10_STARReadsPerGene.out.tab", header=FALSE, sep="\t")
GAR1216<-GAR1216[c(5:(nrow(GAR1216))),]
GAR1218<-read.table("./liver_Dex_E47KO_545_GAR1218_mm10_STARReadsPerGene.out.tab", header=FALSE, sep="\t")
GAR1218<-GAR1218[c(5:(nrow(GAR1218))),]

#2. generate one table for raw countdata; We used the counts ignoring strandedness as not all RNA-Seq samples of this study were stranded.

counts<-cbind(GAR1216[,c(1,2)],GAR1217$V2,GAR1220$V2,GAR1218$V2,GAR1219$V2,GAR1221$V2)
rownames(counts)<-counts[,1]
counts<-counts[,-1]
colnames(counts)<-c("GAR1216","GAR1217","GAR1220","GAR1218","GAR1219","GAR1221")
#filter lowest 25% counts (rowmean below )
#remove below 6
#counts<-counts[rowMeans(counts)>5,]
#quantile(counts$GAR1221,.25)
#1216: 36; 1217: 40; 1218:38; 1219: 33; 1220: 44; 1221: 42 => take lowest as threshold: 33
counts<-counts[rowMeans(counts)>33,]

#3. prepare metadata

genotype<-c("wt","wt","wt","E47KO","E47KO","E47KO")
samples<-colnames(counts)
pdata<-as.data.frame(cbind(samples,genotype))

#4. Differential expression analysis using DESeq2

cds<-DESeqDataSetFromMatrix(countData=counts,colData=pdata,design=~genotype)
#set wildtype as basal level
cds$genotype <- relevel(cds$genotype, "wt")
cds<-DESeq(cds)

# get differential expression wt versus E47KO
res<-as.data.frame(results(cds,contrast=c("genotype","E47KO","wt")))

#5. annotated gene IDs

#load annotation DB
mart <- useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")
#annotate results object
res$ensembl<-rownames(res)
symbol<-getBM(c("mgi_symbol","entrezgene","ensembl_gene_id"),"ensembl_gene_id", res$ensembl, mart)
res<-merge(res,symbol,by.x="ensembl",by.y="ensembl_gene_id")
#add normalized counts to results table
counts2<-DESeq2::counts(cds,normalized=TRUE)
counts2<-cbind(counts2,ensembl=rownames(counts2))
res<-merge(counts2,res,by="ensembl")
#remove genes that have NA for log2FC and adjP
res_f<-res[!is.na(res$log2FoldChange)&!is.na(res$padj),]
#sort differential expressed genes by log2FC in descending order
res_f<-arrange(res_f,desc(abs(log2FoldChange)))
#remove duplicated entries (By sorting data first, we remove the once with lower fold change)
res_f<-res_f[!duplicated(res_f$mgi_symbol),]

# write results to table. These data is published as processed data along with the raw counts on the NCBI Gene Omnibus with the accession number GSE.

write.table(res_f, file="./2018-3_liver_E47Kovswt_shorttermDex_unfiltered.txt",sep="\t",row.names = FALSE)

#6. Vulcano plot to display differential expression

#add color flag for genes significantly (P<0.05) changing more the 1.3x
res_f$color_flag <- factor(ifelse(res_f$log2FoldChange>0.38&res_f$pvalue<0.05, 1, ifelse(res_f$log2FoldChange<(-0.38)&res_f$pvalue<0.05,2,0)))
#plot; labels are added for genes showing significant (adjP<0.05) fold-change of minimal 1.5-fold
ggplot(data=res_f, aes(x=log2FoldChange, y=-log10(pvalue),col=color_flag)) +
  geom_point(alpha=0.7, size=5) +
  xlim(c(-4, 4)) + ylim(c(0, 10)) +
  xlab("log2 fold change") + ylab("-log10 p-value")+
  geom_text_repel(aes(label=ifelse(abs(log2FoldChange)>0.58&padj<0.05, mgi_symbol,"")), size=3)+
  theme_bw()+
  scale_color_manual(values=c("grey","red","blue"))+
  theme(legend.position = "none",panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
