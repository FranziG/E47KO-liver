library("DESeq2")
library("AnnotationDbi")
library("org.Mm.eg.db")
library("biomaRt")
library("ggrepel")
library("ggplot2")
library("dplyr")

# data preprocessing: make count table. The raw count data is diposited at GEO under the accession GSE.

# 1. read STAR count files: 1st column: gene ID (UCSC), 2nd column: counts for both strands, 3rd: 5´-->3´counts; 4th: 3`-->5`counts
GAR293<-read.table("./WAT_Dex_E47KO_694_GAR0293_mm10_STAR_count.txt", header=FALSE, sep="\t")
GAR293<-GAR293[c(1:(nrow(GAR293)-5)),]
GAR295<-read.table("./WAT_Dex_E47KO_755_GAR0295_mm10_STAR_count.txt", header=FALSE, sep="\t")
GAR295<-GAR295[c(1:(nrow(GAR295)-5)),]
GAR292<-read.table("./WAT_Dex_wt_693_GAR0292_mm10_STAR_count.txt", header=FALSE, sep="\t")
GAR292<-GAR292[c(1:(nrow(GAR292)-5)),]
GAR294<-read.table("./WAT_Dex_wt_753_GAR0294_mm10_STAR_count.txt", header=FALSE, sep="\t")
GAR294<-GAR294[c(1:(nrow(GAR294)-5)),]

#2. generate one table for raw countdata; We used the counts ignoring strandedness as not all RNA-Seq samples of this study were stranded.

counts<-cbind(GAR293,GAR295$V2,GAR292$V2,GAR294$V2)
rownames(counts)<-counts[,1]
counts<-counts[,-1]
colnames(counts)<-c("GAR293","GAR295","GAR292","GAR294")
#filter lowest 25% counts (rowmean below )
#counts<-counts[rowMeans(counts)>5,]
#quantile(counts$GAR296,.25)
# 292: 18; 293:20; 294: 18; 295: 25=> take lowest as threshold: 18
counts<-counts[rowMeans(counts)>18,]

#3. prepare metadata

genotype<-c("E47KO","E47KO","wt","wt")
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
counts2<-DESeq2::counts(cds,normalized=TRUE)
counts2<-cbind(counts2,ensembl=rownames(counts2))
res<-merge(counts2,res,by="ensembl")
res_f<-res[!is.na(res$log2FoldChange)&!is.na(res$padj),]
res_f<-arrange(res_f,desc(abs(log2FoldChange)))
res_f<-res_f[!duplicated(res_f$mgi_symbol),]

# write results to table. These data is published as processed data along with the raw counts on the NCBI Gene Omnibus with the accession number GSE.

write.table(res_f, file="./2018_3_WAT_E47Kovswt_Dex_DE_unfiltered.txt",sep="\t",row.names = FALSE)

#6. Vulcano plot to display differential expression

#add color flag for genes significantly (P<0.05) changing more the 1.3-fold
res_f$color_flag <- factor(ifelse(res_f$log2FoldChange>0.38&res_f$pvalue<0.05, 1, ifelse(res_f$log2FoldChange<(-0.38)&res_f$pvalue<0.05,2,0)))
#plot; labels are added for genes showing significant (adjP<0.05) fold-change of minimal 1.5-fold
ggplot(data=res_f, aes(x=log2FoldChange, y=-log10(pvalue),col=color_flag)) +
  geom_point(alpha=0.7, size=4) +
  xlim(c(-5.2, 9.1)) + ylim(c(0, 17)) +
  xlab("log2 fold change") + ylab("-log10 p-value")+
  geom_text_repel(aes(label=ifelse(abs(log2FoldChange)>0.58&padj<0.05, mgi_symbol,"")), size=3)+
  theme_bw()+
  scale_color_manual(values=c("grey","red","blue"))+
  theme(legend.position = "none",panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

