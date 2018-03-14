library("DESeq2")
library("AnnotationDbi")
library("org.Mm.eg.db")
library("biomaRt")
library("ggrepel")
library("ggplot2")
library("dplyr")

# data preprocessing: make count table. The raw count data is diposited at GEO under the accession GSE.

# 1. read STAR count files: 1st column: gene ID (UCSC), 2nd column: counts for both strands, 3rd: 5´-->3´counts; 4th: 3`-->5`counts
GAR1479<-read.table("./liver_E47KO_untreated_rep4_1630_GAR1479_mm10_STARReadsPerGene.out.tab", header=FALSE, sep="\t")
GAR1479<-GAR1479[c(5:(nrow(GAR1479))),]
GAR1480<-read.table("./liver_E4wt_untreated_rep4_1632_GAR1480_mm10_STARReadsPerGene.out.tab", header=FALSE, sep="\t")
GAR1480<-GAR1480[c(5:(nrow(GAR1480))),]
GAR1482<-read.table("./liver_E47wt_untreated_rep2_111_GAR1482_mm10_STARReadsPerGene.out.tab", header=FALSE, sep="\t")
GAR1482<-GAR1482[c(5:(nrow(GAR1482))),]
GAR1568<-read.table("./liver_E47KO_untreated_No55_GAR1568_mm10_STARReadsPerGene.out.tab", header=FALSE, sep="\t")
GAR1568<-GAR1568[c(5:(nrow(GAR1568))),]
GAR1571<-read.table("./liver_E47KO_untreated_No672_GAR1571_mm10_STARReadsPerGene.out.tab", header=FALSE, sep="\t")
GAR1571<-GAR1571[c(5:(nrow(GAR1571))),]
GAR1572<-read.table("./liver_E47KO_untreated_No673_GAR1572_mm10_STARReadsPerGene.out.tab", header=FALSE, sep="\t")
GAR1572<-GAR1572[c(5:(nrow(GAR1572))),]
GAR1573<-read.table("./liver_wt_untreated_No807_GAR1573_mm10_STARReadsPerGene.out.tab", header=FALSE, sep="\t")
GAR1573<-GAR1573[c(5:(nrow(GAR1573))),]
GAR1569<-read.table("./liver_wt_untreated_No58_GAR1569_mm10_STARReadsPerGene.out.tab", header=FALSE, sep="\t")
GAR1569<-GAR1569[c(5:(nrow(GAR1569))),]
GAR1570<-read.table("./liver_wt_untreated_No670_GAR1570_mm10_STARReadsPerGene.out.tab", header=FALSE, sep="\t")
GAR1570<-GAR1570[c(5:(nrow(GAR1570))),]
GAR1574<-read.table("./liver_E47KO_untreated_No808_GAR1574_mm10_STARReadsPerGene.out.tab", header=FALSE, sep="\t")
GAR1574<-GAR1574[c(5:(nrow(GAR1574))),]
GAR1575<-read.table("./liver_E47KO_untreated_No825_GAR1575_mm10_STARReadsPerGene.out.tab", header=FALSE, sep="\t")
GAR1575<-GAR1575[c(5:(nrow(GAR1575))),]

#2. generate one table for raw countdata; We used the counts ignoring strandedness as not all RNA-Seq samples of this study were stranded.
counts_un<-cbind(GAR1479[,c(1,2)],GAR1480$V2,GAR1482$V2,GAR1568$V2,GAR1569$V2,GAR1570$V2,GAR1571$V2,GAR1572$V2,GAR1573$V2,GAR1574$V2,GAR1575$V2)
rownames(counts_un)<-counts_un[,1]
counts_un<-counts_un[,-1]
colnames(counts_un)<-c("GAR1479","GAR1480","GAR1482","GAR1568","GAR1569","GAR1570","GAR1571","GAR1572","GAR1573","GAR1574","GAR1575")
#filter lowest 25% counts (rowmean below )
#a remove below 6
#counts_un<-counts_un[rowMeans(counts_un)>5,]
#quantile(counts_un$GAR1575,.25)
#GAR1479:46 GAR1480:47 GAR1482:47 GAR1568:50 GAR1569:51 GAR1570:55 GAR1571:58 GAR1572:51 GAR1573:47 GAR1574:47  GAR1575:44
counts_un<-counts_un[rowMeans(counts_un)>44,]

#3. prepare metadata
genotype_un<-c("E47KO","wt","wt","E47KO","wt","wt","E47KO","E47KO","wt","E47KO","E47KO")
samples_un<-colnames(counts_un)
pdata_un<-as.data.frame(cbind(samples_un,genotype_un))

#4. Differential expression analysis using DESeq2

cds_un<-DESeqDataSetFromMatrix(countData=counts_un,colData=pdata_un,design=~genotype_un)
#set wildtype as basal level
cds_un$genotype_un <- relevel(cds_un$genotype_un, "wt")
cds_un<-DESeq(cds_un)

# get differential expression wt versus E47KO

res_un<-as.data.frame(results(cds_un,contrast=c("genotype_un","E47KO","wt")))

#5. annotated gene IDs

#load annotation DB
mart <- useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")
#annotate results object
res_un$ensembl<-rownames(res_un)
symbol<-getBM(c("mgi_symbol","entrezgene","ensembl_gene_id"),"ensembl_gene_id", res_un$ensembl, mart)
res_un<-merge(res_un,symbol,by.x="ensembl",by.y="ensembl_gene_id")
#add normalized counts to results table
counts2_un<-DESeq2::counts(cds_un,normalized=TRUE)
counts2_un<-cbind(counts2_un,ensembl=rownames(counts2_un))
res_un<-merge(counts2_un,res_un,by="ensembl")
#remove genes that have NA for log2FC and adjP
res_un_f<-res_un[!is.na(res_un$log2FoldChange)&!is.na(res_un$padj),]
#sort differential expressed genes by log2FC in descending order
res_un_f<-arrange(res_un_f,desc(abs(log2FoldChange)))
#remove duplicated entries (By sorting data first, we remove the once with lower fold change)
res_un_f<-res_un_f[!duplicated(res_un_f$mgi_symbol),]

# write results to table. These data is published as processed data along with the raw counts on the NCBI Gene Omnibus with the accession number GSE.

write.table(res_un_f, file="./2018-3_liver_E47Kovswt_untreated_unfiltered.txt",sep="\t",row.names = FALSE)

#6. Vulcano plot to display differential expression

#add color flag for genes significantly (P<0.05) changing more the 1.3x
res_un_f$color_flag <- factor(ifelse(res_un_f$log2FoldChange>0.38&res_un_f$pvalue<0.05, 1, ifelse(res_un_f$log2FoldChange<(-0.38)&res_un_f$pvalue<0.05,2,0)))
#plot
ggplot(data=res_un_f, aes(x=log2FoldChange, y=-log10(pvalue),col=color_flag)) +
  geom_point(alpha=0.7, size=5) +
  xlim(c(-5, 5)) + ylim(c(0, 15)) +
  xlab("log2 fold change") + ylab("-log10 p-value")+
  geom_text_repel(aes(label=ifelse(abs(log2FoldChange)>0.58&padj<0.05, mgi_symbol,"")), size=3)+
  theme_bw()+
  scale_color_manual(values=c("grey","red","blue"))+
  theme(legend.position = "none",panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
