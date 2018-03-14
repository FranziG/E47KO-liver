library("DESeq2")
library("AnnotationDbi")
library("org.Mm.eg.db")
library("biomaRt")
library("gplots")
library("ggplot2")
library("dplyr")

# data preprocessing: make count table. The raw count data is diposited at GEO under the accession GSE.

# 1. read STAR count files: 1st column: gene ID (UCSC), 2nd column: counts for both strands, 3rd: 5´-->3´counts; 4th: 3`-->5`counts

#Samples for the cort anaylsis were taken from mice in different labs resulting in a batch effect that was included as cofounder in the analysis.
#old samples
GAR0886<-read.table("./GAR0886_liver_cort_E47wt_No2_mm10_STARReadsPerGene.out.tab", header=FALSE, sep="\t",stringsAsFactors = FALSE)
GAR0886<-GAR0886[c(5:(nrow(GAR0886))),]
GAR0886<-arrange(GAR0886,V1)
GAR0887<-read.table("./GAR0887_liver_cort_E47KO_No3_mm10_STARReadsPerGene.out.tab", header=FALSE, sep="\t",stringsAsFactors = FALSE)
GAR0887<-GAR0887[c(5:(nrow(GAR0887))),]
GAR0887<-arrange(GAR0887,V1)
GAR0888<-read.table("./GAR0888_liver_cort_E47KO_No15_mm10_STARReadsPerGene.out.tab", header=FALSE, sep="\t",stringsAsFactors = FALSE)
GAR0888<-GAR0888[c(5:(nrow(GAR0888))),]
GAR0888<-arrange(GAR0888,V1)
#newer samples
GAR0890<-read.table("./GAR0890_liver_cort_E47wt_No195_mm10_STARReadsPerGene.out.tab", header=FALSE, sep="\t",stringsAsFactors = FALSE)
GAR0890<-GAR0890[c(5:(nrow(GAR0890))),]
GAR0890<-arrange(GAR0890,V1)
GAR0891<-read.table("./GAR0891_liver_cort_E47KO_No197_mm10_STARReadsPerGene.out.tab", header=FALSE, sep="\t",stringsAsFactors = FALSE)
GAR0891<-GAR0891[c(5:(nrow(GAR0891))),]
GAR0891<-arrange(GAR0891,V1)
GAR0892<-read.table("./GAR0892_liver_cort_E47KO_No199_mm10_STARReadsPerGene.out.tab", header=FALSE, sep="\t",stringsAsFactors = FALSE)
GAR0892<-GAR0892[c(5:(nrow(GAR0892))),]
GAR0892<-arrange(GAR0892,V1)
GAR0893<-read.table("./STAR/GAR0893_liver_cort_E47wt_No202_mm10_STARReadsPerGene.out.tab", header=FALSE, sep="\t",stringsAsFactors = FALSE)
GAR0893<-GAR0893[c(5:(nrow(GAR0893))),]
GAR0893<-arrange(GAR0893,V1)
GAR0894<-read.table("./GAR0894_liver_cort_E47wt_No203_mm10_STARReadsPerGene.out.tab", header=FALSE, sep="\t",stringsAsFactors = FALSE)
GAR0894<-GAR0894[c(5:(nrow(GAR0894))),]
GAR0894<-arrange(GAR0894,V1)

#2. generate one table for raw countdata; We used the counts ignoring strandedness as not all RNA-Seq samples of this study were stranded.

counts_cortOnly<-cbind(GAR0886[,c(1,2)],GAR0890$V2,GAR0893$V2,GAR0894$V2,GAR0891$V2,GAR0892$V2,GAR0887$V2,GAR0888$V2)
rownames(counts_cortOnly)<-counts_cortOnly[,1]
counts_cortOnly<-counts_cortOnly[,-1]
colnames(counts_cortOnly)<-c("GAR0886","GAR0890","GAR0893","GAR0894","GAR0891","GAR0892","GAR0887","GAR0888")
counts_cortOnly[is.na(counts_cortOnly)]<-0
#counts_all<-counts_all[rowMeans(counts_all)>5,]
#filter lowest 25% counts (rowmean below )
#quantile(counts_all$E47wt_Dex_GAR0886,.25)
#886: 57; 887:51; 888: 44; 890:63 891:50; 892:42 ; 893:51 ; 894:44 => lowest is 42
counts_cortOnly<-counts_cortOnly[rowMeans(counts_cortOnly)>42,]

#3. prepare metadata

genotype_cortOnly<-c("wt","wt","wt","wt","KO","KO","KO","KO")
samples_cortOnly<-colnames(counts_cortOnly)
#We introduce an additional cofounding variable here, which we called batch effect as we observed differences according to the mouse house. 1-old samples; 2-new samples
batch_cortOnly<-c(1,2,2,2,2,2,1,1)
pdata_cortOnly<-as.data.frame(cbind(samples_cortOnly,genotype_cortOnly,batch_cortOnly))

#4. Differential expression analysis using DESeq2

cds_cortOnly<-DESeqDataSetFromMatrix(countData=counts_cortOnly,colData=pdata_cortOnly,design=~batch_cortOnly+genotype_cortOnly)
#set wildtype as basal level
cds_cortOnly$genotype_cortOnly <- relevel(cds_cortOnly$genotype_cortOnly, "wt")
cds_cortOnly<-DESeq(cds_cortOnly)

# get differential expression wt versus E47KO
res_cortOnly_indOnBatch<-data.frame(results(cds_cortOnly,name = "genotype_cortOnly_KO_vs_wt"),stringsAsFactors = FALSE)

#5. annotated gene IDs

#load annotation DB
mart <- useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")
#annotate results object
res_cortOnly_indOnBatch$ensembl<-rownames(res_cortOnly_indOnBatch)
symbol<-getBM(c("mgi_symbol","ensembl_gene_id","entrezgene"),"ensembl_gene_id",res_cortOnly_indOnBatch$ensembl, mart)
res_cortOnly_indOnBatch<-merge(res_cortOnly_indOnBatch,symbol,by.x="ensembl",by.y="ensembl_gene_id")
#add normalized counts to results table
counts_cortOnly2<-data.frame(DESeq2::counts(cds_cortOnly,normalized=TRUE),stringsAsFactors = FALSE)
counts_cortOnly2<-cbind(counts_cortOnly2,ensembl=rownames(counts_cortOnly2))
res_cortOnly_indOnBatch<-merge(counts_cortOnly2,res_cortOnly_indOnBatch,by="ensembl")
#remove genes that have NA for log2FC and adjP
res_cortOnly_indOnBatch_f<-res_cortOnly_indOnBatch[!is.na(res_cortOnly_indOnBatch$log2FoldChange)&!is.na(res_cortOnly_indOnBatch$padj),]
#sort differential expressed genes by log2FC in descending order
res_cortOnly_indOnBatch_f<-arrange(res_cortOnly_indOnBatch_f,desc(abs(log2FoldChange)))
#remove duplicated entries (By sorting data first, we remove the once with lower fold change)
res_cortOnly_indOnBatch_f<-res_cortOnly_indOnBatch_f[!duplicated(res_cortOnly_indOnBatch_f$mgi_symbol),]

#write results to table. These data is published as processed data along with the raw counts on the NCBI Gene Omnibus with the accession number GSE.

write.table(res_cortOnly_indOnBatch_f, file="./2018-3_liver_E47Kovswt_longTermCort_unfiltered.txt",sep="\t",row.names = FALSE)

#6. Heatmap to display differential expression

counts_cortOnly2<-DESeq2::counts(cds_cortOnly,normalized=TRUE)
colfunc<-colorRampPalette(colors=c("DodgerBlue2","black","red"))
col_breaks<-c(seq(min(log2(counts_cortOnly2+1)),quantile(log2(counts_cortOnly2+1),.45),length=50),seq(quantile(log2(counts_cortOnly2+1),.46),quantile(log2(counts_cortOnly2+1),.55),length=10),seq(quantile(log2(counts_cortOnly2+1),.56),max(log2(counts_cortOnly2+1)),length=50))
#Only use differential expressed genes with a baseMean expression above 200 counts and a significant (adjP<0.01) fold-change of 1.5.
res_cortOnly_indOnBatch_f<-res_cortOnly_indOnBatch_f[abs(res_cortOnly_indOnBatch_f$log2FoldChange)>0.58&!duplicated(res_cortOnly_indOnBatch_f$mgi_symbol)&res_cortOnly_indOnBatch_f$padj<0.05&res_cortOnly_indOnBatch_f$baseMean>200,]

heatmap.2(log2(counts_cortOnly2[rownames(counts_cortOnly2)%in%res_cortOnly_indOnBatch_f$ensembl,]+1),na.rm=TRUE,scale="row",distfun = function(x) as.dist(1-cor(t(x))),trace="none",dendrogram = "none",Colv = FALSE,col=colfunc(length(col_breaks)-1), main="DE_genes\n in liver lacking E47KO after cort",cexRow=0.6,cexCol=0.6,labRow = FALSE,labCol=FALSE,margins = c(10,10))
