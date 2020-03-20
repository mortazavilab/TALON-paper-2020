setwd("~/Dropbox/Gaby/TALON-2020/Figure5/")

abundance_Brain <- read.delim("Brain_talon_abundance_filtered.tsv",sep="\t",header = T)

Brain_PB <- data.frame(abundance_Brain$annot_transcript_name,abundance_Brain$PacBio_Cortex_Rep1,abundance_Brain$PacBio_Cortex_Rep1,
           abundance_Brain$PacBio_Hippocampus_Rep1, abundance_Brain$PacBio_Hippocampus_Rep2)

annot_names <- data.frame(abundance_Brain$annot_transcript_name,abundance_Brain$annot_gene_name) 
names(annot_names)[1]="transcript"
names(annot_names)[2]="gene_name"

novelty_transcript_names <- data.frame(abundance_Brain$annot_transcript_name,abundance_Brain$annot_gene_name,abundance_Brain$transcript_novelty) 
names(novelty_transcript_names)[1]="annot_transcript_name"
names(novelty_transcript_names)[2]="gene_name"
names(novelty_transcript_names)[3]="Novelty"

#Rename columns
names(Brain_PB)[1]="transcript_id"
names(Brain_PB)[2]="Cortex_rep1"
names(Brain_PB)[3]="Cortex_rep2"
names(Brain_PB)[4]="Hippocampus_rep1"
names(Brain_PB)[5]="Hippocampus_rep2"

#get rid of duplicates, by aggregation of data based on gene_id column
#install.packages("dplyr")
#library(dplyr)
#library(tibble)
#Brain_PB_agg <- Brain_PB %>% group_by(gene_id) %>% summarize_all(sum)


#Column1 to rownames
#install.packages("tidyverse")

rownames(Brain_PB) = Brain_PB$transcript_id

#get rid of Column1
Brain_PB$transcript_id=NULL

#Round numbers for EdgeR
Brain_PB <- round(Brain_PB)

#Save matrix
write.table(Brain_PB, "CountTranscripts_matrix_PacBio_4EdgeR_pvalue01.txt")

#EdgeR
library(edgeR)

#in case I am working from the table i save rpeviously
group <- factor(c("Cortex","Cortex","Hippocampus","Hippocampus"))
y <- DGEList(counts=Brain_PB,group=group)
y <- calcNormFactors(y)
design <- model.matrix(~group)
y <- estimateDisp(y,design)

#genes are dropped if they can’t possibly be expressed in all the samples for any of the conditions
keep <- rowSums(cpm(y)>1) >= 2
y <- y[keep, , keep.lib.sizes=FALSE]

#Testing for DE Genes
et <- exactTest(y)
topTags(et,n = Inf)

#
fit <- glmFit(y, design)
lrt <- glmLRT(fit)
topTags(lrt)
o <- order(lrt$table$PValue)
cpm(y)[o[1:10],]

#DE genes summary
summary(decideTests(lrt))
plotMD(lrt)
abline(h=c(-1, 1), col="blue")

#Extract et table for VolcanoPlot

Brain_PB_out <- (topTags(et,n=Inf,adjust.method = "bonferroni"))
Brain_PB_et <- Brain_PB_out$table
Brain_PB_et$annot_transcript_name <- rownames(Brain_PB_et)

View(Brain_PB_et)


##Add noveltyType 
Brain_PB_et_novelty <- merge(Brain_PB_et,novelty_transcript_names,by="annot_transcript_name",all=F)
View(Brain_PB_et_novelty)
##save table without the DE column to submit the paper
#first change the FWER for adj.pValue
names(Brain_PB_et_novelty)[5]="adj.pValue"
write.table(Brain_PB_et_novelty,"Brain_PB_et_novelty.txt",sep="\t")


#label DE genes
Brain_PB_et_novelty$DE <- ""
Brain_PB_et_novelty$DE[ Brain_PB_et_novelty$logFC > 1 &
                 Brain_PB_et_novelty$adj.pValue < 0.01 ] = "Up"

Brain_PB_et_novelty$DE[ Brain_PB_et_novelty$logFC < -1 &
                  Brain_PB_et_novelty$adj.pValue < 0.01 ] = "Down"




View(Brain_PB_et_novelty)
##add a number to sort the dataframe
Brain_PB_et_novelty$order=0
Brain_PB_et_novelty$order[Brain_PB_et_novelty$Novelty == "Known"]= 1
Brain_PB_et_novelty$order[Brain_PB_et_novelty$Novelty == "ISM"] = 5
Brain_PB_et_novelty$order[Brain_PB_et_novelty$Novelty == "NIC"] = 6
Brain_PB_et_novelty$order[Brain_PB_et_novelty$Novelty == "NNC"] = 7
Brain_PB_et_novelty$order[Brain_PB_et_novelty$Novelty == "Antisense"] = 4
Brain_PB_et_novelty$order[Brain_PB_et_novelty$Novelty == "Intergenic"] = 3
Brain_PB_et_novelty$order[Brain_PB_et_novelty$Novelty == "Genomic"] = 2

#Brain_PB_DElist <- data.frame(Brain_PB_et$transcript[Brain_PB_et$DE != ""])
#Brain_PB_DElist$DE <-Brain_PB_et$DE[Brain_PB_et$DE != ""]
#names(Brain_PB_DElist)[1]="Transcript"

#merge with annotnames so I don't have to split the transcript column
#Brain_PB_et_transgenes <- merge(Brain_PB_et,annot_names,by="transcript",all = F)

#select for genes that have a transcript with higher FC in cortex
Cortex_FC_genes <- unique(subset(Brain_PB_et_novelty, logFC <-1 & adj.pValue < 0.01 )$gene_name)
Hippo_FC_genes <- unique(subset(Brain_PB_et_novelty, logFC > 1 & adj.pValue < 0.01 )$gene_name)

logFC <- intersect(Cortex_FC_genes,Hippo_FC_genes)
View(Brain_PB_et_novelty)

library(ggrepel)

#pretty volcano plot
#colors
Brain_PB_et_novelty$color <- ""
Brain_PB_et_novelty$color[Brain_PB_et_novelty$Novelty == "Known" & Brain_PB_et_novelty$adj.pValue < 0.01] = "#009E73"
Brain_PB_et_novelty$color[Brain_PB_et_novelty$Novelty == "ISM" & Brain_PB_et_novelty$adj.pValue < 0.01] = "#0072B2"
Brain_PB_et_novelty$color[Brain_PB_et_novelty$Novelty == "NIC" & Brain_PB_et_novelty$adj.pValue < 0.01] = "#D55E00"
Brain_PB_et_novelty$color[Brain_PB_et_novelty$Novelty == "NNC" & Brain_PB_et_novelty$adj.pValue < 0.01] = "#E69F00"
Brain_PB_et_novelty$color[Brain_PB_et_novelty$Novelty == "Antisense" & Brain_PB_et_novelty$adj.pValue < 0.01] = "#000000"
Brain_PB_et_novelty$color[Brain_PB_et_novelty$Novelty == "Intergenic" & Brain_PB_et_novelty$adj.pValue < 0.01] = "#CC79A7"
Brain_PB_et_novelty$color[Brain_PB_et_novelty$Novelty == "Genomic" & Brain_PB_et_novelty$adj.pValue < 0.01] = "#F0E442"
Brain_PB_et_novelty$color[Brain_PB_et_novelty$color == ""] = "grey"
table(Brain_PB_et_novelty$color,Brain_PB_et_novelty$Novelty)

#add the genes with the highest fold change to the NovGenes(because we use this for label)
Brain_PB_et_novelty$Label=""
Brain_PB_et_novelty$Label[Brain_PB_et_novelty$annot_transcript_name == "Pnisr-201"]="Pnisr-201"
Brain_PB_et_novelty$Label[Brain_PB_et_novelty$annot_transcript_name == "Pnisr-203"]="Pnisr-203"
Brain_PB_et_novelty$Label[Brain_PB_et_novelty$annot_transcript_name == "ENCODEMT000648241"]="Meg3"
Brain_PB_et_novelty$Label[Brain_PB_et_novelty$annot_transcript_name == "ENCODEMT000638681"]="Srsf5"
Brain_PB_et_novelty$Label[Brain_PB_et_novelty$annot_transcript_name == "ENCODEMT000361240"]="Amigo2"
Brain_PB_et_novelty$Label[Brain_PB_et_novelty$annot_transcript_name == "ENCODEMT000184694"]="2010300C02Rik"
Brain_PB_et_novelty$Label[Brain_PB_et_novelty$annot_transcript_name == "ENCODEMT000648240"]="Meg3"
Brain_PB_et_novelty$Label[Brain_PB_et_novelty$annot_transcript_name == "ENCODEMT000352747"]="Nell2"


ENCODEMT000648241 - Meg3
ENCODEMT000638681 - Srf5
ENCODEMT000361240 - Amigo2
ENCODEMT000633920 - Map4k5
ENCODEMT000648240 - Meg3 
ENCODEMT000352747 - Nell2

View(Brain_PB_et)

#set transparency
Brain_PB_et_novelty$alpha <- 0
Brain_PB_et_novelty$alpha[Brain_PB_et_novelty$Novelty == "Known"] = 0.5
Brain_PB_et_novelty$alpha[Brain_PB_et_novelty$Novelty != "Known"] = 1


#sort the data frame by the order column
Brain_PB_et_novelty_order <- Brain_PB_et_novelty[order(Brain_PB_et_novelty$order),] 
View(Brain_PB_et_novelty_order)


#original plot
PB <- ggplot(Brain_PB_et_novelty_order,aes(logFC, -log10(adj.pValue))) + geom_point(color=Brain_PB_et_novelty_order$color,
                                                                       alpha=Brain_PB_et_novelty_order$alpha) 
PB <- PB + theme_bw() + geom_label_repel(aes(label=Label),size=7)
PB <- PB + theme(axis.text=element_text(size=20),
            axis.title=element_text(size=22)) + 
           xlab("log2(Fold Change Hippocampus/Cortex)") +
           ylab("-log10(adjusted pValue)")

PB



colored_points <- subset(Brain_PB_et_novelty, color != "grey" )
table(colored_points$Novelty)
unique(colored_points$gene_name)
write.table(colored_points,"colored_dots_table", sep="\t")

view(colored_points)

#make a bar plot for the isoforms in the UCSC genBrow
#subset the Psnir reads
Psnir_abundance <- subset(abundance_Brain, annot_gene_name=="Pnisr" &  transcript_novelty != "ISM" )
Psnir_abundance2 <- data.frame(Psnir_abundance$annot_transcript_name, Psnir_abundance$PacBio_Cortex_Rep1,Psnir_abundance$PacBio_Cortex_Rep2,
                               Psnir_abundance$PacBio_Hippocampus_Rep1,Psnir_abundance$PacBio_Hippocampus_Rep2)
#rownames(Psnir_abundance2)=Psnir_abundance2$Psnir_abundance.annot_transcript_name

#Psnir_abundance2$Psnir_abundance.annot_transcript_name=NULL

Psnir_abundance3 <- data.frame(t(Psnir_abundance2))
names(Psnir_abundance3)<- c("Pnisr-201","Pnisr-211","Pnisr-203","Pnisr-206","ENCODEMT000409037")
Psnir_abundance3 <- Psnir_abundance3[-1, ]
Psnir_abundance3$Tissue<- c("Cortex","Cortex","Hippocampus", "Hippocampus")
View(Psnir_abundance3)

library(reshape2)
library(ggplot2)
Psnir_abundance2$Transcript <- factor(c("Pnisr-201","Pnisr-211","Pnisr-203","Pnisr-206","ENCODEMT000409037"), 
                  levels=c("Pnisr-201","Pnisr-211","Pnisr-203","Pnisr-206","ENCODEMT000409037"))


mdat <- melt(Psnir_abundance2, id.vars="Transcript")
head(mdat)
mdat$Tissue <- ""
mdat$Tissue[mdat$variable == "Psnir_abundance.PacBio_Cortex_Rep1"] = "Cortex"
mdat$Tissue[mdat$variable == "Psnir_abundance.PacBio_Cortex_Rep2"] = "Cortex"
mdat$Tissue[mdat$variable == "Psnir_abundance.PacBio_Hippocampus_Rep1"] = "Hippocampus"
mdat$Tissue[mdat$variable == "Psnir_abundance.PacBio_Hippocampus_Rep2"] = "Hippocampus"

mdat$order <- 0
mdat$order[mdat$Transcript == "Pnisr-201"] = 1
mdat$order[mdat$Transcript == "Pnisr-203"] = 2
mdat$order[mdat$Transcript == "Pnisr-206"] = 3
mdat$order[mdat$Transcript == "Pnisr-211"] = 4
mdat$order[mdat$Transcript == "ENCODEMT000409037	"] = 5

mdat_order <- mdat[order(mdat$order),] 

View(mdat_order)
p <- ggplot(mdat_order, aes(Transcript, value, fill=Tissue)) + 
  geom_bar(stat="identity", position="dodge") + theme_classic() 

p <- p +scale_fill_manual(values=c("lightpink1","paleturquoise2"))

p <- p + ylab("average TPM per Tissue") + theme(axis.text=element_text(size=10),
                                              axis.title=element_text(size=22))

p
