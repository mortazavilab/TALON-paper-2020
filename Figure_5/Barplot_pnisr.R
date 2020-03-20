abundance_Brain <- read.delim("Brain_talon_abundance_filtered.tsv",sep="\t",header = T)

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
