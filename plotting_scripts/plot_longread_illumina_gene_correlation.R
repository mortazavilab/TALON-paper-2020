main <- function() {

    load_packages()
    opt <- parse_options()

    # Get colors
    if (opt$color_scheme == "red") {
        color <- "red2"
    } else if (opt$color_scheme == "blue") {
        color <- "navy"
    } else if (opt$color_scheme == "green") {
        color <- "yellowgreen"
    }

    # Average the gene TPMs of the two long read replicates
    longread_data <- process_longread_data(opt)
    lr_abundance <- longread_data$abundance

    # Filter Illumina genes
    shortread_data <- process_shortread_data(opt)   
    sr_abundance <- shortread_data$abundance

    print(nrow(lr_abundance))
    print(nrow(sr_abundance))

    # Merge 
    merged_data <- merge(sr_abundance, lr_abundance, by = "gene_ID",
                         all.x = T, all.y = T)
    merged_data[is.na(merged_data)] <- 0
    colnames(merged_data) <- c("gene_ID", "sr_expression", "lr_expression")

    # Setup for output
    celltype = opt$celltype
    joined_names <- paste(opt$outdir, "/", opt$longread_type, "-", 
                          opt$shortread_type, sep = "")

    if (opt$byCount == F) {
        outprefix.1 <- paste(joined_names, celltype, "gene", "TPM", "allGenes", sep="_")
        outprefix.2 <- paste(joined_names, celltype, "gene", "TPM", "mutualGenes", sep="_")
    } else {
        outprefix.1 <- paste(joined_names, celltype, "gene", "count", "allGenes", sep="_")
        outprefix.2 <- paste(joined_names, celltype, "gene", "count", "mutualGenes", sep="_")
    } 
    plot_expression(merged_data, opt, color, outprefix.1)

    mutual <- subset(merged_data, sr_expression > 0 & lr_expression > 0)
    print(nrow(mutual))
    plot_expression(mutual, opt, color, outprefix.2)   
}

process_shortread_data <- function(opt) {

    # Get genes and transcripts expressed in the Illumina data from the Kallisto
    # abundance file
    illumina_1 <- filter_kallisto_illumina_genes(opt$illumina_kallisto_1)
    illumina_2 <- filter_kallisto_illumina_genes(opt$illumina_kallisto_2)
    colnames(illumina_1) <- c("gene_ID", "illumina_TPM_1")
    colnames(illumina_2) <- c("gene_ID", "illumina_TPM_2")
    illumina_gene_table <- merge(illumina_1, illumina_2, by = "gene_ID",
                                 all.x = T, all.y = T)

    # Average the TPMs
    # Genes must be detected in both reps
    illumina_gene_table[is.na(illumina_gene_table)] <- 0
    illumina_gene_table$tpm <- (illumina_gene_table$illumina_TPM_1 +
                                illumina_gene_table$illumina_TPM_2) / 2

    illumina_gene_table <- illumina_gene_table[,c("gene_ID", "tpm")]

    if (opt$byCount == T) {
        illumina_gene_table$tpm <- NULL
    } else {
        illumina_gene_table$count <- NULL
    }

    data <- list("abundance" = illumina_gene_table, "dname" = opt$shortread_type)
    return(data)
}

filter_kallisto_illumina_genes <- function(kallisto_file) {
    # This function takes a Kallisto abundance file and filters the genes
    # based on criteria designed to make the gene set comparable to what can
    # be detected using PacBio

    gencode.quantitation <- as.data.frame(read_delim(kallisto_file, "\t", escape_double = FALSE,
                                  col_names = TRUE, trim_ws = TRUE, na = "NA"))

    # Split GENCODE transcript multi-id by '|'
    extraCols =str_split_fixed(gencode.quantitation$target_id, "\\|",9)[,c(5,6,8,1,2)]
    colnames(extraCols) = c("transcript", "gene", "class", "t_ID", "g_ID")
    gencode.quantitation = cbind(extraCols, gencode.quantitation)

    # Aggregate on genes
    gene_gencode <- gencode.quantitation %>%
                                 dplyr::group_by(g_ID) %>%
                                 dplyr::summarize(tpm = sum(tpm))
    gene_gencode <- subset(gene_gencode, tpm > 0)
    colnames(gene_gencode) <- c("gene_ID", "tpm")

    return(gene_gencode)
}

process_longread_data <- function(opt) {
    # Read abundance table
    abundance_table <- as.data.frame(read_delim(opt$infile, delim = "\t",
                                  col_names = TRUE, trim_ws = TRUE, na = "NA"))

    # Remove genomic
    abundance_table <- subset(abundance_table, transcript_novelty != "Genomic")

    # Get dataset name and check that it exists in the file
    d1 <- opt$r1
    d2 <- opt$r2
    longread_type <- opt$longread_type
    if (d1 %in% colnames(abundance_table) == F |
        d2 %in% colnames(abundance_table) == F) {
        print("One of the provided dataset names is not in the abundance file provided. Exiting...")
        quit()
    }

    abundance_table <- subset(abundance_table, gene_novelty == "Known")

    # Restrict to genes that were observed in these datasets
    abundance_table <- abundance_table[abundance_table[,d1] + abundance_table[,d2] > 0, ]

    d1_abundance <- abundance_table[,c("annot_gene_id", d1)]
    d2_abundance <- abundance_table[,c("annot_gene_id", d2)]
    colnames(d1_abundance) <- c("gene_ID", "count")
    colnames(d2_abundance) <- c("gene_ID", "count")

    if (opt$byCount == F) {
        # Compute gene TPMs for each gene in each dataset
        d1_abundance$CPM <- (d1_abundance$count)*1000000/sum(d1_abundance$count)
        d2_abundance$CPM <- (d2_abundance$count)*1000000/sum(d2_abundance$count)
 
        d1_abundance <- aggregate(d1_abundance$CPM, 
                                      by=list(d1_abundance[, "gene_ID"]), FUN=sum)
        colnames(d1_abundance) <- c( "gene_ID", "CPM")
        d2_abundance <- aggregate(d2_abundance$CPM, 
                                      by=list(d2_abundance[, "gene_ID"]), FUN=sum)
        colnames(d2_abundance) <- c("gene_ID", "CPM")
    }
     
    merged_abundances <- merge(d1_abundance, d2_abundance, by = "gene_ID", all.x = T, all.y = T)
    colnames(merged_abundances) <- c("gene_ID", "data1", "data2")

    # Remove genes only detected in only one rep, then average the TPMS
    merged_abundances[is.na(merged_abundances)] <- 0
    merged_abundances$expression <- (merged_abundances$data1 +
                                     merged_abundances$data2) / 2
    merged_abundances <- merged_abundances[, c("gene_ID", "expression")]
    merged_abundances <- subset(merged_abundances, expression > 0)

    data <- list("abundance" = merged_abundances, "dname" = longread_type)
 
    return(data)

}

plot_expression <- function(merged_abundances, options, color, outprefix) {

    # Take log2(TPM + 0.1)
    merged_abundances$lr_log_expression = log(merged_abundances$lr_expression + 0.1, base=10)
    merged_abundances$sr_log_expression = log(merged_abundances$sr_expression + 0.1, base=10)

    # Compute correlations
    pearsonCorr = cor.test(~lr_expression + sr_expression,
                           data=merged_abundances, method = "pearson",
                           continuity = FALSE, conf.level = 0.95)$estimate
    spearmanCorr = cor.test(~lr_expression + sr_expression,
                             data=merged_abundances, method = "spearman",
                             continuity = FALSE, conf.level = 0.95)$estimate

    # Least-Square Regression Line
    mod<-lm(sr_expression~lr_expression, data=merged_abundances)

    # Output names and labels
    celltype = options$celltype
    fname <- paste(outprefix, "correlationPlot.png", sep="_")
    corr_fname <- paste(outprefix, "correlations.txt", sep="_")
    if (options$byCount == F) {
        xlabel <- paste("CPM+0.1 in ", celltype, " ", options$longread_type, sep="")
        ylabel <- paste("TPM+0.1 in ", celltype, " ", options$shortread_type, sep="")
    } else {
        xlabel <- paste("Count+0.1 in ", celltype, " ", options$longread_type, sep="")
        ylabel <- paste("Count+0.1 in ", celltype, " ", options$shortread_type, sep="")
    }

    corr_label <- paste("Pearson r: ",
                            round(pearsonCorr, 2), "\nSpearman rho: ",
                            round(spearmanCorr, 2), "\nLSR slope: ",
                            round(mod$coefficients[2], 2), sep="")
    if (options$lsr == T) {
        plot_label <- paste("Pearson r: ",
                            round(pearsonCorr, 2), "\nSpearman rho: ",
                            round(spearmanCorr, 2), "\nLSR slope: ",
                            round(mod$coefficients[2], 2), sep="")
    } else {
        plot_label <- paste("Pearson r: ",
                    round(pearsonCorr, 2), "\nSpearman rho: ",
                    round(spearmanCorr, 2), "\nLSR slope: ", sep="")
    }

    # write correlation numbers to outfile
    write(corr_label, corr_fname)

    png(filename = fname,
        width = 2700, height = 2500, units = "px",
        bg = "white",  res = 300)

    # Main scatterplot
    scatterplot = ggplot(merged_abundances, aes(x = lr_expression + 0.1, 
                                                y = sr_expression + 0.1)) +
                         geom_jitter(alpha = 0.5, color = color) + theme_bw() +
                         xlab(xlabel)  + ylab(ylabel) +
                         theme(text= element_text(size=24)) +
                         theme(axis.text.x = element_text(color = "black", size=24),
                               axis.text.y = element_text(color = "black", size=24)) +
                         scale_x_continuous(trans=log10_trans(), limits=c(0.1,1000000))+
                         scale_y_continuous(trans=log10_trans(), limits=c(0.1,1000000))+

    # add regression line
    if (options$regression_line) {
        scatterplot <- scatterplot+
              geom_abline(slope = mod$coefficients[2], intercept = mod$coefficients[1],
                   color = "gray", lwd=1, lty=2)
    }

    if (options$corr_labs) {
        scatterplot <- scatterplot+
              annotate("text", x = 5, y = 14, label = plot_label,
                   color="black", size = 10)
    }

    print(plot_label)
    print(scatterplot)
    dev.off()

}

compute_max_density <- function(data, dataset) {
    # Compute the max density y value
    data <- data[,dataset]
    return(max(density(data)$y))
}


load_packages <- function() {
    suppressPackageStartupMessages(library("tidyverse"))
    suppressPackageStartupMessages(library("ggplot2"))
    suppressPackageStartupMessages(library("optparse"))
    suppressPackageStartupMessages(library("gridExtra"))
    suppressPackageStartupMessages(library("cowplot"))
    suppressPackageStartupMessages(library("scales"))

    return
}

parse_options <- function() {

    option_list <- list(
        make_option(c("--f"), action = "store", dest = "infile",
                    default = NULL, help = "UNFILTERED TALON abundance output file"),
        make_option(c("--ik1"), action = "store", dest = "illumina_kallisto_1",
                    default = NULL, help = "Rep1 Illumina Kallisto file."),
        make_option(c("--ik2"), action = "store", dest = "illumina_kallisto_2",
                    default = NULL, help = "Rep2 Illumina Kallisto file."),
        make_option(c("--color"), action = "store", dest = "color_scheme",
                    default = NULL, help = "blue, red, or green"),
        make_option("--r1", action = "store", dest = "r1",
                    default = NULL, help = "Dataset name of long read Rep 1"),
        make_option("--lrtype", action="store", dest="longread_type",
                    help="Data type of long reads (i.e. PacBio)"),
        make_option("--r2", action = "store", dest = "r2",
                    default = NULL, help = "Dataset name of long read Rep 2"),
        make_option("--srtype", action="store", dest="shortread_type",
                    help="Datatype of short read dataset ie 'Illumina'"),
        make_option(c("--lsr"), action="store_true", dest="lsr",
              help="Include this option if you want the LSR label on", default = F),
        make_option("--celltype", action = "store", dest = "celltype",
                    default = NULL, help = "Celltype to use in plot labels"),
        make_option(c("-o","--outdir"), action = "store", dest = "outdir",
                    default = NULL, help = "Output directory for plots and outfiles"),
        make_option("--correlations", action = "store_true", dest = "corr_labs",
              help="Add correlation labels to plot", default = F),
        make_option("--regression_line", action = "store_true", dest = "regression_line",
              help="Add regression line to plot", default = F),
        make_option("--byCount", action = "store_true", dest = "byCount",
              help="If false (default), correlation is computed on gene TPMs. If true, correlation is computed on read counts instead.", default = F)
        )

    opt <- parse_args(OptionParser(option_list=option_list))
    opt$r1 <- as.character(opt$r1)
    opt$r2 <- as.character(opt$r2)
    return(opt)
}

main()
