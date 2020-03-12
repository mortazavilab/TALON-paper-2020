main <-function() {

    set.seed(100)
    load_packages()
    opt <- parse_options()

    # Get colors
    if (opt$color_scheme == "red") {
        fill_color <- "red2"
    } else if (opt$color_scheme == "blue") {
        fill_color <- "navy"
    } else if (opt$color_scheme == "green") {
        fill_color <- "yellowgreen"
    }

    # Get the names of the first and second dataset that we will be working with
    data_names <- str_split(opt$datasets, ",")[[1]]
    dataset1 <- data_names[1]
    dataset2 <- data_names[2]

    # datatype label (ONT or PacBio)
    dtype <- opt$dtype

    # Get genes expressed in the Illumina data from the Kallisto
    # abundance files
    illumina_1 <- filter_kallisto_illumina_genes(opt$illumina_kallisto_1)
    illumina_2 <- filter_kallisto_illumina_genes(opt$illumina_kallisto_2)
    colnames(illumina_1) <- c("annot_gene_id", "illumina_counts_1")
    colnames(illumina_2) <- c("annot_gene_id", "illumina_counts_2")
    illumina_gene_table <- merge(illumina_1, illumina_2, by = "annot_gene_id",
                                 all.x = T, all.y = T)
    print(paste0("Illumina genes: ", nrow(illumina_gene_table)))
    illumina_gene_table[is.na(illumina_gene_table)] <- 0
    illumina_gene_table$illumina_counts_1 <- round(illumina_gene_table$illumina_counts_1)
    illumina_gene_table$illumina_counts_2 <- round(illumina_gene_table$illumina_counts_2)

    # Read PacBio abundance file
    pb_abundance_orig <- as.data.frame(read_delim(opt$infile, "\t", escape_double = FALSE,
                                  col_names = TRUE, trim_ws = TRUE, na = "NA"))

    # Keep known genes only
    pb_abundance <- subset(pb_abundance_orig, gene_novelty == "Known")

    # Remove genomic transcripts
    pb_abundance <- subset(pb_abundance, transcript_novelty != "Genomic")

    # Cut out unnecessary cols
    pb_abundance <- pb_abundance[, c("annot_gene_id", dataset1, dataset2)]

    # Aggregate PacBio by gene ID to get gene-wise counts
    pb_gene_abundance <- ddply(pb_abundance, c("annot_gene_id"), 
                               function(x) colSums(x[c(dataset1, dataset2)]))

    # Merge PacBio with Illumina on annot_gene_name
    merged_illumina_pacbio <- merge(illumina_gene_table, pb_gene_abundance, by = "annot_gene_id",
                                    all.x = T, all.y = T)

    merged_illumina_pacbio[is.na(merged_illumina_pacbio)] <- 0
    merged_illumina_pacbio <- merged_illumina_pacbio[, c("annot_gene_id", 
                                                         "illumina_counts_1", 
                                                         "illumina_counts_2",  
                                                         dataset1, dataset2)]

    # Now, format table for edgeR by setting the rownames to the gene names
    rownames(merged_illumina_pacbio) <- merged_illumina_pacbio$annot_gene_id
    merged_illumina_pacbio$annot_gene_id <- NULL

    # Remove rows that only contain zeros
    merged_illumina_pacbio <- merged_illumina_pacbio[rowSums(merged_illumina_pacbio) > 0, ]
    print(head(merged_illumina_pacbio))

    # edgeR basics: 
    group <- factor(c("Illumina","Illumina","PacBio","PacBio")) # Indicate which group each col belongs to
    y <- DGEList(counts=merged_illumina_pacbio, group = group) # Create a DGEList object
    design <- model.matrix(~group)

    # Filter out lowly expressed
    keep <- filterByExpr(y, design)
    y <- y[keep, , keep.lib.sizes=FALSE]

    y <- calcNormFactors(y) # Normalize counts in the object
    y <- estimateDisp(y,design)

    # Pairwise testing approach for DE Genes. "classic" edgeR
    et <- exactTest(y, pair=c("Illumina","PacBio"))

    # Extract exact test table for plotting
    illumina_PB_et <- et$table
    illumina_PB_et$gene_id <- rownames(illumina_PB_et)

    # Adjust p-values
    illumina_PB_et$adj_pval <- p.adjust(illumina_PB_et$PValue, method = "bonferroni")
    illumina_PB_et$status <- as.factor(ifelse(abs(illumina_PB_et$logFC) > 1 & illumina_PB_et$adj_pval <= 0.01,
                             "significant", "not_sig"))

    # MA plot
    ma_plot(illumina_PB_et, fill_color, opt$outdir, dtype, opt$xmax, opt$ymax)

    # Merge the EdgeR table with the other information
    illumina_PB_et <- merge(illumina_PB_et, merged_illumina_pacbio,
                            by.x = "gene_id", by.y = "row.names",
                            all.x = T, all.y = F)

    # Merge in human-readable gene names
    gene_names <- get_gene_names(opt$illumina_kallisto_1)
    illumina_PB_et <- merge(illumina_PB_et, gene_names, by.x = "gene_id",
                            by.y = "g_ID", all.x = T, all.y = F)

    # Merge in length info
    lengths <- get_gene_lengths(opt$illumina_kallisto_1)
    illumina_PB_et <- merge(illumina_PB_et, lengths, by.x = "gene_id",
                            by.y = "g_ID", all.x = T, all.y = F)

    illumina_PB_et <- illumina_PB_et[order(illumina_PB_et$adj_pval),]
    write.table(illumina_PB_et, 
                paste(opt$outdir, "/edgeR_", dtype, "_illumina_gene_counts.tsv", sep=""),
                row.names=F, col.names=T, quote=F, sep = "\t")
}

ma_plot <- function(data, fillcolor, outdir, dtype, xmax, ymax) {

    n_sig <- length(data$status[data$status == "significant"])
    n_no_sig <- length(data$status[data$status == "not_sig"])

    fname <- paste(outdir, "/edgeR_", dtype, "_illumina_gene_counts_MA_plot.png", sep="")
    xlabel <- "log2(Counts per million)"
    ylabel <- paste0(dtype, " to Illumina log2-fold change")

    png(filename = fname,
        width = 2500, height = 2500, units = "px",
        bg = "white",  res = 300)

    print(head(data))
    data$status <- factor(data$status, levels = c("significant", "not_sig"))

    g <- ggplot(data, aes(x=logCPM, y=logFC, color = status)) +
         geom_point(alpha = 0.4, size = 2) +
         xlab(xlabel) + ylab(ylabel) + theme_bw() +
         coord_cartesian(xlim=c(-5,xmax), ylim = c(-1*ymax,ymax)) +
         scale_color_manual(values = c("orange", fillcolor),
                                  labels = c(paste0("Significant (n=", n_sig, ")"),
                                             paste0("Not sig. (n=", n_no_sig, ")"))) +
         theme(axis.text.x = element_text(color="black", size=24),
               axis.text.y = element_text(color="black", size=24),
               axis.title.x = element_text(color="black", size=24),
               axis.title.y = element_text(color="black", size=24)) +
         guides(colour = guide_legend(override.aes = list(alpha=1,size=3))) +
         theme(legend.position=c(0.27,0.08),
                     legend.title = element_blank(),
                     legend.background = element_rect(fill="white", color = "black"),
                     legend.key = element_rect(fill="transparent"),
                     legend.text = element_text(colour = 'black', size = 24))

    print(g)
    dev.off()
}

filter_kallisto_illumina_genes <- function(kallisto_file) {
    # This function takes a Kallisto abundance file and filters the genes
    # based on criteria designed to make the gene set comparable to what can
    # be detected using PacBio

    gencode.quantitation <- as.data.frame(read_delim(kallisto_file, "\t", escape_double = FALSE,
                                  col_names = TRUE, trim_ws = TRUE, na = "NA"))

    # Split GENCODE transcript multi-id by '|'
    extraCols <- str_split_fixed(gencode.quantitation$target_id, "\\|",9)[,c(5,6,8,1,2)]
    colnames(extraCols) <- c("transcript", "gene", "class", "t_ID", "g_ID")
    gencode.quantitation <- cbind(extraCols, gencode.quantitation)

    # Aggregate by gene
    gene_gencode <- gencode.quantitation %>%
                                 dplyr::group_by(g_ID) %>%
                                 dplyr::summarize(counts = sum(est_counts), tpm = sum(tpm))
    colnames(gene_gencode) <- c("g_ID", "count", "tpm")

    gene_gencode <- subset(gene_gencode, tpm > 0)

    return(gene_gencode[,c("g_ID", "count")])
}

get_gene_names <- function(kallisto_file) {
    data <- as.data.frame(read_delim(kallisto_file, "\t", escape_double = FALSE,
                                  col_names = TRUE, trim_ws = TRUE, na = "NA"))

    extraCols <- str_split_fixed(data$target_id, "\\|",9)[,c(5,6,8,1,2)]
    colnames(extraCols) <- c("transcript", "gene", "class", "t_ID", "g_ID")

    return(unique(extraCols[,c("g_ID", "gene")]))

}

get_gene_lengths <- function(kallisto_file) {
    data <- as.data.frame(read_delim(kallisto_file, "\t", escape_double = FALSE,
                                  col_names = TRUE, trim_ws = TRUE, na = "NA"))

    extraCols <- str_split_fixed(data$target_id, "\\|",9)[,c(5,6,8,1,2)]
    colnames(extraCols) <- c("transcript", "gene", "class", "t_ID", "g_ID")
    data <- cbind(extraCols, data)

    gene_gencode <- data %>%
                                 dplyr::group_by(g_ID) %>%
                                 dplyr::summarize(min_length = min(length), 
                                                  max_length = max(length),
                                                  median_length = median(length),
                                                  mean_length = mean(length))   

    return(gene_gencode)
}

load_packages <- function() {
    suppressPackageStartupMessages(library("ggplot2"))
    suppressPackageStartupMessages(library("dplyr"))
    suppressPackageStartupMessages(library("plyr"))
    suppressPackageStartupMessages(library("Hmisc"))
    suppressPackageStartupMessages(library("optparse"))
    suppressPackageStartupMessages(library("readr"))
    suppressPackageStartupMessages(library("reshape"))
    suppressPackageStartupMessages(library("stringr"))
    suppressPackageStartupMessages(library("data.table"))
    suppressPackageStartupMessages(library("preprocessCore"))
    suppressPackageStartupMessages(library("edgeR"))

    return
}

parse_options <- function() {

    option_list <- list(
    make_option(c("--f"), action = "store", dest = "infile",
                    default = NULL, help = "TALON abundance file (not filtered)"),
    make_option(c("--datasets"), action = "store", dest = "datasets",
                    default = NULL, help = "Comma-delimited list of two dataset names to include in the analysis."),
    make_option(c("--ik1"), action = "store", dest = "illumina_kallisto_1",
                    default = NULL, help = "Rep1 Illumina Kallisto file."),
    make_option(c("--ik2"), action = "store", dest = "illumina_kallisto_2",
                    default = NULL, help = "Rep2 Illumina Kallisto file."),
    make_option(c("--color"), action = "store", dest = "color_scheme",
                    default = NULL, help = "blue, red, or green"),
    make_option(c("--xmax"), action = "store", dest = "xmax",
                default = 16, help = "Max x-value for plot (default = 16)"),
    make_option(c("--ymax"), action = "store", dest = "ymax",
                default = 20, help = "Max y-value for plot (default = 20)"),
    make_option(c("-o","--outdir"), action = "store", dest = "outdir",
                default = NULL, help = "Output directory for plots and outfiles"),
    make_option(c("--dtype"), action = "store", dest = "dtype",
                default = "PacBio", help = "Datatype label to display on plot"))
    

    opt <- parse_args(OptionParser(option_list=option_list))
    return(opt)
}


main()
