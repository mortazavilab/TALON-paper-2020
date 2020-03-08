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
        fill_color <- "#009E73"
    }

    # Get the names of the first and second dataset that we will be working with
    data_names <- str_split(opt$datasets, ",")[[1]]
    dataset1 <- data_names[1]
    dataset2 <- data_names[2]

    # datatype label (ONT or PacBio)
    dtype <- opt$dtype

    # Get transcripts expressed in the Illumina data from the Kallisto
    # abundance files
    illumina_1 <- filter_kallisto_illumina_transcripts(opt$illumina_kallisto_1)
    illumina_2 <- filter_kallisto_illumina_transcripts(opt$illumina_kallisto_2)
    colnames(illumina_1) <- c("annot_transcript_id", "illumina_counts_1")
    colnames(illumina_2) <- c("annot_transcript_id", "illumina_counts_2")
    illumina_table <- merge(illumina_1, illumina_2, by = "annot_transcript_id",
                                 all.x = T, all.y = T)
    print(paste0("Illumina transcripts: ", nrow(illumina_table)))
    illumina_table[is.na(illumina_table)] <- 0
    illumina_table$illumina_counts_1 <- round(illumina_table$illumina_counts_1)
    illumina_table$illumina_counts_2 <- round(illumina_table$illumina_counts_2)

    # Read PacBio abundance file
    pb_abundance_orig <- as.data.frame(read_delim(opt$infile, "\t", escape_double = FALSE,
                                  col_names = TRUE, trim_ws = TRUE, na = "NA"))

    # Keep known transcripts only
    pb_abundance <- subset(pb_abundance_orig, transcript_novelty == "Known")

    # Cut out unnecessary cols
    pb_abundance <- pb_abundance[, c("annot_transcript_id", dataset1, dataset2)]


    # Merge PacBio with Illumina on annot_transcript_name
    merged_illumina_pacbio <- merge(illumina_table, pb_abundance, by = "annot_transcript_id",
                                    all.x = T, all.y = T)

    merged_illumina_pacbio[is.na(merged_illumina_pacbio)] <- 0
    merged_illumina_pacbio <- merged_illumina_pacbio[, c("annot_transcript_id", 
                                                         "illumina_counts_1", 
                                                         "illumina_counts_2",  
                                                         dataset1, dataset2)]

    # Now, format table for edgeR by setting the rownames to the gene names
    rownames(merged_illumina_pacbio) <- merged_illumina_pacbio$annot_transcript_id
    merged_illumina_pacbio$annot_transcript_id <- NULL

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
    illumina_PB_et$transcript_id <- rownames(illumina_PB_et)

    # Adjust p-values
    illumina_PB_et$adj_pval <- p.adjust(illumina_PB_et$PValue, method = "bonferroni")
    illumina_PB_et$status <- as.factor(ifelse(abs(illumina_PB_et$logFC) > 1 & illumina_PB_et$adj_pval <= 0.01,
                             "significant", "not_sig"))

    # MA plot
    ma_plot(illumina_PB_et, fill_color, opt$outdir, dtype, opt$xmax, opt$ymax)

    # Merge the EdgeR table with the other information
    illumina_PB_et <- merge(illumina_PB_et, merged_illumina_pacbio, 
                            by.x = "transcript_id", by.y = "row.names",
                            all.x = T, all.y = F)

    # Merge in human-readable gene names
    transcript_names <- get_transcript_names(opt$illumina_kallisto_1)
    illumina_PB_et <- merge(illumina_PB_et, transcript_names, by.x = "transcript_id",
                            by.y = "t_ID", all.x = T, all.y = F)

    # Merge in transcript lengths
    transcript_lengths <- get_transcript_lengths(opt$illumina_kallisto_1)
    illumina_PB_et <- merge(illumina_PB_et, transcript_lengths, by.x = "transcript_id",
                            by.y = "t_ID", all.x = T, all.y = F)

    illumina_PB_et <- illumina_PB_et[order(illumina_PB_et$adj_pval),]
    write.table(illumina_PB_et, 
                paste(opt$outdir, "/edgeR_", dtype, "_illumina_transcript_counts.tsv", sep=""),
                row.names=F, col.names=T, quote=F, sep = "\t")
}

ma_plot <- function(data, fillcolor, outdir, dtype, xmax, ymax) {

    n_sig <- length(data$status[data$status == "significant"])
    n_no_sig <- length(data$status[data$status == "not_sig"])

    fname <- paste(outdir, "/edgeR_", dtype, "_illumina_transcript_counts_MA_plot.png", sep="")
    xlabel <- "log2(Counts per million)"
    ylabel <- paste0(dtype, " to Illumina log2-fold change")

    png(filename = fname,
        width = 2500, height = 2500, units = "px",
        bg = "white",  res = 300)

    print(head(data))

    g <- ggplot(data, aes(x=logCPM, y=logFC, color = status)) +
         geom_point(alpha = 0.4, size = 2) +
         xlab(xlabel) + ylab(ylabel) + theme_bw() +
         coord_cartesian(xlim=c(-5,xmax), ylim = c(-1*ymax,ymax)) +
         scale_color_manual(values = c("orange", fillcolor),
                                  labels = c(paste0("Significant (n = ", n_sig, ")"),
                                             paste0("Not significant (n = ", n_no_sig, ")"))) +
         theme(axis.text.x = element_text(color="black", size=22),
               axis.text.y = element_text(color="black", size=22),
               axis.title.x = element_text(color="black", size=22),
               axis.title.y = element_text(color="black", size=22)) +
         guides(colour = guide_legend(override.aes = list(alpha=1,size=2.5))) +
         theme(legend.position=c(0.25,0.1),
                     legend.title = element_blank(),
                     legend.background = element_rect(fill="white", color = "black"),
                     legend.key = element_rect(fill="transparent"),
                     legend.text = element_text(colour = 'black', size = 18))

    print(g)
    dev.off()
}

filter_kallisto_illumina_transcripts <- function(kallisto_file) {
    # This function takes a Kallisto abundance file and filters the transcripts
    # based on criteria designed to make the transcript set comparable to what can
    # be detected using PacBio

    gencode.quantitation <- as.data.frame(read_delim(kallisto_file, "\t", escape_double = FALSE,
                                  col_names = TRUE, trim_ws = TRUE, na = "NA"))

    # Split GENCODE transcript multi-id by '|'
    extraCols <- str_split_fixed(gencode.quantitation$target_id, "\\|",9)[,c(5,6,8,1,2)]
    colnames(extraCols) <- c("transcript", "gene_name", "class", "t_ID", "g_ID")
    gencode.quantitation <- cbind(extraCols, gencode.quantitation)

    gencode.quantitation <- subset(gencode.quantitation, tpm > 0)

    return(gencode.quantitation[,c("t_ID", "est_counts")])
}

get_transcript_names <- function(kallisto_file) {
    data <- as.data.frame(read_delim(kallisto_file, "\t", escape_double = FALSE,
                                  col_names = TRUE, trim_ws = TRUE, na = "NA"))

    extraCols <- str_split_fixed(data$target_id, "\\|",9)[,c(5,6,8,1,2)]
    colnames(extraCols) <- c("transcript", "gene", "class", "t_ID", "g_ID")

    return(unique(extraCols[,c("t_ID", "transcript")]))
}

get_transcript_lengths <- function(kallisto_file) {
    data <- as.data.frame(read_delim(kallisto_file, "\t", escape_double = FALSE,
                                  col_names = TRUE, trim_ws = TRUE, na = "NA"))

    extraCols <- str_split_fixed(data$target_id, "\\|",9)[,c(5,6,8,1,2)]
    colnames(extraCols) <- c("transcript", "gene", "class", "t_ID", "g_ID")
    data <- cbind(extraCols, data)

    return(data[,c("t_ID", "length")])

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
