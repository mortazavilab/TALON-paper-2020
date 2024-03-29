
main <-function() {

    set.seed(100)
    load_packages()
    opt <- parse_options()

    # Get colors
    if (opt$color_scheme == "red") {
        color_vec <- c("white", "orange", "red2")
    } else if (opt$color_scheme == "blue") {
        color_vec <- c("white", "skyblue", "navy")
    } else if (opt$color_scheme == "green") {
        color_vec <- c("white", "yellowgreen", "springgreen4")
    }

    # Get genes and transcripts expressed in the Illumina data from the Kallisto
    # abundance file
    illumina_1 <- filter_kallisto_illumina_genes(opt$illumina_kallisto_1)
    illumina_2 <- filter_kallisto_illumina_genes(opt$illumina_kallisto_2)
    colnames(illumina_1) <- c("gene", "illumina_TPM_1")
    colnames(illumina_2) <- c("gene", "illumina_TPM_2")
    illumina_gene_table <- merge(illumina_1, illumina_2, by = "gene",
                                 all.x = T, all.y = T)

    # Remove genes only detected in one Illumina rep, then average the TPMS
    illumina_gene_table <- illumina_gene_table[complete.cases(illumina_gene_table), ]
    illumina_gene_table$tpm <- (illumina_gene_table$illumina_TPM_1 +
                                illumina_gene_table$illumina_TPM_2) / 2
    print(paste0("Illumina genes: ", nrow(illumina_gene_table)))

    
    # Get the names of the first and second dataset that we will be working with
    data_names <- str_split(opt$datasets, ",")[[1]]
    dataset_1 <- data_names[1]
    dataset_2 <- data_names[2]

    # Read in abundance file
    abundance_table <- as.data.frame(read_delim(opt$infile, delim = "\t",
                                  col_names = TRUE, trim_ws = TRUE, na = "NA"))
    # Remove genomic transcripts
    abundance_table <- subset(abundance_table, transcript_novelty != "Genomic")

    # Get genes in each dataset
    d1_genes <- get_detected_genes_for_dataset(abundance_table, dataset_1)
    d2_genes <- get_detected_genes_for_dataset(abundance_table, dataset_2)
    

    # Combine the Illumina tables with information about which genes/transcripts are observed in Pacbio
    illumina_gene_detection <- get_detection(illumina_gene_table, d1_genes, d2_genes, "gene", opt$dtype)

    # Group into buckets
    illumina_gene_detection_buckets <- get_buckets(illumina_gene_detection)    

    # Plot detection by TPM
    plot_detection(illumina_gene_detection_buckets$illumina, illumina_gene_detection_buckets$interval_labels, "gene", color_vec, opt$outdir, opt$dtype)
    print ("--------------------------")

}

get_buckets <- function(illumina) {
    intervals = c(0, 1, 2, 5, 10, 50, 100, 500)
    illumina$group = cut(illumina$tpm, c(intervals, 100000000000))   

    cIntervals = as.character(intervals)
    interval_starts_with_dash = paste0(cIntervals, c(rep("-", length(cIntervals)-1), "+"))
    full_intervals = paste0(interval_starts_with_dash, c(as.character(intervals)[2:length(intervals)], ""))

    return(list("illumina" = illumina, "interval_labels" = full_intervals))
}


get_detection <- function(illumina, d1_items, d2_items, cat_type, dtype) {

    shared <- intersect(d1_items, d2_items)
    d1_d2_union <- union(d1_items, d2_items) 
    not_shared <- d1_d2_union[(d1_d2_union %in% shared) == F]

    illumina$detection = ifelse(illumina[,cat_type] %in% not_shared, paste0("Detected in one ", dtype, " rep"),
                                             ifelse(illumina[,cat_type] %in% shared, paste0("Detected in both ", dtype, " reps"), 
                                                                                     paste0("Not detected in ", dtype)))

    return(illumina)
}

plot_detection <- function(illumina, cIntervals, cat_type, color_vec, outdir, dtype) {

    # Compute number of genes per bin
    totals <- illumina %>%
        dplyr::group_by(group) %>%
        dplyr::tally()
    print(totals)

    # Merge in to detection table
    illumina <- merge(illumina, totals, by = "group", all.x = T, all.y = T)

    # Plot 
    fname <- paste(outdir, "/", cat_type, "_detection_by_TPM.png", sep="")
    xlabel <- paste(capitalize(cat_type), "expression level in Illumina data (TPM)")
    ylabel <- paste("Fraction of ", cat_type, "s", sep="")

    png(filename = fname,
        width = 2500, height = 2500, units = "px",
        bg = "white",  res = 300)
    g = ggplot(illumina, aes(x = group, fill = factor(detection, levels = c(paste0("Not detected in ", dtype),
                                                                            paste0("Detected in one ", dtype, " rep"), 
                                                                            paste0("Detected in both ", dtype, " reps"))))) +
       geom_bar(position = "fill", col = "black") + 
       scale_fill_manual("",values=color_vec)  + 
       theme_bw() +
       coord_cartesian(ylim = c(0, 1.05)) + 
       scale_y_continuous(breaks=c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)) + 
       theme(axis.text.x = element_text(angle = 30, hjust = 1, color = "black", size=26)) + 
       xlab(xlabel)  + ylab(ylabel) + theme(text= element_text(size=26)) + 
       theme(axis.text.x = element_text(color = "black", size=24),
             axis.text.y = element_text(color = "black", size=24)) + 
       scale_x_discrete(labels=cIntervals) + 
       theme(legend.position=c(0.6,0.2),
              legend.title = element_blank(),
              legend.background = element_rect(fill="white", color = "black"),
              legend.key = element_rect(fill="transparent"),
              legend.text = element_text(colour = 'black', size = 24),
              plot.margin = margin(t = 0.5, r = 1.5, l = 0.5, b = 0.5, "cm")) +
        geom_text(aes(x= group, y = 1, label=n), vjust=-1, size = 7)
    print(g)
    dev.off()
}


get_detected_genes_for_dataset <- function(abundance, dataset) {
    data <- abundance[,c("annot_gene_id", dataset)]
    observed_genes <- unique(data[data[,dataset] > 0, "annot_gene_id"])

    return(observed_genes)
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

    # Aggregate on genes
    gene_gencode <- gencode.quantitation %>%
                                 dplyr::group_by(g_ID) %>%
                                 dplyr::summarize(tpm = sum(tpm))
    print(nrow(gene_gencode))
    gene_gencode <- subset(gene_gencode, tpm > 0)
    print(nrow(gene_gencode))

    return(gene_gencode)

}

load_packages <- function() {
    suppressPackageStartupMessages(library("DBI"))
    suppressPackageStartupMessages(library("ggplot2"))
    suppressPackageStartupMessages(library("dplyr"))
    suppressPackageStartupMessages(library("plyr"))
    suppressPackageStartupMessages(library("Hmisc"))
    suppressPackageStartupMessages(library("optparse"))
    #suppressPackageStartupMessages(library("RSQLite"))
    suppressPackageStartupMessages(library("tidyverse"))
    suppressPackageStartupMessages(library("reshape"))
    suppressPackageStartupMessages(library("stringr"))
    suppressPackageStartupMessages(library("data.table"))

    return
}


parse_options <- function() {

    option_list <- list(
        make_option(c("--f"), action = "store", dest = "infile",
                    default = NULL, help = "Unfiltered TALON abundance file"),
        make_option(c("--datasets"), action = "store", dest = "datasets",
                    default = NULL, help = "Comma-delimited list of two dataset names to include in the analysis."),
        make_option(c("--ik1"), action = "store", dest = "illumina_kallisto_1",
                    default = NULL, help = "Rep1 Illumina Kallisto file."),
        make_option(c("--ik2"), action = "store", dest = "illumina_kallisto_2",
                    default = NULL, help = "Rep2 Illumina Kallisto file."),
        make_option(c("--color"), action = "store", dest = "color_scheme",
                    default = NULL, help = "blue, red, or green"),
        make_option(c("--dtype"), action="store", dest='dtype',
                    help = 'Platform used to generate the data ie PacBio or ONT'),
        make_option(c("-o","--outdir"), action = "store", dest = "outdir",
                    default = NULL, help = "Output directory for plots and outfiles")
        )

    opt <- parse_args(OptionParser(option_list=option_list))
    return(opt)
}


main()
