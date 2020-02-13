main <-function() {

    load_packages()
    opt <- parse_options()
    datasets <- opt$datasets
    outdir <- opt$outdir

    # Read abundance table
    abundance_table <- as.data.frame(read_delim(opt$infile, delim = "\t",
                                  col_names = TRUE, trim_ws = TRUE, na = "NA"))

    # Filter the results to limit them to selected datasets
    datasets <- unlist(strsplit(datasets, ","))
    print(datasets)
    abundance_table <- abundance_table[,c("transcript_novelty", datasets)]
    
    # Remove transcripts if they do not appear in any of the datasets
    if(length(datasets) > 1) {
        abundance_table$total <- rowSums(abundance_table[,datasets])
    } else {
        abundance_table$total <- abundance_table[,datasets]
    }
    abundance_table <- subset(abundance_table, total > 0) 

    # Plot
    plot_distinct_novelty(abundance_table, outdir, datasets, opt$ymax)
}

plot_distinct_novelty <- function(distinct_transcripts, outdir, datasets, ymax){
    # This function plots the number of transcripts that belong to each 
    # transcript novelty class. So it amounts to the number of unique transcripts
    # of each type that were identified in the datasets together

    distinct_transcripts$novelty <- factor(distinct_transcripts$transcript_novelty,
                                           levels = c("Known", "ISM", "NIC", "NNC",
                                                      "Antisense", "Intergenic",
                                                      "Genomic"))

    total <- nrow(distinct_transcripts)

    # Compute percentages
    freqs_by_novelty <- distinct_transcripts %>% group_by(novelty) %>% count() %>% 
                        mutate(percent = round(100*n/total,1))
    print(freqs_by_novelty)
    distinct_transcripts <- merge(distinct_transcripts, freqs_by_novelty,
                                  by = c("novelty"), all.x = T, all.y = F)

    # Plotting
    str_datasets <- paste(datasets, collapse='-')
    fname <- paste(outdir, "/", str_datasets, "_distinct_isoforms_by_category.png", sep="")
    xlabel <- "Category"
    ylabel <- "Number of distinct transcript models"
    #ymax = 10^4.6
    if (is.null(ymax)) {
        ymax <- 1.35*(max((distinct_transcripts %>% group_by(novelty) %>% count())$n))
    } else {
        ymax <- as.numeric(ymax)
    }

    colors <- c("Known" = "#009E73","ISM" = "#0072B2", "NIC" = "#D55E00", 
                "NNC" = "#E69F00", "Antisense" = "#000000", 
                "Intergenic" = "#CC79A7", "Genomic" = "#F0E442")

    png(filename = fname,
        width = 3100, height = 3500, units = "px",
        bg = "white",  res = 300)
    total_isoforms = nrow(distinct_transcripts)
    g = ggplot(distinct_transcripts, aes(x = novelty, width=.6,
               fill = novelty)) + 
               geom_bar() + 
               ylab(ylabel) +
               scale_fill_manual("Transcript Type", values = colors) +
               theme_bw(base_family = "Helvetica", base_size = 20) +
               theme(axis.line.x = element_line(color="black", size = 0.5),
                     axis.line.y = element_line(color="black", size = 0.5),
                     axis.text.x = element_text(color="black", size = rel(2.5),
                                             angle = 25, vjust = 1, hjust=1),
                     axis.text.y = element_text(color="black", size = rel(2.5)),
                     axis.title.x = element_blank(),
                     axis.title.y = element_text(color="black", size=rel(2))) +
               guides(fill = FALSE) + 
               coord_cartesian(ylim = c(1, ymax)) +
               yscale("log10", .format = TRUE) + 
               geom_text(aes(y = ..count..,
                  label = paste0(percent, '%')),
                  stat = 'count',
                  position = position_dodge(.9),
                  size = rel(11), vjust=-0.25) +
               annotate("text", x = 3.5, y = ymax*0.9,
                          label = paste("Total isoforms:", total_isoforms), 
                          color="black", size = 10) 
                

    print(g)
    dev.off()
    print(total_isoforms)
}


load_packages <- function() {
    suppressPackageStartupMessages(library("readr"))
    suppressPackageStartupMessages(library("ggplot2"))
    #suppressPackageStartupMessages(library("plyr"))
    suppressPackageStartupMessages(library("optparse"))
    suppressPackageStartupMessages(library("ggpubr"))
    suppressPackageStartupMessages(library("dplyr"))
  
    return
}

parse_options <- function() {

    option_list <- list(
        make_option(c("--f"), action = "store", dest = "infile",
                    default = NULL, help = "Filtered TALON abundance file"),
        make_option(c("--datasets"), action = "store", dest = "datasets",
                    default = NULL, help = "Comma-separated list of datasets to include"),
        make_option(c("--ymax"), action = "store", dest = "ymax",
                    default = NULL, help = "Optional: max value for y axis."),
        make_option(c("-o","--outdir"), action = "store", dest = "outdir",
                    default = NULL, help = "Output directory for plots and outfiles")
        )

    opt <- parse_args(OptionParser(option_list=option_list))
    return(opt)
}

main()
