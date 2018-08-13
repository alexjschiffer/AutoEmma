##################################################
#   Function Name  : Make Input Graphs           #
#   Program Author : Alex Schiffer               #
#   Last Updated   : July 25, 2018               #
#                                                #
#   D'Amato Lab, Boston Children's Hospital      #
##################################################

ae.input.graphs <- function(gfile, pfile, stdev = FALSE, labelsize = 10) {
  # Announce start of function
  cat("* Creating diagnostic graphs for input files.\n")

  # Load required packages
  ae.check.required()

  # Import phenotypes as data frame
  pheno_df <- ae.import(pfile)

  # Import genotypes as data frame
  geno_df <- ae.import(gfile)

  # Check input files
  if (is.null(pheno_df)){
    print("Invalid Phenotype File.")
    return()
  }
  if (is.null(geno_df)){
    print("Invalid Genotype File.")
    return()
  }

  # Extract strain names and values
  colnames(pheno_df) <- toupper(colnames(pheno_df))
  strains_table <- pheno_df %>%
    dplyr::select("STRAIN") %>%
    dplyr::distinct()
  strains_vector <- toupper(strains_table[["STRAIN"]])
  values <- pheno_df[["VALUE"]]

  # Extract SNP information from genotype data frame
  colnames(geno_df) <- toupper(colnames(geno_df))
  geno_matrix <- geno_df[,4:dim(geno_df)[2]]
  geno_matrix <- dplyr::select(geno_matrix, one_of(strains_vector))
  geno_matrix <- data.matrix(geno_matrix)

  # Create transposed SNP matrix (for gene distance calculation)
  trans_geno_matrix <- t(geno_matrix)

  # Create kinship matrix
  k <- emma.kinship(geno_matrix)
  colnames(k) <- strains_vector
  rownames(k) <- strains_vector

  # Generate tree based on gene distance
  stree <- trans_geno_matrix %>%
    dist.gene(method = "pairwise") %>%
    nj()

  # Create heatmap
  heatmap <- suppressMessages((ggcorrplot(k) + ggtitle("Kinship Matrix") + coord_fixed(ratio = 1) +
                                 scale_fill_continuous(low = "yellow", high = "red",name = "Relatedness\n",limits= c(0,1)) +
                                 theme(plot.title = element_text(hjust = 0.5), legend.title = element_text(size = 12),
                                       legend.text = element_text(size = labelsize),
                                       legend.key.height=unit(2,"line"), axis.text.x=element_blank())))

  # Create cladogram
  cladogram <- (ggtree(stree, layout = "circular", branch.length = "none") + ggtitle("Population Structure") +
                  theme(plot.title = element_text(hjust = 0.5)) +
                  geom_tiplab(size = 2.75, aes(angle=angle)))

  # Summarise averages for phenotype data
  new_data <- pheno_df %>%
    group_by(STRAIN) %>%
    summarise(mean = mean(VALUE), stdev = sd(VALUE), mad = mad(VALUE), sem = sd(VALUE)/sqrt(length(VALUE)), count = n())

  # Create results folder or use an existing one
  cat("\n** Checking for existing results folder.")
  if(!file.exists("results")) {
    cat("\n** -- Creating results folder.")
    dir.create("results")
  } else {
    cat("\n** -- Using existing results folder.")
  }

  # Save averages spreadsheet
  write_tsv(x = new_data, path = (paste("results/", tools::file_path_sans_ext(pfile), "-averages.txt", sep = "")), append = FALSE, col_names = TRUE)

  # Create bar graph of average phenotype per strain
  if(stdev){
    average <- ggplot(data = new_data, mapping = aes(x = reorder(STRAIN,mean), y = mean)) +
      geom_col(width = 1, color = "white", fill = "royalblue3") +
      theme(axis.text.x = element_text(angle = 45, hjust = 0.95, vjust = 0.95, size = labelsize)) +
      ggtitle("Average Phenotype by Strain") +
      labs(x = "Strain", y = "Average") +
      geom_errorbar(aes(ymin = mean-stdev, ymax = mean+stdev))
  } else {
    average <- ggplot(data = new_data, mapping = aes(x = reorder(STRAIN,mean), y = mean)) +
      geom_col(width = 1, color = "white", fill = "royalblue3") +
      theme(axis.text.x = element_text(angle = 45, hjust = 0.95, vjust = 0.95, size = labelsize)) +
      ggtitle("Average Phenotype by Strain") +
      labs(x = "Strain", y = "Average") +
      geom_errorbar(aes(ymin = mean-sem, ymax = mean+sem)) +
      scale_y_continuous(expand = c(0,0))
  }

  # Create histogram of all values of phenotype
  histogram <- ggplot(pheno_df, aes(x = pheno_df$VALUE)) +
    geom_histogram(bins = 20, color = "white", fill = "royalblue3") +
    labs(x = "Phenotype Values", y = "Count") +
    ggtitle("Distribution of Phenotypes") +
    scale_y_continuous(expand = c(0,0))

  # Save all plots
  ggsave(filename = (paste("results/", tools::file_path_sans_ext(pfile), "-heatmap.png", sep = "")), plot = heatmap)
  ggsave(filename = (paste("results/", tools::file_path_sans_ext(pfile), "-cladogram.png", sep = "")), plot = cladogram)
  ggsave(filename = (paste("results/", tools::file_path_sans_ext(pfile), "-histogram.png", sep = "")), plot = histogram)
  ggsave(filename = (paste("results/", tools::file_path_sans_ext(pfile), "-averages.png", sep = "")), plot = average)

}
