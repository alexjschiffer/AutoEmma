##################################################
#   Function Name  : Make Input Graphs           #
#   Program Author : Alex Schiffer               #
#   Last Updated   : July 25, 2018               #
#                                                #
#   D'Amato Lab, Boston Children's Hospital      #
##################################################

ae.input.graphs <- function(gfile, pfile, labelsize = 5) {
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
  heatmap <- suppressMessages((ggcorrplot(k) + ggtitle("Kinship Matix") + coord_fixed(ratio = 1)
     + scale_fill_continuous(low = "white", high = "blue",name = "Relatedness\n",limits= c(0,1))
     + theme(plot.title = element_text(hjust = 0.5),
             axis.text.x = element_text(size = labelsize),
             axis.text.y = element_text(size = labelsize),
             legend.title = element_text(size =12),
             legend.text = element_text(size = labelsize),
             legend.margin = margin(t = 1, r = 0, b = 0, l = 0, unit = "pt"),
             legend.key.height=unit(2,"line"))))

  # Create cladogram
  cladogram <- (ggtree(stree, layout = "circular", branch.length = "none") + ggtitle("Population Structure")
    + theme(plot.title = element_text(hjust = 0.5), plot.margin = unit(c(0.5,0.25,0.25,0.25),"cm"))
    + geom_tiplab(size = 2.75, aes(angle=angle)))

  # Initialize vectors for each statistic to be calculated
  count = c()
  avg = c()
  sd = c()
  mad = c()
  counter <- 1

  # Calculate Average, standard deviation, and median absolute deviation for each strain
  for (i in seq_along(strains_vector)) {
    temp <- dplyr::filter(pheno_df, toupper(STRAIN) == strains_vector[i])
    count[counter] <- nrow(temp)
    avg[counter] <- mean(temp$VALUE)
    sd[counter] <- sd(temp$VALUE)
    mad[counter] <- mad(temp$VALUE)
    counter <- counter + 1
  }

  # Assemble averages data and save
  new_data <- dplyr::data_frame(Count = count, Strain = strains_vector, Average = avg, Stdev = sd, Mad = mad)
  write_tsv(x = new_data, path = (paste("results/", tools::file_path_sans_ext(pfile), "-averages.txt", sep = "")), append = FALSE, col_names = TRUE)

  # Create bar graph of average phenotype per strain
  arranged <- arrange(new_data, Average)
  order <- arranged[["Strain"]]
  average <- (ggplot(arranged, aes(x = Strain, y = Average, width = 1)) + geom_bar(stat = "identity", color = "white")
    + theme(axis.text.x = element_text(angle = 45, hjust = 0.95, vjust = 0.95, size = 11))
    + scale_x_discrete(limits = order) + ggtitle("Average Phenotype by Strain"))

  # Create histogram of all values of phenotype
  histogram <- (ggplot(pheno_df, aes(x = pheno_df$VALUE)) + geom_histogram(bins = 20, color = "white")
    + labs(x = "Phenotype Values", y = "Count") + ggtitle("Distribution of Phenotypes"))

  # Create results folder or use an existing one
  cat("\n** Checking for existing results folder.")
  if(!file.exists("results")) {
    cat("\n** -- Creating results folder.")
    dir.create("results")
  } else {
    cat("\n** -- Using existing results folder.")
  }

  # Save all plots
  ggsave(filename = (paste("results/", tools::file_path_sans_ext(pfile), "-heatmap.png", sep = "")), plot = heatmap)
  ggsave(filename = (paste("results/", tools::file_path_sans_ext(pfile), "-cladogram.png", sep = "")), plot = cladogram)
  ggsave(filename = (paste("results/", tools::file_path_sans_ext(pfile), "-histogram.png", sep = "")), plot = histogram)
  ggsave(filename = (paste("results/", tools::file_path_sans_ext(pfile), "-averages.png", sep = "")), plot = average)

  # Save grid
  average2 <- (ggplot(arranged, aes(x = Strain, y = Average, width = 1)) + geom_bar(stat = "identity", color = "white")
    + theme(axis.text.x = element_text(angle = 45, hjust = 0.95, vjust = 0.95, size = 6))
    + scale_x_discrete(limits = order) + ggtitle("Average Phenotype by Strain"))
  cladogram2 <- (ggtree(stree, layout = "circular", branch.length = "none") + ggtitle("Population Structure")
    + theme(plot.title = element_text(hjust = 0.5), plot.margin = unit(c(0.5,0.25,0.25,0.25),"cm"))
    + geom_tiplab(size = 2, aes(angle=angle)))
  grid <- plot_grid(heatmap, cladogram2, histogram, average2)
  ggsave(filename = (paste("results/", tools::file_path_sans_ext(pfile), "-grid.png", sep = "")), plot = grid)
}

  # For Future Development
  #if (color) {
    # Create categories
    #pheno_data <- dplyr::mutate(pheno_data, cat = 1)
    #strain_data$cat["Average" < 1] <- 2
    #strain_data$cat[Average >= 1] <- 2
    #pheno_cats <- strain_data[["cat"]]

    # Create color spectrum
    #cs <- colorRamps::blue2red(length(pheno_cats))

    # Map color spectrum onto phenotype data
    # colors_vecs <- c()
    #for (x in c(1:length(stree$tip.label))) {
    #  category <- pheno_cats[x]
    #  colors_vecs[x] <- colors[category]
    #}
  #}
