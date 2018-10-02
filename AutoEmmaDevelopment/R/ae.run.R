##################################################
#   Function Name  : Run AutoEmma                #
#   Program Author : Alex Schiffer               #
#   Last Updated   : June 22, 2018               #
#                                                #
#   D'Amato Lab, Boston Children's Hospital      #
##################################################

ae.run <- function(genotype = "", phenotype = "", kinship = NULL, map = NULL, method = "REML", use_snps = "all") {

  # Check if required packages are installed and load them
  cat("** Starting AutoEmma run.\n** -- Ensuring all required packages are installed and loaded.")
  ae.check.required()

  # Check if valid method is selected
  if (method != "REML" & method != "ML") {
    stop("INVALID METHOD SELECTED.")
  }

  # Check for valid phenotype file
  cat("\n** -- Attempting to import phenotype file.")
  phenotype_data <- ae.import(phenotype)

  if (!is.null(phenotype_data)) {
    colnames(phenotype_data) <- toupper(colnames(phenotype_data))
    groups <- toupper(dplyr::distinct(dplyr::select(phenotype_data, "STRAIN"))[["STRAIN"]])
    phenotype_matrix <- phenotype_data[["VALUE"]]
    phenotype_matrix <- t(phenotype_matrix)
    phenotype_matrix <- data.matrix(phenotype_matrix)
    cat("\n** -- Successfully imported phenotype file.")
    cat("\n** --", length(groups), "distinct strains detected.")
  } else {
    stop("UNABLE TO IMPORT PHENOTYPE FILE.")
  }

  # Check for valid genotype file
  cat("\n** -- Attempting to import genotype file.")
  genotype_data <- ae.import(genotype)

  if (!is.null(genotype_data)) {

    # Capitalize all column headers
    colnames(genotype_data) <- toupper(colnames(genotype_data))

    # Extract out snp information from genotype data
    columns <- dplyr::select(genotype_data, one_of(c("SNP", "CHR", "BP")))

    # Extract out snps by strain and convert to numerical matrix
    # snps <- genotype_data[,4:dim(genotype_data)[2]]
    snps <- dplyr::select(genotype_data, one_of(groups))
    genotype_matrix <- data.matrix(snps, rownames.force = NA)
    cat("\n** -- Successfully imported genotype file.")
    cat("\n** --", dim(genotype_matrix)[1], "snps for", dim(genotype_matrix)[2], "strains loaded.")
  } else {
    stop("UNABLE TO IMPORT GENOTYPE FILE.")
  }

  # Generate IBS kinship matrix
  if (is.null(kinship)){
    cat("\n** Creating IBS kinship matrix.\n")
    kinship_matrix <- emma.kinship(genotype_matrix, "additive", use = use_snps)
    cat("\n")
  } else {
    cat("\n** Using provided kinship matrix\n")
    kinship_matrix <- kinship
  }

  # Run GWAS EMMA
  cat("\n** Starting EMMA run.\n")
  start <- as.numeric(Sys.time())
  if(is.null(map)){
    if (method == "REML") {
      results <- emma.REML.t(ys = phenotype_matrix, xs = genotype_matrix, K = kinship_matrix)
    } else if (method == "ML") {
      results <- emma.ML.LRT(ys = phenotype_matrix, xs = genotype_matrix, K = kinship_matrix)
    }
  } else {
    if (method == "REML") {
      results <- emma.REML.t(phenotype_matrix, genotype_matrix, kinship_matrix, map)
    } else if (method == "ML") {
      results <- emma.ML.LRT(phenotype_matrix, genotype_matrix, kinship_matrix, map)
    }
  }
  end <- as.numeric(Sys.time())
  cat("\n** Run finished successfully.")
  cat("\n** --GWAS EMMA took ", round(((end - start)/60), 3), " minutes.", sep = "")

  # Extract results data
  pvalues <- as.data.frame(results$ps)
  var_gene <- as.data.frame(results$vgs)
  var_rand <- as.data.frame(results$ves)

  # Label results headers
  colnames(pvalues) <- c("P")
  colnames(var_gene) <- c("SNP Genetic Variance")
  colnames(var_rand) <- c("SNP Environmental Variance")

  # Assemble results data
  manhattan_data <- dplyr::bind_cols(columns, pvalues, var_gene, var_rand)

  # Create results folder or use an existing one
  cat("\n** Checking for existing results folder.")
  if(!file.exists("results")) {
    cat("\n** -- Creating results folder.")
    dir.create("results")
  } else {
    cat("\n** -- Using existing results folder.")
  }
  cat("\n** Saving results.")

  ################################################################################
  ##               Create and save all plots and spreadsheets                   ##
  ################################################################################

  # Save results spreadsheet
  write_csv(manhattan_data, (paste("results/", tools::file_path_sans_ext(phenotype), "-results.csv", sep = "")))

  # Make Manhattan, QQ, and Regional Plots
  cat("\n** -- Saving manhattan plot.")
  ae.manhattan(manhattan_data,(paste("results/", tools::file_path_sans_ext(phenotype), "manhattan.png", sep = "")),
               title = (paste(tools::file_path_sans_ext(phenotype), "Manhattan Plot", sep = " ")))
  cat("\n** -- Saving Q-Q plot.")
  ae.qq(manhattan_data, (paste("results/", tools::file_path_sans_ext(phenotype), "qq.png", sep = "")),
        title = (paste(tools::file_path_sans_ext(phenotype), "Q-Q Plot", sep = " ")))
  cat("\n** -- Saving regional manhattan plots.")
  ae.make.regional(manhattan_data, (paste("results/", tools::file_path_sans_ext(phenotype), "regional.pdf", sep = "")))

  # Make Input Diagnostic Graphs
  cat("\n** -- Saving input diagnostic graphs")
  ae.input.graphs(genotype, phenotype)

  return(manhattan_data)
}
