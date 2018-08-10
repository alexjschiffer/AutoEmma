##################################################
#   Function Name  : Convert From atcg to binary #
#   Program Author : Alex Schiffer               #
#   Last Updated   : June 19, 2018               #
#                                                #
#   D'Amato Lab, Boston Children's Hospital      #
##################################################

ae.recode <- function(input){
  # Check if required packages are installed and load them
  ae.check.required()

  # Convert All Lowercase to Uppercase
  cat("Converting lowercase SNPs to uppercase.\n")
  pb1 <- txtProgressBar(min = 0, initial = 0, max = 6, style = 3)
  input[input == "a"] <- "A"
  Sys.sleep(0.001)
  setTxtProgressBar(pb1, 1)
  input[input == "t"] <- "T"
  Sys.sleep(0.001)
  setTxtProgressBar(pb1, 2)
  input[input == "c"] <- "C"
  Sys.sleep(0.001)
  setTxtProgressBar(pb1, 3)
  input[input == "g"] <- "G"
  Sys.sleep(0.001)
  setTxtProgressBar(pb1, 4)
  input[input == "n"] <- "N"
  Sys.sleep(0.001)
  setTxtProgressBar(pb1, 5)
  input[input == "h"] <- "H"
  Sys.sleep(0.001)
  setTxtProgressBar(pb1, 6)

  # Extract out SNP information
  snps <- input[["ALLELES"]]
  input <- dplyr::select(input, -one_of(c("ALLELES")))
  ALLELE1 <- c(1:length(snps))
  ALLELE2 <- c(1:length(snps))
  for (i in seq_along(snps)) {
    ALLELE1[i] <- unlist(strsplit(snps[i],""))[1]
    ALLELE2[i] <- unlist(strsplit(snps[i],""))[3]
  }
  alleledf <- data.frame(ALLELE1, ALLELE2)

  # Recode snps
  cat("\nRecoding SNPS from Letters to Binary.\n")
  pb <- txtProgressBar(min = 0, initial = 0, max = length(input), style = 3)
  for (i in seq_along(input)){

    # Update Progress Bar
    if(i %% 10 == 0) {
      Sys.sleep(0.001)
      setTxtProgressBar(pb, i)
    }

    input[i][input[i] == alleledf["ALLELE1"]] <- 0
    input[i][input[i] == alleledf["ALLELE2"]] <- 1
  }

  input[input == "N"] <- "NA"

  readr::write_tsv(input, "recoded.txt")
  return(input)
}
