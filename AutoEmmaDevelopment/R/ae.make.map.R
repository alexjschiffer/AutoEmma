##################################################
#   Function Name  : Run AutoEmma                #
#   Program Author : Alex Schiffer               #
#   Last Updated   : June 28, 2018               #
#                                                #
#   D'Amato Lab, Boston Children's Hospital      #
##################################################

ae.make.map <- function(file){

  # Check if required packages are installed and load them
  ae.check.required()

  # Load in list of strains and extract unique strains
  strains <- ae.import(file)
  if (!is.null(strains)){
    colnames(strains) <- toupper(colnames(strains))
    strains <- dplyr::select(strains, "STRAIN")
    unique <- dplyr::distinct(strains)
    unique <- unique[[1]]
  } else {
    message("Unable to import phenotype file.")
  }

  # Create new table
  new_table <- matrix(1:((length(unique))+1),1:((length(unique))+1))
  colnames(new_table) <- c("STRAIN", unique)

  # Expand table to include all unique strains as column names
  for (i in 1:(length(unique))) {
    strains <- dplyr::mutate(strains, !!unique[i])
  }
  colnames(strains) <- colnames(new_table)

  # Create Map
  for (i in 1:(dim(strains)[1])){
    row <- strains[i,]
    value <- as.character(row[1])
    row[colnames(row) == value] <- 1
    new_table <- rbind(new_table, row)
  }

  # Clean and format new table
  new_table[new_table != 1] <- 0
  new_table <- new_table[2:dim(new_table)[1],1:length(unique)+1]
  new_table <- data.matrix(new_table, rownames.force = NA)

  return(new_table)
}
