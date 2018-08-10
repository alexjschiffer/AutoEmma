##################################################
#   Function Name  : Import Data                 #
#   Program Author : Alex Schiffer               #
#   Last Updated   : June 28, 2018               #
#                                                #
#   D'Amato Lab, Boston Children's Hospital      #
##################################################

ae.import <- function(input, unique = FALSE) {
  # Check if required packages are installed and load them
  ae.check.required()

  # Hide messages about column interpretation
  options(readr.num_columns = 0)

  # Check if file exists and is either TSV of CSV
  if (!file.exists(input)){
    stop("No file found.")
  } else {
    # Load in file as data table
    if (tools::file_ext(input) == "txt" | tools::file_ext(input) == "tab" | tools::file_ext(input) == "tsv"){
      suppressWarnings(data <- readr::read_tsv(input, col_names = TRUE, trim_ws = TRUE))
      if(unique){
        data <- dplyr::distinct(data)
      }
      data <- dplyr::filter(data, !is.na(data[1]) & !is.na(data[2]))
      colnames(data) <- toupper(colnames(data))
      return(data)
    } else if (tools::file_ext(input) == "csv") {
      suppressWarnings(data <- readr::read_csv(input, col_names = TRUE, trim_ws = TRUE))
      if(unique){
        data <- dplyr::distinct(data)
      }
      data <- dplyr::filter(data, !is.na(data[1]) & !is.na(data[2]))
      colnames(data) <- toupper(colnames(data))
      return(data)
    } else {
      stop("Invalid file extension")
    }
  }
}
