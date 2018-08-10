##################################################
#   Function Name  : Check Required Packages     #
#   Program Author : Alex Schiffer               #
#   Last Updated   : June 28, 2018               #
#                                                #
#   D'Amato Lab, Boston Children's Hospital      #
##################################################

ae.check.required <- function() {
  if(!suppressPackageStartupMessages(library(cowplot, quietly = TRUE, logical.return = TRUE, warn.conflicts = FALSE))) {
    suppressPackageStartupMessages(install.packages("cowplot", dependencies = TRUE))
  }
    if(!suppressPackageStartupMessages(library(stringi, quietly = TRUE, logical.return = TRUE, warn.conflicts = FALSE))) {
      suppressPackageStartupMessages(install.packages("stringi", dependencies = TRUE))
    }
    if(!suppressPackageStartupMessages(library(readr, quietly = TRUE, logical.return = TRUE, warn.conflicts = FALSE))) {
      suppressPackageStartupMessages(install.packages("readr", dependencies = TRUE))
    }
    if(!suppressPackageStartupMessages(library(dplyr, quietly = TRUE, logical.return = TRUE, warn.conflicts = FALSE))) {
      suppressPackageStartupMessages(install.packages("dplyr", dependencies = TRUE))
    }
    if(!suppressPackageStartupMessages(library(qqman, quietly = TRUE, logical.return = TRUE, warn.conflicts = FALSE))) {
      suppressPackageStartupMessages(install.packages("qqman", dependencies = TRUE))
    }
    if(!suppressPackageStartupMessages(library(ape, quietly = TRUE, logical.return = TRUE, warn.conflicts = FALSE))) {
      suppressPackageStartupMessages(install.packages("ape", dependencies = TRUE))
    }
    if(!suppressPackageStartupMessages(library(colorRamps, quietly = TRUE, logical.return = TRUE, warn.conflicts = FALSE))) {
      install.packages("colorRamps", dependencies = TRUE)
    }
    if(!suppressPackageStartupMessages(library(ggplot2, quietly = TRUE, logical.return = TRUE, warn.conflicts = FALSE))) {
      install.packages("ggplot2", dependencies = TRUE)
    }
    if(!suppressPackageStartupMessages(library(ggcorrplot, quietly = TRUE, logical.return = TRUE, warn.conflicts = FALSE))) {
      install.packages("ggcorrplot", dependencies = TRUE)
    }
    if(!suppressPackageStartupMessages(library(ggtree, quietly = TRUE, logical.return = TRUE, warn.conflicts = FALSE))) {
      source("https://bioconductor.org/biocLite.R")
      suppressPackageStartupMessages(biocLite("ggtree", suppressUpdates = TRUE))
    }
    if(!suppressPackageStartupMessages(library(treeio, quietly = TRUE, logical.return = TRUE, warn.conflicts = FALSE))) {
      source("https://bioconductor.org/biocLite.R")
      suppressPackageStartupMessages(biocLite("treeio", suppressUpdates = TRUE))
    }
}
