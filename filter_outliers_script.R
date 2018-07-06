##################################################
#   Program Name   : Filter Outliers             #
#   Program Author : Alex Schiffer               #
#                                                #
#   Date Created   : June 1, 2018                #
#   Last Update    : June 14, 2018               #
#                                                #
#   D'Amato Lab, Boston Children's Hospital      #
##################################################

setwd("~/Desktop")
STATS <- "new_data.csv" # Output file from Average By Group
IND <- "ind_data.csv"   # Individual phenotype data

# Select Conservativeness of Outlier Removal
## Options: 2, 2.5, 3
NUM <- 2

##################################################
##      Only change values above this box       ##
##   Data must be tab separated values (TSV)    ##
##################################################

library(dplyr)
library(readr)

# Load in Summary Statistics
data <- readr::read_tsv(STATS, col_names = TRUE)
data <- dplyr::select(data, one_of(c("Count", "Strain", "Average", "Stdev", "Mad",
                                     "Min3", "Min25", "Min2", "Max2", "Max25", "Max3")))
# Load in raw, individual data
individuals <- readr::read_tsv(IND, col_names = TRUE)
values <- individuals[["Pheno"]]

# Create variables necessary for expansion
nstrain <- c()
navg = c()
nsd = c()
nmad = c()
nmin3 = c()
nmin25 = c()
nmin2 = c()
npos2 = c()
npos25 = c()
npos3 = c()
counter <- 1

# Extract out necessary information from average by group data
count <- data[["Count"]]
groups <- data[["Strain"]]
avg <- data[["Average"]]
sd <- data[["Stdev"]]
mad <- data[["Mad"]]
min3 <- data[["Min3"]]
min25 <- data[["Min25"]]
min2 <- data[["Min2"]]
pos2 <- data[["Max2"]]
pos25 <- data[["Max25"]]
pos3 <- data[["Max3"]]

# Expand to original number of individuals
for (row in count) {
  nstrain <- append(nstrain, rep(groups[counter], count[counter]))
  navg <- append(navg, rep(avg[counter], count[counter]))
  nsd <- append(nsd, rep(sd[counter], count[counter]))
  nmad <- append(nmad, rep(mad[counter], count[counter]))
  nmin3 <- append(nmin3, rep(min3[counter], count[counter]))
  nmin25 <- append(nmin25, rep(min25[counter], count[counter]))
  nmin2 <- append(nmin2, rep(min2[counter], count[counter]))
  npos2 <- append(npos2, rep(pos2[counter], count[counter]))
  npos25 <- append(npos25, rep(pos25[counter], count[counter]))
  npos3 <- append(npos3, rep(pos3[counter], count[counter]))
  counter <- counter + 1
}

expanded_data <- dplyr::data_frame(Strain = nstrain, Pheno = values, Average = navg, Stdev = nsd, Mad = nmad,
                              Min3 = nmin3, Min25 = nmin25, Min2 = nmin2, Pos2 = npos2, Pos25 = npos25, Pos3 = npos3)

# Filter Outliers
if (NUM == 2) {
  filtered_data <- dplyr::filter(expanded_data, Mad == 0 | Pheno > Min2 & Pheno < Pos2)
} else if (NUM == 2.5) {
  filtered_data <- dplyr::filter(expanded_data, Mad == 0 | Pheno > Min25 & Pheno < Pos25)
} else if (NUM == 3) {
  filtered_data <- dplyr::filter(expanded_data, Mad == 0 | Pheno > Min3 & Pheno < Pos3)
} else {
  print ("Invalid NUM Selected.")
}

# Save Filtered Data
write_tsv(filtered_data, "filtered_data.txt", append = FALSE, col_names = TRUE)

##################################################
##                  End of File                 ##
##################################################