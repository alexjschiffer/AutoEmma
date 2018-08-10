##################################################
#   Function Name  : Convert From atcg to binary #
#   Program Author : Alex Schiffer               #
#   Last Updated   : June 19, 2018               #
#                                                #
#   D'Amato Lab, Boston Children's Hospital      #
##################################################

ae.gwl <- function(pval, snps) {
  return(-log10(pval/snps))
}
