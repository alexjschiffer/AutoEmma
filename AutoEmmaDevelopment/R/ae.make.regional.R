##################################################
#   Function Name  : Make Regional Manhattan Plot#
#   Program Author : Alex Schiffer               #
#   Last Updated   : June 19, 2018               #
#                                                #
#   D'Amato Lab, Boston Children's Hospital      #
##################################################


ae.make.regional <- function(resultsdata, filename = "filename.pdf", chromosomes = c(1:20), sl = 4.5, gwl = 7, title = "Regional Manhattan Plot"){

  ae.check.required()

  pdf(filename, width = 11, height = 8.5)
  for(i in chromosomes){
  par(mar=c(5,6,5,3)) # bottom, left, top, right
  manhattan(subset(resultsdata, CHR == i), main = title, ylim = c(0, 8), cex = 1, cex.axis = 1.2,
            cex.lab = 1.7, cex.main = 1.6, suggestiveline = sl, genomewideline = gwl,col = "royalblue3")
  }
  dev.off()
}
