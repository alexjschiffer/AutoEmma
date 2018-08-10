##################################################
#   Function Name  : Make Manhattan Plot         #
#   Program Author : Alex Schiffer               #
#   Last Updated   : June 19, 2018               #
#                                                #
#   D'Amato Lab, Boston Children's Hospital      #
##################################################

ae.manhattan <- function(resultsdata, filename, colors = c("royalblue3", "gray18"), sl = 4.5, gwl = 7, title = "Manhattan Plot"){

  ae.check.required()

  png(filename, width = 1280, height = 720)
  par(mar=c(5,6,5,3)) # bottom, left, top, right
  manhattan(resultsdata, main = title,
            ylim = c(0, 8), cex = 1, cex.axis = 1.2, chrlabs = c(1:19, "X"),
            cex.lab = 1.7, cex.main = 1.6, suggestiveline = sl, genomewideline = gwl,
            col = colors)
  dev.off()
}
