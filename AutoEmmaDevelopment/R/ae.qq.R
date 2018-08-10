##################################################
#   Function Name  : Make Q-Q Plot               #
#   Program Author : Alex Schiffer               #
#   Last Updated   : August 10, 2018             #
#                                                #
#   D'Amato Lab, Boston Children's Hospital      #
##################################################

ae.qq <- function(df, file_name, title = "Q-Q Plot"){

  ae.check.required()

  n <- df[["P"]]
  end_df <- df %>%
    mutate(expected = sort(-log10(ppoints(length(n))))) %>%
    mutate(observed = sort(-log10(P))) %>%
    filter(-log10(P) > 0.25)

  ggplot(end_df) + geom_point(aes(expected, observed), color = "royalblue3") +
    geom_abline(intercept = 0, slope = 1, alpha = 0.5) +
    scale_x_continuous(expand = c(0.01,0)) +
    scale_y_continuous(expand = c(0.01,0)) +
    ggtitle(title) +
    labs(x = "Expected -log10(p)", y = "Observed -log10(p)") +
    ggsave(file_name)
}
