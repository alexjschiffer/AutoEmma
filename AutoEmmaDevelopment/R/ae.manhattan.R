##################################################
#   Function Name  : Make Manhattan Plot         #
#   Program Author : Alex Schiffer               #
#   Last Updated   : August 10, 2018             #
#                                                #
#   D'Amato Lab, Boston Children's Hospital      #
##################################################

ae.manhattan <- function(df, file_name, point_size = 1.2, colors = c("royalblue3", "gray18"), title = "Manhattan Plot", sl = 4.5, gwl = 7, exclude = 0.5){

  ae.check.required()

  # Calculate cummulative position of each chromosome
  end <- df %>%
    group_by(CHR) %>%
    summarize(chr_len = max(BP)) %>%
    mutate(total = cumsum(chr_len) - chr_len) %>%
    select(-chr_len) %>%
    left_join(df, ., by = c("CHR" = "CHR")) %>%
    arrange(CHR, BP) %>%  # Sort by CHR, then BP
    mutate(BPcum = BP + total)
  # Determine Where to place the chromosome labels on the x-axis
  axis_df <- end %>%
    group_by(CHR) %>%
    summarize(center = (max(BPcum) + min(BPcum)) / 2)
  # Remove data points below 1
  end <- filter(end, -log10(P) > exclude)
  # Annotate data points above suggestive line
  end <- end %>%
    mutate(is_annotate = ifelse(-log10(P) > sl, "yes", "no"))
  # Create manhattan plot
  ggplot(data = end, mapping = aes(x = BPcum, y = -log10(P))) +
    geom_point(mapping = aes(color = as.factor(CHR)), size = point_size, alpha = 0.9) +
    scale_color_manual(values = rep(colors, dim(axis_df)[1])) +
    scale_x_continuous(label = axis_df$CHR, breaks = axis_df$center, expand = c(0.01,0)) +
    scale_y_continuous(limits = c(exclude,8), expand = c(0,0), breaks = seq(0,8)) +
    theme(legend.position = "none") +
    ggtitle(title) +
    labs(x = "Chromosome", y = "-log10(p)") +
    geom_line(y = sl, color = "blue", alpha = 0.8) +
    geom_line(y = gwl, color = "red", alpha = 0.8) +
    geom_label_repel(data = subset(end, is_annotate == "yes"), mapping = aes(label = SNP), size = 3) +
    ggsave(file_name)
}
