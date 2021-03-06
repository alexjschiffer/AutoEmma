##################################################
#   Function Name  : Make Manhattan Plot         #
#   Program Author : Alex Schiffer               #
#   Last Updated   : September 18, 2018          #
#                                                #
#   D'Amato Lab, Boston Children's Hospital      #
##################################################

ae.manhattan <- function(df, file_name, point_size = 1.2, colors = c("royalblue3", "gray18"), title = "Manhattan Plot",
                         sl = 4.5, pval = 0.05, snps = 100000, exclude = 0.5, annotate = TRUE, wide = TRUE){

  ae.check.required()
  gwl <- -log10(pval/snps)

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
  if(annotate){
    manhattan <- ggplot(data = end, mapping = aes(x = BPcum, y = -log10(P))) +
      geom_point(mapping = aes(color = as.factor(CHR)), size = point_size, alpha = 0.9) +
      scale_color_manual(values = rep(colors, dim(axis_df)[1])) +
      scale_x_continuous(label = axis_df$CHR, breaks = axis_df$center, expand = c(0.01,0)) +
      scale_y_continuous(limits = c(exclude,gwl+0.5), expand = c(0,0), breaks = seq(0,gwl+0.5)) +
      theme(legend.position = "none") +
      ggtitle(title) +
      labs(x = "Chromosome", y = "-log10(p)") +
      geom_line(y = sl, color = "blue", alpha = 0.8) +
      geom_line(y = gwl, color = "red", alpha = 0.8) +
      geom_label_repel(data = subset(end, is_annotate == "yes"), mapping = aes(label = paste(SNP,"\np =",formatC(P,format = "e",digits = 2))), size = 3)
  } else {
    manhattan <- ggplot(data = end, mapping = aes(x = BPcum, y = -log10(P))) +
      geom_point(mapping = aes(color = as.factor(CHR)), size = point_size, alpha = 0.9) +
      scale_color_manual(values = rep(colors, dim(axis_df)[1])) +
      scale_x_continuous(label = axis_df$CHR, breaks = axis_df$center, expand = c(0.01,0)) +
      scale_y_continuous(limits = c(exclude,gwl+0.5), expand = c(0,0), breaks = seq(0,gwl+0.5)) +
      theme(legend.position = "none") +
      ggtitle(title) +
      labs(x = "Chromosome", y = "-log10(p)") +
      geom_line(y = sl, color = "blue", alpha = 0.8) +
      geom_line(y = gwl, color = "red", alpha = 0.8)
  }

  if(wide){
    ggsave(filename = file_name, plot = manhattan, width = 16, height = 9)
  }else{
    ggsave(filename = file_name, plot = manhattan, width = 9, height = 9)
  }
}
