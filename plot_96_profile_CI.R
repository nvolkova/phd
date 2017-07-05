plot_96_profile_CI <- function (types,mut_matrix, mut_matrix_lower = NULL, 
                                mut_matrix_upper = NULL, CI=FALSE, 
                                colors=c("deepskyblue2","black","red3","grey","olivedrab3","pink"), ymax = 0.2, size=9) 
{
  norm_mut_matrix = mut_matrix
  mult <- colSums(norm_mut_matrix)
  norm_mut_matrix = apply(mut_matrix, 2, function(x) x/sum(x))
  context = types
  substitution = rep(c("C>A","C>G","C>T","T>A","T>C","T>G"), each = 16)
  substring(context, 2, 2) = "."
  df = data.frame(substitution = substitution, context = context)
  rownames(norm_mut_matrix) = NULL
  df2 = cbind(df, as.data.frame(norm_mut_matrix))
  df3 = melt(df2, id.vars = c("substitution", "context"))
  value = NULL
  plot = ggplot(data = df3, aes(x = context, y = value, fill = substitution, width = 0.6)) + 
         geom_bar(stat = "identity", colour = "black", size = 0.2) + 
         scale_fill_manual(values = colors) + 
         facet_grid(variable ~ substitution) + 
         coord_cartesian(ylim = c(0,ymax)) + 
         scale_y_continuous(breaks = seq(0, ymax, 0.1)) +
         guides(fill = FALSE) + theme_bw() + theme(axis.title.y = element_text(size = 16, vjust = 1), axis.title.x = element_text(size = 16),
                                                   axis.text.y = element_text(size = size), 
                                                   axis.text.x = element_text(size = size, angle = 90, vjust = 0.4), 
                                                   strip.text.x = element_text(size = size), strip.text.y = element_text(size = size), 
                                                   panel.grid.major.x = element_blank()) +
         labs(y="Relative contribution", x="Context")
  
  if (CI) { # add confidence intervals
    mut_matrix_lower <- mut_matrix_lower[row.names(mut_matrix),]
    mut_matrix_upper <- mut_matrix_upper[row.names(mut_matrix),]
    norm_mut_matrix_lower <- mut_matrix_lower
    norm_mut_matrix_upper <- mut_matrix_upper
    for (i in 1:ncol(norm_mut_matrix_lower)) {
      norm_mut_matrix_lower[,i] = norm_mut_matrix_lower[,i] / mult[colnames(norm_mut_matrix_lower)[i]]
      norm_mut_matrix_upper[,i] = norm_mut_matrix_upper[,i] / mult[colnames(norm_mut_matrix_lower)[i]]
    }
    rownames(norm_mut_matrix_upper) = rownames(norm_mut_matrix_lower) = NULL
    df_lower = cbind(df, as.data.frame(norm_mut_matrix_lower))
    df_upper = cbind(df, as.data.frame(norm_mut_matrix_upper))
    df_lower = melt(df_lower, id.vars = c("substitution", "context"))
    df_upper = melt(df_upper, id.vars = c("substitution", "context"))
    df2 = cbind(df, as.data.frame(norm_mut_matrix[,colnames(norm_mut_matrix_lower)]))
    df3 = melt(df2, id.vars = c("substitution", "context"))
    df3 <- cbind(df3, value_min = df_lower$value, value_max = df_upper$value)
    plot = plot + geom_pointrange(data=df3, aes(ymin=value,ymax=value_max,colour = substitution), size=1.5, fatten = 0.001,show.legend = F) +
                  scale_color_manual(values=c(colors,"white")) +
                  geom_pointrange(data=df3, aes(ymin=value_min,ymax=value,colour="white"), size=1.5, fatten = 0.001,show.legend = F)
  
  }

  return(plot)
}
