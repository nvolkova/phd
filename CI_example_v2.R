plot_96_profile_CI <- function (mut_matrix, mut_matrix_lower = NULL, 
                                mut_matrix_upper = NULL, CI=FALSE, 
                                colors=c("deepskyblue2","black","red3","grey","olivedrab3","pink"), ymax=NA, size=9) 
{
  if(is.na(ymax)) ymax = max(mut_matrix_upper)
  if (ymax>1) step = 10
  if (ymax<=1) step = 0.1
  C_TRIPLETS = c(
    "ACA", "ACC", "ACG", "ACT",
    "CCA", "CCC", "CCG", "CCT",
    "GCA", "GCC", "GCG", "GCT",
    "TCA", "TCC", "TCG", "TCT")
  T_TRIPLETS = c(
    "ATA", "ATC", "ATG", "ATT",
    "CTA", "CTC", "CTG", "CTT",
    "GTA", "GTC", "GTG", "GTT",
    "TTA", "TTC", "TTG", "TTT")
  types <- c(rep(C_TRIPLETS,3), rep(T_TRIPLETS,3))
  context = types
  substitution = rep(c("C>A","C>G","C>T","T>A","T>C","T>G"), each = 16)
  substring(context, 2, 2) = "*"
  df = data.frame(substitution = substitution, context = context)
  df2 = cbind(df, as.data.frame(mut_matrix))
  df3 = melt(df2, id.vars = c("substitution", "context"))
  value = NULL
  plot = ggplot(data = df3, aes(x = context, y = value, fill = substitution, width = 0.6)) + 
    geom_bar(stat = "identity", colour = "black", size = 0.2) + 
    scale_fill_manual(values = colors) + 
    facet_grid(variable ~ substitution) + 
    coord_cartesian(ylim = c(0,ymax)) + 
    scale_y_continuous(breaks = seq(0, ymax, step)) +
    guides(fill = FALSE) + theme_bw() + 
    theme(text = element_text(family='ArialMT'),
          axis.title.y = element_text(size = 14,vjust = 1), 
          axis.text.y = element_text(size = 8), 
          axis.title.x = element_text(size = 14), 
          axis.text.x = element_text(size = 8, angle = 90, vjust = 0.4), 
          strip.text.x = element_text(size = 16), 
          strip.text.y = element_text(size = size), 
          panel.grid.major.x = element_blank(), 
          strip.background = element_blank(), 
          panel.border = element_rect(colour="white"),
          panel.spacing = unit(0.1,'lines')) +
    labs(y="Number of mutations", x="Context")
  
  if (CI) { # add confidence intervals
    rownames(mut_matrix_upper) = rownames(mut_matrix_lower) = NULL
    df_lower = cbind(df, as.data.frame(mut_matrix_lower))
    df_upper = cbind(df, as.data.frame(mut_matrix_upper))
    df_lower = melt(df_lower, id.vars = c("substitution", "context"))
    df_upper = melt(df_upper, id.vars = c("substitution", "context"))
    df2 = cbind(df, as.data.frame(mut_matrix[,colnames(mut_matrix_lower)]))
    df3 = melt(df2, id.vars = c("substitution", "context"))
    df3 <- cbind(df3, value_min = df_lower$value, value_max = df_upper$value)
    plot = plot + geom_pointrange(data=df3, aes(ymin=value,ymax=value_max,colour = substitution), size=1.5, fatten = 0.001,show.legend = F) +
      scale_color_manual(values=c(colors,"white")) +
      geom_pointrange(data=df3, aes(ymin=value_min,ymax=value,colour="white"), size=1.5, fatten = 0.001,show.legend = F)
  }
  return(plot)
}

library(ggplot2)
library(reshape2)
load('Learned.sigs.RData')
mmr <- mut_mat[1:96,c('CD0134a','CD0134c','CD0134d','CD0135a','CD0135c','CD0135d')]
polepms <- mut_mat[1:96,c('CD0246d','CD0246e')]
low <- cbind(rowMeans(mmr) - apply(mmr,1,sd)/sqrt(6), rowMeans(polepms) - apply(polepms,1,sd)/sqrt(2) )
high <- cbind(rowMeans(mmr) + apply(mmr,1,sd)/sqrt(6) , rowMeans(polepms) + apply(polepms,1,sd)/sqrt(2) )
colnames(low) = colnames(high) = c('mmr','polepms')
ts <- cbind(rowMeans(mmr), rowMeans(polepms))
colnames(ts) <- c('mmr','polepms')
p <- plot_96_profile_CI(mut_matrix = ts,
                        mut_matrix_lower = low, mut_matrix_upper = high, CI=T)
ggsave(plot=p, file='~/1E.pdf', width=12, height=6)

mult <- colSums(ts)
ts <- apply(ts, 2, function (x) x/sum(x))
low[,1] <- low[,1]/mult[1]; low[,2] <- low[,2]/mult[2]
high[,1] <- high[,1]/mult[1]; high[,2] <- high[,2]/mult[2]
p <- plot_96_profile_CI(mut_matrix = ts,
                        mut_matrix_lower = low, mut_matrix_upper = high, CI=T, ymax=0.15)
ggsave(plot=p, file='~/1E_rel.pdf', width=12, height=6)



  }