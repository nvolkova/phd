plot_96_profile_CI <- function (mut_matrix, mut_matrix_lower = NULL, 
                                mut_matrix_upper = NULL, CI=FALSE, 
                                colors=c("deepskyblue2","black","red3","grey","olivedrab3","pink"), ymax=NA, size=9) 
{
  if(is.na(ymax)) ymax = max(mut_matrix_upper)
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
         scale_y_continuous(breaks = seq(0, ymax, 10)) +
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

plot_104_profile_CI <- function (mut_matrix, mut_matrix_lower = NULL, 
                                 mut_matrix_upper = NULL,  
                                 colors=c("deepskyblue2","black","red3","orange","purple","grey","olivedrab3","pink"),size=6,ymax=0.25) 
{
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
  
  
  norm_mut_matrix = mut_matrix
  mult <- colSums(norm_mut_matrix)
  norm_mut_matrix = apply(mut_matrix, 2, function(x) x/sum(x))
  types <- c(rep(C_TRIPLETS,3), rep(T_TRIPLETS,3))
  df = as.data.frame(norm_mut_matrix)
  df$Type = c(rep(c("C>A","C>G","C>T","T>A","T>C","T>G"),each=16),rep("INS",4),rep("DEL",4))
  df$Base = c(types,rep(c("A","C","G","T"),2))
  df2 = melt(df,id.vars = c("Type","Base"))
  df2$Type_f = factor(df2$Type, levels=c("C>A","C>G","C>T","T>A","T>C","T>G","DEL","INS"))
  plot = ggplot(data = df2, aes(x = Base, y = value, fill = Type, width = 0.6)) +
    geom_bar(stat = "identity", colour = "black",size = 0.2) +
    scale_fill_manual(values = colors) + 
    facet_grid(variable ~ Type_f,scales = "free_x") +
    ylab("Mutation counts") + coord_cartesian(ylim = c(0,ymax)) + 
    guides(fill = FALSE) + theme_bw() +
    theme(axis.title.y = element_text(size = 16,vjust = 1), axis.text.y = element_text(size = size), axis.title.x = element_text(size = 16), 
          axis.text.x = element_text(size = size, angle = 90, vjust = 0.4), 
          strip.text.x = element_text(size = size), strip.text.y = element_text(size = size), 
          panel.grid.major.x = element_blank())
  
  mut_matrix_lower <- mut_matrix_lower[row.names(mut_matrix),]
  mut_matrix_upper <- mut_matrix_upper[row.names(mut_matrix),]
  norm_mut_matrix_lower <- mut_matrix_lower
  norm_mut_matrix_upper <- mut_matrix_upper
  for (i in 1:ncol(norm_mut_matrix_lower)) {
    norm_mut_matrix_lower[,i] = norm_mut_matrix_lower[,i] / mult[colnames(norm_mut_matrix_lower)[i]]
    norm_mut_matrix_upper[,i] = norm_mut_matrix_upper[,i] / mult[colnames(norm_mut_matrix_lower)[i]]
  }
  rownames(norm_mut_matrix_upper) = rownames(norm_mut_matrix_lower) = NULL
  df_lower = as.data.frame(norm_mut_matrix_lower)
  df_upper = as.data.frame(norm_mut_matrix_upper)
  df_lower$Type = c(rep(c("C>A","C>G","C>T","T>A","T>C","T>G"),each=16),rep("INS",4),rep("DEL",4))
  df_lower$Base = c(types,rep(c("A","C","G","T"),2))
  df_upper$Type = c(rep(c("C>A","C>G","C>T","T>A","T>C","T>G"),each=16),rep("INS",4),rep("DEL",4))
  df_upper$Base = c(types,rep(c("A","C","G","T"),2))
  df_lower = melt(df_lower, id.vars = c("Type","Base"))
  df_upper = melt(df_upper, id.vars = c("Type","Base"))
  df3 <- cbind(df2, value_min = df_lower$value, value_max = df_upper$value)
  plot2 = plot + geom_pointrange(data=df3, aes(ymin=value,ymax=value_max,colour = Type_f), fatten = 0.01, size=1, show.legend = F) +
    scale_color_manual(values=c(colors,"white")) +
    geom_pointrange(data=df3, aes(ymin=value_min,ymax=value,colour="white"), size=1, fatten = 0.01,show.legend = F)
  
  return(plot2)
}
