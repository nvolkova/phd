plot_sig_wb <- function (mut_matrix, colors=NA, ymin=NA, ymax = NA, flip=F, ind=F, font=11, size=12, norm=T,
                         CI = F, low = NULL, high = NULL) # plotting 96-signatures / profiles with no lines with/out CIs
{
  
  if (is.na(colors)) {
    library(RColorBrewer)
    colors <- c("#2EBAED", "#000000", "#DE1C14","#D4D2D2", "#ADCC54", "#F0D0CE")
  }
  
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
  
  IND_TRIPLETS = c(
    "A*A", "A*C", "A*G", "A*T",
    "C*A", "C*C", "C*G", "C*T",
    "G*A", "G*C", "G*G", "G*T",
    "T*A", "T*C", "T*G", "T*T")
  
  ylab = 'Number of mutations'
  if (norm) {
    mult <- colSums(mut_matrix)
    mut_matrix <- apply(mut_matrix, 2, function(x) x/sum(x))
    ylab = 'Relative contribution'
  }
  if (is.na(ymax)) ymax = max(mut_matrix)
  if (is.na(ymin)) ymin = 0
  
  TRIPLETS_96 = c(rep(C_TRIPLETS, 3), rep(T_TRIPLETS, 3))
  TRIPLETS_112 = c(rep(C_TRIPLETS, 3), rep(T_TRIPLETS, 3),IND_TRIPLETS)
  if (norm) norm_mut_matrix = apply(mut_matrix, 2, function(x) x/sum(x)) else norm_mut_matrix <- mut_matrix
  if (missing(colors)) {
    if (!ind)  colors = c("#2EBAED", "#000000", "#DE1C14",
                          "#D4D2D2", "#ADCC54", "#F0D0CE")
    else colors = c("#2EBAED", "#000000", "#DE1C14", "yellow",
                    "#D4D2D2", "#ADCC54", "#F0D0CE")
  }
  context = TRIPLETS_96
  substitution = rep(c("C>A","C>G","C>T","T>A","T>C","T>G"), each = 16)
  if (ind) {
    context = TRIPLETS_112
    substitution = rep(c("C>A","C>G","C>T","T>A","T>C","T>G","INDEL"), each = 16)
  }
  substring(context, 2, 2) = "*"
  df = data.frame(substitution = substitution, context = context)
  rownames(norm_mut_matrix) = NULL
  if (flip==T) norm_mut_matrix = -norm_mut_matrix
  df2 = cbind(df, as.data.frame(norm_mut_matrix))
  df3 = melt(df2, id.vars = c("substitution", "context"))
  value = NULL
  if (ymax>0) breaks = c(0,0.1) else breaks = c(-0.1,0)
  if (ymax>0) lbls = c('0','0.1') else lbls = c('0.1','0')
  plot = ggplot(data = df3, aes(x = context, y = value, fill = substitution, width = 0.6)) +
    geom_bar(stat = "identity", colour = "black",size = NA) +
    scale_fill_manual(values = colors) + facet_grid(variable ~ substitution) +
    ylab("Relative contribution") + coord_cartesian(ylim = c(ymin,ymax)) + 
    guides(fill = FALSE) + theme_bw() + scale_y_continuous(breaks = breaks, labels = lbls) +
    theme(axis.title.y = element_text(size = 14,vjust = 1), 
          axis.text.y = element_text(size = 8), 
          axis.title.x = element_blank(), 
          axis.text.x = element_text(size = 8, angle = 90, vjust = 0.4), 
          strip.text.x = element_text(size = 16), 
          strip.text.y = element_text(size = size), 
          panel.grid.major.x = element_blank(), 
          strip.background = element_blank(), 
          panel.border = element_rect(colour="white"),
          panel.spacing = unit(0.1,'lines'))
  
  if (CI) {
    norm_mut_matrix_lower <- as.matrix(low)
    norm_mut_matrix_upper <- as.matrix(high)
    if (norm)
      for (i in 1:ncol(norm_mut_matrix_lower)) {
        norm_mut_matrix_lower[,i] = norm_mut_matrix_lower[,i] / mult[colnames(norm_mut_matrix_lower)[i]]
        norm_mut_matrix_upper[,i] = norm_mut_matrix_upper[,i] / mult[colnames(norm_mut_matrix_upper)[i]]
      }
    rownames(norm_mut_matrix_upper) = rownames(norm_mut_matrix_lower) = NULL
    df_lower = cbind(df,as.data.frame(norm_mut_matrix_lower))
    df_upper = cbind(df,as.data.frame(norm_mut_matrix_upper))
    df_lower = melt(df_lower, id.vars = c("substitution", "context"))
    df_upper = melt(df_upper, id.vars = c("substitution", "context"))
    df3$value_min = df_lower$value
    df3$value_max = df_upper$value
    plot <- plot + geom_errorbar(data=df3, aes(ymin = value_min, ymax = value_max), colour = 'darkgrey', width = 0.05)
  }
  
  return(plot)
}

plot_bigsig_wb <- function(mut_matrix, colors=NA, norm=T, ymax=NA) # plotting 112-long signatures / profiles (Indels unresolved)
{
  
  if (is.na(colors)) {
    library(RColorBrewer)
    colors <- c(brewer.pal(7, "Set1"), brewer.pal(8,'Dark2'), brewer.pal(7,"Pastel1"))
  }
  
  TRIPLETS = c(
    "A*A", "A*C", "A*G", "A*T",
    "C*A", "C*C", "C*G", "C*T",
    "G*A", "G*C", "G*G", "G*T",
    "T*A", "T*C", "T*G", "T*T")

  DINUCL = 'DNV'
  
  INDELS = c(
    'D:1-4', "D:5-49","D:50-400",
    'DI:1-49','DI:50-400',
    'I:1-4', 'I:5-49','D:50-400')
  
  SVS = c("TD","DEL","INV","COMPLEX","TRSL","INTCHR","FOLDBACK")
  
  if (norm) norm_mut_matrix <- apply(mut_matrix, 2, function(x) x/sum(x)) else norm_mut_matrix <- mut_matrix
  
  if (is.na(ymax)) ymax <- max(norm_mut_matrix)
  
  TRIPLETS_96 = rep(TRIPLETS, 6)
  context = c(TRIPLETS_96,DINUCL,INDELS,SVS)
  substitution = c(rep(c("C>A","C>G","C>T","T>A","T>C","T>G"), each = 16),DINUCL,rep('INDEL',8),rep('SV',7))
  df = data.frame(substitution = substitution, context = context)
  df2 = cbind(df, as.data.frame(norm_mut_matrix))
  df3 = melt(df2, id.vars = c("substitution", "context"))
  df3$substitution = factor(df3$substitution, levels=c("C>A","C>G","C>T","T>A","T>C","T>G","DNV","INDEL","SV"))
  plot = ggplot(data = df3, aes(x = context, y = value, fill = substitution, width = 0.6)) +
    geom_bar(stat = "identity", colour = "black",size = 0.2) +
    scale_fill_manual(values = colors) + facet_grid(variable ~ substitution, scales = 'free_x') +
    ylab("Number of mutations") + coord_cartesian(ylim = c(0,ymax)) + 
    guides(fill = FALSE) + theme_bw() +
    theme(text = element_text(family='ArialMT'),
          axis.title.y = element_text(size = 12,vjust = 1), 
          axis.text.y = element_text(size = 8), 
          axis.title.x = element_blank(), 
          axis.text.x = element_text(size = 8, angle = 90, vjust = 0.4), 
          strip.text.x = element_text(size = 12), 
          strip.text.y = element_text(size = 12), 
          panel.grid.major.x = element_blank(), 
          strip.background = element_blank(), 
          panel.border = element_rect(colour="white"),
          panel.spacing = unit(0.1,'lines'))
  return(plot)
  
}

plot_subindel_wb <- function(mut_matrix, colors=NA, norm=T, ymax = NA, ymin = NA, CI = F, low = NULL, high = NULL, DNV = FALSE) # plotting 104 signatures / profiles with/without CI
{
  
  if (is.na(colors)) {
    library(RColorBrewer)
    colors <- c("#2EBAED", "#000000", "#DE1C14","#D4D2D2", "#ADCC54", "#F0D0CE","orange",'darkmagenta',"purple",'brown')
  }
  
  TRIPLETS = c(
    "A*A", "A*C", "A*G", "A*T",
    "C*A", "C*C", "C*G", "C*T",
    "G*A", "G*C", "G*G", "G*T",
    "T*A", "T*C", "T*G", "T*T")
  
  INDELS = c(
    '1-4 bp', "5-49 bp","50-400 bp",
    '1-49 bp','50-400 bp',
    '1-4 bp', '5-49 bp','50-400 bp')
  
  ylab = 'Number of mutations'
  if (norm) {
    mult <- colSums(mut_matrix)
    mut_matrix <- apply(mut_matrix, 2, function(x) x/sum(x))
    ylab = 'Relative contribution'
  }
  
  if (is.na(ymax)) ymax = max(mut_matrix)
  if (is.na(ymin)) ymin = 0
  
  TRIPLETS_96 = rep(TRIPLETS, 6)
  context = c(TRIPLETS_96,INDELS)
  if (DNV) context <- c(context,'DNV')
  substitution = c(rep(c("C>A","C>G","C>T","T>A","T>C","T>G"), each = 16),
                   rep('Deletions',3),rep('Complex\n indels',2),rep('Insertions',3))
  if (DNV) substitution <- c(substitution,'DNV')
  df = data.frame(substitution = substitution, context = context)
  df2 = cbind(df, as.data.frame(mut_matrix))
  df3 = melt(df2, id.vars = c("substitution", "context"))
  if (!DNV) df3$substitution = factor(df3$substitution, levels=c("C>A","C>G","C>T","T>A","T>C","T>G",
                                                       "Deletions","Complex\n indels","Insertions"))
  if (DNV) df3$substitution = factor(df3$substitution, levels=c("C>A","C>G","C>T","T>A","T>C","T>G",
                                                               "Deletions","Complex\n indels","Insertions",'DNV'))
  df3$width <- 0.6
  df3$width[df3$substitution=='Deletions' | df3$substitution=='Insertions'] <- 0.3;
  if (DNV) df3$width[df3$substitution == 'DNV'] <- 0.1
  df3$width[df3$substitution=='Complex\n indels'] <- 0.2
  plot = ggplot(data = df3, aes(x = context, y = value, fill = substitution, width = width)) +
    geom_bar(stat = "identity", colour = "black", size = NA) +
    scale_fill_manual(values = colors) + facet_grid(variable ~ substitution, scales = 'free_x') +
    ylab(ylab) + coord_cartesian(ylim = c(ymin,ymax)) + 
    guides(fill = FALSE) + theme_bw() +
    theme(axis.title.y = element_text(size = 14,vjust = 1), 
          axis.text.y = element_text(size = 8), 
          axis.title.x = element_blank(), 
          axis.text.x = element_text(size = 8, angle = 90, vjust = 0.4), 
          strip.text.x = element_text(size = 14), 
          strip.text.y = element_text(size = 14), 
          panel.grid.major.x = element_blank(), 
          strip.background = element_blank(), 
          panel.border = element_rect(colour="white"),
          panel.spacing = unit(0.1,'lines'))
  
  if (CI) {
    norm_mut_matrix_lower <- as.matrix(low)
    norm_mut_matrix_upper <- as.matrix(high)
    if (norm)
      for (i in 1:ncol(norm_mut_matrix_lower)) {
        norm_mut_matrix_lower[,i] = norm_mut_matrix_lower[,i] / mult[colnames(norm_mut_matrix_lower)[i]]
        norm_mut_matrix_upper[,i] = norm_mut_matrix_upper[,i] / mult[colnames(norm_mut_matrix_upper)[i]]
      }
    rownames(norm_mut_matrix_upper) = rownames(norm_mut_matrix_lower) = NULL
    df_lower = cbind(df,as.data.frame(norm_mut_matrix_lower))
    df_upper = cbind(df,as.data.frame(norm_mut_matrix_upper))
    df_lower = melt(df_lower, id.vars = c("substitution", "context"))
    df_upper = melt(df_upper, id.vars = c("substitution", "context"))
    df3$value_min = df_lower$value
    df3$value_max = df_upper$value
    plot <- plot + geom_errorbar(data=df3, aes(ymin = value_min, ymax = value_max), colour = 'darkgrey', width = 0.05)
  }
  
  return(plot)
  
}


############################################################################################################

# visualisation for the huge profile! (indels and DNV resolved)
  
plot_hugesig_wb <- function(mut_matrix, colors=NA, norm=T, ymax=NA, ytitle = 'Number of mutations')
  {
    types.full <- paste(rep(rep(c('A','C','G','T'), each = 4), 6), '[',
                        rep(c('C','T'), each = 48), '>', rep(c('A','G','T','A','C','G'), each=16), ']',
                        rep(c('A','C','G','T'), 24), sep='')
    
    DNV <- c('AC','AT','CC','CG','CT','GC','TA','TC','TG','TT')
    
    dnv.types.full <- c(paste0(rep('AC',9),'>',unlist(lapply(c('C','G','T'), function(x) paste0(rep(x,3),c('A','G','T'))))),
                        paste0(rep('AT',6),'>',c('CA','CC','CG','GA','GC','TA')),
                        paste0(rep('CC',9),'>',unlist(lapply(c('A','G','T'), function(x) paste0(rep(x,3),c('A','G','T'))))),
                        paste0(rep('CG',6),'>',c('AT','GC','GT','TA','TC','TT')),
                        paste0(rep('CT',9),'>',unlist(lapply(c('A','G','T'), function(x) paste0(rep(x,3),c('A','C','G'))))),
                        paste0(rep('GC',6),'>',c('AA','AG','AT','CA','CG','TA')),
                        paste0(rep('TA',6),'>',c('AT','CG','CT','GC','GG','GT')),
                        paste0(rep('TC',9),'>',unlist(lapply(c('A','C','G'), function(x) paste0(rep(x,3),c('A','G','T'))))),
                        paste0(rep('TG',9),'>',unlist(lapply(c('A','C','G'), function(x) paste0(rep(x,3),c('A','C','T'))))),
                        paste0(rep('TT',9),'>',unlist(lapply(c('A','C','G'), function(x) paste0(rep(x,3),c('A','C','G'))))))
    
    ID <- c('D','DI','I')
    
    if (is.na(colors)) {
      library(RColorBrewer)
      colors <- c(brewer.pal(6, "Set1"), brewer.pal(10,'Set3'), 'orange', brewer.pal(3,"Set2"), 'darkmagenta')
    }
    
    TRIPLETS = c(
      "A*A", "A*C", "A*G", "A*T",
      "C*A", "C*C", "C*G", "C*T",
      "G*A", "G*C", "G*G", "G*T",
      "T*A", "T*C", "T*G", "T*T")
    
    INDELS = c(
      'D:1:rep', 'D:1:nonrep', 'D:2-5:rep','D:2-5:nonrep',
      "D:6-50","D:51-400",
      'DI:1-49','DI:50-400',
      'I:1:rep','I:1:nonrep','I:2-5:rep','I:2-5:nonrep', 
      'I:6-50','I:51-400')
    
    SVS = c("TD","DEL","INV","COMPLEX","TRSL","INTCHR","FOLDBACK")
    
    if (norm) {
      norm_mut_matrix <- apply(mut_matrix, 2, function(x) x/sum(x))
      ytitle <- 'Relative contribution'
    } else norm_mut_matrix <- mut_matrix
    
    
    if (is.na(ymax)) ymax <- max(norm_mut_matrix)
    
    TRIPLETS_96 = rep(TRIPLETS, 6)
    context = c(TRIPLETS_96,dnv.types.full,'MNV',INDELS,SVS)
    substitution = c(rep(c("C>A","C>G","C>T","T>A","T>C","T>G"), each = 16),
                     rep(DNV, c(9,6,9,6,9,6,6,9,9,9)),
                     'MNV',
                     rep(ID,c(6,2,6)),
                     rep('SV',7))
    df = data.frame(substitution = substitution, context = context)
    df2 = cbind(df, as.data.frame(norm_mut_matrix))
    df3 = melt(df2, id.vars = c("substitution", "context"))
    df3$substitution = factor(df3$substitution, levels=c("C>A","C>G","C>T","T>A","T>C","T>G",DNV,'MNV',ID,"SV"))
    plot = ggplot(data = df3, aes(x = context, y = value, fill = substitution, width = 0.6)) +
      geom_bar(stat = "identity", colour = "black",size = 0.2) +
      scale_fill_manual(values = colors) + facet_grid(variable ~ substitution, scales = 'free_x') +
      ylab(ytitle) + coord_cartesian(ylim = c(0,ymax)) + 
      guides(fill = FALSE) + theme_bw() +
      theme(text = element_text(family='ArialMT'),
            axis.title.y = element_text(size = 12,vjust = 1), 
            axis.text.y = element_text(size = 8), 
            axis.title.x = element_blank(), 
            axis.text.x = element_text(size = 6, angle = 90, vjust = 0.4), 
            strip.text.x = element_text(size = 12), 
            strip.text.y = element_text(size = 12), 
            panel.grid.major.x = element_blank(), 
            strip.background = element_blank(), 
            panel.border = element_rect(colour="white"),
            panel.spacing = unit(0.1,'lines'))
    return(plot)
    
  }

plot_dnv_only <- function(mut_matrix, colors=NA, norm=T, ymax=NA, ytitle = 'Number of mutations') {

  DNV <- c('AC','AT','CC','CG','CT','GC','TA','TC','TG','TT')
  
  dnv.types.full <- c(paste0(rep('AC',9),'>',unlist(lapply(c('C','G','T'), function(x) paste0(rep(x,3),c('A','G','T'))))),
                      paste0(rep('AT',6),'>',c('CA','CC','CG','GA','GC','TA')),
                      paste0(rep('CC',9),'>',unlist(lapply(c('A','G','T'), function(x) paste0(rep(x,3),c('A','G','T'))))),
                      paste0(rep('CG',6),'>',c('AT','GC','GT','TA','TC','TT')),
                      paste0(rep('CT',9),'>',unlist(lapply(c('A','G','T'), function(x) paste0(rep(x,3),c('A','C','G'))))),
                      paste0(rep('GC',6),'>',c('AA','AG','AT','CA','CG','TA')),
                      paste0(rep('TA',6),'>',c('AT','CG','CT','GC','GG','GT')),
                      paste0(rep('TC',9),'>',unlist(lapply(c('A','C','G'), function(x) paste0(rep(x,3),c('A','G','T'))))),
                      paste0(rep('TG',9),'>',unlist(lapply(c('A','C','G'), function(x) paste0(rep(x,3),c('A','C','T'))))),
                      paste0(rep('TT',9),'>',unlist(lapply(c('A','C','G'), function(x) paste0(rep(x,3),c('A','C','G'))))))
  
  if (is.na(colors)) {
    library(RColorBrewer)
    colors <- c( brewer.pal(10,'Set3'), 'orange')
  }

  if (norm) {
    norm_mut_matrix <- apply(mut_matrix, 2, function(x) x/sum(x))
    ytitle <- 'Relative contribution'
  } else norm_mut_matrix <- mut_matrix
  
  if (is.na(ymax)) ymax <- max(norm_mut_matrix)
  
  context = c(dnv.types.full,'MNV')
  substitution = c(rep(DNV, c(9,6,9,6,9,6,6,9,9,9)),
                   'MNV')
  df = data.frame(substitution = substitution, context = context)
  df2 = cbind(df, as.data.frame(norm_mut_matrix))
  df3 = melt(df2, id.vars = c("substitution", "context"))
  df3$substitution = factor(df3$substitution, levels=c(DNV,'MNV'))
  plot = ggplot(data = df3, aes(x = context, y = value, fill = substitution, width = 0.6)) +
    geom_bar(stat = "identity", colour = "black",size = 0.2) +
    scale_fill_manual(values = colors) + facet_grid(variable ~ substitution, scales = 'free_x') +
    ylab(ytitle) + coord_cartesian(ylim = c(0,ymax)) + 
    guides(fill = FALSE) + theme_bw() +
    theme(text = element_text(family='ArialMT'),
          axis.title.y = element_text(size = 12,vjust = 1), 
          axis.text.y = element_text(size = 8), 
          axis.title.x = element_blank(), 
          axis.text.x = element_text(size = 6, angle = 90, vjust = 0.4), 
          strip.text.x = element_text(size = 12), 
          strip.text.y = element_text(size = 12), 
          panel.grid.major.x = element_blank(), 
          strip.background = element_blank(), 
          panel.border = element_rect(colour="white"),
          panel.spacing = unit(0.1,'lines'))
  return(plot)
  
}

TRIPLETS = c(
  "A*A", "A*C", "A*G", "A*T",
  "C*A", "C*C", "C*G", "C*T",
  "G*A", "G*C", "G*G", "G*T",
  "T*A", "T*C", "T*G", "T*T")
SVS = c("Tandem dup.","Deletions","Inversions","Complex events","Translocations","Interchr. events","Foldbacks")
INDELS = c(
  '1 bp (in repeats)', '1 bp (non-rep.)', '2-5 bp (in repeats)','2-5 bp (non-rep.)',
  "6-50 bp","over 50 bp",
  '1-49 bp','50-400 bp',
  '1 bp (in repeats)', '1 bp (non-rep.)', '2-5 bp (in repeats)','2-5 bp (non-rep.)',
  '6-50 bp','over 50 bp')
library(grid)
library(ggplot2)
library(reshape2)

plot_dnvhugesig_wb <- function(mut_matrix, colors=NA, norm=T, ymax=NA, ytitle = 'Number of mutations', # signatures for 119 profile, DNV collapsed, indel resolved
                               CI = F, low = NA, high = NA, diff_scale = F)
{
  types.full <- paste(rep(rep(c('A','C','G','T'), each = 4), 6), '[',
                      rep(c('C','T'), each = 48), '>', rep(c('A','G','T','A','C','G'), each=16), ']',
                      rep(c('A','C','G','T'), 24), sep='')
  
  ID <- c('Deletions','Complex\n indels','Insertions')
  
  if (is.na(colors)) {
    library(RColorBrewer)
    colors <- c(brewer.pal(6, "Set1"), 'brown', brewer.pal(3,"Set2"), 'darkmagenta')
  }
  
  TRIPLETS = c(
    "A*A", "A*C", "A*G", "A*T",
    "C*A", "C*C", "C*G", "C*T",
    "G*A", "G*C", "G*G", "G*T",
    "T*A", "T*C", "T*G", "T*T")
  
  INDELS = c(
    '1 bp (in repeats)', '1 bp (non-rep.)', '2-5 bp (in repeats)','2-5 bp (non-rep.)',
    "6-50 bp","over 50 bp",
    '1-49 bp','50-400 bp',
    '1 bp (in repeats)', '1 bp (non-rep.)', '2-5 bp (in repeats)','2-5 bp (non-rep.)',
    '6-50 bp','over 50 bp')
  
  SVS = c("Tandem dup.","Deletions","Inversions","Complex events","Translocations","Interchr. events","Foldbacks")
  
  if (norm) {
    mult <- colSums(mut_matrix)
    norm_mut_matrix <- apply(mut_matrix, 2, function(x) x/sum(x))
    ytitle <- 'Relative contribution'
  } else norm_mut_matrix <- mut_matrix
  
  
  if (is.na(ymax)) ymax <- max(norm_mut_matrix)
  
  TRIPLETS_96 = rep(TRIPLETS, 6)
  context = c(TRIPLETS_96,'2 bp','over 2 bp',INDELS,SVS)
  substitution = c(rep(c("C>A","C>G","C>T","T>A","T>C","T>G"), each = 16),
                   'MNV',
                   'MNV',
                   rep(ID,c(6,2,6)),
                   rep('Structural\n Variants',7))
  df = data.frame(substitution = substitution, context = context)
  df2 = cbind(df, as.data.frame(norm_mut_matrix))
  df3 = melt(df2, id.vars = c("substitution", "context"))
  df3$substitution = factor(df3$substitution, levels=c("C>A","C>G","C>T","T>A","T>C","T>G",'MNV',ID,"Structural\n Variants"))
  df3$width = rep(0.7,nrow(df3));
  df3$width[df3$substitution == 'MNV' | df3$substitution == 'Complex\n indels'] <- 0.2
  if (diff_scale)
    scales <- 'free'
  else scales <- 'free_x'
  p = ggplot(data = df3, aes(x = context, y = value, fill = substitution, width = width)) +
    geom_bar(stat = "identity", colour = "black", size = NA) +
    scale_fill_manual(values = colors) + facet_grid(variable ~ substitution, scales = scales) +
    ylab(ytitle) + 
    guides(fill = FALSE) + theme_bw() +
    theme(axis.title.y = element_text(size = 10,vjust = 1), 
          axis.text.y = element_text(size = 10), 
          axis.title.x = element_blank(), 
          axis.text.x = element_text(size = 10, angle = 90, vjust = 0.4), 
          strip.text.x = element_text(size = 10), 
          strip.text.y = element_text(size = 10),
          panel.grid.major.x = element_blank(), 
          strip.background = element_blank(), 
          panel.border = element_rect(colour="white"),
          panel.spacing = unit(0.1,'lines'))
  if (!diff_scale) p <- p + coord_cartesian(ylim = c(0,ymax)) 
  
  
  if (CI) {
    norm_mut_matrix_lower <- as.matrix(low)
    norm_mut_matrix_upper <- as.matrix(high)
    if (norm)
      for (i in 1:ncol(norm_mut_matrix_lower)) {
        norm_mut_matrix_lower[,i] = norm_mut_matrix_lower[,i] / mult[colnames(norm_mut_matrix_lower)[i]]
        norm_mut_matrix_upper[,i] = norm_mut_matrix_upper[,i] / mult[colnames(norm_mut_matrix_upper)[i]]
      }
    rownames(norm_mut_matrix_upper) = rownames(norm_mut_matrix_lower) = NULL
    df_lower = cbind(df,as.data.frame(norm_mut_matrix_lower))
    df_upper = cbind(df,as.data.frame(norm_mut_matrix_upper))
    df_lower = melt(df_lower, id.vars = c("substitution", "context"))
    df_upper = melt(df_upper, id.vars = c("substitution", "context"))
    df3$value_min = df_lower$value
    df3$value_max = df_upper$value
    p <- p + geom_errorbar(data=df3, aes(ymin = value_min, ymax = value_max), colour = 'darkgrey', width = 0.02)
#    plot <- plot + geom_pointrange(data=df3, aes(ymin=value, ymax=value_max, colour=substitution),
#                                  fatten = 0.001, size=0.5, show.legend = F) +
#      scale_color_manual(values=c(colors,"white")) +
#      geom_pointrange(data=df3, aes(ymin=value_min,ymax=value),
#                      size=0.5, fatten = 0.001,show.legend = F, colour="white")
  }

  return(p)
  
}

intersperse_119 <- function(y) {
  return(c(y[1:16],NA,
           y[17:32],NA,
           y[33:48],NA,
           y[49:64],NA,
           y[65:80],NA,
           y[81:96],NA,
           y[97],NA,
           y[98],NA,
           y[99:104],NA,
           y[105],NA,
           y[106],NA,
           y[107:112],NA,
           y[113:119]))
}

plot_huge_barplot <- function(y, at=c(-1,0,1), plot_main = '', log = T,CI=F,low,high,
                              col = c("#2EBAED","#000000","#DE1C14","#D4D2D2","#ADCC54","#F0D0CE","brown","#8DD3C7","goldenrod","#BEBADA","darkmagenta")) {
  cl <- rep(col,c(rep(16,6),2,6,2,6,7))
  cl1 <- intersperse_119(cl)
  fff <- data.frame(bar = intersperse_119(y))
  if (log) fff$bar <- fff$bar / log(10)
  #                  low = Wsig_all_low[match(y,colnames(W_all)),],
  #                  high = Wsig_all_up[match(y,colnames(W_all)),])
  foo <- barplot(fff$bar, col = cl1,
                 width = c(rep(1,101),3,6,3,6,4,rep(3,6),4,6,3,6,4,rep(3,6),2,rep(2.5,7)),
                 border=NA, main=NA,
                 ylim = c(min(c(at,min(fff$bar)),na.rm=T),max(c(at,max(fff$bar)),na.rm=T)),
                 yaxt = 'n')
  title(main=plot_main)
  if (log) axis(side = 2, at = at, labels = 10**at, las = 2) else axis(side = 2, at = at, labels = at, las = 2)
  if (CI) {
    fff$low <- intersperse_119(low)
    fff$high <- intersperse_119(high)
    if (log) {
      fff$low <- fff$low / log(10)
      fff$high <- fff$high / log(10)
    }
    arrows(x0 = foo, y0 = fff$low, y1 = fff$high, code=3, length=0.01, lwd=1, angle=90, col='darkgrey')
  }
  return(foo)
}

lighten <- function(color, factor=1.4){
  col <- col2rgb(color)
  col <- col*factor
  col[col>255] <- 255
  col <- rgb(t(col), maxColorValue=255)
  col
}

plot_huge_line <- function(y, at=c(-1,0,1), plot_main = '', log = T,CI=F,low,high, cex = 3, lwd = 4, labels=NA,
                              col = c("#2EBAED","#000000","#DE1C14","#D4D2D2","#ADCC54","#F0D0CE","brown","#8DD3C7","goldenrod","#BEBADA","darkmagenta")) {
  fff <- data.frame(bar = intersperse_119(y))
  if (log) fff$bar <- fff$bar / log(10)
  #                  low = Wsig_all_low[match(y,colnames(W_all)),],
  #                  high = Wsig_all_up[match(y,colnames(W_all)),])
  line_X_axis <- cumsum(c(rep(1,101),2,4,3,4,2,rep(3,6),2,4,3,4,2,rep(3,6),2,rep(2.5,7)))
  plot(x = line_X_axis, y = rep(NA,length(line_X_axis)),
       #col = col[1], type = 'l', lwd = lwd,
       xlim = c(0,max(line_X_axis)),
       ylim = c(min(c(at,min(fff$bar,na.rm=T)),na.rm=T),max(c(at,max(fff$bar,na.rm=T)),na.rm=T)),
       yaxt = 'n', xaxt = 'n', xlab = '', ylab = '', bty='n')

  fff[fff < -1] <- -1
  
  if (CI) {
    fff$low <- intersperse_119(low)
    fff$high <- intersperse_119(high)
    if (log) {
        fff$low <- fff$low / log(10)
        fff$high <- fff$high / log(10)
    }
    
    fff[fff < -1] <- -1
    
    col1 <- sapply(col, function(x) {
      tmpf <- colorRampPalette(c(x, "white"))
      return(tmpf(5)[3])
    })
    
    #pos = 1
    #for (j in 1:6) {
    #  arrows(x0 = line_X_axis[pos:(pos+15)], y0 = fff$low[pos:(pos+15)], y1 = fff$high[pos:(pos+15)], col = col1[j], code=3, length=0.05, lwd=1, angle=90)
    #  pos = pos + 17
    #}
    #polygon(x = c(line_X_axis[1:16],line_X_axis[16:1]), y = c(fff$high[1:16],fff$low[16:1]), col = col1[1], border = NA)
    #arrows(x0 = line_X_axis[c(103,105)], y0 = fff$low[c(103,105)], y1 = fff$high[c(103,105)], col = col1[7], code=3, length=0.05, lwd=1, angle=90)
    #arrows(x0 = line_X_axis[107:112], y0 = fff$low[107:112], y1 = fff$high[107:112], col = col1[8], code=3, length=0.05, lwd=1, angle=90)
    #arrows(x0 = line_X_axis[c(114,116)], y0 = fff$low[c(114,116)], y1 = fff$high[c(114,116)], col = col1[9], code=3, length=0.05, lwd=1, angle=90)
    #arrows(x0 = line_X_axis[118:123], y0 = fff$low[118:123], y1 = fff$high[118:123], col = col1[10], code=3, length=0.05, lwd=1, angle=90)
    #arrows(x0 = line_X_axis[125:131], y0 = fff$low[125:131], y1 = fff$high[125:131], col = col1[11], code=3, length=0.05, lwd=1, angle=90)
    
    pos = 1
    for (j in 1:6) {
      for ( k in pos:(pos+15)) {
        lines(x = c(line_X_axis[k] - 0.5,line_X_axis[k] + 0.5,line_X_axis[k] + 0.5),
              y = c(rep(fff$low[k],2),fff$low[k+1]), col = col1[j], lwd = lwd)
        lines(x = c(line_X_axis[k] - 0.5,line_X_axis[k] + 0.5,line_X_axis[k] + 0.5),
              y = c(rep(fff$high[k],2),fff$high[k+1]), col = col1[j], lwd = lwd)
      }
      pos = pos + 17
    }

    lines(x = c(line_X_axis[103] - 0.5,line_X_axis[104],line_X_axis[104],line_X_axis[105] + 0.5),
          y = c(rep(fff$low[103],2),rep(fff$low[105],2)), col = col1[7], lwd = lwd)
    lines(x = c(line_X_axis[103] - 0.5,line_X_axis[104],line_X_axis[104],line_X_axis[105] + 0.5),
          y = c(rep(fff$high[103],2),rep(fff$high[105],2)), col = col1[7], lwd = lwd)
    
    for ( k in 107:112) {
      lines(x = c(line_X_axis[k] - 1.5,line_X_axis[k] + 1.5,line_X_axis[k] + 1.5),
            y = c(rep(fff$low[k],2),fff$low[k+1]), col = col1[8], lwd = lwd)
      lines(x = c(line_X_axis[k] - 1.5,line_X_axis[k] + 1.5,line_X_axis[k] + 1.5),
            y = c(rep(fff$high[k],2),fff$high[k+1]), col = col1[8], lwd = lwd)
    }
    
    lines(x = c(line_X_axis[114] - 1.5,line_X_axis[115],line_X_axis[115],line_X_axis[116] + 1.5),
          y = c(rep(fff$low[114],2),rep(fff$low[116],2)), col = col1[9], lwd = lwd)
    lines(x = c(line_X_axis[114] - 1.5,line_X_axis[115],line_X_axis[115],line_X_axis[116] + 1.5),
          y = c(rep(fff$high[114],2),rep(fff$high[116],2)), col = col1[9], lwd = lwd)
    
    for ( k in 118:123) {
      lines(x = c(line_X_axis[k] - 1.5,line_X_axis[k] + 1.5,line_X_axis[k] + 1.5),
            y = c(rep(fff$low[k],2),fff$low[k+1]), col = col1[10], lwd = lwd)
      lines(x = c(line_X_axis[k] - 1.5,line_X_axis[k] + 1.5,line_X_axis[k] + 1.5),
            y = c(rep(fff$high[k],2),fff$high[k+1]), col = col1[10], lwd = lwd)
    }
    
    for ( k in 125:131) {
      lines(x = c(line_X_axis[k] - 1.25,line_X_axis[k] + 1.25,line_X_axis[k] + 1.25),
            y = c(rep(fff$low[k],2),fff$low[k+1]), col = col1[11], lwd = lwd)
      lines(x = c(line_X_axis[k] - 1.25,line_X_axis[k] + 1.25,line_X_axis[k] + 1.25),
            y = c(rep(fff$high[k],2),fff$high[k+1]), col = col1[11], lwd = lwd)
    }
    #lines(x = line_X_axis[c(103,105)], y = fff$bar[c(103,105)], col = col[7], lwd = lwd)
    #lines(x = line_X_axis[107:112], y = fff$bar[107:112], col = col[8], lwd = lwd)
    #lines(x = line_X_axis[c(114,116)], y = fff$bar[c(114,116)], col = col[9], lwd = lwd)
    #lines(x = line_X_axis[118:123], y = fff$bar[118:123], col = col[10], lwd = lwd)
    #lines(x = line_X_axis[125:131], y = fff$bar[125:131], col = col[11], lwd = lwd)
    
    #lines(x = line_X_axis[c(103,105)], y = fff$low[c(103,105)], col = col[7], lwd = lwd/2)
    #lines(x = line_X_axis[107:112], y = fff$low[107:112], col = col[8], lwd = lwd/2)
    #lines(x = line_X_axis[c(114,116)], y = fff$low[c(114,116)], col = col[9], lwd = lwd/2)
    #lines(x = line_X_axis[118:123], y = fff$low[118:123], col = col[10], lwd = lwd/2)
    #lines(x = line_X_axis[125:131], y = fff$low[125:131], col = col[11], lwd = lwd/2)
    
    #lines(x = line_X_axis[c(103,105)], y = fff$high[c(103,105)], col = col[7], lwd = lwd/2)
    #lines(x = line_X_axis[107:112], y = fff$high[107:112], col = col[8], lwd = lwd/2)
    #lines(x = line_X_axis[c(114,116)], y = fff$high[c(114,116)], col = col[9], lwd = lwd/2)
    #lines(x = line_X_axis[118:123], y = fff$high[118:123], col = col[10], lwd = lwd/2)
    #lines(x = line_X_axis[125:131], y = fff$high[125:131], col = col[11], lwd = lwd/2)
  }
  
  pos = 1
  for (j in 1:6) {
    points(x = line_X_axis[pos:(pos+15)], y = fff$bar[pos:(pos+15)], col = col[j], pch = '-', cex = cex)
    pos = pos + 17
  }
  points(x = line_X_axis[c(103,105)], y = fff$bar[c(103,105)], col = col[7], pch = '-', cex = cex)
  points(x = line_X_axis[107:112], y = fff$bar[107:112], col = col[8], pch = '-', cex = cex)
  points(x = line_X_axis[c(114,116)], y = fff$bar[c(114,116)], col = col[9], pch = '-', cex = cex)
  points(x = line_X_axis[118:123], y = fff$bar[118:123], col = col[10], pch = '-', cex = cex)
  points(x = line_X_axis[125:131], y = fff$bar[125:131], col = col[11], pch = '-', cex = cex)
  
  col2 <- sapply(col, function(x) {
    tmpf <- colorRampPalette(c(x, "black"))
    return(tmpf(5)[3])
  })
  
  pos = 1
  for (j in 1:6) {
    lines(x = line_X_axis[pos:(15+pos)], y = rep(mean(fff$bar[pos:(15+pos)]),16), col = col2[j], lwd = 2*lwd, lty = 3)
    pos = pos + 17
  }
  lines(x = line_X_axis[c(103,105)], y = rep(mean(fff$bar[c(103,105)]),2), col = col2[7], lwd = 2*lwd,  lty = 3)
  lines(x = line_X_axis[107:112], y = rep(mean(fff$bar[107:112]),6), col = col2[8], lwd = 2*lwd, lty = 3)
  lines(x = line_X_axis[c(114,116)], y = rep(mean(fff$bar[c(114,116)]),2), col = col2[9], lwd = 2*lwd, lty = 3)
  lines(x = line_X_axis[118:123], y = rep(mean(fff$bar[118:123]),6), col = col2[10], lwd = 2*lwd, lty = 3)
  lines(x = line_X_axis[125:131], y = rep(mean(fff$bar[125:131]),7), col = col2[11], lwd = 2*lwd, lty = 3)
  
  
  
  title(main=plot_main)
  if (is.na(labels))
    if (log) labels = 10**at else labels = at
  if (log) {
    axis(side = 2, at = at, labels = labels, las = 2)
    abline(h = 0, lty = 2, col = 'red')
    abline(h = -1, lty = 2)
  }  else {
    axis(side = 2, at = at, labels = labels, las = 2)
    abline(h = 1,lty = 2, col = 'red')
    abline(h = 0.1,lty = 2)
  }
  #if (CI) {
  #  fff$low <- intersperse_119(low)
  #  fff$high <- intersperse_119(high)
  #  if (log) {
  #    fff$low <- fff$low / log(10)
  #    fff$high <- fff$high / log(10)
  #  }
  #  arrows(x0 = cumsum(c(rep(1,101),3,6,3,6,4,rep(3,6),4,6,3,6,4,rep(3,6),2,rep(2.5,7))), 
  #         y0 = fff$low, y1 = fff$high, code=3, length=0.01, lwd=1, angle=90, col='darkgrey', lty = 3)
  #}
}



plot_small_barplot <- function(y, at=c(-1,0,1), plot_main = '', log = T, CI = F, low = NA, high = NA, indel = T,
                              col = c("#2EBAED","#000000","#DE1C14","#D4D2D2","#ADCC54","#F0D0CE",'orange','darkmagenta','purple')) {
  if (indel) {
    cl <- rep(col,c(rep(16,6),3,2,3))
    cl1 <- c(cl[1:16],NA,cl[17:32],NA,cl[33:48],NA,cl[49:64],NA,cl[65:80],NA,cl[81:96],NA,
             cl[97],NA,cl[98],NA,cl[99],NA,cl[100],NA,cl[101],NA,cl[102],NA,cl[103],NA,cl[104])
    fff <- data.frame(bar = c(y[1:16],NA,
                              y[17:32],NA,
                              y[33:48],NA,
                              y[49:64],NA,
                              y[65:80],NA,
                              y[81:96],NA,
                              y[97],NA,y[98],NA, y[99],NA,
                              y[100],NA,y[101],NA,
                              y[102],NA,y[103],NA,y[104]))
    if (log) fff$bar <- fff$bar / log(10)
  #                  low = Wsig_all_low[match(y,colnames(W_all)),],
  #                  high = Wsig_all_up[match(y,colnames(W_all)),])
    foo <- barplot(fff$bar, col = cl1,
                   width = c(rep(1,101),2,
                             4,2,4,2,4,5,
                             6,3,6,5,
                             4,2,4,2,4),
                   border=NA, main=NA,
                   ylim = c(min(c(at,min(fff$bar)),na.rm=T),max(c(at,max(fff$bar)),na.rm=T)),
                   yaxt = 'n')
    title(main=plot_main)
    if (log) axis(side = 2, at = at, labels = 10**at, las = 2) else axis(side = 2, at = at, labels = at, las = 2)
    if (CI) {
      fff$low <- c(low[1:16],NA,low[17:32],NA,low[33:48],NA,
                   low[49:64],NA,low[65:80],NA,low[81:96],NA,
                   low[97],NA,low[98],NA, low[99],NA,
                   low[100],NA,low[101],NA,
                   low[102],NA,low[103],NA,low[104])
      fff$high = c(high[1:16],NA,
                   high[17:32],NA,
                   high[33:48],NA,
                   high[49:64],NA,
                   high[65:80],NA,
                   high[81:96],NA,
                   high[97],NA,high[98],NA, high[99],NA,
                   high[100],NA,high[101],NA,
                   high[102],NA,high[103],NA,high[104])
      if (log) {
        fff$low <- fff$low / log(10)
        fff$high <- fff$high / log(10)
      }
      arrows(x0 = foo, y0 = fff$low, y1 = fff$high, code=3, length=0.01, lwd=1, angle=90, col='darkgrey')
    }
  }
  else {
    cl <- rep(col,c(rep(16,6)))
    cl1 <- c(cl[1:16],NA,cl[17:32],NA,cl[33:48],NA,cl[49:64],NA,cl[65:80],NA,cl[81:96])
    fff <- data.frame(bar = c(y[1:16],NA,
                              y[17:32],NA,
                              y[33:48],NA,
                              y[49:64],NA,
                              y[65:80],NA,
                              y[81:96]))
    if (log) fff$bar <- fff$bar / log(10)
    #                  low = Wsig_all_low[match(y,colnames(W_all)),],
    #                  high = Wsig_all_up[match(y,colnames(W_all)),])
    foo <- barplot(fff$bar, col = cl1,
                   width = c(rep(1,101)),
                   border=NA, main=NA,
                   ylim = c(min(c(at,min(fff$bar)),na.rm=T),max(c(at,max(fff$bar)),na.rm=T)),
                   yaxt = 'n')
    title(main=plot_main)
    if (log) axis(side = 2, at = at, labels = 10**at, las = 2) else axis(side = 2, at = at, labels = at, las = 2)
    if (CI) {
      fff$low <- c(low[1:16],NA,low[17:32],NA,low[33:48],NA,
                   low[49:64],NA,low[65:80],NA,low[81:96])
      fff$high = c(high[1:16],NA,high[17:32],NA,
                   high[33:48],NA,high[49:64],NA,
                   high[65:80],NA,high[81:96])
      if (log) {
        fff$low <- fff$low / log(10)
        fff$high <- fff$high / log(10)
      }
      arrows(x0 = foo, y0 = fff$low, y1 = fff$high, code=3, length=0.01, lwd=1, angle=90, col='darkgrey')
    }
  }
}

intersperse_96 <- function(y) {
  return(c(y[1:16],NA,
           y[17:32],NA,
           y[33:48],NA,
           y[49:64],NA,
           y[65:80],NA,
           y[81:96]))
}
intersperse_104 <- function(y) {
  return(c(y[1:16],NA,
           y[17:32],NA,
           y[33:48],NA,
           y[49:64],NA,
           y[65:80],NA,
           y[81:96],NA,
           y[97],NA,
           y[98],NA,
           y[99],NA,
           y[100],NA,
           y[101],NA,
           y[102],NA,
           y[103],NA,
           y[104]))
}

intersperse_105 <- function(y) {
  return(c(y[1:16],NA,
           y[17:32],NA,
           y[33:48],NA,
           y[49:64],NA,
           y[65:80],NA,
           y[81:96],NA,
           y[97],NA,
           y[98],NA,
           y[99],NA,
           y[100],NA,
           y[101],NA,
           y[102],NA,
           y[103],NA,
           y[104],NA,
           y[105]))
}


plot_small_line <- function(y, at=c(-1,0,1), plot_main = '', log = T,CI=F,low,high, cex = 3, lwd = 4, labels=NA, indels = F,
                           col = c("#2EBAED", "#000000", "#DE1C14","#D4D2D2", "#ADCC54", "#F0D0CE","orange",'darkmagenta',"purple")) {
  if (indels) {
    print('no indels yet!')
    return()
  }
  
  fff <- data.frame(bar = intersperse_96(y[1:96]))
  if (log) fff$bar <- fff$bar / log(10)
  line_X_axis <- cumsum(c(rep(1,101)))
  plot(x = line_X_axis[1:16], y = fff$bar[1:16],
       col = col[1], type = 'l', lwd = lwd,
       xlim = c(0,max(line_X_axis)),
       ylim = c(min(c(at,min(fff$bar,na.rm=T)),na.rm=T),max(c(at,max(fff$bar,na.rm=T)),na.rm=T)),
       yaxt = 'n', xaxt = 'n', xlab = '', ylab = '', bty='n')
  if (CI) {
    fff$low <- intersperse_96(low[1:96])
    fff$high <- intersperse_96(high[1:96])
    if (log) {
      fff$low <- fff$low / log(10)
      fff$high <- fff$high / log(10)
    }
    col1 <- sapply(col, function(x) {
      tmpf <- colorRampPalette(c(x, "white"))
      return(tmpf(5)[3])
    })
    polygon(x = c(line_X_axis[1:16],line_X_axis[16:1]), y = c(fff$high[1:16],fff$low[16:1]), col = col1[1], border = NA)
    polygon(x = c(line_X_axis[18:33],rev(line_X_axis[18:33])), y = c(fff$high[18:33],fff$low[33:18]), col = col1[2],border = NA)
    polygon(x = c(line_X_axis[35:50],line_X_axis[50:35]), y =c(fff$high[35:50],fff$low[50:35]), col = col1[3], border = NA)
    polygon(x = c(line_X_axis[52:67],line_X_axis[67:52]), y = c(fff$high[52:67],fff$low[67:52]), col = col1[4], border = NA)
    polygon(x = c(line_X_axis[69:84],line_X_axis[84:69]), y = c(fff$high[69:84],fff$low[84:69]), col = col1[5], border = NA)
    polygon(x = c(line_X_axis[86:101],line_X_axis[101:86]), y = c(fff$high[86:101],fff$low[101:86]), col = col1[6],border = NA)
  }
  lines(x = line_X_axis[1:16], y = fff$bar[1:16], col = col[1], lwd = lwd)
  lines(x = line_X_axis[18:33], y = fff$bar[18:33], col = col[2], lwd = lwd)
  lines(x = line_X_axis[35:50], y = fff$bar[35:50], col = col[3], lwd = lwd)
  lines(x = line_X_axis[52:67], y = fff$bar[52:67], col = col[4], lwd = lwd)
  lines(x = line_X_axis[69:84], y = fff$bar[69:84], col = col[5], lwd = lwd)
  lines(x = line_X_axis[86:101], y = fff$bar[86:101], col = col[6], lwd = lwd)

  title(main=plot_main)
  if (is.na(labels))
    if (log) labels = 10**at else labels = at
  if (log) {
    axis(side = 2, at = at, labels = labels, las = 2)
    abline(h = 0, lty = 2, col = 'red')
    abline(h = -1, lty = 2)
  }  else {
    axis(side = 2, at = at, labels = labels, las = 2)
    abline(h = 1,lty = 2, col = 'red')
    abline(h = 0.1,lty = 2)
  }

}

##################

interaction_effect_plot <- function(y, at=c(-1,0,1), plot_main = '', log = T,CI=F,low,high, cex = 3,  lwd = 4, labels=NA, lwd.means = 4,
                           col = c("#2EBAED","#000000","#DE1C14","#D4D2D2","#ADCC54","#F0D0CE","brown","#8DD3C7","goldenrod","#BEBADA","darkmagenta")) {
  if (length(y) != 119) return('Working progress!')
  fff <- data.frame(bar = intersperse_119(y))
  if (log) fff$bar <- fff$bar / log(10)
  line_X_axis <- cumsum(c(rep(1,101),4,2,6,2,4,rep(2.5,6),5,2,5,2,5,rep(2.5,6),2,rep(2,7)))
  plot(x = line_X_axis, y = rep(NA,length(line_X_axis)),
       xlim = c(0,max(line_X_axis)),
       ylim = c(min(c(at,min(fff$bar,na.rm=T)),na.rm=T),max(c(at,max(fff$bar,na.rm=T)),na.rm=T)),
       yaxt = 'n', xaxt = 'n', xlab = '', ylab = '', bty='n')
  
  fff[fff < -1] <- -1
  
  if (CI) {
    fff$low <- intersperse_119(low)
    fff$high <- intersperse_119(high)
    if (log) {
      fff$low <- fff$low / log(10)
      fff$high <- fff$high / log(10)
    }
    
    fff[fff < -1] <- -1
    
    col1 <- sapply(col, function(x) {
      tmpf <- colorRampPalette(c(x, "white"))
      return(tmpf(7)[6])
    })
    
    pos = 1
    for (j in 1:6) {
      for (k in pos:(pos+15))
        polygon(x = c(line_X_axis[k] - 1, rep(line_X_axis[k], 2),
                      rep(line_X_axis[k] - 1,2)),
                y = c(rep(fff$high[k],2),
                      rep(fff$low[k],2), fff$high[k]),
                col = col1[j], border = NA)
      pos = pos + 17
    }
    
    for (k in c(103,105))
      polygon(x = c(line_X_axis[k] - 2, rep(line_X_axis[k],2), 
                    rep(line_X_axis[k] - 2,2)),
              y = c(rep(fff$high[k],2),
                    rep(fff$low[k],2), fff$high[k]),
              col = col1[7], border = NA)
      
    for (k in c(114,116))
      polygon(x = c(line_X_axis[k] - 2, rep(line_X_axis[k],2), 
                    rep(line_X_axis[k] - 2,2)),
              y = c(rep(fff$high[k],2),
                    rep(fff$low[k],2), fff$high[k]),
              col = col1[9], border = NA)
    
    for (k in c(107:112))
      polygon(x = c(line_X_axis[k] - 2, rep(line_X_axis[k],2), 
                    rep(line_X_axis[k] - 2,2)),
              y = c(rep(fff$high[k],2),
                    rep(fff$low[k],2), fff$high[k]),
              col = col1[8], border = NA)
    
    for (k in c(118:123))
      polygon(x = c(line_X_axis[k] - 2, rep(line_X_axis[k],2), 
                    rep(line_X_axis[k] - 2,2)),
              y = c(rep(fff$high[k],2),
                    rep(fff$low[k],2), fff$high[k]),
              col = col1[10], border = NA)

    for (k in c(125:131))
      polygon(x = c(line_X_axis[k] - 2, rep(line_X_axis[k],2), 
                    rep(line_X_axis[k] - 2,2)),
              y = c(rep(fff$high[k],2),
                    rep(fff$low[k],2), fff$high[k]),
              col = col1[11], border = NA)

  }
  # All means
  col1 <- sapply(col, function(x) {
    tmpf <- colorRampPalette(c(x, "white"))
    return(tmpf(7)[4])
  })
  pos = 1
  for (j in 1:6) {
    for (k in pos:(pos+15))
      lines(x = c(line_X_axis[k] - 1, line_X_axis[k]),
            y = rep(fff$bar[k],2), col = col1[j], lwd = lwd.means/2)
    pos = pos + 17
  }
  for (k in c(103,114))
    lines(x = c(line_X_axis[k] - 2, line_X_axis[k]), y = rep(fff$bar[k],2), col = col1[ifelse(k==103,7,9)], lwd = lwd.means/2)
  for (k in c(105,116))
    lines(x = c(line_X_axis[k] - 2, line_X_axis[k]), y = rep(fff$bar[k],2), col = col1[ifelse(k==105,7,9)], lwd = lwd.means/2)
  for (k in c(107:112,118:123))
    lines(x = c(line_X_axis[k] - 2, line_X_axis[k]), y = rep(fff$bar[k],2), col = col1[ifelse(k<117,8,10)], lwd = lwd.means/2)
  for (k in c(125:131))
    lines(x = c(line_X_axis[k] - 2, line_X_axis[k]), y = rep(fff$bar[k],2), col = col1[11], lwd = lwd.means/2)
  
  # Significant means
  col2 <- sapply(col, function(x) {
    tmpf <- colorRampPalette(c(x, "black"))
    return(tmpf(5)[2])
  })
  pos = 1
  for (j in 1:6) {
    pos1 <- c(pos:(pos+15))[which(fff$high[pos:(pos+15)] < 0 | fff$low[pos:(pos+15)] > 0)]
    for (k in pos1)
      lines(x = c(line_X_axis[k] - 1, line_X_axis[k]),
            y = rep(fff$bar[k],2), col = col2[j], lwd = lwd.means)
    pos = pos + 17
  }
  pos1 <- which(fff$high < 0 | fff$low > 0)
  pos1 <- pos1[pos1>101]
  
  for (k in pos1) {
    if (k %in% c(103,105))
      lines(x = c(line_X_axis[k] - 2, line_X_axis[k]), y = rep(fff$bar[k],2), col = col2[7], lwd = lwd.means)
    if (k %in% c(107:112))
      lines(x = c(line_X_axis[k] - 2, line_X_axis[k]), y = rep(fff$bar[k],2), col = col2[8], lwd = lwd.means)
    else if (k %in% c(114,116))
      lines(x = c(line_X_axis[k] -2, line_X_axis[k]), y = rep(fff$bar[k],2), col = col2[9], lwd = lwd.means)
    else if (k %in% c(118:123))
      lines(x = c(line_X_axis[k] - 2, line_X_axis[k]), y = rep(fff$bar[k],2), col = col2[10], lwd = lwd.means)
    else 
      lines(x = c(line_X_axis[k] - 2, line_X_axis[k]), y = rep(fff$bar[k],2), col = col2[11], lwd = lwd.means)
  }
  
  
  # Class means
  pos = 1
  for (j in 1:6) {
    lines(x = c(line_X_axis[pos]-1,line_X_axis[pos+15]), y = rep(median(fff$bar[pos:(15+pos)]),2))
    pos = pos + 17
  }
  lines(x = c(line_X_axis[103]-2,line_X_axis[105]), y = rep(median(fff$bar[c(103,105)]),2))
  lines(x = c(line_X_axis[107]-2,line_X_axis[112]), y = rep(median(fff$bar[107:112]),2))
  lines(x = c(line_X_axis[114]-2,line_X_axis[116]), y = rep(median(fff$bar[c(114,116)]),2))
  lines(x = c(line_X_axis[118]-2,line_X_axis[123]), y = rep(median(fff$bar[118:123]),2))
  lines(x = c(line_X_axis[125]-2,line_X_axis[131]), y = rep(median(fff$bar[125:131]),2))
  
  
  
  title(main=plot_main)
  if (is.na(labels))
    if (log) labels = 10**at else labels = at
    axis(side = 2, at = c(log10(do.call('c',lapply(c((-1):(max(at)-1)), function(i) 
      c(2:9) * 10**i)))), labels = c(rep('',8*(max(at) - min(at)))), las = 2, pos = -0.5, tck = -0.02)
    axis(side = 2, at = c(at), labels = c(as.character(labels)), las = 2, pos = -0.5, tck = -0.03)
    #axis(side = 2, col = 'white', tcl = 0, labels = NA, lwd = 3, pos = -0.5)
    lines(x = c(-0.5,185.5), y = c(0,0),lty = 2, col = 'red')
  }  else {
    axis(side = 2, at = at, labels = labels, las = 2, pos = -0.5)
    #axis(side = 2, col = 'white', lwd = 3, tcl = 0, labels = NA, pos = -0.5)
    lines(x = c(-0.5,187), y = c(1,1),lty = 2, col = 'red')
  }
  
}

interaction_effect_plot_human <- function(y, at=c(-1,0,1), plot_main = '', log = T,CI=F,low,high, cex = 3, lwd = 4, labels=NA, lwd.means = 4,
                                    col = c("#2EBAED","#000000","#DE1C14","#D4D2D2","#ADCC54","#F0D0CE","orange",'darkmagenta',"purple",'brown')) {
  if (length(y) != 105) return('Working progress!')
  fff <- data.frame(bar = intersperse_105(y))
  if (log) fff$bar <- fff$bar / log(10)
  line_X_axis <- cumsum(c(rep(1,101),3,2,3,2,3,2,6,2,5,2,6,2,3,2,3,2,3,2)) - 1
  plot(x = line_X_axis, y = rep(NA,length(line_X_axis)),
       xlim = c(0,max(line_X_axis)),
       ylim = c(min(at),max(at)),
       yaxt = 'n', xaxt = 'n', xlab = '', ylab = '', bty='n')
  fff[fff < -1] <- -1
  if (CI) {
    fff$low <- intersperse_105(low)
    fff$high <- intersperse_105(high)
    if (log) {
      fff$low <- fff$low / log(10)
      fff$high <- fff$high / log(10)
    }
    fff[fff < -1] <- -1
    col1 <- sapply(col, function(x) {
      tmpf <- colorRampPalette(c(x, "white"))
      return(tmpf(7)[6])
    })
    pos = 1
    for (j in 1:6) {
      for (k in pos:(pos+15))
        polygon(x = c(line_X_axis[k] - 1, rep(line_X_axis[k], 2),
                      rep(line_X_axis[k] - 1,2)),
                y = c(rep(fff$high[k],2),
                      rep(fff$low[k],2), fff$high[k]),
                col = col1[j], border = NA)
      pos = pos + 17
    }
    for (k in c(103,105,107))
      polygon(x = c(line_X_axis[k] - 2, rep(line_X_axis[k],2), 
                    rep(line_X_axis[k] - 2,2)),
              y = c(rep(fff$high[k],2),
                    rep(fff$low[k],2), fff$high[k]),
              col = col1[7], border = NA)
    
    for (k in c(109,111))
      polygon(x = c(line_X_axis[k] - 2, rep(line_X_axis[k],2), 
                    rep(line_X_axis[k] - 2,2)),
              y = c(rep(fff$high[k],2),
                    rep(fff$low[k],2), fff$high[k]),
              col = col1[8], border = NA)
    
    for (k in c(113,115,117))
      polygon(x = c(line_X_axis[k] - 2, rep(line_X_axis[k],2), 
                    rep(line_X_axis[k] - 2,2)),
              y = c(rep(fff$high[k],2),
                    rep(fff$low[k],2), fff$high[k]),
              col = col1[9], border = NA)
    
    polygon(x = c(line_X_axis[119] - 2, rep(line_X_axis[119],2), 
                  rep(line_X_axis[119] - 2,2)),
            y = c(rep(fff$high[119],2),
                  rep(fff$low[119],2), fff$high[119]),
            col = col1[10], border = NA)
  }
  
  col1 <- sapply(col, function(x) {
    tmpf <- colorRampPalette(c(x, "white"))
    return(tmpf(7)[4])
  })
  pos = 1
  for (j in 1:6) {
    for (k in pos:(pos+15))
      lines(x = c(line_X_axis[k] - 1, line_X_axis[k]),
            y = rep(fff$bar[k],2), col = col1[j], lwd = lwd.means/2)
    pos = pos + 17
  }
  for (k in c(103,105,107))
    lines(x = c(line_X_axis[k] - 2, line_X_axis[k]), y = rep(fff$bar[k],2), col = col1[7], lwd = lwd.means/2)
  for (k in c(109,111))
    lines(x = c(line_X_axis[k] - 2, line_X_axis[k]), y = rep(fff$bar[k],2), col = col1[8], lwd = lwd.means/2)
  for (k in c(113,115,117))
    lines(x = c(line_X_axis[k] - 2, line_X_axis[k]), y = rep(fff$bar[k],2), col = col1[9], lwd = lwd.means/2)
  lines(x = c(line_X_axis[119] - 2, line_X_axis[119]), y = rep(fff$bar[119],2), col = col1[10], lwd = lwd.means/2)
  
  col2 <- sapply(col, function(x) {
    tmpf <- colorRampPalette(c(x, "black"))
    return(tmpf(5)[2])
  })
  pos = 1
  for (j in 1:6) {
    pos1 <- c(pos:(pos+15))[which(fff$high[pos:(pos+15)] < 0 | fff$low[pos:(pos+15)] > 0)]
    for (k in pos1)
      lines(x = c(line_X_axis[k] - 1, line_X_axis[k]),
            y = rep(fff$bar[k],2), col = col2[j], lwd = lwd.means)
    pos = pos + 17
  }
  pos1 <- which(fff$high < 0 | fff$low > 0)
  pos1 <- pos1[pos1>101]
  for (k in pos1) {
    if (k %in% c(103,105,107))
      lines(x = c(line_X_axis[k] - 2, line_X_axis[k]), y = rep(fff$bar[k],2), col = col2[7], lwd = lwd.means)
    else if (k %in% c(109,111))
      lines(x = c(line_X_axis[k] - 2, line_X_axis[k]), y = rep(fff$bar[k],2), col = col2[8], lwd = lwd.means)
    else if (k %in% c(113,115,117))
      lines(x = c(line_X_axis[k] - 2, line_X_axis[k]), y = rep(fff$bar[k],2), col = col2[9], lwd = lwd.means)
    else if (k == 119)
      lines(x = c(line_X_axis[k] - 2, line_X_axis[k]), y = rep(fff$bar[k],2), col = col2[10], lwd = lwd.means)
  }
  
  pos = 1
  for (j in 1:6) {
    lines(x = c(line_X_axis[pos]-1,line_X_axis[pos+15]), y = rep(median(fff$bar[pos:(15+pos)]),2))
    pos = pos + 17
  }
  lines(x = c(line_X_axis[103]-2,line_X_axis[107]), y = rep(median(fff$bar[c(103,105,107)]),2))
  lines(x = c(line_X_axis[109]-2,line_X_axis[111]), y = rep(median(fff$bar[c(109,111)]),2))
  lines(x = c(line_X_axis[113]-2,line_X_axis[117]), y = rep(median(fff$bar[c(113,115,117)]),2))

  title(main=plot_main)
  if (is.na(labels))
    if (log) labels = 10**at else labels = at
  if (log) {
    axis(side = 2, at = at, labels = labels, las = 2, pos = -1.5)
    axis(side = 2, lwd = 3, tcl = 0, labels = NA, pos = -1.5)
    lines(x = c(-0.5,148), y = c(0,0),lty = 2, col = 'red')
  }  else {
    axis(side = 2, at = at, labels = labels, las = 2, pos = -1.5)
    axis(side = 2, lwd = 3, tcl = 0, labels = NA, pos = -1.5)
    lines(x = c(-0.5,148), y = c(1,1),lty = 2, col = 'red')
  }
  
}


#  to.show <- data.frame(t(S_est))
#  for (j in 1:nrow(beta_est)) {
#    to.show <- cbind(to.show, S_est[nrow(S_est),] * exp(beta_est[j,]))
#  }
#  colnames(to.show) <- c(paste0('Sig',1:nrow(S_est)), paste0('Sig*',rownames(beta_est)))
#  plot_fc(to.show, norm = n)
#}

# devtools::install_github("thomasp85/patchwork") for combining ggplots plots