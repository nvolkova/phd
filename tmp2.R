# PolK mutants in ICGC

t <- read.delim('~/Downloads/ICGC-POLK/simple_somatic_mutation.open.tsv')


mutmat <- matrix(0,nrow=length( unique(as.character(t$icgc_donor_id))),ncol=6,dimnames=list( unique(as.character(t$icgc_donor_id)), c('C>A','C>G','C>T','T>A','T>C','T>G')))
for (donor in unique(as.character(t$icgc_donor_id))) {
  tmp <- unique(t[t$icgc_donor_id==donor,c(1:3,9:11,14:17)])
  muts <- paste(tmp[tmp$mutation_type=='single base substitution','reference_genome_allele'],'>',tmp[tmp$mutation_type=='single base substitution','mutated_to_allele'],sep='')
  muts[muts=='A>C'] <- 'T>G'
  muts[muts=='A>G'] <- 'T>C'
  muts[muts=='A>T'] <- 'T>A'
  muts[muts=='G>A'] <- 'C>T'
  muts[muts=='G>C'] <- 'C>G'
  muts[muts=='G>T'] <- 'C>A'
  mutmat[donor,] <- table(muts)[colnames(mutmat)]
}

barplot(t(mutmat),col=c('blue','black','red','grey','green','pink'))







to.show <- cbind(colMeans(big.spectrum[c('CD0926b','CD0926c','CD0926d'),], na.rm=T),
                 colMeans(big.spectrum[c('CD0986b','CD0986c','CD0986d'),], na.rm=T),
                 colMeans(big.spectrum[c('CD0989b','CD0989c','CD0989d'),], na.rm=T),
                 colMeans(big.spectrum[c('CD0920b','CD0920c','CD0920d'),], na.rm=T),
                 colMeans(big.spectrum[c('CD0992b','CD0992c','CD0992d'),], na.rm=T),
                 colMeans(big.spectrum[c('CD0957b','CD0957c','CD0957d'),], na.rm=T),
                 colMeans(big.spectrum[c('CD0944b','CD0944c','CD0944d'),], na.rm=T),
                 colMeans(big.spectrum[c('CD0950b','CD0950c','CD0950d'),], na.rm=T),
                 colMeans(big.spectrum[c('CD0974b','CD0974c','CD0974d'),], na.rm=T),
                 colMeans(big.spectrum[c('CD0970b','CD0970c','CD0970d'),], na.rm=T))
colnames(to.show) = c('N2 25', 'xpa-1 5', 'xpc-1 5', 'csb-1 25', 'xpf-1 2', 'polk-1 25', 'pole-4 25', 'polh-1 25', 'rev-3 25', 'rev-1 10')
p <- plot_prof_wb(to.show)
ggsave('~/Desktop/AAprof.jpg', device = 'jpeg', plot = p, width=12, height = 8)
