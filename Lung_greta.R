##################################################################
## MCMC sampling model for interaction effects in lung cancers  ##
## N. Volkova, EBI, November 2018                               ##
##################################################################

library(VariantAnnotation)
library(MASS)
source('~/Desktop/Git/phd/useful_functions.R')
library(ggplot2)
library(rstan)
library(reshape2)
source('~/Desktop/Git/phd/plot_sigs.R')
library(greta)
library(openxlsx)

cosine <- function(x,y) {
  x %*% y / sqrt(sum(x**2)) / sqrt(sum(y**2))
}
sp_url <- paste("http://cancer.sanger.ac.uk/cancergenome/assets/",
                "signatures_probabilities.txt", sep = "")
cancer_signatures = read.table(sp_url, sep = "\t", header = TRUE)
cancer_signatures = cancer_signatures[order(cancer_signatures[,1]),]
types <- as.character(cancer_signatures$Trinucleotide) # trinucleotide classes
types.full <- as.character(cancer_signatures$Somatic.Mutation.Type) # substitution types
row.names(cancer_signatures) <- types.full
cancer_signatures = as.matrix(cancer_signatures[,4:33])

bigmat <- read.table('/Volumes/r/nadia/TCGA.caveman.matrix.dat',sep='\t')
metadata <- read.table('/Volumes/r/nadia/TCGA_caveman_patients_to_files.dat', sep='\t', header=T)

lung <- metadata$tumour[metadata$project %in% c('LUAD','LUSC')]
load('~/lung.donor.mut.mat.RData')

# Labels
ner.genes <- c('XAB2','XPC','XPA','ERCC1','ERCC2','ERCC3','ERCC4','ERCC5','ERCC6','ERCC8','RAD23A','DDB1','DDB2','RAD23B')
tcr <- c('XAB2','ERCC6','ERCC8')
ggr <- c('XPC', 'RAD23A','DDB1','DDB2','RAD23B')
tls.genes <- c('POLH','POLK','REV1','REV3','POLI','POLB')

# Expression
check.exp <- function(z, tr, m = metadata) {
  return(which(m[,match(z,colnames(m))] < tr))
}
for (i in 6:ncol(metadata)) {
  metadata[,i] <- metadata[,i] / median(metadata[,i], na.rm = T)
  #  metadata[,i] <- rank(metadata[,i], na.last='keep') / length(metadata[,i])
}
exp.ner <- metadata$tumour[sort(unique(unlist(sapply(ner.genes, check.exp, tr = 0.3))))]
exp.tcr <- metadata$tumour[sort(unique(unlist(sapply(tcr, check.exp, tr = 0.3))))]
exp.ggr <- metadata$tumour[sort(unique(unlist(sapply(ggr, check.exp, tr = 0.3))))]
exp.polh <- metadata$tumour[sort(check.exp('POLH',0.3))]
exp.polk <- metadata$tumour[sort(check.exp('POLK',0.3))]
exp.rev3 <- metadata$tumour[sort(check.exp('REV3L',0.3))]

ids <- read.table('~/Downloads/ICGCtoTCGA.tsv', sep = '\t', header = T)
samples <- list()
for (x in c(ner.genes,'POLH','REV3L','POLK')) {
  tmp <- read.delim(paste0('~/TCGAmutations/results/',x,'.txt'), sep='\t', header=T)
  tmp <- tmp[tmp$IMPACT %in% c('HIGH') & grepl(x, tmp$SYMBOL),]
  if (nrow(tmp)==0) next
  samples[[x]] <- unique(sapply(as.character(tmp$Uploaded_variation), 
                                function(y) substr(unlist(strsplit(y,split='[:]'))[1],1,12)))
}

# Germline

SNPs <- list()
tmp <- read.table('~/TCGAmutations/results/germline_SKCM_NER_38.txt', sep = '\t', header = T)
for (x in ner.genes) {
  SNPs[[x]] <- unique(as.character(tmp$Uploaded_variation[tmp$IMPACT %in% c('HIGH') & grepl(x, tmp$SYMBOL)]))
}
SNPs <- SNPs[sapply(SNPs,length)>0]

cancername = 'SKCM'
snps_per_sample <- matrix(0,nrow = length(skcm), ncol = length(SNPs), dimnames = list(skcm, names(SNPs)))
for (genename in names(SNPs)) {
  f <- list.files(paste0('/Volumes/SNP/',genename,'/',cancername))
  for (i in 1:length(skcm)) {
    file <- f[grep(skcm[i],f)]
    if (length(file)==0) next
    if (file.info(paste0('/Volumes/SNP/',genename,'/',cancername,'/',file))$size==0) next
    vcf <- read.delim(paste0('/Volumes/SNP/',genename,'/',cancername,'/',file), header = F)
    identifier <- paste0(vcf$V1,':',vcf$V2,'_',vcf$V4,'/',vcf$V5) 
    het <- substr(vcf$V10,1,3)[identifier %in% SNPs[[genename]]]
    if (length(het)>0) print(c(skcm[i],het,genename))
    snps_per_sample[skcm[i], genename] <- length(intersect(identifier, SNPs[[genename]]))
  }
  print(genename)
}
colSums(snps_per_sample>0)


# Damaging mutations

copynumber <- rbind(read.delim('LUAD.ascat.segments.txt', sep = '\t', header = F),read.delim('LUSC.ascat.segments.txt', sep = '\t', header = F))
colnames(copynumber) <- c('Sample','Chromosome','Start','End','MajorCN','MinorCN')
copynumber <- GRanges(copynumber)
ids <- read.table('~/Downloads/ICGCtoTCGA.tsv', sep = '\t', header = T)
samples <- list()
for (x in ner.genes) {
  tmp <- read.delim(paste0('~/TCGAmutations/results/',x,'.txt'), sep='\t', header=T)
  tmp <- tmp[tmp$IMPACT %in% c('HIGH','MODERATE') & grepl(x, tmp$SYMBOL),]
  if (nrow(tmp)==0) next
  uv_mutations <- read.delim(paste0('/Volumes/r/nadia/tcga_melanoma/',x,'.dat'), sep='\t', header = T)
  identifier <- paste0(uv_mutations$Sample,'-',uv_mutations$Chrom,':',uv_mutations$Pos,'_',uv_mutations$Ref,'/',uv_mutations$Alt)
  samples[[x]] <- unique(sapply(as.character(tmp$Uploaded_variation), function(y) substr(unlist(strsplit(y,split='[:]'))[1],1,12)))
  ind <- match(unique(as.character(tmp$Uploaded_variation)), identifier)
  ind <- ind[!is.na(ind)]
  print(c(x,as.character(uv_mutations$GT.Tum[ind])))
  samples[[x]] <- unlist(sapply(samples[[x]], function(y) skcm[grep(y,skcm)]))
  cn <- copynumber[queryHits(findOverlaps(copynumber, genes_of_interest[genes_of_interest$name==x]))]
  cn <- cn[cn$Sample %in% names(samples[[x]])]
  samples[[paste0(x,'LOH')]] <- samples[[x]][as.character(cn$Sample[cn$MinorCN==0])]
}


data <- list()

#NER <- sapply(ner.genes, function(x) unique(tab$Sample))

data$NER <- (substr(rownames(lung.mut.mat),6,12) %in% substr(unlist(samples[ner.genes]),6,12)) | (rownames(lung.mut.mat) %in% exp.ner)
data$TCR <- (substr(rownames(lung.mut.mat),6,12) %in% substr(unlist(samples[tcr]),6,12)) | (rownames(lung.mut.mat) %in% exp.tcr)
data$GGR <- (substr(rownames(lung.mut.mat),6,12) %in% substr(unlist(samples[ggr]),6,12)) | (rownames(lung.mut.mat) %in% exp.ggr)
data$POLH <- (substr(rownames(lung.mut.mat),6,12) %in% substr(samples[['POLH']],6,12)) | (rownames(lung.mut.mat) %in% exp.polh)
data$POLK <- (substr(rownames(lung.mut.mat),6,12) %in% substr(samples[['POLK']],6,12)) | (rownames(lung.mut.mat) %in% exp.polk)
data$REV3 <- (substr(rownames(lung.mut.mat),6,12) %in% substr(samples[['REV3L']],6,12)) | (rownames(lung.mut.mat) %in% exp.rev3)

colSums(as.data.frame(data))

lung_data <- list(
  N = nrow(lung.mut.mat),
  R = ncol(lung.mut.mat),
  K = 4,
  X = model.matrix( ~ ., data = as.data.frame(data))[,-c(1,5:6)],
  y = as.matrix(lung.mut.mat),
  constr = max(rowSums(lung.mut.mat))
)