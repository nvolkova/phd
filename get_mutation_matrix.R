get_mutation_matrix <- function(vlist,ref_genome) {
   TYPES = c("C>A","C>G","C>T","T>A","T>C","T>G")
   TYPES.FULL = c(paste0(rep(rep(c("A","C","G","T"),each=4),6),"[",rep(TYPES,each=16),"]",rep(c("A","C","G","T"),24)),paste0("DEL_",c("A","C","G","T")),paste0("INS_",c("A","C","G","T")))
   mm <- matrix(0,nrow=length(vlist),ncol=104,dimnames=list(names(vlist),TYPES.FULL))
   sub_list <- sapply(vlist, function(vcf) {
     vcf[width(vcf$REF)==1 & width(unlist(vcf$ALT))==1,]
   })
   del_list <- sapply(vlist, function(vcf) {
     vcf[width(vcf$REF)==2 & width(unlist(vcf$ALT))==1,]
   })
   ins_list <- sapply(vlist, function(vcf) {
     vcf[width(vcf$REF)==1 & width(unlist(vcf$ALT))==2,]
   })
   for (n in 1:length(vlist)) {
     counts <- table(type_context(sub_list[[n]], ref_genome))
     for (a in rownames(counts)) {
       tmp = unlist(strsplit(a,split="[>]"))
       inds <- colnames(counts)[counts[a,]>0]
       columns = as.vector(sapply(inds, function(x) paste(substr(x,1,1),"[",a,"]",substr(x,nchar(x),nchar(x)),sep="")))
       mm[n,columns] = counts[a,inds]
     }
     mm[n,97:100] <- table(substr(unlist(ins_list[[n]]$ALT),2,2))[c("A","C","G","T")]
     mm[n,101:104] <- table(substr(del_list[[n]]$REF,2,2))[c("A","C","G","T")]
   }
   mm[is.na(mm)]=0
   return(mm)
}