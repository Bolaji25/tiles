library(rtracklayer)
#to split genome-wide to autosomes and X chromosome
dsmf<- import.bw("../Bigwiggle/sevinc/sm_allCs_combined_M_ucsc.bw")
a<- which(seqnames(dsmf)=="chrX")
Xchr<-dsmf[a]
b<-which(seqnames(dsmf)!="chrX")
Achr<-dsmf[b]

export(Xchr, "../Bigwiggle/sevinc/sm_Xchr_M_ucsc.bw", "bw")
export(Achr, "../Bigwiggle/sevinc/sm_Autosomes_M_ucsc.bw", "bw")
