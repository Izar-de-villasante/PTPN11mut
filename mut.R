#library(Biostrings)
library(data.table)
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[2] = "gene.fasta"
} else if (length(args)==2){
  
}else stop("more than 2 args is not allowed. arg1 = prefix, arg2= gene fasta")


prefix <- args[1]
readCounts<- paste0(prefix,".readCounts")
codonCounts <- paste0(prefix,".codonCounts")
variantCounts <- paste0(prefix,".variantCounts")
var_counts<-fread(variantCounts,fill=T)

# reads

n_reads<-round(read.delim(readCounts,header = F)[1,2]/2,0)

# Wrong call variants
untrusted <- var_counts[ as.numeric(V4) > 3  | V6 > V4 |  V9 %like% "insdel" | V8 %like% "`fs*`" ,]
sum (untrusted$V1)
(sum (untrusted$V1)/sum(var_counts$V1) )*100

# Stop codons:
stopcodons <- var_counts[  endsWith(V8,"X"),]
(sum(stopcodons$V1) / sum(var_counts$V1) )*100

# Clean variants:
clean <- var_counts[!( as.numeric(V4) > 3  | V6 > V4 | V9 %like% "ins"| V9 %like% "del" | V8 %like% "`fs*`" ),]
codonlist <-tstrsplit(clean$V7,",") 

# wt gene seq:
gene_fasta <- args[2]
#gene <- Biostrings::readDNAStringSet(gene_fasta)
#wt_seq<-as.character(gene)
wt_seq<- paste0(readLines(gene_fasta)[-1],collapse="")

wt_codons<-sapply(1:594,function(x){
  end<- x*3     #start + 2
  start<- end-2      #1+(x-1)*3
  unname(wt_codon<-substring(wt_seq,start,end))
})

# Load dictionary:
ct<-readRDS("ct.RDS")
setDT(ct)
setkeyv(ct,"GeneticCode")
wt_aa <- ct[wt_codons,AminoAcids]

# Generate count table:
a <- lapply(codonlist,function(x){ stringi::stri_sub_all(x,c(1,-3),c(-9,-1))})
b<-lapply(a, function(x) cbind(clean$V1,data.table(do.call(rbind,x))))
ce<-rbindlist(b)
names(ce)<-c("count","pos","subs")
d<-ce[,sum(count),by=c("pos","subs")]
codons_clean<-dcast(d,as.numeric(pos)~subs,value.var = "V1",fun.aggregate = sum)
codons_clean<-codons_clean[-1,]#Remove first row that belongs to NA
codons_clean<-codons_clean[,.SD,.SDcols=intersect(names(codons_clean),ct$GeneticCode)]


# dict:
setkeyv(ct,"GeneticCode")
nct<-data.table(ct[,sapply(GeneticCode, function(x){
  cod1<-unlist(ct[x,1:3])
  apply(ct[,1:3],1,function(y){
    sum(cod1!=y)
    # print(cod1)
    # print(y)
    # print(sum(cod1!=y))
  })
})])

nct$GeneticCode<-colnames(nct)
nct[,AminoAcids:=ct[nct$GeneticCode,AminoAcids]]
dict<-melt(nct,id.vars = c("GeneticCode","AminoAcids"),variable.name = "codon_sub",value.name = "dist")
dict$aa_sub<-ct[as.character(dict$codon_sub),"AminoAcids"]

# Main:

codons_main<-codons_clean
codons_main$wt_codons <- wt_codons[-length(wt_codons)]
codons_main$wt_aa <- wt_aa[-length(wt_aa)]
codons_main[,pos:=1:.N]

codons_long<-melt(codons_main,
                  id.vars=c("pos","wt_codons","wt_aa"),
                  variable.name = "codon_sub",
                  value.name = "codon_counts"
)
setkey(ct,"AminoAcids")
# Amino Acids counts:
aa_main<-as.data.table(sapply(unique(ct$AminoAcids), function(aa){
  cols <- names(codons_clean) %in% ct[aa,GeneticCode]
  rowSums(codons_clean[,
                       names(codons_clean) %in% ct[aa,GeneticCode]
                       ,with=FALSE
  ]
  )
})
)
aa_main$wt_codons <- wt_codons[-length(wt_codons)]
aa_main$wt_aa <- wt_aa[-length(wt_aa)]
aa_main[,pos:=1:.N]

aa_long<-melt(aa_main,
              id.vars=c("pos","wt_codons","wt_aa"),
              variable.name = "aa_sub",
              value.name = "aa_counts"
)

setkeyv(codons_long,c("wt_aa","wt_codons","pos"))
setkeyv(aa_long,c("wt_aa","wt_codons","pos"))

## Merge:
long_codons_aa <- merge(codons_long,aa_long,allow.cartesian=T)
names(dict)[1:2]<-c("wt_codons","wt_aa")
main <- merge(dict,long_codons_aa,by=c("wt_aa","wt_codons","codon_sub","aa_sub"))
main[,N_reads:=n_reads]

#Overall codon transformation by number of consecutive bases changed:
library_eff<-main[,.(observed=sum(codon_counts>0),technical=length(codon_counts>0)),by=c("pos","dist","wt_codons","codon_counts")]
library_eff[,codon_eff:=round((observed/technical)*100,2),by=c("pos","dist")]
library_eff[,N_observed_c:=sum(observed),by=c("dist")]
library_eff[,N_technical_c:=sum(technical),by=c("dist")]
library_eff[,Total_codon_eff:=mean(codon_eff),by=c("dist")]


## Overall AminoAcid changes by number of consecutive bases changed:
aa_library_eff<-unique(main[,.SD,.SDcols=c("pos","wt_codons","wt_aa","aa_sub","aa_counts","dist")])
aa_library_eff[,(c("observed_aa","technical_aa")):=
                 .(observed_aa=sum(aa_counts>0),technical_aa=length(aa_counts>0)),
               by=c("pos","dist","wt_aa")]
aa_library_eff[,N_observed_aa:=sum(observed_aa),by=c("dist")]
aa_library_eff[,N_technical_aa:=sum(technical_aa),by=c("dist")]
aa_library_eff[,aa_eff:=round((observed_aa/technical_aa)*100,2),by=c("pos","dist")]
aa_library_eff[,Total_aa_eff:=mean(aa_eff),by=c("dist")]
setkeyv(library_eff,c("pos","wt_codons","dist"))
setkeyv(aa_library_eff,c("pos","wt_codons","dist"))
eff<-merge(library_eff,aa_library_eff,allow.cartesian=T)
eff$N_reads<-n_reads
unique(eff[,.SD,.SDcols=c("N_reads","dist","N_observed_c","N_technical_c","Total_codon_eff","N_observed_aa","N_technical_aa","Total_aa_eff")])

