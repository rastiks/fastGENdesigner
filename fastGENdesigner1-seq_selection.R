
BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
setwd("~/primer3/src") # set Working DIRECTORY where u have input.txt!)
#####################################################333
library(EnsDb.Hsapiens.v86)
library(ensembldb)
library("BSgenome.Hsapiens.UCSC.hg38")

edb <- EnsDb.Hsapiens.v86

protein <- proteins(edb, filter = ~ genename == "PIK3CA")
protein #hladanie spravneho transkriptu / produktu transkriptu
protein@listData[["tx_id"]]

#nacitanie ROI
d <- read.delim("input.txt", header=T)
prt <- IRanges(start=d$start, end=d$stop, names=d$gene) 
gnm <- proteinToGenome(prt, edb)
gr_total = NULL
gr_total = gnm[[1]]
for(i in 2:length(gnm)) {
  gr_total<-c(gr_total,gnm[[i]]) 
} 
gr_total
   

#hladanie region ROI +-200bp
chrom = paste("chr", gr_total@seqnames,sep="")
start=gr_total@ranges@start-200
middle = gr_total@ranges@start
width<-gr_total@ranges@width
stop = gr_total@ranges@start+width+200

#ziskavanie sekvencii - seq_complete
seqA<-getSeq(Hsapiens,chrom,start = start, end = middle-1)
seqA_ch<-getSeq(Hsapiens,chrom,start = start, end = middle-1, as.character = TRUE)
seqB<-getSeq(Hsapiens,chrom,start = middle, end = middle+width-1)
seqC<-getSeq(Hsapiens,chrom,start = stop-200, end = stop)
seqC_ch<-getSeq(Hsapiens,chrom,start = stop-200, end = stop, as.character = TRUE)
seq_complete<-seqA<-getSeq(Hsapiens,chrom,start = start, end = stop)

# ulozenie bedu ROI 
df <- data.frame(seqnames=seqnames(gr_total),
  starts=start(gr_total)-1,
  ends=end(gr_total),
  names=paste(gr_total$exon_id,gr_total$protein_start,gr_total$protein_end,sep='_'),
  scores=c(rep(".", length(gr_total))),
  strands=strand(gr_total))

write.table(df, file="ROI.bed", quote=F, sep="\t", row.names=F, col.names=F)

df_seq <- data.frame(seqnames=seqnames(gr_total),
  starts=start,
  ends=stop,
  names=paste(gr_total$exon_id,gr_total$protein_start,gr_total$protein_end,sep='_'),
  scores=c(rep(".", length(gr_total))),
  strands=strand(gr_total))
 
write.table(df_seq, file="full_seq_ROI.bed", quote=F, sep="\t", row.names=F, col.names=F)

# ulozenie seq_Complete
library(seqinr)
file.remove('sequences.fasta')
write.fasta(sequences = as.list(seq_complete), names = seq(1, 12, 1),file.out = 'sequences.fasta',open = "a")

gr_total
