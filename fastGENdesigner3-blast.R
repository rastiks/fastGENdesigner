# INSTALLATION
# install.packages('rBLAST', repos = 'https://mhahsler.r-universe.dev')
# download.file("https://ftp.ncbi.nlm.nih.gov/blast/db/human_genome.00.tar.gz","human_genome.00.tar.gz",  mode='wb')
# download.file("https://ftp.ncbi.nlm.nih.gov/blast/db/human_genome.01.tar.gz","human_genome.01.tar.gz",  mode='wb')
# untar("human_genome.00.tar.gz", exdir="human_genome")
# untar("human_genome.01.tar.gz", exdir="human_genome")

library('rBLAST')

# DB CONFIGURATION!!!!
bl <-blast(db="/home/ppola/ncbi-blast-2.14.0+/blast/db/GCF_000001405.39_top_level")
seq_for_blast<-readDNAStringSet("primers.fasta")
cl <- predict(bl, seq_for_blast, BLAST_args = '-task blastn-short')
write.csv(cl, file= "blast.out.csv")

#Treba spravit skript aby identifikoval fastGEN designer off-targety (v tomto priklade offtarget je 596 a 678) ??????