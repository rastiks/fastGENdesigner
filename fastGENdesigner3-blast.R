**************************************************************** 
Instalace, download.file sa nejak sekal, asi cez wget to bolo myslim OK  
***************************************************************

install.packages('rBLAST', repos = 'https://mhahsler.r-universe.dev')

download.file("https://ftp.ncbi.nlm.nih.gov/blast/db/human_genome.00.tar.gz","human_genome.00.tar.gz",  mode='wb')
download.file("https://ftp.ncbi.nlm.nih.gov/blast/db/human_genome.01.tar.gz","human_genome.01.tar.gz",  mode='wb')
untar("human_genome.00.tar.gz", exdir="human_genome")
untar("human_genome.01.tar.gz", exdir="human_genome")


**********************************************************************************
setwd("~/primer3/src") 
library('rBLAST')
bl <-blast(db="./human_genome/GCF_000001405.39_top_level")
seq_for_blast<-readDNAStringSet("poolfile1-final.txt.fasta")
cl <- predict(bl, seq_for_blast, BLAST_args = '-task blastn-short')
View(cl)
write.csv(cl, file= "blast.out.csv")

****************************************************************************************
> cl
                      QueryID      SubjectID Perc.Ident Alignment.Length Mismatches Gap.Openings Q.start Q.end   S.start
1  LRG_310t1e21_1004_1025_1_R   NC_000003.12        100               23          0            0       1    23 179234258
2  LRG_310t1e21_1004_1025_1_R   NC_000003.12         95               20          1            0       1    20 168682636
3  LRG_310t1e21_1004_1025_1_R   NC_000003.12        100               16          0            0       3    18 175965820
4  LRG_310t1e21_1004_1025_1_R   NC_000011.10        100               20          0            0       1    20  55010619
5  LRG_310t1e21_1004_1025_1_R   NC_000011.10        100               16          0            0       2    17 134731360

Treba spravit skript aby identifikoval fastGEN designer off-targety (v tomto priklade offtarget je 596 a 678)
 
	QueryID	SubjectID	Perc.Ident	Alignment.Length	Mismatches	Gap.Openings	Q.start	Q.end	S.start	S.end	E	Bits
592	LRG_310t1e10_539_546_4_F	NC_000003.12	100	21	0	0	1	21	179218240	179218260	0.003	42.1
596	LRG_310t1e10_539_546_4_F	NC_000022.11	100	21	0	0	1	21	16572058	16572078	0.003	42.1
675	LRG_310t1e10_539_546_4_R	NC_000003.12	100	20	0	0	1	20	179218353	179218334	0.008	40.1
678	LRG_310t1e10_539_546_4_R	NC_000022.11	100	20	0	0	1	20	16572170	16572151	0.008	40.1
