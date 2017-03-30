library(stringr)
library(dplyr)
# Convert GENCODE transcripts into sequences ------------------------------

gencode = read.csv("./Data/All_gencode_transcripts.csv")

lncRNAs = unlist(lapply(str_split(gencode$ID, "\\."), "[[", 1))
gencode$ID = lncRNAs


dups = gencode[duplicated(gencode$ID), ]
dups = dups[order(dups$ID), ]

df = gencode %>% 
    group_by(ID) %>% 
    summarise(RPKM1.RPKM2 = mean(X.RPKM1.RPKM2.))

library(biomaRt)
mart = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
# map = getBM(mart = mart, attributes = c("ensembl_gene_id","percentage_gc_content","gene_biotype","start_position","end_position","chromosome_name", "external_gene_name","transcript_length"),
#             filters = "ensembl_gene_id", values=rownames(res05))

seq = getSequence(id = df$ID, type="ensembl_transcript_id", seqType = "cdna", mart = mart)


df = cbind(df, seq$cdna[match(df$ID, seq$ensembl_transcript_id)] )
df = df[!is.na(df$`seq$cdna[match(df$ID, seq$ensembl_transcript_id)]`), ]

df = df[order(df$RPKM1.RPKM2, decreasing = TRUE), ]


write.csv(df, "./Data/gencode_seqs.csv")
