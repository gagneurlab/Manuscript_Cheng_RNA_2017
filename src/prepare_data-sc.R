## prepare the csv script for modelling
#dt <- fread("./data-offline/jun-yeast/jun-complete.csv")
library(data.table)
dt <- readRDS('./data/Sun_mutation_UTRs.rds')
dt$num_uAUG <- dt$num_uAUG > 0
motifs <- c("TGTAAATA", "TGCAT", "ATATTC", "AAACAAA", "TTTTTTA")
#motifs <- c("TGTAAATA", "TGCAT", "ATATTC", "AAACAAA", "ATATATGC")
char_columns <- c("genename", "UTR3_seq", "UTR5_seq")
length <- c("UTR5_length", "UTR3_length", "CDS_length")
GC <- c("GC_content_UTR5", "GC_content_UTR3", "GC_content_CDS")
deltaG <- c("UTR5_deltaG", "UTR3_deltaG", "CDS_deltaG")
codons <- c("TTC", "TTA", "TTG", "TCT", "TCC", "TCA", "TCG",
            "TAT", "TAC", "TGT", "TGC", "TGG", "CTT", "CTC", "CTA", "CTG",
            "CCT", "CCC", "CCA", "CCG", "CAT", "CAC", "CAA", "CAG", "CGT",
            "CGC", "CGA", "CGG", "ATT", "ATC", "ATA", "ATG", "ACT", "ACC",
            "ACA", "ACG", "AAT", "AAC", "AAA", "AAG", "AGT", "AGC", "AGA",
            "AGG", "GTT", "GTC", "GTA", "GTG", "GCT", "GCC", "GCA", "GCG",
            "GAT", "GAC", "GAA", "GAG", "GGT", "GGC", "GGA", "GGG")
response <- "hlt.wt"
model_data <- dt[, c(length, response), with = F]
#model_data$Free.energy.of.ensemble <- -model_data$Free.energy.of.ensemble

not_log <- c(codons, "stopPlus1", "num_uAUG", "startMinus3", motifs, GC, deltaG)
## not_log <- c("stopPlus1", "num_uAUG", "startMinus3", motifs,
##              "GC_content_UTR5", "GC_content_UTR3")

model_data <- cbind(log(model_data), dt[, c(not_log, char_columns), with = F])
# normalize codon counts with CDS_length ---> #codon / CDS_length
for (cod in codons){
  model_data[, (cod) := get(cod) / exp(CDS_length)] # CDS_length was taken a log before
}
model_data <- model_data[complete.cases(model_data)]

model_data <- data.table(model.matrix(~., model_data[, -char_columns, with = F]),
                         model_data[, char_columns, with = F])
#model_data <- data.table(model.matrix(~., model_data[, -char_columns, with = F]))
model_data <- model_data[, -"(Intercept)", with = F]
stopifnot(all(!is.na(model_data)))
context_uAUG <- setdiff(names(model_data), c(motifs, GC, deltaG, "hlt.wt", length, codons, "genename",char_columns))

## save to disk
saveRDS(model_data, './data/cerevisiae_model_data_id.rds')

#stopifnot(length(motifs_pad) >1)
info <- list(motifs = motifs,
             length = length,
             context_uAUG = context_uAUG,
             codons = codons,
             GC = GC,
             deltaG = deltaG,
             response = "hlt.wt",
             data_path = "./data/cerevisiae_model_data_id.rds",
             id = "genename",
             n_folds = 10
             )
write(jsonlite::toJSON(info, pretty = TRUE, auto_unbox = TRUE), "./data/benchmark_metadata.json")
