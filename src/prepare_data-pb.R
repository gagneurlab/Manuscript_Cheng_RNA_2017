## prepare the csv script for modelling
#dt <- fread("./data-offline/jun-yeast/jun-complete.csv")
library(data.table)
dt <- readRDS('./data/TUs_Pombe.rds')
dt <- plyr::rename(dt, c("Txn.unit" = "genename",
                         "UTR5.length" = "UTR5_length",
                         "UTR3.length" = "UTR3_length",
                         "cds_seqs" = "CDS_seq",
                         "PositionM3" = "startMinus3"))

dt$num_uAUG <- dt$num_uAUG > 0
motifs <- c( "CAACCA",
            "ACCAAC",
            #"AACCAC",
            #"CCAACA",
            "TATTTAT",
            "TTAATGA",
            "ACTAAT")
char_columns <- c("genename", "UTR3_seq", "UTR5_seq")
length <- c("UTR5_length", "UTR3_length", 
            "CDS_length")
codons <- c( "TTT", "TTC", "TTA", "TTG", "TCT", "TCC", "TCA", "TCG", 
            "TAT", "TAC", "TGT", "TGC", "TGG", "CTT", "CTC", "CTA", "CTG", 
            "CCT", "CCC", "CCA", "CCG", "CAT", "CAC", "CAA", "CAG", "CGT", 
            "CGC", "CGA", "CGG", "ATT", "ATC", "ATA", "ATG", "ACT", "ACC", 
            "ACA", "ACG", "AAT", "AAC", "AAA", "AAG", "AGT", "AGC", "AGA", 
            "AGG", "GTT", "GTC", "GTA", "GTG", "GCT", "GCC", "GCA", "GCG", 
            "GAT", "GAC", "GAA", "GAG", "GGT", "GGC", "GGA", "GGG")
response <- "half.life"

model_data <- dt[, c(length, codons, response), with = F]

not_log <- c("stopPlus1", "num_uAUG", "startMinus3", motifs,
             "GC_content_UTR5", "GC_content_UTR3", "GC_content_CDS",
             "UTR5_deltaG", "UTR3_deltaG", "CDS_deltaG")

model_data <- cbind(log(model_data), dt[, c(not_log, char_columns), with = F])
#model_data <- subset(model_data, half.life>0)
model_data <- model_data[complete.cases(model_data)]                    

model_data <- data.table(model.matrix(~., model_data[, -char_columns, with = F]),
                          model_data[, char_columns, with = F])
#model_data <- data.table(model.matrix(~., model_data[, -char_columns, with = F]))
model_data <- model_data[, -"(Intercept)", with = F]
stopifnot(all(!is.na(model_data)))
GC_context_uAUG <- setdiff(names(model_data), c(motifs, "half.life", length_energy, codons, "genename",char_columns))

## save to disk
saveRDS(model_data, './data/pombe_model_data_id.rds')

#stopifnot(length(motifs_pad) >1)
info <- list(motifs = motifs,
             length = length,
             GC_context_uAUG = GC_context_uAUG,
             codons = codons,
             response = "half.life",
             data_path = "./data/pombe_model_data_id.rds",
             id = "genename",
             n_folds = 10
             )
write(jsonlite::toJSON(info, pretty = TRUE, auto_unbox = TRUE), "./data/benchmark_metadata_pombe.json")

