## compute the snv perturbation for the codons
library(dplyr)
library(data.table)
codons <- c("TTC", "TTA", "TTG", "TCT", "TCC", "TCA", "TCG",
            "TAT", "TAC", "TGT", "TGC", "TGG", "CTT", "CTC", "CTA", "CTG",
            "CCT", "CCC", "CCA", "CCG", "CAT", "CAC", "CAA", "CAG", "CGT",
            "CGC", "CGA", "CGG", "ATT", "ATC", "ATA", "ATG", "ACT", "ACC",
            "ACA", "ACG", "AAT", "AAC", "AAA", "AAG", "AGT", "AGC", "AGA",
            "AGG", "GTT", "GTC", "GTA", "GTG", "GCT", "GCC", "GCA", "GCG",
            "GAT", "GAC", "GAA", "GAG", "GGT", "GGC", "GGA", "GGG")
#stopifnot(length(codons) == 61)
## get all the information you need to fit the data
info <- jsonlite::fromJSON("./data/benchmark_metadata.json")
## Delete GC_context_CDS, because it is already explained by codon usage
##info$GC_context_uAUG <- info$GC_context_uAUG[info$GC_context_uAUG != "GC_content_CDS"]
## Delete not significant variants
# GC_content_UTR3, stopPlus1, startMinus3, AAACAAA
info$context_uAUG <- c('num_uAUGTRUE')
info$motifs <- info$motifs[info$motifs != "AAACAAA"]
info$GC <- c("GC_content_UTR5")
info$deltaG <- c("CDS_deltaG")
info$length <- info$length[info$length != "CDS_length"]

##' ## Fit the LM model with all features + motif countsand no sequence information
dt <- readRDS(info$data_path)
f <- paste0(info$response, " ~ ",
            paste0(info$length, collapse = "+"),
            " + ",
            paste0(info$motifs, collapse = "+"),
            " + ",
            paste0(info$context_uAUG, collapse = "+"),
            " + ",
            paste0(info$GC, collapse = "+"),
            " + ",
            paste0(info$deltaG, collapse = "+"),
            " + ",
            paste0(info$codons, collapse = "+")
)

fit <- lm(f, data = dt)
y_ref_predict <- predict(fit, newdata = dt)

############################################
## useful formulae:
##
## csv_length_nat = exp(CDS_length)
## codon_n = codon / CDS_length_nat
############################################


dt[, y_ref_predict := y_ref_predict]
#dt[, CDS_length_nat := exp(CDS_length) - 0.1]

###################################################################################################
## GC + length
## mutate
## log -> natural
#new_data <- copy(dt)
##
#new_data %>% names %>% grep("GC_content", . , value = T)
## mutate
vars <- list(c("GC_content_CDS", "CDS_length"),
             c("GC_content_UTR5", "UTR5_length"),
             c("GC_content_UTR3", "UTR3_length"))

dt_result <- lapply(vars, function(var_comb) {
  gc_var <- var_comb[1]
  length_var <- var_comb[2]
  new_data <- copy(dt)
  new_data2 <- copy(dt)
  ##
  new_data[, gc_var_mutate := get(gc_var) - 1/ (exp(get(length_var)))]
  new_data[, (gc_var) := gc_var_mutate]
  y_perturbed <- predict(fit, newdata = new_data)
  dtp_gc <- data.table(id = new_data[[info$id]],
                           y_measured = new_data[, get(info$response)],
                           y_ref_predict = new_data[, y_ref_predict],
                           y_perturbed_predict = y_perturbed,
                           variable = gc_var
                           )
  ##
  ## length
  if (!is.na(length_var)){
    new_data2[, length_nat_mutate := exp(get(length_var)) - 1]
    ## natural -> log
    new_data2[, (length_var) := log(length_nat_mutate)]
    y_perturbed <- predict(fit, newdata = new_data2)

    dtp_length <- data.table(id = new_data2[[info$id]],
                             y_measured = new_data2[, get(info$response)],
                             y_ref_predict = new_data2[, y_ref_predict],
                             y_perturbed_predict = y_perturbed,
                             variable = length_var
    )
    return(rbind(dtp_length, dtp_gc))
  }else{
    return(dtp_gc)
  }

}) %>% rbindlist

dt_result <- dt_result[variable != "CDS_length"] # CDS_length was no longer significant

###################################################################################################

results <- list(
  info = info,
  lm_model = fit,
  dt_result = dt_result
)

saveRDS(results, "./data/gc_and_length_perturbation.rds")
