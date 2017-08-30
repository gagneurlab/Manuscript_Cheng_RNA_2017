## compute the snv perturbation for the codons
## get the codon counts:
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
y_ref_predict <- predict(fit)

## SNV occurs in the cds -> we change the codon
## hence, I have to sample through all the codon combinations

library(stringdist)
## Synonymous codons
degenerate_codons <- list(c("GCT", "GCC", "GCA", "GCG"),
                          c("CGT", "CGC", "CGA", "CGG", "AGA", "AGG"),
                          c("AAT", "AAC"),
                          c("GAT", "GAC"),
                          c("TGT", "TGC"),
                          c("CAA", "CAG"),
                          c("GAA", "GAG"),
                          c("GGT", "GGC", "GGA", "GGG"),
                          c("CAT", "CAC"),
                          c("ATT", "ATC", "ATA"),
                          c("TTA", "TTG", "CTT", "CTC", "CTA", "CTG"),
                          c("AAA", "AAG"),
                          #c("TTT", "TTC"),
                          c("CCT", "CCC", "CCA", "CCG"),
                          c("TCT", "TCC", "TCA", "TCG", "AGT", "AGC"),
                          c("ACT", "ACC", "ACA", "ACG"),
                          c("TAT", "TAC"))
degenerate_pairs <- lapply(degenerate_codons, function(i) expand.grid(from_codon = i,
                                                                      to_codon = i)) %>% rbindlist
degenerate_pairs <- degenerate_pairs[stringdist(from_codon, to_codon, method = "hamming") == 1]
codon_combinations <- expand.grid(from_codon = codons,
                                  to_codon = codons,
                                  stringsAsFactors = FALSE)  %>% as.data.table
codon_combinations <- codon_combinations[stringdist(from_codon, to_codon, method = "hamming") == 1]

## remove non-synonymous?
## codon_combinations <- rbind(codon_combinations, degenerate_pairs)
## codon_combinations <- codon_combinations[(duplicated(codon_combinations) | duplicated(codon_combinations, fromLast = TRUE))]
## Only consider Synonymous
codon_combinations <- degenerate_pairs

codon_combinations <- codon_combinations  %>%
  t %>% as.data.table(stringsAsFactors = F) %>% as.list

############################################
## useful formulae:
##
## csv_length_nat = exp(cds_length)
## codon_n = codon*CDS_length_nat
############################################

dt[, CDS_length_nat := exp(CDS_length)]
dt[, y_ref_predict := y_ref_predict]

# epsilon
eps <- 1e-8

dt_result <- lapply(codon_combinations, function(codon_pair) {
  from_codon <- codon_pair[1]
  to_codon <- codon_pair[2]

  ## if I change from one codon to another,
  ## what is the effect on half-life for all the genes?

  ## 1. convert codon variables into actual codon counts
  new_data <- copy(dt)
  new_data[, from_codon_n := get(from_codon) * CDS_length_nat]
  new_data[, to_codon_n := get(to_codon)  * CDS_length_nat]

  # adjust for machine epsilon
  new_data[, from_codon_n := ifelse(abs(from_codon_n) > eps, from_codon_n, 0)]
  new_data[, to_codon_n := ifelse(abs(to_codon_n) > eps, to_codon_n, 0)]

  ## mutate:
  new_data[, from_codon_n := ifelse(from_codon_n, from_codon_n - 1, 0)]
  new_data[, to_codon_n := ifelse(to_codon_n, to_codon_n + 1, 0)]
  # adjust for machine epsilon again
  new_data[, from_codon_n := ifelse(abs(from_codon_n) > eps, from_codon_n, 0)]
  new_data[, to_codon_n := ifelse(abs(to_codon_n) > eps, to_codon_n, 0)]

  ## consider only valid mutations
  #new_data <- new_data[from_codon_n >= 0 & to_codon_n >= 0]
  stopifnot(dim(new_data[from_codon_n >= 0 & to_codon_n >= 0]) == dim(new_data))

  ## convert actual codon counts into codon variables
  new_data[, (from_codon) := from_codon_n / CDS_length_nat]
  new_data[, (to_codon) := to_codon_n/ CDS_length_nat]

  y_perturbed <- predict(fit, newdata = new_data)

  data.table(id = new_data[[info$id]],
             y_measured = new_data[, get(info$response)],
             y_ref_predict = new_data[, y_ref_predict],
             y_perturbed_predict = y_perturbed,
             from_codon = from_codon,
             to_codon = to_codon
             )
}) %>% rbindlist

setnames(dt_result , "id", info$id)
dt_result[, abs_dy := abs(y_ref_predict - y_perturbed_predict)]

results <- list(
  info = info,
  lm_model = fit,
  dt_result = dt_result
)

saveRDS(results, "./data/codon_perturbation.rds")
