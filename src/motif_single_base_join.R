## get all the information you need to fit the data
source('./src/utilities.R')
UTRs <- readRDS('./data/cerevisiae_model_data_id.rds')
UTRs <- data.frame(UTRs)
info <- jsonlite::fromJSON("./data/benchmark_metadata.json")

## Delete GC_context_CDS, because it is already explained by codon usage
##info$GC_context_uAUG <- info$GC_context_uAUG[info$GC_context_uAUG != "GC_content_CDS"]
## Delete not significant variants
# GC_content_UTR3, stopPlus1, startMinus3, AAACAAA
info$GC_context_uAUG <- c('num_uAUGTRUE', "GC_content_UTR5")
info$motifs <- info$motifs[info$motifs != "AAACAAA"]
#info$length <- info$length[info$length != "Free.energy.of.ensemble"]

info$GC_context_uAUG <- info$GC_context_uAUG[info$GC_context_uAUG != "GC_content_CDS"]
## Delete not significant variants
# GC_content_UTR3, stopPlus1, startMinus3, AAACAAA
info$GC_context_uAUG <- c('num_uAUGTRUE', "GC_content_UTR5")
info$motifs <- info$motifs[info$motifs != "AAACAAA"]
add <- c(info$length, info$motifs, info$GC_context_uAUG, info$codons)
# remove the motif we are fitting single base
motifs <- c("TGTAAATA", "TGCAT", "ATATTC",  "TTTTTTA") # "ATATATGC"

motif_single_base_effect_plot <- lapply(motifs, function(motif){
  fit <- lm_motif(motif, n=1, extend=0, UTRs, ratecol='hlt.wt', seqcol='UTR3_seq', add_features = add[add != motif],log=F)
  motif_single_base_effect <- plot.coef(fit=fit$fit,X=fit$X,main=paste(fit$X$consensus),rates = F,ylim=c(0.75,1.35))
})

# 5'UTR motif, but not significant in joint model
fit_AAACAAA <- lm_motif('AAACAAA', n=1, extend=0, UTRs, ratecol='hlt.wt', seqcol='UTR5_seq', add_features = add, log=F)
p_AAACAAA <- plot.coef(fit=fit_AAACAAA$fit,X=fit_AAACAAA$X,main=paste(fit_AAACAAA$X$consensus),rates = F,ylim=c(0.75,1.3))
saveRDS(p_AAACAAA, './data/5UTR_motif_coef_plot.rds')

names(motif_single_base_effect_plot) <- c('p_ATATTC', 'p_TGTAAATA', 'p_TGCAT', 'p_TTTTTTA')
p2 <- motif_single_base_effect_plot[[3]]
p1 <- motif_single_base_effect_plot[[1]]
p3 <- motif_single_base_effect_plot[[2]]
p4 <- motif_single_base_effect_plot[[4]]
#p5 <- motif_single_base_effect_plot[[5]]
motif_coef_data <- rbindlist(list(TGTAAATA=p1$data,ATATTC=p2$data,TGCAT=p3$data,TTTTTTA=p4$data),idcol=T)
motif_coef_data[, motif_length := nchar(as.character(.id))]
motif_coef_data$positions <- as.numeric(motif_coef_data$positions)
motif_coef_data[, pos_label := ifelse(positions>motif_length+2, paste0('+',positions-motif_length-2), # right outside of motif
                            ifelse(positions < 3, as.character(-3+positions), as.character(positions-2))), # left outside of motif
                            by = .(.id)]
motif_coef_data$pos_label <- factor(motif_coef_data$pos_label, c('-2','-1','1','2','3','4','5','6','7','8','+1','+2'))
motif_order <- c('ATATTC','TGCAT','TGTAAATA', 'TTTTTTA')
motif_coef_data$.id <- factor(motif_coef_data$.id, levels = motif_order)

saveRDS(motif_coef_data, './data/motif_coef_data.rds')
