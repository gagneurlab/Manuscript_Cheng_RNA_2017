library(ggplot2)
library(data.table)
library(ggrepel)
res_list <- readRDS("./data/codon_perturbation.rds")
dtc <- res_list$dt_result
dtc[, summary(abs_dy)]
#dtc[, y_ref_predict - y_perturbed_predict] %>%  abs %>% log10 %>% hist
#dtc[, y_ref_predict - y_perturbed_predict]  %>% hist(breaks = 100)
dtc[, dy := y_perturbed_predict - y_ref_predict ]
dtc <- dtc[order(-abs_dy)]
od <- dtc[, .(abs_dy = median(abs_dy)), by = .(from_codon, to_codon)][, .(abs_dy = mean(abs_dy)), by = .(from_codon)][order(-abs_dy)][,from_codon]
## sort by median
dt_codon <- dtc[, .(abs_dy = median(abs_dy)), by = .(from_codon, to_codon)]
dt_codon[, variable := from_codon]
dt_codon[, Features := 'codon']

## sTAI <- read.csv("./data/tAI_per_codon_SC_SP.csv")
## sTAI <- data.table(sTAI[,c('Codon', 'S_cerevisiae')])
## dt_codon <- merge(x=dt_codon, y=sTAI, by.x='to_codon', by.y='Codon')

## codon_opt <- fread('./data/classicalcodonoptimality.txt')
## # Change opt non-opt to colors
## codon_opt <- codon_opt[, .(codons, Scer)]
## codon_opt$Scer <- ifelse(codon_opt$Scer=='N', '#FF0000', '#0000FF')
## dt_codon <- merge(x=dt_codon, y=codon_opt, by.x='to_codon', by.y='codons')
## dt_codon$to_codon_col <- dt_codon$Scer
dt_codon$from_codon <- factor(dt_codon$from_codon, levels=od)

## codon_col <- codon_opt$Scer
## names(codon_col) <- codon_opt$codons
## codon_col <- codon_col[od]

cairo_pdf('./fig/Paper/Supplements/codon_snv_effect.pdf', width=6, height=8)
codon_single_base <- ggplot(aes(from_codon, exp(abs_dy)),
                            data = dt_codon) +
  #geom_label(aes(label = to_codon, fill = S_cerevisiae), size=2.5) +
  geom_point() +
  geom_text_repel(aes(label = to_codon), size=2.5) +
  #geom_jitter() +
  theme_bw() +
  labs(x = "codons", y = "Expected half-life fold-change") +
  theme(axis.title = element_text(size=rel(1.3))) +
  #theme(axis.text.y = element_text(color=codon_col)) +
  coord_flip()
print(codon_single_base)
dev.off()

saveRDS(codon_single_base, './fig/ggplot_obj/codon_single_base.rds')

## ggplot(dtc, aes(from_codon, log10(abs_dy))) +
##   geom_violin() +
##   coord_flip(expand=F, ylim=c(-4,0))

# from src/snv-perturbation-gc-length.R
res_list2 <- readRDS("./data/gc_and_length_perturbation.rds")
dtc2 <- res_list2$dt_result
dtc2[, abs_dy := abs(y_ref_predict - y_perturbed_predict)]
dtc2[, Features := variable]
## Drop out GC_content_CDS
dtc2 <- dtc2[!(dtc2$variable %in% c("GC_content_CDS", "GC_content_UTR3"))]
## Summarize with median
dtc2 <- dtc2[, .(abs_dy = median(abs_dy), variable), by = "Features"]
dtc2 <- unique(dtc2)

## ggplot(aes(variable, log10(abs_dy)),
##        data = dtc2
##        ) + geom_violin() +
##   coord_flip() +
##   labs(x = "Features")

## Get coefficients of start, stop and uAUG
coef_fit <- coef(res_list$lm_model)
startMinus3 <- data.table(coef_fit[grepl('start', names(coef_fit))])
stopPlus1 <- data.table(coef_fit[grepl('stop', names(coef_fit))])
uAUG <- data.table(coef_fit[grepl('uAUG', names(coef_fit))])
dtc3 <- rbindlist(list(startMinus3=startMinus3, stopPlus1=stopPlus1, uAUG=uAUG), idcol = 'variable')
dtc3[, abs_dy := abs(V1)]
dtc3[, Features := variable]

## Add motif single base effect
motifs <- c("TGTAAATA", "TGCAT", "ATATTC", "AAACAAA", "TTTTTTA")
motif_coef_data <- readRDS('./data/motif_coef_data.rds')
## Select only motif region
motif_coef_data[, keep := (positions > 2) & (positions<(motif_length+2)), by = .id]
motif_coef_data <- motif_coef_data[motif_coef_data[, keep]]
motif_coef_data[, Features := .id]
motif_coef_data[, variable := .id]
motif_coef_data[, abs_dy := abs(log(coef))]

## Combine codon with start stop uAUG
betas <- rbind(motif_coef_data[, .(variable, Features, abs_dy)], dtc2[, .(variable, Features, abs_dy)], dtc3[, .(variable, Features, abs_dy)], dt_codon[, .(variable, Features, abs_dy)])
## Only show points for those possible to show
betas[, show_points := ifelse(Features %in% c('codon','uAUG','startMinus3','stopPlus1', motifs), abs_dy, NA)]
## single point feature, we don't use jitter for them
single_point_features <- betas[, .N, by="Features"][N==1, Features]
betas[, single_points := ifelse(Features %in% single_point_features, abs_dy, NA)]
## multiple point features, use jitter
jitter_features <- betas[, .N, by="Features"][N>1, Features]
betas[, jitter_points := ifelse(Features %in% jitter_features, abs_dy, NA)]

# Make label publication quality
betas[, Features := ifelse(Features == "UTR3_length", "3\' UTR length",
                           ifelse(Features == "UTR5_length", "5\' UTR length",
                                  ifelse(Features == "GC_content_UTR5", "5\' UTR GC content", as.character(Features))))]

## Feature mean, order them by feature mean
betas[, Feature_mean := mean(abs_dy, na.rm=TRUE), by = Features]
betas$Features <- factor(betas$Features, levels = unique(betas[order(Feature_mean), Features]))


## single_snv_plot <- ggplot(betas, aes(Features, log10(abs_dy))) +
##   geom_violin() +
##   coord_flip() +
##   labs(x = "Features",
##        y = "Expected half-life log10(log(fold-change))") +
##   #geom_jitter(aes(Features, log10(show_points))) +
##   geom_jitter(aes(Features, log10(jitter_points))) +
##   geom_point(aes(Features, log10(single_points)), size=2) +
##   theme(axis.title=element_text(size=rel(1.3)))
## print(single_snv_plot)

# separate codons abs_dy
betas[, codon_abs_dy := ifelse(Features == 'codon', abs_dy, NA)]
betas[, abs_dy := ifelse(Features == 'codon', NA, abs_dy)]

cairo_pdf('./fig/Paper/Figure6_single_nucleotide.pdf', width=6, height=5)
single_snv_plot <- ggplot(betas, aes(Features, exp(abs_dy))) +
  geom_violin(fill = 'skyblue3', color = "skyblue3") +
  #geom_violin() +
  geom_violin(aes(Features, exp(codon_abs_dy)), fill = 'skyblue3', color = "skyblue3") +
  theme_bw() +
  #coord_flip() +
  labs(#x = "Features",
    title='Effect of single-nucleotide variations', x = '', y = "Expected half-life fold-change") +
  #geom_jitter(aes(Features, exp(show_points))) +
  geom_jitter(aes(Features, exp(jitter_points)), width=0.2, alpha=1) +
  geom_point(aes(Features, exp(single_points))) +
  theme(axis.title=element_text(size=rel(1.3))) +
  theme(axis.text.x = element_text(angle=45, hjust = 1)) +
  theme(plot.title = element_text(rel(1.3), face='bold', hjust = 0.5)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(single_snv_plot)
dev.off()

saveRDS(single_snv_plot, './fig/ggplot_obj/single_snv_plot.rds')
