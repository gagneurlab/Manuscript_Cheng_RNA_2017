##' ---
##' title: Motif validation qPCR results
##' author: Jun Cheng
##' output:
##'   html_document:
##'     pandoc_args: ["+RTS", "-K64m","-RTS"]
##'     toc: yes
##'---
##'
##+ set_workdir, echo=F
library(knitr)
library(rmarkdown)
library(data.table)
#library(tidyr)
library(ggplot2)
library(cowplot)
#library(plotly)
library(ggpval)
opts_knit$set(root.dir = getwd())
opts_chunk$set(echo=TRUE, cache=F, message=FALSE, warning=FALSE, fig.width = 12)
options(width=90)
##+

##' ## load data
#dt <- fread('./data/qPCR_motif_validation.csv')
dt <- fread('./data/170712_qPCR_tidy.csv')
dt <- melt(dt, measure.vars = paste0("rep", 1:9), value.name = "Ct", variable.name = "Bio_rep")

##' ## Normalize Ct by primer efficiency
primer_eff <- data.table(gene = c("sfg1", "tub2", "act1", "YLR093C"), primer_eff = c(0.97, 0.96, 0.98, 0.95))
primer_eff[, primer_eff := log2(1 + primer_eff)]
dt <- merge(dt, primer_eff, by = "gene")
dt[, Ct := Ct * primer_eff]
dt[, sample := paste(strain, Bio_rep, design, sep = '.')]
## y (i.e. -Ct) relate with log2(concentration) of target gene
dt[, y := -Ct]
dt[, sample := paste(strain, Bio_rep, design, sep = '.')]
dt[, act1 := gene == 'act1']
dt[, tub2 := gene == 'tub2']
dt[, sfg1 := gene == "sfg1"]
dt[, YLR093C := gene == "YLR093C"]

get_betas <- function(g = 'sfg1TRUE', beta = beta){
  wh_sfg1 <- grep(g, names(beta))
  # add act1 and tub2 control
  beta = data.table(beta=beta[wh_sfg1], condition=names(beta)[wh_sfg1])
  beta[, condition := substr(condition, nchar("sample")+1, nchar(condition)-nchar(g)-1)]
  beta[, Bio_rep := sapply(condition, function(x) strsplit(x, '.', fixed = T)[[1]][2])]
  beta[, design := sapply(condition, function(x) strsplit(x, '.', fixed = T)[[1]][3])]
  beta[, beta := 2^beta]
}

##' ## sfg1 gene
dt_sfg1 <- dt[strain == 'sfg1']
fit <- lm(y ~ sample + act1 + tub2 + sfg1:sample, data=dt_sfg1)
# The reference of this model is "samplesfg1.rep1.motif2", i.e. the sample size of the first replicate of sample motif2
betas_sfg1 <- get_betas(beta = coefficients(fit))
betas_sfg1 <- betas_sfg1[design %in% c('motif2', 'scram_motif2')]
p_sfg1 <- ggplot(betas_sfg1, aes(design, beta)) +
  geom_jitter(width=0) +
  ggtitle("sfg1")
p_sfg1 <- add_pval(p_sfg1,
                   list(c(1,2)),
                   heights = c(1.3),
                   log = F,
                   textsize = 5)
p_sfg1
#ggplotly(p_sfg1)

##' ## YLR093C gene
dt_YLR093C <- dt[strain == 'YLR093C']
fit <- lm(y ~ sample + act1 + tub2 + YLR093C:sample, data=dt_YLR093C)
betas_YLR093C <- get_betas(g = "YLR093CTRUE", beta = coefficients(fit))
betas_YLR093C <- betas_YLR093C[design %in% c('motif2', 'scram_motif2')]
p_YLR093C <- ggplot(betas_YLR093C, aes(design, beta)) +
  geom_jitter(width=0) +
  ggtitle("NYV1")
p_YLR093C <- add_pval(p_YLR093C,
                      list(c(1,2)),
                      heights = c(1.9),
                      log = F,
                      textsize = 5)
p_YLR093C
#ggplotly(p_YLR093C)

##' ## Plot fitted Betas
dt_2_test <- rbind(betas_sfg1, betas_YLR093C)
dt_2_test[, gene := tstrsplit(condition, "\\.")[[1]]]
dt_2_test[, design := ifelse(design == 'motif2', 'motifs', 'scrambled')]
dt_2_test[, gene := ifelse(gene == 'sfg1', 'SFG1', gene)]
dt_2_test[, gene := ifelse(gene == 'YLR093C', 'NYV1', gene)]
dt_2_test[, mu := mean(beta), by = c('design', 'gene')]
dt_2_test[, md := median(beta), by = c('design', 'gene')]
dt_2_test[, design := factor(design, levels = c('scrambled', 'motifs'))]
dt_2_test[, se := sd(beta) / sqrt(.N), by = c('design', 'gene')]

motif_effect <- ggplot(dt_2_test, aes(design, mu)) +
  geom_bar(stat = "identity", position = 'dodge', fill = 'skyblue3') +
  geom_errorbar(aes(ymax = mu + se, ymin=mu-se), width = .2) +
  geom_jitter(aes(design, beta), width=0.1) +
  theme(strip.text.x = element_text(face="bold.italic")) +
  facet_wrap(~ gene) +
  labs(title = "Effects of embedding 2 ATATTC motifs",
       x = NULL,
       y = 'Normalized expression level')

# Plot and add test pvalue
pdf('./fig/Paper/validation_motif2.pdf', width = 6, height = 4)
print(add_pval(motif_effect,
               barheight = 0.04,
               pairs = list(c(1,2)),
               #annotation = c('FC = ...'),
               textsize = 4,
               response = 'beta'))
dev.off()
saveRDS(motif_effect, './fig/ggplot_obj/validation_motif_effect.rds')

##' ## Join betas from two gene and test
# fit <- lm(log(beta) ~ gene + design, data = dt_2_test)
# pval <- summary(fit)$coefficients['motif2', "Pr(>|t|)"]
# pval
# resid_plot <- ggplot(data.frame(res = fit$residuals),
#                      aes(sample = res)) +
#               geom_point(stat = 'qq') +
#               labs(title = "Residual qqnorm plot")
#
# plot_grid(motif_effect, resid_plot)
#
# ##' ## Possible to do with one linear model?
# dt <- dt[design %in% c("motif2", "scram_motif2")]
# dt[, motif2 := ifelse(gene %in% c('sfg1', 'YLR093C') & design == 'motif2', 1, 0)]
# fit <- lm(y ~ sample + act1 + tub2 + sfg1 + YLR093C + motif2, data = dt)
# summary(fit)

##' ## Combine p value approach
# dt_2_test[gene=="sfg1"]
# pv_1 =dt_2_test[gene=="sfg1",wilcox.test(beta ~ design, alternative="less")$p.value]
# pv_2 =dt_2_test[gene=="YLR093C",wilcox.test(beta ~ design, alternative="less")$p.value]
#
# require(metap)
# pv = sumlog(c(pv_1,pv_2))


