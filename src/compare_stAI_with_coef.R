library(data.table)
library(cowplot)

## Correlation between codon coefficients with tAI
tAIs_coef_scer <- readRDS('./data/Scer_coef_tAI.rds')
tAIs_coef_scer <- tAIs_coef_scer[, .(Codon, beta_scer, log(S_cerevisiae))]
setnames(tAIs_coef_scer, c("Codon", "beta", "tAI"))

tAIs_coef_spom <- readRDS('./data/Spom_coef_tAI.rds')
tAIs_coef_spom <- tAIs_coef_spom[, .(Codon, beta_spom, log(S_pombe))]
setnames(tAIs_coef_spom, c("Codon", "beta", "tAI"))

tAIs_coef <- rbindlist(list(S.cerevisiae = tAIs_coef_scer, S.pombe = tAIs_coef_spom), idcol = T)

# Inf to NA
tAIs_coef_no_inf <- tAIs_coef[!is.infinite(tAI)]

correlation <- tAIs_coef_no_inf[, .(cor(beta, tAI, use = "pairwise.complete.obs"),
                             cor.test(beta, tAI)$p.value), by = '.id']
correlation[, x := c(-45)]
correlation[, y := c(1)]

corr_eqn <- function(corr_coef, digits = 2) {
  corr_coef <- round(corr_coef, digits)
  paste("italic(r) == ", corr_coef)
}
corr_eqn <- Vectorize(corr_eqn)

format_pval <- function(pval){
  pval <- format.pval(pval, digits = 2)
  paste("italic(P) == ", pval)
}
format_pval <- Vectorize(format_pval)

correlation[, correlation := corr_eqn(V1)]
correlation[, pval := format_pval(V2)]

pdf("./fig/Compare_sTAI_with_coef.pdf", width = 8, height = 4)
ggplot(tAIs_coef, aes(beta, exp(tAI))) +
  geom_text(aes(label = Codon), size = rel(2)) +
  facet_wrap(~.id) +
  labs(title = "Compare fitted codon coefficients with sTAI",
       x = "Fitted Coefficients",
       y = "sTAI") +
  scale_y_continuous(trans = 'log', breaks = c(0, 0.1, 0.2, 0.3, 0.5, 1)) +
  theme(strip.text = element_text(face = "italic")) +
  geom_text(data = correlation, aes(x=x, y=y, label = correlation), parse = T) +
  geom_text(data = correlation, aes(x=x, y=y-0.2, label = pval), parse = T)
dev.off()
