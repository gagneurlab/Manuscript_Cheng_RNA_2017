###################################################################################################
# Figure 3
###################################################################################################
source('./src/Paper/ini.R')

codon_counts <- Scer[ ,c(codons, "CDS_length", 'hlt.wt'), with=F]
# normalize codon counts with CDS_length ---> #codon  / CDS_length
for (cod in codons){
  codon_counts[, (cod) := get(cod) / CDS_length]
}
codon_counts[, CDS_length := NULL]
# Delete TTT, use as reference. The beta for it will be 0
codon_counts[, TTT := NULL]
# log on response is essential, otherwise fit much worse
codon_counts[, hlt.wt := log(hlt.wt)]

fit <- lm(data = codon_counts, hlt.wt ~ .)
cv_rsquare2(fit, return_predict = FALSE)
# This explains 0.5504367 of the variance
predictions <- cv_rsquare2(fit, return_predict = TRUE)
predictions <- predictions[, .(hlt.wt = exp(y),
                                prediction = exp(y_hat))]

##' ## Translation related
cairo_pdf('./fig/Paper/Figure3_translation_features_ceveriviae.pdf',width=11, height=3.7)
plot.new()
pushViewport(viewport(width=1, height=1,layout=grid.layout(1, 11)))
pushViewport(viewport(layout.pos.col = c(1:4), layout.pos.row = 1))
par(fig = gridFIG(),mar=c(5,6,4,4),new = TRUE)
pvalheatscatter(predictions$prediction, predictions$hlt.wt,
                xlab="", ylab="", main='Codon-explained variance 55%',
                log='xy', cex=0.3, s_cor = F,
                xlim = c(2,200), ylim = c(2,200))
mtext(side=2,'Measured half-life [min]',line=2, font=2)
mtext(side=1,'Predicted half-life [min]',line=2, font=2)
grid.text(x=0.05,y=0.95,'A',gp=gpar(cex=1.5,font=2))
popViewport()
pushViewport(viewport(layout.pos.col = c(5:11), layout.pos.row = 1))
# codon_sig_knockout.rds come from src/sun_mutation_analysis.R
translation_features_sig_knockout <- readRDS('./fig/ggplot_obj/codon_sig_knockout.rds')
translation_features_sig_knockout <- translation_features_sig_knockout + theme(plot.margin=unit(c(1,0,0,0.3),"cm")) +
  guides(fill = guide_legend(keywidth = 0.7, keyheight = 0.8))
print(translation_features_sig_knockout,newpage=FALSE)
grid.text(x=0,y=0.95,'B', gp=gpar(cex=1.5,font=2))
popViewport(2)
dev.off()
