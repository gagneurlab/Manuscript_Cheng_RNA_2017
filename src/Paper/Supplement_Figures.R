###################################################################################################
# All supplementary figures
###################################################################################################
# Load data and functions
source('./src/Paper/ini.R')

##' ### Length
cairo_pdf('./fig/Paper/Supplements/Fig_S1_Length.pdf', width=6, height=9)
par(mfrow=c(3,2))
par(mgp = c(2.2, 1, 0))
##' #### UTR5 length
##' Long 5' prime UTR length can increase the dificulties of ribosome to find the correct start codon. We also need more literatuer search
pvalheatscatter(Scer$UTR5_length, Scer$hlt.wt, s_cor = T,ylab='half-life [min]',xlab='5\'UTR length', main='S. cerevisiae',log='xy')
label_plot(2.5, 800, 'A', cex=2, font=2)
pvalheatscatter(Pombe$UTR5.length, Pombe$half.life, s_cor=T,ylab='half-life [min]',xlab='5\'UTR length', main='S. pombe',log='xy',ylim=c(1,200))
label_plot(2.5, 600, 'B', cex=2, font=2)
##' #### CDS length
pvalheatscatter(Scer$CDS_length, Scer$hlt.wt, s_cor=T,ylab='half-life [min]',xlab='CDS length', main='S. cerevisiae',log='xy')
label_plot(40, 800, 'C', cex=2, font=2)
pvalheatscatter(Pombe$CDS_length, Pombe$half.life, s_cor=T,ylab='half-life [min]',xlab='CDS length', main='S. pombe',log='xy',ylim=c(1,200))
label_plot(40, 600, 'D', cex=2, font=2)
##' **Note:** 3' UTR length also positively correlate with GC content
pvalheatscatter(Scer$UTR3_length, Scer$hlt.wt, s_cor=T,ylab='half-life [min]',xlab='3\'UTR length', main='S. cerevisiae',log='xy')
label_plot(0.3, 800, 'E', cex=2, font=2)
pvalheatscatter(Pombe$UTR3.length, Pombe$half.life, s_cor=T,ylab='half-life [min]',xlab='3\'UTR length', main='S. pombe',log='xy',ylim=c(1,200))
label_plot(2, 600, 'F', cex=2, font=2)
dev.off()

##' ### Folding energy
cairo_pdf('./fig/Paper/Supplements/Fig_S1_Free_energy.pdf', width=6, height=9)
par(mfrow=c(3,2))
par(mgp = c(2.2, 1, 0))
deltaG_label <- c(-500, -100, -10, -1)
##' #### Folding
pvalheatscatter(-log(-Scer$UTR5_deltaG), Scer$hlt.wt, s_cor=T,ylab='half-life [min]',xlab='5\' Free Energy', main='S. cerevisiae', log='y',xaxt='n')
axis(1, at = -log(-deltaG_label), labels=deltaG_label,cex.axis=1 )
label_plot(-7.5, 550, 'A', cex=2, font=2)
pvalheatscatter(-log(-Pombe$UTR5_deltaG), Pombe$half.life, s_cor=T,ylab='half-life [min]',xlab='5\' Free Energy', main='S. pombe',log='y',ylim=c(1,200),xaxt='n')
axis(1, at = -log(-deltaG_label), labels=deltaG_label,cex.axis=1 )
label_plot(-7.5, 600, 'B', cex=2, font=2)
pvalheatscatter(-log(-Scer$CDS_deltaG), Scer$hlt.wt, s_cor=T,ylab='half-life [min]',xlab='CDS Free Energy', main='S. cerevisiae', log='y',xaxt='n')
axis(1, at = -log(-deltaG_label), labels=deltaG_label,cex.axis=1 )
label_plot(-7.5, 800, 'C', cex=2, font=2)
pvalheatscatter(-log(-Pombe$CDS_deltaG), Pombe$half.life, s_cor=T,ylab='half-life [min]',xlab='CDS Free Energy', main='S. pombe',log='y',ylim=c(1,200),xaxt='n')
axis(1, at = -log(-deltaG_label), labels=deltaG_label,cex.axis=1 )
label_plot(-7.5, 600, 'D', cex=2, font=2)
pvalheatscatter(-log(-Scer$UTR3_deltaG), Scer$hlt.wt, s_cor=T,ylab='half-life [min]',xlab='3\' Free Energy', main='S. cerevisiae', log='y',xaxt='n')
axis(1, at = -log(-deltaG_label), labels=deltaG_label,cex.axis=1 )
label_plot(-7.5, 800, 'E', cex=2, font=2)
pvalheatscatter(-log(-Pombe$UTR3_deltaG), Pombe$half.life, s_cor=T,ylab='half-life [min]',xlab='3\' Free Energy', main='S. pombe',log='y',ylim=c(1,200),xaxt='n')
axis(1, at = -log(-deltaG_label), labels=deltaG_label,cex.axis=1 )
label_plot(-7.5, 600, 'F', cex=2, font=2)
dev.off()

##' ### GC content
cairo_pdf('./fig/Paper/Supplements/Fig_S2_GC_content.pdf',height=9,width=6)
par(mfrow=c(3,2))
par(mgp = c(2.2, 1, 0))
pvalheatscatter(Scer$GC_content_UTR5, Scer$hlt.wt, s_cor=T,ylab='half-life [min]',xlab='GC content 5\'UTR', main='S. cerevisiae',log='xy',xlim=c(0.1,0.6))
label_plot(0.08, 800, 'A', cex=2, font=2)
pvalheatscatter(Pombe$GC_content_UTR5, Pombe$half.life, s_cor=T,ylab='half-life [min]',xlab='GC content 5\'UTR', main='S. pombe',log='xy',ylim=c(1,200), xlim=c(0.1, 0.7))
label_plot(0.08, 600, 'B', cex=2, font=2)
pvalheatscatter(Scer$GC_content_CDS, Scer$hlt.wt, s_cor=T,ylab='half-life [min]',xlab='GC content CDS', main='S. cerevisiae',log='xy')
label_plot(0.25, 800, 'C', cex=2, font=2)
pvalheatscatter(Pombe$GC_content_CDS, Pombe$half.life, s_cor=T,ylab='half-life [min]',xlab='GC content CDS', main='S. pombe',log='xy',ylim=c(1,200), xlim=c(0.3,0.55))
label_plot(0.27, 600, 'D', cex=2, font=2)
pvalheatscatter(Scer$GC_content_UTR3, Scer$hlt.wt, s_cor=T,ylab='half-life [min]',xlab='GC content 3\'UTR', main='S. cerevisiae',log='xy')
label_plot(0.08, 800, 'E', cex=2, font=2)
pvalheatscatter(Pombe$GC_content_UTR3, Pombe$half.life, s_cor=T,ylab='half-life [min]',xlab='GC content 3\'UTR', main='S. pombe',log='xy',xlim=c(0.1,0.45),ylim=c(1,200))
label_plot(0.08, 600, 'F', cex=2, font=2)
dev.off()


## S.cer start stop
start_scer <- boxplotOmmitN(Scer[ ,hlt.wt],Scer[ ,startMinus3],breaks=y_breaks,
                    main='S. cerevisiae',xlab='start codon\\-3',ylim=c(1,300),give.n.size=3,n.vjust=0,axis.title.size=1.2,plot.title.size=1.2,axis.text.size=1,fill='white')
start_scer <- start_scer + scale_fill_identity() +
  #theme(plot.margin = unit(c(0,1,0,0), "cm")) +
  theme(plot.title = element_text(face='bold.italic', hjust = 0.5))
start_scer <- add_pval_ggplot(start_scer,pairs = list(c(1,2),c(1,3),c(1,4)),heights = c(100,160,240),size=3,pval_text_adj=0.15)
print(start_scer,newpage=FALSE)
grid.text(x=0.05,y=0.95,'A',gp=gpar(cex=1.5,font=2))

stop_scer <- ggplot(Scer, aes(stopPlus1, hlt.wt, fill=Stop_codon)) +
  geom_boxplot(outlier.size = 0.5) +
  scale_y_continuous('Half-life [min]', trans='log10', breaks=y_breaks) +
  theme_bw() +
  theme(axis.title = element_text(size=rel(1.2), face='bold')) +
  theme(axis.text.x = element_text(size=rel(1.2), angle=90, face='bold')) +
  theme(axis.text.y = element_text(size=rel(1.2), face='bold')) +
  theme(plot.title = element_text(size=rel(1.2), face = 'bold.italic')) +
  labs(x='stop codon\\+1', title='S. cerevisiae') +
  stat_summary(fun.data = give.n, geom = "text", fun.y = median, position = position_dodge(width = 0.75),size=rel(3))
stop_scer <- add_pval_ggplot(stop_scer,pairs = list(c(1,4),c(5,8),c(9,12)),heights = c(240,240,240),size=3,pval_text_adj=0.15) +
  theme(plot.title = element_text(face='bold.italic', hjust = 0.5))


##' ## Figure Supplement 4
##' ## Plot start stops for pombe
start_spom <- boxplotOmmitN(Pombe[ ,half.life],Pombe[ ,PositionM3],breaks=y_breaks,main='S. Pombe',xlab='start codon\\+3',ylim=c(1,250), give.n.size=3,n.vjust=0,axis.title.size=1.2,plot.title.size=1.2,axis.text.size=1, fill='white')
start_spom <- start_spom + scale_fill_identity() + theme(plot.title = element_text(face='bold.italic', hjust = 0.5))
start_spom <- add_pval_ggplot(start_spom,pairs = list(c(1,2),c(1,3),c(1,4)),heights = c(100,140,200),pval_text_adj=0.15,size=3)

stop_spom <- ggplot(Pombe, aes(stopPlus1, half.life, fill=Stop_codon)) +
  geom_boxplot(outlier.size = 0.5) +
  scale_y_continuous('Half-life [min]', trans='log10', breaks=y_breaks) +
  theme_bw() +
  theme(axis.title = element_text(size=rel(1.2), face='bold')) +
  theme(axis.text.x = element_text(size=rel(1.2), angle=90, face='bold')) +
  theme(axis.text.y = element_text(size=rel(1.2), face='bold')) +
  theme(plot.title = element_text(size=rel(1.2), face = 'bold.italic')) +
  labs(x='stop codon\\+1', title='S. pombe') +
  stat_summary(fun.data = give.n, geom = "text", fun.y = median, position = position_dodge(width = 0.75),size=rel(3))
stop_spom <- add_pval_ggplot(stop_spom, pairs = list(c(1,4),c(5,8),c(9,12)),heights = c(240,240,240),size=3,pval_text_adj=0.15) +
  theme(plot.title = element_text(face='bold.italic', hjust = 0.5))


## Calculate seqlogo matrix
Cerevisiae_start <- consensusMatrix(start_context_seq(Scer,UTR5seqcol = 'UTR5_seq',CDSseqcol = 'CDS_seq',l=6,r=6))
Cerevisiae_start <- Cerevisiae_start/colSums(Cerevisiae_start)
Pombe_start <- consensusMatrix(start_context_seq(Pombe,UTR5seqcol = 'UTR5_seq',CDSseqcol = 'cds_seqs',l=6,r=6))
Pombe_start <- Pombe_start/colSums(Pombe_start)

## Pombe tAI
## Calculate tAI
tAI_codons <- fread('./data/tAI_per_codon_SC_SP.csv')
codons = as.character(tAI_codons$Codon)
tAIs = tAI_codons$S_pombe
names(tAIs) = codons

cod2tAI = function(cod,tAIs){
  return(as.numeric(tAIs[cod]))
}
# function to remove stop codon from CDS sequence
remove_stop <- function(s){
  n <- nchar(s)
  substr(s, 1, n - 3)
}
Pombe[, CDS_seq := sapply(cds_seqs, remove_stop)]
tAI_Spom = tAImatrix(CDSseq=Pombe$CDS_seq,tAIs=tAIs,percentage = 10)
Pombe$tAI <- tAI_Spom$avg_gene_tAI


pdf('./fig/Paper/Supplements/Sup_Pombe_start_stop1.pdf', width=10, height=8)
plot.new()
pushViewport(viewport(width=1, height=1,layout=grid.layout(2, 7)))
pushViewport(viewport(layout.pos.col = c(1:4), layout.pos.row = 1))
seqLogo(Cerevisiae_start,xaxislabel = c(-6:0,0,0,1:6), fold = 3)
grid.text(x=0.05,y=0.95,'A',gp=gpar(cex=1.5,font=2))
popViewport()
pushViewport(viewport(layout.pos.col = c(1:4), layout.pos.row = c(2)))
seqLogo(Pombe_start,xaxislabel = c(-6:0,0,0,1:6),fold = 3)
grid.text(x=0.05,y=0.95,'C',gp=gpar(cex=1.5,font=2))
popViewport()
pushViewport(viewport(layout.pos.col = c(5:7), layout.pos.row = 1))
print(start_scer, newpage=FALSE)
grid.text(x=0.05,y=0.95,'B',gp=gpar(cex=1.5,font=2))
popViewport()
pushViewport(viewport(layout.pos.col = c(5:7), layout.pos.row = 2))
print(start_spom, newpage=FALSE)
grid.text(x=0.05,y=0.95,'D',gp=gpar(cex=1.5,font=2))
popViewport(2)
dev.off()

pdf('./fig/Paper/Supplements/Sup_Pombe_start_stop2.pdf', width=11, height=4)
plot.new()
pushViewport(viewport(width=1, height=1, layout=grid.layout(1, 8)))
pushViewport(viewport(layout.pos.col = c(1:4), layout.pos.row = 1))
print(stop_scer, newpage = FALSE)
grid.text(x=0.05,y=0.95,'E',gp=gpar(cex=1.5,font=2))
popViewport()
pushViewport(viewport(layout.pos.col = c(5:8), layout.pos.row = 1))
print(stop_spom, newpage = FALSE)
grid.text(x=0.05,y=0.95,'F',gp=gpar(cex=1.5,font=2))
popViewport(2)
dev.off()


##' ## Compare sTAI with codon regression coeficients
# Please refer to script src/compare_stAI_with_coef.R

## tAI half-life
pdf('./fig/Paper/Supplements/Sup_tAI.pdf', width=10, height=5)
par(mfrow = c(1, 2))
Scer_cor <- cor(log(Scer$hlt.wt), log(Scer$tAI),m='s',use='pairwise.complete.obs')
pvalheatscatter(Scer$tAI, Scer$hlt.wt, xlab="", ylab="", main='', xlim=c(0.35,0.7), log='xy', cex=0.3)
mtext(side=2,'Half-life [min]',line=2, font=2)
mtext(side=1, expression(paste(bolditalic('S. cerevisiae'), bold('-specific tAI'))), line=2)
label_plot(0.31, 600, 'G', cex=2, font=2)
Spom_cor <- cor(log(Pombe$half.life), log(Pombe$tAI),m='s',use='pairwise.complete.obs')
title <- paste0('Spearman rho=',round(Spom_cor,2))
pvalheatscatter(Pombe$tAI, Pombe$half.life, xlab="", ylab="", main='', xlim=c(0.38,0.55), log='xy', cex=0.3)
mtext(side=2,'Half-life [min]',line=2, font=2)
mtext(side=1, expression(paste(bolditalic('S. pombe'), bold('-specific tAI'))), line=2)
label_plot(0.36, 400, 'H', cex=2, font=2)
dev.off()


##' ### 5'UTR motif
p_AAACAAA = boxplot_kmer('AAACAAA',0,seqcol='UTR5_seq',ratecol='hlt.wt',dat=Scer,ylab='half-life (mins)',breaks=y_breaks, omitlessthan=30)
tmp <- data.frame(Scer)
fit1 = lm_motif(motif='AAACAAA',UTRs=tmp,n = 1,extend = 2,ratecol = 'hlt.wt',seqcol = 'UTR5_seq')
#p1 <- plot.coef(fit=fit1$fit,X=fit1$X,main=paste(fit1$X$consensus,'Puf3'),rates = F,ylim=c(0.8,1.1))
p_AAACAAA_coef <- readRDS('./data/5UTR_motif_coef_plot.rds')
# Conservation plot
dtpml <- readRDS("./data/motif_conservation_base_UTR5.rds")
dtpml <- dtpml[motif == 'AAACAAA']
dtpml[, base_position := pos - motif_start + 1, by = .(motif_i, ID, motif_start, motif, max_missmatch)]
#dtpml <- dtpml[, base_position := pos - motif_start + 1]
dtpml <- dtpml[, motif_length := nchar(motif)]
#dtpml <- dtpml[, base_position := ifelse(base_position > motif_length, paste('+', as.character(base_position-motif_length)), base_position)]
dtpml[, outside_motif := base_position < 1 | base_position > nchar(motif)]
dtpml[, base := substr(motif, base_position, base_position), by = .(motif_i, ID, motif_start, motif, max_missmatch)]
# Make pos_base the correct x-tick label, with the factor order same as base_position
dtpml[, pos_label := ifelse(base_position>motif_length, paste0('+',base_position-motif_length), # right outside of motif
                            ifelse(base_position < 1, as.character(base_position-1), as.character(base_position))), # left outside of motif
                            by = .(motif_i, ID, motif_start, motif, max_missmatch)]
dtpml$pos_label <- factor(dtpml$pos_label, c('-2','-1','1','2','3','4','5','6','7','8','+1','+2'))

motif_conservation_UTR5 <- ggplot(dtpml, aes(x = pos_label, y = conservation)) +
  #geom_boxplot() +
  stat_summary(fun.y='median', geom = 'point', size=2) +
  stat_summary(fun.y='median', geom = 'line', aes(group=1), size=1) +
  facet_wrap(~motif, scale='free', nrow = 1) +
  labs(x = 'Position') +
  theme_bw() + theme(strip.text=element_text(size=rel(1.1),face='bold'))

cairo_pdf('./fig/Paper/Supplements/Fig_S4_AAACAAA.pdf', width=4, height=10)
plot.new()
pushViewport(viewport(width=1, height=1,layout=grid.layout(5, 1)))
pushViewport(viewport(layout.pos.row = c(1:2), layout.pos.col = 1))
print(p_AAACAAA,newpage=FALSE)
grid.text(x=0.03,y=0.95,'A',gp=gpar(cex=1,font=2))
popViewport()
pushViewport(viewport(layout.pos.row = 3, layout.pos.col = 1))
print(p_AAACAAA_coef,newpage=FALSE)
grid.text(x=0.03,y=0.95,'B',gp=gpar(cex=1,font=2))
popViewport()
pushViewport(viewport(layout.pos.row = 4, layout.pos.col = 1))
seqlog.motif2(fit1$X$pwm,ic.scale = F,xaxis=F,yaxis=F)
grid.text(x=0.03,y=0.95,'C',gp=gpar(cex=1,font=2))
popViewport()
pushViewport(viewport(layout.pos.row = 5, layout.pos.col = 1))
print(motif_conservation_UTR5, newpage=FALSE)
grid.text(x=0.03,y=0.95,'D',gp=gpar(cex=1,font=2))
popViewport(2)
dev.off()

## Compare datasets
half_life_5_data <- readRDS("./data/half_life_5_data.rds")
half_life_5_data[, ATATTC := cut_lessthan(findmatch2(UTR3_seq, "ATATTC"), n = 2)]
half_life_5_data[, TGTAAATA := cut_lessthan(findmatch2(UTR3_seq, "TGTAAATA"), n = 2)]
half_life_5_data[, TGCAT := cut_lessthan(findmatch2(UTR3_seq, "TGCAT"), n = 2)]
half_life_5_data[, TTTTTTA := cut_lessthan(findmatch2(UTR3_seq, "TTTTTTA"), n = 2)]

lev <- c('Sun', "Miller", "Neymotin", "Gupta", "Presnyak", "Geisberg")
half_life_5_data$Geisberg <- as.numeric(half_life_5_data$Geisberg)
half_life_5_data <- melt(half_life_5_data, 
                         measure.vars = lev,
                         value.name = 'half_life',
                         variable.name = 'Study')

half_life_5_data <- melt(half_life_5_data,
                         measure.vars = c("ATATTC", "TGTAAATA", "TGCAT", "TTTTTTA"),
                         value.name = 'counts',
                         variable.name = 'motif')

half_life_5_data[, counts := factor(counts, levels = c('0', '1', ">=2"))]
half_life_5_data <- half_life_5_data[, Study := factor(Study, levels = lev)]

## Test for significance
pval <- half_life_5_data[, wilcox.test(half_life ~ counts == '0')$p.value,
                         by = c('motif', 'Study')]

pval <- sapply(pval$V1, function(i) deparse(format_pval(i)))

pdf('./fig/Paper/Supplements/Sun_Miller_Presnyak_motifs.pdf',width = 10,height = 7)
ggplot(half_life_5_data, aes(counts, half_life)) +
  geom_boxplot(outlier.size = 0.3) +
  annotate('text', x = ">=2", y = 150, label = pval, parse = TRUE, size = 2.5) +
  scale_y_continuous(breaks=c(5, 10, 20, 40, 150),trans='log2') +
  facet_grid(motif ~ Study) +
  theme_bw() +
  stat_summary(fun.data = give.n, geom = "text", fun.y = median, position = position_dodge(width = 0.75),size=rel(2)) +
  labs(title = "Functional motifs recovered in different datasets",
       x = 'Motif counts',
       y = 'Half-life [min]') +
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_cartesian(ylim = c(1, 150))
dev.off()

## S. cerevisiae UTR length problem
UTR_len <- UTRs[, .(genename, UTR5_length, UTR3_length)]
UTR_len <- melt(UTR_len, id.vars = 'genename')

cairo_pdf('./fig/Paper/UTR_Length.pdf', width = 6, height = 3)
ggplot(UTR_len, aes(log2(value))) +
  geom_histogram() +
  facet_wrap(~ variable) +
  labs(title = 'UTR length distribution of S. cerevisiae',
       x = 'log2 UTR length') +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = 'bold'))
dev.off()

UTR_len[, value < 6] %>%
  sum

UTR_len[, value == 1] %>%
  sum
