###################################################################################################
# Figure 5 motifs
###################################################################################################
# Load data and functions
source('./src/Paper/ini.R')

##' ## Half-life determinant motifs in Cerevisiae
motif_order <- c('ATATTC','TGCAT','TGTAAATA','TTTTTTA')
p_ATATTC = boxplot_kmer('ATATTC',0,seqcol='UTR3_seq',ratecol='hlt.wt',dat=Scer,ylab='half-life (mins)',breaks = y_breaks)
p_TGTAAATA = boxplot_kmer('TGTAAATA',0,seqcol='UTR3_seq',ratecol='hlt.wt',dat=Scer,ylab='half-life (mins)',breaks = y_breaks)
p_TGCAT = boxplot_kmer('TGCAT',0,seqcol='UTR3_seq',ratecol='hlt.wt',dat=Scer,ylab='half-life (mins)',breaks = y_breaks)
p_AAACAAA <- boxplot_kmer('AAACAAA',0,seqcol='UTR5_seq',ratecol='hlt.wt',dat=Scer,ylab='half-life (mins)',breaks = y_breaks)
p_TTTTTTA <- boxplot_kmer('TTTTTTA',0,seqcol='UTR3_seq',ratecol='hlt.wt',dat=Scer,ylab='half-life (mins)',breaks = y_breaks,omitlessthan = 20)
## Assemble them into one
motif_data <- rbindlist(list(ATATTC=p_ATATTC$data,
                             TGTAAATA=p_TGTAAATA$data,
                             TGCAT=p_TGCAT$data,
                             TTTTTTA=p_TTTTTTA$data), idcol=TRUE)
a <- motif_data$group
b <- ifelse(a=='(-Inf,0]', '0',
            ifelse(a=='(0,1]', '1',
                   ifelse(a=='(1,2]', '2',
                          ifelse(a=='(1, Inf]', '>=2',
                                 ifelse(a=='(0, Inf]', '>=1',
                                        ifelse(a=='(2, Inf]', '>=3',NA))))))
motif_data$group <- factor(b, levels=c('0', '1', '>=1', '2', '>=2', '>=3'))

motif_data$.id <- factor(motif_data$.id, levels = motif_order)

motif_labeller <- list('ATATTC'='ATATTC\nFDR=1.96e−05',
                       'TGCAT'='TGCAT\nFDR=6.94e−04',
                       'TGTAAATA'='TGTAAATA\nFDR=3.2e−05',
                       'TTTTTTA'='TTTTTTA\nFDR=0.086')
plot_labeller <- function(variable,value){
  return(motif_labeller[value])
}
motif_boxplot <- ggplot(motif_data, aes(group, response)) +
  geom_boxplot(fill='white', outlier.size = 0.5) +
  theme_bw() +
  scale_y_log10(breaks=y_breaks) +
  facet_grid(. ~ .id, scale='free', labeller=plot_labeller) +
  labs(x='Occurence per transcript', y='Half-life [min]') +
  stat_summary(fun.data = give.n, geom = "text", fun.y = median, position = position_dodge(width = 0.75),size=rel(3)) +
  #theme(axis.title=element_text(size=rel(1.2), face='bold'),
  #      axis.text=element_text(size=rel(1.2), face='bold')) +
  theme(strip.text=element_text(size=rel(1.1),face='bold')) +
  theme(panel.spacing=unit(0.5, "lines"))

##' Single-nucleotide contribution
UTRs <- data.frame(Scer)
fit1 = lm_motif(motif='TGTAAATA',UTRs=UTRs,n = 1,extend = 2,ratecol = 'hlt.wt',seqcol = 'UTR3_seq')
## p1 <- plot.coef(fit=fit1$fit,X=fit1$X,main=paste(fit1$X$consensus,'Puf3'),rates = F,ylim=c(0.75,1.7))
fit2 = lm_motif(motif='ATATTC',UTRs=UTRs,n = 1,extend = 2,ratecol = 'hlt.wt',seqcol = 'UTR3_seq')
## p2 <- plot.coef(fit=fit2$fit,X=fit2$X,main=paste(fit2$X$consensus,''),rates = F,ylim=c(0.75,1.7))
fit3 = lm_motif(motif='TGCAT',UTRs=UTRs,n = 1,extend = 2,ratecol = 'hlt.wt',seqcol = 'UTR3_seq')
## p3 <- plot.coef(fit=fit3$fit,X=fit3$X,main=paste(fit3$X$consensus,'Whi3'),rates = F,ylim=c(0.75,1.7))
fit5 <- lm_motif(motif='TTTTTTA',UTRs=UTRs,n = 1,extend = 2,ratecol = 'hlt.wt',seqcol = 'UTR3_seq')
## coef_data <- rbindlist(list(TGTAAATA=p1$data,ATATTC=p2$data,TGCAT=p3$data,ATATATGC=p4$data),idcol=T)

## coef_data fit from ./src/analyze-snv-perturbation.R
coef_data <- readRDS('./data/motif_coef_data.rds')

coef_plot <- ggplot(data=coef_data,aes(y=coef,x=pos_label,group=base,colour=base))  +
  geom_point(size=1) +
  geom_line(size=0.7) +
  facet_wrap(~.id,scales = "free_x",nrow=1)+
  theme_bw() +
  labs(x='Position',y='Effects on half-life') +
  geom_hline(yintercept=1) +
  theme(strip.text=element_text(size=rel(1.1),face='bold')) +
  theme(panel.spacing=unit(2, "lines")) +
  #theme(legend.position='top')
  theme(legend.position='none')

## This is conservation score analysis from Ziga
## deepcis/src/r/conservation/motif-analysis.R
dtpml <- readRDS("./data/motif_conservation_base_UTR3.rds")
dtpml[, base_position := pos - motif_start + 1, by = .(motif_i, ID, motif_start, motif, max_missmatch)]
#dtpml <- dtpml[, base_position := pos - motif_start + 1]
dtpml <- dtpml[, motif_length := nchar(motif)]
#dtpml <- dtpml[, base_position := ifelse(base_position > motif_length, paste('+', as.character(base_position-motif_length)), base_position)]
dtpml[, outside_motif := base_position < 1 | base_position > nchar(motif)]
dtpml[, base := substr(motif, base_position, base_position), by = .(motif_i, ID, motif_start, motif, max_missmatch)]
dtpml$motif <- factor(dtpml$motif, motif_order)
# Make pos_base the correct x-tick label, with the factor order same as base_position
dtpml[, pos_label := ifelse(base_position>motif_length, paste0('+',base_position-motif_length), # right outside of motif
                            ifelse(base_position < 1, as.character(base_position-1), as.character(base_position))), # left outside of motif
                            by = .(motif_i, ID, motif_start, motif, max_missmatch)]
dtpml$pos_label <- factor(dtpml$pos_label, c('-2','-1','1','2','3','4','5','6','7','8','+1','+2'))

motif_conservation <- ggplot(dtpml, aes(x = pos_label, y = conservation)) +
  #geom_boxplot() +
  stat_summary(fun.y='median', geom = 'point', size=2) +
  stat_summary(fun.y='median', geom = 'line', aes(group=1), size=1) +
  facet_wrap(~motif, scale='free', nrow = 1) +
  labs(x = 'Position') +
  theme_bw() + theme(strip.text=element_text(size=rel(1.1),face='bold'))

# Motif positinal bias  colors=c('#d62f0c','#8C96C6','#F781BF','#4DAF4A')
motif_density <- density_kmer_ggplot(kmers=c('ATATTC','TGTAAATA','TGCAT','TTTTTTA'),colors=rep("#000000",4),
                          dat=Scer,seqcol='UTR3_seq',ratecol='hlt.wt',window=20,max.length=200,
                          align.right=T,promoter=F)

motif_density$data$motif <- factor(motif_density$data$motif, motif_order)
motif_density <- motif_density +
  facet_wrap(~motif, nrow = 1) +
  theme_bw() + theme(strip.text=element_text(size=rel(1.1),face='bold')) +
  theme(panel.margin=unit(0.5, "lines")) + scale_x_continuous(labels=c(150,100,50,0))

## A position holder for seqlogo
coef_data$Frequency <- seq(0,1,0.1)
place_holder <- ggplot(coef_data, aes(pos_label, Frequency)) +
  facet_wrap(~.id, nrow=1, scales = 'free_x') +
  theme_classic() +
  theme(axis.line.x=element_line(),
        axis.line.y=element_line()) +
  #theme(strip.text=element_blank(), strip.background=element_blank()) +
  theme(panel.spacing=unit(2, "lines")) +
  xlab('Position')


## Read validation plot
motif_effect <- readRDS('./fig/ggplot_obj/validation_motif_effect.rds')

cairo_pdf('./fig/Paper/Figure5_motifs.pdf',width = 10,height=12)
plot.new()
pushViewport(viewport(layout = grid.layout(8, 4)))
pushViewport(viewport(layout.pos.col = c(1,4),layout.pos.row=c(1,2)))
print(motif_boxplot, newpage = FALSE)
grid.text(x=0.03,y=0.95,'A',gp=gpar(cex=1,font=2))
popViewport()
pushViewport(viewport(layout.pos.col=c(1,4),layout.pos.row=3)) # The last row
print(motif_density, newpage = FALSE)
grid.text(x=0.03,y=0.95,'B',gp=gpar(cex=1,font=2))
popViewport()
pushViewport(viewport(layout.pos.col = c(1,4),layout.pos.row=4))
print(coef_plot, newpage = FALSE)
grid.text(x=0.03,y=0.95,'C',gp=gpar(cex=1,font=2))
popViewport()
pushViewport(viewport(layout.pos.col=c(1,4),layout.pos.row=5))
print(place_holder, newpage=FALSE)
grid.text(x=0.03,y=0.95,'D',gp=gpar(cex=1,font=2))
popViewport()
pushViewport(viewport(layout.pos.col = c(1,4), layout.pos.row=6))
print(motif_conservation, newpage = FALSE)
grid.text(x=0.03,y=0.95,'E',gp=gpar(cex=1,font=2))
popViewport()
pushViewport(viewport(layout.pos.col = c(1,2), layout.pos.row=c(7,8)))
print(motif_effect, newpage = FALSE)
grid.text(x=0.03,y=0.95,'F',gp=gpar(cex=1,font=2))
popViewport()
pushViewport(viewport(layout.pos.col = c(4), layout.pos.row=c(7,8)))
popViewport(2) # Pop out the last row
dev.off()


cairo_pdf('./fig/Paper/Figure5_motifs2.pdf',width = 10,height=5/3)
pushViewport(viewport(layout = grid.layout(1,4))) # define a 1X4 grid on the last row
pushViewport(viewport(layout.pos.col=1,layout.pos.row=1))
seqlog.motif2(fit2$X$pwm,ic.scale = F,xaxis=F,yaxis=F)
popViewport() # pop out from this viewport
pushViewport(viewport(layout.pos.col=2,layout.pos.row=1))
seqlog.motif2(fit3$X$pwm,ic.scale = F,xaxis=F,yaxis=F)
popViewport()
pushViewport(viewport(layout.pos.col=3,layout.pos.row=1))
seqlog.motif2(fit1$X$pwm,ic.scale = F,xaxis=F,yaxis=F)
popViewport()
pushViewport(viewport(layout.pos.col=4,layout.pos.row=1))
seqlog.motif2(fit5$X$pwm,ic.scale = F,xaxis=F,yaxis=F)
popViewport(2)
dev.off()
