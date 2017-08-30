###################################################################################################
# Figure2 uAUG
###################################################################################################
# Load data and functions
source('./src/Paper/ini.R')

uAUG_counts <- lapply(UTRs$UTR5_seq, findStartCodon)
uAUG_counts <- lengths(uAUG_counts)
UTRs$uAUG_counts <- ifelse(uAUG_counts > 0, 'with uAUG', 'no uAUG')
#UTRs$uAUG_counts <- factor(uAUG_counts)
uAUG_to_plot <- melt(UTRs, id.vars=c('uAUG_counts'),measure.vars=c('WILD','XRN1','UPF2','UPF3'))
uAUG_to_plot$variable <- ifelse(uAUG_to_plot$variable=='XRN1','xrn1\u0394',ifelse(uAUG_to_plot$variable=='UPF2','upf2\u0394',ifelse(uAUG_to_plot$variable=='UPF3','upf3\u0394','WT')))
uAUG_to_plot$variable <- factor(uAUG_to_plot$variable,levels=c('WT','xrn1\u0394','upf2\u0394','upf3\u0394'))

Figure1_uAUG_plot <- ggplot(uAUG_to_plot, aes(x=uAUG_counts,y=value)) + 
  geom_boxplot(outlier.size = 0.5) + 
  scale_y_continuous(trans='log10',breaks=y_breaks, limits=c(1,200)) + 
  facet_grid(~variable,scales='free') + xlab('') + 
  ylab('Half-life [min]') +
  theme_bw() +
  theme(axis.title=element_text(size=rel(1.2),face='bold')) + 
  theme(axis.text=element_text(size=rel(1),face='bold')) +
  theme(strip.text.x=element_text(size=rel(1),face='bold.italic')) +
  stat_summary(fun.data = give.n, geom = "text", fun.y = median, position = position_dodge(width = 0.75),size=rel(3)) +
  theme(panel.margin=unit(0.5, "lines"), panel.grid.minor=element_blank())
## Do the testing
Figure1_uAUG_plot <- add_pval(Figure1_uAUG_plot, pairs = list(c(1,2)),heights=1.7,barheight=-0.1,pval_text_adj=-0.4, textsize=3)
# Compute fold change
test_dt <- data.table(Figure1_uAUG_plot$data)
fold_changes <- test_dt[ ,median(response, na.rm=T), by=c('variable', 'group')][,.SD[1]/.SD[2],by='variable',.SDcols='V1'][,V1]
fold_changes <- paste0('Median FC = ', round(fold_changes, 2))
Figure1_uAUG_plot <- Figure1_uAUG_plot + annotate('text', x=1.5, y=1, label=fold_changes, size=3)

##' ### Does uORF matters for uAUG?
##' #### Cerevisiae
uAUG_status <- get_uAUG_status(UTRs)
UTRs$uAUG_status <- factor(uAUG_status, labels=c('no uAUG', 'in-frame\nno PTC', 'out-frame\nPTC in CDS', 'uORF'))

# Number of uAUG...
uAUG_counts <- lapply(UTRs$UTR5_seq, findStartCodon)
uAUG_counts <- lengths(uAUG_counts)
uAUG_numbers <- boxplotOmmitN(UTRs$hlt.wt,uAUG_counts,breaks=y_breaks, main='S. cerevisiae', xlab='Number of uAUGs',ylim = c(1,300),
                    give.n.size=3,n.vjust=0,axis.title.size=1.2,plot.title.size=1.2,axis.text.size=0.8,fill='white',n=100)
uAUG_numbers <- uAUG_numbers + scale_fill_identity() + theme(plot.title=element_text(face='bold.italic'), panel.grid.minor=element_blank()) 
uAUG_numbers <- add_pval(uAUG_numbers, pairs = list(c(1,2),c(1,3),c(2,3)),textsize=3,
                             heights=c(150,230,100),barheight=0.5,pval_text_adj=rep(0.13,3))
uAUG_numbers <- uAUG_numbers + theme(plot.title = element_text(hjust = 0.5))

# Dissect different categories of uAUGs
uAUG_uORF <- ggplot(UTRs, aes(uAUG_status, hlt.wt)) +
  geom_boxplot(outlier.size = 0.5) +
  scale_y_continuous(breaks=y_breaks, trans='log10') +
  labs(x = '', y = 'Half-life [min]', title='S. cerevisiae') +
  theme_bw() +
  theme(axis.title=element_text(size=rel(1.2),face='bold')) + 
  theme(axis.text=element_text(size=rel(0.8),face='bold')) +
  theme(plot.title=element_text(face='bold.italic'), panel.grid.minor=element_blank()) +
  stat_summary(fun.data = give.n, geom = "text", fun.y = median, position = position_dodge(width = 0.75),size=rel(3)) 
uAUG_uORF <- add_pval(uAUG_uORF, pairs = list(c(1,2),c(1,3),c(1,4),c(3,4)),textsize=3,
                             heights=c(80,130,200,80),barheight=rep(0.05,4),pval_text_adj=rep(0.13,4))
uAUG_uORF <- add_pval(uAUG_uORF, pairs = list(c(1,2), c(3,4)), textsize=3,
                heights=c(1,1), barheight=rep(-0.1,2), pval_text_adj=rep(-0.2,2),
                annotation=c('no PTC', 'with PTC'))
uAUG_uORF <- uAUG_uORF + theme(plot.title = element_text(hjust = 0.5))

##' #### Pombe
uAUG_status <- get_uAUG_status(Pombe, CDSseqcol='cds_seqs')
Pombe$uAUG_status <- factor(uAUG_status, labels=c('no uAUG', 'in-frame\nno PTC', 'out-frame\nPTC in CDS', 'uORF'),
                                      levels=c('no ATG', 'Only_inframe_ATG', 'Out_frame_PTC_in_CDS', 'uORF'))

uAUG_uORF_pombe <- ggplot(Pombe, aes(uAUG_status, half.life)) +
  geom_boxplot(outlier.size = 0.5) +
  scale_y_continuous(breaks=y_breaks, trans='log10') +
  labs(x = '', y = 'Half-life [min]', title='S. pombe') +
  theme_bw() +
  theme(axis.title=element_text(size=rel(1.2),face='bold')) + 
  theme(axis.text=element_text(size=rel(0.8),face='bold')) +
  theme(plot.title=element_text(face='bold.italic'), panel.grid.minor=element_blank()) +
  coord_cartesian(ylim=c(1,250)) + 
  stat_summary(fun.data = give.n, geom = "text", fun.y = median, position = position_dodge(width = 0.75),size=rel(3))
uAUG_uORF_pombe <- add_pval(uAUG_uORF_pombe, pairs = list(c(1,2),c(1,3),c(1,4),c(3,4)),textsize=3,
                             heights=c(85,130,200,85),barheight=rep(0.05,4),pval_text_adj=rep(0.13,4))
uAUG_uORF_pombe <- add_pval(uAUG_uORF_pombe, pairs = list(c(1,2), c(3,4)), textsize=3,
                heights=c(1.5,1.5), barheight=rep(-0.1,2), pval_text_adj=rep(-0.4,2),
                annotation=c('no PTC', 'with PTC'))
uAUG_uORF_pombe <- uAUG_uORF_pombe + theme(plot.title = element_text(hjust = 0.5))

cairo_pdf('./fig/Paper/Figure2_uAUG.pdf', width=12, height=8)
plot.new()
pushViewport(viewport(width=1, height=1,layout=grid.layout(2, 3)))
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = c(1:3)))
print(Figure1_uAUG_plot, newpage=FALSE)
grid.text(x=0.02,y=0.95,'A',gp=gpar(cex=1.5,font=2))
popViewport()
pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
print(uAUG_numbers, newpage=FALSE)
grid.text(x=0.05,y=0.95,'B',gp=gpar(cex=1.5,font=2))
popViewport()
pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 2))
print(uAUG_uORF, newpage=FALSE)
grid.text(x=0.05,y=0.95,'C',gp=gpar(cex=1.5,font=2))
popViewport()
pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 3))
print(uAUG_uORF_pombe, newpage=FALSE)
grid.text(x=0.05,y=0.95,'D',gp=gpar(cex=1.5,font=2))
popViewport(2)
#grid.arrange(uAUG_numbers, uAUG_uORF, uAUG_uORF_pombe, ncol=3)
dev.off()
