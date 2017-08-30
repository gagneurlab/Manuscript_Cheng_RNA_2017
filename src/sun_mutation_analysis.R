##' ---
##' title: Analysis Sun Mutation data
##' author: Jun Cheng
##' ---
##'

#############################
# Analysis Sun Mutation data
#############################
library(gplots)
library(multtest)
source('./src/Paper/ini.R')

mutants = names(sun_mutants_hlt)
mutants <- sapply(mutants, tolower)
mutants[which(mutants=='wild')] = 'WT'
# Add Delta
mutants[which(mutants!='WT')] <- paste0(mutants[which(mutants!='WT')], '\u0394')
names(sun_mutants_hlt) = mutants
names(mutants) <- mutants

apply_to_all_mutants <- function(sun_mutants_hlt,FUN='getMotif_FC',nametoExtract='motif_effect_relative',label=group,id='.id',adj_m='BH',...){
  all_mutants <- lapply(sun_mutants_hlt,function(i) get(FUN)(i,...))
  all_mutants_stat <- lapply(all_mutants, function(i) data.table(fc=i[[nametoExtract]]))  # can be fold change or correlation
  all_mutants_stat <- rbindlist(all_mutants_stat,idcol = T)
  pvars <- cross_test_dt(all_mutants_stat,measurevar = 'fc',ref='WT',adjust = F)
  adjusted <- mt.rawp2adjp(pvars, proc=adj_m)
  pvars_adj <- adjusted$adj[order(adjusted$index), ][ ,2]
  pvars_stars <- data.table(.id = names(pvars), sig = fdr2star(pvars_adj), pval = pvars, pval_adj = pvars_adj)
  all_mutants_stat_summary <- summarySE(all_mutants_stat,measurevar = 'fc',groupvars = id)
  all_mutants_stat_summary <- merge(all_mutants_stat_summary,pvars_stars,by=id)
  all_mutants_stat_summary$Group <- label
  all_mutants_stat_summary <- data.table(all_mutants_stat_summary)
  all_mutants_stat_summary <- all_mutants_stat_summary[order(group,fc)]
  # Make WILD the first one
  all_mutants_stat_summary$.id <- factor(all_mutants_stat_summary$.id, levels=c('WT',all_mutants_stat_summary$.id[all_mutants_stat_summary$.id!='WT']))
  ## Extract plots
  # check whether plot exist
  if('plots' %in% names(all_mutants[[1]])){
    plots <- lapply(all_mutants, function(i) i[['plots']])
  }else{
    plots <- NULL
  }
  return(list(summary=all_mutants_stat_summary,stat=all_mutants_stat,plots=plots))
}


plot_ko <- function(rt_obj, ylab=''){
  ggplot(rt_obj$summary,aes(x=.id,y=fc)) +
    geom_bar(stat="identity",aes(fill=Group)) +
    #geom_errorbar(aes(ymin=fc-se,ymax=fc+se),width=.2)+
    theme(axis.text.x=element_text(angle = -90, size=20,vjust=0.5)) +
    theme(axis.title.y = element_text(size = 20, angle = 90)) +
    xlab('') + ylab(ylab)+
    geom_text(aes(.id,fc,label=sig),vjust=0.8,hjust=-1,size=10,color='red') +
    geom_point(data=rt_obj$stat,aes(.id,fc),position = position_jitter(width=.15),color='blue') +
    geom_hline(yintercept=rt_obj$summary[rt_obj$summary$.id=='WT',fc]) +
    #theme(axis.text.x=element_text(colour=c(rep('black',length(sun_mutants_hlt)-2),"red",'black'))) +
    coord_flip()+theme_bw()
}


#' ## Group the knockout genes
## CCR4-NOT comple
CCR4_NOT <- c('CCR4',
              'CAF40',
              'NOT3',
              'PAN2',
              'PAN3',
              'POP2')
CCR4_NOT <- sapply(CCR4_NOT, tolower)
CCR4_NOT <- sapply(CCR4_NOT, function(i) paste0(i, '\u0394'))
## Decapping
Decapping <- c( 'DCS1',
                'DCS2',
                'DHH1',
                'EDC2',
                'EDC3',
                'LSM1',
                'LSM6',
                'LSM7',
                'PAT1',
                'SCD6')
Decapping <- sapply(Decapping, tolower)
Decapping <- sapply(Decapping, function(i) paste0(i, '\u0394'))
## RNA surveillance
Surveillance <- c('DOM34',
                  'HBS1',
                  'UPF2',
                  'UPF3')
Surveillance <- sapply(Surveillance, tolower)
Surveillance <- sapply(Surveillance, function(i) paste0(i, '\u0394'))
## PUF
RBP <- c(       'PUB1',
                'PUF1',
                'PUF2',
                'PUF3',
                'PUF4',
                'PUF5',
                'PUF6')
RBP <- sapply(RBP, tolower)
RBP <- sapply(RBP, function(i) paste0(i, '\u0394'))
## Exosome
Exosome <- c(   'SKI2',
                'SKI3',
                'SKI7',
                'SKI8',
                'RRP47',
                'RRP6')
Exosome <- sapply(Exosome, tolower)
Exosome <- sapply(Exosome, function(i) paste0(i, '\u0394'))
group <- ifelse(mutants %in% CCR4_NOT, 'Deadenylation',
                ifelse(mutants %in% Decapping, 'Decapping',
                       ifelse(mutants %in% Surveillance, 'Surveillance',
                              ifelse(mutants %in% RBP,'RBP',
                                     ifelse(mutants %in% Exosome,'Exosome',mutants)))))
group[which(group=='xrn1\u0394')] <- 'XRN1'
group <- data.table(mutants=mutants,group=group)
group <- group[order(tolower(mutants))]
group <- group$group
names(group) <- mutants

#' ## Start codon context
startMinus3 = as.character(sapply(UTRs$UTR5_seq, function(x) substr(x,nchar(x)-2,nchar(x)-2)))
startMinus3[startMinus3==''] = NA
names(startMinus3) = UTRs$genename
startMinus3 <- startMinus3[!is.na(startMinus3)]
startMinus3 <- as.data.frame(startMinus3)
getstartMinus3 <- function(hltdmat,startMinus3,explained_variance=FALSE){
  hltdmat <- merge(hltdmat,startMinus3,by='row.names')
  hltdmat <- data.table(hltdmat)
  hltdmat <- hltdmat[complete.cases(hltdmat)]
  replicates <- colnames(hltdmat)
  replicates <- replicates[2:(length(replicates)-2)]
  C2A <- c()
  C2A_relative <- c() # Fold change relative to C, from C to A
  plots <- list()
  for(i in replicates){
    p = boxplotOmmitN(hltdmat[,get(i)], hltdmat[,startMinus3],n = 30,main=i,ylim = c(1,150))
    if(explained_variance){
      ##fit a linear model, and get the explained variance
      tmp <- hltdmat[,.(get(i),startMinus3)]
      setnames(tmp,c(i,'startMinus3'))
      fit <- lm(data=tmp, get(i) ~ startMinus3)
      differ <- summary(fit)$r.squared  # explained variance
      C2A_relative <- c(C2A_relative,differ)
    }else{
      ## Use median(A) - median(C)
      tmp = hltdmat[,.(hlt=median(get(i))),by=startMinus3]
      differ = tmp[startMinus3=='A',hlt] - tmp[startMinus3=='C',hlt]
      C2A_relative <- c(C2A_relative,differ/tmp[startMinus3=='C',hlt])
    }
    C2A = c(C2A,differ)
    plots <- list(plots,p)
  }
  return(list(plots=plots,C2A=C2A,C2A_relative=C2A_relative))
}

startMinus3_mutants_fc <- apply_to_all_mutants(sun_mutants_hlt,FUN = 'getstartMinus3',nametoExtract = 'C2A_relative',startMinus3=startMinus3)
cairo_pdf('./fig/startMinus3.pdf',width=7,height=9)
plot_ko(startMinus3_mutants_fc, ylab='Fold of half-life decrease from C to A')
dev.off()


#' ## Codon usage, variance explained
codon_counts <- data.frame(UTRs[ ,c(codons, "CDS_length"), with=F])
rownames(codon_counts) <- UTRs$genename
getcodon_exp_var <- function(hltdmat,codontable){
  replicates <- colnames(hltdmat)
  # Find out which are the mutant replicate half-life
  replicates <- replicates[1:(length(replicates)-1)]
  hltdmat <- merge(hltdmat, codontable,by='row.names')
  hltdmat <- data.table(hltdmat)
  hltdmat <- hltdmat[complete.cases(hltdmat)]
  codon_exp_var <- c()
  plots <- list()
  for(i in replicates){
    # Use regression
    dat <- hltdmat[, c(codons, "CDS_length"), with=F]
    hlt <- hltdmat[, log(get(i))]
    dat <- data.table(hlt=hlt, dat)
    for (cod in codons){
      dat[, (cod) := get(cod) / CDS_length]
    }
    dat <- dat[, TTT := NULL] # The reference
    fit <- lm(data=dat, hlt ~ .)
    s_cor <- cv_rsquare2(fit)
    codon_exp_var <- c( codon_exp_var, s_cor)
    p <- ggplot(data.frame(prediction=predict(fit), measured=dat$hlt),aes(x=prediction,y=measured))
    p <- p + geom_point(alpha=0.3,stat='identity')
    plots <- list(plots,p)
  }
  return(list(codon_exp_var=codon_exp_var,plots=plots))
}

codon_var_explained <- apply_to_all_mutants(sun_mutants_hlt,FUN = 'getcodon_exp_var',nametoExtract = 'codon_exp_var',codontable=codon_counts)
cairo_pdf('./fig/Codon_regression.pdf',width=7,height=9)
codon_mutants_plot <- plot_ko(codon_var_explained, ylab='Correlation of tRNA adaptation index and mRNA half-life')
print(codon_mutants_plot)
dev.off()

# Same plot horizontal direction
cairo_pdf('./fig/Codon_regression2.pdf',width=8,height=4)
p_codons <- ggplot(codon_var_explained$summary,aes(x=.id,y=fc)) +
  geom_bar(stat="identity",aes(fill=Group)) +
  xlab('') + ylab('Codon-explained-variance (%)')+
  geom_text(aes(.id,fc,label=sig),vjust=-0.2,hjust=0.5,size=7,color='red') +
  geom_point(data=codon_var_explained$stat,aes(.id,fc),position = position_jitter(width=.15),color='blue',alpha=0.7) +
  geom_hline(yintercept=codon_var_explained$summary[codon_var_explained$summary$.id=='WT',fc]) + theme_bw() +
  theme(axis.text.x=element_text(size = rel(1.3), angle = 90, face='italic', hjust=1, vjust=0.5)) +
  theme(axis.title.x=element_text(size = rel(1.2), vjust=-0.5), axis.title.y = element_text(size=rel(1.2), hjust = 0.8)) +
  scale_y_continuous(labels=c(0,20,40,60)) + guides(fill=guide_legend(title='Group'))
print(p_codons)
dev.off()

gene_order <- levels(codon_var_explained$summary$.id)
saveRDS(p_codons,file='./fig/ggplot_obj/codon_sig_knockout.rds')

# Get the color for later usage
g <- ggplot_build(p_codons)
pathway_color <- c(g$data[[1]]["fill"])[[1]]
#pathway_color <- pathway_color[1:(length(pathway_color) / 3)]
names(pathway_color) <- levels(codon_var_explained$summary$.id)

#' ## Stop + 1, compare TGAC and TGAG, but plot all
stopPlus1 <- UTRs[, .(genename, stopPlus1)]
getstopPlus1 <- function(hltdmat,stopPlus1,cod1='TGAG',cod2='TGAC'){
  genename <- rownames(hltdmat)
  hltdmat <- data.table(hltdmat)
  hltdmat$genename <- genename
  hltdmat <- merge(hltdmat,stopPlus1,by='genename')
  hltdmat <- hltdmat[complete.cases(hltdmat)]
  replicates <- colnames(hltdmat)
  replicates <- replicates[2:(length(replicates)-2)]
  TGAC2TGAG <- c()
  TGAC2TGAG_relative <- c() # Fold change relative to TGAC
  plots <- list()
  for(i in replicates){
    p = boxplotOmmitN(hltdmat[,get(i)], hltdmat[,stopPlus1],n = 30,main=i)
    tmp = hltdmat[,.(hlt=median(get(i))),by=stopPlus1]
    differ = tmp[stopPlus1==cod2,hlt] - tmp[stopPlus1==cod1,hlt]
    TGAC2TGAG = c(TGAC2TGAG,differ)
    TGAC2TGAG_relative <- c(TGAC2TGAG_relative,differ/tmp[stopPlus1==cod1,hlt]) # Fold increase from cod1 to cod2
    plots <- list(plots,p)
  }
  return(list(plots=plots,TGAC2TGAG=TGAC2TGAG,TGAC2TGAG_relative=TGAC2TGAG_relative))
}

TGAC2TGAG_fc <- apply_to_all_mutants(sun_mutants_hlt,FUN = 'getstopPlus1',nametoExtract = 'TGAC2TGAG_relative',stopPlus1=stopPlus1,cod1='TGAC',cod2='TGAG')
cairo_pdf('./fig/TGAC2TGAG_fc.pdf',width=7,height=9)
plot_ko(TGAC2TGAG_fc, ylab='Fold decrease from TGAG to TGAC')
dev.off()

TGAC2TAAG_fc <- apply_to_all_mutants(sun_mutants_hlt,FUN = 'getstopPlus1',nametoExtract = 'TGAC2TGAG_relative',stopPlus1=stopPlus1,cod1='TGAC',cod2='TAAG')
cairo_pdf('./fig/TGAC2TAAG_fc.pdf',width=7,height=9)
plot_ko(TGAC2TAAG_fc, ylab='Fold increase from TGAC to TAAG')
dev.off()


#' ## uAUG
##' ### In and out of frame uAUG
##'
getcod <- function(UTR5,cod='ATG',reverse=T){
  # reverse is true for UTR5, not for cds or UTR3
  if(reverse){
    cod_rev <- Biostrings::reverse(cod)
  }else{
    cod_rev <- cod
  }
  num_cod <- as.numeric(sapply(UTR5, function(x) count_triplet(cod,x)$num))
  num_cod_inframe <- sapply(UTR5, function(x) count_triplet_inframe(cod_rev,x,rev = reverse))
  cod_status <- ifelse(num_cod>num_cod_inframe,'Out-frame',num_cod)
  cod_status <- ifelse(cod_status==0,'No uAUG',cod_status)
  cod_status <- ifelse(cod_status==num_cod_inframe,'In-frame',cod_status)
  cod_status <- factor(cod_status,levels=c('No uAUG','In-frame','Out-frame'))
  return(data.table(num_cod=num_cod,status=cod_status))
}

## subset no NA UTR5 sequence
UTRs <- subset(UTRs, !is.na(UTR5_seq))
uAUG_status <- getcod(UTRs$UTR5_seq,cod='ATG')
uAUG_status[,genename := UTRs$genename]
# Just distinguish whether have or not have uAUG
uAUG_status$status <- ifelse(uAUG_status$num_cod > 0, 'with', 'without')

getuAUG <- function(hltdmat,uAUG_status){
  genename <- rownames(hltdmat)
  hltdmat <- data.table(hltdmat)
  hltdmat$genename <- genename
  hltdmat <- merge(hltdmat,uAUG_status,by='genename')
  hltdmat <- hltdmat[complete.cases(hltdmat)]
  replicates <- colnames(hltdmat)
  replicates <- replicates[2:(length(replicates)-3)]
  uAUG_effect <- c()
  uAUG_effect_relative <- c() # Fold change relative to with uAUG
  plots <- list()
  for(i in replicates){
    p = boxplotOmmitN(hltdmat[,get(i)], hltdmat[,status],n = 30,main=i,ylim=c(1,150))
    tmp = hltdmat[,.(hlt=median(get(i))),by=status]
    differ = tmp[status=='with',hlt] - tmp[status=='without',hlt]
    uAUG_effect = c(uAUG_effect,differ)
    uAUG_effect_relative <- c(uAUG_effect_relative,differ/tmp[status=='without',hlt])
    plots <- list(plots,p)
  }
  return(list(plots=plots,uAUG_effect=uAUG_effect,uAUG_effect_relative=uAUG_effect_relative))
}

uAUG_mutants_fc <- apply_to_all_mutants(sun_mutants_hlt,FUN = 'getuAUG',nametoExtract = 'uAUG_effect_relative',uAUG_status=uAUG_status)
cairo_pdf('./fig/uAUG_mutants_fc.pdf',width=7,height=9)
plot_ko(uAUG_mutants_fc, ylab='Fold increase without uAUG')
dev.off()

#' ## 5' UTR folding energy
foldingenergy5 <- UTRs[,.(genename,Free.energy.of.ensemble)]
# Remove NA
foldingenergy5 <- foldingenergy5[complete.cases(foldingenergy5)]
foldingenergy5$Free.energy.of.ensemble = -foldingenergy5$Free.energy.of.ensemble
geteng_cor <- function(hltdmat,foldingenergy5){
  genename <- rownames(hltdmat)
  hltdmat <- data.table(hltdmat)
  hltdmat$genename <- genename
  hltdmat <- merge(hltdmat,foldingenergy5,by='genename')
  hltdmat <- hltdmat[complete.cases(hltdmat)]
  replicates <- colnames(hltdmat)
  replicates <- replicates[2:(length(replicates)-2)]
  eng_cor <- c()
  plots <- list()
  for(i in replicates){
    eng_cor <- c( eng_cor,-cor(hltdmat[,get(i),with=T],hltdmat$Free.energy.of.ensemble,m='s') )
    p <- ggplot(hltdmat,aes(x=get(i),y=Free.energy.of.ensemble))
    p <- p + geom_point(alpha=0.3,stat='identity')
    p <- p + scale_x_continuous(trans='log')
    p <- p + scale_y_continuous(trans='log')
    plots <- list(plots,p)
  }
  return(list(eng_cor=eng_cor,plots=plots))
}

eng5_mutants_cor <- apply_to_all_mutants(sun_mutants_hlt,FUN = 'geteng_cor',nametoExtract = 'eng_cor',foldingenergy5=foldingenergy5)
cairo_pdf('./fig/eng5_cor.pdf',width=7,height=9)
plot_ko(eng5_mutants_cor, ylab='correlation of 5\' folding energy and half-life')
dev.off()

#' ## Get motif fold change
UTRseqDT <- UTRs[,.(genename, UTR3_seq, UTR5_seq)]
UTRseqDT <- UTRseqDT[nchar(UTR3_seq)>10 & nchar(UTR5_seq) > 10]
getMoif_FC <- function(hltdmat,motif,UTRseqDT,seqcol='UTR3_seq'){
  UTRseqDT_fun <- copy(UTRseqDT)
  genename <- rownames(hltdmat)
  hltdmat <- data.table(hltdmat)
  hltdmat$genename <- genename
  # Grep number of motifs
  counts = as.numeric(sapply(UTRseqDT_fun[,get(seqcol)],function(s) findmatch(s,motif)$count ))
  UTRseqDT_fun[,c(seqcol):=NULL][,c('counts'):=counts]
  hltdmat <- merge(hltdmat,UTRseqDT_fun,by='genename')
  hltdmat <- hltdmat[complete.cases(hltdmat)]
  replicates <- colnames(hltdmat)
  replicates <- replicates[2:(length(replicates)-3)]
  motif_effect <- c() ## Compare with genes without motif, the effect of one motif
  motif_effect_relative <- c() # Fold change relative to no motif
  plots <- list()
  for(i in replicates){
    p = boxplotOmmitN(hltdmat[,get(i)], hltdmat[,counts],n = 30,main=i, ylim=c(1,150))
    tmp = hltdmat[,.(hlt=median(get(i))),by=counts]
    differ = tmp[counts==1,hlt] - tmp[counts==0,hlt]
    motif_effect = c(motif_effect,differ)
    motif_effect_relative <- c(motif_effect_relative,differ/tmp[counts==0,hlt])
    plots <- list(plots,p)
  }
  return(list(plots=plots,motif_effect=motif_effect,motif_effect_relative=motif_effect_relative))
}


Puf3 <- apply_to_all_mutants(sun_mutants_hlt,FUN = 'getMoif_FC',nametoExtract = 'motif_effect_relative',motif='TGTAAATA',UTRseqDT=UTRseqDT,seqcol='UTR3_seq')
ATATTC <- apply_to_all_mutants(sun_mutants_hlt,FUN = 'getMoif_FC',nametoExtract = 'motif_effect_relative',motif='ATATTC',UTRseqDT=UTRseqDT,seqcol='UTR3_seq')
Whi3 <- apply_to_all_mutants(sun_mutants_hlt,FUN = 'getMoif_FC',nametoExtract = 'motif_effect_relative',motif='TGCAT',UTRseqDT=UTRseqDT,seqcol='UTR3_seq')
TTTTTTA <- apply_to_all_mutants(sun_mutants_hlt,FUN = 'getMoif_FC',nametoExtract = 'motif_effect_relative',motif='TTTTTTA',UTRseqDT=UTRseqDT,seqcol='UTR3_seq')
AAACAAA <- apply_to_all_mutants(sun_mutants_hlt,FUN = 'getMoif_FC',nametoExtract = 'motif_effect_relative',motif='AAACAAA',UTRseqDT=UTRseqDT,seqcol='UTR5_seq')

pdf('./fig/motifs_ko.pdf',width=7,height=9)
plot_ko(Puf3)
plot_ko(ATATTC)
plot_ko(Whi3)
plot_ko(TTTTTTA)
plot_ko(AAACAAA)
dev.off()

#### Draw a heatmap
data_plot  <- rbindlist(list(startMinus3_mutants_fc=startMinus3_mutants_fc$summary,
                             codon_var_explained=codon_var_explained$summary,
                             TGAC2TGAG_fc=TGAC2TGAG_fc$summary,
                             uAUG_mutants_fc=uAUG_mutants_fc$summary,
                             eng5_mutants_cor=eng5_mutants_cor$summary,
                             Puf3=Puf3$summary,
                             ATATTC=ATATTC$summary,
                             Whi3=Whi3$summary,
                             TTTTTTA=TTTTTTA$summary,
                             AAACAAA=AAACAAA$summary),idcol = 'feature')
data_plot$.id <- as.character(data_plot$.id)
data_plot_wide <- dcast(data_plot[,.(feature,fc,.id)], feature ~.id, value.var = 'fc')
data_plot_wide <- t(data.frame(data_plot_wide))
colnames(data_plot_wide) <- data_plot_wide[1,]
data_plot_wide <- data_plot_wide[-1,]
hc <- hclust(dist(data_plot_wide, method = "euclidean"), method = "complete" )
ord <- hc$order

## Order by a clustering
dend_order <- rownames(data_plot_wide)[ord]
# Make wide type first
dend_order <- c('WT',dend_order[dend_order!='WT'])
data_plot$.id <- factor(data_plot$.id, levels=dend_order)
data_plot$color <- pathway_color[dend_order]

## Order by pathway
data_plot$.id <- factor(data_plot$.id, levels=gene_order)
#levels(data_plot$.id) <- c('WT',data_plot$.id[data_plot$.id!='WT'])
data_plot$color <- pathway_color[gene_order]

# Order feature from 5' to 3'
data_plot$feature <- factor(data_plot$feature,
                            levels = c('uAUG_mutants_fc',
                                       'AAACAAA',
                                       'eng5_mutants_cor',
                                       'startMinus3_mutants_fc',
                                       'codon_var_explained',
                                       'TGAC2TGAG_fc',
                                       'Puf3',
                                       'ATATTC',
                                       'Whi3',
                                       'TTTTTTA'))

cairo_pdf('./fig/FigEV5_summary_heatmap.pdf',width = 5,height = 10)
p <- ggplot(aes(x=.id, y=feature, fill=fc), data=data_plot)
p <- p + geom_tile() + scale_fill_gradient2(low="red", mid="white", high="blue") +
  #   geom_text(aes(label=stars, color=value), size=8) + scale_colour_gradient(low="grey30", high="white", guide="none") +
  geom_text(aes(label=sig), color="black", size=rel(5),vjust=0.8) +
  labs(y=NULL, x=NULL, fill="Effects") +
  theme_bw() + theme(axis.text=element_text(hjust = 0, size=rel(1.5)),axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(axis.text.x=element_text(colour=data_plot$color, face='italic')) +
  scale_y_discrete(labels = c("uAUG","AAACAAA","5\' energy","St-3 C-A","codon usage","TGAC-TGAG","TGTAAATA","ATATTC","TGCAT","TTTTTTA"))
print(p)
dev.off()
