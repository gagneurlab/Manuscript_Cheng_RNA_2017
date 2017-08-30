# utilities, usefull functions

###################################################
# Find the number of matches. and return position
###################################################

# A function to permutate given string, allow n mismatches, return a regular expression ready to mach
seqdic <- function(motif,n){
  seqlength = nchar(motif)
  if(n>seqlength){
    stop('Number of mismatch cannot larger than string length')
  }
  subseq = strsplit(motif,NULL)[[1]]
  #permute = function(subseq,seqlength,n){
  allComb = lapply( 1:n, function(i) combn(seqlength,i,simplify = F) ) # All permutation combinations, from 1 mismatch to n mismatches
  dic =c( paste(subseq,collapse='') )  # all permutated sequences
  for (i in seq(n)){
    # i mismatches
    tmp = allComb[[i]] # length choose(seqlength,i) possible mismatch combinations
    listlength = choose(seqlength,i)
    for (j in seq(listlength)){
      # permutate
      runsubseq = subseq
      runsubseq[tmp[[j]]] = '[ATGCU]'
      runsubseq = paste(runsubseq,collapse='')
      dic = c(dic,runsubseq)
    }
  }
  #}
  return(dic)
}

# This function exhaust all posible motifs with mismatches, but not regular expression form ([ATGCU]CAC) like seqdic.
# Permute in ATGC order, the returned first one is the consensus motif
seqdic2 <- function(motif,n){
  seqlength = nchar(motif)
  if(n>seqlength){
    stop('Number of mismatch cannot larger than string length')
  }
  subseq = strsplit(motif,NULL)[[1]]
  #permute = function(subseq,seqlength,n){
  #allComb = lapply( 1:n, function(i) combn(seqlength,i,simplify = F) ) # All permutation combinations, from 1 mismatch to n mismatches
  dic =c()  # all permutated sequences
  if (n == 1){
    dic = as.vector(
      sapply( 1:nchar(motif),
              function(i) sapply(c('A','T','G','C'),
                                 function(s)
                                 {runsubseq = subseq;runsubseq[i] = s;runsubseq=paste(runsubseq,collapse='');return(runsubseq)}
              )
      )
    )
  }else{
    if(n == 0){
      dic = c(motif)
    }else{
      stop('Not support n!=1 yet')
    }
  }
  #}
  dic=dic[!dic==motif]
  dic = c(paste(subseq,collapse=''),unique(dic))
  return(dic)
}

#findmatch, give a sequence find the number and positions of mathes
findmatch <- function(seq,motif,n=0){
  if (n == 0){
    pois = unlist(as.vector(gregexpr(motif,seq)[[1]]))
  }else{
    # Get regular expression dictionary
    dic = seqdic(motif,n=n)
    pois = unique(unlist( lapply( dic,function(m) as.vector(gregexpr(m,seq)[[1]]) ) )) # match positions
  }
  if (-1 %in% pois){
    pois = sort(pois[-which(pois==-1)])
  }else{
    pois = sort(pois)
  }
  if (length(pois)==0){
    pois = ''
    count = 0
  }else{
    pois = sapply(pois,toString)
    count=toString(length(pois))
  }
  return( list(count=count,positions = pois ))
}

findmatch2 <- function(seq, motif, n=0, bidir = FALSE){
  library(Biostrings)
  motif <- DNAString(motif)
  seq <- DNAStringSet(seq)
  count <- vcountPattern(motif, seq, fixed = F)
  if (bidir){
    # birectional motif
    count <- count + vcountPattern(reverseComplement(motif), seq)
  }
  return(count)
}

# Compare two sequence, return base substitution (can understand as difference or motation) and substitution position
seq.compare <- function(seq1,seq2,subsplit=NULL){ # seq1 compared to seq2
  # subsplit. Provide if seq1 is shorter than seq2, used to truncate seq2 in order to compare with seq2. if provided, should be a region like c(1:3)
  seq1 = strsplit(seq1,NULL)[[1]]
  seq2 = strsplit(seq2,NULL)[[1]]
  # Remove . from seq1
  if('.' %in% seq1){
    subsplit = seq(length(seq2))[-which(seq1=='.')]
    if (1 %in% subsplit){
      shift = 0
    }else{
      shift = length(which(seq1=='.'))
    }
    seq1 = seq1[which(seq1!='.')]
    seq2 = seq2[subsplit]
  }else{
    shift = 0
  }
  #   if (length(seq1)<length(seq2)){
  #     if(!is.null(subsplit)){
  #       seq2 = seq2[subsplit]
  #       # check again
  #       if(length(seq1)!=length(seq2)){
  #         stop('Please provide correct subsplit length.')
  #       }
  #     }else{
  #       stop('seq1 is not equal to seq2, please provide subsplit to substr seq2 for comparison')
  #     }
  #   }else{
  #     if(length(seq1)>length(seq2)){
  #       stop('seq1 cannot be longer than seq2.')
  #     }
  #   }
  posi = c()
  substitution = c()
  for (i in seq(length(seq1))){
    if(seq1[i]==seq2[i]){
      next
    }else{
      substitution = c(substitution,seq1[i])
      posi = c(posi,i+shift)
    }
  }
  res = data.frame(substitution,posi)
  #names(res) = posi
  return(res)
  #return(list(positions = posi,substitution=substitution))
}

## Site matrix. return a matrix with each sites per row, and count of appearence for each motif allowing mismatches at column
sitematrix <- function(seqs,motif,n=1,extend=2){
  seqs = as.character(seqs)
  dic = seqdic(motif,n=n)
  positions = sapply( seqs, function(s) {r=unlist(sapply(dic, function(m) findmatch(s,m,n=0)$positions));r=unique(as.numeric(r[r!='']))-extend;return(r)} ) #yes n=0
  # get matched sequence
  matches = sapply( seq(length(positions)),function(p) unlist(sapply(positions[[p]],function(i) getseq(names(positions[p]),i,nchar(motif)+2*extend))) )
  sites = sapply(matches,length)
  # Find consensus sequence
  matches.as.matrix = matrix(unlist(sapply(unlist(matches),function(i) strsplit(i,NULL))), ncol=nchar(motif)+2*extend, byrow = TRUE)
  pwm = apply(matches.as.matrix,2,table)
  if(class(pwm)=='matrix'){
    consensus = paste( apply(pwm,2,function(i) names(which(i==max(i)))),collapse='' )
  }else{
    consensus = paste(sapply(pwm,function(i)names(i[which(i==max(i))])),collapse='')
  }
  # compare to consensus and get the sitematrix
  positions.X = rep(seq(nchar(consensus)),each=3) # positons of consensus kmer
  base.X = c(sapply(strsplit(consensus,NULL)[[1]],function(b) c('A','T','G','C')[which(c('A','T','G','C')!=b)]))
  run.compare = function(g){
    if(is.null(g)){
      return(c(rep(0,nchar(consensus)*3)))
    }else{
      res=data.frame()
      for(s in g){
        res=rbind(res,seq.compare(s,consensus))
      }
      res.v = sapply(seq(nchar(consensus)*3), function(i) sum(res$posi==positions.X[i] & res$substitution==base.X[i]) )
      return(res.v)
    }
  }
  X = t(sapply( matches, run.compare ))
  X = cbind(sites,X)
  return(list(X=X,pwm=pwm,base.X=base.X,positions.X=positions.X,consensus=consensus))
}


getseq <- function(seq,start,length,keeplength=T){
  if(start<=0){
    if(keeplength){
      append = paste(rep('.',1-start),collapse='')
      res = paste(append,substr(seq,start,start+length-1),sep='')
    }else{
      res = substr(seq,start,start+length-1)
    }
  }else{
    if((start+length-1)>nchar(seq)){
      if(keeplength){
        append = paste(rep('.',start+length-1-nchar(seq)),collapse='')
        res = paste(substr(seq,start,start+length-1),append,sep='')
      }else{
        res = substr(seq,start,start+length-1)
      }
    }else{
      res = substr(seq,start,start+length-1)
    }
  }
  return(res)
}

################################################
# Boxplot motif with rates, mismatches possible
################################################

boxplot_kmer <- function(kmer,n,seqcol='UTR3_seq',ratecol='dr',dat=UTRs,log=T,rate2time=FALSE,breaks=seq(1,30,5),
                         ylab='Degradation rates',xlab=NULL,main=NULL,origscale=T,omitlessthan=10, bidir = FALSE){
  require(ggplot2)
  # kill factor
  #dat <- data.frame(lapply(dat, as.character), stringsAsFactors=FALSE)
  dat <- data.frame(dat)
  # get counts
  #counts = as.numeric(unlist(sapply( dat[,seqcol], function(s) findmatch(s,kmer,n)$count )))
  counts <- findmatch2(dat[, seqcol], kmer, n, bidir = bidir)
  # combine boxes have less than omitlessthan members
  lessThanN = as.numeric(names(table(counts))[table(counts)<omitlessthan])[1]-2
  #lessThanN = c(lessThanN-1,lessThanN)
  #delIndx = which(counts %in% lessThanN)
  #if(length(delIndx)>0){
  #counts = counts[-delIndx]
  #}
  motifs <- cut(counts,c(-Inf,c(0:lessThanN),Inf))
  #tab = names(table(motifs))

  if (!is.null(main)){
    main = main
  }else{
    main=kmer #paste(gsub("T", "U", kmer),sep='\n') # paste(n,'mistmates',sep=' '),
  }
  if(rate2time){
    toplot = 1/as.numeric(dat[,ratecol])
  }else{
    toplot = as.numeric(dat[,ratecol])
  }
  #if(length(delIndx)>0){
  #toplot = toplot[-delIndx]
  #}

  # Significance test, count 0 vs others
  p=format.pval(wilcox.test(toplot ~ cut(counts,c(-Inf,0,Inf)))$p.value,digit=3)
  p = paste('p=',p,sep='')
  main = paste(main,p,sep='\n')

  # ggplot2
  g <- ggplot(data.frame(response=toplot, group=motifs), aes(factor(group),response))
  #g <- g + geom_boxplot(aes(fill=factor(group))) + labs(title=main,x=xlab,y=ylab)
  g <- g + geom_boxplot() + theme_bw() + labs(title=main,x=xlab,y=ylab)
  if(log){
    g <- g + scale_y_continuous(breaks=breaks,trans='log2')
  }
  g <- g + stat_summary(fun.data = give.n, geom = "text", fun.y = median,position = position_dodge(width = 0.75),size=rel(3))
  g <- g + theme(legend.position = "none",
                 axis.title = element_text(face="bold",size=rel(1.5)),
                 axis.text = element_text(face="bold",size=rel(1.2)),
                 plot.title = element_text(face="bold",size=rel(1.2), hjust = 0.5))
  g <- g + scale_x_discrete(labels=c(0:lessThanN,paste('>=',lessThanN+1)))
  return(g)
}

#################################################################
#  Compute fold change, between with and without certain motif
#################################################################
FC.motif <- function(kmer,n,seqcol='UTR3_seq',ratecol='dr',dat=UTRs,rate2time=FALSE,relative=T){
  # kill factor
  dat <- data.frame(lapply(dat, as.character), stringsAsFactors=FALSE)
  # get counts
  counts = as.numeric(unlist(sapply( dat[,seqcol], function(s) findmatch(s,kmer,n)$count )))
  tab = cut(counts,c(-Inf,0,Inf))

  if(rate2time){
    obj = 1/as.numeric(dat[,ratecol])
  }else{
    obj = as.numeric(dat[,ratecol])
  }

  p=wilcox.test(obj ~ cut(counts,c(-Inf,0,Inf)))$p.value
  with_motif = median(obj[tab=='(0, Inf]' & obj!='Inf'],na.rm=T)
  without_motif = median(obj[tab=='(-Inf,0]' & obj!='Inf'],na.rm=T)
  log2.fc = log2(with_motif/without_motif)
  relative_change = (with_motif - without_motif) / without_motif
  return(list(p=p,fc=log2.fc,relative_change=relative_change))
}



################################################
# Density plot of a given kmer
################################################
# kmers supply as a vector of kmers
# Need source function findmatch
# source('~/Documents/data/DTA_S.cerevisiae/Sc.Motifs/utilities.R')
density_kmer <- function(kmers,colors=NULL,dat=UTRs,seqcol='UTR3_seq',ratecol='dr',genenamecol='genename',
                         seqlengthcol='UTR3_length',window=20,max.length=200,align.right=T, promoter=F){
  # Check colors
  if (is.null(colors)){
    colors <- kmers
  }else{
    if (length(kmers)!=length(colors)){
      stop('Length of kmers not equal to length of color vector.')
    }
  }
  # kill factor
  dat <- data.frame(lapply(dat, as.character), stringsAsFactors=FALSE)
  rownames(dat) = dat[,genenamecol]
  # Compute statistics
  k_n_list = list()
  for (k in seq(length(kmers))){
    motif <- kmers[k]
    ## Create data frame with motif, gene, rates, UTR.length, dist
    # Find motif positions at each gene
    positions = sapply( 1:nrow(dat),function(i) findmatch(dat[,seqcol][i],n=0,motif)$positions )
    names(positions) = rownames(dat)
    # Find the index of genes which have target motif
    indx =c()
    for(i in 1:length(positions)){
      if (positions[[i]][1] != ''){
        indx = c(indx,i)
      }
    }
    tab = t(as.data.frame(sapply(rownames(dat)[indx],function(i)  # apply to all genes have the target motif
      sapply(positions[[i]], # For each gene, loop over all positions of kmers (some gene has more than one kmer)
             function(j)
               # position and other information of the motif, each motif instance one row
               c(i,motif,dat[i,ratecol],dat[i,seqlengthcol],
                 as.numeric(j)+nchar(motif)/2)  # center position of the motif
      )))
    )
    colnames(tab) = c('genename','motif','rates','length','dist')  # dist: is the distance w.r.t start of the sequence
    tab = as.data.frame(tab,stringsAsFactors=F)
    if(align.right){
      tab$dist = as.numeric(tab$length)-as.numeric(tab$dist)
    }

    ## Compute number of motifs at each positon and number of sequence coverage
    # normalization vector
    n=round(nchar(motif)/2,0)
    theo = c(rep(0,n),rep(1,max.length-2*n),rep(0,n)) # Theoretical posibale kmer numbers at each position
    theo.w = c()  # Theoretical kmer numbers at each window
    for(i in 1:(max.length)){
      theo.w = c(theo.w,sum(theo[max((i-(window/2)+1),1):(min((i+window/2),max.length))]))
    }
    # k: number of motifs, n: sequence coverage
    if(promoter){
      k_n = t(sapply( # number of motifs in the window
        1:max.length,
        function(i){
          c(sum(as.numeric(tab$dist) >= (i-window/2) & as.numeric(tab$dist) <= i+(window/2)),nrow(tab))
        }
      ))
    }else{
      k_n = t(sapply( # number of motifs in the window
        1:max.length, # slice one bp step
        function(i){
          c(sum(as.numeric(tab$dist) >= (i-window/2) & as.numeric(tab$dist) <= i+(window/2)), # fall into the 20bp window
            sum(as.numeric(tab$length) >= i+(window/2)))
        }
      ))
    }
    k_n = data.table(k_n,theo.w,1:length(theo.w),colors[k]) # the columns are: oberseved # of motifs in the window, number of genes that UTR/promter longer than the extend of the window, theoretical number of possible motifs
    colnames(k_n) <- c('k','n','theo.w','position','color')
    if (align.right){
      k_n$plot_position <- rev(k_n$position)
    }else{
      k_n$plot_position <- k_n$position
    }
    k_n_list[[motif]] = k_n
  }
  return(k_n_list)
}

## ggplot2 plot k-mer density
## Use custom colors
density_kmer_ggplot <- function(align.right=TRUE, ylab='Motif density', xlab='Distance from poly(A) site', window=20, max=150,...){
  require(data.table)
  k_n_list <- density_kmer(...)
  dat <- rbindlist(k_n_list, idcol='motif')
  dat[ ,density := window*(k/n/theo.w)]
  dat[ ,sd := sqrt(density*(1-density)/n)]
  dat <- dat[,.SD[1:max],by=motif]
  #dat[ ,position := rev(position)]
  density_plot <- ggplot(dat, aes(plot_position, density, ymax=density+sd, ymin=density-sd, color=color, fill=color)) +
    geom_smooth(stat='identity', alpha=0.2)
  if (!isTRUE(all.equal(dat$motif, dat$color))){
    density_plot <- density_plot + scale_fill_identity() + scale_color_identity()
  }
  if (align.right){
    breaks <- seq(0,max(dat$position),50)
    density_plot <- density_plot + scale_x_continuous(breaks=breaks,labels=rev(breaks))
  }
  density_plot <- density_plot +
    labs(x=xlab, y=ylab) +
    theme(axis.title=element_text(size=1.2, face='bold'))
  return(density_plot)
}

## Find percentage of two overlaping motifs co-appear
# motif2 should be the smaller motif which contained in motif1
# For example, find the overlap percentage of GCATAT and TTGCAT, motif1=TTGCATAT, motif2=GCATAT
coappear <- function(motif1,motif2,UTRs,seqcol='UTR3_seq'){
  a = sum(as.numeric(sapply(UTRs[,seqcol],function(s) findmatch(s,motif1,n=0)$count)))
  b = sum(as.numeric(sapply(UTRs[,seqcol],function(s) findmatch(s,motif2,n=0)$count)))
  print(a)
  print(b)
  percentage = a/b
  return(percentage)
}

# Two distinct motif, explore how many times they coappear
coappear2<- function(motif1,motif2,UTRs,seqcol1='UTR3_seq',seqcol2='UTR3_seq'){
  UTRs <- data.frame(UTRs)
  # Subset UTRs make sure no NA in the sequence column
  UTRs <- subset(UTRs, !is.na(UTRs[ ,seqcol1]) & !is.na(UTRs[ ,seqcol2]))
  a = as.numeric(sapply(UTRs[,seqcol1],function(s) findmatch(s,motif1,n=0)$count))>0
  b = as.numeric(sapply(UTRs[,seqcol2],function(s) findmatch(s,motif2,n=0)$count))>0
  co = sum(a&b)
  a.only = sum(a)-co
  b.only = sum(b)-co
  # Test of enrichment
  total <- dim(UTRs)[1]
  contingency_table <- matrix(c(co, b.only, a.only, total-a.only-b.only-co), ncol=2)
  test <- fisher.test(contingency_table)
  return(list(co=co,a_only=a.only,b_only=b.only,a=sum(a),b=sum(b),t=test))
}
##################################################
# Fit linear model
##################################################
lm_motif <- function(motif,n = 1,extend = 2,UTRs,ratecol='sr',seqcol='promoterseq',add_features=NULL,log=T){
  X=sitematrix(UTRs[,seqcol],motif=motif, n=1, extend=2)
  y=as.numeric(as.character(UTRs[,ratecol]))
  X.matrix = X$X[!is.na(y),]
  if(log){
    y = log(y[!is.na(y)])
  }else{
    y = y[!is.na(y)]
  }
  if (is.null(add_features)){
   fit = lm(y ~ X.matrix)
  }else{
    fit <- lm(y ~ cbind(X.matrix, data.matrix(UTRs[, add_features])))
  }
  return(list(fit=fit,X=X))
}

#################################################
# Plot coef of linear regression of motifs
#################################################
# input a coef vector, including intercept and site, they will be removed by the function
plot.coef <- function(fit,X=X,rates=T,main=X$consensus,ylim=c(0.8,2),ylab='effects on half-life'){
  require(ggplot2)
  # number of coef for motif single base
  n <- dim(X$X)[2]
  if(rates){
    coef = 1/exp(as.numeric(coef(fit))[3:(n+1)])
  }else{
    coef = exp(as.numeric(coef(fit))[3:(n+1)])
  }
  plotdat = data.frame(coef=coef,base=X$base.X,positions=X$positions.X,stringsAsFactors=F)
  c = strsplit(X$consensus,NULL)[[1]]
  t= data.frame(t(sapply(seq(length(c)),function(i) c(NA,c[i],i))),stringsAsFactors=F)
  colnames(t) = colnames(plotdat)
  t$coef=as.numeric(t$coef)
  t = rbind(t,plotdat)
  t$base = as.character(t$base)
  t$positions = as.numeric(t$positions)
  t$positions = as.factor(t$positions)
  d <- ggplot(data=t,aes(y=coef,x=positions,group=base,colour=base))  +
    geom_point(size=3) +
    geom_line(size=1.5) +
    ggtitle(main) +
    ylim(ylim) +
    theme_bw() +
    labs(y=ylab) +
    geom_hline(yintercept=1)
  return(d)
  #substitution = sapply(m, function(b) s=c('A','T','G','C') return(s[s!=b]) )
}


## Plot codon bias with decay
# useCount: use absolute count of codon used to calculate correlation, If FALSE, use Count/CDS_length to calculate correlation
# In my matrix, codons named like CGA.1 are count, names like CGA are relative percentage. This is the difault working mechnism of this function.
plot_certain_codon_rates <- function(codon,UTRs,ratecol='dr',ylim=c(0,0.3),n=T,log=T,ylab='Degradation rate',useCount=TRUE,mainadd=NULL){
  if(useCount){
    cod = codon
  }else{
    cod = substr(codon,1,3)
  }

  if(log){
    response = log(as.numeric(as.character(UTRs[,ratecol])))
  }else{
    response = as.numeric(as.character(UTRs[,ratecol]))
  }


  correlation = paste('cor=',as.character(round(cor(response,log(UTRs[,cod]+1),use='pairwise.complete.obs',m='s'),2)))
  sig = paste('pvalue=',format.pval(cor.test(response,log(UTRs[,cod]+1))$p.value,digit=3))
  main = paste( substr(codon,1,3),correlation,sig,sep=' ' )
  if(!is.null(mainadd)){
    main = paste(main,mainadd,sep='\n')
  }

  if(n){
    count = table(UTRs[,codon])
    boxplot(response ~ as.factor(UTRs[,codon]),ylim=ylim,xaxt='n',main=main,ylab=ylab)
    axis(1, at=1:length(levels(as.factor(UTRs[,codon]))), labels=paste(paste('n=',count,sep=''),0:(length(levels(as.factor(UTRs[,codon])))-1),sep=','),las=2 )
  }else{
    boxplot(response ~ as.factor(UTRs[,codon]),ylim=ylim)
  }
}


# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }

  if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


### get the matched sequence of findmatch function. This does not usefull for n=0 (no mismatch case), because the return always the input motif itself
motifMatch = function(seq,motif,n=1){
  indices = as.numeric(findmatch(seq,motif,n)$positions) # if no match, this return NA
  matches = sapply( indices, function(i) substr(seq,i,i+nchar(motif)-1) )
  return(matches)
}


## Get PWM matrix of a motif with mismatch in the given sequences
PWM.motif = function(seq,motif,n=1){
  matches = as.character(unlist(sapply(seq,function(i) motifMatch(seq=i,motif=motif,n=n))))
  # remove NA
  matches = matches[!is.na(matches)]
  # split into matrix
  matches.split = sapply(matches, function(i) strsplit(i,split=""))
  matches.matrix = matrix(unlist(matches.split),ncol=nchar(motif),byrow=T)
  pwm = sapply(1:nchar(motif),function(i) table(matches.matrix[,i]))
  return(pwm)
}


## generate seqlog for motif
seqlog.motif = function(seq,motif,n=1,ic.scale=T){
  require(seqLogo)
  pwm = PWM.motif(seq=seq,motif=motif,n=n)
  seqLogo(makePWM(pwm/colSums(pwm)),ic.scale=ic.scale)
}

## generate seqlog from motif pwm result from sitematrix function
seqlog.motif2 = function(sitematrix.pwm,ic.scale=T,...){
  #require(myseqLogo)
  # remove "."
  pwm.adjust = sapply( sitematrix.pwm, function(i) i=i[names(i)!="."] )
  seqLogo(pwm.adjust/colSums(pwm.adjust),ic.scale=ic.scale,...)
}


### Compute explained variance of linear model from the fit result
# method: "var", "mad"
variance_explainded_lm_fit <- function(fit, method='var'){
  y <- residuals(fit) + predict(fit)
  if(method == 'var'){
    # residual variance
    var_resid <- var(residuals(fit))
    # y response variance
    var_y <- var(y)
  }else{
    if(method == 'mad'){
      var_resid <- mad(residuals(fit))^2
      var_y <- mad(y)^2
    }
  }
  explained <- 1 - var_resid / var_y
  return(explained)
}


### Explained variance of PUF3
fit.motif = function(UTRs,ratecol='hlt.wt',seqcol='UTR3_seq',motif='TGTAAATA',log=T,rate=F,bootstrap=1000,...){
  UTRs <- data.frame(UTRs)
  motif.count = as.numeric(sapply(UTRs[,seqcol], function(i) findmatch(i,motif,n=0)$count))
  if(rate){
    response = 1/as.numeric(UTRs[,ratecol])
  }else{
    response = as.numeric(UTRs[,ratecol])
  }
  response = data.frame(response,motif.count)
  # if bootstraping
  if(bootstrap){
    rsquare = c()
    i=0
    while(i<=bootstrap){
      response.bootstrap = response[sample(nrow(response),replace=T),]
      if(log){
        fit = lm(log(as.numeric(response.bootstrap$response)) ~ response.bootstrap$motif.count)
      }else{
        fit = lm(as.numeric(response.bootstrap$response) ~ response.bootstrap$motif.count)
      }
      rsquare = c(rsquare, variance_explainded_lm_fit(fit, ...))
      i=i+1
    }
  }else{
    if(log){
      fit = lm(log(as.numeric(response$response)) ~ response$motif.count)
    }else{
      fit = lm(as.numeric(response$response) ~ response$motif.count)
    }
    rsquare = variance_explainded_lm_fit(fit, ...)
  }
  return(rsquare)
}

## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)

  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }

  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )

  # Rename the "mean" column
  datac <- rename(datac, c("mean" = measurevar))

  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean

  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval:
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult

  return(datac)
}


# Box plot with given order or with frequency
boxplotwithfreq = function(dat,plotc,lev=NULL,main='boxplotwithfreq',ylab='log2 half lives',ylim=c(0,8)){
  lv = dat[,plotc[2]]
  if(is.null(lev)){
    a=sort(table(lv))
    lv = factor(lv,names(a))
    boxplot(log2(dat[,plotc[1]]) ~ lv,las=2,xaxt='n',main=main,ylab=ylab,ylim=ylim)
    axis(1, at=1:length(a), labels=paste(names(a),a),las=2 )
  }else{
    lv = factor(lv,lev)
    boxplot(log2(dat[,plotc[1]]) ~ lv,las=2,main=main,ylab=ylab,ylim=ylim)
  }
}


# boxplot, ommit boxes with less than n instances
boxplotOmmitN <- function(x,y,n=10,log2x=TRUE,...){
  if(length(x) != length(y)){
    stop('x, y length differ!')
  }
  counts = table(y)
  lessThanN = names(table(y))[table(y)<n]
  delIndx = which(as.character(y) %in% lessThanN)
  if (log2x){
    boxplot(log2(x[-delIndx]) ~ y[-delIndx],... )
  }else{
    boxplot(x[-delIndx] ~ y[-delIndx],... )
  }
}

# Count codon number, without in frame
count_triplet <- function(triplet,seq){
  posi = unlist(gregexpr(triplet,seq))
  if (-1 %in% posi){
    posi = sort(posi[-which(posi==-1)])
  }
  num = length(posi)
  return(list(posi=posi,num=num))
}


# codon usage vector. This function return a vector of codon usage counts
# '~/Documents/data/DTA_S.cerevisiae/Model/Features/Codon/codons.RData'
# function width and step defined similar to oligonucleotidefrequency function from Biostrings
codonusagevector <- function(sequence,codons2count,normalize=T,width=3,step=3){
  if(is.na(sequence)){
    countcodon=rep(0,length(codons2count))
    names(countcodon) = codons2count
  }else{
    # split sequence into triplet
    #n <- nchar(codons2count[1])
    starts <- seq(1, nchar(sequence)-width+1, step)
    ends <- starts + width - 1
    usedcodons <- substring(sequence, starts, ends)
    usedcodons <- table(usedcodons)
    countcodon <- sapply(codons2count, function(x) {
      counts = usedcodons[x]
      if(is.na(counts)) counts = 0;
      return(counts)
    }
    )
    names(countcodon) = codons2count
    if(normalize){
      countcodon = countcodon / sum(countcodon)
    }
  }
  return(countcodon)
}


##' Concatenate the spliced sequence of a single transcript to one single sequence
##' Need a data.frame(data.table) with two columns, the id column and the sequence column
concatenate_seq_per_transcript <- function(dat,idcol='transcript',seqcol='CDS_seq'){
  if(is.data.table(dat)){
    select_dat <- dat[,c(idcol,seqcol),with=FALSE]
  }else{
    select_dat <- data.table(dat[,c(idcol,seqcol)])
  }
  select_dat <- select_dat[,paste(get(seqcol),collapse=''),by=c(idcol)]
  setnames(select_dat, c(idcol,seqcol))
  return(select_dat)
}

## ggplot2 related
give.n <- function(x){
  return(c(y = median(x)*1.05, label = length(x)))
  # experiment with the multiplier to find the perfect position
}


boxplotOmmitN <- function(x,y,n=10,breaks=seq(1,30,10),log=TRUE,main='title',fill=NULL,
                          xlab='group',ylab='half-life [min]',las=1,n.vjust=2,ylim=NULL,
                          axis.title.size=2,axis.text.size=2,plot.title.size=2,give.n.size=3,give_n=TRUE,violin=FALSE){
  require(ggplot2)
  give.n <- function(x){
    return(c(y = median(x)*1.05, label = length(x)))
    # experiment with the multiplier to find the perfect position
  }

  if(length(x) != length(y)){
    stop('x, y length differ!')
  }
  if(breaks[1]==0 & log){
    stop('If log transform, breaks should not contain 0.')
  }

  if(is.numeric(y)){
    lessThanN = as.numeric(names(table(y))[table(y)<n])[1]
    if(!is.na(lessThanN)){
      groups <- cut(y,c(-Inf,c(0:(lessThanN-1)),Inf))
    }else{
      groups <- y
    }
  }else{
    counts = table(y)
    lessThanN = names(table(y))[table(y)<n]
    delIndx = which(as.character(y) %in% lessThanN)
    if(length(delIndx)!=0){
      groups = y[-delIndx]
      x = x[-delIndx]
    }else{
      groups = y
    }
  }

  # ggplot2
  if (is.null(fill)){
    df <- data.frame(response=x, group=groups, fill=groups)
  }else{
    df <- data.frame(response=x, group=groups, fill=fill)
  }

  df <- df[complete.cases(df),]
  g <- ggplot(df, aes(x=factor(group),y=response,fill=fill)) + theme_bw()
  if (violin){
    g <- g + geom_violin() +
      geom_boxplot(width=0.08, fill='black', outlier.shape = NA) +
      stat_summary(fun.y='median', color = 'white', size=1, geom='point') +
      labs(title=main,x=xlab,y=ylab)
  }else{
    g <- g + geom_boxplot(outlier.size = 0.5) + labs(title=main,x=xlab,y=ylab)
  }

  if(log){
    if(!is.null(ylim)){
      g <- g + scale_y_continuous(breaks=breaks,trans='log10',limits=ylim)
    }else{
      g <- g + scale_y_continuous(breaks=breaks,trans='log10',limits=NA)
    }
  }
  if(give_n){
    g <- g + stat_summary(fun.data = give.n, geom = "text", vjust=n.vjust, size=rel(give.n.size)) # position = position_dodge(width = 0.75),
  }
  g <- g + theme(legend.position = "none",
                 axis.title = element_text(face="bold",size=rel(axis.title.size)),
                 axis.text = element_text(face="bold",size=rel(axis.text.size)),
                 plot.title = element_text(face="bold",size=rel(plot.title.size)),
                 panel.grid.minor = element_blank())
  if(las==2){
    g <- g + theme(axis.text.x=element_text(angle=90,hjust=0))
  }
  #   if(!is.null(ylim)){
  #     g <- g + coord_cartesian(ylim = ylim)
  #   }
  if(is.numeric(y)){
    if(!is.na(lessThanN))
      #g <- g + scale_x_discrete(labels=c(0:(lessThanN-1),paste0('\u2265',lessThanN)))
      g <- g + scale_x_discrete(labels=c(0:(lessThanN-1),paste0('>=',lessThanN)))
  }
  return(g)
}

#### Codon and tAIs
#load('~/Documents/data/DTA_S.cerevisiae/Model/Features/Codon/tAI_function.RData')
#source('/data/ouga/home/ag_gagneur/chengju/Documents/data/DTA_S.cerevisiae/Model/Features/Codon/tAI_function.R')

count_triplet_inframe <- function(triplet,sequence,rev=FALSE,selectfirst=NULL){
  require('IRanges')
  if(rev){
    sequence = reverse(sequence)
  }
  if(!is.null(selectfirst)) sequence = substr(sequence,1,selectfirst)
  if(is.na(sequence)){
    count = 0
  }else{
    if (nchar(sequence) > 2){
      sequence = sapply(seq(1,nchar(sequence),3), function(i) substr(sequence, i, i+2))
      count = table(sequence)
      count = count[triplet]
      if(is.na(count)){
        count=0
      }
    }else{
      count = 0
    }
  }
  return(count)
}

count_triplet <- function(triplet,seq,selectlast=NULL){
  if(!is.null(selectlast)){
    if(nchar(seq) > selectlast) seq = substr(seq,nchar(seq)-selectlast+1,nchar(seq))
  }
  posi = unlist(gregexpr(triplet,seq))
  if (-1 %in% posi){
    posi = sort(posi[-which(posi==-1)])
  }
  # remove NA
  posi = sort(posi)
  num = length(posi)
  return(list(posi=posi,num=num))
}


# Calculate GC content from given sequence vectors
# Only consider sequences longer than min_length
GC_content <- function(sequences,min_length=20){
  #
  require(Biostrings)
  # Need to remove NAs
  gc_content <- sapply(sequences, function(s) {
    if(is.na(s) | nchar(s)<min_length){
      return(NA)
    }else{
      return(letterFrequency(DNAString(s),'GC'))
    }
  })
  gc_content <- as.numeric(gc_content) / nchar(sequences)
  return(gc_content)
}

# GC content of third codon positon
GC3_content <- function(sequence){
  if(is.na(sequence)){
    return(NA)
  }else{
    # check whether is a ORF
    if(nchar(sequence) %% 3 != 0) {
      stop("Sequence not a open reading frame")
    }else{
      # Extract third codon poisition
      third_position <- paste0(strsplit(sequence,'')[[1]][(seq(nchar(sequence))%%3)==0],collapse='')
      gc_content <- letterFrequency(DNAString(third_position),'GC') / nchar(third_position)
      return(gc_content)
    }
  }
}

## Variance explained by each predictor from multivariate regression model
variance_explained <- function(model_fit){
  variable_summary <- anova(model_fit)
  variable_ss <- variable_summary$"Sum Sq"
  return(cbind(variable_summary,PctVarExp=variable_ss/sum(variable_ss)*100))
}

concatenate_seq_per_transcript <- function(dat,idcol='transcript',seqcol='CDS_seq'){
  if(is.data.table(dat)){
    select_dat <- dat[,c(idcol,seqcol),with=FALSE]
  }else{
    select_dat <- data.table(dat[,c(idcol,seqcol)])
  }
  select_dat <- select_dat[,paste(get(seqcol),collapse=''),by=c(idcol)]
  setnames(select_dat, c(idcol,seqcol))
  return(select_dat)
}


## Function calculate uAUG status
getcod <- function(UTR5,cod='ATG',reverse=T){
  # reverse is true for UTR5, but not for cds or UTR3
  if(reverse){
    cod_rev <- Biostrings::reverse(cod)
  }else{
    cod_rev <- cod
  }
  num_cod <- as.numeric(sapply(UTR5, function(x) count_triplet(cod,x)$num))
  num_cod_inframe <- sapply(UTR5, function(x) count_triplet_inframe(cod_rev,x,rev = T))
  cod_status <- ifelse(num_cod,num_cod,paste0('no ',cod))
  cod_status <- ifelse(cod_status==num_cod_inframe,paste0('Only_inframe_',cod),cod_status)
  cod_status <- ifelse(grepl(cod,cod_status),cod_status,paste0('w.o.frame ',cod))
  return(data.table(num_cod=num_cod,status=cod_status))
}

get_uAUG_status <- function(UTRs, UTR5seqcol='UTR5_seq', CDSseqcol='CDS_seq', cod='ATG'){
  # Four categories:
  # No PTC
  # 1. no uAUG
  # 2. In-frame uAUG w/o PTC
  # 3. Out-frame uAUG w/o PTC
  # 3. uORF
  # 4. Out-frame uAUG w PTC
  # 5. In-frame uAUG with PTC in CDS
  require(Biostrings)
  UTRs <- data.frame(UTRs)
  UTRs <- subset(UTRs, !is.na(UTRs[ ,UTR5seqcol]))
  UTR5_seq <- UTRs[ ,UTR5seqcol]
  CDSseq <- UTRs[ ,CDSseqcol]
  all_uORF <- lapply(UTR5_seq, findORFsinSeq)
  uORF_counts <- sapply(all_uORF, function(i) length(i[[1]]))
  num_uAUG_inframe <- sapply(UTR5_seq, function(x) count_triplet_inframe(Biostrings::reverse(cod),x,rev = T))
  in_out_frame <- getcod(UTR5_seq, cod=cod)
  uAUG_counts <- in_out_frame[ ,num_cod]
  in_out_frame <- in_out_frame[ ,status]
  uAUG_status <- ifelse(in_out_frame=='no ATG', 'no ATG', ifelse(uORF_counts>0, 'uORF', in_out_frame))
  # next dissect with out-frame uAUG to PTC_in_CDS and no PTC
  # take out those transcripts
  out_frame_uAUG_transcripts <- paste0(UTR5_seq, CDSseq)[uAUG_status=='w.o.frame ATG']
  with_PTC_in_CDS <- lapply(out_frame_uAUG_transcripts, findORFsinSeq)
  # check whether the starts of this ORFs are in the 5'UTR
  # some transcripts may have out-frame uAUG, but at the same time have in-frame uAUG form a ORF with stop codon in CDS (Not true, cannot have in-frame stop in CDS)
  # check whether positions of with_PTC_in_CDS starts are in the 5'UTR
  out_frame_uAUG_transcripts_5UTR_length <- nchar(UTR5_seq)[uAUG_status=='w.o.frame ATG']
  with_PTC_in_CDS <- sapply(seq(length(with_PTC_in_CDS)), function(i) any(with_PTC_in_CDS[[i]][['starts']] < out_frame_uAUG_transcripts_5UTR_length[i]))
  with_PTC_in_CDS <- which(uAUG_status=='w.o.frame ATG')[with_PTC_in_CDS]
  uAUG_status[with_PTC_in_CDS] = 'Out_frame_PTC_in_CDS'
  uAUG_status
}

## Perform hypothesis test of each element agaist the reference, input a data table
cross_test_dt <- function(dt,test='wilcox.test',idvar='.id',ref='WILD',measurevar='measure',adjust=T,adj_m='BH'){
  if(!is.data.frame(dt)){
    stop('Need input a data.table')
  }
  dt <- data.table(dt)
  setkeyv(dt,idvar)
  ids <- unique(dt[,get(idvar)])
  pvars <- c()
  pvars <- sapply(ids, function(s) get(test)(dt[ref,get(measurevar)],dt[s,get(measurevar)])$p.value )
  pvars <- pvars[names(pvars)!=ref]
  if(adjust){
    pvars <- p.adjust(pvars,adj_m)
  }
  pvars[ref]=1
  return(pvars)
}


## Convert a vector of pvalues to characters, p<0.1 ., p<0.05 **, p<0.01 ***, the others are empty string
pvars2star <- function(pvars){
  pvars <- ifelse(pvars<0.001,'***',ifelse(pvars<0.01,'**',ifelse(pvars<0.05,'*',ifelse(pvars<0.1,'',''))))
  return(pvars)
}

fdr2star <- function(fdrs, alpha=0.1){
  fdrs <- ifelse(fdrs<alpha, "*", "")
  fdrs
}

#############################################################################
#
# Analyse variance explained by each predictor
# select variable in increase oder, determine Rsquare, adjusted Rsquare and BIC at each step
# will take all columns except the response as features
#
#############################################################################

### Rsqaure from cross validation
### Calculate by square of pearson correlation
cv_rsquare <- function(data,response,formular,fold=10){
  data <- data.table(data)
  fold <- fold
  n <- nrow(data)
  myFolds <- split(sample(n),rep(1:fold,length=n))
  predicted <- c()
  for (ii in myFolds) {
    testset <- data[ii,]
    trainset <- data[-ii,]
    fit <- lm(data=trainset, as.formula(formular))
    y_hat <- predict(fit,testset)
    names(y_hat) <- ii
    predicted <- c(predicted,y_hat)
  }
  y <- data[,get(response)]
  names(y) <- 1:length(y)
  predicted <- predicted[names(y)]
  rsqaure <- cor(predicted,y)
  rsqaure <- rsqaure^2
  return(rsqaure)
}


variance_explainded_lm <- function(data,response){
  if(!is.data.frame(data)){
    stop('Input should be data frame!')
  }
  features <- colnames(data)
  features <- features[-which(features==response)]
  generate_formular <- function(response, predictors){
    paste(response, '~', paste(predictors,collapse = '+'))
  }
  ## First fit a model with each predictor, rank the AIC
  fit_indv <- sapply(features, function(f) {
    fited <- lm(data=data,as.formula(generate_formular(response,f)))
    aic <- extractAIC(fited)[2]
    rsquare <- summary(fited)$r.squared
    adj_rsquare <- summary(fited)$adj.r.squared
    return(c(AIC=aic,rsquare=rsquare,adj_rsquare=adj_rsquare))
  })
  fit_indv <- data.table(t(fit_indv),keep.rownames = T)[order(AIC)]
  fit_indv_predictors <- fit_indv[,rn]

  ## Evaluate them in a forward session
  #bAIC <- fit_indv[1,]
  initial_model <- lm(data=data,as.formula(generate_formular(response,fit_indv_predictors[1])))
  aic <- c(fit_indv[1,AIC])
  rsquare <- c(summary(initial_model)$r.squared)
  adj_rsquare <- c(summary(initial_model)$adj.r.squared)
  form <- paste('~',fit_indv_predictors[1])
  brsquare_cv <- cv_rsquare(data, response,paste(response,form)) # first rsquare_cv
  rsquare_cv <- c(brsquare_cv)
  rsquare_cv_gain <- c(brsquare_cv)
  formular <- c(form)
  feature_add <- c(fit_indv_predictors[1])
  for(f in fit_indv_predictors[2:length(fit_indv_predictors)]){
    form <- paste(form,'+',f)
    arsquare_cv <- cv_rsquare(data,response,paste(response,form))
    rsquare_cv <- c(rsquare_cv, arsquare_cv)
    rsquare_cv_gain <- c(rsquare_cv_gain,arsquare_cv-brsquare_cv)
    fited <- update(initial_model, as.formula(form))
    formular <- c(formular,form)
    aic <- c(aic,extractAIC(fited)[2])
    rsquare <- c(rsquare, summary(fited)$r.squared)
    adj_rsquare <- c(adj_rsquare,summary(fited)$adj.r.squared)
    feature_add <- c(feature_add,f)
    brsquare_cv <- arsquare_cv # arsquare_cv become brsquare_cv
  }

  stat <- data.table(formular=formular,
                     feature_add=feature_add,
                     AIC=aic,rsquare=rsquare,
                     adj_rsquare=adj_rsquare,
                     rsquare_cv=rsquare_cv,
                     rsquare_cv_gain=rsquare_cv_gain)
  stat

}



analysis_model <- function(data,response){
  require(ggplot2)
  variance_explain = variance_explainded_lm(data,response)
  plotdata <- melt(variance_explain,id.vars = c('formular','feature_add'))
  plotdata$feature_add <- factor(plotdata$feature_add, levels = variance_explain[,feature_add])
  p1 <- ggplot(plotdata[variable!='AIC'],aes(x=feature_add,y=value,group=variable,color=variable)) + geom_line() + theme(axis.text.x=element_text(angle=-45,hjust=0,size=15))
  p2 <- ggplot(plotdata[variable=='AIC'],aes(x=feature_add,y=value,group=variable,color=variable)) + geom_line() + theme(axis.text.x=element_text(angle=-45,hjust=0,size=15))
  print(p1)
  print(p2)
}

##########################################################################################
## This function intend to analysis the variance_explained with backward selection procedure
##########################################################################################
# function to update model by drop one variable
dropvar <- function(model, dropvar){
  droped_model <- update(model, as.formula(paste(".~.-",dropvar)))
  return(droped_model)
}

generate_formular <- function(response, predictors){
  paste(response, '~', paste(predictors,collapse = '+'))
}

## This function calculate cross validated rsquare with a given lm object
cv_rsquare2 <- function(lmfit, fold=10, return_predict = FALSE){
  data <- lmfit$model
  fold <- fold
  n <- nrow(data)
  myFolds <- split(sample(n),rep(1:fold,length=n))
  predicted <- c()
  for (ii in myFolds) {
    testset <- data[ii,]
    trainset <- data[-ii,]
    fit <- lm(data=trainset, as.formula(lmfit$call))
    y_hat <- predict(fit,testset)
    names(y_hat) <- ii
    predicted <- c(predicted,y_hat)
  }
  y <- data[,1]
  names(y) <- 1:length(y)
  predicted <- predicted[names(y)]
  rsqaure <- cor(predicted,y)
  rsqaure <- rsqaure^2
  if (return_predict){
    return(data.table(y = y, y_hat = predicted))
  }
  return(rsqaure)
}


variance_explainded_lm_backward <- function(model_data,response){
  if(!is.data.frame(model_data)){
    stop('Input should be data frame!')
  }
  features <- colnames(model_data)
  features <- features[-which(features==response)]
  ## First fit a model with all predictors
  full_model_formular <- as.formula(generate_formular(response,features))
  full_model <- lm(data=model_data,full_model_formular)

  ## Evaluate them in a forward session
  # test to drop each variable, see which one have the largest influence on rsquare
  model <- full_model
  feature_drop <- c('FullModel')
  aic <- c(AIC(full_model))
  rsquare <- c(summary(full_model)$r.squared)
  adj_rsquare <- c(summary(full_model)$adj.r.squared)
  rsquare_cv <- c(cv_rsquare(data=model_data,response=response,formular=full_model_formular))
  rsquare_cv_drop <- c(0)
  while(length(features)>1){
    reduced_model <- most_important_var(model)
    model <- reduced_model$model
    rsquare_cv_new <- cv_rsquare2(model)
    rsquare_cv_drop <- c(rsquare_cv_drop,rsquare_cv[length(rsquare_cv)]-rsquare_cv_new)
    rsquare_cv <- c(rsquare_cv,rsquare_cv_new)
    rsquare <- c(rsquare,summary(model)$r.squared)
    adj_rsquare <- c(adj_rsquare,summary(model)$adj.r.squared)
    aic <- c(aic,AIC(model))
    feature_drop <- c(feature_drop,reduced_model$feature)
    features <- attr(model$terms,'term.labels')
  }

  stat <- data.table(feature_drop=feature_drop,
                     AIC=aic,rsquare=rsquare,
                     adj_rsquare=adj_rsquare,
                     rsquare_cv=rsquare_cv,
                     rsquare_cv_drop=rsquare_cv_drop)
  stat
}

most_important_var <- function(model){
  features <- attr(model$terms,'term.labels')
  if (length(features)>1){
    droped_models <- lapply(features, function(v) dropvar(model,v) )
    rsquared_droped_models <- sapply(droped_models, function(m) cv_rsquare2(m))
    names(rsquared_droped_models) <- features
    most_drop <- which(rsquared_droped_models == min(rsquared_droped_models))
    #sprintf('droped %s', features[most_drop])
    # myfunc(droped_models[[most_drop]])
    return(list(feature=features[most_drop],rsquare=rsquared_droped_models[most_drop],model=droped_models[[most_drop]]))
  }
}

analysis_model_backward <- function(model_data,response){
  require(ggplot2)
  variance_explain = variance_explainded_lm_backward(model_data,response)
  plotdata <- melt(variance_explain,id.vars = 'feature_drop')
  plotdata$feature_drop <- factor(plotdata$feature_drop, levels = variance_explain[,feature_drop])
  p1 <- ggplot(plotdata[variable!='AIC'],aes(x=feature_drop,y=value,group=variable,color=variable)) + geom_line() + theme(axis.text.x=element_text(angle=-45,hjust=0,size=15))
  p2 <- ggplot(plotdata[variable=='AIC'],aes(x=feature_drop,y=value,group=variable,color=variable)) + geom_line() + theme(axis.text.x=element_text(angle=-45,hjust=0,size=15))
  print(p1)
  print(p2)
}

###########################
# Functions related to tAI
###########################
cod2tAI = function(cod,tAIs){
  return(tAIs[cod])
}
#cod2tAI('TTT',tAIs)


# Calculate geometric mean
gm <- function(x,na.rm=T){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

# This function return: 1) a matrix of tAI of all genes, all codons
#                       2) a matrix of tAI of all genes, first and last p percent of codons
#                       3) a vector of average tAI of all genes averaged from all codons
#                       4) a vector of average tAI of all genes averaged from first and last p percent of codons
tAImatrix = function(CDSseq,tAIs,percentage=10){
  m = sapply(CDSseq, function(i) substring(i,seq(1,nchar(i)-2,3),seq(3,nchar(i),3)))
  tai = sapply(m, function(i) as.numeric(sapply(i,function(k) cod2tAI(k,tAIs))) )
  avg_gene_tAI = as.numeric(sapply(tai, gm))
  #names(avg_gene_tAI) = UTRs$genename
  # Select first or last p percent of codons
  ncod = as.numeric(sapply(tai,length))
  ncod2select = round(as.numeric(sapply(tai,length))*(percentage/100))
  selected_start = sapply(seq(length(CDSseq)), function(i) tai[[i]][1:ncod2select[i]])
  avg_selected_start_tAI = as.numeric(sapply(selected_start,gm))
  #names(avg_selected_start_tAI) = UTRs$genename
  selected_end = sapply(seq(length(CDSseq)), function(i) tai[[i]][(ncod[i]-ncod2select[i]-1):ncod[i]])
  avg_selected_end_tAI = as.numeric(sapply(selected_end,gm))
  #names(avg_selected_end_tAI) = UTRs$genename
  return(list(tAI=tai,
              avg_gene_tAI=avg_gene_tAI,
              selected_start = selected_start,
              avg_selected_start_tAI=avg_selected_start_tAI,
              selected_end = selected_end,
              avg_selected_end_tAI=avg_selected_end_tAI,
              ncod2select=ncod2select))
}



################################
# start codon context sequence
################################
start_context_seq <- function(UTRs,UTR5seqcol='UTR5_seq',CDSseqcol='CDSseq',l=6,r=6){
  # get sequence
  require(data.table)
  UTRs <- data.table(UTRs)
  seqs <- UTRs[,c(UTR5seqcol,CDSseqcol),with=F]
  seqs <- seqs[nchar(get(UTR5seqcol))>l]
  kozaks <- paste0(seqs[,substr(get(UTR5seqcol),nchar(get(UTR5seqcol))-l+1,nchar(get(UTR5seqcol)))],seqs[,substr(get(CDSseqcol),1,l+3)])
  return(kozaks)
}


################################
# stop codon context sequence
################################
stop_context_seq <- function(UTRs,UTR3seqcol='UTR3_seq',CDSseqcol='CDSseq',stopplus1='stop1',l=6,r=6){
  # get sequence
  require(data.table)
  UTRs <- data.table(UTRs)
  seqs <- UTRs[,c(UTR3seqcol,CDSseqcol,stopplus1),with=F]
  seqs <- seqs[nchar(get(UTR3seqcol))>r & nchar(get(CDSseqcol))>l]
  stopcontext <- paste0(seqs[,substr(get(CDSseqcol),nchar(get(CDSseqcol))-l+1,nchar(get(CDSseqcol)))],
                        seqs[,substr(get(stopplus1),1,3)],
                        seqs[,substr(get(UTR3seqcol),1,l)])
  stopcodon <- seqs[,substr(get(stopplus1),1,3)]
  return(data.table(stopcodon=stopcodon,stopcontext=stopcontext))
}

stop_context_seq2 <- function(UTRs,UTR3seqcol='UTR3_seq',CDSseqcol='CDSseq',l=6,r=6){
  # get sequence
  require(data.table)
  UTRs <- data.table(UTRs)
  seqs <- UTRs[,c(UTR3seqcol,CDSseqcol),with=F]
  seqs <- seqs[nchar(get(UTR3seqcol))>r & nchar(get(CDSseqcol))>l]
  stopcontext <- paste0(seqs[,substr(get(CDSseqcol),nchar(get(CDSseqcol))-l-2,nchar(get(CDSseqcol)))],
                        seqs[,substr(get(UTR3seqcol),1,l)])
  stopcodon <- seqs[,substr(get(CDSseqcol),nchar(get(CDSseqcol))-2,nchar(get(CDSseqcol)))]
  return(data.table(stopcodon=stopcodon,stopcontext=stopcontext))
}



########################
#  glmnet functions
########################
cv_glmnet_prediction <- function(X,y,fold=10,alpha=0,lambda=0){
  fold = 10
  n = length(y)
  myFolds  = split(sample(n),rep(1:fold,length=n))
  predicted = c()
  for (ii in myFolds) {
    testset <- X[ii,]
    trainset <- X[-ii,]
    y_testset <- y[ii]
    y_trainset <- y[-ii]
    fit = glmnet(trainset,y_trainset,alpha=alpha,lambda=lambda)

    y_hat = predict(fit,newx=testset,s=lambda)
    names(y_hat) = ii
    predicted = c(predicted,y_hat)
  }
  predicted <- predicted[as.character(1:length(y))]
}

cv_glmnet_rquared <- function(X,y,fold=10,alpha=0,lambda=0,repeats=1){
  i <- 0
  collection <- c()
  while(i<repeats){
    predicted <- cv_glmnet_prediction(X,y,fold,alpha,lambda)
    rsquared <- cor(y,predicted)^2
    collection <- c(collection,rsquared)
    i <- i+1
  }
  collection
}

plot_coef_glmnet <- function(glmnetfit,s='lambda.min'){
  barplot(as.data.frame(as.matrix(coef(glmnetfit,s=s)))[,1][-1],names=rownames(as.data.frame(as.matrix(coef(glmnetfit))))[-1],las=2)
}

# list(plot, row(s), column(s))
lay_out = function(...) {
  x <- list(...)
  n <- max(sapply(x, function(x) max(x[[2]])))
  p <- max(sapply(x, function(x) max(x[[3]])))
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(n, p)))
  for (i in seq_len(length(x))) {
    print(x[[i]][[1]], vp = grid::viewport(layout.pos.row = x[[i]][[2]],
                                           layout.pos.col = x[[i]][[3]]))
  }
}


## Find ORFs
## From  Avril Coghlan, Wellcome Trust Sanger Institute, Cambridge, U.K
## http://a-little-book-of-r-for-bioinformatics.readthedocs.io/en/latest/src/chapter7.html
## Under a Creative Commons Attribution 3.0 Licensehttps://creativecommons.org/licenses/by/3.0/
findPotentialStartsAndStops <- function(sequence)
{
  # Define a vector with the sequences of potential start and stop codons
  codons            <- c("ATG", "TAA", "TAG", "TGA")
  # Find the number of occurrences of each type of potential start or stop codon
  for (i in 1:4)
  {
    codon <- codons[i]
    # Find all occurrences of codon "codon" in sequence "sequence"
    occurrences <- matchPattern(codon, sequence)
    # Find the start positions of all occurrences of "codon" in sequence "sequence"
    codonpositions <- start(occurrences) # attr(occurrences,"start")
    # Find the total number of potential start and stop codons in sequence "sequence"
    numoccurrences <- length(codonpositions)
    if (i == 1)
    {
      # Make a copy of vector "codonpositions" called "positions"
      positions <- codonpositions
      # Make a vector "types" containing "numoccurrences" copies of "codon"
      types <- rep(codon, numoccurrences)
    }
    else
    {
      # Add the vector "codonpositions" to the end of vector "positions":
      positions   <- append(positions, codonpositions, after=length(positions))
      # Add the vector "rep(codon, numoccurrences)" to the end of vector "types":
      types       <- append(types, rep(codon, numoccurrences), after=length(types))
    }
  }
  # Sort the vectors "positions" and "types" in order of position along the input sequence:
  indices <- order(positions)
  positions <- positions[indices]
  types <- types[indices]
  # Return a list variable including vectors "positions" and "types":
  mylist <- list(positions,types)
  return(mylist)
}

findORFsinSeq <- function(sequence)
{
  require(Biostrings)
  if (is.na(sequence)){
    orfstarts = NA
    orfstops = NA
    orflengths = NA
  }else{
    # Make vectors "positions" and "types" containing information on the positions of ATGs in the sequence:
    mylist <- findPotentialStartsAndStops(sequence)
    positions <- mylist[[1]]
    types <- mylist[[2]]
    # Make vectors "orfstarts" and "orfstops" to store the predicted start and stop codons of ORFs
    orfstarts <- numeric()
    orfstops <- numeric()
    # Make a vector "orflengths" to store the lengths of the ORFs
    orflengths <- numeric()
    # Print out the positions of ORFs in the sequence:
    # Find the length of vector "positions"
    numpositions <- length(positions)
    # There must be at least one start codon and one stop codon to have an ORF.
    if (numpositions >= 2)
    {
      for (i in 1:(numpositions-1))
      {
        posi <- positions[i]
        typei <- types[i]
        found <- 0
        while (found == 0)
        {
          for (j in (i+1):numpositions)
          {
            posj  <- positions[j]
            typej <- types[j]
            posdiff <- posj - posi
            posdiffmod3 <- posdiff %% 3
            # Add in the length of the stop codon
            orflength <- posj - posi + 3
            if (typei == "ATG" && (typej == "TAA" || typej == "TAG" || typej == "TGA") && posdiffmod3 == 0)
            {
              # Check if we have already used the stop codon at posj+2 in an ORF
              numorfs <- length(orfstops)
              usedstop <- -1
              if (numorfs > 0)
              {
                for (k in 1:numorfs)
                {
                  orfstopk <- orfstops[k]
                  if (orfstopk == (posj + 2)) { usedstop <- 1 }
                }
              }
              if (usedstop == -1)
              {
                orfstarts <- append(orfstarts, posi, after=length(orfstarts))
                orfstops <- append(orfstops, posj+2, after=length(orfstops)) # Including the stop codon.
                orflengths <- append(orflengths, orflength, after=length(orflengths))
              }
              found <- 1
              break
            }
            if (j == numpositions) { found <- 1 }
          }
        }
      }
    }
    # Sort the final ORFs by start position:
    indices <- order(orfstarts)
    orfstarts <- orfstarts[indices]
    orfstops <- orfstops[indices]
    # Find the lengths of the ORFs that we have
    orflengths <- numeric()
    numorfs <- length(orfstarts)
    for (i in 1:numorfs)
    {
      orfstart <- orfstarts[i]
      orfstop <- orfstops[i]
      orflength <- orfstop - orfstart + 1
      orflengths <- append(orflengths,orflength,after=length(orflengths))
    }
  }
  mylist <- list(starts=orfstarts, stops=orfstops, lengths=orflengths)
  return(mylist)
}

findStartCodon <- function(sequence){
  require(Biostrings)
  if(is.na(sequence)){
    starts <- NA
  }else{
    all_atgs <- matchPattern('ATG',sequence)
    starts <- start(all_atgs)
  }
  starts
}


###############################################
# My format pval function, return a bquote expression
##################################################################################################
format_pval <- function(pval){
  if (is.character(pval)){
    # pval contains fold change
    pval <- strsplit(pval, ' ')[[1]]
    fc <- pval[2]
    pval <- as.numeric(pval[1])
  }else{
    fc <- ""
  }
  pval <- format.pval(pval, digits = 2)
  if (grepl("<", pval)){
    pval <- gsub("< ?", "", pval)
    pval <- bquote(italic(P) < .(paste(pval, fc)))
  }else{
    pval <- bquote(italic(P) == .(paste(pval, fc)))
  }
  pval
}


##########################################################
# Correlation test, return a formated pvalue string
##########################################################
cor_pval <- function(...){
  test <- cor.test(...)
  pval <- format_pval(test$p.value)
  # expression(italic('P'), pval)
  # main=expression( italic(p~value) == 0.01 )
  return(pval)
}


#############################
# Clean working environment
#############################
cl <- function(){
  rm(list=ls())
}


##############################
# myheatscatter, add pval
###############################
pvalheatscatter <- function(x, y, pval = TRUE, s_cor = TRUE, pval_size=0.7, main='', main_font=4, cex.lab=1.5,...){
  heatscatter(x, y, cex.lab=cex.lab, main='', ...)
  pval <- cor_pval(x, y, exact = FALSE, m = 's')
  mtext(pval, side = 3, cex=pval_size)
  if (s_cor){
    main <- paste0(main, "\n", "Spearman rho=", round(cor(x, y, m='s', use='pairwise.complete.obs'),2))
  }else{
    main <- main
  }
  title(main, font=main_font)
}


########################################################################################
# Add plot label
########################################################################################
label_plot <- function(...){
  par(xpd=TRUE)
  text(...)
  par(xpd=FALSE)
}


## Cut values to groups
cut_lessthan <- function(counts, n){
  counts <- ifelse(counts >= n,
                   paste0('>=', n),
                   counts)
  # set levels
  counts <- factor(counts, levels = c(as.character(c(0 : (n-1))), paste0('>=', n)))
}
