###################################################################################################
# Figure 6 predictions
###################################################################################################
# Load data and functions
source('./src/Paper/ini.R')

##' ## Prediction of mRNA half-life from sequencing features, quantify the contribution of every single features. 
outofsamplePrediction <- function(selected, y_col='hlt.wt', fold=10, species='S.cerevisiae', ...){
  fold = fold
  n = nrow(selected)
  #set.seed(1)
  myFolds  = split(sample(n),rep(1:fold,length=n))
  predicted = c()
  selected_features <- list()
  for (ii in myFolds) {
    testset <- selected[ii,]
    #testset <- data.frame(model.matrix(formula(paste0(y_col, "~.")), testset)[, -1]) # uncomment this if use lm_refit
    trainset <- selected[-ii,]
    fit = lm(as.matrix(paste0(y_col, "~.")), data=trainset)
    ## select dataset
    # fit <- lm_refit(fit, y = y_col) # lm_refit will just throw one iterate non-significant ones and refit
    fit <- step(fit, trace = 0, steps = 10, k = log(n))
    selected_features <- c(selected_features, list(names(coef(fit)[-1])))
    y_hat = predict(fit, testset)
    names(y_hat) = ii
    predicted = c(predicted,y_hat)
  }
  y=selected[,get(y_col)]
  names(y) = c(1:n)
  fit2 = lm(predicted ~ y[names(predicted)])
  heatscatter(exp(predicted), exp(y[names(predicted)]), log='xy', main='',xlab='',ylab='', ...)
  rsq <- round(summary(fit2)$r.squared,digit=2)
  title(main=paste0(species, "\n", paste('R-squared=',rsq)),
        cex.main=1.5,cex.lab=1.5,
        xlab='Predicted half-life [min]',
        ylab='Measured half-life [min]')
  abline(0,1,lwd=2)
  return(list(
    measured=exp(y),
    selected_features = selected_features, 
    predicted=exp(predicted[order(as.numeric(names(predicted)))])))
}

# Load feature matrix
Cerevisiae_model_data <- readRDS('./data/Cerevisiae_model_data.rds') # generated with regression script
Pombe_model_data <- readRDS('./data/Pombe_model_data.rds')

codons_no_TTT <- codons[codons != "TTT"]
feature_groups <- list(codons_no_TTT,
                       c('TGTAAATA','TGCAT','ATATTC','TTTTTTA'),
                       'num_uAUG',
                       'UTR5_length',
                       'CDS_deltaG',
                       'GC_content_UTR5',
                       'UTR3_length')
# Only significant ones
Cerevisiae_model_data <- Cerevisiae_model_data[, c("half.life", unlist(feature_groups)), with = F]

fit_cerevisiae <- lm(data=Cerevisiae_model_data, half.life ~ .)

variance_explained_analysis <- function(model_data, feature_groups, full_fit, N=100){
  # N <- 100  # Repeat 100 times, take the mean
  Resquare_combine_matrix <- matrix(NA, nrow=N, ncol=length(feature_groups))
  Resquare_alone_matrix <- matrix(NA, nrow=N, ncol=length(feature_groups))
  Resquare_add_matrix <- matrix(NA, nrow=N, ncol=length(feature_groups))
  Resquare_drop_matrix <- matrix(NA, nrow=N, ncol=length(feature_groups))
  i = 1
  while(i < (N+1) ){
    Resquare_combine <- c()
    Resquare_alone <- c()
    Resquare_drop <- c() # Droped rsquare if one feature(group) left out
    model_feature <- c('half.life')
    for (feature in feature_groups){
      model_feature <- c(model_feature, feature)
      data_tofit <- model_data[ ,model_feature, with=F]
      fit_feature_combine <- lm(data = data_tofit, half.life ~ .)
      Resquare_combine <- c(Resquare_combine,cv_rsquare2(fit_feature_combine)) # c(summary(fit_feature)$r.squared)
      # Feature alone
      data_tofit <- model_data[,c('half.life',feature),with=F]
      fit_feature_alone <- lm(data = data_tofit, half.life ~ .)
      Resquare_alone <- c(Resquare_alone, cv_rsquare2(fit_feature_alone)) # c(summary(fit_feature)$r.squared)
      # Feature add to the combined model
      droped_fit <- update(full_fit, paste(c('. ~ .', feature), collapse='-'))
      drop <- cv_rsquare2(full_fit) - cv_rsquare2(droped_fit)
      Resquare_drop <- c(Resquare_drop, drop)
    }
    Resquare_add <- c(Resquare_combine[1], diff(Resquare_combine))
    ## Fill into matrix
    Resquare_combine_matrix[i, ] <- Resquare_combine
    Resquare_alone_matrix[i, ] <- Resquare_alone
    Resquare_add_matrix[i, ] <- Resquare_add
    Resquare_drop_matrix[i, ] <- Resquare_drop
    i <- i + 1
    print(i)
  }
  Resquare_combine <- colMeans(Resquare_combine_matrix)
  Resquare_alone <- colMeans(Resquare_alone_matrix)
  Resquare_add <- colMeans(Resquare_add_matrix)
  Resquare_drop <- colMeans(Resquare_drop_matrix)
  data.frame(Individual = round(Resquare_alone, 4) * 100,
             Cumulative = round(Resquare_combine, 4) * 100,
             Add = round(Resquare_add, 4) * 100,
             Drop = round(Resquare_drop, 4) * 100
             )
}

# This step is time consuming, load pre-computed results
#variance_explain <- variance_explained_analysis(Cerevisiae_model_data, feature_groups, fit_cerevisiae, N=100)
variance_explain <- readRDS('./fig/ggplot_obj/variance_explain_cerevisiae.rds')
#saveRDS(variance_explain, './fig/ggplot_obj/variance_explain_cerevisiae.rds')
# 
# feature_groups_pombe <- list('num_uAUG', "UTR5.length", "GC_content_UTR5", "Free_energy_ensemble",
#                           codons, 'CDS_length', 'stopPlus1',
#                           c("CAACCA", "ACCAAC", "TATTTAT", "TTAATGA", "ACTAAT"), 'UTR3.length')
#variance_explain_pombe <- variance_explained_analysis(Pombe_model_data, feature_groups_pombe, N=100)
#saveRDS(variance_explain_pombe, './fig/ggplot_obj/variance_explain_pombe.rds')

library(gridBase)
library(gridExtra)
# Take only Drop
variance_explain <- data.frame(CREs = c('Codon usage',
                                        '3\' UTR motifs',
                                        'uAUG',
                                        '5\' UTR length',
                                        'CDS free energy',
                                        '5\' UTR GC content',
                                        '3\' UTR length'),
                               variance_explain[ ,c('Individual', 'Cumulative', 'Drop')])

theme1 <- ttheme_default(core=list(fg_params=list(hjust=0, x=0.1)))
g2 <- tableGrob(variance_explain, rows=NULL, theme=theme1)
grid.arrange(g2)
save(g2, file = './fig/ggplot_obj/variance_explained.RData')

cairo_pdf('./fig/Paper/Variance_explained.pdf')
grid.arrange(gridExtra::combine(g2,along=1))
dev.off()

# This code that does complete neasted cross-validation takes around 10mins to run.
cairo_pdf('./fig/Paper/Figure5_prediction.pdf', width=9, height=10)
plot.new()
pushViewport(viewport(width=1, height=1,layout=grid.layout(2, 2))) # a 1X3 panel
pushViewport(viewport(layout.pos.col=1,layout.pos.row=1)) # a viewport for one of Cerevisiae prediction
par(fig = gridFIG(), new = TRUE)
grid.text(x=0.05,y=0.95,'A',gp=gpar(cex=2,font=2))
outofsamplePrediction(Cerevisiae_model_data,fold=10,y_col='half.life',ylim=c(2,150),xlim=c(2,150),species='S. cerevisiae')
popViewport()
pushViewport(viewport(layout.pos.col=2, layout.pos.row=1)) 
par(fig = gridFIG(), new = TRUE)
grid.text(x=0,y=0.95,'B',gp=gpar(cex=2,font=2))
outofsamplePrediction(Pombe_model_data,fold=10,y_col='half.life',ylim=c(1,150),xlim=c(1,150),species='S. pombe')
popViewport()
pushViewport(viewport(layout.pos.col=1,layout.pos.row=2)) # a viewport for one of Pombe prediction
grid.text(x=0.5,y=0.8,'Contribution of each sequence features (%)',gp=gpar(cex=1.2,font=2))
load('./fig/ggplot_obj/variance_explained.RData')
grid.arrange(gridExtra::combine(g2,along=1), newpage=F)
grid.text(x=0,y=1,'C',gp=gpar(cex=2,font=2))
popViewport()
pushViewport(viewport(layout.pos.col=2,layout.pos.row=2)) # a viewport for one of Pombe prediction
# single_snv_plot was generated in src/analyze-snv-perturbation.R
single_snv_plot <- readRDS('./fig/ggplot_obj/single_snv_plot.rds')
print(single_snv_plot, newpage = F)
grid.text(x=0.05,y=1,'D',gp=gpar(cex=2,font=2))
popViewport(2)
dev.off()
