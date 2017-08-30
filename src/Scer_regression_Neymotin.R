library(data.table)
library(LSD)
library(ggplot2)
library(Biostrings)
source('./src/utilities.R')

## A function to refit lm by throw out non-significant covariates
lm_refit <- function(fit, y = 'hlt.wt'){
  dt <- model.matrix(formula(paste0(y, "~.")), fit$model)
  dt <- cbind(y = fit$model[, y], dt[, -1])
  use_variables <- summary(fit)$coeff[-1,4] < 0.05
  use_variables <- names(use_variables)[use_variables]
  dt <- data.frame(dt[, c('y', use_variables)])
  refit <- lm(y ~ ., data = dt)
  refit
}

outofsamplePrediction_feat_select <- function(selected, y_col='hlt.wt', fold=10, ...){
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
  title(main=paste('R-squared=',round(summary(fit2)$r.squared,digit=2)),cex.main=2,cex.lab=2,xlab='Predicted half-life [min]',ylab='Measured half-life [min]')
  abline(0,1,lwd=2)
  return(list(
    measured=exp(y),
    selected_features = selected_features,
    predicted=exp(predicted[order(as.numeric(names(predicted)))])))
}

UTRs <- readRDS('./data/Sun_mutation_UTRs.rds')
Neymotin <- fread('./data/Neymotin_TableS5.csv')
setnames(Neymotin, c("Syst"), "genename")
UTRs <- merge(UTRs, Neymotin, by = 'genename')
codons <- readRDS('./data/codons.rds')
UTRs[, hlt.wt := thalf]
# normalize codon counts with CDS_length ---> #codon / CDS_length
for (cod in codons){
  #UTRs[, (cod) := (get(cod) + 1) / (CDS_length + 1)]
  UTRs[, (cod) := get(cod) / CDS_length]
}
UTRs[, TTT := NULL]
codon_no_TTT <- codons[codons != "TTT"]

### upstream AUG, make into a indicator variable
UTRs$num_uAUG <- UTRs$num_uAUG > 0


##' ## Model with tAI instead of each codon as independent metrics
model_data <- UTRs[, c("hlt.wt", 'tAI')]
model_data <- model_data[, c("hlt.wt", "tAI") := .(log(hlt.wt), log(tAI))]
fit_tAI <- lm(data = model_data, hlt.wt ~ tAI)
# Variance explained in cross-validation
cv_rsquare2(fit_tAI)
# This gives 0.09


##' ## Model with codon usage and all other variables
##'
model_data = UTRs[,c('hlt.wt',
                     'UTR5_length',
                     'UTR3_length',
                     'CDS_length'
                     ),with=F]

model_data <- cbind(log(model_data),
                    UTRs[ , c(codon_no_TTT,
                              "GC_content_UTR5",
                              "GC_content_UTR3",
                              "GC_content_CDS",
                              'UTR5_deltaG',
                              'UTR3_deltaG',
                              'CDS_deltaG'
                              ),
                         with=F])

model_data <- cbind(model_data,
                    UTRs[,
                         c("stopPlus1",
                           "num_uAUG",
                           'startMinus3',
                           "TGTAAATA",
                           "TGCAT",
                           "ATATTC",
                           "TTTTTTA",
                           "AAACAAA"
                           ),
                         with=F])

selected_genenames <- UTRs$genename[complete.cases(model_data)]
model_data = model_data[complete.cases(model_data),]

## Fit Neymotin data with the same covariates (same model) as with Sun 2013 data
fit_codons_init <- lm(data = model_data, hlt.wt ~ ., na.action=na.exclude)
write.csv(coef(summary(fit_codons_init)), './fig/Paper/Supplements/Scer_regression_Neymotin_coef.csv')
pdf('./fig/Paper/Supplements/Neymotin_prediction.pdf', width = 5, height = 5)
prediction_res <- outofsamplePrediction_feat_select(model_data, fold=10, y_col = 'hlt.wt')
dev.off()
prediction_res$selected_features %>% unlist %>% table %>% sort
prediction_res$genename <- selected_genenames
prediction_res[['selected_features']] <- NULL
write.csv(data.frame(prediction_res), './fig/Paper/Supplements/Scer_regression_predictions_Neymotin.csv')


## Fit Neymotin data without codons
model_data[, (codons) := NULL]
fit_no_codons <- lm(data = model_data, hlt.wt ~ ., na.action=na.exclude)

summary(fit_codons_init)$r.squared - summary(fit_no_codons)$r.squared
