##'---
##' title: Cerevisiae Regression Models
##' author: Jun Cheng

library(data.table)
library(LSD)
library(ggplot2)
library(Biostrings)
library(tidyr)
library(glmnet)
source('./src/utilities.R')

codons <- readRDS('./data/codons.rds')
tAIs <- fread('./data/tAI_per_codon_SC_SP.csv')

## A function to refit lm by throw out non-significant covariates
lm_refit <- function(fit, y = 'hlt.wt', except = codons){
  dt <- model.matrix(formula(paste0(y, "~.")), fit$model)
  dt <- cbind(y = fit$model[, y], dt[, -1])
  coefs <- summary(fit)$coeff
  use_variables <- (coefs[-1,4] < 0.05) | (rownames(coefs)[-1] %in% except)
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
  title(main=paste('R-squared=',round(summary(fit2)$r.squared,digit=3)),cex.main=2,cex.lab=2,xlab='Predicted half-life [min]',ylab='Measured half-life [min]')
  abline(0,1,lwd=2)
  return(list(
    measured=exp(y),
    selected_features = selected_features,
    predicted=exp(predicted[order(as.numeric(names(predicted)))])))
}

UTRs <- readRDS('./data/Sun_mutation_UTRs.rds')
## #' ### upstream AUG, make into a indicator variable
UTRs$num_uAUG <- UTRs$num_uAUG > 0

##'  ## Fit codon models
codon_counts <- UTRs[ ,c(codons, "CDS_length", 'hlt.wt'), with=F]
for (cod in codons){
  codon_counts[, (cod) := get(cod) / (CDS_length)]
  #codon_counts[, (cod) := log((get(cod) + 1) / (CDS_length + 1))]
}
codon_counts <- codon_counts[, CDS_length := NULL]
codon_counts[, hlt.wt := log(hlt.wt)]
codon_counts[, TTT := NULL]
fit <- lm(data = codon_counts, hlt.wt ~ .)
# Fit with cross-validation
cv_rsquare2(fit, return_predict = FALSE)

fitted_coef <- as.data.table(coef(fit), keep.rownames = T)
setnames(fitted_coef, c("Codon", "beta_scer"))
fitted_coef <- fitted_coef[Codon %in% codons]
tAIs_coef <- merge(fitted_coef, tAIs, by = 'Codon')
#tAIs_coef_scer <- tAIs_coef[, .(Codon, beta_scer, S_cerevisiae)]
saveRDS(tAIs_coef, "./data/Scer_coef_tAI.rds")

library(cowplot)
ggplot(tAIs_coef, aes(beta_scer, log(S_cerevisiae))) +
  geom_text(aes(label = Codon)) +
  labs(title = "Compare fitted codon coefficients with sTAI",
       x = "Fitted Coefficients",
       y = "sTAI")
cor(tAIs_coef$beta_scer, log(tAIs_coef$S_cerevisiae))

##' ## Model with tAI instead of each codon as independent metrics
model_data <- UTRs[, c("hlt.wt", 'tAI')]
model_data <- model_data[, c("hlt.wt", "tAI") := .(log(hlt.wt), log(tAI))]
fit_tAI <- lm(data = model_data, hlt.wt ~ tAI)
# Variance explained in cross-validation
cv_rsquare2(fit_tAI)
# This gives 0.4


##' ## Model with codon usage and all other variables
# normalize codon counts with CDS_length ---> (#codon + 1) / (CDS_length+1)
for (cod in codons){
  UTRs[, (cod) := get(cod) / CDS_length]
  #UTRs[, (cod) := log((get(cod) + 1) / (CDS_length + 1))]
}

UTRs[, TTT := NULL]
codon_no_TTT <- codons[codons != "TTT"]

model_data = UTRs[,c('hlt.wt',
                     'UTR5_length',
                     'UTR3_length',
                     'CDS_length'
                     ),with=F]

model_data <- cbind(log(model_data),
                    UTRs[ , c(
                      codon_no_TTT,
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

fit_codons_init <- lm(data = model_data, hlt.wt ~ .)
fit_codons_refit <- lm_refit(fit_codons_init, y = "hlt.wt")
write.csv(coef(summary(fit_codons_init)), './fig/Paper/Supplements/Scer_regression_coef_init.csv')
write.csv(coef(summary(fit_codons_refit)), './fig/Paper/Supplements/Scer_regression_coef.csv')

prediction_res <- outofsamplePrediction_feat_select(model_data, fold=10, y_col = 'hlt.wt', xlim=c(3,150), ylim=c(3,150))
prediction_res$genename <- selected_genenames
prediction_res[['selected_features']] <- NULL
write.csv(data.frame(prediction_res), './fig/Paper/Supplements/Scer_regression_predictions.csv')

tmp <- colnames(model_data)
tmp[1] <- 'half.life'
colnames(model_data) <- tmp
saveRDS(model_data,'./data/Cerevisiae_model_data.rds')

##' ## glmnet model selection?
y <- model_data$half.life
X <- model.matrix(half.life ~., model_data)
X <- X[,-1]
cvfit <- cv.glmnet(X,y,nfolds=10,alpha=0)
plot(cvfit)

## Do not really helped
