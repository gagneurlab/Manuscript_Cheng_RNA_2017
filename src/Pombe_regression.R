##'---
##' title: Pombe Regression Models
##' author: Jun Cheng

#' ## Pombe regression model
library(rmarkdown)
library(data.table)
library(Biostrings)
library(LSD)
library(glmnet)
source('./src/utilities.R')
tAIs <- fread('./data/tAI_per_codon_SC_SP.csv')

TUs <- readRDS('./data/TUs_Pombe.rds')
TUs <- TUs[!is.na(UTR5_seq) & !is.na(UTR3_seq)]
TUs <- TUs[UTR5.length > 6]
TUs <- subset(TUs, half.life > 0.5 & half.life < 250)
TUs <- TUs[half.life > 1]
TUs <- TUs[half.life < 250]
TUs$num_uAUG <- TUs$num_uAUG > 0
codons <- readRDS('./data/codons.rds')


##'  ## First with do not log transform
codon_counts <- TUs[ ,c(codons, "CDS_length", 'half.life'), with=F]
for (cod in codons){
  codon_counts[, (cod) := get(cod) / (CDS_length)]
  #codon_counts[, (cod) := log((get(cod) + 1) / (CDS_length + 1))]
}
codon_counts <- codon_counts[, CDS_length := NULL]
codon_counts[, half.life := log(half.life)]
codon_counts[, TTT := NULL]
fit <- lm(data = codon_counts, half.life ~ .)
# Fit with cross-validation
cv_rsquare2(fit, return_predict = FALSE)

fitted_coef <- as.data.table(coef(fit), keep.rownames = T)
setnames(fitted_coef, c("Codon", "beta_spom"))
fitted_coef <- fitted_coef[Codon %in% codons]
tAIs_coef <- merge(fitted_coef, tAIs, by = 'Codon')
#tAIs_coef_scer <- tAIs_coef[, .(Codon, beta_scer, S_cerevisiae)]
saveRDS(tAIs_coef, "./data/Spom_coef_tAI.rds")

library(cowplot)
ggplot(tAIs_coef, aes(beta_spom, log(S_pombe))) +
  geom_text(aes(label = Codon)) +
  labs(title = "Compare fitted codon coefficients with sTAI",
       x = "Fitted Coefficients",
       y = "sTAI")
cor(tAIs_coef$beta_spom, log(tAIs_coef$S_pombe))

##' ## Fit complete model
# normalize codon counts with CDS_length ---> #codon / CDS_length
for (cod in codons){
  TUs[, (cod) := get(cod) / CDS_length]
  #TUs[, (cod) := log((get(cod) + 1) / (CDS_length + 1))]
}

TUs[, TTT := NULL]
codon_no_TTT <- codons[codons != "TTT"]

train_test_split <- function(dt, fold = 5){
  n <- dim(dt)[1]
  myFolds  = split(sample(n),rep(1:fold,length=n))
  # make folder for each fold
  for(i in seq(fold)){
    dir.create(paste0("fold", i))
    # store train, test at each folder
    trainset <- dt[-myFolds[[i]],]
    testset <- dt[myFolds[[i]],]
    saveRDS(trainset, paste0("fold", i, "/train.rds"))
    saveRDS(testset, paste0("fold", i, "/test.rds"))
  }
  return(myFolds)
}

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

##' ## Include other features
model_data <- TUs[, c(
                      'half.life',
                      'UTR3.length',
                      "UTR5.length",
                      "CDS_length"),
                  with=F]

model_data <- cbind(log(model_data),
                TUs[,
                  c(codon_no_TTT,
                    "stopPlus1",
                    "num_uAUG",
                    'PositionM3',
                    "CAACCA",
                    "ACCAAC",
                    "TATTTAT",
                    "TTAATGA",
                    "ACTAAT",
                    "GC_content_UTR5",
                    "GC_content_UTR3",
                    "UTR5_deltaG",
                    "UTR3_deltaG",
                    "CDS_deltaG",
                    "GC_content_CDS"
                    ),
                  with=F])

selected_genenames <- TUs$Txn.unit[complete.cases(model_data)]
model_data <- model_data[complete.cases(model_data)]
saveRDS(model_data,'./data/Pombe_model_data.rds')

fit_codons_init <- lm(data = model_data, half.life ~ .)
fit_codons_refit <- lm_refit(fit_codons_init, y = "half.life")
write.csv(coef(summary(fit_codons_init)), './fig/Paper/Supplements/Pombe_regression_coef_init.csv')
write.csv(coef(summary(fit_codons_refit)), './fig/Paper/Supplements/Pombe_regression_coef.csv')

prediction_res <- outofsamplePrediction_feat_select(model_data, fold=10, y_col='half.life', xlim = c(1, 100), ylim = c(1, 100))
prediction_res$selected_features %>% unlist %>% table %>% sort
prediction_res$genename <- selected_genenames
prediction_res[['selected_features']] <- NULL
write.csv(data.frame(prediction_res), './fig/Paper/Supplements/Spom_regression_predictions.csv')

## Try with regularization
X <- model.matrix(half.life ~., model_data)
X <- X[, -1]
y <- model_data$half.life
cvfit <- cv.glmnet(X, y, nfolds=10, alpha=0)
plot(cvfit)
## Do not really helped
