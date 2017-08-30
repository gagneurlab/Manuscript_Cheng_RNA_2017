library(knitr)
source('./src/utilities.R')
library(data.table)
library(LSD)
library(gridExtra)
library(ggplot2)
library(Biostrings)
library(grid)
library(gridBase)
library(magrittr)
library(ggpval)

Scer <- readRDS('./data/Sun_mutation_UTRs.rds')
Scer <- Scer[!is.na(UTR5_seq) & !is.na(UTR3_seq)]
Scer <- Scer[UTR5_length > 6]
Scer_syn <- readRDS('./data/sun_mutants_synt.rds') %>% as.data.table
Scer_syn <- merge(Scer_syn, Scer, by = 'genename')
Pombe <- readRDS('./data/TUs_Pombe.rds')
Pombe <- Pombe[!is.na(UTR5_seq) & !is.na(UTR3_seq)]
Pombe <- Pombe[UTR5.length > 6]
Pombe <- subset(Pombe, half.life > 0.5 & half.life < 250)
# Pombe has some extream low half-life transcripts, remove them
Pombe <- Pombe[half.life > 1]
Pombe <- Pombe[half.life < 250]
codons = readRDS('./data/codons.rds')
codon1 <- sapply(codons, function(i) paste0(i, '.1'))
sun_mutants_hlt <- readRDS('./data/sun_mutants_hlt.rds')
UTRs <- readRDS('./data/sun_mutants_hlt_medians.rds')
UTRs <- UTRs[!is.na(UTR5_seq) & !is.na(UTR3_seq)]
UTRs <- UTRs[UTR5_length > 6]
mutants = names(sun_mutants_hlt)
tAIs <- fread('./data/tAI_per_codon_SC_SP.csv')

source('./src/seqLogo.R')
source('./src/pwm.R')

### y axis breaks
y_breaks <- c(1, 5, 10, 20, 40, 150)

# Levels
Scer$startMinus3 <- factor(Scer$startMinus3, levels=c('A','G','T','C'))
stopcodon_level <- c('TAAG','TAAA','TAAT','TAAC','TAGG','TAGA','TAGT','TAGC','TGAG','TGAA','TGAT','TGAC')
Scer$Stop_codon <- sapply(Scer$stopPlus1, function(i) substring(i, 1, 3))
Scer$stopPlus1 <- factor(Scer$stopPlus1,levels = stopcodon_level)

Pombe$PositionM3 <- factor(Pombe$PositionM3, levels=c('A','G','T','C'))
stopcodon_level <- c('TAAG','TAAA','TAAT','TAAC','TAGG','TAGA','TAGT','TAGC','TGAG','TGAA','TGAT','TGAC')
Pombe$stopPlus1 <- factor(Pombe$stopPlus1, levels=stopcodon_level)
Pombe$Stop_codon <- sapply(Pombe$stopPlus1, function(i) substring(i, 1, 3))

# write my session info
# writeLines(capture.output(sessionInfo()), "src/sessionInfo.txt")
