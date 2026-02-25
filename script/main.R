# -- head -- #

setwd('~/Github/KRR')

library(tidyverse)

# -- read -- #

source("script/functions.R")

train = read_tsv('dat/lakok_train.tsv')
test = read_tsv('dat/lakok_test.tsv')
dist = read_tsv('dat/word_distances.tsv.gz')

# -- setup --#

train = train |> 
  group_by(base) |> 
  arrange(-lemma_freq_corrected) |> 
  slice(1) |> 
  ungroup()

# -- model -- #

# Train (LOO tuning on corpus data)
m <- train_krr(train, dist, word_col = "lemma", outcome_col = "p", link = "logit", criterion = 'rmse')
m$sigma   # best sigma
m$alpha   # best alpha
m$best_score
m$predictions

# Predict on nonwords
preds <- predict_krr(train, test, dist, word_col = "lemma", outcome_col = "p",
                     sigma = m$sigma, alpha = m$alpha, link = "logit")

# Multiple outcomes
# results <- map(outcome_cols, ~ {
#   m <- train_krr(c, dist, "lemma", .x, link = "logit")
#   predict_krr(c, d, dist, "lemma", .x, sigma = m$sigma, alpha = m$alpha, link = "logit")
# }) |> list_rbind()

# -- viz -- #

m$predictions |> 
  ggplot(aes(observed,predicted_loo)) +
  geom_point() +
  geom_smooth()

preds |> 
  ggplot(aes(observed,predicted)) +
  geom_point() +
  geom_smooth()

with(preds, cor.test(observed,predicted))
