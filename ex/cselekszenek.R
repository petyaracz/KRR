# -- head -- #

setwd('~/Github/KRR')

library(tidyverse)

# -- read -- #

source("script/functions.R")

train = read_tsv('~/Github/cselekszenek/dat/train.tsv')
test = read_tsv('~/Github/cselekszenek/dat/test.tsv')
dist = read_tsv('~/Github/cselekszenek/dat/word_distances.tsv.gz')

# -- setup --#

# limit n training forms
train = train |> 
  filter(lemma_freq >= 10)

# unique test
test = test |> 
  summarise(
    p_cc = mean(p_cc),
    .by = base_ipa
  )

# -- hmm -- #

# train_words = train$lemma
# dist_words = unique(c(dist$word1, dist$word2))
# missing = setdiff(train_words, dist_words)
# length(missing) # no
# 
# nrow(train)
# length(unique(train$lemma)) # no
# 
# n = length(train_words)
# 
# # should be n^2
# dist |> 
#   filter(word1 %in% train_words, word2 %in% train_words) |> 
#   nrow()
# 
# n^2

# test$base_ipa %in% train$base_ipa
train |> filter(duplicated(base_ipa) | duplicated(base_ipa, fromLast = TRUE))
test |> filter(duplicated(base_ipa) | duplicated(base_ipa, fromLast = TRUE))

# -- model -- #

train$p_cc = pmin(pmax(train$p_cc, 0.001), 0.999)

# Train KRR with LOO hyperparameter tuning on corpus data.
# link = "logit" transforms probabilities to log-odds before fitting and
# back-transforms predictions; criterion = "rmse" selects the (sigma, alpha)
# pair with the lowest LOO root-mean-squared error in the transformed space.
m <- train_krr(train, dist, word_col = "base_ipa", outcome_col = "p_cc", link = "logit", criterion = 'rmse')
m$sigma   # best sigma
m$alpha   # best alpha
m$best_score
m$predictions

# Predict on nonwords using the best hyperparameters found above
preds <- predict_krr(train, test, dist, word_col = "base_ipa", outcome_col = "p_cc",
                     sigma = m$sigma, alpha = m$alpha, link = "logit")

# Multiple outcomes
# Example of looping over multiple outcome columns, e.g. different response types:
# results <- map(outcome_cols, ~ {
#   m <- train_krr(c, dist, "lemma", .x, link = "logit")
#   predict_krr(c, d, dist, "lemma", .x, sigma = m$sigma, alpha = m$alpha, link = "logit")
# }) |> list_rbind()

# -- viz -- #

# Scatter plot of observed vs LOO-predicted values on training data with smoother
m$predictions |> 
  ggplot(aes(observed,predicted_loo)) +
  geom_point() +
  geom_smooth()

# Scatter plot of observed vs predicted values on test nonwords with smoother
preds |> 
  ggplot(aes(observed,predicted)) +
  geom_point() +
  geom_smooth()

# Pearson correlation between observed and predicted values on test data
with(preds, cor.test(observed,predicted))
# upszika