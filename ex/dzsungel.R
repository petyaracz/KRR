# -- head -- #

setwd('~/Github/KRR')

library(tidyverse)

# -- read -- #

source("script/functions.R")

train = read_tsv('~/Documents/Misc/habilitation/dzsungel/dat/dzsungel.tsv')
dist = read_tsv('~/Documents/Misc/habilitation/dzsungel/dat/word_distances.tsv.gz')

# -- setup --#

# unique training forms w/ overall rate of back/front
nrow(distinct(train,stem)) == nrow(train)

# -- check -- #

check_krr_inputs(train = train, dist = dist, word_col = 'transcribed', outcome_col = 'p_back')

# -- model -- #

# Train KRR with LOO hyperparameter tuning on corpus data.
# link = "logit" transforms probabilities to log-odds before fitting and
# back-transforms predictions; criterion = "rmse" selects the (sigma, alpha)
# pair with the lowest LOO root-mean-squared error in the transformed space.
m <- train_krr(train, dist, word_col = "transcribed", outcome_col = "log_odds_back", link = "identity", criterion = 'rmse')
m$sigma   # best sigma
m$alpha   # best alpha
m$best_score # logit: 2.28 # identity: 2.39
m$predictions

# Predict on nonwords using the best hyperparameters found above
# preds <- predict_krr(train, test, dist, word_col = "base_ipa", outcome_col = "p_cc",
                     # sigma = m$sigma, alpha = m$alpha, link = "logit")

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
# preds |> 
#   ggplot(aes(observed,predicted)) +
#   geom_point() +
#   geom_smooth()

# Pearson correlation between observed and predicted values on test data
with(m$predictions, cor.test(observed,predicted_loo))
# upszika

# -- alright -- #

pred = m$predictions |> 
  mutate(
    sigma = m$sigma,   # best sigma
    alpha = m$alpha,   # best alpha
    rmse = m$best_score
  )

write_tsv(pred, '~/Documents/Misc/habilitation/dzsungel/dat/dzsungel_pred.tsv')
