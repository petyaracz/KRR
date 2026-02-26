# -- head -- #

setwd('~/Github/KRR')

library(tidyverse)

# -- read -- #

source("script/functions.R")

train = read_tsv('dat/lakok_train.tsv')
test = read_tsv('dat/lakok_test.tsv')
dist = read_tsv('dat/word_distances.tsv.gz')

# -- setup --#

# Deduplicate: keep only the highest-frequency form per base lemma so that
# each base contributes exactly one training observation
train = train |> 
  group_by(base) |> 
  arrange(-lemma_freq_corrected) |> 
  slice(1) |> 
  ungroup()

# -- model -- #

# Train KRR with LOO hyperparameter tuning on corpus data.
# link = "logit" transforms probabilities to log-odds before fitting and
# back-transforms predictions; criterion = "rmse" selects the (sigma, alpha)
# pair with the lowest LOO root-mean-squared error in the transformed space.
m <- train_krr(train, dist, word_col = "lemma", outcome_col = "p", link = "logit", criterion = 'rmse')
m$sigma   # best sigma
m$alpha   # best alpha
m$best_score
m$predictions

# Predict on nonwords using the best hyperparameters found above
preds <- predict_krr(train, test, dist, word_col = "lemma", outcome_col = "p",
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
