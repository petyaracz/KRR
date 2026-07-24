# Manual smoke test for link = "binary" in train_krr()/predict_krr()/check_krr_inputs().
# Run this after reinstalling the package (devtools::install() or pak::pak("petyaracz/KRR", upgrade = TRUE)).

library(krr)      # the package itself
library(dplyr)    # for mutate() below

set.seed(1337)    # reproducible if any randomness creeps in later

# -- toy data: synthetic 4-letter "words" over {a, b} --------------------------
# label = 1 if a word has 3+ b's, 0 if it has <=1 b -- words with exactly 2 b's
# are left out so the two classes have a clean margin (min Hamming distance 2)
# instead of sitting right on the decision boundary. Swap this out for real
# phonological distances (e.g. via JANET) and a real binary outcome column
# for actual use.

train = tibble::tibble(
  word = c("aaaa", "aaab", "aaba", "abaa",   # <=1 b
           "bbbb", "abbb", "babb", "bbab"),  # >=3 b
  y    = c(0,      0,      0,      0,
           1,      1,      1,      1)
)

test = tibble::tibble(
  word = c("baaa", "bbba"),   # one held-out word per class, same margin
  y    = c(0,      1)
)

# -- build a distance table from character (Hamming) distance -----------------

all_words = union(train$word, test$word)                     # every word we need distances for
dist_mat  = adist(all_words)                                 # equal-length words -> adist == Hamming distance
dimnames(dist_mat) = list(all_words, all_words)               # name rows/cols for the pivot below

dist_df = as.data.frame(dist_mat) |>                          # matrix -> data frame
  tibble::rownames_to_column("word1") |>                      # row names -> word1 column
  tidyr::pivot_longer(-word1, names_to = "word2", values_to = "phon_dist")  # wide -> long

# -- 1. validate inputs before fitting -----------------------------------------

check_krr_inputs(train, test, dist_df,
                 word_col = "word", outcome_col = "y", link = "binary")

# sanity check: validation should FAIL if the outcome isn't exactly 0/1
bad_train = train |> mutate(y = c(0.1, 0, 0.2, 0, 1, 0.9, 1, 0.8))  # proportions, not labels
check_krr_inputs(bad_train, test, dist_df,
                 word_col = "word", outcome_col = "y", link = "binary")

# -- 2. train with default criterion (rmse) -> expect a warning ----------------

m_rmse = train_krr(train, dist_df, word_col = "word", outcome_col = "y",
                   link = "binary")   # criterion defaults to "rmse"; should warn to use "accuracy"

print(m_rmse$predictions)             # predicted_loo (probability) and predicted_class (0/1)

# -- 3. train tuned directly on classification accuracy ------------------------

m_acc = train_krr(train, dist_df, word_col = "word", outcome_col = "y",
                  link = "binary", criterion = "accuracy")   # no warning this time

print(m_acc$best_score)               # best LOO accuracy found during tuning
print(m_acc$predictions)

# -- 4. predict on held-out test words ------------------------------------------

preds = predict_krr(train, test, dist_df,
                    word_col = "word", outcome_col = "y",
                    sigma = m_acc$sigma, alpha = m_acc$alpha,
                    link = "binary")   # threshold defaults to 0.5

print(preds)                          # predicted (probability), observed, predicted_class

# -- 5. try a stricter threshold -------------------------------------------------

preds_strict = predict_krr(train, test, dist_df,
                           word_col = "word", outcome_col = "y",
                           sigma = m_acc$sigma, alpha = m_acc$alpha,
                           link = "binary", threshold = 0.7)  # only call it "1" if p >= 0.7

print(preds_strict)
