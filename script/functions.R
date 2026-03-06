# krr.R
# Kernel ridge regression using phonological distance matrices.
# Two entry points:
#   train_krr()   — fit on training data with LOO tuning, return best hyperparameters
#   predict_krr() — predict on new data using supplied sigma and alpha
#
# Kernel: RBF / Gaussian kernel  K(x,y) = exp(-d(x,y)^2 / (2σ^2))
#   where d(x,y) is the precomputed aligned phonological distance.
#
# Hyperparameter tuning: grid search over sigma (RBF bandwidth) and alpha
#   (ridge regularisation strength).  The search minimises LOO RMSE or
#   maximises LOO Pearson r, both evaluated in the (optionally link-transformed)
#   target space.
#
# LOO efficiency: leave-one-out predictions are computed via the PRESS /
#   hat-matrix formula without refitting for each fold:
#   loo_i = y_i - (y_i - ŷ_i) / (1 - H_ii)
#
# Logit link note: unlike GLM binomial, KRR is not likelihood-based.
# With link = "logit", targets are logit-transformed before fitting and
# predictions are inverse-logit-transformed back. Exact 0/1 targets are
# clamped to [epsilon, 1 - epsilon] before transformation.
#
# Dependency: tidyverse

library(tidyverse)

# -- internal helpers --------------------------------------------------------

#' Build a distance matrix from long-format distance data.
#' Returns a matrix with rows ordered by row_words, columns by col_words.
.build_dist_matrix <- function(dist_df, row_words, col_words) {
  mat <- dist_df |>
    filter(word1 %in% row_words, word2 %in% col_words) |>
    pivot_wider(names_from = word2, values_from = phon_dist) |>
    arrange(factor(word1, levels = row_words)) |>
    select(all_of(col_words)) |>
    as.matrix()
  
  n_na <- sum(is.na(mat))
  if (n_na > 0) {
    warning(sprintf("%d missing distances set to 0.", n_na))
    mat[is.na(mat)] <- 0
  }
  
  mat
}

#' LOO predictions via the hat-matrix shortcut (PRESS formula).
#' Takes a raw distance matrix (not a kernel), computes kernel internally.
.loo_krr <- function(sigma, alpha, dist_matrix, target) {
  K <- exp(-dist_matrix^2 / (2 * sigma^2))
  n <- nrow(K)
  K_inv <- solve(K + alpha * diag(n))
  H <- K %*% K_inv
  fitted <- H %*% target
  residuals <- target - fitted
  loo_pred <- target - residuals / (1 - diag(H))
  as.vector(loo_pred)
}

#' Apply link transformation to target vector.
.apply_link <- function(target, link, epsilon) {
  if (link == "logit") {
    target <- pmax(pmin(target, 1 - epsilon), epsilon)
    target <- qlogis(target)
  }
  target
}

#' Apply inverse link to predictions.
.apply_inverse_link <- function(predictions, link) {
  if (link == "logit") {
    predictions <- plogis(predictions)
  }
  predictions
}

# -- public functions --------------------------------------------------------

#' Train KRR with LOO cross-validation for hyperparameter selection.
#'
#' @param data       Data frame containing words and outcome.
#' @param dist_df    Long-format distance data frame with columns
#'                   word1, word2, phon_dist. Must contain all pairwise
#'                   distances (including self-pairs) for words in data.
#' @param word_col   Name of the word/lemma column (string).
#' @param outcome_col Name of the outcome column (string).
#' @param link       "identity" (default) or "logit".
#' @param criterion  "rmse" (default) or "r". Tuning metric for LOO.
#'                   RMSE is computed in the link-transformed space.
#'                   "r" uses Pearson correlation in link-transformed space.
#' @param epsilon    Clamping bound for logit link. Default 0.001.
#' @param sigma_grid Numeric vector of sigma values to try.
#' @param alpha_grid Numeric vector of alpha (regularisation) values to try.
#'
#' @return A list with components:
#'   - sigma, alpha: best hyperparameters
#'   - best_score: best tuning metric value (RMSE or r, in transformed space)
#'   - tuning: tibble of full tuning grid with score values
#'   - predictions: tibble with word, observed (original scale),
#'                  and predicted_loo (original scale)
#'   - link, criterion, epsilon: stored for reference
train_krr <- function(data, dist_df, word_col, outcome_col,
                      link = c("identity", "logit"),
                      criterion = c("rmse", "r"),
                      epsilon = 0.001,
                      sigma_grid = c(1, 2, 3, 4, 5, 8, 16, 32, 64),
                      alpha_grid = c(1, 10, 100, 1000)) {
  
  link <- match.arg(link)
  criterion <- match.arg(criterion)
  
  words <- data[[word_col]]
  target_raw <- data[[outcome_col]]
  target <- .apply_link(target_raw, link, epsilon)
  
  dist_matrix <- .build_dist_matrix(dist_df, words, words)
  
  tuning <- crossing(sigma = sigma_grid, alpha = alpha_grid) |>
    mutate(
      score = map2_dbl(sigma, alpha, ~ {
        loo_pred <- .loo_krr(.x, .y, dist_matrix, target)
        if (criterion == "rmse") {
          sqrt(mean((target - loo_pred)^2))
        } else {
          cor(target, loo_pred, method = "pearson")
        }
      })
    )
  
  best <- if (criterion == "rmse") {
    filter(tuning, score == min(score)) |> slice(1)
  } else {
    filter(tuning, score == max(score)) |> slice(1)
  }
  
  loo_pred <- .loo_krr(best$sigma, best$alpha, dist_matrix, target)
  loo_pred_out <- .apply_inverse_link(loo_pred, link)
  
  list(
    sigma      = best$sigma,
    alpha      = best$alpha,
    best_score = best$score,
    tuning     = tuning,
    predictions = tibble(
      !!word_col  := words,
      observed    = target_raw,
      predicted_loo = loo_pred_out
    ),
    link      = link,
    criterion = criterion,
    epsilon   = epsilon
  )
}

#' Predict on new data using a trained KRR specification.
#'
#' @param train_data  Training data frame (words + outcome).
#' @param test_data   Test data frame (words, and optionally outcome).
#' @param dist_df     Long-format distance data frame. Must cover all
#'                    train-train and test-train word pairs.
#' @param word_col    Name of the word/lemma column (string).
#' @param outcome_col Name of the outcome column (string).
#' @param sigma       RBF kernel bandwidth.
#' @param alpha       Ridge regularisation parameter.
#' @param link        "identity" (default) or "logit".
#' @param epsilon     Clamping bound for logit link. Default 0.001.
#'
#' @return A tibble with word, predicted (original scale), and observed
#'         (if outcome_col exists in test_data).
predict_krr <- function(train_data, test_data, dist_df,
                        word_col, outcome_col,
                        sigma, alpha,
                        link = c("identity", "logit"),
                        epsilon = 0.001) {
  
  link <- match.arg(link)
  
  train_words <- train_data[[word_col]]
  test_words  <- test_data[[word_col]]
  
  train_target <- .apply_link(train_data[[outcome_col]], link, epsilon)
  
  # Build distance matrices, then kernelise
  train_dist <- .build_dist_matrix(dist_df, train_words, train_words)
  test_dist  <- .build_dist_matrix(dist_df, test_words, train_words)
  
  train_kernel <- exp(-train_dist^2 / (2 * sigma^2))
  test_kernel  <- exp(-test_dist^2 / (2 * sigma^2))
  
  n <- nrow(train_kernel)
  coefficients <- solve(train_kernel + alpha * diag(n), train_target)
  predictions <- as.vector(test_kernel %*% coefficients)
  
  predictions <- .apply_inverse_link(predictions, link)
  
  out <- tibble(
    !!word_col := test_words,
    predicted  = predictions
  )
  
  # Attach observed values if available
  if (outcome_col %in% names(test_data)) {
    out <- out |> mutate(observed = test_data[[outcome_col]])
  }
  
  out
}

# nice checking function for various failure modes of the data
check_krr_inputs = function(train, test = NULL, dist, 
                            word_col = "transcription2", 
                            outcome_col = "prop_a",
                            link = "identity") {
  
  ok = TRUE
  
  train_words = train[[word_col]]
  dist_words = unique(c(dist$word1, dist$word2))
  
  # 1. Duplicate words in train
  dupes_train = train_words[duplicated(train_words)]
  if (length(dupes_train) > 0) {
    ok = FALSE
    cat(sprintf("FAIL: %d duplicate values in train$%s\n", length(dupes_train), word_col))
    cat("  ", head(dupes_train, 5), "\n")
  }
  
  # 2. Duplicate words in test
  if (!is.null(test)) {
    test_words = test[[word_col]]
    dupes_test = test_words[duplicated(test_words)]
    if (length(dupes_test) > 0) {
      ok = FALSE
      cat(sprintf("FAIL: %d duplicate values in test$%s\n", length(dupes_test), word_col))
      cat("  ", head(dupes_test, 5), "\n")
    }
  }
  
  # 3. Train words missing from distance table
  missing_train = setdiff(train_words, dist_words)
  if (length(missing_train) > 0) {
    ok = FALSE
    cat(sprintf("FAIL: %d train words missing from dist\n", length(missing_train)))
    cat("  ", head(missing_train, 5), "\n")
  }
  
  # 4. Test words missing from distance table
  if (!is.null(test)) {
    missing_test = setdiff(test_words, dist_words)
    if (length(missing_test) > 0) {
      ok = FALSE
      cat(sprintf("FAIL: %d test words missing from dist\n", length(missing_test)))
      cat("  ", head(missing_test, 5), "\n")
    }
  }
  
  # 5. Incomplete pairwise distances (train)
  n_train = length(train_words)
  n_pairs_train = dist |> 
    filter(word1 %in% train_words, word2 %in% train_words) |> 
    nrow()
  if (n_pairs_train < n_train^2) {
    ok = FALSE
    cat(sprintf("FAIL: train distance pairs incomplete (%d of %d expected)\n", 
                n_pairs_train, n_train^2))
    # check self-pairs specifically
    n_self = dist |> 
      filter(word1 %in% train_words, word1 == word2) |> 
      nrow()
    if (n_self < n_train) {
      cat(sprintf("  -> %d of %d self-pairs (diagonal) missing\n", 
                  n_train - n_self, n_train))
    }
  }
  
  # 6. Incomplete cross-distances (test x train)
  if (!is.null(test)) {
    n_test = length(test_words)
    n_pairs_cross = dist |> 
      filter(word1 %in% test_words, word2 %in% train_words) |> 
      nrow()
    if (n_pairs_cross < n_test * n_train) {
      ok = FALSE
      cat(sprintf("FAIL: cross distance pairs incomplete (%d of %d expected)\n", 
                  n_pairs_cross, n_test * n_train))
    }
  }
  
  # 7. Outcome values incompatible with logit link
  y = train[[outcome_col]]
  if (link == "logit") {
    n_zero = sum(y == 0, na.rm = TRUE)
    n_one = sum(y == 1, na.rm = TRUE)
    if (n_zero + n_one > 0) {
      ok = FALSE
      cat(sprintf("FAIL: logit link but outcome has %d zeros and %d ones (-> Inf)\n", 
                  n_zero, n_one))
    }
  }
  
  # 8. NAs in outcome
  n_na = sum(is.na(y))
  if (n_na > 0) {
    ok = FALSE
    cat(sprintf("FAIL: %d NAs in train$%s\n", n_na, outcome_col))
  }
  
  if (ok) cat("All checks passed.\n")
  
  invisible(ok)
}