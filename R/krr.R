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
# Binary link note: for a genuine 0/1 label (as opposed to a proportion),
# link = "binary" reuses the logit transform to fit in log-odds space, but
# check_krr_inputs() requires outcomes to be exactly 0/1 (the reverse of the
# logit check), and predictions gain a predicted_class column thresholded
# at `threshold` (default 0.5).

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
  if (link == "logit" || link == "binary") {
    target <- pmax(pmin(target, 1 - epsilon), epsilon)
    target <- qlogis(target)
  }
  target
}

#' Apply inverse link to predictions.
.apply_inverse_link <- function(predictions, link) {
  if (link == "logit" || link == "binary") {
    predictions <- plogis(predictions)
  }
  predictions
}

# -- public functions --------------------------------------------------------

#' Train KRR with LOO cross-validation for hyperparameter selection.
#'
#' @param data       Data frame containing words and outcome.
#' @param dist_df    Long-format distance data frame with columns
#'   `word1`, `word2`, `phon_dist`. Must contain all pairwise distances
#'   (including self-pairs) for words in `data`.
#' @param word_col   Name of the word/lemma column (string).
#' @param outcome_col Name of the outcome column (string).
#' @param link       `"identity"` (default), `"logit"`, or `"binary"`.
#'   `"binary"` is for a genuine 0/1 label rather than a proportion; it fits
#'   in log-odds space like `"logit"` but expects exact 0/1 outcomes and adds
#'   a thresholded `predicted_class` column to `predictions`.
#' @param criterion  `"rmse"` (default), `"r"`, or `"accuracy"`. Tuning metric
#'   for LOO. RMSE and `"r"` (Pearson correlation) are computed in the
#'   link-transformed space. `"accuracy"` compares LOO predictions
#'   thresholded at `threshold` against the raw outcome, and is only
#'   meaningful for `link = "binary"`.
#' @param epsilon    Clamping bound for logit/binary link. Default 0.001.
#' @param threshold  Cutoff applied to probability predictions to obtain
#'   `predicted_class` and to compute the `"accuracy"` criterion when
#'   `link = "binary"`. Default 0.5.
#' @param sigma_grid Numeric vector of sigma values to try.
#' @param alpha_grid Numeric vector of alpha (regularisation) values to try.
#'
#' @return A list with components:
#'   - `sigma`, `alpha`: best hyperparameters
#'   - `best_score`: best tuning metric value (RMSE, r, or accuracy)
#'   - `tuning`: tibble of full tuning grid with score values
#'   - `predictions`: tibble with word, observed (original scale),
#'     predicted_loo (original scale), and predicted_class (if `link = "binary"`)
#'   - `link`, `criterion`, `epsilon`, `threshold`: stored for reference
#'
#' @export
#' @examples
#' \dontrun{
#' train <- readr::read_tsv(system.file("extdata", "lakok_train.tsv", package = "krr"))
#' dist  <- readr::read_tsv(system.file("extdata", "word_distances.tsv.gz", package = "krr"))
#' m <- train_krr(train, dist, word_col = "lemma", outcome_col = "p", link = "logit")
#' }
train_krr <- function(data, dist_df, word_col, outcome_col,
                      link = c("identity", "logit", "binary"),
                      criterion = c("rmse", "r", "accuracy"),
                      epsilon = 0.001,
                      threshold = 0.5,
                      sigma_grid = c(1, 2, 3, 4, 5, 8, 16, 32, 64),
                      alpha_grid = c(1, 10, 100, 1000)) {

  if (!outcome_col %in% names(data)) stop(sprintf("Column '%s' not found in data.", outcome_col))
  if (!word_col %in% names(data)) stop(sprintf("Column '%s' not found in data.", word_col))

  link <- match.arg(link)
  criterion <- match.arg(criterion)

  if (link == "binary" && criterion != "accuracy") {
    warning("link = \"binary\" with criterion = \"", criterion,
             "\": consider criterion = \"accuracy\" to tune directly on classification performance.")
  }

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
        } else if (criterion == "r") {
          cor(target, loo_pred, method = "pearson")
        } else {
          loo_pred_class <- as.numeric(.apply_inverse_link(loo_pred, link) >= threshold)
          mean(loo_pred_class == target_raw)
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

  predictions <- tibble(
    !!word_col  := words,
    observed    = target_raw,
    predicted_loo = loo_pred_out
  )
  if (link == "binary") {
    predictions <- mutate(predictions, predicted_class = as.numeric(predicted_loo >= threshold))
  }

  list(
    sigma      = best$sigma,
    alpha      = best$alpha,
    best_score = best$score,
    tuning     = tuning,
    predictions = predictions,
    link      = link,
    criterion = criterion,
    epsilon   = epsilon,
    threshold = threshold
  )
}

#' Predict on new data using a trained KRR specification.
#'
#' @param train_data  Training data frame (words + outcome).
#' @param test_data   Test data frame (words, and optionally outcome).
#' @param dist_df     Long-format distance data frame. Must cover all
#'   train–train and test–train word pairs.
#' @param word_col    Name of the word/lemma column (string).
#' @param outcome_col Name of the outcome column (string).
#' @param sigma       RBF kernel bandwidth.
#' @param alpha       Ridge regularisation parameter.
#' @param link        `"identity"` (default), `"logit"`, or `"binary"`.
#' @param epsilon     Clamping bound for logit/binary link. Default 0.001.
#' @param threshold   Cutoff applied to probability predictions to obtain
#'   `predicted_class` when `link = "binary"`. Default 0.5.
#'
#' @return A tibble with `word`, `predicted` (original scale), `observed`
#'   (if `outcome_col` exists in `test_data`), and `predicted_class`
#'   (if `link = "binary"`).
#'
#' @export
#' @examples
#' \dontrun{
#' preds <- predict_krr(train, test, dist,
#'                      word_col = "lemma", outcome_col = "p",
#'                      sigma = m$sigma, alpha = m$alpha, link = "logit")
#' }
predict_krr <- function(train_data, test_data, dist_df,
                        word_col, outcome_col,
                        sigma, alpha,
                        link = c("identity", "logit", "binary"),
                        epsilon = 0.001,
                        threshold = 0.5) {

  if (!outcome_col %in% names(train_data)) stop(sprintf("Column '%s' not found in data.", outcome_col))
  if (!word_col %in% names(train_data)) stop(sprintf("Column '%s' not found in data.", word_col))

  link <- match.arg(link)

  train_words <- train_data[[word_col]]
  test_words  <- test_data[[word_col]]

  train_target <- .apply_link(train_data[[outcome_col]], link, epsilon)

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

  if (outcome_col %in% names(test_data)) {
    out <- out |> mutate(observed = test_data[[outcome_col]])
  }

  if (link == "binary") {
    out <- out |> mutate(predicted_class = as.numeric(predicted >= threshold))
  }

  out
}

#' Validate inputs before fitting KRR.
#'
#' Checks for duplicate words, missing distance pairs, incomplete pairwise
#' coverage, and outcome values incompatible with the chosen link. Prints a
#' summary of any issues found.
#'
#' @param train      Training data frame.
#' @param test       Test data frame, or `NULL`.
#' @param dist       Long-format distance data frame (`word1`, `word2`, `phon_dist`).
#' @param word_col   Name of the word column (string).
#' @param outcome_col Name of the outcome column (string).
#' @param link       `"identity"`, `"logit"`, or `"binary"`.
#'
#' @return `TRUE` invisibly if all checks pass, `FALSE` otherwise.
#'
#' @export
check_krr_inputs <- function(train, test = NULL, dist,
                             word_col = "transcription2",
                             outcome_col = "prop_a",
                             link = "identity") {

  ok <- TRUE

  train_words <- train[[word_col]]
  dist_words  <- unique(c(dist$word1, dist$word2))

  dupes_train <- train_words[duplicated(train_words)]
  if (length(dupes_train) > 0) {
    ok <- FALSE
    cat(sprintf("FAIL: %d duplicate values in train$%s\n", length(dupes_train), word_col))
    cat("  ", head(dupes_train, 5), "\n")
  }

  if (!is.null(test)) {
    test_words  <- test[[word_col]]
    dupes_test  <- test_words[duplicated(test_words)]
    if (length(dupes_test) > 0) {
      ok <- FALSE
      cat(sprintf("FAIL: %d duplicate values in test$%s\n", length(dupes_test), word_col))
      cat("  ", head(dupes_test, 5), "\n")
    }
  }

  missing_train <- setdiff(train_words, dist_words)
  if (length(missing_train) > 0) {
    ok <- FALSE
    cat(sprintf("FAIL: %d train words missing from dist\n", length(missing_train)))
    cat("  ", head(missing_train, 5), "\n")
  }

  if (!is.null(test)) {
    missing_test <- setdiff(test_words, dist_words)
    if (length(missing_test) > 0) {
      ok <- FALSE
      cat(sprintf("FAIL: %d test words missing from dist\n", length(missing_test)))
      cat("  ", head(missing_test, 5), "\n")
    }
  }

  n_train       <- length(train_words)
  n_pairs_train <- dist |>
    filter(word1 %in% train_words, word2 %in% train_words) |>
    nrow()
  if (n_pairs_train < n_train^2) {
    ok <- FALSE
    cat(sprintf("FAIL: train distance pairs incomplete (%d of %d expected)\n",
                n_pairs_train, n_train^2))
    n_self <- dist |>
      filter(word1 %in% train_words, word1 == word2) |>
      nrow()
    if (n_self < n_train) {
      cat(sprintf("  -> %d of %d self-pairs (diagonal) missing\n",
                  n_train - n_self, n_train))
    }
  }

  if (!is.null(test)) {
    n_test        <- length(test_words)
    n_pairs_cross <- dist |>
      filter(word1 %in% test_words, word2 %in% train_words) |>
      nrow()
    if (n_pairs_cross < n_test * n_train) {
      ok <- FALSE
      cat(sprintf("FAIL: cross distance pairs incomplete (%d of %d expected)\n",
                  n_pairs_cross, n_test * n_train))
    }
  }

  y      <- train[[outcome_col]]
  if (link == "logit") {
    n_zero <- sum(y == 0, na.rm = TRUE)
    n_one  <- sum(y == 1, na.rm = TRUE)
    if (n_zero + n_one > 0) {
      ok <- FALSE
      cat(sprintf("WARNING: logit link but outcome has %d zeros and %d ones (-> clamped by epsilon during fitting)\n",
                  n_zero, n_one))
    }
  }

  if (link == "binary") {
    n_other <- sum(!y %in% c(0, 1), na.rm = TRUE)
    if (n_other > 0) {
      ok <- FALSE
      cat(sprintf("FAIL: binary link but %d outcome values are not exactly 0 or 1\n", n_other))
    }
  }

  n_na <- sum(is.na(y))
  if (n_na > 0) {
    ok <- FALSE
    cat(sprintf("FAIL: %d NAs in train$%s\n", n_na, outcome_col))
  }

  if (ok) cat("All checks passed.\n")

  invisible(ok)
}
