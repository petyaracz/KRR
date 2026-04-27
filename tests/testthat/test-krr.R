test_that("train_krr returns expected list structure", {
  m <- train_krr(.train, .train_dist, word_col = "word", outcome_col = "p", link = "logit")
  expect_named(m, c("sigma", "alpha", "best_score", "tuning", "predictions",
                    "link", "criterion", "epsilon"))
  expect_gt(m$sigma, 0)
  expect_gt(m$alpha, 0)
})

test_that("train_krr predictions are in [0, 1] with logit link", {
  m <- train_krr(.train, .train_dist, word_col = "word", outcome_col = "p", link = "logit")
  expect_true(all(m$predictions$predicted_loo >= 0 & m$predictions$predicted_loo <= 1))
})

test_that("train_krr works with identity link", {
  m <- train_krr(.train, .train_dist, word_col = "word", outcome_col = "p", link = "identity")
  expect_named(m, c("sigma", "alpha", "best_score", "tuning", "predictions",
                    "link", "criterion", "epsilon"))
})

test_that("train_krr works with criterion = 'r'", {
  m <- train_krr(.train, .train_dist, word_col = "word", outcome_col = "p",
                 link = "logit", criterion = "r")
  expect_gte(m$best_score, -1)
  expect_lte(m$best_score,  1)
})

test_that("train_krr errors on missing columns", {
  expect_error(train_krr(.train, .train_dist, word_col = "word", outcome_col = "missing"))
  expect_error(train_krr(.train, .train_dist, word_col = "missing", outcome_col = "p"))
})

test_that("predict_krr returns tibble with word, predicted, observed columns", {
  m     <- train_krr(.train, .train_dist, word_col = "word", outcome_col = "p", link = "logit")
  preds <- predict_krr(.train, .test, .dist_full,
                       word_col = "word", outcome_col = "p",
                       sigma = m$sigma, alpha = m$alpha, link = "logit")
  expect_true(tibble::is_tibble(preds))
  expect_true(all(c("word", "predicted", "observed") %in% names(preds)))
  expect_equal(nrow(preds), nrow(.test))
})

test_that("predict_krr predictions are in [0, 1] with logit link", {
  m     <- train_krr(.train, .train_dist, word_col = "word", outcome_col = "p", link = "logit")
  preds <- predict_krr(.train, .test, .dist_full,
                       word_col = "word", outcome_col = "p",
                       sigma = m$sigma, alpha = m$alpha, link = "logit")
  expect_true(all(preds$predicted >= 0 & preds$predicted <= 1))
})

test_that("predict_krr omits observed when test has no outcome column", {
  test_no_p <- dplyr::select(.test, -p)
  m     <- train_krr(.train, .train_dist, word_col = "word", outcome_col = "p", link = "logit")
  preds <- predict_krr(.train, test_no_p, .dist_full,
                       word_col = "word", outcome_col = "p",
                       sigma = m$sigma, alpha = m$alpha, link = "logit")
  expect_false("observed" %in% names(preds))
  expect_true("predicted" %in% names(preds))
})

test_that("check_krr_inputs returns TRUE on clean data", {
  result <- capture.output(
    ok <- check_krr_inputs(.train, .test, .dist_full,
                           word_col = "word", outcome_col = "p", link = "logit")
  )
  expect_true(ok)
})

test_that("check_krr_inputs returns FALSE with duplicate train words", {
  train_dupes <- dplyr::bind_rows(.train, .train[1, ])
  capture.output(
    ok <- check_krr_inputs(train_dupes, NULL, .train_dist,
                           word_col = "word", outcome_col = "p")
  )
  expect_false(ok)
})

test_that("check_krr_inputs returns FALSE with NAs in outcome", {
  train_na    <- .train
  train_na$p[1] <- NA
  capture.output(
    ok <- check_krr_inputs(train_na, NULL, .train_dist,
                           word_col = "word", outcome_col = "p")
  )
  expect_false(ok)
})

test_that("check_krr_inputs returns FALSE when logit link sees exact 0 or 1", {
  train_boundary    <- .train
  train_boundary$p  <- c(0, 0.3, 0.7, 1)
  capture.output(
    ok <- check_krr_inputs(train_boundary, NULL, .train_dist,
                           word_col = "word", outcome_col = "p", link = "logit")
  )
  expect_false(ok)
})

test_that("logit link and inverse-logit are exact inverses", {
  x           <- c(0.1, 0.3, 0.5, 0.7, 0.9)
  transformed <- krr:::.apply_link(x, "logit", epsilon = 0.001)
  recovered   <- krr:::.apply_inverse_link(transformed, "logit")
  expect_equal(recovered, x, tolerance = 1e-10)
})

test_that("identity link is a no-op", {
  x           <- c(0.1, 0.5, 0.9)
  transformed <- krr:::.apply_link(x, "identity", epsilon = 0.001)
  recovered   <- krr:::.apply_inverse_link(transformed, "identity")
  expect_equal(recovered, x)
})
