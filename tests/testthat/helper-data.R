# Synthetic fixtures shared across test files.
# Four 2-character train words, two test words.
# Distances are Hamming distances between characters.

.train <- tibble::tibble(
  word = c("aa", "ab", "ba", "bb"),
  p    = c(0.1,  0.3,  0.7,  0.9)
)

# Complete train x train distance table (symmetric, zero diagonal)
.train_dist <- tibble::tibble(
  word1     = c("aa","aa","aa","aa", "ab","ab","ab","ab", "ba","ba","ba","ba", "bb","bb","bb","bb"),
  word2     = c("aa","ab","ba","bb", "aa","ab","ba","bb", "aa","ab","ba","bb", "aa","ab","ba","bb"),
  phon_dist = c(  0,   1,   1,   2,    1,   0,   2,   1,    1,   2,   0,   1,    2,   1,   1,   0)
)

.test <- tibble::tibble(
  word = c("ac", "bc"),
  p    = c(0.2,  0.8)
)

# test x train cross-distances (word1 = test, word2 = train)
.cross_dist <- tibble::tibble(
  word1     = c("ac","ac","ac","ac", "bc","bc","bc","bc"),
  word2     = c("aa","ab","ba","bb", "aa","ab","ba","bb"),
  phon_dist = c(  1,   1,   2,   2,    2,   2,   1,   1)
)

.dist_full <- dplyr::bind_rows(.train_dist, .cross_dist)

# Binary-label fixtures (same words/distances, 0/1 outcome instead of a proportion)
.train_bin <- tibble::tibble(
  word = c("aa", "ab", "ba", "bb"),
  y    = c(0,    0,    1,    1)
)

.test_bin <- tibble::tibble(
  word = c("ac", "bc"),
  y    = c(0,    1)
)
