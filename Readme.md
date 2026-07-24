## Kernel ridge regression

- take word distances (word1-word2) based on aligned phonological distance (see [here](https://github.com/petyaracz/JANET))
- take training data with words and continuous behaviour (e.g. English past tense forms and p(regular))
- use hyperparameter tuning with LOO to fit model that predicts word behaviour from behaviour of similar words
- use best model to predict test words

## Installation

```r
# install.packages("pak")
pak::pak("petyaracz/KRR")
```

## Repository structure

```
KRR/
├── DESCRIPTION
├── NAMESPACE
├── LICENSE
├── Readme.md
├── R/                — package source (krr.R, krr-package.R)
├── inst/
│   ├── extdata/      — bundled data files (lakok_train, lakok_test, word_distances)
│   └── scripts/      — setup script for regenerating data from source
├── tests/testthat/   — unit tests
├── dat/              — raw data files (inputs and precomputed distances)
├── script/           — standalone R scripts (main.R, functions.R)
└── ex/               — examples from other projects
```

## Data

The package ships bundled data accessible via `system.file("extdata", ..., package = "krr")`:

- `lakok_train.tsv` — training data: real Hungarian verbs from WebCorpus2, with log-odds converted to probability (`p`) for the _lakok/lakom_ alternation
- `lakok_test.tsv` — test data: nonword experimental ratings for the _lakok/lakom_ variation
- `word_distances.tsv.gz` — precomputed pairwise aligned phonological distances between word forms (columns: word1, word2, phon_dist)

Additional raw files in `dat/`:

- `lakok_forms.tsv` — unique lemma forms (training + test) used for phonological distance computation via JANET
- `siptar_torkenczy_toth_racz_hungarian.tsv` — Hungarian phonological feature matrix (segments × features such as cons, son, cont, labial, coronal, etc.)

## Outcome types (`link`)

- `"identity"` (default) — continuous outcome, fit directly (like `lm`)
- `"logit"` — outcome is a proportion in `[0, 1]` (e.g. p(regular) across a word's occurrences); fit in log-odds space, predictions are back-transformed probabilities
- `"binary"` — outcome is a genuine 0/1 label for a single item (e.g. "is this a proper noun"), not a proportion; fits in log-odds space like `"logit"`, but `check_krr_inputs()` requires outcomes to be exactly 0/1, and `train_krr()`/`predict_krr()` add a `predicted_class` column thresholded at `threshold` (default 0.5). Hyperparameter tuning defaults to `criterion = "rmse"` on the log-odds scale; pass `criterion = "accuracy"` to tune directly on classification accuracy instead (a warning is printed otherwise).

## Functions

- `train_krr()` — LOO-tuned training: grid search over sigma (RBF bandwidth) and alpha (ridge regularisation), returns best hyperparameters, full tuning grid, and LOO predictions
- `predict_krr()` — predict on new words using supplied sigma and alpha
- `check_krr_inputs()` — validate inputs before fitting: checks for duplicates, missing distance pairs, incomplete pairwise coverage, and outcome values incompatible with the chosen link

## Usage

```r
library(krr)
library(readr)

train <- read_tsv(system.file("extdata", "lakok_train.tsv", package = "krr"))
test  <- read_tsv(system.file("extdata", "lakok_test.tsv",  package = "krr"))
dist  <- read_tsv(system.file("extdata", "word_distances.tsv.gz", package = "krr"))

# Optional: validate inputs
check_krr_inputs(train, test, dist, word_col = "lemma", outcome_col = "p", link = "logit")

# Fit model with LOO hyperparameter tuning
m <- train_krr(train, dist, word_col = "lemma", outcome_col = "p", link = "logit")

# Predict on new (test) words
preds <- predict_krr(train, test, dist,
                     word_col = "lemma", outcome_col = "p",
                     sigma = m$sigma, alpha = m$alpha, link = "logit")
```

See `inst/scripts/setup.R` for data preparation and `script/main.R` for the full analysis pipeline including visualisation.

## Dependencies

Imports: `dplyr`, `tidyr`, `purrr`, `tibble`

Suggests: `testthat`, `readr`, `ggplot2`, `knitr`, `rmarkdown`

## License

MIT — see `LICENSE`.

## Related

- [JANET](https://github.com/petyaracz/JANET) — aligned phonological distance computation
- [Racz2024](https://github.com/petyaracz/Racz2024) — source data for the Hungarian _lakok/lakom_ variation
- ex: various examples from other projects, with slightly adapted scripts, for reference.