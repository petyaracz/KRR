## Kernel ridge regression

- take word distances (word1-word2) based on aligned phonological distance (see [here](https://github.com/petyaracz/JANET))
- take training data with words and continuous behaviour (e.g. English past tense forms and p(regular))
- use hyperparameter tuning with LOO to fit model that predicts word behaviour from behaviour of similar words
- use best model to predict test words

## Repository structure

```
KRR/
├── LICENSE
├── Readme.md
├── dat/          — data files (inputs and precomputed distances)
└── script/       — R scripts for setup, modelling, and analysis
```

## Data

- `lakok_forms.tsv` — unique lemma forms (training + test) used for phonological distance computation via JANET
- `lakok_train.tsv` — training data: real Hungarian verbs from WebCorpus2, with log-odds converted to probability (`p`) for the _lakok/lakom_ alternation
- `lakok_test.tsv` — test data: nonword experimental ratings for the _lakok/lakom_ variation
- `siptar_torkenczy_toth_racz_hungarian.tsv` — Hungarian phonological feature matrix (segments × features such as cons, son, cont, labial, coronal, etc.)
- `word_distances.tsv.gz` — precomputed pairwise aligned phonological distances between word forms (columns: word1, word2, phon_dist)

## Scripts

- `setup.R` — downloads raw data from the Racz2024 repository, applies Hungarian orthography-to-phonological transcription, filters to the _lakok/lakom_ variation, computes outcome probabilities, and writes local TSV files
- `functions.R` — core KRR implementation with two public entry points (`train_krr` for LOO-tuned training, `predict_krr` for prediction on new data), plus internal helpers for building distance matrices, LOO via the PRESS/hat-matrix formula, and logit link transformations
- `main.R` — end-to-end analysis pipeline: reads data, deduplicates training set by highest frequency per base, trains KRR with logit link, predicts on test nonwords, and visualises results

## Usage

```r
source("script/functions.R")

train <- read_tsv("dat/lakok_train.tsv")
test  <- read_tsv("dat/lakok_test.tsv")
dist  <- read_tsv("dat/word_distances.tsv.gz")

# Fit model with LOO hyperparameter tuning
m <- train_krr(train, dist, word_col = "lemma", outcome_col = "p", link = "logit")

# Predict on new (test) words
preds <- predict_krr(train, test, dist,
                     word_col = "lemma", outcome_col = "p",
                     sigma = m$sigma, alpha = m$alpha, link = "logit")
```

See `script/main.R` for the full analysis pipeline including visualisation.

## Dependencies

- R package: `tidyverse`

## License

MIT — see `LICENSE`.

## Related

- [JANET](https://github.com/petyaracz/JANET) — aligned phonological distance computation
- [Racz2024](https://github.com/petyaracz/Racz2024) — source data for the Hungarian _lakok/lakom_ variation