# -- head -- #

setwd('~/Github/KRR')

library(tidyverse)

# -- read -- #

# Download corpus data (real Hungarian ik-verbs) from the Racz2024 repository
train = read_tsv('https://raw.githubusercontent.com/petyaracz/Racz2024/refs/heads/main/resource/real_words/ik_verbs/ikes_pairs_webcorpus2.tsv')
# Download experimental data (nonword ratings) from the Racz2024 repository
test = read_tsv('https://raw.githubusercontent.com/petyaracz/Racz2024/refs/heads/main/exp_data/baseline/baseline_tidy_proc.tsv')

# -- transcribe -- #

# Map Hungarian orthographic digraphs to single phonological symbols so that
# each character corresponds to exactly one segment:
#   cs → č, sz → ß (temp), zs → ž, ty → ṯ, gy → ḏ, ny → ṉ, ly → j
#   s  → š, ß  → s  (completes the two-step sz → s substitution)
train = train |> 
  mutate(
    lemma = str_replace_all(base, c('cs' = 'č', 'sz' = 'ß', 'zs' = 'ž', 'ty' = 'ṯ', 'gy' = 'ḏ', 'ny' = 'ṉ', 'ly' = 'j', 's' = 'š', 'ß' = 's')),
  )

test = test |> 
  mutate(
    lemma = str_replace_all(base, c('cs' = 'č', 'sz' = 'ß', 'zs' = 'ž', 'ty' = 'ṯ', 'gy' = 'ḏ', 'ny' = 'ṉ', 'ly' = 'j', 's' = 'š', 'ß' = 's')),
  )

# -- setup --#

# Keep only the lakok/lakom alternation (one of several ik-verb variation types)
test = test[test$variation == 'lakok/lakom',]

# Convert log-odds to probability: p = 1 / (1 + exp(-log_odds))
train$p = plogis(train$log_odds)
test$p = plogis(test$log_odds)


# Collect unique lemma forms from both sets for phonological distance
# computation via JANET (see dat/word_distances.tsv.gz)
out1 = train |> 
  select(lemma)

out1 = test |> 
  select(lemma) |> 
  bind_rows(out1) |> 
  distinct(lemma)

# length(unique(out1$lemma))
# length(unique(test$base_tr))

# -- write -- #

write_tsv(out1, 'dat/lakok_forms.tsv')
write_tsv(train, 'dat/lakok_train.tsv')
write_tsv(test, 'dat/lakok_test.tsv')
