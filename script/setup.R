# -- head -- #

setwd('~/Github/KRR')

library(tidyverse)

# -- read -- #

train = read_tsv('https://raw.githubusercontent.com/petyaracz/Racz2024/refs/heads/main/resource/real_words/ik_verbs/ikes_pairs_webcorpus2.tsv')
test = read_tsv('https://raw.githubusercontent.com/petyaracz/Racz2024/refs/heads/main/exp_data/baseline/baseline_tidy_proc.tsv')

# -- transcribe -- #

train = train |> 
  mutate(
    lemma = str_replace_all(base, c('cs' = 'č', 'sz' = 'ß', 'zs' = 'ž', 'ty' = 'ṯ', 'gy' = 'ḏ', 'ny' = 'ṉ', 'ly' = 'j', 's' = 'š', 'ß' = 's')),
  )

test = test |> 
  mutate(
    lemma = str_replace_all(base, c('cs' = 'č', 'sz' = 'ß', 'zs' = 'ž', 'ty' = 'ṯ', 'gy' = 'ḏ', 'ny' = 'ṉ', 'ly' = 'j', 's' = 'š', 'ß' = 's')),
  )

# -- setup --#

# narrow
test = test[test$variation == 'lakok/lakom',]

# outcome
train$p = plogis(train$log_odds)
test$p = plogis(test$log_odds)


# forms
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
