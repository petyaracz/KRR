# -- head -- #

setwd('~/Github/KRR')

library(tidyverse)

# -- fun -- #

# ipa to orth
transcribe3 = function(transcription){
  str_replace_all(transcription,
                  c(
                    'A' = 'ā',
                    'E' = 'ē',
                    'I' = 'ī',
                    'O' = 'ō',
                    'U' = 'ū',
                    'f' = 'wh',
                    'N' = 'ng'
                  )
  )
}

# orth to moraic """IPA"""
transcribe2 = function(orthography){
  str_replace_all(orthography,
                  c(
                    'ā' = 'aa',
                    'ē' = 'ee',
                    'ī' = 'ii',
                    'ō' = 'oo',
                    'ū' = 'uu',
                    'wh' = 'f',
                    'ng' = 'N'
                  )
  )
}

# -- read -- #

source("script/functions.R")

train = read_csv('~/Github/HarrisRacz2026maori/raw/mbc_passives_all.csv')
test = read_csv('~/Github/HarrisRacz2026maori/raw/maonze_passives_anonymised.csv')
dist = read_tsv('~/Github/HarrisRacz2026maori/dat/word_distances.tsv.gz')

# -- setup --#

train = train |> 
  rename(ipa_root = root) |> 
  mutate(
    orth_root = transcribe3(ipa_root),
    dist_root = transcribe2(orth_root)
  ) |> 
  count(ipa_root,orth_root,dist_root,suffix) |> 
  mutate(
    root_fq = sum(n),
    .by = c(ipa_root,orth_root,dist_root)
  ) |> 
  filter(root_fq > 10, suffix %in% c('a', 'Nia')) |> 
  mutate(
    p_suffix = n / root_fq
  ) |> 
  select(-n) |> # !!!
  pivot_wider(names_from = suffix, values_from = p_suffix, names_prefix = 'observed_', values_fill = 0)

test = test |> 
  rename(ipa_root = root) |> 
  mutate(
    orth_root = transcribe3(ipa_root),
    dist_root = transcribe2(orth_root),
    speaker_group = case_when(
      str_starts(konae, "MU") | str_starts(konae, "WU") ~ "historical elder",
      str_starts(konae, "L1") ~ "L1",
      str_starts(konae, "L2") ~ "L2",
      str_starts(konae, "K") | str_starts(konae, "R") ~ "present elder"
    )
  ) |> 
  filter(!is.na(speaker_group)) |> 
  count(ipa_root,orth_root,dist_root,suffix,speaker_group) |> 
  mutate(
    root_fq = sum(n),
    .by = c(ipa_root,orth_root,dist_root)
  ) |> 
  filter(suffix %in% c('a', 'Nia')) |> 
  mutate(
    p_suffix = n / root_fq
  ) |> 
  filter(dist_root %in% dist$word1)

test_words_only = test |> distinct(dist_root)

# -- check -- #

sample_n(dist, 10)
train |> select(dist_root) |> sample_n(10)
test |> select(dist_root) |> sample_n(10)
summary(train$dist_root %in% dist$word1)
summary(train$dist_root %in% dist$word2)
summary(test$dist_root %in% dist$word1)
summary(test$dist_root %in% dist$word2)

# -- outcomes -- #

# -- model -- #

# a
m1 <- train_krr(train, dist, word_col = "dist_root", outcome_col = "observed_a", link = "logit", criterion = 'rmse')
m1$sigma   # best sigma
m1$alpha   # best alpha
m1$best_score

m1$predictions |> 
  ggplot(aes(observed,predicted_loo)) +
  geom_point() +
  geom_smooth()

with(m1$predictions, cor.test(observed,predicted_loo))

trainb = train[train$observed_a > 0,]

m1b <- train_krr(trainb, dist, word_col = "dist_root", outcome_col = "observed_a", link = "logit", criterion = 'rmse')
m1b$sigma   # best sigma
m1b$alpha   # best alpha
m1b$best_score

m1b$predictions |> 
  ggplot(aes(observed,predicted_loo)) +
  geom_point() +
  geom_smooth()

with(m1b$predictions, cor.test(observed,predicted_loo))

# maybe we're actually doing classification of "this is a". worth looking into.

# tia
m2 <- train_krr(train, dist, word_col = "dist_root", outcome_col = "observed_Nia", link = "logit", criterion = 'rmse')
m2$sigma   # best sigma
m2$alpha   # best alpha
m2$best_score

m2$predictions |> 
  ggplot(aes(observed,predicted_loo)) +
  geom_point() +
  geom_smooth(method = 'lm')

with(m2$predictions, cor.test(observed,predicted_loo))

trainc = train[train$observed_Nia > 0,]
m2b <- train_krr(trainc, dist, word_col = "dist_root", outcome_col = "observed_Nia", link = "logit", criterion = 'rmse')
m2b$sigma   # best sigma
m2b$alpha   # best alpha
m2b$best_score

m2b$predictions |> 
  ggplot(aes(observed,predicted_loo)) +
  geom_point() +
  geom_smooth(method = 'lm')

with(m2$predictions, cor.test(observed,predicted_loo))

# doesn't matter

# let's predict on test

# -- predict -- #

predict_a = predict_krr(train, test_words_only, dist, word_col = 'dist_root', outcome_col = 'observed_a', sigma = m1$sigma, alpha = m1$alpha, link = 'logit')

predict_Nia = predict_krr(train, test_words_only, dist, word_col = 'dist_root', outcome_col = 'observed_Nia', sigma = m2$sigma, alpha = m2$alpha, link = 'logit')

# -- combine -- #

predict_a = predict_a |> 
  rename(predict_a = predicted)

predict_Nia = predict_Nia |> 
  rename(predict_Nia = predicted)

predictions_a_Nia = left_join(predict_a,predict_Nia) |> 
  pivot_longer(-dist_root, names_to = 'suffix', values_to = 'predicted')
  
# -- write -- #

predictions_a_Nia |> 
  write_tsv('~/Github/HarrisRacz2026maori/dat/predictions_a_Nia_from_mbc_to_maonze.tsv')
