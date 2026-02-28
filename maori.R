# -- head -- #

setwd('~/Github/KRR')

library(tidyverse)

# -- fun -- #

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

train = read_tsv('~/Github/HarrisRacz2026maori/dat/maori_roots_suffixes_corpus.tsv')
dist = read_tsv('~/Github/HarrisRacz2026maori/dat/word_distances.tsv.gz')

# -- setup --#

train = train |> 
  filter(root_fq > 10) |> 
  group_by(root) |> 
  arrange(-root_fq) |> 
  slice(1) |> 
  ungroup() |> 
  mutate(transcription = transcribe2(root))

# -- check -- #

sample_n(dist, 10)
train |> select(root,transcription) |> sample_n(10)
summary(train$transcription %in% dist$word1)
summary(train$transcription %in% dist$word2)

# erm
train = train |> 
  filter(transcription %in% dist$word1)

# -- outcomes -- #

train = train |> 
  mutate(
    p_a = a / root_fq,
    p_tia = tia / root_fq,
    p_hia = hia / root_fq,
    p_ngia = ngia / root_fq
  )

# -- model -- #

# a
m1 <- train_krr(train, dist, word_col = "transcription", outcome_col = "p_a", link = "logit", criterion = 'rmse')
m1$sigma   # best sigma
m1$alpha   # best alpha
m1$best_score

m1$predictions |> 
  ggplot(aes(observed,predicted_loo)) +
  geom_point() +
  geom_smooth()

with(m1$predictions, cor.test(observed,predicted_loo))

# tia
m2 <- train_krr(train, dist, word_col = "transcription", outcome_col = "p_tia", link = "logit", criterion = 'rmse')
m2$sigma   # best sigma
m2$alpha   # best alpha
m2$best_score

m2$predictions |> 
  ggplot(aes(observed,predicted_loo)) +
  geom_point() +
  geom_smooth()

with(m2$predictions, cor.test(observed,predicted_loo))

# hia
m3 <- train_krr(train, dist, word_col = "transcription", outcome_col = "p_hia", link = "logit", criterion = 'rmse')
m3$sigma   # best sigma
m3$alpha   # best alpha
m3$best_score

m3$predictions |> 
  ggplot(aes(observed,predicted_loo)) +
  geom_point() +
  geom_smooth()

with(m3$predictions, cor.test(observed,predicted_loo))

# ngia
m4 <- train_krr(train, dist, word_col = "transcription", outcome_col = "p_ngia", link = "logit", criterion = 'rmse')
m4$sigma   # best sigma
m4$alpha   # best alpha
m4$best_score

m4$predictions |> 
  ggplot(aes(observed,predicted_loo)) +
  geom_point() +
  geom_smooth()

with(m4$predictions, cor.test(observed,predicted_loo))
