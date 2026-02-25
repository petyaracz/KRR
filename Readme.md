## Kernel ridge regression

- take word distances (word1-word2) based on aligned phonological distance (see [here](github.com/petyaracz/JANET))
- take training data with words and continuous behaviour (e.g. English past tense forms and p(regular))
- use hyperparameter tuning with LOO to fit model that predicts word behaviour from behaviour of similar words
- use best model to predict test words