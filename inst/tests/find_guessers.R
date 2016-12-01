library(bp)
## library(runjags)
## library(coda)
library(ggplot2)

score = c(21,17,21,18,22,31,31,34,34,35,35,36,39,36,35)
n = rep(40, length(score))
res = find_guessers(score, n)

res$fit

plot_guessers(res)
