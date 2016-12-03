library(bp)
library(lme4)

df = expand.grid(id = 1:40, trial = 1, x = 0:1)
df$gr = 1
df$gr[df$id > 20] = 2
df$acc = NA
x.eff = ct(rnorm(max(df$id), sd = 1))
intercept.eff = ct(rnorm(max(df$id), sd = 0.5))
n = 20
for(i in 1:nrow(df))df$acc[i] = rbinom(1, size = n, binomial()$linkinv(c(-.5, .5)[df$gr[i]] + intercept.eff[df$id[i]] +
                                                                             (df$gr[i] + x.eff[df$id[i]]) * df$x[i]))
df$gr = as.factor(df$gr)

## Testujemy zwykle najbardziej aktualną wersję
source('~/cs/code/r/bp/R/robust_mixed.R')

## nu = 30, żeby przybliżyć zwykłą mieszaną regresję logistyczną
fit = robust_mixed(acc ~ -1 + gr / x, id ~ x, df, n = 20, type = 'robit',
                   pars = 'y_new', y_nu = 30, ranef_nu = 30)

round(fit$summary, 2)

## Porównujemy efekty ustalone z wartościami prawdziwymi
round(rbind(apply(fit$samples[,fit$fixef], 2, mean), c(-.5, .5, 1, 2)), 2)

## Współczynniki bardzo podobne
m <- glmer(cbind(acc,n-acc) ~ -1 + gr / x + (x|id), df, family = 'binomial')
round(rbind(coef(summary(m))[fit$fixef, 1],
            apply(fit$samples[,fit$fixef], 2, mean)), 2)

## Błędy standardowe bardzo podobne
round(coef(summary(m))[fit$fixef, 2] / apply(fit$samples[,fit$fixef], 2, sd), 2)
