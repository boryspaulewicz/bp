library(bp)
library(lme4)
library(ggplot2)
data(sleepstudy)

fit = robust_mixed(Reaction ~ Days, Subject ~ Days, sleepstudy, beta_sigma = 1000, family = 'normal',
                   pars = c('ranef', 'y_new'))

## Rhat should be <= 1.01, neff >= 1000 for 95% CI
round(fit$summary, 2)

m <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy)
print(summary(m), corr = F)

## Bardzo podobne współczynniki efektów ustalonych
round(data.frame(lmer = fixef(m), robust = apply(fit$s[,fit$fixef], 2, mean)), 2)
##               lmer robust
## (Intercept) 251.41 252.65
## Days         10.47  10.75

ranef.1 = data.frame(x = ranef(m)$Subject[,1], y = apply(fit$s[, paste(unique(sleepstudy$Subject), '(Intercept)', sep = '.')], 2, mean))
ranef.2 = data.frame(x = ranef(m)$Subject[,2], y = apply(fit$s[, paste(unique(sleepstudy$Subject), 'Days', sep = '.')], 2, mean))
lm(y ~ x, ranef.1)
lm(y ~ x, ranef.2)
## Nachylenia bliskie 1, ale widać odstające przypadki
ggplot(ranef.1, aes(x, y)) + geom_point() + geom_smooth(method = 'lm') + geom_text(nudge_y = 1, size = 1.5, aes(label = unique(sleepstudy$Subject)))
## 308, 332, 335, 352
ggplot(ranef.2, aes(x, y)) + geom_point() + geom_smooth(method = 'lm') + geom_text(nudge_y = 1, size = 1.5, aes(label = unique(sleepstudy$Subject)))
## 332, 352, 308

## Większy poziom niepewności z bayesa
tbl = lmer_sig(m, T)
round(data.frame(lmer = as.numeric(tbl[,3]), robust = apply(fit$s[,fit$fixef], 2, sd)), 2)
##             lmer robust
## (Intercept) 6.82   8.10
## Days        1.55   1.74

## Większa wariancja efektów losowych z bayesa, korelacja z odwrotnym
## znakiem (sic!), ale jest "nieistotna"
round(apply(fit$s, 2, mean)[c(3:4, 6)], 2)
## (Intercept) SD        Days SD          C.2.1
##          27.93           6.04          -0.11
##
## Random effects:
##  Groups   Name        Std.Dev. Corr
##  Subject  (Intercept) 24.740
##           Days         5.922   0.07

## Dopasowanie na podstawie posteriora predykcyjnego
y_new = t(apply(fit$s[,rmatch('y_new', names(fit$s))], 2, function(x)c(quantile(x, .025), mean(x), quantile(x, .975))))
colnames(y_new) = c('lo', 'y', 'hi')
ggplot(cbind(sleepstudy, y_new), aes(x = Days, y = Reaction)) + geom_point() +
    geom_line(color = 'red', aes(y = fitted(m))) + geom_line(aes(y = y)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = .1) + facet_wrap(~Subject)
## Pięknie olewa obserwacje odstające: 308, 332, 335, 352

######################################################################
## Robit

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

## nu = 30, żeby przybliżyć zwykłą mieszaną regresję logistyczną
fit = robust_mixed(acc ~ -1 + gr / x, id ~ x, df, n = 20, family = 'binomial',
                   pars = 'y_new', y_nu = 30, ranef_nu = 30)

## Rhat <= 1.01, neff >= 1000 for 95% CI
round(fit$summary, 2)

## Porównujemy efekty ustalone z wartościami prawdziwymi
round(rbind(apply(fit$s[,fit$fixef], 2, mean), c(-.5, .5, 1, 2)), 2)

## Współczynniki bardzo podobne
m <- glmer(cbind(acc,n-acc) ~ -1 + gr / x + (x|id), df, family = 'binomial')
round(rbind(coef(summary(m))[fit$fixef, 1],
            apply(fit$s[,fit$fixef], 2, mean)), 2)

## Błędy standardowe bardzo podobne
round(coef(summary(m))[fit$fixef, 2] / apply(fit$s[,fit$fixef], 2, sd), 2)
