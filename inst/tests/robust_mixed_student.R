library(bp)
library(lme4)
library(ggplot2)
data(sleepstudy)

## Testujemy zwykle najbardziej aktualną wersję
source('~/cs/code/r/bp/R/bayes.R')

fit = robust_mixed(Reaction ~ Days, Subject ~ Days, sleepstudy, beta_sigma = 1000, type = 'student',
                   pars = c('ranef', 'y_new'))

## Rhat should be <= 1.01, neff >= 1000 for 95% CI
round(fit$summary, 2)

plot(fit$s$y_nu)
## Pięknie się próbkuje
hist(fit$s$y_nu)

m <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy)
print(summary(m), corr = F)

## Bardzo podobne współczynniki efektów ustalonych
round(data.frame(lmer = fixef(m), robust = apply(fit$s[,fit$fixef], 2, mean)), 2)
              lmer robust
(Intercept) 251.41 252.67
Days         10.47  10.74

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
## (Intercept) 6.82   8.42
## Days        1.55   1.73

## Większa wariancja efektów losowych z bayesa, korelacja z odwrotnym
## znakiem (sic!), ale jest "nieistotna"
round(apply(fit$s, 2, mean)[c(3:4, 7)], 2)
## (Intercept) SD        Days SD          C.2.1 
##          28.58           6.07          -0.11 
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
