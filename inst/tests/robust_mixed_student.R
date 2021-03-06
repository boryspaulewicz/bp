library(bp)
library(lme4)
library(ggplot2)
data(sleepstudy)

## Testujemy zwykle najbardziej aktualną wersję
source('~/cs/code/r/bp/R/robust_mixed.R')

fit = robust_mixed(Reaction ~ Days, Subject ~ Days, sleepstudy, beta_sigma = 1000, type = 'student',
                   pars = c('ranef', 'y_new'), y_nu_rate = 1/29, ranef_nu_rate = 1/29)

round(fit$summary, 2)
## ranef_nu próbkuje się znacznie gorzej, niż pozostałe parametry i
## nic dziwnego (mało osób badanych, ciężar ogonów odpowiada rzadkim i
## odległym obserwacjom). Ponieważ jego średnia posterioryczna > 30,
## można założyć, że efekty losowe mają rozkład normalny.

fit$summary[fit$ranef, 'mean'] / apply(fit$samples[,fit$ranef], 2, mean)
## Nazwy efektów losowych w ramce próbek i summary są zgodne.

m <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy)
print(summary(m), corr = F)

## Bardzo podobne współczynniki efektów ustalonych
round(data.frame(lmer = fixef(m), robust = apply(fit$samples[,fit$fixef], 2, mean)), 2)

## Porównujemy efekty losowe z lmer i robit
ranef.1 = data.frame(x = ranef(m)$Subject[,1], y = apply(fit$samples[, paste(unique(sleepstudy$Subject), '(Intercept)', sep = '.')], 2, mean))
ranef.2 = data.frame(x = ranef(m)$Subject[,2], y = apply(fit$samples[, paste(unique(sleepstudy$Subject), 'Days', sep = '.')], 2, mean))
lm(y ~ x, ranef.1)
lm(y ~ x, ranef.2)
## Nachylenia bliskie 1, ale widać odstające przypadki
ggplot(ranef.1, aes(x, y)) + geom_point() + geom_smooth(method = 'lm') + geom_text(nudge_y = 1, size = 1.5, aes(label = unique(sleepstudy$Subject)))
## 308, 332, 335, 352
ggplot(ranef.2, aes(x, y)) + geom_point() + geom_smooth(method = 'lm') + geom_text(nudge_y = 1, size = 1.5, aes(label = unique(sleepstudy$Subject)))
## 332, 352, 308

## Większy poziom niepewności z bayesa
tbl = lmer_sig(m, T)
round(data.frame(lmer = as.numeric(tbl[,3]), robust = apply(fit$samples[,fit$fixef], 2, sd)), 2)

## Większa wariancja efektów losowych z bayesa, korelacja z odwrotnym
## znakiem (sic!), ale jest "nieistotna"
round(apply(fit$samples, 2, mean)[c(3:4, 7)], 2)
## Random effects:
##  Groups   Name        Std.Dev. Corr
##  Subject  (Intercept) 24.740
##           Days         5.922   0.07

## Dopasowanie na podstawie posteriora predykcyjnego
y_new = t(apply(fit$samples[,rmatch('y_new', names(fit$samples))], 2, function(x)c(quantile(x, .025), mean(x), quantile(x, .975))))
colnames(y_new) = c('lo', 'y', 'hi')
ggplot(cbind(sleepstudy, y_new), aes(x = Days, y = Reaction)) + geom_point() +
    geom_line(color = 'red', aes(y = fitted(m))) + geom_line(aes(y = y)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = .1) + facet_wrap(~Subject)
## Wydaje się dobrze sobie radzić z odstającymi: 308, 332, 335, 352
