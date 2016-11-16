#' Poziomy istotności dla mieszanego modelu liniowego
#' @param m Model dopasowany za pomocą funkcji lmer.
#' @param return Czy zwracać tabelkę regresji jako obiekt.
#' @param short Czy drukować tabelkę w wersji skróconej.
#' @param show Czy drukować tabelkę.
#' @return Zwraca tabelkę z istotnościami, jeżeli return ==
#' TRUE. Drukuje tabelkę z istotnościami na wyjściu, jeżeli show ==
#' TRUE.
#' @export
lmer.sig = function(m, return = FALSE, short = FALSE, show = TRUE){
  ms = summary(m)
  ## fef = as.data.frame(ms@coefs)
  fef = as.data.frame(coef(ms))
  rownames(fef) = gsub('\\(Intercept\\)', 'Intercept', rownames(fef))
  names(fef) = c('Coef.', 'Coef. SE', 't')
  ## tdf = nrow(m@X) - rankMatrix(cbind(m@X, (as.matrix(t(m@Zt)))))
  ## tdf = nrow(model.matrix(m)) - rankMatrix(cbind(m@X, (as.matrix(t(m@Zt)))))
  tdf = nrow(model.matrix(m)) - nrow(fef) - ms$devcomp$dims['q']
  ## To co ponizej jest raczej liberalne:
  ## tdf = nrow(m@X) - ncol(m@X):
  ##
  ## If you would like to use a t-distribution to calculate a confidence
  ## interval, I would argue that the degrees of freedom for the t should
  ## be somewhere between n - p = 178 in this case (n = number of
  ## observations, p = number of coefficients in the fixed effects or
  ## rank(X) where X is the model matrix for the fixed effects) and n -
  ## rank([Z:X]) = 144.  (there are 36 random effects and 2 fixed effects
  ## but the rank of [Z:X] = 36)
  ## If we use the lower value of 144 we produce a confidence interval of
  ##
  ## > 10.4673 + c(-1,1) * qt(0.975, df = 144) * 1.5458
  ## [1]  7.41191 13.52269
  ##
  ## Notice that this interval is contained in the HPD interval and in the
  ## interval obtained from the 2.5% and 97.5% quantiles of the empirical
  ## distribution.  I attribute this to the fact that the standard errors
  ## are calculated conditional on the variance of the random effects.
  ## Thus the t-based confidence interval takes into account imprecision of
  ## the estimate of \sigma^2 but it assumes that the variance of the
  ## random effects is fixed at the estimated values.
  ##
  ## Doug Bates w odpowiedzi na post na forum:
  ## https://stat.ethz.ch/pipermail/r-help/2006-August/110736.html
  pval = 2*pt(-abs(fef$t), tdf)
  ## 2*pnorm(-abs(fef[,1] / fef[,2]))
  fef[,1:3] = apply(fef[,1:3], 2, function(x)sprintf('%.2f', x))
  fef$p = sprintf('%.3f', pval)
  stars = rowSums(cbind(pval < .05, pval < .01, pval < .001))
  fef$sig = ''
  fef$sig[stars == 3] = '***'
  fef$sig[stars == 2] = '**'
  fef$sig[stars == 1] = '*'
  ## Formatujemy macierz efektoww losowych
  ## remat = ms@REmat ## juz nie dziala
  ## remat = ms$varcor
  ## Przerabiamy na dokladnosc do 2 miejsc po przecinku, zmieniamy nazwy
  ## kolumn etc.
  ##
  ## if(ncol(remat) == 4){ colnames(remat)[3:4] = c('Variance', 'SD')
  ## }else if(ncol(remat) == 3){
  ##     colnames(remat)[3] = 'SD'
  ## }else{ colnames(remat)[3:5] = c('Variance', 'SD', 'Correlations') }
  ## for(r in 1:nrow(remat)){
  ##   for(c in 1:ncol(remat)){
  ##     if(rmatch('[0-9]+\\.[0-9]+', remat[r,c])){
  ##       remat[r,c] = sprintf("%.2f", as.numeric(remat[r,c]))
  ##     }}}
  ## remat[,2] = gsub('\\(Intercept\\)', 'Intercept', remat[,2])
  ## remat[,1:2] = apply(remat[,1:2], 2, function(x)format(x, justify = 'left'))
  ## ## Pokazujemy wszystko
  ## if(show){
  ##   cat('\nRandom effects\n\n')
  ##   print.table(remat, justify = 'right')
  ##   cat('\nFixed effects\n\n')
  ## }
  if(show){
      if(short){
          print(fef[,-c(2:4)])
      }else{
          print(fef)
      }
      cat('\n\nt df: ', tdf[1], '\n')
  }
  if(return)cbind(rownames(fef), fef) ## list(fef = fef, ref = remat)
}

#' Funkcja zwracająca ładną prostą tabelkę z gwiazdkami dla modelu
#' liniowego
#' @export
lms = function(..., digits = 3, return.table = FALSE){
  sm = summary(lm(...))
  ## tylko coef, se i p
  res = as.data.frame(round(sm$coef,digits)[,c(1,2,4)])
  names(res) = c('Współcz.', 'Błąd std.', 'p')
  p = summary(lm(...))$coef[,4]
  res$sig = ifelse(p <= 0.001, '***', ifelse(p <= 0.01, '**', ifelse(p <= 0.05, '*', '')))
  print(formula(lm(...)))
  cat('\n')
  print(res)
  cat(paste('\nk =', sm$df[1],
            'n = ', sum(sm$df[1:2]),
            'błąd std. reszt =', round(sqrt(sum((sm$residuals - mean(sm$residuals))^2)/sm$df[2]),digits),
            '\n'))
  if(return.table)res
}

#' Oblicza macierz kontrastów dla tabelki (summary) modelu
#'
#' Zakładamy, że druga kolumna ramki (m) zawiera błędy standardowe i
#' dla każdego porównania współczynników stosujemy średnią
#' geometryczną ich błędów standardowych jako błąd standardowy
#' różnicy.
#' @export
contrast.matrix = function(m, rnd = 2, draw.stars = TRUE, fun = I){
  m = as.data.frame(m)
  nms = rownames(m)
  res = sigtbl = matrix(nrow = nrow(m), ncol = nrow(m))
  rownames(res) = nms
  colnames(res) = nms
  for(n in 1:length(nms)){
    dif = m[n,1] - m[,1]
    se = sqrt(m[n,2]^2 + m[,2]^2)
    sig = 2 * pnorm(-abs(dif / se))
    stars = rep('', nrow(m))
    stars[sig < 0.05] = '*'
    stars[sig < 0.01] = '**'
    stars[sig < 0.001] = '***'
    ## W ten sposób można np. zastosować inv.logit na interceptach
    val = round(fun(m[n,1]) - fun(m[,1]), rnd)
    if(draw.stars){
      res[,n] = paste(val, stars)
    }else{
      res[,n] = paste(val, 'p = ', round(sig, 3))
    }
    sigtbl[,n] = sig
  }
  ## Wymazujemy przekątną
  for(n in 1:nrow(res)){
    res[n:nrow(res),n] = ''
  }
  list(tbl = res, sig = sigtbl)
}

#' Przybliżone przedziały ufności dla proporcji
#'
#' Rozkład normalny może być kiepskim przybliżeniem przedziałów
#' ufności (w wikipedii jest odnośnik do artykułu opisującego kiedy ta
#' reguła nie działa). Generalnie proporcja nie może być zbyt blisko 0
#' ani 1. Istnieje wiele innych metod liczenia interwałów, wiki sporo
#' o tym mówi.
#' @param fits Wektor wartości oczekiwanych (prawdopodobieństw)
#' @param n Wektor wielkości prób
#' @param alpha 1 - prawdopodobieństwo pokrycia
#' @export
binomial.pred.ci = function(fits, n, alpha = .05){
  if(sum(n * fits <= 5 | n * (1 - fits) <= 5) != 0){
    warning("Rule of thumb for normal approximation violated, using Wilson score interval")
    z = qnorm(1-alpha/2)
    est = (fits + z^2/(2*n)) / (1 + z^2 / n)
    ci.hi = z * sqrt(fits * (1 - fits) / n + z^2 / (4 * n^2)) / (1 + z^2 / n)
    ci.lo = -1 * ci.hi
    ci.lo = est + ci.lo - fits
    ci.hi = est + ci.hi - fits
  }else{
    ci.hi = qnorm(1-alpha/2) * sqrt((fits*(1-fits))/n)
    ci.lo = -1 * ci.hi
  }
  data.frame(ci.hi = ci.hi, ci.lo = ci.lo, fits = fits, n = n)
}
