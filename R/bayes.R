#' Dopasowuje odporny mieszany model liniowy lub logistyczny
#'
#' Implementacja modelu jest oparta na publikacji Sorensen'a,
#' Hohenstain'a i Vasishth'a (2015): Bayesian linear mixed models
#' using Stan: A tutorial for psychologists, linguists, and cognitive
#' scientists. Warto pamiętać, że prior dla efektów ustalonych nie
#' jest płaski (domyślnie beta_sigma = sd prioru normalnego == 20).
#' @param fixed Formuła modelu dla efektów ustalonych.
#' @param random Formuła modelu dla efektów losowych, postaci g ~ x1 +
#'     x2 ...
#' @param d Zbiór danych.
#' @param chains Liczba równoległych próbników - domyślnie liczba
#'     rdzeni.
#' @param pars Wektor tekstowy nazw parametrów, dla których chcemy
#'     próbki, poza beta, ranef_sigma i C. Może zawierać 'y_new'
#'     (posterior predykcyjny) i 'ranef' (efekty losowe).
#' @param y_nu Liczba stopni swobody rozkładu t modelującego zmienną
#'     zależną.
#' @param y_sigma (Uwzględniane tylko dla modelu logistycznego)
#'     Odchylenie standardowe rozkładu t modelującego zmienną
#'     zależną. To jest wartość arbitralna - domyślnie przyjęto
#'     wartość dającą współczynniki porównywalne do zwykłej regresji
#'     logistycznej.
#' @param ranef_nu Liczba stopni swobody rozkładu t modelującego
#'     rozkład efektów losowych.
#' @param beta_sigma Odchylenie standardowe prioru dla efektów
#'     ustalonych. UWAGA Jeżeli predyktory mają bardzo małe odchylenie
#'     standardowe, albo oczekiwane prawdopodobieństwa dla punktów
#'     przecięcia są skrajne, domyślna wartość może być
#'     nieodpowiednia.
#' @return Wynik próbkowania zwrócony przez funkcję stan z dodanym
#'     atrybutem fixed.effects, zawierającym nazwy kolumn macierzy
#'     modelu.
#' @export
robust.mixed = function(fixed, random, d,
                        y_nu = 4, y_sigma = 1.548435, ranef_nu = 4, beta_sigma = 20,
                        chains = parallel::detectCores(), pars = NULL, family = 'binomial', ...){
    require(rstan)
    if(family == 'binomial'){
        if(!('n' %in% names(d)))stop("Number of observations per data point variable (n) missing")
        model.path = 'stan_models/robust_mixed_logistic.stan'
    }
    if(family == 'normal'){
        if(beta_sigma == 20)warning('SD of normal prior for beta = 20 in the linear model')
        model.path = 'stan_models/robust_mixed_linear.stan'
    }
    pars = c('beta', 'ranef_sigma', 'C', pars)
    if(!all(pars %in% c('beta', 'ranef_sigma', 'C', 'ranef', 'y_new'))){
        error('Podano niewłaściwe nazwy dodatkowych parametrów: dopuszczalne wartości to \'ranef\' i \'y_new\'.')
    }
    rstan_options(auto_write = TRUE)
    options(mc.cores = parallel::detectCores())
    id = as.numeric(as.factor(as.character(d[[as.character(random[2])]])))
    y = d[,as.character(fixed[2])]
    X = model.matrix(fixed, d)
    Z = model.matrix(random, d)
    data = list(D = ncol(X), R = ncol(Z), N = nrow(X), I = max(id),
                X = X, Z = Z, y = y, id = id,
                y_nu = y_nu, ranef_nu = ranef_nu, beta_sigma = beta_sigma)
    if(family == 'binomial'){
        data$n = d$n
        data$y_sigma = y_sigma
    }
    fit = stan(paste(path.package('bp'), model.path, sep = '/'), data = data,
               chains = chains, pars = pars, ...)
    ## Zwracamy próbki z czytelnymi nazwami parametrów
    s = as.data.frame(extract(fit))
    names(s)[rmatch('beta', names(s))] = colnames(X)
    names(s)[rmatch('ranef_sigma', names(s))] = paste(colnames(Z), 'SD')
    ranef_names = apply(expand.grid(unique(as.character(d[[as.character(random[2])]])), colnames(Z)), 1, function(x)paste(x, collapse = '.'))
    ## names(ranef_names) = names(s)[rmatch('ranef.', names(s))]
    if('ranef' %in% pars)names(s)[rmatch('ranef.', names(s))] = ranef_names
    stan.sum = rstan::summary(fit)$summary
    ## rownames(stan.sum) = names(s) ## nazwy efektów losowych w summary są w innej kolejności
    list(s = s, summary = stan.sum, ranef = ranef_names, fixef = colnames(X))
}

#' Test funkcji robust.mixed.logistic
#' @export
test.robust.mixed.logistic = function(){
    require(lme4)
    df = expand.grid(id = 1:40, trial = 1, x = 0:1)
    df$gr = 1
    df$gr[df$id > 20] = 2
    df$n = 20
    df$acc = NA
    x.eff = rnorm(max(df$id), sd = 1)
    x.eff = x.eff - mean(x.eff)
    intercept.eff = rnorm(max(df$id), sd = 0.5)
    intercept.eff = intercept.eff - mean(intercept.eff)
    for(i in 1:nrow(df))df$acc[i] = rbinom(1, size = df$n[i], binomial()$linkinv(c(-.8, .8)[df$gr[i]] + intercept.eff[df$id[i]] +
                                                                        (df$gr[i] + x.eff[df$id[i]]) * df$x[i]))
    df$gr = as.factor(df$gr)
    fit = robust.mixed.logistic(acc ~ -1 + gr / x, id ~ x, df,
        y_nu = 100, ranef_nu = 100)
    sm = round(summary(fit)$summary, 3)
    ## Porównujemy współczynniki dostarczone przez glmer i te bayesowskie
    round(rbind(coef(summary(glmer(cbind(acc,n-acc) ~ -1 + gr / x + (x|id), df, family = 'binomial')))[1:4, 1],
                sm[1:4, 'mean']), 3)
    ## Porównujemy błędy standardowe do sd posteriorów
    round(coef(summary(glmer(cbind(acc,n-acc) ~ -1 + as.factor(gr) / x + (x|id), df, family = 'binomial')))[1:4, 2] /
              sm[1:4, 'sd'], 1)
    ## W symulacjach niepewność była wyższa w przypadku odpornego
    ## modelu bayesowskiego gdy nu = 4, a gdy nu = 100 (czyli
    ## zakładamy w przybliżeniu model prawdziwy), była praktycznie
    ## taka sama.
    fit
}
