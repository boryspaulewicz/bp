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
#' @param n Liczba obserwacji na punkt danych (robit)
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
#' @param family (='binomial') Rozkład zmiennej zależnej: binomial
#'     oznacza odporną regresję logistyczną, normal oznacza odporny
#'     model liniowy.
#' @param return.stanfit (=F) Czy zwracać obiekt zwracany przez Stana,
#'     czy listę z ramką z próbkami i summary
#' @return Obiekt zwracany przez rstan lub lista złożona z elementów
#'     s: ramka próbek, summary: podsumowanie wyników STAN'a.
#' @export
robust.mixed = function(fixed, random, d, n = NULL,
                        y_nu = 4, y_sigma = 1.548435,
                        beta_sigma = NULL, ranef_nu = 4,
                        chains = parallel::detectCores() - 1,
                        pars = NULL, family = 'binomial',
                        return.stanfit = F, ...){
    require(rstan)
    if(family == 'binomial'){
        if(is.null(n))stop("Number of observations per data point (n) missing.")
        model.path = 'stan_models/robust_mixed_logistic.stan'
        if(is.null(beta_sigma)){
            warning('Using defaulf SD=20 for fixed effects priors. This could be inappropriate for unstandardized predictors.')
            beta_sigma = 20
        }
    }
    if(family == 'normal'){
        if(beta_sigma == 20)warning('SD of normal prior for beta = 20 in the linear model')
        if(is.null(beta_sigma))stop('beta_sigma (SD of fixed effects priors) not set')
        model.path = 'stan_models/robust_mixed_linear.stan'
    }
    if(!all(pars %in% c('ranef', 'y_new'))){
        error('Wrong parameter names. Valid values are \'ranef\' and \'y_new\'.')
    }
    pars = c(c('beta', 'ranef_sigma', 'C'), pars)
    rstan_options(auto_write = TRUE)
    options(mc.cores = parallel::detectCores())
    id = as.numeric(as.factor(as.character(d[[as.character(random[2])]])))
    y = d[,as.character(fixed[2])]
    X = model.matrix(fixed, d)
    Z = model.matrix(random, d)
    data = list(D = ncol(X), R = ncol(Z), N = nrow(X), I = max(id),
                X = X, Z = Z, y = y, id = id,
                y_nu = y_nu, ranef_nu = ranef_nu,
                beta_sigma = beta_sigma)
    if(family == 'binomial'){
        d$n = n ## Zapewnia odpowiednią długość, gdy n to skalar
        data$n = d$n
        data$y_sigma = y_sigma
    }
    fit = stan(paste(path.package('bp'), model.path, sep = '/'), data = data,
               chains = chains, pars = pars, ...)
    if(!return.stanfit){
        ## Zwracamy próbki z czytelnymi nazwami parametrów
        s = as.data.frame(extract(fit))
        names(s)[rmatch('beta', names(s))] = colnames(X)
        names(s)[rmatch('ranef_sigma', names(s))] = paste(colnames(Z), 'SD')
        ranef_names = apply(expand.grid(unique(as.character(d[[as.character(random[2])]])), colnames(Z)), 1, function(x)paste(x, collapse = '.'))
        ## names(ranef_names) = names(s)[rmatch('ranef.', names(s))]
        if('ranef' %in% pars)names(s)[rmatch('ranef.', names(s))] = ranef_names
        stan.sum = rstan::summary(fit)$summary
        stan.sum = stan.sum[!rmatch('y_new', rownames(stan.sum)),]
        ## rownames(stan.sum) = names(s) ## nazwy efektów losowych w summary są w innej kolejności
        list(s = s, summary = stan.sum, ranef = ranef_names, fixef = colnames(X))
    }else{
        fit
    }
}

#' Dopasowuje model służący do wykrywania osób wykonujących zadanie
#' przypadkowo
#'
#' @param score Wektor sum poprawnych odpowiedzi
#' @param n Wektor liczby prób
#' @param guessing_prob = .5, zakładane prawdopodobieństwo przypadkowego trafienia
#' @param ... Pozostałe argumenty dla funkcji run.jags
#' @return Lista z elementami fit i pred, gdzie fit to obiekt zwracany przez run.jags a pred to ramka z predykcjami posteriorycznymi
#' @export
find_guessers = function(score, n, guessing_prob = .5, n.chains = parallel::detectCores() - 1, sample = 20000,
                         ...){
    ## TODO: eff dla samych prób wybranych ze względu na is_guessing,
    ## dodać wygodną funkcję diagnostyczną
    d = list(score = score,
             N = length(score),
             n = n,
             guessing_prob = guessing_prob)
    ## Bardzo słabo miesza mu i sigma, niektóre not_guessing_prob też tak sobie
    model =
    'model{
    mu ~ dunif(0, 4.5) ## odpowiada .5-.99
    sigma ~ dunif(1.0E-3, 2.25)
    for (s in 1:N){
        is_guessing[s] ~ dbern(0.5)
        logit(not_guessing_prob[s]) <- logit_prob[s]
        logit_prob[s] ~ dnorm(mu, 1/sigma^2)
        theta[s] <- equals(is_guessing[s], 0) * not_guessing_prob[s] + equals(is_guessing[s], 1) * guessing_prob
        score[s] ~ dbin(theta[s], n[s])
        }
    }
    '
    fit = runjags::run.jags(model, c('mu', 'sigma', 'not_guessing_prob', 'is_guessing'), d,
                            n.chains = n.chains, sample = sample,
                            inits = function(chain)list(mu = runif(1, binomial()$linkfun(.6), binomial()$linkfun(.9)),
                                                        sigma = runif(1, .001, 2.2),
                                                        is_guessing = rbinom(length(score), 1, .5)),
                            modules = 'glm',
                            ...)
    s = as.data.frame(as.matrix(coda::as.mcmc.list(fit)))
    ## Zwracamy predykcje posterioryczne
    pred = expand.grid(theta = d$score, lo = NA, hi = NA)
    pred$is_guessing = apply(s[,rmatch('is_guessing', names(s))], 2, mean)
    pred$n = d$n
    pred$i = as.factor(1:nrow(pred))
    for(i in 1:nrow(pred)){
        if(pred$is_guessing[i] > .5){
            pred[i, c('theta', 'lo', 'hi')] = c(.5, qbinom(c(.025, .975), pred$n[i], .5) / pred$n[i])
        }else{
            ## Obliczamy poprawność i niepewność pod warunkiem, że
            ## prawdopodobnie niezgadujący nie zgadują
            samples = s[[sprintf('not_guessing_prob[%d]', i)]][s[[sprintf('is_guessing[%d]', i)]] == 0]
            pred[i, c('theta', 'lo', 'hi')] = c(mean(samples), quantile(samples, c(.025, .975)))
        }
    }
    pred$acc = score / pred$n
    list(fit = fit, pred = pred)
}

#' Wykres dopasowania modelu wykrywającego osoby wykonujące zadanie na poziomie przypadku
#'
#' @param fit Obiekt zwracany przez find_guessers
#' @export
plot_guessers = function(fit){
    fit$pred$is_guessing = c('Above chance', 'Chance')[as.numeric(fit$pred$is_guessing > .5) + 1]
    ggplot2::ggplot(fit$pred, aes(x = i, y = acc)) + ylab('Accuracy') +
        ggplot2::geom_point() + ggplot2::geom_errorbar(aes(ymin = lo, ymax = hi)) + facet_wrap(~is_guessing)
}
