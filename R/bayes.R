## Zmienia wartość fragmentów modelu na podstawie listy POLE =
## 'wartość'
set = function(fields){
    res = paste(readLines(paste(path.package('bp'), '/stan_models/robust_mixed.stan', sep = ''), encoding = 'utf8'),
                collapse = '\n')
    for(f in names(fields))res = gsub(sprintf('// %s', f), fields[[f]], res)
    res
}

create_model = function(type, y_nu_free){
    model = switch(type,
                   robit = list(DATA = 'int<lower=0> y[N];\n real<lower=0> y_sigma;\n int<lower=1> n[N];\n',
                                PARAMETERS = '',
                                TRANSFORMED = 'for(i in 1:N){ eta[i] = student_t_cdf(X[i] * beta + Z[i] * ranef[id[i]], y_nu, 0, y_sigma); }\n',
                                MODEL = 'for(i in 1:N){ y[i] ~ binomial(n[i], eta[i]); }\n',
                                GENERATED = 'for(i in 1:N){ y_new[i] = binomial_rng(n[i], eta[i]); }\n'),
                   student = list(DATA = 'real y[N];\n',
                                PARAMETERS = 'real<lower=0> y_sigma;\n',
                                TRANSFORMED = 'for(i in 1:N){ eta[i] = X[i] * beta + Z[i] * ranef[id[i]]; }\n',
                                MODEL = 'for(i in 1:N){ y[i] ~ student_t(y_nu, eta[i], y_sigma); }\n',
                                GENERATED = 'for(i in 1:N){ y_new[i] = student_t_rng(y_nu, eta[i], y_sigma); }\n'),
                   )
    if(y_nu_free){
        model$DATA = paste(model$DATA, 'real<lower=0> y_nu_rate;\n')
        model$PARAMETERS = paste(model$PARAMETERS, 'real<lower=1> y_nu;\n')
        model$TRANSF_DECL = 'real<lower=0> y_nu_minus_one;\n '
        model$TRANSFORMED = paste(model$TRANSFORMED, 'y_nu_minus_one = y_nu - 1;\n')
        model$MODEL = paste(model$MODEL, 'y_nu_minus_one ~ exponential(y_nu_rate);\n')
    }else{
        model$DATA = paste(model$DATA, 'real<lower=1> y_nu;\n')
    }
    set(model)
}

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
#' @param return_stanfit (=F) Czy zwracać obiekt zwracany przez Stana,
#'     czy listę z ramką z próbkami i summary
#' @return Obiekt zwracany przez rstan lub lista złożona z elementów
#'     s: ramka próbek, summary: podsumowanie wyników STAN'a.
#' @export
robust_mixed = function(fixed, random, d, n = NULL,
                        y_nu_rate = NULL, y_nu = NULL, y_sigma = 1.548435,
                        beta_mu = 0, beta_sigma = NULL,
                        ranef_nu = 4,
                        chains = parallel::detectCores() - 1,
                        pars = NULL, type = 'robit',
                        auto_write = TRUE,
                        return_stanfit = F, ...){
    require(rstan)
    if(!all(pars %in% c('ranef', 'y_new'))){
        stop('Wrong parameter names. Valid values are \'ranef\' and \'y_new\'.')
    }
    if((is.null(y_nu_rate) & is.null(y_nu)) | (is.numeric(y_nu_rate) & is.numeric(y_nu)))stop('Provide either y_nu_rate or y_nu.')
    ## Sprawdzamy poprawność argumentów ze względu na typ modelu
    if(type == 'robit'){
        if(is.null(n))stop("Number of observations per data point (n) missing.")
        if(is.null(beta_sigma)){
            warning('Using defaulf SD=20 for fixed effects priors. This could be inappropriate for unstandardized predictors.')
            beta_sigma = 20
        }
    }
    if(type == 'student'){
        if(is.null(beta_sigma)){
            stop('beta_sigma (SD of fixed effects priors) not set')
        }else if(all(beta_sigma == 20)){
            warning('SD of normal prior for beta = 20 in the linear model')
        }
    }
    ## Uzupełniamy wektor parametrów do monitorowania
    essential.pars = c('beta', 'ranef_sigma', 'C')
    if(is.numeric(y_nu_rate)){
        pars = c(essential.pars, 'y_nu', pars)
    }else{
        pars = c(essential.pars, pars)
    }
    ## Optymalizacja działania STANa
    rstan_options(auto_write = auto_write)
    options(mc.cores = parallel::detectCores())
    ## Przygotowujemy dane
    id = as.numeric(as.factor(as.character(d[[as.character(random[2])]])))
    y = d[,as.character(fixed[2])]
    X = model.matrix(fixed, d)
    ## Ewentualne rozwinięcie wektora priorów
    if(length(beta_sigma) < ncol(X)){
        warning('Using the same SD (beta_sigma[1]) for all fixed effects priors')
        beta_sigma = rep(beta_sigma[1], ncol(X))
    }
    if(length(beta_mu) < ncol(X)){
        warning('Using the same mu (beta_mu[1]) for all fixed effects priors')
        beta_mu = rep(beta_mu[1], ncol(X))
    }
    Z = model.matrix(random, d)
    data = list(D = ncol(X), R = ncol(Z), N = nrow(X), I = max(id),
                X = X, Z = Z, y = y, id = id,
                ranef_nu = ranef_nu,
                beta_mu = beta_mu,
                beta_sigma = beta_sigma)
    if(is.numeric(y_nu_rate)){
        data$y_nu_rate = y_nu_rate
    }else{
        data$y_nu = y_nu
    }
    if(type == 'robit'){
        d$n = n ## Zapewnia odpowiednią długość, gdy n to skalar
        data$n = d$n
        data$y_sigma = y_sigma
    }
    fit = stan(model_code = create_model(type, is.numeric(y_nu_rate)),
               data = data,
               chains = chains, pars = pars, ...)
    if(!return_stanfit){
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
