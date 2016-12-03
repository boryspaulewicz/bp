## Zmienia wartość fragmentów modelu na podstawie listy POLE =
## 'wartość'
set = function(fields){
    res = paste(readLines(paste(path.package('bp'), '/stan_models/robust_mixed.stan', sep = ''), encoding = 'utf8'),
                collapse = '\n')
    for(f in names(fields))res = gsub(sprintf('// %s', f), fields[[f]], res)
    res
}

## Tworzy kod modelu na podstawie szablonu
create_model = function(type, y_nu_rate, ranef_nu_rate){
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

    if(is.numeric(y_nu_rate)){
        model$DATA = paste(model$DATA, 'real<lower=0> y_nu_rate;\n')
        model$PARAMETERS = paste(model$PARAMETERS, 'real<lower=1> y_nu;\n')
        model$TRANSF_DECL = paste(model$TRANSF_DECL, 'real<lower=0> y_nu_minus_one;\n ')
        model$TRANSFORMED = paste(model$TRANSFORMED, 'y_nu_minus_one = y_nu - 1;\n')
        model$MODEL = paste(model$MODEL, 'y_nu_minus_one ~ exponential(y_nu_rate);\n')
    }else{
        model$DATA = paste(model$DATA, 'real<lower=1> y_nu;\n')
    }

    if(is.numeric(ranef_nu_rate)){
        if(length(ranef_nu_rate) > 1){
            model$DATA = paste(model$DATA, 'real<lower=0> ranef_nu_rate[R];\n')
            model$PARAMETERS = paste(model$PARAMETERS, 'vector<lower=1>[R] ranef_nu;\n')
            model$TRANSF_DECL = paste(model$TRANSF_DECL, 'real<lower=0> ranef_nu_minus_one[R];\n ')
            model$TRANSFORMED = paste(model$TRANSFORMED, 'for(i in 1:R){ ranef_nu_minus_one[i] = ranef_nu[i] - 1; }\n')
            model$MODEL = paste(model$MODEL, 'for(i in 1:R){ ranef_nu_minus_one[i] ~ exponential(ranef_nu_rate[i]); }\n')
        }else{
            model$DATA = paste(model$DATA, 'real<lower=0> ranef_nu_rate;\n')
            model$PARAMETERS = paste(model$PARAMETERS, 'real<lower=1> ranef_nu;\n')
            model$TRANSF_DECL = paste(model$TRANSF_DECL, 'real<lower=0> ranef_nu_minus_one;\n ')
            model$TRANSFORMED = paste(model$TRANSFORMED, 'ranef_nu_minus_one = ranef_nu - 1;\n')
            model$MODEL = paste(model$MODEL, 'ranef_nu_minus_one ~ exponential(ranef_nu_rate);\n')
        }
    }else{
        model$DATA = paste(model$DATA, 'real<lower=1> ranef_nu;\n')
    }

    set(model)
}

#' Dopasowanie uogólnionego odpornego mieszanego modelu liniowego
#'
#' Implementacja modelu jest oparta na publikacji Sorensen'a,
#' Hohenstain'a i Vasishth'a (2015): Bayesian linear mixed models
#' using Stan: A tutorial for psychologists, linguists, and cognitive
#' scientists.
#' @param fixed Formuła modelu dla efektów ustalonych.
#' @param random Formuła modelu dla efektów losowych, postaci g ~ x1 +
#'     x2 ...
#' @param d Zbiór danych.
#' @param n Liczba obserwacji na punkt danych (robit)
#' @param y_nu_rate Parametr prioru dla liczby stopni swobody rozkładu
#'     zmiennej zależnej. Jeżeli określony, y_nu będzie zignorowany, a
#'     nu rozkładu zmiennej zależnej będzie wolnym parametrem.
#' @param y_nu Ustalona liczba stopni swobody rozkładu zmiennej
#'     zależnej. Nie można ustawiać jednocześnie z y_nu_rate.
#' @param y_sigma (Uwzględniane tylko dla modelu logistycznego)
#'     Odchylenie standardowe rozkładu t modelującego zmienną
#'     zależną. To jest wartość arbitralna - domyślnie przyjęto
#'     wartość dającą współczynniki porównywalne do zwykłej regresji
#'     logistycznej.
#' @param beta_mu Skalar lub wektor (długości równej liczbie efektów
#'     ustalonych) określający średnią/e prioru normalnego dla efektów
#'     ustalonych.
#' @param beta_sigma Skalar lub wektor (długości równej liczbie
#'     efektów ustalonych) określający odchylenie standardowe prioru
#'     dla efektów ustalonych. W robit domyślnie 20, co daje słabo
#'     informacyjny prior przy założeniu, że predyktory są binarne lub
#'     mają odchylenie standardowe bliskie 1.
#' @param ranef_nu_rate Parametr prioru nu dla efektów
#'     losowych. Skalar albo wektor o długości równej liczbie efektów
#'     losowych. Jeżeli określony, to nu dla efektów losowych będzie
#'     wolnym parametrem, a ranef_nu będzie zignorowany.
#' @param ranef_nu Liczba stopni swobody (nu) rozkładu efektów
#'     losowych. Nie można podawać jednocześnie z ranef_nu_rate.
#' @param chains Liczba równoległych próbników. Domyślnie liczba
#'     rdzeni - 1.
#' @param pars Wektor tekstowy nazw parametrów, dla których chcemy
#'     próbki, poza automatycznie uwzględnianymi beta, ranef_sigma, C
#'     i y_nu, jeżeli y_nu_rate jest ustalone. Dopuszczalne wartości
#'     to 'y_new' (posterior predykcyjny) i 'ranef' (efekty losowe).
#' @param type [student/robit] Określa rozkład zmiennej zależnej.  Typ
#'     robit to rozkład dwumianowy z funkcją łączącą równą
#'     dystrybuancie rozkładu t. Typ student to rozkład t.
#' @param auto_write Parametr przekazywany do funkcji stan.
#' @param return_stanfit (=F) Czy zwracać obiekt zwracany przez Stana,
#'     czy listę z ramką z próbkami i summary
#' @return Obiekt zwracany przez rstan (gdy return_stanfit = T) lub
#'     lista złożona z elementów samples: ramka próbek o nazwach
#'     odpowiadających nazwom efektów w modelu, summary: podsumowanie
#'     wyników próbkowania zwracane przez STAN'a, ranef: wektor nazw
#'     efektów losowych, fixef: wektor nazw efektów ustalonych, model:
#'     kod dopasowanego modelu.
#' @export
robust_mixed = function(fixed, random, d, n = NULL,
                        y_nu_rate = NULL, y_nu = 4,
                        y_sigma = 1.548435,
                        beta_mu = 0, beta_sigma = NULL,
                        ranef_nu_rate = NULL, ranef_nu = 4,
                        chains = parallel::detectCores() - 1,
                        pars = 'ranef', type = 'robit',
                        auto_write = TRUE,
                        return_stanfit = F, ...){

    ## Sprawdzamy poprawność wektora parametrów monitorowanych
    if(!all(pars %in% c('ranef', 'y_new'))){
        stop('Wrong parameter names. Valid values are \'ranef\' and \'y_new\'.')
    }

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
    if(is.numeric(y_nu_rate))essential.pars = c(essential.pars, 'y_nu')
    if(is.numeric(ranef_nu_rate))essential.pars = c(essential.pars, 'ranef_nu')
    pars = unique(c(essential.pars, pars))

    ## Optymalizacja działania STANa
    require(rstan)
    rstan_options(auto_write = auto_write)
    options(mc.cores = parallel::detectCores())

    ## Przygotowujemy dane
    y = d[, as.character(fixed[2])]
    X = model.matrix(fixed, d)
    Z = model.matrix(random, d)
    ## Potrzebujemy całkowitoliczbowej reprezentacji id, ale później
    ## przydadzą nam się oryginalne nazwy grup/osób
    idfactor = as.factor(as.character(d[[as.character(random[2])]]))
    id = as.numeric(idfactor)
    idnames = levels(idfactor)

    ## Ewentualne rozwinięcie wektorów priorów
    if(length(beta_sigma) < ncol(X)){
        warning('Using the same SD (beta_sigma[1]) for all fixed effects priors')
        beta_sigma = rep(beta_sigma[1], ncol(X))
    }
    if(length(beta_mu) < ncol(X)){
        warning('Using the same mu (beta_mu[1]) for all fixed effects priors')
        beta_mu = rep(beta_mu[1], ncol(X))
    }
    if(length(ranef_nu_rate) > 1)stop('More than one free nu parameter for random effects not implemented yet.')
    ## if(!any(length(ranef_nu_rate) == c(0, 1, ncol(Z)))){
    ##     stop('ranef_nu_rate is longer than 1 but shorter than the number of random effects')
    ## }
    
    data = list(D = ncol(X), R = ncol(Z), N = nrow(X), I = max(id),
                X = X, Z = Z, y = y, id = id,
                beta_mu = beta_mu,
                beta_sigma = beta_sigma)

    ## nu wolne czy ustalone
    if(is.numeric(y_nu_rate)){
        data$y_nu_rate = y_nu_rate
    }else{
        data$y_nu = y_nu
    }
    if(is.numeric(ranef_nu_rate)){
        data$ranef_nu_rate = ranef_nu_rate
    }else{
        data$ranef_nu = ranef_nu
    }
    
    if(type == 'robit'){
        d$n = n ## Zapewnia odpowiednią długość, gdy n to skalar
        data$n = d$n
        data$y_sigma = y_sigma
    }

    model = create_model(type, y_nu_rate, ranef_nu_rate)

    fit = stan(model_code = model,
               data = data,
               chains = chains, pars = pars, ...)

    if(!return_stanfit){## Zwracamy próbki z czytelnymi nazwami parametrów
        s = as.data.frame(extract(fit))
        names(s)[rmatch('beta', names(s))] = colnames(X)
        names(s)[rmatch('ranef_sigma', names(s))] = paste(colnames(Z), 'SD')
        ## W nazwach efektów losowych w s najszybciej zmienia się id,
        ranef_names = apply(expand.grid(idnames, colnames(Z)), 1, function(x)paste(x, collapse = '.'))
        if('ranef' %in% pars)names(s)[rmatch('ranef\\.[0-9]+\\.[0-9]+', names(s))] = ranef_names
        stan.sum = rstan::summary(fit)$summary
        ## Usuwamy z summary predykcje posterioryczne
        stan.sum = stan.sum[!rmatch('y_new', rownames(stan.sum)),]
        ## Poprawiamy nazwy parametrów w summary
        rownames(stan.sum)[rmatch('beta', rownames(stan.sum))] = colnames(X)
        rownames(stan.sum)[rmatch('ranef_sigma', rownames(stan.sum))] = paste(colnames(Z), 'SD')
        ## w nazwach efektów losowych w summary STANa najszybciej
        ## zmienia się efekt.
        if('ranef' %in% pars)rownames(stan.sum)[rmatch('ranef\\[', rownames(stan.sum))] =
                                apply(expand.grid(colnames(Z), idnames), 1, function(x)paste(rev(x), collapse = '.'))
        list(samples = s, summary = as.data.frame(stan.sum), ranef = ranef_names, fixef = colnames(X),
             model = model)
    }else{
        fit
    }
}
