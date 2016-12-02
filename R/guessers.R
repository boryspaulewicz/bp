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
                            n.chains = n.chains, sample = sample, method = 'parallel',
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
