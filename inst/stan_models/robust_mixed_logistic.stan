// Model mieszany "t-logistyczny", z jednym czynnikiem losowym i
// efektami losowymi o wielowymiarowym rozk³adzie t. Wielowymiarowy
// rozk³ad t zaimplementowany wed³ug opisu z Wikipedii, wielowymiarowy
// normalny zamodelowany wed³ug artyku³u o modelach mieszanych w
// Stan-ie. (Manuskrypt jest w podkatalogu stan_models).

data{
  // Wymiary macierzy efektów ustalonych
  int<lower=1> D;
  int<lower=1> N;
  row_vector[D] X[N];
  // Macierz efektów losowych
  int<lower=1> R;
  row_vector[R] Z[N];
  int<lower=1> I;
  int<lower=1, upper=I> id[N];
  //#
  int<lower=0> y[N];
  //@
  //
  // Efekty losowe
  real<lower=1> ranef_nu; // Liczba stopni swobody rozk³adu efektów losowych, domy¶lnie 4
  // real<lower=0> ranef_nu_rate; // parametr prioru exp dla nu efektów losowych
  // Prior dla efektów ustalonych
  real beta_mu[D];
  real<lower=0> beta_sigma[D]; // domy¶lnie 20 dla robit, dla stdandardyzowanych predyktorów mo¿na 5
  // Parametry rozk³adu t dla obserwacji
  real<lower=1> y_nu; // domy¶lnie 4
  real<lower=0> y_sigma; // Domy¶lnie 1.548435, ¿eby wspó³czynniki odpowiada³y probitowi
  // Liczba obserwacji na punkt danych
  //#
  int<lower=0> n[N];
  //@
}

parameters{
  vector[D] beta; // Efekty ustalone
  vector<lower=0>[R] ranef_sigma; // sd efektów losowych
  // Parametryzacja skorelowanych efektów losowych
  cholesky_factor_corr[R] L;
  vector[R] z_ranef[I];
  vector<lower=0>[I] u_ranef;
}

transformed parameters{
  vector[R] ranef[I];
  // Macierz korelacji
  matrix[R, R] C;
  //#
  //
  // Warto¶æ oczekiwania zmiennej zale¿nej
  vector[N] eta;
  //@
  for(i in 1:I){
    // Efekty losowe maj± rozk³ad t ze ¶redni± 0
    ranef[i] = sqrt(ranef_nu / u_ranef[i]) * diag_pre_multiply(ranef_sigma, L) * z_ranef[i];
  }
  C = L * L';
  for(i in 1:N){
    //#
    eta[i] = student_t_cdf(X[i] * beta + Z[i] * ranef[id[i]], y_nu, 0, y_sigma);
    //@
  }
}

model{
  for(i in 1:D){
    beta ~ normal(beta_mu[i], beta_sigma[i]);
  }
  L ~ lkj_corr_cholesky(2.0); 
  for(i in 1:I){
    z_ranef[i] ~ normal(0, 1);
    u_ranef[i] ~ chi_square(ranef_nu);
  }
  for(i in 1:N){
    //#
    y[i] ~ binomial(n[i], eta[i]);
    //@
  }
}

// Posterior predykcyjny
generated quantities {
  vector[N] y_new;
  for(i in 1:N){
    //#
    y_new[i] = binomial_rng(n[i], eta[i]);
    //@
  }
}

