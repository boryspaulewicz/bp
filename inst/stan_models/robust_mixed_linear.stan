// Model mieszany liniowy z resztami o rozk³adzie t, jednym czynnikiem
// losowym i efektami losowymi o wielowymiarowym rozk³adzie
// t. Wielowymiarowy rozk³ad t zaimplementowany wed³ug opisu z
// Wikipedii, wielowymiarowy normalny zamodelowany wed³ug artyku³u o
// modelach mieszanych w Stan-ie.

data{
  // Macierz efektów ustalonych
  int<lower=1> D;
  int<lower=1> N;
  row_vector[D] X[N];
  // Macierz efektów losowych
  int<lower=1> R;
  row_vector[R] Z[N];
  int<lower=0> y[N];
  int<lower=1> I;
  int<lower=1, upper=I> id[N];
  // Parametry rozk³adu t dla obserwacji
  real<lower=1> y_nu; // 4
  // Efekty losowe
  real<lower=1> ranef_nu; // 4
  // Prior dla efektów ustalonych
  real<lower=0> beta_sigma; // du¿a liczba
}

parameters{
  vector[D] beta;
  real<lower=0> y_sigma;
  vector<lower=0>[R] ranef_sigma;
  cholesky_factor_corr[R] L;
  // U¿ywane do modelowania efektów losowych
  vector[R] z_ranef[I];
  vector<lower=0>[I] u_ranef;
}

transformed parameters{
  vector[R] ranef[I];
  // Macierz korelacji
  matrix[R, R] C;
  vector[N] fit;
  for(i in 1:I){
    // Efekty losowe maj± rozk³ad t ze ¶redni± 0
    ranef[i] <-sqrt(ranef_nu / u_ranef[i]) * diag_pre_multiply(ranef_sigma, L) * z_ranef[i];
  }
  C <- L * L';
  // Liczymy w ten sposób, bo ponownie u¿ywamy tej warto¶ci do
  // próbkowania predykcji
  for(i in 1:N){
    fit[i] <- X[i] * beta + Z[i] * ranef[id[i]];
  }
}

model{
  beta ~ normal(0, beta_sigma);
  L ~ lkj_corr_cholesky(2.0); 
  for(i in 1:I){
    z_ranef[i] ~ normal(0, 1);
    u_ranef[i] ~ chi_square(ranef_nu);
  }
  for(i in 1:N){
    y[i] ~ student_t(y_nu, fit[i], y_sigma);
  }
}

generated quantities {
  vector[N] y_new;
  for(i in 1:N){
    y_new[i] <- student_t(y_nu, fit[i], y_sigma);
  }
}
