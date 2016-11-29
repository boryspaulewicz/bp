// Model mieszany "t-logistyczny", z jednym czynnikiem losowym i
// efektami losowymi o wielowymiarowym rozk�adzie t. Wielowymiarowy
// rozk�ad t zaimplementowany wed�ug opisu z Wikipedii, wielowymiarowy
// normalny zamodelowany wed�ug artyku�u o modelach mieszanych w
// Stan-ie.

// Elementy modelu specyficzne dla wersji logistycznej s� oznaczone
// tagami //#<n>

data{
  // Wymiary macierzy efekt�w ustalonych
  int<lower=1> D;
  int<lower=1> N;
  row_vector[D] X[N];
  // Macierz efekt�w losowych
  int<lower=1> R;
  row_vector[R] Z[N];
  int<lower=1> I;
  int<lower=1, upper=I> id[N];
  int<lower=0> y[N];
  // Efekty losowe
  // real<lower=1> ranef_nu; // Liczba stopni swobody rozk�adu efekt�w losowych, domy�lnie 4
  real<lower=0> ranef_nu_rate; // parametr prioru exp dla nu efekt�w losowych
  // Prior dla efekt�w ustalonych
  real<lower=0> beta_sigma; // domy�lnie 20, dla stdandardyzowanych predyktor�w mo�na 5
  // Parametry rozk�adu t dla obserwacji
  // real<lower=1> y_nu; // domy�lnie 4
  real<lower=0> residuals_nu_rate; // parametr prioru exp dla nu reszt
  real<lower=0> y_sigma; // Domy�lnie 1.548435, �eby wsp�czynniki odpowiada�y probitowi
  // Liczba obserwacji na punkt danych
  //#1
  int<lower=0> n[N];
  //#1
}

parameters{
  vector[D] beta; // Efekty ustalone
  real<lower=1> y_nu;
  vector<lower=0>[R] ranef_sigma; // sd efekt�w losowych
  real<lower=1>[R] ranef_nu; // nu efekt�w losowych
  // Parametryzacja skorelowanych efekt�w losowych
  cholesky_factor_corr[R] L;
  vector[R] z_ranef[I];
  vector<lower=0>[I] u_ranef;
}

transformed parameters{
  vector[R] ranef[I];
  // Macierz korelacji
  matrix[R, R] C;
  // Warto�� oczekiwania zmiennej zale�nej
  vector[N] eta;
  for(i in 1:I){
    // Efekty losowe maj� rozk�ad t ze �redni� 0
    ranef[i] = sqrt(ranef_nu / u_ranef[i]) * diag_pre_multiply(ranef_sigma, L) * z_ranef[i];
  }
  C = L * L';
  for(i in 1:N){
    // Odporna regresja logistyczna
    //
    //#2
    eta[i] = student_t_cdf(X[i] * beta + Z[i] * ranef[id[i]], y_nu, 0, y_sigma);
    //#2
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
    //#3
    y[i] ~ binomial(n[i], eta[i]);
    //#3
  }
}

// Posterior predykcyjny
generated quantities {
  vector[N] y_new;
  for(i in 1:N){
    //#4
    y_new[i] = binomial_rng(n[i], eta[i]);
    //#4
  }
}

