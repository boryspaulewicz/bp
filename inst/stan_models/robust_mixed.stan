// Szablon dla mieszanego modelu odpornego. W miejsca oznaczone
// komentarzami typu // <SECTION> wstawiane s± fragmenty modelu
// odpowiednio do typu (robit, student) i dodatkowych argumentów
// (np. y_nu_rate).
//
// Model jest oparty na tym opisanym w publikacji Bayesian linear
// mixed models using Stan (Sorensen et al. 2015). Rozk³ad t dla
// efektów losowych zaimplementowany jako

data{
  // Macierz efektów ustalonych
  int<lower=1> D;
  int<lower=1> N;
  row_vector[D] X[N];
  // Macierz efektów losowych, R to liczba efektów, I to liczba
  // poziomów czynnika grupuj±cego
  int<lower=1> R;
  row_vector[R] Z[N];
  int<lower=1> I;
  int<lower=1, upper=I> id[N];
  // Efekty losowe
  real<lower=1> ranef_nu; // domy¶lnie 4
  // Prior dla efektów ustalonych
  real beta_mu[D];
  real<lower=0> beta_sigma[D];
  // DATA
}

parameters{
  vector[D] beta;
  vector<lower=0>[R] ranef_sigma;
  cholesky_factor_corr[R] L;
  // U¿ywane do modelowania skorelowanych efektów losowych
  vector[R] z_ranef[I];
  vector<lower=0>[I] u_ranef;
  // PARAMETERS
}
transformed parameters{
  // TRANSF_DECL
  vector[R] ranef[I];
  // Macierz korelacji
  matrix[R, R] C;
  vector[N] eta;
  for(i in 1:I){
    // Efekty losowe maj± rozk³ad t ze ¶redni± 0, co uzyskujemy mno¿±c
    // zmienn± losow± normaln± przez sqrt(nu / v ~ chiqs(nu)).
    ranef[i] = sqrt(ranef_nu / u_ranef[i]) * diag_pre_multiply(ranef_sigma, L) * z_ranef[i];
  }
  C = L * L';
  // TRANSFORMED
}

model{
  for(i in 1:D){
    beta[i] ~ normal(beta_mu[i], beta_sigma[i]);
  }
  L ~ lkj_corr_cholesky(2.0); 
  for(i in 1:I){
    z_ranef[i] ~ normal(0, 1);
    u_ranef[i] ~ chi_square(ranef_nu);
  }
  // MODEL
}

generated quantities {
  vector[N] y_new;
  // GENERATED
}
