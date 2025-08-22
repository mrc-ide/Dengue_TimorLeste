functions {
// Custom max function in RStan to ensure no negative values
real my_min(real k) {
  
    if (k < 0) {
      return(k);
    } else {
      return(0);
  }
}
}

data {
  
  int<lower=0> EA; //number of EA
  int AgeG; //number of age groups
  array[AgeG] int age; //age
  array[EA, AgeG] int N1; 
  array[EA, AgeG] int NPos1;
  
  array[EA, AgeG] int N0; 
  array[EA, AgeG] int NPos0;
   
}


parameters {
  real<lower=0> b_H; //clustering effect
  array[EA] real<lower=0.0001, upper=1> lambda; //EA specific FOI, i.e. force of infection
}

transformed parameters {
  matrix[EA,AgeG] seroprev1; //seroprevalence for positive
  matrix[EA,AgeG] k; //variable to store lambda age bH
  matrix[EA,AgeG] seroprev0; //seroprevalence for negative
  
  for (i in 1:EA) {for (a in 1:AgeG) {
    k[i,a] = -lambda[i] * (age[a] + b_H) ; 
    seroprev1[i,a] = 1 - exp(my_min(k[i,a])); //Eq 2, has bH
  }}
  
  for (i in 1:EA)  for (a in 1:AgeG) {
    seroprev0[i,a] = 1 - exp(-lambda[i] * age[a]); //Eq 1, NO bH
  }
}


model {
  //Priors
  lambda ~ normal(0,1); //prior on lambda
  b_H ~ normal(1,0.5);  //prior on bH
  
  //Likelihood
  for (i in 1:EA) target += binomial_lpmf(NPos1[i,] | N1[i,], seroprev1[i,]);
  for (i in 1:EA) target += binomial_lpmf(NPos0[i,] | N0[i,], seroprev0[i,]);
  
}

