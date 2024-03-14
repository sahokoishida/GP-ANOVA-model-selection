functions {
 #include GP_helper.stan
}
data {
  int<lower=1> N1 ;
  int<lower=1> N2 ;
  int<lower=1> N3 ;
  matrix[N1, 1] X1 ;
  matrix[N2, 1] X2 ;
  matrix[N3, 1] X3 ;
  vector[N1*N2*N3] y ;
  real<lower=0, upper=1> Hurst1;
  real<lower=0, upper=1> Hurst2;
  real<lower=0, upper=1> Hurst3;

  real<lower=0> alpha0 ;
  real<lower=0> alpha1 ;
  real<lower=0> alpha2 ;
  real<lower=0> alpha3 ;
  real<lower=0> sigma ;
}
transformed data {
  int N = N1*N2*N3;
  real mloglik;
  {
    vector[N] m ;
    vector[N] eval ;
    //vector[N] m_div_eval ;
    matrix[N1, N1] Q1 ;
    matrix[N2, N2] Q2 ;
    matrix[N3, N3] Q3 ;
    vector[N1] l1;
    vector[N2] l2;
    vector[N3] l3;
    {
      matrix[N1,N1] K1 = Gram_fBM_sq_cen(X1, N1, Hurst1);
      {
        matrix[N1, N1+1] R = cen_eigen_decompose(K1, N1);
        Q1 = R[1:N1, 1:N1] ;
        l1 = to_vector(R[1:N1, N1+1]);
      }
    }
    {
      matrix[N2, N2] K2 = Gram_fBM_sq_cen(X2, N2, Hurst2);
      {
        matrix[N2, N2+1] R = cen_eigen_decompose(K2, N2);
        Q2 = R[1:N2, 1:N2] ;
        l2 = to_vector(R[1:N2, N2+1]);
      }
    }
    {
      matrix[N3,N3] K3 = Gram_fBM_sq_cen(X3, N3, Hurst3);
      {
        matrix[N3, N3+1] R = cen_eigen_decompose(K3, N3);
        Q3 = R[1:N3, 1:N3] ;
        l3 = to_vector(R[1:N3, N3+1]);
      }
    }
    // computing m = (\otimes_{l=1}^6 Ql') * y
    {
      vector[N] s = y;
      s = mat_vec_prod(Q3', s, N3, N1*N2);
      s = mat_vec_prod(Q2', s, N2, N1*N3);
      s = mat_vec_prod(Q1', s, N1, N2*N3);
      m = s;
    }
    {
      vector[N1] e1 = square(alpha1) * l1;
      vector[N2] e2 = square(alpha2) * l2;
      vector[N1] e3 = square(alpha3) * l3;
      vector[N1] d1 = rep_vector(0,N1);
      vector[N2] d2 = rep_vector(0,N2);
      vector[N1] d3 = rep_vector(0,N3);
      d1[1] = N1;
      d2[1] = N2;
      d3[1] = N3;
      {
        vector[N] t0 = to_vector(to_vector(d3*d2')*d1');
        vector[N] t1 = to_vector(to_vector(d3*d2')*e1');
        vector[N] t2 = to_vector(to_vector(d3*e2')*d1');
        vector[N] t3 = to_vector(to_vector(e3*d2')*d1');
        vector[N] t12 = to_vector(to_vector(d3*e2')*e1');
        vector[N] t13 = to_vector(to_vector(e3*d2')*e1');
        //vector[N] t23 = to_vector(to_vector(e3*e2')*d1');
        eval = square(alpha0)*(t0 + t1 + t2 + t3 + t12 + t13 ) + square(sigma)*rep_vector(1,N);
       }
    }
    mloglik = -0.5 * sum(square(m)./eval) - 0.5 * sum(log(eval));
  }
}
model{
}
generated quantities {
  // marginal likelihood
  real mllik = mloglik;
}
