functions {
  int n_zero_eval(vector eval, int n){
    // count number of zero eigen values (smaller than threshold)
    // Input: eval - a vector of eigenvalues
    //        n ---- a size of a matrix
    // Output: number of zero eigen values
    real d = eval[1];
    int i = 1 ;
    while (d < 1e-7){
        i += 1 ;
        d = eval[i];
    }
    return i - 1 ;
  }

  vector eval_zero(vector eval, int k, int n){
    // replace eigenvalues with zero
    vector[n] evalz = eval;
    for (i in 1:k)  evalz[i] = 0.0;
    return evalz ;
  }

  matrix GS_complete(matrix Evec, int k, int n){
    // Gram-Schmidt process
    // Input: Evec - original sets of eigenvaectors
    //        k, n --the number of zero eigen values, a size of a matrix
    // Output: sets of eigenvaectors after Gram-Schumidt process
    matrix[n,k] Q;
    matrix[n,k] V ;
    matrix[n,k] X = Evec[,1:k];
    X[,1] =  rep_vector(1/sqrt(n),n) ;
    V = X ;
    Q[,1] =  V[,1];
    for (i in 2:k){
      for (j in 1:(i-1)){
        V[,i] = V[,i] - ((Q[,j]')*X[,i])*Q[,j];
      }
      Q[,i] = V[,i]/sqrt(sum(V[,i].*V[,i]));
    }
    return append_col(Q,Evec[,(k+1):n]);
  }

  matrix cen_eigen_decompose(matrix K, int N){
    // Output: Eigenvalues and eigenvector of a centered Gram matrix
    // Input: K - Centred Gram matrix, N - nrow = ncol of K
    matrix[N, N+1] R;
    {
      matrix[N,N] Q ;
      vector[N] l = eigenvalues_sym(K);
      {
        int k = n_zero_eval(l,N);
        if (k > 0){
          Q = eigenvectors_sym(K);
          l = eval_zero(l, k, N);
          if (k>1){
              Q = GS_complete(Q, k, N);
          } else {
            Q[,1] = rep_vector(1/sqrt(N),N);
          }
        } else {
          reject("k must be positive; found k=", k);
        }
      }
      R[1:N,1:N] = Q ;
      R[1:N, N+1] = l;
    }
    return R ;
  }
  /* Kernel related functions */
  // Gram matrix for different kernels
  real kernel_SE(vector x, vector y, real sq_rho){
    // Output k(x,y) = cov(f(x), f(y)) where k is s.e. kernel
    // Input: x, y, and the length scale (rho) of s.e. kernel
    return exp(-squared_distance(x,y)/(2*sq_rho));
  }
  matrix Gram_SE(matrix X, int N, real rho){
    // Output: n by n Gram matrix with squared exponential kernel
    // Input: X - predictor
    //        n - n_rows of X
    //        rho - length scale of s.e. kernel
    matrix[N,N] K = diag_matrix(rep_vector(1,N));
    real sq_rho = square(rho);
    for (i in 1:(N-1)){
      for (j in (i+1):N){
        K[i,j] = kernel_SE(to_vector(X[i,]),to_vector(X[j,]),sq_rho);
        K[j,i] = K[i,j];
      }
    }
    return K ;
  }

  matrix Gram_fBM(matrix X, int N, real Hurst){
    // Output: Gram matrix with fractional brownian motion
    // Input: X - predictor
    //        N - nrow of X
    //        Hurst - Hurst coefficeint of fBM kernel
    matrix[N, N] K;
    matrix[N, N] E ;
    matrix[N, N] B = rep_matrix(0, N, N);
    vector[N] d;
    matrix[N, N] A = diag_matrix(rep_vector(1,N)) - (1.0/N)*rep_matrix(1, N, N);
    matrix[N, N] Xcp = X * X' ;
    vector[N] dvec = diagonal(Xcp);
    for (i in 1:(N-1)){
      d[i] = pow(fabs(dvec[i]), Hurst);
      for (j in (i+1):N){
        B[i,j] = pow(fabs(dvec[i] + dvec[j] - 2 * Xcp[i,j]), Hurst);
        B[j,i] = B[i,j];
      }
    }
    d[N] = pow(fabs(dvec[N]), Hurst);
    E = rep_matrix(d, N);
    K = 0.5 * (E + E' - B);
    return K ;
  }

  matrix Gram_fBM_sq_cen(matrix X, int N, real Hurst){
    // Output: Gram matrix with square centered fractional brownian motion
    // Input: X - predictor
    //        N - nrow of X
    //        Hurst - Hurst coefficeint of fBM kernel
    matrix[N, N] K = Gram_fBM(X, N, Hurst);
    matrix[N, N] A = diag_matrix(rep_vector(1,N)) - (1.0/N)*rep_matrix(1, N, N);
    K = A * K * A;
    K =  K * K;
    // to ensure that it's symmetric
    K = 0.5 *(K + K');
    return K;
  }

  matrix Gram_centring(matrix K, int N){
      // Output: centered Gram matrix
      // Input: K - Uncentered Gram matrix
      //        N - nrow of K
      matrix[N, N] A = diag_matrix(rep_vector(1,N)) - (1.0/N)*rep_matrix(1, N, N);
      return A * K * A ;
  }

  matrix Gram_square(matrix K, int N){
      // Output: squared Gram matrix
      // Input: K -  Gram matrix
      //        N - nrow of K
      matrix[N,N] Ksq = K * K ;
      Ksq = 0.5 *(Ksq + Ksq');
      return Ksq;
  }

  vector kvec_SE(matrix X, vector x_tes, int N, real rho){
    // Output: a vector of k(x*, x_1),...,k(x*, x_N) for a given test point x*
    //          with squared exponential kernel
    // Input: X - (traing) predictor
    //        x_tes - test point
    //        N - n_rows of X
    //        rho - length scale of s.e. kernel
    vector[N] kvec ;
    real sq_rho = square(rho);
    for (i in 1:N){
      kvec[i] = kernel_SE(x_tes, to_vector(X[i,]), sq_rho);
    }
    return kvec ;
  }

  vector kvec_fBM(matrix X, vector x_new, int N, real Hurst){
    // Output: a vector of k(x*, x_1),...,k(x*, x_N) for a given test point x*
    //          with fBM kernel
    // Input: X - (traing) predictor
    //        x_tes - test point
    //        N - n_rows of X
    //        Hurst - hurst coefficeint of fBM kernel
    vector[N] kvec ;
    real t1 = pow(dot_self(x_new), Hurst) ;
    vector[N] t2 ;
    vector[N] t12 ;
    for (i in 1:N){
      t12[i] = pow(squared_distance(x_new,X[i,]), Hurst);
      t2[i] = pow(dot_self(X[i,]), Hurst);
    }
    kvec = 0.5*(rep_vector(t1,N) + t2 - t12) ;
    return kvec ;
  }

  vector kvec_cen(vector kvec, vector Krowsum, int N){
    return kvec - rep_vector(sum(kvec)/N, N) -(1.0/N)*Krowsum + rep_vector(sum(Krowsum)/square(N), N);
  }

  vector kvec_sq(vector kvec, matrix K){
    return K * kvec ;
  }

  real kstar_cen(real kstar,vector kvec, vector Krowsum,int N){
    return kstar - 2*sum(kvec)/N + sum(Krowsum)/square(N) ;
  }

  vector mat_vec_prod(matrix A, vector s, int Nd, int NNd ){
    //matrix-vector product involving Kronecker product
    matrix[Nd, NNd] S = to_matrix(s, Nd, NNd);
    matrix[NNd, Nd] Z = (A * S)';
    return to_vector(Z);
  }
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
}
transformed data {
  int N = N1*N2*N3;
  vector[N1] l1;
  vector[N2] l2;
  vector[N3] l3;
  vector[N] m ;
  {
    matrix[N1, N1] Q1 ;
    matrix[N2, N2] Q2 ;
    matrix[N3, N3] Q3 ;

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
  }
}
parameters {
  real<lower=0> alpha0 ;
  real<lower=0> alpha1 ;
  real<lower=0> alpha2 ;
  real<lower=0> alpha3 ;
  real<lower=0> sigma ;
}
model{
  vector[N] eval ;
  // computing eigenvalues of the model matrix K + sigma^2*I

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
      //vector[N] t13 = to_vector(to_vector(e3*d2')*e1');
      //vector[N] t23 = to_vector(to_vector(e3*e2')*d1');
      eval = square(alpha0)*(t0 + t1 + t2 + t3 + t12) + square(sigma)*rep_vector(1,N);
     }
  //prior
  target += std_normal_lpdf(alpha0);
  target += std_normal_lpdf(alpha1);
  target += std_normal_lpdf(alpha2);
  target += std_normal_lpdf(alpha3);
  target += std_normal_lpdf(sigma);
  //likelihood
  target += -0.5 * sum(square(m)./eval) - 0.5 * sum(log(eval));
}
