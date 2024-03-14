eucledian_dist = function(x,y){
  return(sqrt(sum((x-y)^2)))
}
fBM_kernel = function(x,y, Hurst=0.5){
  r2 = sum((x-y)^2)
  return(0.5*( (sum(x^2))^Hurst + (sum(y^2))^Hurst - r2^Hurst ))
}
fBM_gram = function(X, Hurst){
  n = dim(X)[1]
  #K = E = matrix(NA,n,n)
  B = matrix(0,n,n)
  A = diag(n) - (1/n)*matrix(1,n,n)
  Xcp = tcrossprod(X)
  dvec = diag(Xcp)
  d = rep(NA, n)
  for (i in 1:(n-1)){
    d[i] = abs(dvec[i])^Hurst
    for (j in (i+1):n){
      B[i,j] = abs(dvec[i]+dvec[j] - 2*Xcp[i,j])^Hurst
      B[j,i] = B[i,j]
    }
  }
  d[n] = abs(dvec[n])^Hurst
  E = matrix(rep(d,times=n),n,n)
  K = 0.5 * (E+t(E) - B)
}

fBM_kvec = function(X, x_new, Hurst=0.5){
  n = dim(X)[1]
  kvec = rep(NA, n)
  for (i in 1:n){
    kvec[i] = fBM_kernel(X[i,], x_new, Hurst)
  }
  return(kvec)
}

sq_cen_fBM_kernel_mat = function(X, Hurst, pred=FALSE){
  n = dim(X)[1]
  #K = E = matrix(NA,n,n)
  B = matrix(0,n,n)
  A = diag(n) - (1/n)*matrix(1,n,n)
  Xcp = tcrossprod(X)
  dvec = diag(Xcp)
  d = rep(NA, n)
  for (i in 1:(n-1)){
    d[i] = abs(dvec[i])^Hurst
    for (j in (i+1):n){
      B[i,j] = abs(dvec[i]+dvec[j] - 2*Xcp[i,j])^Hurst
      B[j,i] = B[i,j]
    }
  }
  d[n] = abs(dvec[n])^Hurst
  E = matrix(rep(d,times=n),n,n)
  K = 0.5 * (E+t(E) - B)
  K = A%*%K%*%A
  if (pred == FALSE){
    return(K%*%K)
  } else {
    return(list('Ksq'=K%*%K, 'K'=K, 'Browsum'=rowSums(B)))
  }
}

sq_cen_kernel_vec = function(K, Browsum, X, x_new, Hurst){
  n = length(Browsum)
  r2vec = rep(NA, n)
  Bsum = sum(Browsum)
  for (i in 1:n){
    r = eucledian_dist(x_new, X[i,])
    r2vec[i] = r^(2*Hurst)
  }
  k = 0.5* (-r2vec + (sum(r2vec)/n)*rep(1,n) + (1.0/n)*Browsum - (Bsum/n^2)*rep(1,n))
  kvar =crossprod(k)
  k = (K%*%k)
  return(list('kvec'=k, 'kvar'=kvar))
}

se_kernel = function(x,y,sq_rho){
  r2 = sum((x-y)^2)
  return(exp(-r2/(2*sq_rho)))
}

se_gram = function(X,sq_rho){
  n = dim(X)[1]
  K = diag(n)
  for (i in 1:(n-1)){
    for (j in (i+1):n){
     K[i,j] =  se_kernel(X[i,],X[j,],sq_rho)
     K[j,i] = K[i,j]
    }
  }
  return(K)
}

se_kvec = function(X, xnew, sq_rho){
  n = dim(X)[1]
  kvec = rep(NA, n)
  for (i in 1:n){
    kvec[i] = se_kernel(X[i,], x_new, sq_rho)
  }
  return(kvec)
}


cen_gram = function(K){
  n = dim(K)[1]
  A = diag(n) - matrix((1/n),n,n)
  return(A%*%K%*%A)
}
sq_gram = function(K){
  return(K%*%K)
}

cen_kvec =  function(kvec, K){
  n = length(kvec)
  Krsum = as.vector(K%*%rep(1,n))
  Ksum = sum(Krsum)
  return(kvec - rep(sum(kvec)/n, n) - ((1.0/n)*Krsum) + rep(Ksum/(n^2), n) )
}

sq_kvec = function(kvec, K){
  return(K%*%kvec)
}

cen_kstar = function(k_star, kvec, K){
  n = length(kvec)
  Ksum = sum(K)
  return(k_star-2*sum(kvec)/n + Ksum/(n^2))
}


