n_zero_eval = function(eval, th = 1e-10){
  # returns number of zero eigen values
  return(sum(eval < th))
}
eval_zero = function(eval, k){
  eval_res = eval;
  eval_res[1:k] = 0.0
  return(eval_res)
}

GS_complete = function(Evec, k){
  n = dim(Evec)[1]
  if (k < n){
    Q = matrix(NA, n, k)
    X = Evec[,1:k]
    X[,1] = rep(1/sqrt(n), n)
    Q[,1] = rep(1/sqrt(n), n)
    V = X 
    for (i in 2:k){
      for (j in 1:(i-1)){
        V[,i] = V[,i] -((t(Q[,j])%*%X[,i])%*%Q[,j])
      }
      Q[,i] = V[,i]/sqrt(sum(V[,i]*V[,i]))
    }
    return(cbind(Q, Evec[,(k+1):n]))
  }else {
    stop("the number of zero eigenvalue is", k, "out of", n)
  }
}

cen_eigen = function(K){
  n = dim(K)[1]
  E = eigen(K)
  l = E$values[n:1]
  Q = E$vectors[,n:1]
  k = n_zero_eval(l)
  if(k>0){
    l = eval_zero(l,k)
    if(k>1){
      Q = GS_complete(Q,k)
    } else {
      Q[,1] = rep(1/sqrt(n), n) 
    }
    return(list('values'= l,'vectors'= Q)) 
  } else{
    stop('k must be positive; found k=' , as.character(k))
  }
}

mat_vec_prod = function(A,s, nd, nnd){
  S = matrix(s, nd,nnd)
  return(c(t(A%*%S)))
}

kron_mat_vec_2d = function(A1, A2, x){
  n1 = dim(A1)[1]
  n2 = dim(A2)[1]
  X = matrix(x,n2,n1)
  return(c(A2%*%(X%*%t(A1))))
}

pcg = function(Q1,Q2,eval,d,y,initx, tol = 1e-3, itmax=1000){
  # output : x = (K + D)^{-1} y  with preconditioning matrix C = D^{-1/2} and D = diag_matrix(d)
  #          note K = Q %*% diag(eval) %*% Q' where Q = Q1 \otime Q2
  n1 = dim(Q1)[1] 
  n2 = dim(Q2)[1]
  n = n1*n2
  #dinv = (1/d) 
  #initial values
  x = initx
  v = eval*kron_mat_vec_2d(t(Q1), t(Q2), x) # diag_matrix(eval)%*%(Q)%*%x
  Kx = kron_mat_vec_2d(Q1,Q2,v)             # Q%*%v
  r = y - Kx - (d*x)                        # y- (K+D)x = y - Kx - Dx
  z = r/d                                   # D^{-1}%*%r , preconditioning
  p = z
  k = 0
  while (k<itmax){
    rz_old = c(crossprod(r,z))
    #  p'.(K+D).p = p'.(K.p + D.p)
    v = eval*kron_mat_vec_2d(t(Q1), t(Q2), p)
    KaddD.p = kron_mat_vec_2d(Q1,Q2,v) + d*p # (K+D)p = Kp + Dp
    alpha = rz_old/c(crossprod(p, KaddD.p))
    x = x + alpha*p
    r = r - alpha*(KaddD.p)
    if ( (sqrt(c(crossprod(r))) < tol) == TRUE){
      break
    } 
    z = r/d  #D^{-1}%*%r 
    beta = c(crossprod(r,z))/rz_old
    p = z + beta*p
    k = k + 1
  }
  return(list('sol'= x, 'iteration'=k))
}
