ArithmeticMean = function(X)                
{                                           
  N = nrow(X)                               
  return(c(sum(X[,1])/N,sum(X[,2])/N))      
}

does.intersect = function(X,m)
{
  if(X[2] <= m[2] & X[2] >= m[1])
  {
    return(T)
  }
  else if(X[1] <= m[2] & X[1] >= m[1])
  {
    return(T)
  }
  else
  {
    return(F)
  }
}

small = function(i,X)
{
  return(c(X[i],X[i + 1]))
}

sigma2k = function(X,x)
{
  Xl = X[,1]
  Xh = X[,2]
  xk = x[1]
  xk1 = x[2]
  Skl = Xl[Xl >= xk1]
  Skh = Xh[Xh <= xk]
  Nk = length(Skl) + length(Skh)
  Sk = (sum(Skl) + sum(Skh))
  
  if(Nk==0)
  {
   sigma = 0
  }
  else
  {

    rk = Sk/Nk
    
    if(rk <= xk1 & rk >= xk)
    {
      
      sigma = (1/nrow(X))*(sum((Skl - rk)^2) + sum((Skh - rk)^2))
    }
    else
    {
      sigma = NA
    }
    
  }
  return(sigma)
}

small = function(i,X)
{
  return(c(X[i],X[i + 1]))
}

l_var = function(X)
{
  A = as.vector(X)
  ordered_A = A[order(A)]
  mean = ArithmeticMean(X)

  smallints = lapply(1:(length(A)-1),function(i,X = ordered_A) small(i,X))
  ints_intersect = smallints[unlist(lapply(smallints, function(x) does.intersect(x,mean)))]
  sigmas = lapply(ints_intersect,function(x) sigma2k(X = X,x))
  sigmas = unlist(sigmas)
  sigmas = na.omit(sigmas)
  if(length(sigmas)==0)
  {
    return(0)
  }
  else
  {
    return(min(sigmas))
  }
}



definitely_different = function(X)
{
  Xmid = (X[,2] + X[,1])/2
  
  Xdiff = (X[,2] - X[,1])/2
  n = nrow(X)
  Xn = cbind(Xmid - Xdiff/n, Xmid + Xdiff/n)
  
  A = as.vector(Xn)
  A = A[order(A)]
  B = as.vector(t(Xn))
  

  
  return(all(A==B))
  
}

u_var = function(X)
{
  
    Xmid = (X[,2] + X[,1])/2
    
    Xdiff = (X[,2] - X[,1])/2
    n = nrow(X)
    Xn = cbind(Xmid - Xdiff/n, Xmid + Xdiff/n)
    
    A = as.vector(Xn)
    Xr = A[order(A)]
    
    mean = ArithmeticMean(X)
    smallints = lapply(1:(length(A)-1),function(i,X = Xr) small(i,X))
    I = unlist(lapply(smallints, function(x) does.intersect(x,mean)))
    
    ints_intersect = smallints[I]
    ints_left = smallints[-I]
    
    nl = length(ints_left)
    
    if(nl == 0)
    {
      return(0)
    }
    else
    {
      xi_list = list()
      for(i in 1:nl)
      {
        x = ints_left[[i]]
        xi = c()
        for(j in 1:n)
        {
          if(x[2] < Xn[j,1])
          {
            xi = c(xi,X[j,2])
          }
          else if(x[1] > Xn[j,2])
          {
            xi = c(xi,X[j,1])
          }
          else
          {
            xi = c(xi,X[j,1],X[j,2])
          }
        }
        
        xi_list[[i]] = xi
      }
      
      M = unlist(lapply(xi_list,function(X) mean(X)))
      vars = c()
      for(i in 1:nl)
      {
        x = ints_left[[i]]
        if(!is.numeric(M[i]))
        {
          vars = c(vars,NA)
        }
        else if(M[i] >= x[1] & M[i] <= x[2])
        {
          vars = c(vars,var(xi_list[[i]]))
        }
      }
      
      vars = na.omit(vars)
      
      if(length(vars) == 0)
      {
        return(0)
      }
      else
      {
        return(max(vars))
      }
      
    }
    
  
}

u_var_exp = function(X)
{
  M = expand.grid(as.list(as.data.frame(t(X))))
  L = as.list(as.data.frame(t(M)))
  VAR = lapply(L,function(x) var(x))
  return(max(unlist(VAR)))
  
}

int_var = function(X)
{
  return(c(l_var(X),u_var(X)))
}

l_order = function(XI,XJ)
{
  if(XI[1] < XJ[1])
  {
    return(T)
  }
  else if(XI[1] == XJ[1] & XI[2] <= XJ[2])
  {
    return(T)
  }
  else
  {
    return(F)
  }
}

swap_if_larger = function(pair) {
  if(larger(pair)) {
    return(rev(pair)) 
  } else {
    return(pair)
  }
}



lex_sort = function(X) {
  no_passes = 0
  
  while(1) 
  {
    no_swaps = 0

    for (j in 1 : (nrow(X) - 1 - no_passes)) 
    {

      if (l_order(X[j,],X[j + 1,])) 
      {
   
        temp = X[j,]
        X[j,] = X[j+1,]
        X[j+1,] = temp
        no_swaps = no_swaps + 1
      }
    }
    no_passes = no_passes + 1
    if(no_passes == nrow(X) - 1)
    {
      break
    }
    if(no_swaps == 0) 
    {
      break
    }
  }
  X[rev(1:nrow(X)),]
}

lexicographic = function(X)
{
  Y = data.frame(t(rep(NA,2)))
  flag = T
  
  
  
  while(flag)
  {
    flag = F
    for(j in 1:(nrow(X)-1))
    {

      if(l_order(X[j,],X[j+1,]))
      {
        temp = X[j,]
        X[j,] = X[j+1,]
        X[j+1,] = temp
        flag = T
      }
    }
  }
  
  return(X[rev(1:nrow(X)),])
}

nonest_test = function(X)
{
  val = c()
  for(i in 1:(nrow(X)-1))
  {
    if(X[i + 1,2] < X[i,2])
    {
      val[i] = 1
      
    }
    else
    {
      val[i] = 0
    }
  }
  
  S = sum(val)
  
  if(S == 0)
  {
    return(T)
  }
  else
  {
    return(F)
  }
  
  
}

u_var_nonest = function(X)
{
  N = nrow(X)

  X = lex_sort(X)

  
  if(nonest_test(X))
  {
    M = 0
    E = 0
    v = c()
    for(k in 2:(N-1))
    {

      v = c(v,var(c(X[1:k,1],X[(k+1):N,2])))
    }
    
    return(max(v))
  }
  else
  {
    X = un_nest(X)
    M = 0
    E = 0
    v = c()
    for(k in 2:(N-1))
    {
      
      v = c(v,var(c(X[1:k,1],X[(k+1):N,2])))
    }
    
    return(max(v))
  }
  
}



final_var = function(X)
{
  return(c(l_var(X),u_var_nonest(X)))
}

un_nest = function(X)
{
  # THIS ASSUMES X IS LEXICOGRAPHICALLY ORDERED
  
  for(i in 1:(nrow(X)-1))
  {
    if(X[i + 1,2] < X[i,2])
    {
      X[i + 1,2] = X[i,2] + 0.0000001
      
    }
    else
    {
      next
    }
  }
  
  return(X)
}

UminusV_sq_INT = function(u,v)
{
  xl = u[1] - v[2]
  xh = u[2] - v[1]
  
  if(xl <= 0)
  {
    fl = 0
  }
  else
  {
    fl = xl^2
  }
  
  fu = max(xl^2,xh^2)
  
  return(c(fl,fu))
}

UminusV_abs_INT = function(u,v)
{
  xl = u[1] - v[2]
  xh = u[2] - v[1]
  
  if(xl <= 0)
  {
    fl = 0
  
  }
  else
  {
    fl = xl
  }
  
  if(xl <= 0)
  {
    if(xh >= 0)
    {
      fu = max(-xl,xh)
    }
    else
    {
      fu = -xl
    }
  }
  else
  {
    fu = xh
  }
  

  
  return(c(fl,fu))

  
  
}

Envelope = function(X)
{
  return(c(min(X[,1]),max(X[,2])))
}

##################################################
# Function implements subinterval reconstitution #
##################################################

split_int = function(X,n)
{
  L = lapply(1:n,function(x) {
    return(c(X[1] + ((x-1)/n)*(X[2] - X[1]),X[1] + ((x)/n)*(X[2] - X[1])))
  })
  
  return(L)
}

subint_recon = function(FUN,vars,n)
{
  library(gtools)


  sub_ints = lapply(vars,function(X){
    return(split_int(X,n))
  })
  

  
  C = permutations(n = n,r = length(vars),repeats.allowed = T)
  
  
  for(i in 1:nrow(C))
  {
    A = cbind(1:length(vars),C[i,])

    list_v = list()
    
    for(j in 1:nrow(A))
    {
      list_v[[j]] = sub_ints[[A[j,1]]][[A[j,2]]]
    }
    

    X = do.call(FUN,list_v)

    
    if(i == 1)
    {
      xl = X[1]
      xu = X[2]
    }
    else
    {
      if(X[1] < xl)
      {
        xl = X[1]
      }
      
      if(X[2] > xu)
      {
        xu = X[2]
      }
    }
  }
  
  return(c(xl,xu))
}


