uncertainty_handler_sqrt = function(X,p,addval,addname)
{
  nr = nrow(X)
  nc = ncol(X)
  if(nr == 0 | nc == 0)
  {
    return(NULL)
  }
  Y = array(0,c(nr,nc,2))
  Y[,,1] = as.matrix(X) - p*as.matrix(sqrt(X))
  for(i in 1:nr)
  {
    for(j in 1:nc)
    {
      if(Y[i,j,1] <= 0)
      {
        Y[i,j,1] = 1
      }
    }
  }
  
  Y[,,2] = as.matrix(X) + p*as.matrix(sqrt(X))
  
  for(j in 1:nc)
  {
    Y[,j,1] = round(Y[,j,1],3)
    Y[,j,2] = round(Y[,j,2],3)
  }
  
  if(!is.null(addval))
  {
    d = dim(Y)
    
    Z = array(0,c(d[1],d[2] + 1,d[3]))
    Z[,,1] = cbind(Y[,,1],addval)
    Z[,,2] = cbind(Y[,,2],addval)
    dimnames(Z)[[2]] = c(colnames(X),addname)
    Y = Z
  }
  else
  {
    dimnames(Y)[[2]] = colnames(X)
  }
  
  
  
  return(Y)
}

uncertainty_handler_poisson = function(X,p,addval,addvalname)
{
  p. = p
  nr = nrow(X)
  nc = ncol(X)
  
  # Scott's interval functions
  Lpoisson_interval = function(x, conflevel) {
    0.5*qchisq((1-conflevel)/2, 2*x)
  }
  Rpoisson_interval = function(x, conflevel, ruleofthree=FALSE) {
    if (ruleofthree) return(0.5*qchisq(ifelse(x==0,conflevel,1-(1-conflevel)/2), 2*x+2)) else
      return(0.5*qchisq(1-(1-conflevel)/2, 2*x+2))
  }
  
  
  if(nr == 0 | nc == 0)
  {
    return(NULL)
  }
  
  Y = array(0,c(nr,nc,2))
  Y[,,1] = apply(t(as.matrix(X)),1,function(x) Lpoisson_interval(x,conflevel = p.))
  
  for(i in 1:nr)
  {
    for(j in 1:nc)
    {
      if(Y[i,j,1] <= 0)
      {
        Y[i,j,1] = 1
      }
    }
  }
  
  Y[,,2] = apply(t(as.matrix(X)),1,function(x) Rpoisson_interval(x,conflevel = p.))
  

  
  for(j in 1:nc)
  {
    Y[,j,1] = round(Y[,j,1],3)
    Y[,j,2] = round(Y[,j,2],3)
  }
  
  if(!is.null(addval))
  {
    d = dim(Y)
    
    Z = array(0,c(d[1],d[2] + 1,d[3]))
    Z[,,1] = cbind(Y[,,1],addval)
    Z[,,2] = cbind(Y[,,2],addval)
    dimnames(Z)[[2]] = c(colnames(X),addvalname)
    Y = Z
  }
  else
  {
    dimnames(Y)[[2]] = colnames(X)
  }
  
  return(Y)
}

as.interval.frame = function(data,exclude = NULL,method,param)
{
  if(is.null(exclude))
  {
    if(method == "sqrt")
    {
      data = uncertainty_handler_sqrt(data,param)
    }
    else if(method == "poisson")
    {
      data = uncertainty_handler_poisson(data,param)
    }
    else
    {
      show("ERROR: invalid method chosen.")
    }
  }
  else if(exclude %in% colnames(data))
  {
    M = match(exclude,colnames(data))
    exval= data[,M]
    data = data[,-M]
    

    if(method == "sqrt")
    {
      data = uncertainty_handler_sqrt(data,param,exval,exclude)
    }
    else if(method == "poisson")
    {
      data = uncertainty_handler_poisson(data,param,exval,exclude)
    }
    else
    {
      show("ERROR: invalid method chosen.")
    }
  }
}