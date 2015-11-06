bindit = function(X,Y)
{
  A = array(0,c(dim(X),2))
  A[,,1] = X
  A[,,2] = Y
  return(A)
}

prop_growth_rate_INT = function(dataL,dataH,Time)
{
  # Data passed to this must be separated from the time variable

  nr = length(dataL)
  
  rL = c()
  rU = c()
  
  xL = dataL
  xU = dataH
  for(i in 2:nr)
  {
    rU[i-1] = (log(xU[i]) - log(xL[i - 1]))/(Time[i] - Time[i - 1])
    rL[i-1] = (log(xL[i]) - log(xU[i - 1]))/(Time[i] - Time[i - 1])
  }
  
  return(list(rL,rU))
  
}

pgr_frame_INT = function(data,t)
{
  
  df = data[1:(nrow(data)-1),,]
  
  for(j in 1:ncol(df))
  {
    
    pgr = prop_growth_rate_INT(data[,j,1],data[,j,2],t)
    df[,j,1] = pgr[[1]]
    df[,j,2] = pgr[[2]]
  }
  return(df)
}

rel_abundance_INT = function(data)
{
  # Data passed to this must be separated from the time variable
  
  data.layer.1 = as.data.frame(data[,,1])
  data.layer.2 = as.data.frame(data[,,2])
  
  DF_L = data.frame(t(rep(NA,ncol(data))))
  DF_U = data.frame(t(rep(NA,ncol(data))))
  
  
  
  for(t in 1:nrow(data))
  {
    for(i in 1:ncol(data))
    {
      DF_L[t,i] = 1/(1 + (sum(data.layer.2[t,]) - data.layer.2[t,i])/data.layer.1[t,i])
      DF_U[t,i] = 1/(1 + (sum(data.layer.1[t,])- data.layer.2[t,i])/data.layer.2[t,i])
    }
  }
  
  
  colnames(DF_L) = colnames(DF_U) = paste(colnames(data),"relative abundance")
  return(list(DF_L,DF_U))
}

unit_sphere_INT = function(data)
{
  # Data passed to this must be separated from the time variable
  data.layer.1 = data[,,1]
  data.layer.2 = data[,,2]
  
  DF_L = data.frame(t(rep(NA,ncol(data))))
  DF_U = data.frame(t(rep(NA,ncol(data))))
  
  for(t in 1:nrow(data))
  {
    for(i in 1:ncol(data))
    {
      DF_L[t,i] = 1/(1 + sum(data.layer.2[t,-i]^2)/(data.layer.1[t,i])^2)
      DF_U[t,i] = 1/(1 + sum(data.layer.1[t,-i]^2)/(data.layer.2[t,i])^2)
    }
  }
  colnames(DF_L) = colnames(DF_U) = paste(colnames(data),"unit sphere")
  return(list(DF_L,DF_U))
}

sizechange_INT = function(data,time)
{
  # Creates a list to store the proportional growth rates for each species
  rates = list()
  
  # For now I'm using for loops to ensure accuracy, later on we can vectorize these
  for(i in 1:ncol(data))
  {
    rates[[i]] = prop_growth_rate_INT(data[,i,1],data[,i,2],time)
  }
  
  # These are vectors that store the mean growth rate (Lower and upper bounds)
  # across all species for each time interval.
  # Formula used is given in Eq. 6 on pg. 5
  mean_rate_L = c()
  mean_rate_U = c()
  
  # If there are nrow(data) measurements, there are *nrow(data) - 1* intervals
  
  for(t in 1:(nrow(data)-1))
  {
    # Take the average
    S_L = 0
    S_U = 0
    
    for(i in 1:length(rates))
    {
      # rates[[i]] is a list, with two vectors in it, the first one is the lower bound, the other is upper
      S_L = S_L + rates[[i]][[1]][t]
      S_U = S_U + rates[[i]][[2]][t]
    }
    
    mean_rate_L[t] = S_L/length(rates)
    mean_rate_U[t] = S_U/length(rates)
  }
  return(list(mean_rate_L,mean_rate_U))
}

scalar_shapechange_Spencer_INT = function(data,time)
{
  # Creates a list to store the proportional growth rates for each species
  rates = list()
  
  # For now I'm using for loops to ensure accuracy, later on we can vectorize these
  
  nm = pgr_frame_INT(data,time)

  
  
  shape = data.frame(t(rep(NA,3)))
  colnames(shape) = c("time","LS","US")
  

  
  for(t in 1:(nrow(data)-1))
  {
    shape[t,1] = time[t]

    s = final_var(nm[t,,])

    shape[t,2] = s[1]
    shape[t,3] = s[2]
    
  }
  A = bindit(cbind(time[2:length(time)],sqrt(shape[,2])),cbind(time[2:length(time)],sqrt(shape[,3])))
  dimnames(A) = list(2:length(time),c("time","shape"))
  return(A)
}

shape_change_fostertilman_INT = function(data,time)
{
  r.m = rel_abundance_INT(data)
  
  
  L = r.m[[1]]
  U = r.m[[2]]
  
  
  shapeL = c()
  shapeU = c()
  
  for(t in 1:(nrow(data)-1))
  {
    SL = 0
    SU = 0
    
    for(j in 1:ncol(data))
    {
      u = c(L[t+1,j],U[t+1,j])
      v = c(L[t,j],U[t,j])
      add = UminusV_sq_INT(u,v)
      SL = SL + add[1]
      SU = SU + add[2]
    }
    
    delta_t = time[t + 1] - time[t]
    shapeL[t] = (1/delta_t)*sqrt(SL)
    shapeU[t] = (1/delta_t)*sqrt(SU)
  }
  A = bindit(cbind(time[2:length(time)],sqrt(shapeL)),cbind(time[2:length(time)],sqrt(shapeU)))
  dimnames(A) = list(2:length(time),c("time","shape"))
  return(A)
}

shape_change_jassbygoldman_INT = function(data,time)
{
  uns = unit_sphere_INT(data)
  
  L = uns[[1]]
  U = uns[[2]]
  
  
  shapeL = c()
  shapeU = c()
  
  for(t in 1:(nrow(data)-1))
  {
    SL = 0
    SU = 0
    
    for(j in 1:ncol(data))
    {
      u = c(L[t+1,j],U[t+1,j])
      v = c(L[t,j],U[t,j])
      add = UminusV_sq_INT(u,v)
      SL = SL + add[1]
      SU = SU + add[2]
    }
    
    delta_t = time[t + 1] - time[t]
    shapeL[t] = (1/delta_t)*sqrt(SL)
    shapeU[t] = (1/delta_t)*sqrt(SU)
  }
  
  A = bindit(cbind(time[2:length(time)],sqrt(shapeL)),cbind(time[2:length(time)],sqrt(shapeU)))
  dimnames(A) = list(2:length(time),c("time","shape"))
  return(A)
}

shape_change_lewis_INT = function(data,time)
{
  r.m = rel_abundance_INT(data)
  
  L = r.m[[1]]
  U = r.m[[2]]
  
  
  shapeL = c()
  shapeU = c()
  
  for(t in 1:(nrow(data)-1))
  {
    SL = 0
    SU = 0
    
    for(j in 1:ncol(data))
    {
      u = c(L[t+1,j],U[t+1,j])
      v = c(L[t,j],U[t,j])
      add = UminusV_abs_INT(u,v)
      SL = SL + add[1]
      SU = SU + add[2]
    }
    
    delta_t = time[t + 1] - time[t]
    shapeL[t] = (1/delta_t)*SL
    shapeU[t] = (1/delta_t)*SU
  }
  
  A = bindit(cbind(time[2:length(time)],sqrt(shapeL)),cbind(time[2:length(time)],sqrt(shapeU)))
  dimnames(A) = list(2:length(time),c("time","shape"))
  return(A)
}

shape_change_braycurtis_Field_INT = function(data,time)
{
  
  r.m = rel_abundance_INT(data)
  
  L = r.m[[1]]
  U = r.m[[2]]
  
  f = function(a,b)
  {
    n = nrow(a)
    
    A1 = A2 = B1 = B2 = rep(0,n)
    
    for (i in 1:n) 
    {
      if (max(a[i,1],b[i,1]) <= min(a[i,2], b[i,2])) 
      {
        A1[[i]] = B1[[i]] = max(a[i,1],b[i,1])
      }else if (abs(a[i,1]-b[i,2]) < abs(a[i,2]-b[i,1])) 
      {
        A1[[i]] = a[i,1]
        B1[[i]] = b[i,2]
      } else 
      {
        A1[[i]] = a[i,2]
        B1[[i]] = b[i,1]
      }
    }
      
    
    for (i in 1:n)
    {
      if (abs(a[i,1]-b[i,2]) < abs(a[i,2]-b[i,1])) 
      {
        A2[[i]] = a[i,2]
        B2[[i]] = b[i,1]
      } else 
      {
        A2[[i]] = a[i,1]
        B2[[i]] = b[i,2]
      }
      
      
    }
    
    exact = c(sum(abs(A1 - B1))/sum(A1+B1), sum(abs(A2 - B2))/sum(A2+B2))
    return(exact)
    
  }
  
  shapeL = c()
  shapeU = c()
  
  for(t in 1:(nrow(data)-1))
  {
    a = cbind(t(L[t + 1,]),t(U[t + 1,]))
    b = cbind(t(L[t,]),t(U[t,]))

    
    I = f(a,b)
    
    delta = time[t + 1] - time[t]
    I = (1/delta)*I
    
    shapeL[t] = I[1]
    shapeU[t] = I[2]
  }
  
  A = bindit(cbind(time[2:length(time)],sqrt(shapeL)),cbind(time[2:length(time)],sqrt(shapeU)))
  dimnames(A) = list(2:length(time),c("time","shape"))
  return(A)
}
