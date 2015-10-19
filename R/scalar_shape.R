# Function that gives the proportional growth rate for a species.
# Input in the format of a 2-column data frame - first column is
# time, second column is abundance

# Formula is that given in Eq. 4 on pg. 4

prop_growth_rate = function(data,Time)
{
  X = data
  R = diff(log(X))
  R = R/diff(Time)
  return(R)
}

# Implementation of prop_growth_rate over entire data frame

pgr_frame = function(data,Time)
{
  df = data[1:(nrow(data)-1),]
  for(j in 1:ncol(df))
  {
    df[,j] = prop_growth_rate(data[,j],Time)
  }
  return(df)
}

#
# p sub i
#
rel_abundance = function(data)
{
  DF = data.frame(t(rep(NA,ncol(data))))
  for(i in 1:nrow(data))
  {
    for(j in 1:ncol(data))
    {
      DF[i,j] = data[i,j]/sum(data[i,])
    }
  }
  colnames(DF) = paste(colnames(data),"relative abundance")
  return(DF)
}

#
# e sub i
#

unit_sphere = function(data)
{
  DF = data.frame(t(rep(NA,ncol(data))))
  
  for(i in 1:nrow(data))
  {
    for(j in 1:ncol(data))
    {
      DF[i,j] = data[i,j]/sqrt(sum((data[i,])^2))
    }
  }
  colnames(DF) = paste(colnames(data),"unit sphere")
  return(DF)
}



sizechange = function(data,time)
{
  # Creates a list to store the proportional growth rates for each species
  rates = list()
  # For now I'm using for loops to ensure accuracy, later on we can vectorize these
  for(i in 1:ncol(data))
  {
    rates[[i]] = prop_growth_rate(data[,i],time)
  }
  # This is a vector that stores the mean growth rate 
  # across all species for each time interval.
  # Formula used is given in Eq. 6 on pg. 5
  mean_rate = c()
  # If there are nrow(data) measurements, there are *nrow(data) - 1* intervals
  for(t in 1:(nrow(data)-1))
  {
    # Take the average
    S = 0
    
    for(i in 1:length(rates))
    {
      S = S + rates[[i]][t]
    }
    
    mean_rate[t] = S/length(rates)
  }
  return(mean_rate)
}

# Function that gives the scalar rate of shape change (s sub r or s_r)
# Input is in the form of two items: a data frame where each column is
# abundance data for each species, and a vector of time values

scalar_shapechange_Spencer = function(data,time)
{
  # Creates a list to store the proportional growth rates for each species
  rates = pgr_frame(data,time)
  # This stores the scalar rate of shape change for each interval
  shape = c()
  for(t in 1:(nrow(data) - 1))
  {
    M = sum(rates[t,])/length(rates[t,])
    S = 0
    d = ncol(data) - 1
    for(j in 1:ncol(data))
    {
      S = S + (rates[t,j] - M)^2
    }
    shape[t] = sqrt(S/d)
  }
  # Creates a 2-column data frame, first column is the time (removing the first time instance)
  frame = cbind(time[2:length(time)],shape)
  colnames(frame) = c("time","shape")
  return(frame)
}

cos_theta = function(data,time)
{
  # Creates a list to store the proportional growth rates for each species
  rates = list()
  
  # For now I'm using for loops to ensure accuracy, later on we can vectorize these
  
  for(i in 1:ncol(data))
  {
    rates[[i]] = prop_growth_rate(data[,i],time)
  }
  
  # This is a vector that stores the mean growth rate 
  # across all species for each time interval.
  # Formula used is given in Eq. 6 on pg. 5
  mean_rate = c()
  
  # If there are nrow(data) measurements, there are *nrow(data) - 1* intervals
  for(t in 1:(nrow(data)-1))
  {
    # Take the average
    S = 0
    
    for(i in 1:length(rates))
    {
      S = S + rates[[i]][t]
    }
    
    mean_rate[t] = S/length(rates)
  }
  
  # This stores the scalar rate of shape change for each interval
  theta = c()
  
  for(t in 1:(nrow(data) - 1))
  {
    n = length(rates)
    u = rep(1/sqrt(n),n)
    r = c()
    for(i in 1:length(rates))
    {
      r[i] = rates[[i]][t]
    }
    
    theta[t] = 1/sum( r^2 )
    
  }
  
  # Creates a 2-column data frame, first column is the time (removing the first time instance)
  
  frame = cbind(time[2:length(time)],theta)
  
  colnames(frame) = c("time","shape")
  
  return(frame)
}

shape_change_fostertilman = function(data,time)
{
  sc = c()
  
  r.df = rel_abundance(data)
  
  for(i in 1:(length(time)-1))
  {
    sc[i] = (1/(time[i + 1] - time[i]))*sqrt(sum((r.df[i + 1,] - r.df[i,])^2))
  }
  
  frame = cbind(time[2:length(time)],sc)
  colnames(frame) = c("time","shape")
  
  return(frame)
}

shape_change_jassbygoldman = function(data,time)
{
  s.df = unit_sphere(data)
  
  sc = c()
  
  for(i in 1:(length(time)-1))
  {
    sc[i] = (1/(time[i + 1] - time[i]))*sqrt( sum(  (s.df[i + 1,] - s.df[i,])^2   )  )
  }
  
  frame = cbind(time[2:length(time)],sc)
  colnames(frame) = c("time","shape")
  
  return(frame)
}

shape_change_lewis = function(data,time)
{
  sc = c()
  
  r.df = rel_abundance(data)
  for(i in 1:(length(time)-1))
  {
    sc[i] = (1/(time[i + 1] - time[i]))*sum(abs(  r.df[i + 1,] - r.df[i,]  ))
  }
  
  frame = cbind(time[2:length(time)],sc)
  colnames(frame) = c("time","shape")
  
  return(frame)
}

shape_change_braycurtis_Field = function(data,time)
{
  r.df = rel_abundance(data)
  kappa = r.df
  
  for(i in 1:nrow(r.df))
  {
    kappa[i,] = (r.df[i,])^(1/4)
  }
  sc = c()
  for(i in 1:(length(time)-1))
  {
    dti = (1/(time[i + 1] - time[i]))
    sc[i] = sum(  abs( kappa[i + 1,] - kappa[i,]) )/(dti*sum(kappa[i + 1,] + kappa[i,]))
  }
  
  frame = cbind(time[2:length(time)],sc)
  colnames(frame) = c("time","shape")
  
  return(frame)
}

shape_change_chi = function(data,time)
{
  x.. = sum(data)
  
  x.i = c()
  
  r.df = rel_abundance(data)
  
  for(j in 1:ncol(data))
  {
    x.i = sum(data[,j])
  }
  
  sc = c()
  
  for(t in 1:(length(time)-1))
  {
    S = 0
    dti = (1/(time[t + 1] - time[t]))
    
    for(j in 1:ncol(data))
    {
      S = S + (1/x.i[j])*((r.df[t + 1,j] - r.df[t,j])^2)
    }
    
    sc[t] = (sqrt(S)*x..)/dti
  }
  
  frame = cbind(time[2:length(time)],sc)
  
  colnames(frame) = c("time","shape")
  
  return(frame)
  
  
  
}