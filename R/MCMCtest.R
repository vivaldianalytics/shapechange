library(coda)
library(abind)

progbar = function(x,Y,pb)
{
  div = Y%/%10
  num = x%/%div
  if(num !=0)
  {
    for(i in 1:num)
    {
      setTxtProgressBar(pb,i)
    }
  }
  close(pb)
}
mc_cent = function(data,time,name)
{
  d = dim(data)
  n = d[3]
  load("MCMC")
  N = names(MCMC)
  if(name %in% N)
  {
    C_final = MCMC[[name]]
    d2 = dim(C_final[,,1])
  }
  else
  {
    c = list()
    for(i in 1:20)
    {
      c[[i]] = cent_mod(data[,,i],time)
    }
    d2 = dim(c[[1]])

    C_final = abind(c,along = 3)
    
    MCMC[[name]] = C_final
    save(MCMC,file = "MCMC")
    
  }
  
  fin = array(0,dim = c(d2[1],d2[2],2))
  
  for(t in 1:d2[1])
  {
    z95 = HPDinterval(as.mcmc(t(C_final[t,,])))
    
    fin[t,,] = z95 
  }
  colnames(fin) = colnames(c[[1]])
  return(list(fin[,,1],fin[,,2]))
}

cent_to_int = function(C_final,p)
{
  d2 = dim(C_final)
  show(d2)
  fin = array(0,dim = c(d2[1],2,2))
  
  for(t in 1:d2[1])
  {
    z95 = HPDinterval(as.mcmc(t(C_final[t,,])),prob = p)
    show(z95[2,])
    show(fin[t,,])
    fin[t,,] = z95
  }
  dimnames(fin) = list(dimnames(C_final)[[1]],dimnames(C_final)[[2]],1:2)
  return(fin)
}

mc_cent_all = function(data,time,name)
{
    d = dim(data)
    n = d[3]
    c1 = cent_mod(data[,,1],time)
    dc = dim(c1)
    c = array(0,c(dc[1],dc[2],n))
    for(i in 1:n)
    {
      c[,,i] = cent_mod(data[,,i],time)
    }
    A = c
  return(A)
}

mc_SP = function(data,time,maxit,prog,int,p)
{
  d = dim(data)
  if(is.null(maxit))
  {
    n = d[3]
  }
  else
  {
    n = min(maxit,d[3])
    if(n <= 0 )
    {
      n = 2
    }
  }
  c1 = cent_mod_SP(data[,,1],time)
  dc = dim(c1)
  c = array(0,c(dc[1],dc[2],n))
  if(prog)
  {
    print("Running Monte Carlo Simulations...",quote = F)
  }
  pb = txtProgressBar(min = 0,max = 1,style = 3)
  for(i in 1:n)
  {
    
    c[,,i] = cent_mod_SP(data[,,i],time)
 
    if(prog)
    {
      setTxtProgressBar(pb,i/n)
    }
    

  }
  close(pb)
  A = c
  dimnames(A) = list(2:length(time),c("time","shape"),1:dim(A)[3])
  
  if(int)
  {
    A = cent_to_int(A,p)
  }
  return(A)
}

mc_FT = function(data,time,maxit = NULL,prog,int,p)
{
  d = dim(data)
  if(is.null(maxit))
  {
    n = d[3]
  }
  else
  {
    n = min(maxit,d[3])
    if(n <= 0 )
    {
      n = 2
    }
  }
  c1 = cent_mod_FT(data[,,1],time)
  dc = dim(c1)
  c = array(0,c(dc[1],dc[2],n))
  if(prog)
  {
    print("Running Monte Carlo Simulations...",quote = F)
  }
  pb = txtProgressBar(min = 0,max = 1,style = 3)
  for(i in 1:n)
  {
    
    c[,,i] = cent_mod_FT(data[,,i],time)
    
    if(prog)
    {
      setTxtProgressBar(pb,i/n)
    }
    
    
  }
  close(pb)
  A = c
  dimnames(A) = list(2:length(time),c("time","shape"),1:dim(A)[3])
  if(int)
  {
    A = cent_to_int(A,p)
  }
  return(A)
}

mc_JG = function(data,time,maxit = NULL,prog,int,p)
{
  d = dim(data)
  if(is.null(maxit))
  {
    n = d[3]
  }
  else
  {
    n = min(maxit,d[3])
    if(n <= 0 )
    {
      n = 2
    }
  }
  c1 = cent_mod_JG(data[,,1],time)
  dc = dim(c1)
  c = array(0,c(dc[1],dc[2],n))
  if(prog)
  {
    print("Running Monte Carlo Simulations...",quote = F)
  }
  pb = txtProgressBar(min = 0,max = 1,style = 3)
  for(i in 1:n)
  {
    
    c[,,i] = cent_mod_JG(data[,,i],time)
    
    if(prog)
    {
      setTxtProgressBar(pb,i/n)
    }
    
    
  }
  close(pb)
  A = c
  dimnames(A) = list(2:length(time),c("time","shape"),1:dim(A)[3])
  if(int)
  {
    A = cent_to_int(A,p)
  }
  return(A)
}


mc_L = function(data,time,maxit = NULL,prog,int,p)
{
  d = dim(data)
  if(is.null(maxit))
  {
    n = d[3]
  }
  else
  {
    n = min(maxit,d[3])
    if(n <= 0 )
    {
      n = 2
    }
  }
  c1 = cent_mod_L(data[,,1],time)
  dc = dim(c1)
  c = array(0,c(dc[1],dc[2],n))
  if(prog)
  {
    print("Running Monte Carlo Simulations...",quote = F)
  }
  pb = txtProgressBar(min = 0,max = 1,style = 3)
  for(i in 1:n)
  {
    
    c[,,i] = cent_mod_L(data[,,i],time)
    
    if(prog)
    {
      setTxtProgressBar(pb,i/n)
    }
    
    
  }
  close(pb)
  A = c
  dimnames(A) = list(2:length(time),c("time","shape"),1:dim(A)[3])
  if(int)
  {
    A = cent_to_int(A,p)
  }
  return(A)
}

mc_BCF = function(data,time,maxit = NULL,prog,int,p)
{
  d = dim(data)

  if(is.null(maxit))
  {
    n = d[3]
  }
  else
  {
    n = min(maxit,d[3])
    if(n <= 0 )
    {
      n = 2
    }
  }
  c1 = cent_mod_BCF(data[,,1],time)
  dc = dim(c1)
  c = array(0,c(dc[1],dc[2],n))
  if(prog)
  {
    print("Running Monte Carlo Simulations...",quote = F)
  }
  pb = txtProgressBar(min = 0,max = 1,style = 3)
  for(i in 1:n)
  {
    
    c[,,i] = cent_mod_BCF(data[,,i],time)
    if(prog)
    {
      setTxtProgressBar(pb,i/n)
    }
    
    
    
  }
  close(pb)
  A = c
  dimnames(A) = list(2:length(time),c("time","shape"),1:dim(A)[3])
  if(int)
  {
    A = cent_to_int(A,p)
  }
  return(A)
}




cent_mod = function(data,time)
{

  
  pgrowth = pgr_frame(data,time)
  rab = rel_abundance(data)
  uns = unit_sphere(data)
  schange = sizechange(data,time)
  sc_spencer = scalar_shapechange_Spencer(data,time)
  
  sc_FT = shape_change_fostertilman(data,time)
  sc_JG = shape_change_jassbygoldman(data,time)
  sc_lew = shape_change_lewis(data,time)
  sc_BCF = shape_change_braycurtis_Field(data,time)
  
  
  
  final = cbind(time[2:length(time)],
                data[2:(nrow(data)),],
                pgrowth,
                rab[2:(nrow(rab)),],
                uns[2:(nrow(uns)),],
                schange,
                sc_spencer[,2],
                sc_FT[,2],
                sc_JG[,2],
                sc_lew[,2],
                sc_BCF[,2])
  colnames(final) = c("Time",
                      colnames(data),
                      paste(colnames(data),"prop"),
                      colnames(rab),
                      colnames(uns),
                      "Rate of Size Change",
                      "Population Shape Change: Spencer",
                      "Population Shape Change: Foster and Tilman",
                      "Population Shape Change: Jassby and Goldman",
                      "Population Shape Change: Lewis",
                      "Population Shape Change: Bray-Curtis/Field et al"
  )
  return(final)
}

cent_mod_SP = function(data,time)
{
  sc_spencer = scalar_shapechange_Spencer(data,time)
  final = cbind(time[2:length(time)],sc_spencer[,2])
  colnames(final) = c("Time","Population Shape Change: Spencer")
  return(final)
}

cent_mod_FT = function(data,time)
{
  sc_FT = shape_change_fostertilman(data,time)
  final = cbind(time[2:length(time)],sc_FT[,2])
  colnames(final) = c("Time","Population Shape Change: Foster and Tilman")
  return(final)
}

cent_mod_JG = function(data,time)
{
  sc_JG = shape_change_jassbygoldman(data,time)
  final = cbind(time[2:length(time)],sc_JG[,2])
  colnames(final) = c("Time","Population Shape Change: Jassby and Goldman")
  return(final)
}

cent_mod_L = function(data,time)
{
  sc_L = shape_change_lewis(data,time)
  final = cbind(time[2:length(time)],sc_L[,2])
  colnames(final) = c("Time","Population Shape Change: Lewis")
  return(final)
}

cent_mod_BCF = function(data,time)
{
  sc_L = shape_change_braycurtis_Field(data,time)
  final = cbind(time[2:length(time)],sc_L[,2])
  colnames(final) = c("Time","Population Shape Change: Bray-Curtis/Field et al")
  return(final)
}