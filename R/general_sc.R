source("R/scalar_shape_int.R")
source("R/MCMCtest.R")
source("R/successionmeasures.surreal2.R")
source("R/interval_variance.R")
source("R/scalar_shape.R")
detect_type = function(data)
{
  type = typeof(data)
  if(type == "list")
  {
    type = "scalar"
  }
  else if(type == "double")
  {
    d = dim(data)
    if(length(d) == 2)
    {
      type = "scalar"
    }
    else if(length(d) > 2)
    {
      if(d[3] == 2)
      {
        type = "interval"
      }
      else if(d[3] > 2)
      {
        type = "mcmc"
      }
      else
      {
        return(NULL)
      }
    }
    else
    {
      return(NULL)
    }
  }
  else
  {
    return(NULL)
  }
  return(type)
}
shape_change = function(data,time_val,
                        method,surr = F,
                        iter = 40,log_transform = F,
                        mc_progress = F,mc_int = F,mc_prob = 0.95)
{
  type = detect_type(data)
  if(type == "scalar")
  {
    time = data[,time_val]
    data = data[,-match(time_val,colnames(data))]
  }
  else if(type == "interval")
  {
    time = data[,time_val,1]
    data = data[,-match(time_val,colnames(data)),]
  }
  else if(type == "mcmc")
  {
    time = data[,time_val,1]
    data = data[,-match(time_val,colnames(data)),]
  }
  
  if(log_transform)
  {
    data = exp(data)
  }
  
  if(is.null(type))
  {
    show("Wrong data format.")
    return(NULL)
  }
  
  if(method == "Spencer")
  {
    if(type == "scalar")
    {
      if(surr)
      {
        SC = getsurr(data,time)
        SC = SC$sr
        SC = cbind(time,SC)
        SC = as.data.frame(SC[2:nrow(SC),])
        colnames(SC) = c("time","real","infv")
      }
      else
      {
        data[data < 1] = 1
        SC = scalar_shapechange_Spencer(data,time)
      }
    }
    else if(type == "interval")
    {
      data[data < 1] = 1
      SC = scalar_shapechange_Spencer_INT(data,time)
    }
    else if(type == "mcmc")
    {
      
      
      
      data[data < 1] = 1
      SC = mc_SP(data = data,time =time,maxit =iter,prog = mc_progress,mc_int,mc_prob)
    }
    else
    {
      show("ERROR - data must take the form of a data frame, matrix, or 3-dimensional array")
      return(NULL)  
    }
  }
  else if(method == "Foster and Tilman")
  {
    if(type == "scalar")
    {
      data[data < 1] = 1
      SC = shape_change_fostertilman(data,time)
    }
    else if(type == "interval")
    {
      data[data < 1] = 1
      SC = shape_change_fostertilman_INT(data,time)
    }
    else if(type == "mcmc")
    {
      data[data < 1] = 1
      
      SC = mc_FT(data,time,maxit =iter,prog = mc_progress,mc_int,mc_prob)
    }
    else
    {
      show("ERROR - data must take the form of a data frame, matrix, or 3-dimensional array")
      return(NULL)
    }
  }
  else if(method == "Jassby and Goldman")
  {
    if(type == "scalar")
    {
      data[data < 1] = 1
      SC = shape_change_jassbygoldman(data,time)
    }
    else if(type == "interval")
    {
      data[data < 1] = 1
      SC = shape_change_jassbygoldman_INT(data,time)
    }
    else if(type == "mcmc")
    {
      data[data < 1] = 1
      SC = mc_JG(data,time,maxit =iter,prog = mc_progress,mc_int,mc_prob)
    }
    else
    {
      show("ERROR - data must take the form of a data frame, matrix, or 3-dimensional array")
      return(NULL)
    }
  }
  else if(method == "Lewis")
  {
    if(type == "scalar")
    {
      data[data < 1] = 1
      SC = shape_change_lewis(data,time)
    }
    else if(type == "interval")
    {
      data[data < 1] = 1
      SC = shape_change_lewis_INT(data,time)
    }
    else if(type == "mcmc")
    {
      data[data < 1] = 1
      SC = mc_L(data,time,maxit =iter,prog = mc_progress,mc_int,mc_prob)
    }
    else
    {
      show("ERROR - data must take the form of a data frame, matrix, or 3-dimensional array")
      return(NULL)
    }
  }
  else if(method == "Bray Curtis Field")
  {
    if(type == "scalar")
    {
      data[data < 1] = 1
      SC = shape_change_braycurtis_Field(data,time)
    }
    else if(type == "interval")
    {
      data[data < 1] = 1
      SC = shape_change_braycurtis_Field_INT(data,time)
    }
    else if(type == "mcmc")
    {
      data[data < 1] = 1
      SC = mc_BCF(data,time,maxit =iter,prog = mc_progress,mc_int,mc_prob)
    }
    else
    {
      show("ERROR - data must take the form of a data frame, matrix, or 3-dimensional array")
      return(NULL)
    }
  }
  else
  {
    show("Invalid method chosen. Please choose one of the following:")
    show("- Spencer")
    show("- Foster and Tilman")
    show("- Jassby and Goldman")
    show("- Lewis")
    show("- Bray Curtis Field")
    return(NULL)
    
  }
  return(SC)
}