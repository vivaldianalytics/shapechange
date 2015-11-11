# A = get(load("data/z.Rdata"))
# load("cnames")
# load("time")
# B = array(0,c(30,15,1000))
# 
# for(j in 1:1000)
# {
#   Ax = cbind(A[j,,],time)
#   colnames(Ax) = C
#   if(j ==1)
#   {
#     show(Ax)
#   }
#   B[,,j] = Ax
#   
# }
# 
# dimnames(B) = list(1:30,C,1:1000)
# hoverflies_MCMC = B
# 
# show(B[1,,])
# 
# save(hoverflies_MCMC,file = "data/hoverflies_MCMC.rda")