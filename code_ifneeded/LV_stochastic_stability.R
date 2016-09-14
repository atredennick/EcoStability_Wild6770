
# Parameters----------------------------------------------

r=c(1.1,1.1)
Zvar=0.02  # both spp have the same variability in K
Zcov=0.01 # covariance between the two species' K's
alpha=matrix(NA,2,2)
alpha[1,1]=0.1
# alpha[1,2]=0.19 #low synchrony, high CV
# alpha[2,1]=0.09 #low synchrony, high CV
alpha[1,2]=0.01 #high synchrony, low CV
alpha[2,1]=0.01 #high synchrony, low CV
alpha[2,2]=0.2
initialN=c(5,5)
Time=500
tiny=0.001  #extinction threshold

# Function-------------------------------------------------

updateN=function(N,r,alpha,z){
    # N,r,K, and z are a vectors with 2 elements
    # alpha is a 2x2 matrix
    # syntax: out = updateN(N,r,K,alpha)
    # out is a vector with 2 elements (updated population
    #   size of spp 1 and spp2, respectively
    
    spp1 = N[1]+N[1]*(r[1]*(1-alpha[1,1]*N[1]-alpha[1,2]*N[2])+z[1])
    spp2 = N[2]+N[2]*(r[2]*(1-alpha[2,2]*N[2]-alpha[2,1]*N[1])+z[2])
    out = c(spp1,spp2)
    return(out)
}

# Loop --------------------------------------------------

N=matrix(NA,Time,2)
N[1,]=initialN

# make stochastic carrying capacity
library(mvtnorm)
sigma = matrix(data=c(Zvar,Zcov,Zcov,Zvar),nrow=2,ncol=2)
z = rmvnorm(Time,c(0,0),sigma)

for(i in 2:Time){
    tmp=updateN(N[i-1,],r,alpha,z[i,]) #note the indexing on K
    tmp[tmp<tiny] = 0
    N[i,] = tmp
}

# Plot output--------------------------------------------


community <- rowSums(N)
cv <- sd(community) / mean(community) # calculate CV

library(synchrony)
synch <- community.sync(N)[[1]]

N <- cbind(N,rowSums(N))
matplot(N,type="l", main=paste0("CV = ",round(cv,2), "; synch = ", round(synch,2)),
        col = c("dodgerblue", "coral", "grey35"), lty=1)
abline(h = mean(N[,1]), col="dodgerblue", lwd=2)
abline(h = mean(N[,1])+sd(N[,1]), col="dodgerblue", lwd=0.5)
abline(h = mean(N[,1])-sd(N[,1]), col="dodgerblue", lwd=0.5)
abline(h = mean(N[,2]), col="coral", lwd=2)
abline(h = mean(N[,2])+sd(N[,2]), col="coral", lwd=0.5)
abline(h = mean(N[,2])-sd(N[,2]), col="coral", lwd=0.5)
abline(h = mean(N[,3]), col="black", lwd=2)
abline(h = mean(N[,3])+sd(N[,3]), col="black", lwd=0.5)
abline(h = mean(N[,3])-sd(N[,3]), col="black", lwd=0.5)
