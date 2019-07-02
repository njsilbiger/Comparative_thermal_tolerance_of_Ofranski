model{
  #priors
  # variance for observations
  sigma.obs~dunif(0,100)
  tau.obs<-1/sigma.obs^2
  
  # variance between genotypes
  sigma.geno~dunif(0,100)
  tau.geno<-1/sigma.geno^2
  # random intercept for genotype

  # multivariate normal priors
  mu_raw[1]~dnorm(0,1E-2)
  mu_raw[2]~dnorm(320,.01)
  mu_raw[3]~dgamma(1,1)
  mu_raw[4]~dgamma(1,1)
  
  
  
  for (i in 1:4){
    xi[i]~dunif(0,1) 
    mu[i]<-xi[i]*mu_raw[i]
  }
  
  mu.lnc<-mu[1]
  mu.Th<-mu[2]
  mu.Eh<-mu[3] # Eh has to be bigger than E
  mu.E<-mu[4]
  
  #correlated errors for varying slopes using inverse Wishart
  Tau_B_raw[1:4,1:4]~dwish(W[,],5)
  Sigma_B_raw[1:4,1:4]<-inverse(Tau_B_raw[,])
  
  Tc<-300.15 # common temperature 27 degrees
  
  #likelihood
  
  #fine scale for every observation
  for (i in 1:N){
    
    mu.obs[i]<-lnc.rate[rate[i]] + lnc.geno[geno[i]]+ log(exp(E.rate[rate[i]]/8.62e-5*(1/Tc - 1/temp[i]))) + log(1/(1 + exp(Eh.rate[rate[i]]/8.62e-5*(1/Th.rate[rate[i]] - 1/temp[i]))))
    
    y[i]~dnorm(mu.obs[i], tau.obs) # the data
    
    #for posterior predictive checks
    y.new[i]~dnorm(mu.obs[i], tau.obs) # new data for posterior predictive checks
    
  }
  
  # Genotype
  for (g in 1:G){
    lnc.geno[g]~dnorm(lnc.rate[ind_per_group[g]], tau.geno) #between genotypic variance across all genotypes.  Need to change variance some how to get different one per location 
  }
  
  #  broad scale (between the 3 rates) 
  
  #Multivariate distribution for varying coefs
  for (r in 1:R){ #for every Rate
    
    #multivariate normal distribution 
    B.hat[r,1]<-mu_raw[1] 
    B.hat[r,2]<-mu_raw[2]
    B.hat[r,3]<-mu_raw[3]
    B.hat[r,4]<-mu_raw[4]
    
    B[r,1:4]~dmnorm(B.hat[r,],Tau_B_raw[,])
    
    lnc.rate[r]<-xi[1]*B[r,1] 
    Th.rate[r]<-xi[2]*B[r,2]
    Eh.rate[r]<-xi[3]*B[r,3]
    E.rate[r]<-xi[4]*B[r,4]
  }
  #correlated errors using a Scaled inverse-Wisshart Model
  for (i in 1:4){
    sigma[i]<-xi[i]*sqrt(Sigma_B_raw[i,i])
  }
  
  sigma.lnc<-sigma[1]
  sigma.Th<-sigma[2]
  sigma.Eh<-sigma[3]
  sigma.E<-sigma[4]
  
  #correlation coefs
  for (i in 1:4){
    for (j in 1:4){
      rho[i,j]<-Sigma_B_raw[i,j]/sqrt(Sigma_B_raw[i,i]*Sigma_B_raw[j,j])
    }
  }
  
  
  for(g in 1:G) {lnc.pred[g] <-  lnc.rate[ind_per_group[g]]+lnc.geno[g] } # get predicted lnc values
  
  ## for plotting ##something not quite right with the intercepts.. they are all shifted down
  lnc.hat<-c(mean(lnc.pred[c(1,4,7,10)]), mean(lnc.pred[c(2,5,8,11)]),mean(lnc.pred[c(3,6,9,12)])) # to plot one curve per location
  

  
  for(i in 1:length(x.hat)){
    # predictions for figures
    for (r in 1:R){
      y.hat[i,r]<- lnc.hat[r] + E.rate[r]/8.62e-5*(1/Tc - 1/x.hat[i]) + log(1/(1 + exp(Eh.rate[r]/8.62e-5*(1/Th.rate[r] - 1/x.hat[i]))))

    }
    
  }
  
  for (r in 1:R){ # for every location
    # topt for each genotype
    C[r]<-ifelse(Eh.rate[r]/E.rate[r]<=1,1.1,Eh.rate[r]/E.rate[r]) #Eh/E has to be greater than 1. This stops the model from failing
    Topt[r] <- ((Eh.rate[r]*Th.rate[r])/(Eh.rate[r] + (8.62e-05 *Th.rate[r]*log((C[r]) - 1))))
    # calculate the lograte at Topt (maximum rate)
    
    a[r]<-lnc.rate[r] +  E.rate[r]/8.62e-5*(1/Tc - 1/Topt[r]) + log(1/(1 + exp(Eh.rate[r]/8.62e-5*(1/Th.rate[r] - 1/Topt[r]))))
    
    # calculating breadth which is 80% rate at Topt
    # calculate functional dropout at 90% the max rate
    fd[r]<-a[r]*0.9 # rate at functional dropout
    
    
    #temperature at functional dropout (high temp)
    Tmax[r]<- -1*(8.62e-5*(log(ifelse(exp(fd[r])>=1,exp(fd[r]),1.1)-1)/Eh.rate[r])-1/(1/Th.rate[r]))
  }
  
  #posterior predictive checks-- calculate Bayesian p-value
  #pvalue CV
  cv.y <- sd(y[])/mean(y[])
  cv.y.rep <- sd(y.new[])/mean(y.new[])
  pvalue.cv <- step(cv.y.rep-cv.y) # find Bayesian P for CoV
  # #value--the mean of many 0's and 1's returned by
  #the step function, one for each step in the chain
  
  #pvalue mean
  mean.y <-mean(y[])
  mean.y.rep <-mean(y.new[])
  pvalue.mean <-step(mean.y.rep - mean.y)
  # pvalue sd
  
  sigma.y<-sd(y[])
  sigma.y.new<-sd(y.new[])
  pvalue.sd<-step(sigma.y.new - sigma.y)
  
  
}#end model

