model{
  #priors
  # variance for observations
  sigma.obs~dunif(0,100)
  tau.obs<-1/sigma.obs^2
  
  # variance between genotypes
  sigma.geno~dunif(0,100)
  tau.geno<-1/sigma.geno^2
  
  mu_raw[1]~dnorm(0,1E-2)
  mu_raw[2]~dnorm(320,.01)
  mu_raw[3]~dgamma(1,1)
  mu_raw[4]~dgamma(1,1)
  
  
  
  for (i in 1:4){
     xi[i]~dunif(0,1) # is this suppose to be 0,1.  was 0,100?
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
    
    mu.obs[i]<-lnc.loc[loc[i]] + lnc.geno[geno[i]]+ log(exp(E.loc[loc[i]]/8.62e-5*(1/Tc - 1/temp[i]))) + log(1/(1 + exp(Eh.loc[loc[i]]/8.62e-5*(1/Th.loc[loc[i]] - 1/temp[i]))))
 
    y[i]~dnorm(mu.obs[i], tau.obs) # the data
    
    #for posterior predictive checks
    y.new[i]~dnorm(mu.obs[i], tau.obs) # new data for posterior predictive checks
    
  }

   # Genotype
   for (g in 1:G){
       lnc.geno[g]~dnorm(lnc.loc[ind_per_group[g]], tau.geno) #between genotypic variance across all genotypes.  Need to change variance some how to get different one per location 
   }
   
#  broad scale (between locations) 
  
     #Multivariate distribution for varying coefs
     for (l in 1:L){ #for every Location
       
       #multivariate normal distribution 
       B.hat[l,1]<-mu_raw[1] 
       B.hat[l,2]<-mu_raw[2]
       B.hat[l,3]<-mu_raw[3]
       B.hat[l,4]<-mu_raw[4]
        
       B[l,1:4]~dmnorm(B.hat[l,],Tau_B_raw[,])
       
       lnc.loc[l]<-xi[1]*B[l,1] 
       Th.loc[l]<-xi[2]*B[l,2]
       Eh.loc[l]<-xi[3]*B[l,3]
       E.loc[l]<-xi[4]*B[l,4]
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
     
     
     for(g in 1:G) {lnc.pred[g] <-  lnc.loc[ind_per_group[g]]+lnc.geno[g] } # get predicted lnc values
     
     ## for plotting ##something not quite right with the intercepts.. they are all shifted down
     lnc.hat<-c(mean(lnc.pred[c(1,3,5,7)]), mean(lnc.pred[c(2,4,6,8)])) # to plot one curve per location
     lnc.g<-c(mean(lnc.geno[c(1,3,5,7)]), mean(lnc.geno[c(2,4,6,8)])) # to plot one curve per location
     #
     
     
     for(i in 1:length(x.hat)){
       # predictions for figures
       for (l in 1:L){
        y.hat[i,l]<- lnc.hat[l] + E.loc[l]/8.62e-5*(1/Tc - 1/x.hat[i]) + log(1/(1 + exp(Eh.loc[l]/8.62e-5*(1/Th.loc[l] - 1/x.hat[i]))))
        
         }
       
    }
     
    for (l in 1:L){ # for every location
     # topt for each genotype
        C[l]<-ifelse(Eh.loc[l]/E.loc[l]<=1,1.1,Eh.loc[l]/E.loc[l]) #Eh/E has to be greater than 1. This stops the model from failing
     Topt[l] <- ((Eh.loc[l]*Th.loc[l])/(Eh.loc[l] + (8.62e-05 *Th.loc[l]*log((C[l]) - 1))))
     # calculate the lograte at Topt (maximum rate)
          
          a[l]<-lnc.loc[l] +  E.loc[l]/8.62e-5*(1/Tc - 1/Topt[l]) + log(1/(1 + exp(Eh.loc[l]/8.62e-5*(1/Th.loc[l] - 1/Topt[l]))))
          
          # calculating breadth which is 80% rate at Topt
          # calculate functional dropout at 90% the max rate
          fd[l]<-a[l]*0.9 # rate at functional dropout
          
        
            #temperature at functional dropout (high temp)
          Tmax[l]<- -1*(8.62e-5*(log(ifelse(exp(fd[l])>=1,exp(fd[l]),1.1)-1)/Eh.loc[l])-1/(1/Th.loc[l]))
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

