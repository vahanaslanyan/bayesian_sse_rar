#Adapted from https://github.com/timobrakenhoff/BayesianSSRwithPDCCPP/blob/master/Gamma_PDCCPP_procedure.R
library(tidyverse)


set.seed(219)
sample_size_calculation <- function(alpha_prior=1, beta_prior=49, eta=0.95, zeta=0.90, 
                                    xi=0.95,r=c(1/3,1/6,1/6,1/6,1/6),q_prior=c(10,2,2,2,2),delta_star=5) {
  # Counter  
  for (n in 1:100000) {
    # Calculated sample size
    N_treat<-n*r            
    # Posterior alpha
    alpha_posterior<-alpha_prior + 0.5*n
    # eta point on t-distribution
    t_eta<-qt(eta, 2*alpha_posterior)
    # zeta point on t-distribution
    t_zeta<- qt(zeta, 2*alpha_posterior)
    # xi point on beta
    beta_xi<- qbeta(xi, n/2,alpha_prior) 
    # posteriors of q
    q_posterior<-q_prior+N_treat
    # calculation of D
    D<-c()
    for(i in 2:length(q_posterior)){
      D<-c(D,(q_posterior[1]*q_posterior[i])/(q_posterior[1]+q_posterior[i]))  
    }
    D<-min(D)
    # left hand side of equation
    lhs<- D*(alpha_posterior/beta_prior)*(1-beta_xi)
    # right hand side of equation
    rhs<- ((t_eta+t_zeta)/delta_star)^2  
    # Corresponding critical value of xi (see eq. 10 in Brakenhoff et al 2018) 
    xi_crit<-1- (rhs*beta_prior)/(alpha_posterior*D)          
    xi_int<- pbeta(xi_crit,n/2,alpha_prior)                   
    
    #Output only the first row to satisfy the condition, if we reac 10000, output NA and exit
    if (lhs >= rhs) {                                         
      results<-N_treat
      names(results)<-c(paste0("treatment", seq_along(N_treat)))
      results<-as.data.frame(t(results))
      break
      break
    }else if (n==100000){
      results<-rep(NA,length(r))
      names(results)<-c(paste0("treatment", seq_along(N_treat)))
      results<-as.data.frame(t(results))
      break
      break
    }
    }
  return(results)
}

posterior_calculations<-function(alpha_prior=100,beta_prior=100,q_prior=c(10,2,2,2,2),mu_prior=c(0,0,0,0,0),N_treat=c(649, 162,162, 649,324),y_treatment){
  #y treatment should have columns treatment assignment and y
  #Prior nu
  nu_prior=rgamma(1,shape = alpha_prior,rate = beta_prior)
  #Posterior q
  q_posterior<-q_prior+N_treat
  #means per treatment groups
  y_treatment%>%group_by(treatment_assignment)%>%summarize(mean_treatment_group=mean(y))->treatment_effect_means
  #posterior mu
  mu_posterior<-mu_prior*q_prior/q_posterior+N_treat*treatment_effect_means$mean_treatment_group/q_posterior
  #posterior alpha
  alpha_posterior<- alpha_prior + 0.5*(sum(N_treat))
  #sum of squared treatment effects 
  y_treatment%>%group_by(treatment_assignment)%>%summarize(treatment_group_squared=sum(y^2))->treatment_effect_squared
  #H calculation
  H<-sum(treatment_effect_squared$treatment_group_squared+q_prior*mu_prior^2-q_posterior*mu_posterior^2)
  #posterior beta
  beta_posterior=beta_prior+H/2
  #calculation of D
  D<-c()
  for(i in 2:length(q_posterior)){
    D<-c(D,(q_posterior[1]*q_posterior[i])/(q_posterior[1]+q_posterior[i]))  # calculation of D
  }
  #q and mu have the same vector length
  results1<-as.data.frame(t(c(q_posterior,mu_posterior)))
  colnames(results1)<-c(paste0("q", seq_along(q_posterior)),paste0("mu", seq_along(mu_posterior)))
  #alpha posterior, beta prior, beta_posterior, and nu prior have the same length
  results2<-tibble(alpha_posterior,beta_prior,beta_posterior,nu_prior)
  #D has a length of length(N treat)-1
  D<-as.data.frame(D)
  #put results in a list and output it
  final_results<-list(q_andmu_posteriors=results1,alpha_beta_params=results2,D=D)
  return(final_results)
}

get_treatment_difference<-function(q_andmu_posteriors,D,alpha_beta_params){
  #extract relevant parameters from posterior parameters
  means<-q_andmu_posteriors[,grep("mu",colnames(q_andmu_posteriors))]
  means<-as.numeric(means)
  #center the means based on treatment 1, assuming treatment 1 is the placebo
  means<-means-means[1]
  #standard deviation
  D=min(D)
  sd1=1/sqrt(alpha_beta_params$nu_prior*D)
  #sample treatment effect values from their respective posterior distributions
  delta<-c()
  for(i in 1:length(means)){
    aux_delta<-mean(rnorm(1000,means[i],sd=sd1))
    delta<-c(delta,aux_delta)
  }
  #scale them so they follow normal(0,1) distribution
  delta_scaled<-(delta)/sd1
  return(delta_scaled)
}

allocation_calculation<-function(delta_scaled,n=1,N=2){
  #in normal distribution, probability of max is F(y)=P(y>X1,y>X2,y>X3,y>X4)=
  #=P(y>X1)*P(y>X2)*P(y>X2)*P(y>X3)*P(y>X4),and since they all come from the same
  #distribution, this turns into P(y>X_i)^n where n is the number of variables of interest
  #P(y>x_i) in R is pnorm(x_i)
  prob_2_max<-pnorm(delta_scaled)^length(delta_scaled)
  #allocation ratio is scaled based on number of patients recruited at interim analysis, and total number to be recruited
  r=prob_2_max^(n/(2*N))/sum(prob_2_max^(n/(2*N)))
  return(r)
}






 
