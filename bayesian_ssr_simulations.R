set.seed(219)
source("bayesian_ssr.R")

### Scenario 1




alpha<-1:10000
beta<-alpha
df_sample_size<-tibble()
df_posteriors<-tibble()
df_q_andmu_posteriors<-tibble()
df_alpha_beta_posteriors<-tibble()
df_D_posteriors<-tibble()
interim_allocation<-tibble()


for(i in 1:5000){
  alpha1<-sample(alpha,1)
  beta1<-beta[alpha1==alpha]
  aux<-sample_size_calculation(alpha_prior = alpha1,beta_prior = beta1, eta=0.95, zeta=0.90, 
                               xi=0.95,r=c(1/4,1/4,1/4,1/4),q_prior=c(2,2,2,2),delta_star=0.1)
  if(any(is.na(aux)))next
  y1_aux=rnorm(round(aux$treatment1/2),mean=0, sd=1)
  y2_aux=rnorm(round(aux$treatment2/2),mean=0.1, sd=1)
  y3_aux=rnorm(round(aux$treatment3/2),mean=0.2, sd=1)
  y4_aux=rnorm(round(aux$treatment4/2),mean=0.6, sd=1)
  
  y=c(y1_aux,y2_aux,y3_aux,y4_aux)
  treatment_assignment<-c(rep(1,round(aux$treatment1/2)),rep(2,round(aux$treatment2/2)),rep(3,round(aux$treatment3/2)),rep(4,round(aux$treatment4/2)))
  df=tibble(treatment_assignment,y)
  aux_post<-posterior_calculations(alpha_prior=alpha1,beta_prior=beta1,q_prior=c(2,2,2,2),mu_prior=c(0,0,0,0),
                                   N_treat = c(round(aux$treatment1/2),round(aux$treatment2/2),round(aux$treatment3/2),round(aux$treatment4/2)),
                                   y_treatment = df)
  df_q_andmu_posteriors<-rbind(df_q_andmu_posteriors,aux_post$q_andmu_posteriors)
  df_alpha_beta_posteriors<-rbind(df_alpha_beta_posteriors,aux_post$alpha_beta_params)
  df_D_posteriors<-rbind(df_D_posteriors,aux_post$D)
  df_sample_size<-rbind(df_sample_size,aux)
  treatment_differences<-get_treatment_difference(aux_post$q_andmu_posteriors,aux_post$D,aux_post$alpha_beta_params)
  new_r<-allocation_calculation(treatment_differences)
  interim_allocation<-rbind(interim_allocation,new_r)
  #df_posteriors<-rbind(df_posteriors,aux_post)
}
colMeans(df_sample_size)



interim_ss<-tibble()
for(i in 1:nrow(df_sample_size)){
  interim_aux<-sample_size_calculation(alpha_prior =df_alpha_beta_posteriors$alpha_posterior[i],beta_prior = df_alpha_beta_posteriors$beta_posterior[i], eta=0.95, zeta=0.90, 
                                       xi=0.95,r=c(1/4,1/4,1/4,1/4),q_prior =as.numeric( df_q_andmu_posteriors[i,grep("q",colnames(df_q_andmu_posteriors)),]),
                                       delta_star=0.2)
  interim_ss<-rbind(interim_ss,interim_aux)
}

colMeans(interim_ss)

#after interim analysis, we will need 461 per arm

colMeans(interim_allocation)


interim_ss_rar<-tibble()
for(i in 1:nrow(df_sample_size)){
  interim_aux<-sample_size_calculation(alpha_prior =df_alpha_beta_posteriors$alpha_posterior[i],beta_prior = df_alpha_beta_posteriors$beta_posterior[i], eta=0.95, zeta=0.90, 
                                       xi=0.95,r=as.numeric(interim_allocation[i,]),q_prior =as.numeric( df_q_andmu_posteriors[i,grep("q",colnames(df_q_andmu_posteriors)),]),
                                       delta_star=0.2)
  interim_ss_rar<-rbind(interim_ss_rar,interim_aux)
}

colMeans(interim_ss_rar,na.rm = T)


























#Lecanemab
#from sample size rationale:
#sd=2.031 variance is 4.124961
#nu_prior is 0.2424265
#so alpha/beta is 0.2424265

alpha<-1:10000
beta<-round(alpha/0.2424265)
df_sample_size<-tibble()
df_posteriors<-tibble()
df_q_andmu_posteriors<-tibble()
df_alpha_beta_posteriors<-tibble()
df_D_posteriors<-tibble()
interim_allocation<-tibble()
#true data
#se=0.1122449
#N=898
#sd=4.755529 (SE*sqrt(N))
#mean change=0.45
#effect size is 0.0946267


for(i in 1:1000){
  alpha1<-sample(alpha,1)
  beta1<-beta[alpha1==alpha]
  aux<-sample_size_calculation(alpha_prior = alpha1,beta_prior = beta1, eta=0.95, zeta=0.90, 
                               xi=0.95,r=c(0.5,0.5),q_prior=c(1,1),delta_star=0.373)
  if(any(is.na(aux)))next
  yE_aux=rnorm(round(aux$treatment1/2),mean=0.45, sd=4.755529)
  yC_aux=rnorm(round(aux$treatment2/2),mean=0, sd=4.755529)
  y=c(yC_aux,yE_aux)
  treatment_assignment<-c(rep(1,round(aux$treatment1/2)),rep(2,round(aux$treatment2/2)))
  df=tibble(treatment_assignment,y)
  aux_post<-posterior_calculations(alpha_prior=alpha1,beta_prior=beta1,q_prior=c(0.5,0.5),mu_prior =c(0,0),N_treat = c(round(aux$treatment1/2),round(aux$treatment2/2)),y_treatment = df)
  df_q_andmu_posteriors<-rbind(df_q_andmu_posteriors,aux_post$q_andmu_posteriors)
  df_alpha_beta_posteriors<-rbind(df_alpha_beta_posteriors,aux_post$alpha_beta_params)
  df_D_posteriors<-rbind(df_D_posteriors,aux_post$D)
  df_sample_size<-rbind(df_sample_size,aux)
  treatment_differences<-get_treatment_difference(aux_post$q_andmu_posteriors,aux_post$D,aux_post$alpha_beta_params)
  new_r<-allocation_calculation(treatment_differences)
  interim_allocation<-rbind(interim_allocation,new_r)
  #df_posteriors<-rbind(df_posteriors,aux_post)
}

colnames(interim_allocation)<-paste0("r", seq_along(new_r))
mean(df_sample_size$treatment1)

#Estimated sample size per arm 515

interim_ss<-tibble()
for(i in 1:nrow(df_sample_size)){
  interim_aux<-sample_size_calculation(alpha_prior =df_alpha_beta_posteriors$alpha_posterior[i],beta_prior = df_alpha_beta_posteriors$beta_posterior[i], eta=0.95, zeta=0.90, 
                                       xi=0.95,r=c(0.5,0.5),q_prior =as.numeric( df_q_andmu_posteriors[i,grep("q",colnames(df_q_andmu_posteriors)),]),
                                       delta_star=0.373)
  interim_ss<-rbind(interim_ss,interim_aux)
}

mean(interim_ss$treatment1)

#after interim analysis, we will need 461 per arm




interim_ss_rar<-tibble()
for(i in 1:nrow(df_sample_size)){
  interim_aux<-sample_size_calculation(alpha_prior =df_alpha_beta_posteriors$alpha_posterior[i],beta_prior = df_alpha_beta_posteriors$beta_posterior[i], eta=0.95, zeta=0.90, 
                                       xi=0.95,r=as.numeric(interim_allocation[i,]),q_prior =as.numeric( df_q_andmu_posteriors[i,grep("q",colnames(df_q_andmu_posteriors)),]),
                                       delta_star=0.373)
  interim_ss_rar<-rbind(interim_ss_rar,interim_aux)
}




#using unbalanced ratio--458 per arm


mean(interim_ss_rar$treatment1,na.rm = T)




#from  Whitehead et al. 2015
alpha<-1:10000
beta<-round(alpha/0.0204)
df_sample_size<-tibble()
df_posteriors<-tibble()
df_q_andmu_posteriors<-tibble()
df_alpha_beta_posteriors<-tibble()
df_D_posteriors<-tibble()
interim_allocation<-tibble()


for(i in 1:1000){
  alpha1<-sample(alpha,1)
  beta1<-beta[alpha1==alpha]
  aux<-sample_size_calculation(alpha_prior = alpha1,beta_prior = beta1, eta=0.95, zeta=0.90, 
                               xi=0.95,r=c(1/3,1/6,1/6,1/6,1/6),q_prior=c(10,2,2,2,2),delta_star=5)
  if(any(is.na(aux)))next
  y1_aux=rnorm(round(aux$treatment1/2),mean=2.8, sd=12.3)
  y2_aux=rnorm(round(aux$treatment2/2),mean=12.7, sd=14.1)
  y3_aux=rnorm(round(aux$treatment3/2),mean=14.3, sd=11.5)
  y4_aux=rnorm(round(aux$treatment4/2),mean=13.4, sd=14.4)
  y5_aux=rnorm(round(aux$treatment5/2),mean=17, sd=15)
  
  y=c(y1_aux,y2_aux,y3_aux,y4_aux,y5_aux)
  treatment_assignment<-c(rep(1,round(aux$treatment1/2)),rep(2,round(aux$treatment2/2)),rep(3,round(aux$treatment3/2)),rep(4,round(aux$treatment4/2)),
                          rep(5,round(aux$treatment5/2)))
  df=tibble(treatment_assignment,y)
  aux_post<-posterior_calculations(alpha_prior=alpha1,beta_prior=beta1,q_prior=c(10,2,2,2,2),
                                   N_treat = c(round(aux$treatment1/2),round(aux$treatment2/2),round(aux$treatment3/2),round(aux$treatment4/2),round(aux$treatment5/2)),
                                   y_treatment = df)
  df_q_andmu_posteriors<-rbind(df_q_andmu_posteriors,aux_post$q_andmu_posteriors)
  df_alpha_beta_posteriors<-rbind(df_alpha_beta_posteriors,aux_post$alpha_beta_params)
  df_D_posteriors<-rbind(df_D_posteriors,aux_post$D)
  df_sample_size<-rbind(df_sample_size,aux)
  treatment_differences<-get_treatment_difference(aux_post$q_andmu_posteriors,aux_post$D,aux_post$alpha_beta_params)
  new_r<-allocation_calculation(treatment_differences)
  interim_allocation<-rbind(interim_allocation,new_r)
  #df_posteriors<-rbind(df_posteriors,aux_post)
}
colMeans(df_sample_size)



interim_ss<-tibble()
for(i in 1:nrow(df_sample_size)){
  interim_aux<-sample_size_calculation(alpha_prior =df_alpha_beta_posteriors$alpha_posterior[i],beta_prior = df_alpha_beta_posteriors$beta_posterior[i], eta=0.95, zeta=0.90, 
                                       xi=0.95,r=c(1/3,1/6,1/6,1/6,1/6),q_prior =as.numeric( df_q_andmu_posteriors[i,grep("q",colnames(df_q_andmu_posteriors)),]),
                                       delta_star=5)
  interim_ss<-rbind(interim_ss,interim_aux)
}

colMeans(interim_ss)

#after interim analysis, we will need 461 per arm

colMeans(interim_allocation)


interim_ss_rar<-tibble()
for(i in 1:nrow(df_sample_size)){
  interim_aux<-sample_size_calculation(alpha_prior =df_alpha_beta_posteriors$alpha_posterior[i],beta_prior = df_alpha_beta_posteriors$beta_posterior[i], eta=0.95, zeta=0.90, 
                                       xi=0.95,r=as.numeric(interim_allocation[i,]),q_prior =as.numeric( df_q_andmu_posteriors[i,grep("q",colnames(df_q_andmu_posteriors)),]),
                                       delta_star=5)
  interim_ss_rar<-rbind(interim_ss_rar,interim_aux)
}

colMeans(interim_ss_rar,na.rm = T)




alpha<-1:10000
beta<-alpha
df_sample_size<-tibble()
df_posteriors<-tibble()
df_q_andmu_posteriors<-tibble()
df_alpha_beta_posteriors<-tibble()
df_D_posteriors<-tibble()
interim_allocation<-tibble()


for(i in 1:1000){
  alpha1<-sample(alpha,1)
  beta1<-beta[alpha1==alpha]
  aux<-sample_size_calculation(alpha_prior = alpha1,beta_prior = beta1, eta=0.95, zeta=0.90, 
                               xi=0.95,r=c(1/4,1/4,1/4,1/4),q_prior=c(2,2,2,2),delta_star=0.2)
  if(any(is.na(aux)))next
  y1_aux=rnorm(round(aux$treatment1/2),mean=0, sd=1)
  y2_aux=rnorm(round(aux$treatment2/2),mean=0.1, sd=1)
  y3_aux=rnorm(round(aux$treatment3/2),mean=0.2, sd=1)
  y4_aux=rnorm(round(aux$treatment4/2),mean=0.6, sd=1)
  
  y=c(y1_aux,y2_aux,y3_aux,y4_aux)
  treatment_assignment<-c(rep(1,round(aux$treatment1/2)),rep(2,round(aux$treatment2/2)),rep(3,round(aux$treatment3/2)),rep(4,round(aux$treatment4/2)))
  df=tibble(treatment_assignment,y)
  aux_post<-posterior_calculations(alpha_prior=alpha1,beta_prior=beta1,q_prior=c(2,2,2,2),mu_prior=c(0,0,0,0),
                                   N_treat = c(round(aux$treatment1/2),round(aux$treatment2/2),round(aux$treatment3/2),round(aux$treatment4/2)),
                                   y_treatment = df)
  df_q_andmu_posteriors<-rbind(df_q_andmu_posteriors,aux_post$q_andmu_posteriors)
  df_alpha_beta_posteriors<-rbind(df_alpha_beta_posteriors,aux_post$alpha_beta_params)
  df_D_posteriors<-rbind(df_D_posteriors,aux_post$D)
  df_sample_size<-rbind(df_sample_size,aux)
  treatment_differences<-get_treatment_difference(aux_post$q_andmu_posteriors,aux_post$D,aux_post$alpha_beta_params)
  new_r<-allocation_calculation(treatment_differences)
  interim_allocation<-rbind(interim_allocation,new_r)
  #df_posteriors<-rbind(df_posteriors,aux_post)
}
colMeans(df_sample_size)



interim_ss<-tibble()
for(i in 1:nrow(df_sample_size)){
  interim_aux<-sample_size_calculation(alpha_prior =df_alpha_beta_posteriors$alpha_posterior[i],beta_prior = df_alpha_beta_posteriors$beta_posterior[i], eta=0.95, zeta=0.90, 
                                       xi=0.95,r=c(1/4,1/4,1/4,1/4),q_prior =as.numeric( df_q_andmu_posteriors[i,grep("q",colnames(df_q_andmu_posteriors)),]),
                                       delta_star=0.2)
  interim_ss<-rbind(interim_ss,interim_aux)
}

colMeans(interim_ss)

#after interim analysis, we will need 461 per arm

colMeans(interim_allocation)


interim_ss_rar<-tibble()
for(i in 1:nrow(df_sample_size)){
  interim_aux<-sample_size_calculation(alpha_prior =df_alpha_beta_posteriors$alpha_posterior[i],beta_prior = df_alpha_beta_posteriors$beta_posterior[i], eta=0.95, zeta=0.90, 
                                       xi=0.95,r=as.numeric(interim_allocation[i,]),q_prior =as.numeric( df_q_andmu_posteriors[i,grep("q",colnames(df_q_andmu_posteriors)),]),
                                       delta_star=0.2)
  interim_ss_rar<-rbind(interim_ss_rar,interim_aux)
}

colMeans(interim_ss_rar,na.rm = T)






alpha<-1:10000
beta<-alpha
df_sample_size<-tibble()
df_posteriors<-tibble()
df_q_andmu_posteriors<-tibble()
df_alpha_beta_posteriors<-tibble()
df_D_posteriors<-tibble()
interim_allocation<-tibble()


for(i in 1:1000){
  alpha1<-sample(alpha,1)
  beta1<-beta[alpha1==alpha]
  aux<-sample_size_calculation(alpha_prior = alpha1,beta_prior = beta1, eta=0.95, zeta=0.90, 
                               xi=0.95,r=c(1/4,1/4,1/4,1/4),q_prior=c(2,2,2,2),delta_star=0.2)
  if(any(is.na(aux)))next
  y1_aux=rnorm(round(aux$treatment1/4),mean=0, sd=1)
  y2_aux=rnorm(round(aux$treatment2/4),mean=0.1, sd=1)
  y3_aux=rnorm(round(aux$treatment3/4),mean=0.2, sd=1)
  y4_aux=rnorm(round(aux$treatment4/4),mean=0.6, sd=1)
  
  y=c(y1_aux,y2_aux,y3_aux,y4_aux)
  treatment_assignment<-c(rep(1,round(aux$treatment1/4)),rep(2,round(aux$treatment2/4)),rep(3,round(aux$treatment3/4)),rep(4,round(aux$treatment4/4)))
  df=tibble(treatment_assignment,y)
  aux_post<-posterior_calculations(alpha_prior=alpha1,beta_prior=beta1,q_prior=c(2,2,2,2),mu_prior=c(0,0,0,0),
                                   N_treat = c(round(aux$treatment1/4),round(aux$treatment2/4),round(aux$treatment3/4),round(aux$treatment4/4)),
                                   y_treatment = df)
  df_q_andmu_posteriors<-rbind(df_q_andmu_posteriors,aux_post$q_andmu_posteriors)
  df_alpha_beta_posteriors<-rbind(df_alpha_beta_posteriors,aux_post$alpha_beta_params)
  df_D_posteriors<-rbind(df_D_posteriors,aux_post$D)
  df_sample_size<-rbind(df_sample_size,aux)
  treatment_differences<-get_treatment_difference(aux_post$q_andmu_posteriors,aux_post$D,aux_post$alpha_beta_params)
  new_r<-allocation_calculation(treatment_differences)
  interim_allocation<-rbind(interim_allocation,new_r)
  #df_posteriors<-rbind(df_posteriors,aux_post)
}
colMeans(df_sample_size)



interim_ss<-tibble()
for(i in 1:nrow(df_sample_size)){
  interim_aux<-sample_size_calculation(alpha_prior =df_alpha_beta_posteriors$alpha_posterior[i],beta_prior = df_alpha_beta_posteriors$beta_posterior[i], eta=0.95, zeta=0.90, 
                                       xi=0.95,r=c(1/4,1/4,1/4,1/4),q_prior =as.numeric( df_q_andmu_posteriors[i,grep("q",colnames(df_q_andmu_posteriors)),]),
                                       delta_star=0.2)
  interim_ss<-rbind(interim_ss,interim_aux)
}

colMeans(interim_ss)

#after interim analysis, we will need 461 per arm

colMeans(interim_allocation)


interim_ss_rar<-tibble()
for(i in 1:nrow(df_sample_size)){
  interim_aux<-sample_size_calculation(alpha_prior =df_alpha_beta_posteriors$alpha_posterior[i],beta_prior = df_alpha_beta_posteriors$beta_posterior[i], eta=0.95, zeta=0.90, 
                                       xi=0.95,r=as.numeric(interim_allocation[i,]),q_prior =as.numeric( df_q_andmu_posteriors[i,grep("q",colnames(df_q_andmu_posteriors)),]),
                                       delta_star=0.2)
  interim_ss_rar<-rbind(interim_ss_rar,interim_aux)
}

colMeans(interim_ss_rar,na.rm = T)





for(i in 1:1000){
  alpha1<-sample(alpha,1)
  beta1<-beta[alpha1==alpha]
  aux<-sample_size_calculation(alpha_prior = alpha1,beta_prior = beta1, eta=0.95, zeta=0.90, 
                               xi=0.95,r=c(1/4,1/4,1/4,1/4),q_prior=c(2,2,2,2),delta_star=0.2)
  if(any(is.na(aux)))next
  y1_aux=rnorm(round(aux$treatment1/8),mean=0, sd=1)
  y2_aux=rnorm(round(aux$treatment2/8),mean=0.1, sd=1)
  y3_aux=rnorm(round(aux$treatment3/8),mean=0.2, sd=1)
  y4_aux=rnorm(round(aux$treatment4/8),mean=0.6, sd=1)
  
  y=c(y1_aux,y2_aux,y3_aux,y4_aux)
  treatment_assignment<-c(rep(1,round(aux$treatment1/8)),rep(2,round(aux$treatment2/8)),rep(3,round(aux$treatment3/8)),rep(4,round(aux$treatment4/8)))
  df=tibble(treatment_assignment,y)
  aux_post<-posterior_calculations(alpha_prior=alpha1,beta_prior=beta1,q_prior=c(2,2,2,2),mu_prior=c(0,0,0,0),
                                   N_treat = c(round(aux$treatment1/2),round(aux$treatment2/8),round(aux$treatment3/8),round(aux$treatment4/8)),
                                   y_treatment = df)
  df_q_andmu_posteriors<-rbind(df_q_andmu_posteriors,aux_post$q_andmu_posteriors)
  df_alpha_beta_posteriors<-rbind(df_alpha_beta_posteriors,aux_post$alpha_beta_params)
  df_D_posteriors<-rbind(df_D_posteriors,aux_post$D)
  df_sample_size<-rbind(df_sample_size,aux)
  treatment_differences<-get_treatment_difference(aux_post$q_andmu_posteriors,aux_post$D,aux_post$alpha_beta_params)
  new_r<-allocation_calculation(treatment_differences)
  interim_allocation<-rbind(interim_allocation,new_r)
  #df_posteriors<-rbind(df_posteriors,aux_post)
}
colMeans(df_sample_size)



interim_ss<-tibble()
for(i in 1:nrow(df_sample_size)){
  interim_aux<-sample_size_calculation(alpha_prior =df_alpha_beta_posteriors$alpha_posterior[i],beta_prior = df_alpha_beta_posteriors$beta_posterior[i], eta=0.95, zeta=0.90, 
                                       xi=0.95,r=c(1/4,1/4,1/4,1/4),q_prior =as.numeric( df_q_andmu_posteriors[i,grep("q",colnames(df_q_andmu_posteriors)),]),
                                       delta_star=0.2)
  interim_ss<-rbind(interim_ss,interim_aux)
}

colMeans(interim_ss)

#after interim analysis, we will need 461 per arm

colMeans(interim_allocation)


interim_ss_rar<-tibble()
for(i in 1:nrow(df_sample_size)){
  interim_aux<-sample_size_calculation(alpha_prior =df_alpha_beta_posteriors$alpha_posterior[i],beta_prior = df_alpha_beta_posteriors$beta_posterior[i], eta=0.95, zeta=0.90, 
                                       xi=0.95,r=as.numeric(interim_allocation[i,]),q_prior =as.numeric( df_q_andmu_posteriors[i,grep("q",colnames(df_q_andmu_posteriors)),]),
                                       delta_star=0.2)
  interim_ss_rar<-rbind(interim_ss_rar,interim_aux)
}

colMeans(interim_ss_rar,na.rm = T)

