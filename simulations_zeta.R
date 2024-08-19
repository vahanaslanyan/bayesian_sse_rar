set.seed(219)
source("bayesian_ssr.R")
source("helper_xi.R")

alpha<-runif(5000,0,20)
beta<-alpha

df_sample_size<-tibble()
df_q_andmu_posteriors_50<-tibble()
df_alpha_beta_posteriors_50<-tibble()
df_D_posteriors_50<-tibble()
interim_allocation_50<-tibble()
df_q_andmu_posteriors_25<-tibble()
df_alpha_beta_posteriors_25<-tibble()
df_D_posteriors_25<-tibble()
interim_allocation_25<-tibble()
df_q_andmu_posteriors_125<-tibble()
df_alpha_beta_posteriors_125<-tibble()
df_D_posteriors_125<-tibble()
interim_allocation_125<-tibble()
interim_ss50<-tibble()
interim_ss_rar50<-tibble()
interim_ss25<-tibble()
interim_ss_rar25<-tibble()
interim_ss125<-tibble()
interim_ss_rar125<-tibble()

zeta_stars<-seq(0.70,0.99,by=0.02)
for (zeta_value in zeta_stars){
  for(i in 1:length(alpha)){
    aux<-sample_size_calculation(alpha_prior = alpha[i],beta_prior = beta[i], 
                                 eta=0.95, zeta=zeta_value, xi=0.95,
                                 r=c(1/4,1/4,1/4,1/4),q_prior=c(1,1,1,1),
                                 delta_star=0.3)
    aux$zeta<-zeta_value
    if(any(is.na(aux)))next
    df_sample_size<-df_sample_size%>%bind_rows(aux)
    #Scenario 1--50% outcome available
    df_50<-generate_auxiliary_data(percent_outcome_available = 50,
                                   sample_sizes = aux)
    treatment_lengths<-get_number_treated(df=df_50)
    #get posteriors based on generated data
    aux_post_50<-posterior_calculations(alpha_prior=alpha[i],beta_prior=beta[i],
                                        q_prior=c(1,1,1,1),mu_prior=c(0,0,0,0),
                                        N_treat = treatment_lengths,
                                        y_treatment = df_50)
    
    df_q_andmu_posteriors_50<-df_q_andmu_posteriors_50%>%
      bind_rows(aux_post_50$q_andmu_posteriors)
    df_alpha_beta_posteriors_50<-df_alpha_beta_posteriors_50%>%
      bind_rows(aux_post_50$alpha_beta_params)
    df_D_posteriors_50<-df_D_posteriors_50%>%bind_rows(aux_post_50$D)
    #get treatment differences and allocation ratios
    treatment_differences_50<-get_treatment_difference(aux_post_50$q_andmu_posteriors,
                                                       aux_post_50$D,
                                                       aux_post_50$alpha_beta_params)
    new_r<-allocation_calculation(treatment_differences_50)
    names(new_r)<-c("Trt1","Trt2","Trt3","Trt4")
    interim_allocation_50<-interim_allocation_50%>%bind_rows(new_r)
    
    #Scenario 2--25% outcome available
    df_25<-generate_auxiliary_data(percent_outcome_available = 25,
                                   sample_sizes = aux)
    treatment_lengths<-get_number_treated(df=df_25)
    #get posteriors based on generated data
    aux_post_25<-posterior_calculations(alpha_prior=alpha[i],beta_prior=beta[i],
                                        q_prior=c(1,1,1,1),mu_prior=c(0,0,0,0),
                                        N_treat = treatment_lengths,
                                        y_treatment = df_25)
    
    df_q_andmu_posteriors_25<-df_q_andmu_posteriors_25%>%
      bind_rows(aux_post_25$q_andmu_posteriors)
    df_alpha_beta_posteriors_25<-df_alpha_beta_posteriors_25%>%
      bind_rows(aux_post_25$alpha_beta_params)
    df_D_posteriors_25<-df_D_posteriors_25%>%bind_rows(aux_post_25$D)
    
    treatment_differences_25<-get_treatment_difference(aux_post_25$q_andmu_posteriors,
                                                       aux_post_25$D,
                                                       aux_post_25$alpha_beta_params)
    new_r<-allocation_calculation(treatment_differences_25)
    names(new_r)<-c("Trt1","Trt2","Trt3","Trt4")
    interim_allocation_25<-interim_allocation_25%>%bind_rows(new_r)
    #Scenario 3--12.5% outcome available
    df_125<-generate_auxiliary_data(percent_outcome_available = 12.5,
                                    sample_sizes = aux)
    treatment_lengths<-get_number_treated(df=df_125)
    aux_post_125<-posterior_calculations(alpha_prior=alpha[i],beta_prior=beta[i],q_prior=c(1,1,1,1),mu_prior=c(0,0,0,0),
                                         N_treat = treatment_lengths,
                                         y_treatment = df_125)
    df_q_andmu_posteriors_125<-df_q_andmu_posteriors_125%>%
      bind_rows(aux_post_125$q_andmu_posteriors)
    df_alpha_beta_posteriors_125<-df_alpha_beta_posteriors_125%>%
      bind_rows(aux_post_125$alpha_beta_params)
    df_D_posteriors_125<-df_D_posteriors_125%>%bind_rows(aux_post_125$D)
    
    treatment_differences_125<-get_treatment_difference(aux_post_125$q_andmu_posteriors,
                                                        aux_post_125$D,
                                                        aux_post_125$alpha_beta_params)
    new_r<-allocation_calculation(treatment_differences_125)
    names(new_r)<-c("Trt1","Trt2","Trt3","Trt4")
    interim_allocation_125<-rbind(interim_allocation_125,new_r)
    
  }
    
  interim_ss50<-get_interim_ss(alpha_beta_posteriors = df_alpha_beta_posteriors_50,
                               q_and_mu_posteriors = df_q_andmu_posteriors_50,
                               interim_allocation = as.numeric(interim_allocation_50),
                               method ="ER",zeta_value = zeta_value)
  interim_ss_rar50<-get_interim_ss(alpha_beta_posteriors = df_alpha_beta_posteriors_50,
                                   q_and_mu_posteriors = df_q_andmu_posteriors_50,
                                   interim_allocation = as.numeric(interim_allocation_50),
                                   method ="RAR",zeta_value = zeta_value)
  interim_ss25<-get_interim_ss(alpha_beta_posteriors = df_alpha_beta_posteriors_25,
                               q_and_mu_posteriors = df_q_andmu_posteriors_25,
                               interim_allocation = as.numeric(interim_allocation_25),
                               method ="ER",zeta_value = zeta_value)
  interim_ss_rar25<-get_interim_ss(alpha_beta_posteriors = df_alpha_beta_posteriors_25,
                                   q_and_mu_posteriors = df_q_andmu_posteriors_25,
                                   interim_allocation = as.numeric(interim_allocation_25),
                                   method ="RAR",zeta_value = zeta_value)
  interim_ss125<-get_interim_ss(alpha_beta_posteriors = df_alpha_beta_posteriors_125,
                                q_and_mu_posteriors = df_q_andmu_posteriors_125,
                                interim_allocation = as.numeric(interim_allocation_125),
                                method ="ER",zeta_value =zeta_value)
  interim_ss_rar125<-get_interim_ss(alpha_beta_posteriors = df_alpha_beta_posteriors_125,
                                    q_and_mu_posteriors = df_q_andmu_posteriors_125,
                                    interim_allocation = as.numeric(interim_allocation_125),
                                    method ="RAR",zeta_value = zeta_value)
  
  df_q_andmu_posteriors_50<-tibble()
  df_alpha_beta_posteriors_50<-tibble()
  df_D_posteriors_50<-tibble()
  interim_allocation_50<-tibble()
  
  df_q_andmu_posteriors_25<-tibble()
  df_alpha_beta_posteriors_25<-tibble()
  df_D_posteriors_25<-tibble()
  interim_allocation_25<-tibble()
  
  df_q_andmu_posteriors_125<-tibble()
  df_alpha_beta_posteriors_125<-tibble()
  df_D_posteriors_125<-tibble()
  interim_allocation_125<-tibble()
}

write_csv(df_sample_size,"changing_zeta_no_interim_analysis.csv")
write_csv(interim_ss50,"changing_zeta_no_rar_interim_analysis50.csv")
write_csv(interim_ss_rar50,"changing_zeta_rar_interim_analysis50.csv")
write_csv(interim_ss25,"changing_zeta_no_rar_interim_analysis25.csv")
write_csv(interim_ss_rar25,"changing_zeta_rar_interim_analysis25.csv")
write_csv(interim_ss125,"changing_zeta_no_rar_interim_analysis125.csv")
write_csv(interim_ss_rar125,"changing_zeta_rar_interim_analysis125.csv")

interim_ss_rar50<-tibble()

for (zeta_value in zeta_stars){
  interim_ss_rar25<-tibble()
  df_sample_size<-tibble()
  df_q_andmu_posteriors_50<-tibble()
  df_alpha_beta_posteriors_50<-tibble()
  df_D_posteriors_50<-tibble()
  interim_allocation_50<-tibble()
  df_q_andmu_posteriors_25<-tibble()
  df_alpha_beta_posteriors_25<-tibble()
  df_D_posteriors_25<-tibble()
  interim_allocation_25<-tibble()
  for(i in 1:length(alpha)){
    aux<-sample_size_calculation(alpha_prior = alpha[i],beta_prior = beta[i],
                                 eta=0.95, zeta=zeta_value, xi=0.95,
                                 r=c(1/4,1/4,1/4,1/4),
                                 q_prior=c(1,1,1,1),
                                 delta_star=0.3)
    if(any(is.na(aux)))next
    df_sample_size<-rbind(df_sample_size,aux)
    #Scenario 4--25% and 50%
    #generate data from 25% outcomes
    df_25<-generate_auxiliary_data(percent_outcome_available = 25,
                                   sample_sizes = aux)
    treatment_lengths<-get_number_treated(df=df_25)
    #get posterior distributions given the data
    aux_post_25<-posterior_calculations(alpha_prior=alpha[i],beta_prior=beta[i],
                                        q_prior=c(1,1,1,1),mu_prior=c(0,0,0,0),
                                        N_treat = treatment_lengths,
                                        y_treatment = df_25)
    
    df_q_andmu_posteriors_25<-aux_post_25$q_andmu_posteriors
    df_alpha_beta_posteriors_25<-aux_post_25$alpha_beta_params
    df_D_posteriors_25<-aux_post_25$D
    #get the new allocation ratios
    treatment_differences_25<-get_treatment_difference(aux_post_25$q_andmu_posteriors,
                                                       aux_post_25$D,
                                                       aux_post_25$alpha_beta_params)
    new_r<-allocation_calculation(treatment_differences_25)
    names(new_r)<-c("Trt1","Trt2","Trt3","Trt4")
    interim_allocation_25<-interim_allocation_25%>%bind_rows(new_r)
    #calculate sample size again, using posteriors as priors
    
    interim_aux<-sample_size_calculation(alpha_prior=df_alpha_beta_posteriors_25$alpha_posterior,
                                         beta_prior = df_alpha_beta_posteriors_25$beta_posterior,
                                         eta=0.95, zeta=zeta_value, 
                                         xi=0.95,
                                         r=as.numeric(new_r),
                                         q_prior =as.numeric( df_q_andmu_posteriors_25[,grep("q",colnames(df_q_andmu_posteriors_25)),]),
                                         delta_star=0.3)
    if(any(is.na(interim_aux)))next
    interim_aux<-interim_aux+aux/4
    interim_ss_rar25<-interim_ss_rar25%>%bind_rows(interim_aux)
    #generate data for the second interim analysis
    
    df_50=generate_auxiliary_data_two_analyses(percent_outcome_available = 50,
                                               percent_available_at_first_interim = 25,
                                               interim_sample_sizes = interim_aux,
                                               sample_sizes = aux)
    df<-df_25%>%bind_rows(df_50)
    #calculate posteriors with the newly generated data
    aux_post_50<-posterior_calculations(alpha_prior =df_alpha_beta_posteriors_25$alpha_posterior ,
                                        beta_prior = df_alpha_beta_posteriors_25$beta_posterior, 
                                        q_prior=as.numeric( df_q_andmu_posteriors_25[,grep("q",colnames(df_q_andmu_posteriors_25)),]),
                                        mu_prior=as.numeric( df_q_andmu_posteriors_25[,grep("mu",colnames(df_q_andmu_posteriors_25)),]),
                                        N_treat= c(round(interim_aux$treatment1/2),
                                                   round(interim_aux$treatment2/2),
                                                   round(interim_aux$treatment3/2),
                                                   round(interim_aux$treatment4/2)),
                                        y_treatment = df)
    
    #sample size at 50% outcome available
    sample_sizes_before50<-interim_aux/2
    df_q_andmu_posteriors_50<-aux_post_50$q_andmu_posteriors
    df_alpha_beta_posteriors_50<-aux_post_50$alpha_beta_params
    df_D_posteriors_50<-aux_post_50$D
    treatment_differences_50<-get_treatment_difference(aux_post_50$q_andmu_posteriors,
                                                       aux_post_50$D,
                                                       aux_post_50$alpha_beta_params)
    new_r<-allocation_calculation(treatment_differences_50)
    names(new_r)<-c("Trt1","Trt2","Trt3","Trt4")
    interim_allocation_50<-interim_allocation_50%>%bind_rows(new_r)
    interim_aux<-sample_size_calculation(alpha_prior=df_alpha_beta_posteriors_50$alpha_posterior,
                                         beta_prior = df_alpha_beta_posteriors_50$beta_posterior,
                                         eta=0.95, zeta=zeta_value, 
                                         xi=0.95,r=as.numeric(new_r),
                                         q_prior =as.numeric( df_q_andmu_posteriors_50[,grep("q",colnames(df_q_andmu_posteriors_50)),]),
                                         delta_star=0.3)
    if(any(is.na(interim_aux)))next
    interim_aux<-interim_aux+sample_sizes_before50
    interim_aux$zeta<-zeta_value
    interim_ss_rar50<-interim_ss_rar50%>%bind_rows(interim_aux)
    
  }
  
}

write_csv(interim_ss_rar50,"changing_zeta_rar_interim_analyses25_50.csv")

interim_ss_50<-tibble()

for(zeta_value in zeta_stars){
  interim_ss_25<-tibble()
  df_sample_size<-tibble()
  df_q_andmu_posteriors_50<-tibble()
  df_alpha_beta_posteriors_50<-tibble()
  df_D_posteriors_50<-tibble()
  interim_allocation_50<-tibble()
  df_q_andmu_posteriors_25<-tibble()
  df_alpha_beta_posteriors_25<-tibble()
  df_D_posteriors_25<-tibble()
  interim_allocation_25<-tibble()
  for(i in 1:length(alpha)){
    aux<-sample_size_calculation(alpha_prior = alpha[i],beta_prior = beta[i],
                                 eta=0.95, zeta=zeta_value, 
                                 xi=0.95,
                                 r=c(1/4,1/4,1/4,1/4),
                                 q_prior=c(1,1,1,1),
                                 delta_star=0.3)
    if(any(is.na(aux)))next
    df_sample_size<-df_sample_size%>%bind_rows(aux)
    #Scenario 4
    df_25<-generate_auxiliary_data(percent_outcome_available = 25,
                                   sample_sizes = aux)
    treatment_lengths<-get_number_treated(df=df_25)
    aux_post_25<-posterior_calculations(alpha_prior=alpha[i],beta_prior=beta[i],
                                        q_prior=c(1,1,1,1),mu_prior=c(0,0,0,0),
                                        N_treat = treatment_lengths,
                                        y_treatment = df_25)
    
    df_q_andmu_posteriors_25<-aux_post_25$q_andmu_posteriors
    df_alpha_beta_posteriors_25<-aux_post_25$alpha_beta_params
    df_D_posteriors_25<-aux_post_25$D
    interim_aux<-sample_size_calculation(alpha_prior =df_alpha_beta_posteriors_25$alpha_posterior,
                                         beta_prior = df_alpha_beta_posteriors_25$beta_posterior,
                                         eta=0.95, zeta=zeta_value, 
                                         xi=0.95,r=c(1/4,1/4,1/4,1/4),
                                         q_prior =as.numeric( df_q_andmu_posteriors_25[,grep("q",colnames(df_q_andmu_posteriors_25)),]),
                                         delta_star=0.3)
    if(any(is.na(interim_aux)))next
    interim_aux<-interim_aux+aux/4
    interim_ss_rar25<-interim_ss_rar25%>%bind_rows(interim_aux)
    
    
    #Scenario 4
    df_50=generate_auxiliary_data_two_analyses(percent_outcome_available = 50,
                                               percent_available_at_first_interim = 25,
                                               interim_sample_sizes = interim_aux,
                                               sample_sizes = aux)
    df<-df_25%>%bind_rows(df_50)
    aux_post_50<-posterior_calculations(alpha_prior=df_alpha_beta_posteriors_25$alpha_posterior,
                                        beta_prior = df_alpha_beta_posteriors_25$beta_posterior, 
                                        q_prior=as.numeric( df_q_andmu_posteriors_25[,grep("q",colnames(df_q_andmu_posteriors_25)),]),
                                        mu_prior=as.numeric( df_q_andmu_posteriors_25[,grep("mu",colnames(df_q_andmu_posteriors_25)),]),
                                        N_treat= c(round(interim_aux$treatment1/2),
                                                   round(interim_aux$treatment2/2),
                                                   round(interim_aux$treatment3/2),
                                                   round(interim_aux$treatment4/2)),
                                        y_treatment = df)
    sample_sizes_before50<-interim_aux/2
    df_q_andmu_posteriors_50<-aux_post_50$q_andmu_posteriors
    df_alpha_beta_posteriors_50<-aux_post_50$alpha_beta_params
    df_D_posteriors_50<-aux_post_50$D
    interim_aux<-sample_size_calculation(alpha_prior =df_alpha_beta_posteriors_50$alpha_posterior,
                                         beta_prior = df_alpha_beta_posteriors_50$beta_posterior,
                                         eta=0.95, zeta=zeta_value, 
                                         xi=0.95,r=c(1/4,1/4,1/4,1/4),
                                         q_prior =as.numeric( df_q_andmu_posteriors_50[,grep("q",colnames(df_q_andmu_posteriors_50)),]),
                                         delta_star=0.3)
    if(any(is.na(interim_aux)))next
    interim_aux<-interim_aux+sample_sizes_before50
    interim_aux$zeta<-zeta_value
    interim_ss_50<-interim_ss_50%>%bind_rows(interim_aux)
    
  }
}

write_csv(interim_ss_50,"changing_zeta_no_rar_interim_analyses25_50.csv")


interim_ss_rar75<-tibble()

for(zeta_value in zeta_stars){
  interim_ss_rar25<-tibble()
  df_sample_size<-tibble()
  df_q_andmu_posteriors_75<-tibble()
  df_alpha_beta_posteriors_75<-tibble()
  df_D_posteriors_75<-tibble()
  interim_allocation_75<-tibble()
  df_q_andmu_posteriors_25<-tibble()
  df_alpha_beta_posteriors_25<-tibble()
  df_D_posteriors_25<-tibble()
  interim_allocation_25<-tibble()
  
  
  
  for(i in 1:length(alpha)){
    aux<-sample_size_calculation(alpha_prior = alpha[i],beta_prior = beta[i],
                                 eta=0.95, zeta=zeta_value,
                                 xi=0.95,r=c(1/4,1/4,1/4,1/4),
                                 q_prior=c(1,1,1,1),delta_star=0.3)
    if(any(is.na(aux)))next
    df_sample_size<-rbind(df_sample_size,aux)
    #Scenario 5--25% and 75%
    df_25<-generate_auxiliary_data(percent_outcome_available = 25,
                                   sample_sizes = aux)
    treatment_lengths<-get_number_treated(df=df_25)
    
    aux_post_25<-posterior_calculations(alpha_prior=alpha[i],beta_prior=beta[i],
                                        q_prior=c(1,1,1,1),mu_prior=c(0,0,0,0),
                                        N_treat = treatment_lengths,
                                        y_treatment = df_25)
    
    df_q_andmu_posteriors_25<-aux_post_25$q_andmu_posteriors
    df_alpha_beta_posteriors_25<-aux_post_25$alpha_beta_params
    df_D_posteriors_25<-aux_post_25$D
    
    treatment_differences_25<-get_treatment_difference(aux_post_25$q_andmu_posteriors,
                                                       aux_post_25$D,
                                                       aux_post_25$alpha_beta_params)
    new_r<-allocation_calculation(treatment_differences_25)
    names(new_r)<-c("Trt1","Trt2","Trt3","Trt4")
    interim_allocation_25<-interim_allocation_25%>%bind_rows(new_r)
    
    interim_aux<-sample_size_calculation(alpha_prior =df_alpha_beta_posteriors_25$alpha_posterior,
                                         beta_prior = df_alpha_beta_posteriors_25$beta_posterior,
                                         eta=0.95, zeta=zeta_value, 
                                         xi=0.95,r=as.numeric(new_r),
                                         q_prior =as.numeric( df_q_andmu_posteriors_25[,grep("q",colnames(df_q_andmu_posteriors_25)),]),
                                         delta_star=0.3)
    if(any(is.na(interim_aux)))next
    interim_aux<-interim_aux+aux/4
    interim_ss_rar25<-interim_ss_rar25%>%bind_rows(interim_aux)
    
    
    #Scenario 5
    df_75=generate_auxiliary_data_two_analyses(percent_outcome_available = 75,
                                               percent_available_at_first_interim = 25,
                                               interim_sample_sizes = interim_aux,
                                               sample_sizes = aux)
    df<-df_25%>%bind_rows(df_75)
    aux_post_75<-posterior_calculations(alpha_prior=df_alpha_beta_posteriors_25$alpha_posterior,beta_prior=df_alpha_beta_posteriors_25$beta_posterior,q_prior=as.numeric( df_q_andmu_posteriors_25[,grep("q",colnames(df_q_andmu_posteriors_25)),]),
                                        mu_prior=as.numeric( df_q_andmu_posteriors_25[,grep("mu",colnames(df_q_andmu_posteriors_25)),]),
                                        N_treat = c(round(interim_aux$treatment1*0.75),
                                                    round(interim_aux$treatment2*0.75),
                                                    round(interim_aux$treatment3*0.75),
                                                    round(interim_aux$treatment4*0.75)),
                                        y_treatment = df)
    sample_sizes_before75<-interim_aux*0.75
    df_q_andmu_posteriors_75<-aux_post_75$q_andmu_posteriors
    df_alpha_beta_posteriors_75<-aux_post_75$alpha_beta_params
    df_D_posteriors_75<-aux_post_75$D
    treatment_differences_75<-get_treatment_difference(aux_post_75$q_andmu_posteriors,
                                                       aux_post_75$D,
                                                       aux_post_75$alpha_beta_params)
    new_r<-allocation_calculation(treatment_differences_75)
    names(new_r)<-c("Trt1","Trt2","Trt3","Trt4")
    interim_allocation_75<-interim_allocation_75%>%bind_rows(new_r)
    
    interim_aux<-sample_size_calculation(alpha_prior=df_alpha_beta_posteriors_75$alpha_posterior,
                                         beta_prior=df_alpha_beta_posteriors_75$beta_posterior,
                                         eta=0.95, zeta=zeta_value, 
                                         xi=0.95,r=as.numeric(new_r),
                                         q_prior=as.numeric( df_q_andmu_posteriors_75[,grep("q",colnames(df_q_andmu_posteriors_75)),]),
                                         delta_star=0.3)
    
    if(any(is.na(interim_aux)))next
    interim_aux<-interim_aux+sample_sizes_before75
    interim_aux$zeta<-zeta_value
    interim_ss_rar75<-interim_ss_rar75%>%bind_rows(interim_aux)
    
  }
}

write_csv(interim_ss_rar75,"changing_zeta_rar_interim_analyses25_75.csv")

interim_ss_75<-tibble()

for(zeta_value in zeta_stars){

  interim_ss_25<-tibble()
  df_sample_size<-tibble()
  df_q_andmu_posteriors_75<-tibble()
  df_alpha_beta_posteriors_75<-tibble()
  df_D_posteriors_75<-tibble()
  interim_allocation_75<-tibble()
  df_q_andmu_posteriors_25<-tibble()
  df_alpha_beta_posteriors_25<-tibble()
  df_D_posteriors_25<-tibble()
  interim_allocation_25<-tibble()
  for(i in 1:length(alpha)){
    aux<-sample_size_calculation(alpha_prior = alpha[i],beta_prior = beta[i],
                                 eta=0.95, zeta=zeta_value, xi=0.95,
                                 r=c(1/4,1/4,1/4,1/4),q_prior=c(1,1,1,1),
                                 delta_star=0.3)
    if(any(is.na(aux)))next
    df_sample_size<-df_sample_size%>%bind_rows(aux)
    #Scenario 4
    df_25<-generate_auxiliary_data(percent_outcome_available = 25,
                                   sample_sizes = aux)
    treatment_lengths<-get_number_treated(df=df_25)
    aux_post_25<-posterior_calculations(alpha_prior=alpha[i],beta_prior=beta[i],
                                        q_prior=c(1,1,1,1),mu_prior=c(0,0,0,0),
                                        N_treat = treatment_lengths,
                                        y_treatment = df_25)
    
    df_q_andmu_posteriors_25<-aux_post_25$q_andmu_posteriors
    df_alpha_beta_posteriors_25<-aux_post_25$alpha_beta_params
    df_D_posteriors_25<-aux_post_25$D
    interim_aux<-sample_size_calculation(alpha_prior =df_alpha_beta_posteriors_25$alpha_posterior,
                                         beta_prior = df_alpha_beta_posteriors_25$beta_posterior,
                                         eta=0.95, zeta=zeta_value, 
                                         xi=0.95,r=c(1/4,1/4,1/4,1/4),
                                         q_prior =as.numeric( df_q_andmu_posteriors_25[,grep("q",colnames(df_q_andmu_posteriors_25)),]),
                                         delta_star=0.3)
    if(any(is.na(interim_aux)))next
    interim_aux<-interim_aux+aux/4
    interim_ss_rar25<-interim_ss_rar25%>%bind_rows(interim_aux)
    
    
    #Scenario 5
    df_75=generate_auxiliary_data_two_analyses(percent_outcome_available = 75,
                                               percent_available_at_first_interim = 25,
                                               interim_sample_sizes = interim_aux,
                                               sample_sizes = aux)
    df<-df_25%>%bind_rows(df_75)
    aux_post_75<-posterior_calculations(alpha_prior=df_alpha_beta_posteriors_25$alpha_posterior,
                                        beta_prior=df_alpha_beta_posteriors_25$beta_posterior,
                                        q_prior=as.numeric( df_q_andmu_posteriors_25[,grep("q",colnames(df_q_andmu_posteriors_25)),]),
                                        mu_prior=as.numeric( df_q_andmu_posteriors_25[,grep("mu",colnames(df_q_andmu_posteriors_25)),]),
                                        N_treat = c(round(interim_aux$treatment1*0.75),
                                                    round(interim_aux$treatment2*0.75),
                                                    round(interim_aux$treatment3*0.75),
                                                    round(interim_aux$treatment4*0.75)),
                                        y_treatment = df)
    sample_sizes_before75<-interim_aux*0.75
    df_q_andmu_posteriors_75<-aux_post_75$q_andmu_posteriors
    df_alpha_beta_posteriors_75<-aux_post_75$alpha_beta_params
    df_D_posteriors_75<-aux_post_75$D
    
    interim_aux<-sample_size_calculation(alpha_prior=df_alpha_beta_posteriors_75$alpha_posterior,
                                         beta_prior = df_alpha_beta_posteriors_75$beta_posterior,
                                         eta=0.95, zeta=zeta_value, 
                                         xi=0.95,r=c(1/4,1/4,1/4,1/4),
                                         q_prior=as.numeric( df_q_andmu_posteriors_75[,grep("q",colnames(df_q_andmu_posteriors_75)),]),
                                         delta_star=0.3)
    
    if(any(is.na(interim_aux)))next
    interim_aux<-interim_aux+sample_sizes_before75
    interim_aux$zeta<-zeta_value
    interim_ss_75<-interim_ss_75%>%bind_rows(interim_aux)
    
  }
}

write_csv(interim_ss_75,"changing_zeta_no_rar_interim_analyses25_75.csv")
