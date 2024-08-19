# Helper functions for simulating xi


generate_auxiliary_data<-function(percent_outcome_available=50, sample_sizes){
  length1<-round(sample_sizes$treatment1*percent_outcome_available/100)
  y1_aux=rnorm(length1,mean=0, sd=1)
  length2<-round(sample_sizes$treatment2*percent_outcome_available/100)
  y2_aux=rnorm(length2,mean=0.1, sd=1)
  length3<-round(sample_sizes$treatment3*percent_outcome_available/100)
  y3_aux=rnorm(length3,mean=0.2, sd=1)
  length4<-round(sample_sizes$treatment4*percent_outcome_available/100)
  y4_aux=rnorm(length4,mean=0.6, sd=1)
  
  y=c(y1_aux,y2_aux,y3_aux,y4_aux)
  treatment_assignment<-c(rep(1,length1),rep(2,length2),
                          rep(3,length3),rep(4,length4))
  df=tibble(treatment_assignment,y)
  return(df)
}



get_number_treated<-function(df){
  treatment_lengths<-df%>%count(treatment_assignment)%>%select(n)
  treatment_lengths<-treatment_lengths$n
  return(treatment_lengths)
}

get_interim_ss<-function(alpha_beta_posteriors,q_and_mu_posteriors, 
                         interim_allocation, method="ER",xi_value=0.95,eta_value=0.95,
                         zeta_value=0.9,delta_value=0.3){
  outcome<-tibble()
  if(method=="ER"){
    for(i in 1:nrow(alpha_beta_posteriors)){
      interim_aux<-sample_size_calculation(alpha_prior= alpha_beta_posteriors$alpha_posterior[i],
                                           beta_prior = alpha_beta_posteriors$beta_posterior[i],
                                           eta=eta_value, zeta=zeta_value, 
                                           xi=xi_value,r=c(1/4,1/4,1/4,1/4),
                                           q_prior =as.numeric(q_and_mu_posteriors[i,
                                                                                   grep("q",colnames(q_and_mu_posteriors)),]),
                                           delta_star=delta_value)
      interim_aux$xi<-xi_value
      outcome<-outcome%>%bind_rows(interim_aux)
      
    }
  }else if(method=="RAR"){
    for(i in 1:nrow(alpha_beta_posteriors)){
      interim_aux<-sample_size_calculation(alpha_prior= alpha_beta_posteriors$alpha_posterior[i],
                                           beta_prior = alpha_beta_posteriors$beta_posterior[i],
                                           eta=eta_value, zeta=zeta_value, 
                                           xi=xi_value,r=interim_allocation,
                                           q_prior =as.numeric(q_and_mu_posteriors[i,
                                                                                   grep("q",colnames(q_and_mu_posteriors)),]),
                                           delta_star=delta_value)
      interim_aux$xi<-xi_value
      outcome<-outcome%>%bind_rows(interim_aux)
    }
  }
  
  return(outcome)
}


generate_auxiliary_data_two_analyses<-function(percent_outcome_available=50,
                                               percent_available_at_first_interim=25, 
                                               sample_sizes, interim_sample_sizes){
  
  #Scenario 4--25% then 50%
  length1<-max(0,round((interim_sample_sizes$treatment1)*percent_outcome_available/100-
                         sample_sizes$treatment1*percent_available_at_first_interim/100))
  y1_aux=rnorm(length1,mean=0, sd=1)
  length2<-max(0,round((interim_sample_sizes$treatment2)*percent_outcome_available/100-
                         sample_sizes$treatment2*percent_available_at_first_interim/100))
  y2_aux=rnorm(length2,mean=0.1, sd=1)
  length3<-max(0,round((interim_sample_sizes$treatment3)*percent_outcome_available/100-
                         sample_sizes$treatment3*percent_available_at_first_interim/100))
  y3_aux=rnorm(length3,mean=0.2, sd=1)
  length4<-max(0,round((interim_sample_sizes$treatment4)*percent_outcome_available/100-
                         sample_sizes$treatment4*percent_available_at_first_interim/100))
  y4_aux=rnorm(length4,mean=0.6, sd=1)
  
  y=c(y1_aux,y2_aux,y3_aux,y4_aux)
  treatment_assignment<-c(rep(1,length1),rep(2,length2),rep(3,length3),rep(4,length4))
  df=tibble(treatment_assignment,y)
  return(df)
}
