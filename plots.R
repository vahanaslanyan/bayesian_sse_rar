library(ggplot2)
library(tidyverse)
library(ggpubr)
library(ggsci)
setwd("/Users/vahanaslanyan/Documents/bayesian_sse_rar/bayesian_sse_rar/misc files/")
df1<-read_csv("changing_delta_no_interim_analysis.csv")
df1$Total_participants<-df1$treatment1+df1$treatment2+
  df1$treatment3+df1$treatment4
df1%>% 
  group_by(delta_star) %>% 
  summarize(total_participants = mean(Total_participants))->summary_delta_no_interim
summary_delta_no_interim_plot<-summary_delta_no_interim
summary_delta_no_interim$method<-"No"
colnames(summary_delta_no_interim)[1]<-"delta"

ggplot(summary_delta_no_interim_plot,aes(x=delta_star,y=total_participants))+geom_point(size=1.5)+geom_line(size=1.5)+theme_bw()



df2<-read_csv("changing_delta_no_rar_interim_analysis50.csv")
df2$Total_participants<-df2$treatment1+df2$treatment2+
  df2$treatment3+df2$treatment4
df2%>% 
  group_by(delta) %>% 
  summarize(total_participants = mean(Total_participants))->summary_delta_no_rar_50
summary_delta_no_rar_50$method<-"At 50% enrollment"
summary_delta_no_rar_50$total_participants<-summary_delta_no_rar_50$total_participants+
  summary_delta_no_interim$total_participants/2

summary_delta_no_rar50_plot<-summary_delta_no_rar_50


ggplot(summary_delta_no_rar50_plot,aes(x=delta,y=total_participants))+geom_point(size=1.5)+geom_line(size=1.5) 


master_dataset<-rbind(summary_delta_no_interim,summary_delta_no_rar_50)
ggplot(master_dataset,aes(x=delta,y=total_participants,color=method))+geom_point(size=1.5)+geom_line(size=1.5)+theme_classic()+
  ylab("Total Number of Participants")+
  theme(legend.position = "bottom")+xlab("δ*")+theme(text = element_text(size=36))+
  labs(color="Interim Analysis")

df3<-read_csv("changing_delta_no_rar_interim_analysis25.csv")
df3$Total_participants<-df3$treatment1+df3$treatment2+
  df3$treatment3+df3$treatment4
df3%>% 
  group_by(delta) %>% 
  summarize(total_participants = mean(Total_participants))->summary_delta_no_rar_25

summary_delta_no_rar_25$total_participants<-summary_delta_no_rar_25$total_participants+
  summary_delta_no_interim$total_participants/4
summary_delta_no_rar_25$method<-"At 25% enrollment"
ggplot(summary_delta_no_rar_25,aes(x=delta,y=total_participants))+geom_point(size=1.5)+geom_line(size=1.5)  
master_dataset<-rbind(master_dataset,summary_delta_no_rar_25)

ggplot(master_dataset,aes(x=delta,y=total_participants,color=method))+geom_point(size=1.5)+geom_line(size=1.5) 

df4<-read_csv("changing_delta_no_rar_interim_analysis125.csv")
df4$Total_participants<-df4$treatment1+df4$treatment2+
  df4$treatment3+df4$treatment4
df4%>% 
  group_by(delta) %>% 
  summarize(total_participants = mean(Total_participants))->summary_delta_no_rar_125

summary_delta_no_rar_125$total_participants<-summary_delta_no_rar_125$total_participants+
  summary_delta_no_interim$total_participants/8
summary_delta_no_rar_125$method<-"At 12.5% enrollment"

ggplot(summary_delta_no_rar_125,aes(x=delta,y=total_participants))+geom_point(size=1.5)+geom_line(size=1.5) 


master_dataset<-rbind(master_dataset,summary_delta_no_rar_125)

ggplot(master_dataset,aes(x=delta,y=total_participants,color=method))+geom_point(size=1.5)+geom_line(size=1.5) 


df5<-read_csv("changing_delta_no_rar_interim_analyses25_50.csv")
df5$Total_participants<-df5$treatment1+df5$treatment2+
  df5$treatment3+df5$treatment4
df5%>% 
  group_by(delta) %>% 
  summarize(total_participants = mean(Total_participants))->summary_delta_no_rar_25_50
summary_delta_no_rar_25_50$method<-"At 25% and 50% enrollment"
ggplot(summary_delta_no_rar_25_50,aes(x=delta,y=total_participants))+geom_point(size=1.5)+geom_line(size=1.5) 

master_dataset<-rbind(master_dataset,summary_delta_no_rar_25_50)

ggplot(master_dataset,aes(x=delta,y=total_participants,color=method))+geom_point(size=1.5)+geom_line(size=1.5) 


df6<-read_csv("changing_delta_no_rar_interim_analyses25_75.csv")
df6$Total_participants<-df6$treatment1+df6$treatment2+
  df6$treatment3+df6$treatment4
df6%>% 
  group_by(delta) %>% 
  summarize(total_participants = mean(Total_participants))->summary_delta_no_rar_25_75
summary_delta_no_rar_25_75$method<-"At 25% and 75% enrollment"
ggplot(summary_delta_no_rar_25_75,aes(x=delta,y=total_participants))+geom_point(size=1.5)+geom_line(size=1.5) 

master_dataset<-rbind(master_dataset,summary_delta_no_rar_25_75)
df_nhst<-read_csv("NHST.csv")
colnames(df_nhst)[2:3]<-c( "total_participants", "method" )
master_dataset<-rbind(master_dataset,df_nhst)


p1<-ggplot(master_dataset,aes(x=delta,y=total_participants,color=method))+
  geom_point(size=1.5)+geom_line(size=1.5)+theme_classic()+
  ylab("Total Number of Participants")+
  theme(legend.position = "bottom")+xlab("δ*")+theme(text = element_text(size=36))+
  labs(color="Interim Analysis") +scale_color_jama()



master_dataset_small_delta<-master_dataset%>%filter(delta<=0.25)

master_dataset_large_delta<-master_dataset%>%filter(delta>0.25)

p2<-ggplot(master_dataset_small_delta,aes(x=delta,y=total_participants,color=method))+geom_point(size=1.5)+geom_line(size=1.5)+theme_classic()+
  ylab("Total Number of Participants")+
  theme(legend.position = "bottom")+xlab("δ*")+theme(text = element_text(size=36))+
  labs(color="Interim Analysis") +scale_color_jama()

p3<-ggplot(master_dataset_large_delta,aes(x=delta,y=total_participants,color=method))+geom_point(size=1.5)+geom_line(size=1.5)+theme_classic()+
  ylab("Total Number of Participants")+
  theme(legend.position = "bottom")+xlab("δ*")+theme(text = element_text(size=36))+
  labs(color="Interim Analysis") +scale_color_jama()



###RAR


df7<-read_csv("changing_delta_rar_interim_analysis50.csv")
df7$Total_participants<-df7$treatment1+df7$treatment2+
  df7$treatment3+df7$treatment4
df7%>% 
  group_by(delta) %>% 
  summarize(total_participants = mean(Total_participants))->summary_delta_rar_50
summary_delta_rar_50$method<-"At 50% enrollment"
summary_delta_rar_50$total_participants<-summary_delta_rar_50$total_participants+
  summary_delta_no_interim$total_participants/2

summary_delta_rar50_plot<-summary_delta_rar_50


ggplot(summary_delta_rar50_plot,aes(x=delta,y=total_participants))+geom_point(size=1.5)+geom_line(size=1.5) 


master_dataset_rar<-rbind(summary_delta_no_interim,summary_delta_rar_50)
ggplot(master_dataset_rar,aes(x=delta,y=total_participants,color=method))+geom_point(size=1.5)+geom_line(size=1.5) 

df8<-read_csv("changing_delta_rar_interim_analysis25.csv")
df8$Total_participants<-df8$treatment1+df8$treatment2+
  df8$treatment3+df8$treatment4
df8%>% 
  group_by(delta) %>% 
  summarize(total_participants = mean(Total_participants,na.rm = T))->summary_delta_rar_25

summary_delta_rar_25$total_participants<-summary_delta_rar_25$total_participants+
  summary_delta_no_interim$total_participants/4
summary_delta_rar_25$method<-"At 25% enrollment"
ggplot(summary_delta_rar_25,aes(x=delta,y=total_participants))+geom_point(size=1.5)+geom_line(size=1.5)  
master_dataset_rar<-rbind(master_dataset_rar,summary_delta_rar_25)

ggplot(master_dataset_rar,aes(x=delta,y=total_participants,color=method))+geom_point(size=1.5)

df9<-read_csv("changing_delta_rar_interim_analysis125.csv")
df9$Total_participants<-df9$treatment1+df9$treatment2+
  df9$treatment3+df9$treatment4
df9%>% 
  group_by(delta) %>% 
  summarize(total_participants = mean(Total_participants,na.rm = T))->summary_delta_rar_125

summary_delta_rar_125$total_participants<-summary_delta_rar_125$total_participants+
  summary_delta_no_interim$total_participants/8
summary_delta_rar_125$method<-"At 12.5% enrollment"

ggplot(summary_delta_rar_125,aes(x=delta,y=total_participants))+geom_point(size=1.5)+geom_line(size=1.5) 


master_dataset_rar<-rbind(master_dataset_rar,summary_delta_rar_125)

ggplot(master_dataset_rar,aes(x=delta,y=total_participants,color=method))+geom_point(size=1.5)+geom_line(size=1.5) 


df10<-read_csv("changing_delta_rar_interim_analyses25_50.csv")
df10$Total_participants<-df10$treatment1+df10$treatment2+
  df10$treatment3+df10$treatment4
df10%>% 
  group_by(delta) %>% 
  summarize(total_participants = mean(Total_participants))->summary_delta_rar_25_50
summary_delta_rar_25_50$method<-"At 25% and 50% enrollment"
ggplot(summary_delta_rar_25_50,aes(x=delta,y=total_participants))+geom_point(size=1.5)+geom_line(size=1.5) 

master_dataset_rar<-rbind(master_dataset_rar,summary_delta_rar_25_50)

ggplot(master_dataset_rar,aes(x=delta,y=total_participants,color=method))+geom_point(size=1.5)+geom_line(size=1.5) 


df11<-read_csv("changing_delta_rar_interim_analyses25_75.csv")
df11$Total_participants<-df11$treatment1+df11$treatment2+
  df11$treatment3+df11$treatment4
df11%>% 
  group_by(delta) %>% 
  summarize(total_participants = mean(Total_participants))->summary_delta_rar_25_75
summary_delta_rar_25_75$method<-"At 25% and 75% enrollment"
ggplot(summary_delta_rar_25_75,aes(x=delta,y=total_participants))+geom_point(size=1.5)+geom_line(size=1.5) 

master_dataset_rar<-rbind(master_dataset_rar,summary_delta_rar_25_75)

ggplot(master_dataset_rar,aes(x=delta,y=total_participants,color=method))+geom_point(size=1.5)+geom_line(size=1.5) 
master_dataset_rar<-rbind(master_dataset_rar,df_nhst)

p4<-ggplot(master_dataset_rar,aes(x=delta,y=total_participants,color=method))+geom_point(size=1.5)+geom_line(size=1.5)+theme_classic()+
  ylab("Total Number of Participants")+
  theme(legend.position = "bottom")+xlab("δ*")+theme(text = element_text(size=36))+
  labs(color="Interim Analysis") +scale_color_jama()



master_dataset_small_delta_rar<-master_dataset_rar%>%filter(delta<=0.25)

master_dataset_large_delta_rar<-master_dataset_rar%>%filter(delta>0.25)

p5<-ggplot(master_dataset_small_delta_rar,aes(x=delta,y=total_participants,color=method))+geom_point(size=1.5)+geom_line(size=1.5)+theme_classic()+
  ylab("Total Number of Participants")+
  theme(legend.position = "bottom")+xlab("δ*")+theme(text = element_text(size=36))+
  labs(color="Interim Analysis") +scale_color_jama()

p6<-ggplot(master_dataset_large_delta_rar,aes(x=delta,y=total_participants,color=method))+geom_point(size=1.5)+geom_line(size=1.5)+theme_classic()+
  ylab("Total Number of Participants")+
  theme(legend.position = "bottom")+xlab("δ*")+theme(text = element_text(size=36))+
  labs(color="Interim Analysis") +scale_color_jama()


#########Change eta

df13<-read_csv("changing_eta_no_rar_interim_analysis50.csv")
df12<-read_csv("changing_eta_no_interim_analysis.csv")
df12$eta<-df13$eta
df12$Total_participants<-df12$treatment1+df12$treatment2+
  df12$treatment3+df12$treatment4
df12%>% 
  group_by(eta) %>% 
  summarize(total_participants = mean(Total_participants))->summary_eta_no_interim
summary_eta_no_interim_plot<-summary_eta_no_interim
summary_eta_no_interim$method<-"No"
ggplot(summary_eta_no_interim_plot,aes(x=eta,y=total_participants))+geom_point(size=1.5)+geom_line(size=1.5)+theme_bw()



df13<-read_csv("changing_eta_no_rar_interim_analysis50.csv")
df13$Total_participants<-df13$treatment1+df13$treatment2+
  df13$treatment3+df13$treatment4
df13%>% 
  group_by(eta) %>% 
  summarize(total_participants = mean(Total_participants))->summary_eta_no_rar_50
summary_eta_no_rar_50$method<-"At 50% enrollment"
summary_eta_no_rar_50$total_participants<-summary_eta_no_rar_50$total_participants+
  summary_eta_no_interim$total_participants/2

summary_eta_no_rar50_plot<-summary_eta_no_rar_50


ggplot(summary_eta_no_rar50_plot,aes(x=eta,y=total_participants))+geom_point(size=1.5)+geom_line(size=1.5) 


master_dataset_eta<-rbind(summary_eta_no_interim,summary_eta_no_rar_50)
ggplot(master_dataset_eta,aes(x=eta,y=total_participants,color=method))+geom_point(size=1.5)+geom_line(size=1.5)+theme_classic()+
  ylab("Total Number of Participants")+
  theme(legend.position = "bottom")+xlab("η")+theme(text = element_text(size=36))+
  labs(color="Interim Analysis")

df14<-read_csv("changing_eta_no_rar_interim_analysis25.csv")
df14$Total_participants<-df14$treatment1+df14$treatment2+
  df14$treatment3+df14$treatment4
df14%>% 
  group_by(eta) %>% 
  summarize(total_participants = mean(Total_participants))->summary_eta_no_rar_25

summary_eta_no_rar_25$total_participants<-summary_eta_no_rar_25$total_participants+
  summary_eta_no_interim$total_participants/4
summary_eta_no_rar_25$method<-"At 25% enrollment"
ggplot(summary_eta_no_rar_25,aes(x=eta,y=total_participants))+geom_point(size=1.5)+geom_line(size=1.5)  
master_dataset_eta<-rbind(master_dataset_eta,summary_eta_no_rar_25)

ggplot(master_dataset_eta,aes(x=eta,y=total_participants,color=method))+geom_point(size=1.5)+geom_line(size=1.5) 

df15<-read_csv("changing_eta_no_rar_interim_analysis125.csv")
df15$Total_participants<-df15$treatment1+df15$treatment2+
  df15$treatment3+df15$treatment4
df15%>% 
  group_by(eta) %>% 
  summarize(total_participants = mean(Total_participants))->summary_eta_no_rar_125

summary_eta_no_rar_125$total_participants<-summary_eta_no_rar_125$total_participants+
  summary_eta_no_interim$total_participants/8
summary_eta_no_rar_125$method<-"At 12.5% enrollment"

ggplot(summary_eta_no_rar_125,aes(x=eta,y=total_participants))+geom_point(size=1.5)+geom_line(size=1.5) 


master_dataset_eta<-rbind(master_dataset_eta,summary_eta_no_rar_125)

ggplot(master_dataset_eta,aes(x=eta,y=total_participants,color=method))+geom_point(size=1.5)+geom_line(size=1.5) 


df16<-read_csv("changing_eta_no_rar_interim_analyses25_50.csv")
df16$Total_participants<-df16$treatment1+df16$treatment2+
  df16$treatment3+df16$treatment4
df16%>% 
  group_by(eta) %>% 
  summarize(total_participants = mean(Total_participants))->summary_eta_no_rar_25_50
summary_eta_no_rar_25_50$method<-"At 25% and 50% enrollment"
ggplot(summary_eta_no_rar_25_50,aes(x=eta,y=total_participants))+geom_point(size=1.5)+geom_line(size=1.5) 

master_dataset_eta<-rbind(master_dataset_eta,summary_eta_no_rar_25_50)

ggplot(master_dataset_eta,aes(x=eta,y=total_participants,color=method))+geom_point(size=1.5)+geom_line(size=1.5) 


df17<-read_csv("changing_eta_no_rar_interim_analyses25_75.csv")
df17$Total_participants<-df17$treatment1+df17$treatment2+
  df17$treatment3+df17$treatment4
df17%>% 
  group_by(eta) %>% 
  summarize(total_participants = mean(Total_participants))->summary_eta_no_rar_25_75
summary_eta_no_rar_25_75$method<-"At 25% and 75% enrollment"
ggplot(summary_eta_no_rar_25_75,aes(x=eta,y=total_participants))+geom_point(size=1.5)+geom_line(size=1.5) 

master_dataset_eta<-rbind(master_dataset_eta,summary_eta_no_rar_25_75)

p7<-ggplot(master_dataset_eta,aes(x=eta,y=total_participants,color=method))+geom_point(size=1.5)+geom_line(size=1.5)+theme_classic()+
  ylab("Total Number of Participants")+
  theme(legend.position = "bottom")+xlab("η")+theme(text = element_text(size=36))+geom_vline(xintercept =0.95,linetype="dashed",color="gray")+
  labs(color="Interim Analysis") +scale_color_jama()

###RAR


df18<-read_csv("changing_eta_rar_interim_analysis50.csv")
df18$Total_participants<-df18$treatment1+df18$treatment2+
  df18$treatment3+df18$treatment4
df18%>% 
  group_by(eta) %>% 
  summarize(total_participants = mean(Total_participants))->summary_eta_rar_50
summary_eta_rar_50$method<-"At 50% enrollment"
summary_eta_rar_50$total_participants<-summary_eta_rar_50$total_participants+
  summary_eta_no_interim$total_participants/2

summary_eta_rar50_plot<-summary_eta_rar_50


ggplot(summary_eta_rar50_plot,aes(x=eta,y=total_participants))+geom_point(size=1.5)+geom_line(size=1.5) 


master_dataset_eta_rar<-rbind(summary_eta_no_interim,summary_eta_rar_50)
ggplot(master_dataset_eta_rar,aes(x=eta,y=total_participants,color=method))+geom_point(size=1.5)+geom_line(size=1.5) 

df19<-read_csv("changing_eta_rar_interim_analysis25.csv")
df19$Total_participants<-df19$treatment1+df19$treatment2+
  df19$treatment3+df19$treatment4
df19%>% 
  group_by(eta) %>% 
  summarize(total_participants = mean(Total_participants,na.rm = T))->summary_eta_rar_25

summary_eta_rar_25$total_participants<-summary_eta_rar_25$total_participants+
  summary_eta_no_interim$total_participants/4
summary_eta_rar_25$method<-"At 25% enrollment"
ggplot(summary_eta_rar_25,aes(x=eta,y=total_participants))+geom_point(size=1.5)+geom_line(size=1.5)  
master_dataset_eta_rar<-rbind(master_dataset_eta_rar,summary_eta_rar_25)

ggplot(master_dataset_eta_rar,aes(x=eta,y=total_participants,color=method))+geom_point(size=1.5)

df20<-read_csv("changing_eta_rar_interim_analysis125.csv")
df20$Total_participants<-df20$treatment1+df20$treatment2+
  df20$treatment3+df20$treatment4
df20%>% 
  group_by(eta) %>% 
  summarize(total_participants = mean(Total_participants,na.rm = T))->summary_eta_rar_125

summary_eta_rar_125$total_participants<-summary_eta_rar_125$total_participants+
  summary_eta_no_interim$total_participants/8
summary_eta_rar_125$method<-"At 12.5% enrollment"

ggplot(summary_eta_rar_125,aes(x=eta,y=total_participants))+geom_point(size=1.5)+geom_line(size=1.5) 


master_dataset_eta_rar<-rbind(master_dataset_eta_rar,summary_eta_rar_125)

ggplot(master_dataset_eta_rar,aes(x=eta,y=total_participants,color=method))+geom_point(size=1.5)+geom_line(size=1.5) 


df21<-read_csv("changing_eta_rar_interim_analyses25_50.csv")
df21$Total_participants<-df21$treatment1+df21$treatment2+
  df21$treatment3+df21$treatment4
df21%>% 
  group_by(eta) %>% 
  summarize(total_participants = mean(Total_participants))->summary_eta_rar_25_50
summary_eta_rar_25_50$method<-"At 25% and 50% enrollment"
ggplot(summary_eta_rar_25_50,aes(x=eta,y=total_participants))+geom_point(size=1.5)+geom_line(size=1.5) 

master_dataset_eta_rar<-rbind(master_dataset_eta_rar,summary_eta_rar_25_50)

ggplot(master_dataset_eta_rar,aes(x=eta,y=total_participants,color=method))+geom_point(size=1.5)+geom_line(size=1.5) 


df22<-read_csv("changing_eta_rar_interim_analyses25_75.csv")
df22$Total_participants<-df22$treatment1+df22$treatment2+
  df22$treatment3+df22$treatment4
df22%>% 
  group_by(eta) %>% 
  summarize(total_participants = mean(Total_participants))->summary_eta_rar_25_75
summary_eta_rar_25_75$method<-"At 25% and 75% enrollment"
ggplot(summary_eta_rar_25_75,aes(x=eta,y=total_participants))+geom_point(size=1.5)+geom_line(size=1.5) 

master_dataset_eta_rar<-rbind(master_dataset_eta_rar,summary_eta_rar_25_75)

ggplot(master_dataset_eta_rar,aes(x=eta,y=total_participants,color=method))+geom_point(size=1.5)+geom_line(size=1.5) 
library(ggsci)

p8<-ggplot(master_dataset_eta_rar,aes(x=eta,y=total_participants,color=method))+geom_point(size=1.5)+geom_line(size=1.5)+theme_classic()+
  ylab("Total Number of Participants")+
  theme(legend.position = "bottom")+xlab("η")+theme(text = element_text(size=36))+geom_vline(xintercept =0.95,linetype="dashed",color="gray")+
  labs(color="Interim Analysis")+scale_color_jama()


#########Change zeta

df24<-read_csv("changing_zeta_no_rar_interim_analysis50.csv")
df23<-read_csv("changing_zeta_no_interim_analysis.csv")
df23$zeta<-df24$zeta
df23$Total_participants<-df23$treatment1+df23$treatment2+
  df23$treatment3+df23$treatment4
df23%>% 
  group_by(zeta) %>% 
  summarize(total_participants = mean(Total_participants))->summary_zeta_no_interim
summary_zeta_no_interim_plot<-summary_zeta_no_interim
summary_zeta_no_interim$method<-"No"
ggplot(summary_zeta_no_interim_plot,aes(x=zeta,y=total_participants))+geom_point(size=1.5)+geom_line(size=1.5)+theme_bw()



df24<-read_csv("changing_zeta_no_rar_interim_analysis50.csv")
df24$Total_participants<-df24$treatment1+df24$treatment2+
  df24$treatment3+df24$treatment4
df24%>% 
  group_by(zeta) %>% 
  summarize(total_participants = mean(Total_participants))->summary_zeta_no_rar_50
summary_zeta_no_rar_50$method<-"At 50% enrollment"
summary_zeta_no_rar_50$total_participants<-summary_zeta_no_rar_50$total_participants+
  summary_zeta_no_interim$total_participants/2

summary_zeta_no_rar50_plot<-summary_zeta_no_rar_50


ggplot(summary_zeta_no_rar50_plot,aes(x=zeta,y=total_participants))+geom_point(size=1.5)+geom_line(size=1.5) 


master_dataset_zeta<-rbind(summary_zeta_no_interim,summary_zeta_no_rar_50)
ggplot(master_dataset_zeta,aes(x=zeta,y=total_participants,color=method))+geom_point(size=1.5)+geom_line(size=1.5)+theme_classic()+
  ylab("Total Number of Participants")+
  theme(legend.position = "bottom")+xlab("ζ")+theme(text = element_text(size=36))+
  labs(color="Interim Analysis")

df25<-read_csv("changing_zeta_no_rar_interim_analysis25.csv")
df25$Total_participants<-df25$treatment1+df25$treatment2+
  df25$treatment3+df25$treatment4
df25%>% 
  group_by(zeta) %>% 
  summarize(total_participants = mean(Total_participants))->summary_zeta_no_rar_25

summary_zeta_no_rar_25$total_participants<-summary_zeta_no_rar_25$total_participants+
  summary_zeta_no_interim$total_participants/4
summary_zeta_no_rar_25$method<-"At 25% enrollment"
ggplot(summary_zeta_no_rar_25,aes(x=zeta,y=total_participants))+geom_point(size=1.5)+geom_line(size=1.5)  
master_dataset_zeta<-rbind(master_dataset_zeta,summary_zeta_no_rar_25)

ggplot(master_dataset_zeta,aes(x=zeta,y=total_participants,color=method))+geom_point(size=1.5)+geom_line(size=1.5) 

df26<-read_csv("changing_zeta_no_rar_interim_analysis125.csv")
df26$Total_participants<-df26$treatment1+df26$treatment2+
  df26$treatment3+df26$treatment4
df26%>% 
  group_by(zeta) %>% 
  summarize(total_participants = mean(Total_participants))->summary_zeta_no_rar_125

summary_zeta_no_rar_125$total_participants<-summary_zeta_no_rar_125$total_participants+
  summary_zeta_no_interim$total_participants/8
summary_zeta_no_rar_125$method<-"At 12.5% enrollment"

ggplot(summary_zeta_no_rar_125,aes(x=zeta,y=total_participants))+geom_point(size=1.5)+geom_line(size=1.5) 


master_dataset_zeta<-rbind(master_dataset_zeta,summary_zeta_no_rar_125)

ggplot(master_dataset_zeta,aes(x=zeta,y=total_participants,color=method))+geom_point(size=1.5)+geom_line(size=1.5) 


df27<-read_csv("changing_zeta_no_rar_interim_analyses25_50.csv")
df27$Total_participants<-df27$treatment1+df27$treatment2+
  df27$treatment3+df27$treatment4
df27%>% 
  group_by(zeta) %>% 
  summarize(total_participants = mean(Total_participants))->summary_zeta_no_rar_25_50
summary_zeta_no_rar_25_50$method<-"Interims at 25% and 50% enrollment"
ggplot(summary_zeta_no_rar_25_50,aes(x=zeta,y=total_participants))+geom_point(size=1.5)+geom_line(size=1.5) 

master_dataset_zeta<-rbind(master_dataset_zeta,summary_zeta_no_rar_25_50)

ggplot(master_dataset_zeta,aes(x=zeta,y=total_participants,color=method))+geom_point(size=1.5)+geom_line(size=1.5) 


df28<-read_csv("changing_zeta_no_rar_interim_analyses25_75.csv")
df28$Total_participants<-df28$treatment1+df28$treatment2+
  df28$treatment3+df28$treatment4
df28%>% 
  group_by(zeta) %>% 
  summarize(total_participants = mean(Total_participants))->summary_zeta_no_rar_25_75
summary_zeta_no_rar_25_75$method<-"At 25% and 75% enrollment"
ggplot(summary_zeta_no_rar_25_75,aes(x=zeta,y=total_participants))+geom_point(size=1.5)+geom_line(size=1.5) 

master_dataset_zeta<-rbind(master_dataset_zeta,summary_zeta_no_rar_25_75)

p9<-ggplot(master_dataset_zeta,aes(x=zeta,y=total_participants,color=method))+geom_point(size=1.5)+geom_line(size=1.5)+theme_classic()+
  ylab("Total Number of Participants")+
  theme(legend.position = "bottom")+xlab("ζ")+theme(text = element_text(size=36))+geom_vline(xintercept =0.90,linetype="dashed",color="gray")+
  labs(color="Interim Analysis")+scale_color_jama()

###RAR


df29<-read_csv("changing_zeta_rar_interim_analysis50.csv")
df29$Total_participants<-df29$treatment1+df29$treatment2+
  df29$treatment3+df29$treatment4
df29%>% 
  group_by(zeta) %>% 
  summarize(total_participants = mean(Total_participants))->summary_zeta_rar_50
summary_zeta_rar_50$method<-"At 50% enrollment"
summary_zeta_rar_50$total_participants<-summary_zeta_rar_50$total_participants+
  summary_zeta_no_interim$total_participants/2

summary_zeta_rar50_plot<-summary_zeta_rar_50


ggplot(summary_zeta_rar50_plot,aes(x=zeta,y=total_participants))+geom_point(size=1.5)+geom_line(size=1.5) 


master_dataset_zeta_rar<-rbind(summary_zeta_no_interim,summary_zeta_rar_50)
ggplot(master_dataset_zeta_rar,aes(x=zeta,y=total_participants,color=method))+geom_point(size=1.5)+geom_line(size=1.5) 

df30<-read_csv("changing_zeta_rar_interim_analysis25.csv")
df30$Total_participants<-df30$treatment1+df30$treatment2+
  df30$treatment3+df30$treatment4
df30%>% 
  group_by(zeta) %>% 
  summarize(total_participants = mean(Total_participants,na.rm = T))->summary_zeta_rar_25

summary_zeta_rar_25$total_participants<-summary_zeta_rar_25$total_participants+
  summary_zeta_no_interim$total_participants/4
summary_zeta_rar_25$method<-"At 25% enrollment"
ggplot(summary_zeta_rar_25,aes(x=zeta,y=total_participants))+geom_point(size=1.5)+geom_line(size=1.5)  
master_dataset_zeta_rar<-rbind(master_dataset_zeta_rar,summary_zeta_rar_25)

ggplot(master_dataset_zeta_rar,aes(x=zeta,y=total_participants,color=method))+geom_point(size=1.5)

df31<-read_csv("changing_zeta_rar_interim_analysis125.csv")
df31$Total_participants<-df31$treatment1+df31$treatment2+
  df31$treatment3+df31$treatment4
df31%>% 
  group_by(zeta) %>% 
  summarize(total_participants = mean(Total_participants,na.rm = T))->summary_zeta_rar_125

summary_zeta_rar_125$total_participants<-summary_zeta_rar_125$total_participants+
  summary_zeta_no_interim$total_participants/8
summary_zeta_rar_125$method<-"At 12.5% enrollment"

ggplot(summary_zeta_rar_125,aes(x=zeta,y=total_participants))+geom_point(size=1.5)+geom_line(size=1.5) 


master_dataset_zeta_rar<-rbind(master_dataset_zeta_rar,summary_zeta_rar_125)

ggplot(master_dataset_zeta_rar,aes(x=zeta,y=total_participants,color=method))+geom_point(size=1.5)+geom_line(size=1.5) 


df32<-read_csv("changing_zeta_rar_interim_analyses25_50.csv")
df32$Total_participants<-df32$treatment1+df32$treatment2+
  df32$treatment3+df32$treatment4
df32%>% 
  group_by(zeta) %>% 
  summarize(total_participants = mean(Total_participants))->summary_zeta_rar_25_50
summary_zeta_rar_25_50$method<-"At 25% and 50% enrollment"
ggplot(summary_zeta_rar_25_50,aes(x=zeta,y=total_participants))+geom_point(size=1.5)+geom_line(size=1.5) 

master_dataset_zeta_rar<-rbind(master_dataset_zeta_rar,summary_zeta_rar_25_50)

ggplot(master_dataset_zeta_rar,aes(x=zeta,y=total_participants,color=method))+geom_point(size=1.5)+geom_line(size=1.5) 


df33<-read_csv("changing_zeta_rar_interim_analyses25_75.csv")
df33$Total_participants<-df33$treatment1+df33$treatment2+
  df33$treatment3+df33$treatment4
df33%>% 
  group_by(zeta) %>% 
  summarize(total_participants = mean(Total_participants))->summary_zeta_rar_25_75
summary_zeta_rar_25_75$method<-"At 25% and 75% enrollment"
ggplot(summary_zeta_rar_25_75,aes(x=zeta,y=total_participants))+geom_point(size=1.5)+geom_line(size=1.5) 

master_dataset_zeta_rar<-rbind(master_dataset_zeta_rar,summary_zeta_rar_25_75)

ggplot(master_dataset_zeta_rar,aes(x=zeta,y=total_participants,color=method))+geom_point(size=1.5)+geom_line(size=1.5) 


p10<-ggplot(master_dataset_zeta_rar,aes(x=zeta,y=total_participants,color=method))+geom_point(size=1.5)+geom_line(size=1.5)+theme_classic()+
  ylab("Total Number of Participants")+
  theme(legend.position = "bottom")+xlab("ζ")+theme(text = element_text(size=36))+geom_vline(xintercept =0.90,linetype="dashed",color="gray")+
  labs(color="Interim Analysis")+scale_color_jama() 



###Changing xi

df34<-read_csv("changing_xi_no_interim_analysis.csv")
df34$Total_participants<-df34$treatment1+df34$treatment2+
  df34$treatment3+df34$treatment4
df34%>% 
  group_by(xi) %>% 
  summarize(total_participants = mean(Total_participants))->summary_xi_no_interim
summary_xi_no_interim_plot<-summary_xi_no_interim
summary_xi_no_interim$method<-"No"
ggplot(summary_xi_no_interim_plot,aes(x=xi,y=total_participants))+geom_point(size=1.5)+geom_line(size=1.5)+theme_bw()



df35<-read_csv("changing_xi_no_rar_interim_analysis50.csv")
df35$Total_participants<-df35$treatment1+df35$treatment2+
  df35$treatment3+df35$treatment4
df35%>% 
  group_by(xi) %>% 
  summarize(total_participants = mean(Total_participants))->summary_xi_no_rar_50
summary_xi_no_rar_50$method<-"At 50% enrollment"
summary_xi_no_rar_50$total_participants<-summary_xi_no_rar_50$total_participants+
  summary_xi_no_interim$total_participants/2

summary_xi_no_rar50_plot<-summary_xi_no_rar_50


ggplot(summary_xi_no_rar50_plot,aes(x=xi,y=total_participants))+geom_point(size=1.5)+geom_line(size=1.5) 


master_dataset_xi<-rbind(summary_xi_no_interim,summary_xi_no_rar_50)
ggplot(master_dataset_xi,aes(x=xi,y=total_participants,color=method))+geom_point(size=1.5)+geom_line(size=1.5)+theme_classic()+
  ylab("Total Number of Participants")+
  theme(legend.position = "bottom")+xlab("ζ")+theme(text = element_text(size=36))+
  labs(color="Interim Analysis")

df36<-read_csv("changing_xi_no_rar_interim_analysis25.csv")
df36$Total_participants<-df36$treatment1+df36$treatment2+
  df36$treatment3+df36$treatment4
df36%>% 
  group_by(xi) %>% 
  summarize(total_participants = mean(Total_participants))->summary_xi_no_rar_25

summary_xi_no_rar_25$total_participants<-summary_xi_no_rar_25$total_participants+
  summary_xi_no_interim$total_participants/4
summary_xi_no_rar_25$method<-"At 25% enrollment"
ggplot(summary_xi_no_rar_25,aes(x=xi,y=total_participants))+geom_point(size=1.5)+geom_line(size=1.5)  
master_dataset_xi<-rbind(master_dataset_xi,summary_xi_no_rar_25)

ggplot(master_dataset_xi,aes(x=xi,y=total_participants,color=method))+geom_point(size=1.5)+geom_line(size=1.5) 

df37<-read_csv("changing_xi_no_rar_interim_analysis125.csv")
df37$Total_participants<-df37$treatment1+df37$treatment2+
  df37$treatment3+df37$treatment4
df37%>% 
  group_by(xi) %>% 
  summarize(total_participants = mean(Total_participants))->summary_xi_no_rar_125

summary_xi_no_rar_125$total_participants<-summary_xi_no_rar_125$total_participants+
  summary_xi_no_interim$total_participants/8
summary_xi_no_rar_125$method<-"At 12.5% enrollment"

ggplot(summary_xi_no_rar_125,aes(x=xi,y=total_participants))+geom_point(size=1.5)+geom_line(size=1.5) 


master_dataset_xi<-rbind(master_dataset_xi,summary_xi_no_rar_125)

ggplot(master_dataset_xi,aes(x=xi,y=total_participants,color=method))+geom_point(size=1.5)+geom_line(size=1.5) 


df38<-read_csv("changing_xi_no_rar_interim_analyses25_50.csv")
df38$Total_participants<-df38$treatment1+df38$treatment2+
  df38$treatment3+df38$treatment4
df38%>% 
  group_by(xi) %>% 
  summarize(total_participants = mean(Total_participants))->summary_xi_no_rar_25_50
summary_xi_no_rar_25_50$method<-"At 25% and 50% enrollment"
ggplot(summary_xi_no_rar_25_50,aes(x=xi,y=total_participants))+geom_point(size=1.5)+geom_line(size=1.5) 


master_dataset_xi<-rbind(master_dataset_xi,summary_xi_no_rar_25_50)

ggplot(master_dataset_xi,aes(x=xi,y=total_participants,color=method))+geom_point(size=1.5)+geom_line(size=1.5) 


df39<-read_csv("changing_xi_no_rar_interim_analyses25_75.csv")
df39$Total_participants<-df39$treatment1+df39$treatment2+
  df39$treatment3+df39$treatment4
df39%>% 
  group_by(xi) %>% 
  summarize(total_participants = mean(Total_participants))->summary_xi_no_rar_25_75
summary_xi_no_rar_25_75$method<-"At 25% and 75% enrollment"
ggplot(summary_xi_no_rar_25_75,aes(x=xi,y=total_participants))+geom_point(size=1.5)+geom_line(size=1.5) 

master_dataset_xi<-rbind(master_dataset_xi,summary_xi_no_rar_25_75)

p11<-ggplot(master_dataset_xi,aes(x=xi,y=total_participants,color=method))+geom_point(size=1.5)+geom_line(size=1.5)+theme_classic()+
  ylab("Total Number of Participants")+
  theme(legend.position = "bottom")+xlab("ξ")+theme(text = element_text(size=36))+geom_vline(xintercept =0.95,linetype="dashed",color="gray")+
  labs(color="Interim Analysis")+scale_color_jama()

###RAR


df40<-read_csv("changing_xi_rar_interim_analysis50.csv")
df40$Total_participants<-df40$treatment1+df40$treatment2+
  df40$treatment3+df40$treatment4
df40%>% 
  group_by(xi) %>% 
  summarize(total_participants = mean(Total_participants))->summary_xi_rar_50
summary_xi_rar_50$method<-"At 50% enrollment"
summary_xi_rar_50$total_participants<-summary_xi_rar_50$total_participants+
  summary_xi_no_interim$total_participants/2

summary_xi_rar50_plot<-summary_xi_rar_50


ggplot(summary_xi_rar50_plot,aes(x=xi,y=total_participants))+geom_point(size=1.5)+geom_line(size=1.5) 


master_dataset_xi_rar<-rbind(summary_xi_no_interim,summary_xi_rar_50)
ggplot(master_dataset_xi_rar,aes(x=xi,y=total_participants,color=method))+geom_point(size=1.5)+geom_line(size=1.5) 

df41<-read_csv("changing_xi_rar_interim_analysis25.csv")
df41$Total_participants<-df41$treatment1+df41$treatment2+
  df41$treatment3+df41$treatment4
df41%>% 
  group_by(xi) %>% 
  summarize(total_participants = mean(Total_participants,na.rm = T))->summary_xi_rar_25

summary_xi_rar_25$total_participants<-summary_xi_rar_25$total_participants+
  summary_xi_no_interim$total_participants/4
summary_xi_rar_25$method<-"At 25% enrollment"
ggplot(summary_xi_rar_25,aes(x=xi,y=total_participants))+geom_point(size=1.5)+geom_line(size=1.5)  
master_dataset_xi_rar<-rbind(master_dataset_xi_rar,summary_xi_rar_25)

ggplot(master_dataset_xi_rar,aes(x=xi,y=total_participants,color=method))+geom_point(size=1.5)

df42<-read_csv("changing_xi_rar_interim_analysis125.csv")
df42$Total_participants<-df42$treatment1+df42$treatment2+
  df42$treatment3+df42$treatment4
df42%>% 
  group_by(xi) %>% 
  summarize(total_participants = mean(Total_participants,na.rm = T))->summary_xi_rar_125

summary_xi_rar_125$total_participants<-summary_xi_rar_125$total_participants+
  summary_xi_no_interim$total_participants/8
summary_xi_rar_125$method<-"At 12.5% enrollment"

ggplot(summary_xi_rar_125,aes(x=xi,y=total_participants))+geom_point(size=1.5)+geom_line(size=1.5) 


master_dataset_xi_rar<-rbind(master_dataset_xi_rar,summary_xi_rar_125)

ggplot(master_dataset_xi_rar,aes(x=xi,y=total_participants,color=method))+geom_point(size=1.5)+geom_line(size=1.5) 


df43<-read_csv("changing_xi_rar_interim_analyses25_50.csv")
df43$Total_participants<-df43$treatment1+df43$treatment2+
  df43$treatment3+df43$treatment4
df43%>% 
  group_by(xi) %>% 
  summarize(total_participants = mean(Total_participants))->summary_xi_rar_25_50
summary_xi_rar_25_50$method<-"At 25% and 50% enrollment"
ggplot(summary_xi_rar_25_50,aes(x=xi,y=total_participants))+geom_point(size=1.5)+geom_line(size=1.5) 

master_dataset_xi_rar<-rbind(master_dataset_xi_rar,summary_xi_rar_25_50)

ggplot(master_dataset_xi_rar,aes(x=xi,y=total_participants,color=method))+geom_point(size=1.5)+geom_line(size=1.5) 


df44<-read_csv("changing_xi_rar_interim_analyses25_75.csv")
df44$Total_participants<-df44$treatment1+df44$treatment2+
  df44$treatment3+df44$treatment4
df44%>% 
  group_by(xi) %>% 
  summarize(total_participants = mean(Total_participants))->summary_xi_rar_25_75
summary_xi_rar_25_75$method<-"At 25% and 75% enrollment"
ggplot(summary_xi_rar_25_75,aes(x=xi,y=total_participants))+geom_point(size=1.5)+geom_line(size=1.5) 

master_dataset_xi_rar<-rbind(master_dataset_xi_rar,summary_xi_rar_25_75)

ggplot(master_dataset_xi_rar,aes(x=xi,y=total_participants,color=method))+geom_point(size=1.5)+geom_line(size=1.5) 


p12<-ggplot(master_dataset_xi_rar,aes(x=xi,y=total_participants,color=method))+geom_point(size=1.5)+geom_line(size=1.5)+theme_classic()+
  ylab("Total Number of Participants")+
  theme(legend.position = "bottom")+xlab("ξ")+theme(text = element_text(size=36))+geom_vline(xintercept =0.95,linetype="dashed",color="gray")+
  labs(color="Interim Analysis")+scale_color_jama() 





p13<-ggarrange(p1,p2,p3,p4,p5,p6, nrow=2,ncol=3, common.legend = T, legend="bottom",labels = "AUTO",
          font.label=list(color="black",size=36))                 
ggsave("Delta_stars0902.svg",p13,height = 20,width = 36,dpi=600)
#ggsave("Delta_stars0814.png",p13,height = 20,width = 36,dpi=600)

p14<-ggarrange(  p7,p9,p11,p8,p10,p12,nrow=2,ncol=3, common.legend = T, legend="bottom",labels = "AUTO",
               font.label=list(color="black",size=36))             

ggsave("eta_zeta_xi08902.svg",p14,height = 20,width = 36,dpi=600)
#ggsave("eta_zeta_xi0814.png",p14,height = 20,width = 36,dpi=600)

