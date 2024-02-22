

library(survcomp)

rm(list=ls())
setwd("C:/Users/212768837/Downloads/survival/survival")

df_os<-read.csv("train_survival.csv",header = T)

dim(df_os) 

head(df_os)

fit <- coxph(Surv(time, status) ~ RS, data = df_os)
cindex <- concordance.index(predict(fit),surv.time = df_os$time, surv.event = df_os$status,method = "noether")
cindex$c.index; cindex$lower; cindex$upper

fit <- coxph(Surv(time, status) ~ DL, data = df_os)
cindex <- concordance.index(predict(fit),surv.time = df_os$time, surv.event = df_os$status,method = "noether")
cindex$c.index; cindex$lower; cindex$upper

fit <- coxph(Surv(time, status) ~  RS-DL, data = df_os)
cindex <- concordance.index(predict(fit),surv.time = df_os$time, surv.event = df_os$status,method = "noether")
cindex$c.index; cindex$lower; cindex$upper


###HR
###KM

library(survival)
library("survminer")

data.survdiff <- survdiff(Surv(time, status)~Risk1, data=df_os)
p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
HR = (data.survdiff$obs[2]/data.survdiff$exp[2])/(data.survdiff$obs[1]/data.survdiff$exp[1])
up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1])) 


HR;low95;up95

fit <- survfit(Surv(time, status) ~ Risk1, data = df_os)


ggsurvplot(fit, data = df_os,  conf.int = F,pval = T,legend.title="Risk",risk.table = T,
           legend.labs=c("Low","High"),
           #font.legend = 8,
           xlim = c(0,120),        
           break.time.by = 12,    
           ggtheme=theme_survminer(base_size = 20,font.x=18,font.y = 18,font.tickslab = 16),
           tables.theme =theme_survminer(base_size = 20,font.tickslab = 16,font.x=18,font.y = 18,),fontsize=6,pval.size=5)


data.survdiff <- survdiff(Surv(time, status)~Risk2, data=df_os)
p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
HR = (data.survdiff$obs[2]/data.survdiff$exp[2])/(data.survdiff$obs[1]/data.survdiff$exp[1])
up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1])) 


HR;low95;up95
fit <- survfit(Surv(time, status) ~ Risk2, data = df_os)


ggsurvplot(fit, data = df_os,  conf.int = F,pval = T,legend.title="Risk",risk.table = T,
           legend.labs=c("Low","High"),
           #font.legend = 8,
           xlim = c(0,120),        
           break.time.by = 12,    
           ggtheme=theme_survminer(base_size = 20,font.x=18,font.y = 18,font.tickslab = 16),
           tables.theme =theme_survminer(base_size = 20,font.tickslab = 16,font.x=18,font.y = 18,),fontsize=6,pval.size=5)


data.survdiff <- survdiff(Surv(time, status)~Risk3, data=df_os)
p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
HR = (data.survdiff$obs[2]/data.survdiff$exp[2])/(data.survdiff$obs[1]/data.survdiff$exp[1])
up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1])) 


HR;low95;up95

fit <- survfit(Surv(time, status) ~ Risk3, data = df_os)


ggsurvplot(fit, data = df_os,  conf.int = F,pval = T,legend.title="Risk",risk.table = T,
           legend.labs=c("Low","High"),
           #font.legend = 8,
           xlim = c(0,120),        
           break.time.by = 12,    
           ggtheme=theme_survminer(base_size = 20,font.x=18,font.y = 18,font.tickslab = 16),
           tables.theme =theme_survminer(base_size = 20,font.tickslab = 16,font.x=18,font.y = 18,),fontsize=6,pval.size=5)




