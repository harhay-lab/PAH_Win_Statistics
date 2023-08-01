library(WINS)
library(meta)
library(metafor)

######Set working directory
setwd("U:/...")

######read data
win<- read.csv("U:/....csv")

###### 1.0 Setup data
##In our case, we use the following hierarchy  for the events:
#deaths>lung transplant>hospitalization due to PAH>clinical improvement>clinical worsening
win1<-win
win1$arm<-win1$study_arm_assignment
win1$id<-win1$usubjid
win1$Delta_1<-win1$cedeath
win1$Delta_2<-win1$cetxp
win1$Delta_3<-win1$cehosp
win1$Delta_4<-win1$ceimprove
win1$Delta_5<-win1$cedecline
win1$Y_1<-win1$weeknumdeath
win1$Y_2<-win1$weeknumtxp
win1$Y_3<-win1$weeknumhosp
win1$Y_4<-win1$weeknumimprove
win1$Y_5<-win1$weeknumdecline

win2=win1



###### 2.0 Define win strategy to specify the improvement outcome(s)
#The function in the WINS package considers all outcome events as negative events. 
#Since one of our events is clinical improvement, we have to redefine the win strategy so that those with a shorter time to clinical improvement are considered as 'wins' 

######consider Tau as 0 since TTE
win.strategy.default1<-function(trt_con, priority){
  n_ep = length(priority)
  
  #### Obtain the indicator of the first endpoint for treatment and control
  colname.trt_con = colnames(trt_con)
  
  ind.delta1.trt = which(colname.trt_con=="Delta_1_trt")
  ind.delta1.con = which(colname.trt_con=="Delta_1_con")
  
  ind.time1.trt = which(colname.trt_con=="Y_1_trt")
  ind.time1.con = which(colname.trt_con=="Y_1_con")
  
  
  win.status0 = NULL
  
  #### For TTE outcome: Denote the observed survival time as Y_trt and Y_con, and event status
  #### as Delta_trt and Delta_con. There is a win for the treatment group if we have:
  ####                      Delta_con = 1 and Y_trt > Y_con + 0,
  #### where 0 denote the magnitude of difference to determine win/loss/tie.
  
  #### For continuous outcome: Denote the observed value as Y_trt and Y_con. The event status
  #### are Delta_trt = 1 and Delta_con = 1. There is a win for the treatment group if we have:
  ####                      Delta_con = 1 and Y_trt > Y_con + 0,
  #### where 0 denote the magnitude of difference to determine win/loss/tie.
  
  #### For binary outcome (0/1): Denote the observed value as Y_trt and Y_con. The event status
  #### are Delta_trt = 1 and Delta_con = 1. There is a win for the treatment group if we have:
  ####                      Delta_con = 1 and Y_trt > Y_con + 0,
  #### where tau is 0.
  for(l in priority){
    delta_l_trt = trt_con[,ind.delta1.trt+l-1]; delta_l_con = trt_con[,ind.delta1.con+l-1]
    Y_l_trt = trt_con[,ind.time1.trt+l-1]; Y_l_con = trt_con[,ind.time1.con+l-1]
    
    
    
    #user define direction of outcome (specify improvement outcome(s) in 'imp' vector)
    #For the improvement outcome(s) a subject wins if they have a lesser time to the improvement event
    #For the negative outcome(s) a subject wins if they have a longer time to the event
    if(l%in%imp){
      #### Treatment won, per Endpoint l
      win.temp1 = ( (delta_l_trt == 1 & delta_l_con == 1 & Y_l_trt < (Y_l_con + 0)) |
                      (delta_l_trt == 1 & delta_l_con < 1  & Y_l_trt < (Y_l_con + 0)) )
      
      #### Control won, per Endpoint l
      win.temp2 = ( (delta_l_trt == 1 & delta_l_con == 1 & Y_l_trt > (Y_l_con + 0)) |
                      (delta_l_trt < 1  & delta_l_con == 1 & (Y_l_con + 0) < Y_l_trt) )
      
    }else{
      #### Treatment won, per Endpoint l
      win.temp1 = ( (delta_l_trt == 1 & delta_l_con == 1 & Y_l_trt > (Y_l_con + 0)) |
                      (delta_l_trt < 1  & delta_l_con == 1 & (Y_l_con + 0) < Y_l_trt) )
      
      #### Control won, per Endpoint l
      win.temp2 = ( (delta_l_trt == 1 & delta_l_con == 1 & Y_l_trt < (Y_l_con + 0)) |
                      (delta_l_trt == 1 & delta_l_con < 1  & Y_l_trt < (Y_l_con + 0)) )
      
    }
    win.status0 = cbind(win.status0,win.temp1,win.temp2)
  }
  
  #### prioritize: once a winner is determined, then all the subsequent is set to zero
  win_status = t(apply(win.status0, 1, func<-function(x){
    if(sum(x)>1){
      temp = x; temp[min((2*min(ceiling(which(x==1)/2))+1),(2*n_ep-1)):(2*n_ep)] = 0
      return(temp)
    }else{
      return(as.numeric(x))
    }
  }))
  
  colnames(win_status) = paste0(rep(c("Trt","Con"),n_ep),"_Endpoint",rep(priority,each=2))
  win_status = as.data.frame(win_status)
  
  return(win_status)
}


######set vector for the improvement outcome(s)
#In our case, the only improvement outcome is clinical improvement which is outcome 4 in our hierarchy
imp=4






###### 3.0 Compute unstratified win statistics without IPCW
res<-win.stat(data = win2, ep_type = "tte",arm.name = c("Active","Control"), tau = 0, priority = c(1:5), alpha = 0.05, digit = 5,
         weight = "unstratified", censoring_adjust ="No",pvalue = "two-sided",win.strategy = win.strategy.default1,summary.print=F)


###### 3.1 Overall number of wins, losses, and ties
res_all=res
res_all$Win_statistic
##wins
res_all$Win_prop[,2]*nrow(win2[win2$arm=="Active",])*nrow(win2[win2$arm=="Control",])
##losses
res_all$Win_prop[,3]*nrow(win2[win2$arm=="Active",])*nrow(win2[win2$arm=="Control",])
##ties
(1-(res_all$Win_prop[,2]+res_all$Win_prop[,3]))*nrow(win2[win2$arm=="Active",])*nrow(win2[win2$arm=="Control",])

###### 3.2 Number of wins, losses, and ties per outcome event (1 to 5 respectively)
##Wins
res$summary_ep$Trt_Endpoint1$Count
res$summary_ep$Trt_Endpoint2$Count
res$summary_ep$Trt_Endpoint3$Count
res$summary_ep$Trt_Endpoint4$Count
res$summary_ep$Trt_Endpoint5$Count
##losses
res$summary_ep$Con_Endpoint1$Count
res$summary_ep$Con_Endpoint2$Count
res$summary_ep$Con_Endpoint3$Count
res$summary_ep$Con_Endpoint4$Count
res$summary_ep$Con_Endpoint5$Count
##ties
T1=nrow(win2[win2$arm=="Active",])*nrow(win2[win2$arm=="Control",])-(res$summary_ep$Trt_Endpoint1$Count+res$summary_ep$Con_Endpoint1$Count);T1
T2=T1-(res$summary_ep$Trt_Endpoint2$Count+res$summary_ep$Con_Endpoint2$Count);T2
T3=T2-(res$summary_ep$Trt_Endpoint3$Count+res$summary_ep$Con_Endpoint3$Count);T3
T4=T3-(res$summary_ep$Trt_Endpoint4$Count+res$summary_ep$Con_Endpoint4$Count);T4
T5=T4-(res$summary_ep$Trt_Endpoint5$Count+res$summary_ep$Con_Endpoint5$Count);T5
##Win ratio, net benefit, and win odds per outcome event
WR1=res$summary_ep$Trt_Endpoint1$Count/res$summary_ep$Con_Endpoint1$Count;WR1
NB1=res$summary_ep$Trt_Endpoint1$Proportion-res$summary_ep$Con_Endpoint1$Proportion;NB1
WO1=(res$summary_ep$Trt_Endpoint1$Count+(0.5*T1))/(res$summary_ep$Con_Endpoint1$Count+(0.5*T1));WO1
WO1=(res$summary_ep$Trt_Endpoint1$Proportion+(0.5*(1-(res$summary_ep$Trt_Endpoint1$Proportion+res$summary_ep$Con_Endpoint1$Proportion))))/(res$summary_ep$Con_Endpoint1$Proportion+(0.5*(1-(res$summary_ep$Trt_Endpoint1$Proportion+res$summary_ep$Con_Endpoint1$Proportion))));WO1
WR2=res$summary_ep$Trt_Endpoint2$Count/res$summary_ep$Con_Endpoint2$Count;WR2
NB2=res$summary_ep$Trt_Endpoint2$Proportion-res$summary_ep$Con_Endpoint2$Proportion;NB2
WO2=(res$summary_ep$Trt_Endpoint2$Count+(0.5*T2))/(res$summary_ep$Con_Endpoint2$Count+(0.5*T2));WO2
WR3=res$summary_ep$Trt_Endpoint3$Count/res$summary_ep$Con_Endpoint3$Count;WR3
NB3=res$summary_ep$Trt_Endpoint3$Proportion-res$summary_ep$Con_Endpoint3$Proportion;NB3
WO3=(res$summary_ep$Trt_Endpoint3$Count+(0.5*T3))/(res$summary_ep$Con_Endpoint3$Count+(0.5*T3));WO3
WR4=res$summary_ep$Trt_Endpoint4$Count/res$summary_ep$Con_Endpoint4$Count;WR4
NB4=res$summary_ep$Trt_Endpoint4$Proportion-res$summary_ep$Con_Endpoint4$Proportion;NB4
WO4=(res$summary_ep$Trt_Endpoint4$Count+(0.5*T4))/(res$summary_ep$Con_Endpoint4$Count+(0.5*T4));WO4
WR5=res$summary_ep$Trt_Endpoint5$Count/res$summary_ep$Con_Endpoint5$Count;WR5
NB5=res$summary_ep$Trt_Endpoint5$Proportion-res$summary_ep$Con_Endpoint5$Proportion;NB5
WO5=(res$summary_ep$Trt_Endpoint5$Count+(0.5*T5))/(res$summary_ep$Con_Endpoint5$Count+(0.5*T5));WO5




###### 4.0 Compute unstratified win statistics with IPCW
res<-win.stat(data = win2, ep_type = "tte",arm.name = c("Active","Control"), tau = 0, priority = c(1:5), alpha = 0.05, digit = 5,
              weight = "unstratified", censoring_adjust ="IPCW",pvalue = "two-sided",win.strategy = win.strategy.default1,summary.print=F)
res$Win_statistic



###### 4.1 Plot of unstratified win statistics with IPCW over time
stat_t.plot(data=win2, arm.name = c("Active","Control"), priority = c(1:5),
            Ctime = seq(4,192,4),plotTimeUnit = "Weeks", statistic = c("WR"),
            weight = "unstratified", censoring_adjust ="IPCW",win.strategy = win.strategy.default1,
            tau = 0, plot_CI = TRUE)
stat_t.plot(data=win2, arm.name = c("Active","Control"), priority = c(1:5),
            Ctime = seq(4,192,4),plotTimeUnit = "Weeks", statistic = c("WO"),
            weight = "unstratified", censoring_adjust ="IPCW",win.strategy = win.strategy.default1,
            tau = 0, plot_CI = TRUE)
stat_t.plot(data=win2, arm.name = c("Active","Control"), priority = c(1:5),
            Ctime = seq(4,192,4),plotTimeUnit = "Weeks", statistic = c("NB"),
            weight = "unstratified", censoring_adjust ="IPCW",win.strategy = win.strategy.default1,
            tau = 0, plot_CI = TRUE)





###### 5.0 Compute win statistics with IPCW stratified by trial
#Two stage: compute win statistics within each trial, then meta-analyze across trials

####First stage
#Win statistics by trial
fit<-cbind.data.frame(study=NA, wr=NA, "wrul"=NA, "wrll"=NA, nb=NA, "nbul"=NA, "nbll"=NA, wo=NA, "woul"=NA, "woll"=NA)
for (study in unique(win2$studyid)) {
  res<-win.stat(data = win2[win2$studyid==study,], ep_type = "tte",arm.name = c("Active","Control"), tau = 0, priority = c(1:5), alpha = 0.05, digit = 5,
                weight = "unstratified", censoring_adjust ="IPCW",pvalue = "two-sided",win.strategy = win.strategy.default1,summary.print=F)
  
  fit1<-cbind.data.frame(study=study, wr=res$Win_statistic$Win_Ratio[1], "wrul"=res$Win_statistic$Win_Ratio[3], "wrll"=res$Win_statistic$Win_Ratio[2], 
                        nb=res$Win_statistic$Net_Benefit[1], "nbul"=res$Win_statistic$Net_Benefit[3], "nbll"=res$Win_statistic$Net_Benefit[2], 
                        wo=res$Win_statistic$Win_Odds[1], "woul"=res$Win_statistic$Win_Odds[3], "woll"=res$Win_statistic$Win_Odds[2])
  fit<-rbind.data.frame(fit,fit1)
}
fit<-fit[-1,]

fit_all<-fit


#### Second stage
## Inverse variance meta-analysis
fit <- fit[order(fit$study),]

#WR
met<-metagen(TE=wr, lower=`wrll`, upper = `wrul`,sm="OR",
             comb.fixed=F,hakn=T,transf=F, data=fit)
forest(met,
       leftcols = c("studlab","effect.ci"),rightcols = "w.random",
       leftlabs = c("Trial", NA),just = "center",title="yu",
       smlab="Win ratio",text.random="Total",digits=2,lwd=1.5,
       col.square="lightblue",col.diamond.fixed="darkblue",test.overall=T)

#WO
met<-metagen(TE=wo, lower=`woll`, upper = `woul`,sm="OR",
             comb.fixed=F,hakn=T,transf=F, data=fit)
forest(met,
       leftcols = c("studlab","effect.ci"),rightcols = "w.random",
       leftlabs = c("Trial", NA),just = "center",title="yu",
       smlab="Win odds",text.random="Total",digits=2,lwd=1.5,
       col.square="lightblue",col.diamond.fixed="darkblue",test.overall=T)

#NB
met<-metagen(TE=nb, lower=`nbll`, upper = `nbul`,sm="RD",
             comb.fixed=F,hakn=T,transf=F, data=fit)
forest(met,
       leftcols = c("studlab","effect.ci"),rightcols = "w.random",
       leftlabs = c("Trial", NA),just = "center",title="yu",
       smlab="Net benefit",text.random="Total",digits=2,lwd=1.5,
       col.square="lightblue",col.diamond.fixed="darkblue",test.overall=T)




###### 6.0 Compute win statistics with IPCW stratified by trial for each subgroup
##### 6.1 by trial follow-up
win2$trial<-ifelse(win2$studyid%in%c(1004,1010,1013),"event","fixed")
table(win2$trial)
table(win2$trial)/sum(table(win2$trial))

#### first stage
fit<-cbind.data.frame(study=NA, wr=NA, "wrul"=NA, "wrll"=NA, nb=NA, "nbul"=NA, "nbll"=NA, wo=NA, "woul"=NA, "woll"=NA)
for (study in unique(win2[win2$trial=="event",]$studyid)) {
  res<-win.stat(data = win2[win2$trial=="event"&win2$studyid==study,], ep_type = "tte",arm.name = c("Active","Control"), tau = 0, priority = c(1:5), alpha = 0.05, digit = 5,
                weight = "unstratified", censoring_adjust ="IPCW",pvalue = "two-sided",win.strategy = win.strategy.default1,summary.print=F)
  
  fit1<-cbind.data.frame(study=study, wr=res$Win_statistic$Win_Ratio[1], "wrul"=res$Win_statistic$Win_Ratio[3], "wrll"=res$Win_statistic$Win_Ratio[2], 
                         nb=res$Win_statistic$Net_Benefit[1], "nbul"=res$Win_statistic$Net_Benefit[3], "nbll"=res$Win_statistic$Net_Benefit[2], 
                         wo=res$Win_statistic$Win_Odds[1], "woul"=res$Win_statistic$Win_Odds[3], "woll"=res$Win_statistic$Win_Odds[2])
  fit<-rbind.data.frame(fit,fit1)
}
fit<-fit[-1,]
fit_e<-fit

fit<-cbind.data.frame(study=NA, wr=NA, "wrul"=NA, "wrll"=NA, nb=NA, "nbul"=NA, "nbll"=NA, wo=NA, "woul"=NA, "woll"=NA)
for (study in unique(win2[win2$trial=="fixed",]$studyid)) {
  res<-win.stat(data = win2[win2$trial=="fixed"&win2$studyid==study,], ep_type = "tte",arm.name = c("Active","Control"), tau = 0, priority = c(1:5), alpha = 0.05, digit = 5,
                weight = "unstratified", censoring_adjust ="IPCW",pvalue = "two-sided",win.strategy = win.strategy.default1,summary.print=F)
  
  fit1<-cbind.data.frame(study=study, wr=res$Win_statistic$Win_Ratio[1], "wrul"=res$Win_statistic$Win_Ratio[3], "wrll"=res$Win_statistic$Win_Ratio[2], 
                         nb=res$Win_statistic$Net_Benefit[1], "nbul"=res$Win_statistic$Net_Benefit[3], "nbll"=res$Win_statistic$Net_Benefit[2], 
                         wo=res$Win_statistic$Win_Odds[1], "woul"=res$Win_statistic$Win_Odds[3], "woll"=res$Win_statistic$Win_Odds[2])
  fit<-rbind.data.frame(fit,fit1)
}
fit<-fit[-1,]
fit_f<-fit

#### second stage
## Inverse variance meta-analysis
fit_e$trial<-"event"
fit_f$trial<-"fixed"
fit<-rbind.data.frame(fit_e, fit_f)

#WR
met<-metagen(TE=wr, lower=`wrll`, upper = `wrul`,sm="OR",
             comb.fixed=F,hakn=T,transf=F,subgroup=fit$trial,test.subgroup=T, data=fit)
forest(met,
       leftcols = c("studlab","effect.ci"),rightcols = "w.random",
       leftlabs = c("Trial", NA),just = "center",title="yu",
       smlab="Win ratio",text.random="Total",digits=2,lwd=1.5,
       col.square="lightblue",col.diamond.fixed="darkblue",test.overall=T)

#WO
met<-metagen(TE=wo, lower=`woll`, upper = `woul`,sm="OR",
             comb.fixed=F,hakn=T,transf=F,subgroup=fit$trial,test.subgroup=T, data=fit)
forest(met,
       leftcols = c("studlab","effect.ci"),rightcols = "w.random",
       leftlabs = c("Trial", NA),just = "center",title="yu",
       smlab="Win odds",text.random="Total",digits=2,lwd=1.5,
       col.square="lightblue",col.diamond.fixed="darkblue",test.overall=T)

#NB
met<-metagen(TE=nb, lower=`nbll`, upper = `nbul`,sm="RD",
             comb.fixed=F,hakn=T,transf=F,subgroup=fit$trial,test.subgroup=T, data=fit)
forest(met,
       leftcols = c("studlab","effect.ci"),rightcols = "w.random",
       leftlabs = c("Trial", NA),just = "center",title="yu",
       smlab="Win odds",text.random="Total",digits=2,lwd=1.5,
       col.square="lightblue",col.diamond.fixed="darkblue",test.overall=T)


##### 6.2 by age
win2$age1<-ifelse(win2$age<50,"<50",">=50")
table(win2$age1)
table(win2$age1)/sum(table(win2$age1))

#### first stage
fit<-cbind.data.frame(study=NA, wr=NA, "wrul"=NA, "wrll"=NA, nb=NA, "nbul"=NA, "nbll"=NA, wo=NA, "woul"=NA, "woll"=NA)
for (study in unique(win2[win2$age1=="<50"&!is.na(win2$age1),]$studyid)) {
  res<-win.stat(data = win2[win2$age1=="<50"&!is.na(win2$age1)&win2$studyid==study,], ep_type = "tte",arm.name = c("Active","Control"), tau = 0, priority = c(1:5), alpha = 0.05, digit = 5,
                weight = "unstratified", censoring_adjust ="IPCW",pvalue = "two-sided",win.strategy = win.strategy.default1,summary.print=F)
  
  fit1<-cbind.data.frame(study=study, wr=res$Win_statistic$Win_Ratio[1], "wrul"=res$Win_statistic$Win_Ratio[3], "wrll"=res$Win_statistic$Win_Ratio[2], 
                         nb=res$Win_statistic$Net_Benefit[1], "nbul"=res$Win_statistic$Net_Benefit[3], "nbll"=res$Win_statistic$Net_Benefit[2], 
                         wo=res$Win_statistic$Win_Odds[1], "woul"=res$Win_statistic$Win_Odds[3], "woll"=res$Win_statistic$Win_Odds[2])
  fit<-rbind.data.frame(fit,fit1)
}
fit<-fit[-1,]
fit_l<-fit

fit<-cbind.data.frame(study=NA, wr=NA, "wrul"=NA, "wrll"=NA, nb=NA, "nbul"=NA, "nbll"=NA, wo=NA, "woul"=NA, "woll"=NA)
for (study in unique(win2[win2$age1==">=50"&!is.na(win2$age1),]$studyid)) {
  res<-win.stat(data = win2[win2$age1==">=50"&!is.na(win2$age1)&win2$studyid==study,], ep_type = "tte",arm.name = c("Active","Control"), tau = 0, priority = c(1:5), alpha = 0.05, digit = 5,
                weight = "unstratified", censoring_adjust ="IPCW",pvalue = "two-sided",win.strategy = win.strategy.default1,summary.print=F)
  
  fit1<-cbind.data.frame(study=study, wr=res$Win_statistic$Win_Ratio[1], "wrul"=res$Win_statistic$Win_Ratio[3], "wrll"=res$Win_statistic$Win_Ratio[2], 
                         nb=res$Win_statistic$Net_Benefit[1], "nbul"=res$Win_statistic$Net_Benefit[3], "nbll"=res$Win_statistic$Net_Benefit[2], 
                         wo=res$Win_statistic$Win_Odds[1], "woul"=res$Win_statistic$Win_Odds[3], "woll"=res$Win_statistic$Win_Odds[2])
  fit<-rbind.data.frame(fit,fit1)
}
fit<-fit[-1,]
fit_u<-fit

#### second stage
## Inverse variance meta-analysis
fit_l$age<-"<50"
fit_u$age<-">=50"
fit<-rbind.data.frame(fit_l, fit_u)

#WR
met<-metagen(TE=wr, lower=`wrll`, upper = `wrul`,sm="OR",
             comb.fixed=F,hakn=T,transf=F,subgroup=fit$age,test.subgroup=T, data=fit)
forest(met,
       leftcols = c("studlab","effect.ci"),rightcols = "w.random",
       leftlabs = c("Trial", NA),just = "center",title="yu",
       smlab="Win ratio",text.random="Total",digits=2,lwd=1.5,
       col.square="lightblue",col.diamond.fixed="darkblue",test.overall=T)

#WO
met<-metagen(TE=wo, lower=`woll`, upper = `woul`,sm="OR",
             comb.fixed=F,hakn=T,transf=F,subgroup=fit$age,test.subgroup=T, data=fit)
forest(met,
       leftcols = c("studlab","effect.ci"),rightcols = "w.random",
       leftlabs = c("Trial", NA),just = "center",title="yu",
       smlab="Win odds",text.random="Total",digits=2,lwd=1.5,
       col.square="lightblue",col.diamond.fixed="darkblue",test.overall=T)

#NB
met<-metagen(TE=nb, lower=`nbll`, upper = `nbul`,sm="RD",
             comb.fixed=F,hakn=T,transf=F,subgroup=fit$age,test.subgroup=T, data=fit)
forest(met,
       leftcols = c("studlab","effect.ci"),rightcols = "w.random",
       leftlabs = c("Trial", NA),just = "center",title="yu",
       smlab="Win odds",text.random="Total",digits=2,lwd=1.5,
       col.square="lightblue",col.diamond.fixed="darkblue",test.overall=T)


##### 6.3 by etiology
win2$eti<-ifelse(win2$eti_code%in%c(1,2,3),"idio",
                  ifelse(win2$eti_code%in%c(4),"ctd","others"))
table(win2$eti)
table(win2$eti)/sum(table(win2$eti))

#### first stage
fit<-cbind.data.frame(study=NA, wr=NA, "wrul"=NA, "wrll"=NA, nb=NA, "nbul"=NA, "nbll"=NA, wo=NA, "woul"=NA, "woll"=NA)
for (study in unique(win2[win2$eti=="idio"&!is.na(win2$eti_code),]$studyid)) {
  res<-win.stat(data = win2[win2$eti=="idio"&!is.na(win2$eti_code)&win2$studyid==study,], ep_type = "tte",arm.name = c("Active","Control"), tau = 0, priority = c(1:5), alpha = 0.05, digit = 5,
                weight = "unstratified", censoring_adjust ="IPCW",pvalue = "two-sided",win.strategy = win.strategy.default1,summary.print=F)
  
  fit1<-cbind.data.frame(study=study, wr=res$Win_statistic$Win_Ratio[1], "wrul"=res$Win_statistic$Win_Ratio[3], "wrll"=res$Win_statistic$Win_Ratio[2], 
                         nb=res$Win_statistic$Net_Benefit[1], "nbul"=res$Win_statistic$Net_Benefit[3], "nbll"=res$Win_statistic$Net_Benefit[2], 
                         wo=res$Win_statistic$Win_Odds[1], "woul"=res$Win_statistic$Win_Odds[3], "woll"=res$Win_statistic$Win_Odds[2])
  fit<-rbind.data.frame(fit,fit1)
}
fit<-fit[-1,]
fit_i<-fit

fit<-cbind.data.frame(study=NA, wr=NA, "wrul"=NA, "wrll"=NA, nb=NA, "nbul"=NA, "nbll"=NA, wo=NA, "woul"=NA, "woll"=NA)
for (study in unique(win2[win2$eti=="ctd"&!is.na(win2$eti_code),]$studyid)[-2]) {
  res<-win.stat(data = win2[win2$eti=="ctd"&!is.na(win2$eti_code)&win2$studyid==study,], ep_type = "tte",arm.name = c("Active","Control"), tau = 0, priority = c(1:5), alpha = 0.05, digit = 5,
                weight = "unstratified", censoring_adjust ="IPCW",pvalue = "two-sided",win.strategy = win.strategy.default1,summary.print=F)
  
  fit1<-cbind.data.frame(study=study, wr=res$Win_statistic$Win_Ratio[1], "wrul"=res$Win_statistic$Win_Ratio[3], "wrll"=res$Win_statistic$Win_Ratio[2], 
                         nb=res$Win_statistic$Net_Benefit[1], "nbul"=res$Win_statistic$Net_Benefit[3], "nbll"=res$Win_statistic$Net_Benefit[2], 
                         wo=res$Win_statistic$Win_Odds[1], "woul"=res$Win_statistic$Win_Odds[3], "woll"=res$Win_statistic$Win_Odds[2])
  fit<-rbind.data.frame(fit,fit1)
}
fit<-fit[-1,]
fit_c<-fit

fit<-cbind.data.frame(study=NA, wr=NA, "wrul"=NA, "wrll"=NA, nb=NA, "nbul"=NA, "nbll"=NA, wo=NA, "woul"=NA, "woll"=NA)
for (study in unique(win2[win2$eti=="others"&!is.na(win2$eti_code),]$studyid)) {
  res<-win.stat(data = win2[win2$eti=="others"&!is.na(win2$eti_code)&win2$studyid==study,], ep_type = "tte",arm.name = c("Active","Control"), tau = 0, priority = c(1:5), alpha = 0.05, digit = 5,
                weight = "unstratified", censoring_adjust ="IPCW",pvalue = "two-sided",win.strategy = win.strategy.default1,summary.print=F)
  
  fit1<-cbind.data.frame(study=study, wr=res$Win_statistic$Win_Ratio[1], "wrul"=res$Win_statistic$Win_Ratio[3], "wrll"=res$Win_statistic$Win_Ratio[2], 
                         nb=res$Win_statistic$Net_Benefit[1], "nbul"=res$Win_statistic$Net_Benefit[3], "nbll"=res$Win_statistic$Net_Benefit[2], 
                         wo=res$Win_statistic$Win_Odds[1], "woul"=res$Win_statistic$Win_Odds[3], "woll"=res$Win_statistic$Win_Odds[2])
  fit<-rbind.data.frame(fit,fit1)
}
fit<-fit[-1,]
fit_o<-fit


#### second stage
## Inverse variance meta-analysis
fit_i$eti<-"idio"
fit_c$eti<-"ctd"
fit_o$eti<-"others"
fit<-rbind.data.frame(fit_i, fit_c,fit_o)

#WR
met<-metagen(TE=wr, lower=`wrll`, upper = `wrul`,sm="OR",
             comb.fixed=F,hakn=T,transf=F,subgroup=fit$eti,test.subgroup=T, data=fit)
forest(met,
       leftcols = c("studlab","effect.ci"),rightcols = "w.random",
       leftlabs = c("Trial", NA),just = "center",title="yu",
       smlab="Win ratio",text.random="Total",digits=2,lwd=1.5,
       col.square="lightblue",col.diamond.fixed="darkblue",test.overall=T)

#WO
met<-metagen(TE=wo, lower=`woll`, upper = `woul`,sm="OR",
             comb.fixed=F,hakn=T,transf=F,subgroup=fit$eti,test.subgroup=T, data=fit)
forest(met,
       leftcols = c("studlab","effect.ci"),rightcols = "w.random",
       leftlabs = c("Trial", NA),just = "center",title="yu",
       smlab="Win odds",text.random="Total",digits=2,lwd=1.5,
       col.square="lightblue",col.diamond.fixed="darkblue",test.overall=T)

#NB
met<-metagen(TE=nb, lower=`nbll`, upper = `nbul`,sm="RD",
             comb.fixed=T,hakn=T,transf=F,subgroup=fit$eti,test.subgroup=T, data=fit,control=list(stepadj=0.5, maxiter=10000))
forest(met,
       leftcols = c("studlab","effect.ci"),rightcols = "w.random",
       leftlabs = c("Trial", NA),just = "center",title="yu",
       smlab="Win odds",text.random="Total",digits=2,lwd=1.5,
       col.square="lightblue",col.diamond.fixed="darkblue",test.overall=T)


##### 6.4 by baseline 6MWD
win2$wd<-ifelse(win2$walk_dist>=440,">=440","<440")
table(win2$wd)
table(win2$wd)/sum(table(win2$wd))

#### first stage
fit<-cbind.data.frame(study=NA, wr=NA, "wrul"=NA, "wrll"=NA, nb=NA, "nbul"=NA, "nbll"=NA, wo=NA, "woul"=NA, "woll"=NA)
for (study in unique(win2[win2$wd=="<440"&!is.na(win2$walk_dist),]$studyid)) {
  res<-win.stat(data = win2[win2$wd=="<440"&!is.na(win2$walk_dist)&win2$studyid==study,], ep_type = "tte",arm.name = c("Active","Control"), tau = 0, priority = c(1:5), alpha = 0.05, digit = 5,
                weight = "unstratified", censoring_adjust ="IPCW",pvalue = "two-sided",win.strategy = win.strategy.default1,summary.print=F)
  
  fit1<-cbind.data.frame(study=study, wr=res$Win_statistic$Win_Ratio[1], "wrul"=res$Win_statistic$Win_Ratio[3], "wrll"=res$Win_statistic$Win_Ratio[2], 
                         nb=res$Win_statistic$Net_Benefit[1], "nbul"=res$Win_statistic$Net_Benefit[3], "nbll"=res$Win_statistic$Net_Benefit[2], 
                         wo=res$Win_statistic$Win_Odds[1], "woul"=res$Win_statistic$Win_Odds[3], "woll"=res$Win_statistic$Win_Odds[2])
  fit<-rbind.data.frame(fit,fit1)
}
fit<-fit[-1,]
fit_l<-fit

fit<-cbind.data.frame(study=NA, wr=NA, "wrul"=NA, "wrll"=NA, nb=NA, "nbul"=NA, "nbll"=NA, wo=NA, "woul"=NA, "woll"=NA)
for (study in unique(win2[win2$wd==">=440"&!is.na(win2$walk_dist),]$studyid)) {
  res<-win.stat(data = win2[win2$wd==">=440"&!is.na(win2$walk_dist)&win2$studyid==study,], ep_type = "tte",arm.name = c("Active","Control"), tau = 0, priority = c(1:5), alpha = 0.05, digit = 5,
                weight = "unstratified", censoring_adjust ="IPCW",pvalue = "two-sided",win.strategy = win.strategy.default1,summary.print=F)
  
  fit1<-cbind.data.frame(study=study, wr=res$Win_statistic$Win_Ratio[1], "wrul"=res$Win_statistic$Win_Ratio[3], "wrll"=res$Win_statistic$Win_Ratio[2], 
                         nb=res$Win_statistic$Net_Benefit[1], "nbul"=res$Win_statistic$Net_Benefit[3], "nbll"=res$Win_statistic$Net_Benefit[2], 
                         wo=res$Win_statistic$Win_Odds[1], "woul"=res$Win_statistic$Win_Odds[3], "woll"=res$Win_statistic$Win_Odds[2])
  fit<-rbind.data.frame(fit,fit1)
}
fit<-fit[-1,]
fit_u<-fit

#### second stage
## Inverse variance meta-analysis
fit_l$wd<-"<440"
fit_u$wd<-">=440"
fit<-rbind.data.frame(fit_l, fit_u)

#WR
met<-metagen(TE=wr, lower=`wrll`, upper = `wrul`,sm="OR",
             comb.fixed=F,hakn=T,transf=F,subgroup=fit$wd,test.subgroup=T, data=fit)
forest(met,
       leftcols = c("studlab","effect.ci"),rightcols = "w.random",
       leftlabs = c("Trial", NA),just = "center",title="yu",
       smlab="Win ratio",text.random="Total",digits=2,lwd=1.5,
       col.square="lightblue",col.diamond.fixed="darkblue",test.overall=T)

#WO
met<-metagen(TE=wo, lower=`woll`, upper = `woul`,sm="OR",
             comb.fixed=F,hakn=T,transf=F,subgroup=fit$wd,test.subgroup=T, data=fit)
forest(met,
       leftcols = c("studlab","effect.ci"),rightcols = "w.random",
       leftlabs = c("Trial", NA),just = "center",title="yu",
       smlab="Win odds",text.random="Total",digits=2,lwd=1.5,
       col.square="lightblue",col.diamond.fixed="darkblue",test.overall=T)

#NB
met<-metagen(TE=nb, lower=`nbll`, upper = `nbul`,sm="RD",
             comb.fixed=F,hakn=T,transf=F,subgroup=fit$wd,test.subgroup=T, data=fit)
forest(met,
       leftcols = c("studlab","effect.ci"),rightcols = "w.random",
       leftlabs = c("Trial", NA),just = "center",title="yu",
       smlab="Win odds",text.random="Total",digits=2,lwd=1.5,
       col.square="lightblue",col.diamond.fixed="darkblue",test.overall=T)


##### 6.5 by baseline WHO FC
win2$fc1<-ifelse(win2$fc<=2,"2","3")
table(win2$fc1)
table(win2$fc1)/sum(table(win2$fc1))

#### first stage
fit<-cbind.data.frame(study=NA, wr=NA, "wrul"=NA, "wrll"=NA, nb=NA, "nbul"=NA, "nbll"=NA, wo=NA, "woul"=NA, "woll"=NA)
for (study in unique(win2[win2$fc1=="2"&!is.na(win2$fc),]$studyid)) {
  res<-win.stat(data = win2[win2$fc1=="2"&!is.na(win2$fc)&win2$studyid==study,], ep_type = "tte",arm.name = c("Active","Control"), tau = 0, priority = c(1:5), alpha = 0.05, digit = 5,
                weight = "unstratified", censoring_adjust ="IPCW",pvalue = "two-sided",win.strategy = win.strategy.default1,summary.print=F)
  
  fit1<-cbind.data.frame(study=study, wr=res$Win_statistic$Win_Ratio[1], "wrul"=res$Win_statistic$Win_Ratio[3], "wrll"=res$Win_statistic$Win_Ratio[2], 
                         nb=res$Win_statistic$Net_Benefit[1], "nbul"=res$Win_statistic$Net_Benefit[3], "nbll"=res$Win_statistic$Net_Benefit[2], 
                         wo=res$Win_statistic$Win_Odds[1], "woul"=res$Win_statistic$Win_Odds[3], "woll"=res$Win_statistic$Win_Odds[2])
  fit<-rbind.data.frame(fit,fit1)
}
fit<-fit[-1,]
fit_l<-fit

fit<-cbind.data.frame(study=NA, wr=NA, "wrul"=NA, "wrll"=NA, nb=NA, "nbul"=NA, "nbll"=NA, wo=NA, "woul"=NA, "woll"=NA)
for (study in unique(win2[win2$fc1=="3"&!is.na(win2$fc),]$studyid)) {
  res<-win.stat(data = win2[win2$fc1=="3"&!is.na(win2$fc)&win2$studyid==study,], ep_type = "tte",arm.name = c("Active","Control"), tau = 0, priority = c(1:5), alpha = 0.05, digit = 5,
                weight = "unstratified", censoring_adjust ="IPCW",pvalue = "two-sided",win.strategy = win.strategy.default1,summary.print=F)
  
  fit1<-cbind.data.frame(study=study, wr=res$Win_statistic$Win_Ratio[1], "wrul"=res$Win_statistic$Win_Ratio[3], "wrll"=res$Win_statistic$Win_Ratio[2], 
                         nb=res$Win_statistic$Net_Benefit[1], "nbul"=res$Win_statistic$Net_Benefit[3], "nbll"=res$Win_statistic$Net_Benefit[2], 
                         wo=res$Win_statistic$Win_Odds[1], "woul"=res$Win_statistic$Win_Odds[3], "woll"=res$Win_statistic$Win_Odds[2])
  fit<-rbind.data.frame(fit,fit1)
}
fit<-fit[-1,]
fit_u<-fit

#### second stage
## Inverse variance meta-analysis
fit_l$fc<-"<=2"
fit_u$fc<-">=3"
fit<-rbind.data.frame(fit_l, fit_u)

#WR
met<-metagen(TE=wr, lower=`wrll`, upper = `wrul`,sm="OR",
             comb.fixed=F,hakn=T,transf=F,subgroup=fit$fc,test.subgroup=T, data=fit)
forest(met,
       leftcols = c("studlab","effect.ci"),rightcols = "w.random",
       leftlabs = c("Trial", NA),just = "center",title="yu",
       smlab="Win ratio",text.random="Total",digits=2,lwd=1.5,
       col.square="lightblue",col.diamond.fixed="darkblue",test.overall=T)

#WO
met<-metagen(TE=wo, lower=`woll`, upper = `woul`,sm="OR",
             comb.fixed=F,hakn=T,transf=F,subgroup=fit$fc,test.subgroup=T, data=fit)
forest(met,
       leftcols = c("studlab","effect.ci"),rightcols = "w.random",
       leftlabs = c("Trial", NA),just = "center",title="yu",
       smlab="Win odds",text.random="Total",digits=2,lwd=1.5,
       col.square="lightblue",col.diamond.fixed="darkblue",test.overall=T)

#NB
met<-metagen(TE=nb, lower=`nbll`, upper = `nbul`,sm="RD",
             comb.fixed=F,hakn=T,transf=F,subgroup=fit$fc,test.subgroup=T, data=fit)
forest(met,
       leftcols = c("studlab","effect.ci"),rightcols = "w.random",
       leftlabs = c("Trial", NA),just = "center",title="yu",
       smlab="Win odds",text.random="Total",digits=2,lwd=1.5,
       col.square="lightblue",col.diamond.fixed="darkblue",test.overall=T)



##### 6.6 by background PAH therapy
table(win2$background_therapy)
table(win2$background_therapy)/sum(table(win2$background_therapy))

#### first stage
fit<-cbind.data.frame(study=NA, wr=NA, "wrul"=NA, "wrll"=NA, nb=NA, "nbul"=NA, "nbll"=NA, wo=NA, "woul"=NA, "woll"=NA)
for (study in unique(win2[win2$background_therapy=="0",]$studyid)) {
  res<-win.stat(data = win2[win2$background_therapy=="0"&win2$studyid==study,], ep_type = "tte",arm.name = c("Active","Control"), tau = 0, priority = c(1:5), alpha = 0.05, digit = 5,
                weight = "unstratified", censoring_adjust ="IPCW",pvalue = "two-sided",win.strategy = win.strategy.default1,summary.print=F)
  
  fit1<-cbind.data.frame(study=study, wr=res$Win_statistic$Win_Ratio[1], "wrul"=res$Win_statistic$Win_Ratio[3], "wrll"=res$Win_statistic$Win_Ratio[2], 
                         nb=res$Win_statistic$Net_Benefit[1], "nbul"=res$Win_statistic$Net_Benefit[3], "nbll"=res$Win_statistic$Net_Benefit[2], 
                         wo=res$Win_statistic$Win_Odds[1], "woul"=res$Win_statistic$Win_Odds[3], "woll"=res$Win_statistic$Win_Odds[2])
  fit<-rbind.data.frame(fit,fit1)
}
fit<-fit[-1,]
fit_0<-fit

fit<-cbind.data.frame(study=NA, wr=NA, "wrul"=NA, "wrll"=NA, nb=NA, "nbul"=NA, "nbll"=NA, wo=NA, "woul"=NA, "woll"=NA)
for (study in unique(win2[win2$background_therapy=="1",]$studyid)) {
  res<-win.stat(data = win2[win2$background_therapy=="1"&win2$studyid==study,], ep_type = "tte",arm.name = c("Active","Control"), tau = 0, priority = c(1:5), alpha = 0.05, digit = 5,
                weight = "unstratified", censoring_adjust ="No",pvalue = "two-sided",win.strategy = win.strategy.default1,summary.print=F)
  
  fit1<-cbind.data.frame(study=study, wr=res$Win_statistic$Win_Ratio[1], "wrul"=res$Win_statistic$Win_Ratio[3], "wrll"=res$Win_statistic$Win_Ratio[2], 
                         nb=res$Win_statistic$Net_Benefit[1], "nbul"=res$Win_statistic$Net_Benefit[3], "nbll"=res$Win_statistic$Net_Benefit[2], 
                         wo=res$Win_statistic$Win_Odds[1], "woul"=res$Win_statistic$Win_Odds[3], "woll"=res$Win_statistic$Win_Odds[2])
  fit<-rbind.data.frame(fit,fit1)
}
fit<-fit[-1,]
fit_1<-fit


fit<-cbind.data.frame(study=NA, wr=NA, "wrul"=NA, "wrll"=NA, nb=NA, "nbul"=NA, "nbll"=NA, wo=NA, "woul"=NA, "woll"=NA)
for (study in c(1010, 1013, 1026, 1028)) {
  res<-win.stat(data = win2[win2$background_therapy=="2"&win2$studyid==study,], ep_type = "tte",arm.name = c("Active","Control"), tau = 0, priority = c(1:5), alpha = 0.05, digit = 5,
                weight = "unstratified", censoring_adjust ="No",pvalue = "two-sided",win.strategy = win.strategy.default1,summary.print=F)
  
  fit1<-cbind.data.frame(study=study, wr=res$Win_statistic$Win_Ratio[1], "wrul"=res$Win_statistic$Win_Ratio[3], "wrll"=res$Win_statistic$Win_Ratio[2], 
                         nb=res$Win_statistic$Net_Benefit[1], "nbul"=res$Win_statistic$Net_Benefit[3], "nbll"=res$Win_statistic$Net_Benefit[2], 
                         wo=res$Win_statistic$Win_Odds[1], "woul"=res$Win_statistic$Win_Odds[3], "woll"=res$Win_statistic$Win_Odds[2])
  fit<-rbind.data.frame(fit,fit1)
}
fit<-fit[-1,]
fit_2<-fit


fit_0$BacgroundTherapy<-"0"
fit_1$BacgroundTherapy<-"1"
fit_2$BacgroundTherapy<-"2"
fit<-rbind.data.frame(fit_0, fit_1, fit_2)

#### second stage
## Inverse variance meta-analysis

#WR
met<-metagen(TE=wr, lower=`wrll`, upper = `wrul`,sm="OR",
             comb.fixed=F,hakn=T,transf=F,subgroup=fit$BacgroundTherapy,test.subgroup=T, data=fit)
forest(met,
       leftcols = c("studlab","effect.ci"),rightcols = "w.random",
       leftlabs = c("Trial", NA),just = "center",title="yu",
       smlab="Win ratio",text.random="Total",digits=2,lwd=1.5,
       col.square="lightblue",col.diamond.fixed="darkblue",test.overall=T)

#WO
met<-metagen(TE=wo,lower=`woll`, upper = `woul`,sm="OR",
             comb.fixed=F,hakn=T,transf=F,subgroup=fit$BacgroundTherapy,subgroup.name="Background therapy",test.subgroup=T, data=fit)
forest(met,
       leftcols = c("studlab","effect.ci"),rightcols = "w.random",
       leftlabs = c("Trial", NA),just = "center",title="yu",
       smlab="Win odds",text.random="Total",digits=2,lwd=1.5,
       col.square="lightblue",col.diamond.fixed="darkblue",test.overall=T)

#NB
met<-metagen(TE=nb, studlab = study1, lower=`nbll`, upper = `nbul`,sm="RD",
             comb.fixed=F,hakn=T,transf=F,subgroup=fit$BacgroundTherapy,subgroup.name="Background therapy",test.subgroup=T,  data=fit)
forest(met,
       leftcols = c("studlab","effect.ci"),rightcols = "w.random",
       leftlabs = c("Trial", NA),just = "center",title="yu",
       smlab="Net benefit",text.random="Total",digits=2,lwd=1.5,
       col.square="lightblue",col.diamond.fixed="darkblue",test.overall=T)



###### 7.0 median(IQR) for the delay between events for a pair that both had events
win3=win2

##create pairs (unstratified)
win3$m=1
dat<-merge(win3[win3$arm=="Active",],win3[win3$arm=="Control",],by="m")
nrow(dat)

## For each outcome event, determine wins (1), ties(0), and losses(-1) [dat$res1]
## For pairs who both had an event, compute the time difference between events (in weeks) (dat$td1)
dat$res1<-ifelse(dat$Delta_1.x==1&dat$Delta_1.y==1 & dat$Y_1.x > (dat$Y_1.y + 0)|
                   dat$Delta_1.x<1&dat$Delta_1.y==1 & dat$Y_1.y < (dat$Y_1.x + 0),1,
                 ifelse(dat$Delta_1.x==1&dat$Delta_1.y==1 & dat$Y_1.x < (dat$Y_1.y + 0)|
                          dat$Delta_1.x==1&dat$Delta_1.y<1 & dat$Y_1.x < (dat$Y_1.y + 0),-1,0))
dat$td1=ifelse(dat$Delta_1.x==1&dat$Delta_1.y==1, dat$Y_1.x-dat$Y_1.y, NA)
table(dat$res1)

dat$res2<-ifelse(dat$Delta_2.x==1&dat$Delta_2.y==1 & dat$Y_2.x > (dat$Y_2.y + 0)|
                   dat$Delta_2.x<1&dat$Delta_2.y==1 & dat$Y_2.y < (dat$Y_2.x + 0),1,
                 ifelse(dat$Delta_2.x==1&dat$Delta_2.y==1 & dat$Y_2.x < (dat$Y_2.y + 0)|
                          dat$Delta_2.x==1&dat$Delta_2.y<1 & dat$Y_2.x < (dat$Y_2.y + 0),-1,0))
dat$res2<-ifelse(dat$res1==0,dat$res2,NA)
dat$td2=ifelse(dat$Delta_2.x==1&dat$Delta_2.y==1&dat$res1==0, dat$Y_2.x-dat$Y_2.y, NA)
table(dat$res2)

dat$res3<-ifelse(dat$Delta_3.x==1&dat$Delta_3.y==1 & dat$Y_3.x > (dat$Y_3.y + 0)|
                   dat$Delta_3.x<1&dat$Delta_3.y==1 & dat$Y_3.y < (dat$Y_3.x + 0),1,
                 ifelse(dat$Delta_3.x==1&dat$Delta_3.y==1 & dat$Y_3.x < (dat$Y_3.y + 0)|
                          dat$Delta_3.x==1&dat$Delta_3.y<1 & dat$Y_3.x < (dat$Y_3.y + 0),-1,0))
dat$res3<-ifelse(dat$res2==0,dat$res3,NA)
dat$td3=ifelse(dat$Delta_3.x==1&dat$Delta_3.y==1&dat$res2==0, dat$Y_3.x-dat$Y_3.y, NA)
table(dat$res3)

dat$res4<-ifelse(dat$Delta_4.x==1&dat$Delta_4.y==1 & dat$Y_4.x < (dat$Y_4.y + 0)|
                   dat$Delta_4.x==1&dat$Delta_4.y<1 & dat$Y_4.x < (dat$Y_4.y + 0),1,
                 ifelse(dat$Delta_4.x==1&dat$Delta_4.y==1 & dat$Y_4.x > (dat$Y_4.y + 0)|
                          dat$Delta_4.x<1&dat$Delta_4.y==1 & dat$Y_4.y < (dat$Y_4.x + 0),-1,0))
dat$res4<-ifelse(dat$res3==0,dat$res4,NA)
dat$td4=ifelse(dat$Delta_4.x==1&dat$Delta_4.y==1&dat$res3==0, dat$Y_4.x-dat$Y_4.y, NA)
table(dat$res4)

dat$res5<-ifelse(dat$Delta_5.x==1&dat$Delta_5.y==1 & dat$Y_5.x > (dat$Y_5.y + 0)|
                   dat$Delta_5.x<1&dat$Delta_5.y==1 & dat$Y_5.y < (dat$Y_5.x + 0),1,
                 ifelse(dat$Delta_5.x==1&dat$Delta_5.y==1 & dat$Y_5.x < (dat$Y_5.y + 0)|
                          dat$Delta_5.x==1&dat$Delta_5.y<1 & dat$Y_5.x < (dat$Y_5.y + 0),-1,0))
dat$res5<-ifelse(dat$res4==0,dat$res5,NA)
dat$td5=ifelse(dat$Delta_5.x==1&dat$Delta_5.y==1&dat$res4==0, dat$Y_5.x-dat$Y_5.y, NA)
table(dat$res5)

## For each outcome event, compute the min, p25, p50, p75, and max time difference between events where treatment wins (1), treatment loses (-1), and the absolute time difference (1 or -1)
quantile(dat[dat$res1==1,]$td1, probs=c(0,0.25,0.50,0.75,1),na.rm=T)
quantile(dat[dat$res1==-1,]$td1, probs=c(0,0.25,0.50,0.75,1),na.rm=T)
quantile(abs(dat$td1), probs=c(0,0.25,0.50,0.75,1),na.rm=T)

quantile(dat[dat$res2==1,]$td2, probs=c(0,0.25,0.50,0.75,1),na.rm=T)
quantile(dat[dat$res2==-1,]$td2, probs=c(0,0.25,0.50,0.75,1),na.rm=T)
quantile(abs(dat$td2), probs=c(0,0.25,0.50,0.75,1),na.rm=T)

quantile(dat[dat$res3==1,]$td3, probs=c(0,0.25,0.50,0.75,1),na.rm=T)
quantile(dat[dat$res3==-1,]$td3, probs=c(0,0.25,0.50,0.75,1),na.rm=T)
quantile(abs(dat$td3), probs=c(0,0.25,0.50,0.75,1),na.rm=T)

quantile(dat[dat$res4==1,]$td4, probs=c(0,0.25,0.50,0.75,1),na.rm=T)
quantile(dat[dat$res4==-1,]$td4, probs=c(0,0.25,0.50,0.75,1),na.rm=T)
quantile(abs(dat$td4), probs=c(0,0.25,0.50,0.75,1),na.rm=T)


quantile(dat[dat$res5==1,]$td5, probs=c(0,0.25,0.50,0.75,1),na.rm=T)
quantile(dat[dat$res5==-1,]$td5, probs=c(0,0.25,0.50,0.75,1),na.rm=T)
quantile(abs(dat$td5), probs=c(0,0.25,0.50,0.75,1),na.rm=T)



###### 8.0 Sensitivity analysis I: Stratify analysis by baseline WHO FC to account for PAH severity
win2$fc1<-ifelse(win2$fc<=2,2,3)

win2$stratum<-win2$fc1
res<-win.stat(data = win2, ep_type = "tte",arm.name = c("Active","Control"), tau = 0, priority = c(1:5), alpha = 0.05, digit = 5,
              weight = "MH-type", censoring_adjust ="IPCW",pvalue = "two-sided",win.strategy = win.strategy.default1,summary.print=F)

res_all=res
res_all$Win_statistic

#Win statistics by study
# first step
fit<-cbind.data.frame(study=NA, wr=NA, "wrul"=NA, "wrll"=NA, nb=NA, "nbul"=NA, "nbll"=NA, wo=NA, "woul"=NA, "woll"=NA)
for (study in unique(win2$studyid)) {
  res<-win.stat(data = win2[win2$studyid==study,], ep_type = "tte",arm.name = c("Active","Control"), tau = 0, priority = c(1:5), alpha = 0.05, digit = 5,
                weight = "MH-type", censoring_adjust ="IPCW",pvalue = "two-sided",win.strategy = win.strategy.default1,summary.print=F)
  
  fit1<-cbind.data.frame(study=study, wr=res$Win_statistic$Win_Ratio[1], "wrul"=res$Win_statistic$Win_Ratio[3], "wrll"=res$Win_statistic$Win_Ratio[2], 
                         nb=res$Win_statistic$Net_Benefit[1], "nbul"=res$Win_statistic$Net_Benefit[3], "nbll"=res$Win_statistic$Net_Benefit[2], 
                         wo=res$Win_statistic$Win_Odds[1], "woul"=res$Win_statistic$Win_Odds[3], "woll"=res$Win_statistic$Win_Odds[2])
  fit<-rbind.data.frame(fit,fit1)
}
fit<-fit[-1,]

fit_all<-fit

## second step
# Inverse variance meta-analysis
fit <- fit[order(fit$study),]
met<-metagen(TE=wr, lower=`wrll`, upper = `wrul`,sm="OR",
             comb.fixed=F,hakn=T,transf=F, data=fit)
forest(met,
       leftcols = c("studlab","effect.ci"),rightcols = "w.random",
       leftlabs = c("Trial", NA),just = "center",title="yu",
       smlab="Win ratio",text.random="Total",digits=2,lwd=1.5,
       col.square="lightblue",col.diamond.fixed="darkblue",test.overall=T)


met<-metagen(TE=wo, lower=`woll`, upper = `woul`,sm="OR",
             comb.fixed=F,hakn=T,transf=F, data=fit)
forest(met,
       leftcols = c("studlab","effect.ci"),rightcols = "w.random",
       leftlabs = c("Trial", NA),just = "center",title="yu",
       smlab="Win odds",text.random="Total",digits=2,lwd=1.5,
       col.square="lightblue",col.diamond.fixed="darkblue",test.overall=T)


met<-metagen(TE=nb, lower=`nbll`, upper = `nbul`,sm="RD",
             comb.fixed=F,hakn=T,transf=F, data=fit)
forest(met,
       leftcols = c("studlab","effect.ci"),rightcols = "w.random",
       leftlabs = c("Trial", NA),just = "center",title="yu",
       smlab="Net benefit",text.random="Total",digits=2,lwd=1.5,
       col.square="lightblue",col.diamond.fixed="darkblue",test.overall=T)






###### 9.0 Sensitivity analysis I: Rank clinical worsening before clinical improvement in the event hierarchy 
#deaths>lung transplant>hospitalization due to PAH>clinical worsening>clinical improvement
win1<-win
win1$arm<-win1$study_arm_assignment
win1$id<-win1$usubjid
win1$Delta_1<-win1$cedeath
win1$Delta_2<-win1$cetxp
win1$Delta_3<-win1$cehosp
win1$Delta_4<-win1$cedecline
win1$Delta_5<-win1$ceimprove
win1$Y_1<-win1$weeknumdeath
win1$Y_2<-win1$weeknumtxp
win1$Y_3<-win1$weeknumhosp
win1$Y_4<-win1$weeknumdecline
win1$Y_5<-win1$weeknumimprove

win2=win1

######set vector for the improvement outcome(s)
#the only improvement outcome is clinical improvement is now outcome 5 in our hierarchy 
imp=5



###### 9.1 Compute unstratified win statistics without IPCW
res<-win.stat(data = win2, ep_type = "tte",arm.name = c("Active","Control"), tau = 0, priority = c(1:5), alpha = 0.05, digit = 5,
              weight = "unstratified", censoring_adjust ="No",pvalue = "two-sided",win.strategy = win.strategy.default1,summary.print=F)


###### 9.1.1 Overall number of wins, losses, and ties
res_all=res
res_all$Win_statistic
##wins
res_all$Win_prop[,2]*nrow(win2[win2$arm=="Active",])*nrow(win2[win2$arm=="Control",])
##losses
res_all$Win_prop[,3]*nrow(win2[win2$arm=="Active",])*nrow(win2[win2$arm=="Control",])
##ties
(1-(res_all$Win_prop[,2]+res_all$Win_prop[,3]))*nrow(win2[win2$arm=="Active",])*nrow(win2[win2$arm=="Control",])

###### 9.1.2 Number of wins, losses, and ties per outcome event (1 to 5 respectively)
##Wins
res$summary_ep$Trt_Endpoint1$Count
res$summary_ep$Trt_Endpoint2$Count
res$summary_ep$Trt_Endpoint3$Count
res$summary_ep$Trt_Endpoint4$Count
res$summary_ep$Trt_Endpoint5$Count
##losses
res$summary_ep$Con_Endpoint1$Count
res$summary_ep$Con_Endpoint2$Count
res$summary_ep$Con_Endpoint3$Count
res$summary_ep$Con_Endpoint4$Count
res$summary_ep$Con_Endpoint5$Count
##ties
T1=nrow(win2[win2$arm=="Active",])*nrow(win2[win2$arm=="Control",])-(res$summary_ep$Trt_Endpoint1$Count+res$summary_ep$Con_Endpoint1$Count);T1
T2=T1-(res$summary_ep$Trt_Endpoint2$Count+res$summary_ep$Con_Endpoint2$Count);T2
T3=T2-(res$summary_ep$Trt_Endpoint3$Count+res$summary_ep$Con_Endpoint3$Count);T3
T4=T3-(res$summary_ep$Trt_Endpoint4$Count+res$summary_ep$Con_Endpoint4$Count);T4
T5=T4-(res$summary_ep$Trt_Endpoint5$Count+res$summary_ep$Con_Endpoint5$Count);T5
##Win ratio, net benefit, and win odds per outcome event
WR1=res$summary_ep$Trt_Endpoint1$Count/res$summary_ep$Con_Endpoint1$Count;WR1
NB1=res$summary_ep$Trt_Endpoint1$Proportion-res$summary_ep$Con_Endpoint1$Proportion;NB1
WO1=(res$summary_ep$Trt_Endpoint1$Count+(0.5*T1))/(res$summary_ep$Con_Endpoint1$Count+(0.5*T1));WO1
WO1=(res$summary_ep$Trt_Endpoint1$Proportion+(0.5*(1-(res$summary_ep$Trt_Endpoint1$Proportion+res$summary_ep$Con_Endpoint1$Proportion))))/(res$summary_ep$Con_Endpoint1$Proportion+(0.5*(1-(res$summary_ep$Trt_Endpoint1$Proportion+res$summary_ep$Con_Endpoint1$Proportion))));WO1
WR2=res$summary_ep$Trt_Endpoint2$Count/res$summary_ep$Con_Endpoint2$Count;WR2
NB2=res$summary_ep$Trt_Endpoint2$Proportion-res$summary_ep$Con_Endpoint2$Proportion;NB2
WO2=(res$summary_ep$Trt_Endpoint2$Count+(0.5*T2))/(res$summary_ep$Con_Endpoint2$Count+(0.5*T2));WO2
WR3=res$summary_ep$Trt_Endpoint3$Count/res$summary_ep$Con_Endpoint3$Count;WR3
NB3=res$summary_ep$Trt_Endpoint3$Proportion-res$summary_ep$Con_Endpoint3$Proportion;NB3
WO3=(res$summary_ep$Trt_Endpoint3$Count+(0.5*T3))/(res$summary_ep$Con_Endpoint3$Count+(0.5*T3));WO3
WR4=res$summary_ep$Trt_Endpoint4$Count/res$summary_ep$Con_Endpoint4$Count;WR4
NB4=res$summary_ep$Trt_Endpoint4$Proportion-res$summary_ep$Con_Endpoint4$Proportion;NB4
WO4=(res$summary_ep$Trt_Endpoint4$Count+(0.5*T4))/(res$summary_ep$Con_Endpoint4$Count+(0.5*T4));WO4
WR5=res$summary_ep$Trt_Endpoint5$Count/res$summary_ep$Con_Endpoint5$Count;WR5
NB5=res$summary_ep$Trt_Endpoint5$Proportion-res$summary_ep$Con_Endpoint5$Proportion;NB5
WO5=(res$summary_ep$Trt_Endpoint5$Count+(0.5*T5))/(res$summary_ep$Con_Endpoint5$Count+(0.5*T5));WO5




###### 9.2 Compute unstratified win statistics with IPCW
res<-win.stat(data = win2, ep_type = "tte",arm.name = c("Active","Control"), tau = 0, priority = c(1:5), alpha = 0.05, digit = 5,
              weight = "unstratified", censoring_adjust ="IPCW",pvalue = "two-sided",win.strategy = win.strategy.default1,summary.print=F)
res$Win_statistic


###### 9.3 Compute win statistics with IPCW stratified by trial
#Two stage: compute win statistics within each trial, the meta-analyze across trials

####First stage
#Win statistics by trial
fit<-cbind.data.frame(study=NA, wr=NA, "wrul"=NA, "wrll"=NA, nb=NA, "nbul"=NA, "nbll"=NA, wo=NA, "woul"=NA, "woll"=NA)
for (study in unique(win2$studyid)) {
  res<-win.stat(data = win2[win2$studyid==study,], ep_type = "tte",arm.name = c("Active","Control"), tau = 0, priority = c(1:5), alpha = 0.05, digit = 5,
                weight = "unstratified", censoring_adjust ="IPCW",pvalue = "two-sided",win.strategy = win.strategy.default1,summary.print=F)
  
  fit1<-cbind.data.frame(study=study, wr=res$Win_statistic$Win_Ratio[1], "wrul"=res$Win_statistic$Win_Ratio[3], "wrll"=res$Win_statistic$Win_Ratio[2], 
                         nb=res$Win_statistic$Net_Benefit[1], "nbul"=res$Win_statistic$Net_Benefit[3], "nbll"=res$Win_statistic$Net_Benefit[2], 
                         wo=res$Win_statistic$Win_Odds[1], "woul"=res$Win_statistic$Win_Odds[3], "woll"=res$Win_statistic$Win_Odds[2])
  fit<-rbind.data.frame(fit,fit1)
}
fit<-fit[-1,]

fit_all<-fit


#### Second stage
## Inverse variance meta-analysis
fit <- fit[order(fit$study),]

#WR
met<-metagen(TE=wr, lower=`wrll`, upper = `wrul`,sm="OR",
             comb.fixed=F,hakn=T,transf=F, data=fit)
forest(met,
       leftcols = c("studlab","effect.ci"),rightcols = "w.random",
       leftlabs = c("Trial", NA),just = "center",title="yu",
       smlab="Win ratio",text.random="Total",digits=2,lwd=1.5,
       col.square="lightblue",col.diamond.fixed="darkblue",test.overall=T)

#WO
met<-metagen(TE=wo, lower=`woll`, upper = `woul`,sm="OR",
             comb.fixed=F,hakn=T,transf=F, data=fit)
forest(met,
       leftcols = c("studlab","effect.ci"),rightcols = "w.random",
       leftlabs = c("Trial", NA),just = "center",title="yu",
       smlab="Win odds",text.random="Total",digits=2,lwd=1.5,
       col.square="lightblue",col.diamond.fixed="darkblue",test.overall=T)

#NB
met<-metagen(TE=nb, lower=`nbll`, upper = `nbul`,sm="RD",
             comb.fixed=F,hakn=T,transf=F, data=fit)
forest(met,
       leftcols = c("studlab","effect.ci"),rightcols = "w.random",
       leftlabs = c("Trial", NA),just = "center",title="yu",
       smlab="Net benefit",text.random="Total",digits=2,lwd=1.5,
       col.square="lightblue",col.diamond.fixed="darkblue",test.overall=T)





####### End of script