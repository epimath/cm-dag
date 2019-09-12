statistical_analysis<-function(SumRes, Init_states, times)
{

#Calculate all incident events from model run by disease state i.e., just mortality
incEvents<-SumRes[2:length(SumRes[,1]),10:17]-SumRes[1:(length(SumRes[,1])-1),10:17] #first element is mortality from t0 to t1, length t-1 


#calculate person-time from model run by disease state
ptime<-SumRes[1:(length(SumRes[,1])-1),2:9]- (SumRes[1:(length(SumRes[,1])-1),2:9]-
                                                  SumRes[2:length(SumRes[,1]),2:9])/2 #first element is PT of people from t0 to t1, length t-1 

#set up final simulated dataset 1 is 'yes' and 0 is 'no', the final dataset will have person time and deaths by compartment
diabetes = c(1,1,1,1,1,1,1,1)  
smoker =   c(0,0,0,0,1,1,1,1)    
obesity =  c(0,0,1,1,0,0,1,1)    
NW =       c(1,1,0,0,1,1,0,0)     
agegroup = c(1,2,1,2,1,2,1,2)    
deaths_group = colSums(incEvents)
pt_group = colMeans(ptime)
df_wt= data.frame(diabetes,smoker,obesity,NW,agegroup,deaths_group, pt_group)

#split data-set by smoking status strata
df_wtes=df_wt[which(df_wt$smoker==1),] #ever-smokers
df_wtns=df_wt[which(df_wt$smoker==0),] #never smokers

#caluculate mortality rate ratios
Ept<-(sum(df_wtes[df_wtes$NW==1,"deaths_group"])/sum(df_wtes[df_wtes$NW==1,"pt_group"]))/
  (sum(df_wtes[df_wtes$NW==0,"deaths_group"])/sum(df_wtes[df_wtes$NW==0,"pt_group"])) #ever smokers

Npt<-(sum(df_wtns[df_wtns$NW==1,"deaths_group"])/sum(df_wtns[df_wtns$NW==1,"pt_group"]))/
  (sum(df_wtns[df_wtns$NW==0,"deaths_group"])/sum(df_wtns[df_wtns$NW==0,"pt_group"]))  #Never smokers


#now externally standardize to the unexposed (obese) group to adjust for age and calculate mortality rate ratios

#ever-smokers
EptSens<-(((df_wtes[df_wtes$NW==0 & df_wtes$agegroup==1 ,"pt_group"])*(df_wtes[df_wtes$NW==1 & df_wtes$agegroup==1 ,"deaths_group"])/(df_wtes[df_wtes$NW==1 & df_wtes$agegroup==1 ,"pt_group"]))+
            ((df_wtes[df_wtes$NW==0 & df_wtes$agegroup==2 ,"pt_group"])*(df_wtes[df_wtes$NW==1 & df_wtes$agegroup==2 ,"deaths_group"])/(df_wtes[df_wtes$NW==1 & df_wtes$agegroup==2 ,"pt_group"])))/
  (df_wtes[df_wtes$NW==0 & df_wtes$agegroup==1 ,"deaths_group"]+df_wtes[df_wtes$NW==0 & df_wtes$agegroup==2 ,"deaths_group"])
             
#Never-smokers
NptSens<-(((df_wtns[df_wtns$NW==0 & df_wtns$agegroup==1 ,"pt_group"])*(df_wtns[df_wtns$NW==1 & df_wtns$agegroup==1 ,"deaths_group"])/(df_wtns[df_wtns$NW==1 & df_wtns$agegroup==1 ,"pt_group"]))+
            ((df_wtns[df_wtns$NW==0 & df_wtns$agegroup==2 ,"pt_group"])*(df_wtns[df_wtns$NW==1 & df_wtns$agegroup==2 ,"deaths_group"])/(df_wtns[df_wtns$NW==1 & df_wtns$agegroup==2 ,"pt_group"])))/
  (df_wtns[df_wtns$NW==0 & df_wtns$agegroup==1 ,"deaths_group"]+df_wtns[df_wtns$NW==0 & df_wtns$agegroup==2 ,"deaths_group"])


OPmodels<- list("Ept" = Ept, "Npt" = Npt, "NptSens"=NptSens,"EptSens"=EptSens )

return(OPmodels)

}