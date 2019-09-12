rm(list = ls())
#load necessary R packages
library(deSolve)
library(AER)
library(MASS)
library(pscl)
library(pracma)
library(lme4)
library(lhs)
library(plyr)
library(dplyr)
library(ggplot2)
library(RColorBrewer) 
set.seed(12)

# Setup -------------------------------------------------------------------

#set working directory
setwd()

#Read in R Code necessary for simulations
source('age_structured_model.R') 
source('statistical_analysis.R')

#determine the proportion of individuals in each age group
demogroups = c("y","o") #define age groups
StudySize = 1000000 #define total study population
CensusAS<-read.csv("Age2010Census.csv") #read in census data 
Age_wt = CensusAS[,3]/CensusAS[dim(CensusAS)[1],3] #percent of total population in each age group

#percent of individuals in by age group as per census data
agepropy=sum(Age_wt[which(CensusAS[,1]=="40 to 44 years")], Age_wt[which(CensusAS[,1]=="45 to 49 years")],  
             Age_wt[which(CensusAS[,1]=="50 to 54 years")], Age_wt[which(CensusAS[,1]=="55 to 59 years")])
agewpropo=sum(Age_wt[which(CensusAS[,1]=="60 to 64 years")], Age_wt[which(CensusAS[,1]=="65 to 69 years")],
           Age_wt[which(CensusAS[,1]=="70 to 74 years")])

#corresponding individuals in study
agewty= StudySize*agepropy/(agepropy+agewpropo)
agewto= StudySize*agewpropo/(agepropy+agewpropo)
agewts=c(agewty,agewto)

times <- seq(1,2,1/365) #set timesteps: here we are setting 1 year timesteps and breaking into day intervals.
numsamples = 10000 #number of different parameter sets/simulated studies
p=1

#ranges to sample from
paramlower = c(0.01,0.01, 0.01, 0.01, 0.01, 0.01, 0.01,0.01,0.01, 0.01,0, 0,1,0,0,0,0,0,0)
paramupper = c(0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.1,0.1,0.1,2,1,1,1,1,1,1)

paramRangeNames<-c("NWD to ODy", "NWD to ODo","OD to NWDy","OD to NWDo","NWD to NWDSy","NWD to NWDSo","NWDS to ODSy",
                   "NWDS to ODSo","ODS to NWDSy","baseline_mortality_rate","obese_mortality_add_on", "smoking_mortality_add_on", 
                   "age_related_mortality_scaling", "youngIC1","youngIC2","youngIC3","oldIC1","oldIC2","oldIC3") 
#see appendix section: "Initial Condition Calculations" for more informaiton on IC variables

#initialize objects needed to store model inputs and outputs 
index=c() #this will be used to store the runs which result in the obesity paradox
paramValues = matrix(NA,length(paramlower)-1,numsamples)
rownames(paramValues)<-c("NWD to ODy", "NWD to ODo","OD to NWDy","OD to NWDo","NWD to NWDSy","NWD to NWDSo","NWDS to ODSy",
                         "NWDS to ODSo","ODS to NWDSy","ODS to NWDSo","NWDy mortality", "NWDo mortality","ODy mortality", 
                         "ODo mortality","NWDSy mortality", "NWDSo mortality","ODSy mortality", 
                         "ODSo mortality") #'y' denotes young and 'o' denotes old

Init_states = matrix(NA,length(paramlower)-3,numsamples)
rownames(Init_states)<- c("dNWD_y", "dNWD_o", "dOD_y","dOD_o", "dNWDS_y", "dNWDS_o", "dODS_y","dODS_o",
                          "dNWD_yM","dNWD_oM", "dOD_yM", "dOD_oM","dNWDS_yM","dNWDS_oM", "dODS_yM", 
                          "dODS_oM") #'y' denotes young and 'o' denotes old

paramsample = randomLHS(numsamples,(length(paramlower)))
ysample = matrix(NA,length(times),numsamples)
rawLHS= matrix(NA,numsamples, length(paramlower)-6)
Ev_res<-c()
Ne_res<-c()
Npt<-c()
Ept<-c()
NptSens<-c()
EptSens<-c()

# Sample LHS and run model ------------------------------------------------
for(i in 1:numsamples)
{
  # re-scale the parameters to be in the preset upper and lower bounds
  ptemp = paramlower + (paramupper-paramlower)*paramsample[i,] 
  #rawLHS[i,]<-ptemp[1:13] #storing just the lhs'ed parameter values
  
  #determining the distribution of individuals across disease states, i.e., setting initial coniditions (Init_states)
  #see appendix section: "Initial Condition Calculations" for more informaiton on IC variables
  Init_statesy<- diff(c(0, sort(ptemp[14:16]), 1))*agewty
  Init_stateso <- diff(c(0, sort(ptemp[17:19]), 1))*agewto                    
  Init_states[,i]<- c(Init_statesy[1],Init_stateso[1],Init_statesy[2],Init_stateso[2],
                      Init_statesy[3],Init_stateso[3],Init_statesy[4],Init_stateso[4],0,0,0,0,0,0,0,0)                     
  
  #calculating mortlity rates and appending to parameter set
  pmort = c(ptemp[10],ptemp[10]*ptemp[13], #NWD y then o
            ptemp[10]+ptemp[11],(ptemp[10]+ptemp[11])*ptemp[13], #OD y then o 
            ptemp[10]+ptemp[12],(ptemp[10]+ptemp[12])*ptemp[13], #NWS y then 0
            ptemp[10]+ptemp[11]+ptemp[12], (ptemp[10]+ptemp[11]+ptemp[12])*ptemp[13])
  
  ptemp =c(ptemp[1:5], ptemp[5],ptemp[6:9],pmort)
  paramValues[,i] <- ptemp #storing all parameter values
  
  #run ODE model
  xtemp <- ode(Init_states[,i], times, age_structured_model, ptemp, method='ode45') #run the model
  
  #obtain mortality rate ratios
  PtEst= statistical_analysis(xtemp, Init_states[,i], times) 
  
  #save results
  Ept[i] = PtEst$Ept
  Npt[i] = PtEst$Npt
  EptSens[i] = PtEst$EptSens
  NptSens[i] = PtEst$NptSens
  
}


# Save data ---------------------------------------------------------------
Obesity_Paradox<-matrix(0,nrow=numsamples)
Obesity_Paradox[which(Npt<1 & Ept>1)]<-1
Obesity_Paradox_adj<-matrix(0,nrow=numsamples)
Obesity_Paradox_adj[which(NptSens<1 & NptSens>1)]<-1
model2Data<-data.frame(Npt,Ept,NptSens,NptSens,Obesity_Paradox,Obesity_Paradox_adj)   
opSensIndx<-which(NptSens<1 & EptSens>1)
colnames(model2Data)<-c("NeverSmoker", "EverSmoker", "NeverSmoker_Sens", "EverSmoker_Sens","Obesity_Paradox_Unadjusted","Obesity_Paradox_Adjusted")
write.csv(model2Data, "model2data.csv")
write.csv(distData, "model2.csv")

# Generate Results --------------------------------------------------------
length(which(Npt<1 & Ept>1)) #number in unadjusted analysis that acheived obesity paradox
length(which(NptSens<1 & EptSens>1)) #number that acheived obesity paradox after adjusting for age

#Distribution of simulated studies by ever smoker and never smoker mortality rate ratios 
df_rect=data.frame(xmin=1,xmax=3.5, ymin=0, ymax=1, r = "Obesity Paradox")
ggplot() + labs(title="Model 2: Age-Varying Mortality")  +
  xlab("Ever-smoker")+ylab("Never-smoker") +xlim(c(0,3.5)) +ylim(c(0,2))+
  theme(legend.position="none", text = element_text(size=20), axis.text.x = element_text(size=20),axis.title.x = element_text(size=25),
        axis.title.y = element_text(size=25), axis.text.y = element_text(size=20)) +geom_hline(yintercept = 1, colour="red", lty=50) +
  geom_vline(xintercept = 1, colour="red", lty=50) + geom_rect(data=df_rect,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill="black", alpha=0.25) + 
  geom_text(data=df_rect, aes(x=xmin+(xmax-xmin)-0.5, y=ymin+(ymax-ymin)-0.85,label=r ),alpha=1, size = 8) + 
  geom_point(data= model2Data, aes(x=EverSmoker, y=NeverSmoker, fill=""),alpha=0.3)
