age_structured_model<-function(t, states, params)
{


  # Initial states ----------------------------------------------------------
  NWD <-states[1:length(demogroups)]
  OD <-states[(length(demogroups)+1):(2*length(demogroups))]
  NWDS <-states[(2*length(demogroups)+1):(3*length(demogroups))]
  ODS <-states[(3*length(demogroups)+1):(4*length(demogroups))]
  NWD_M <- states[(4*length(demogroups)+1):(5*length(demogroups))]
  OD_M <- states[(5*length(demogroups)+1):(6*length(demogroups))]
  NWDS_M <- states[(6*length(demogroups)+1):(7*length(demogroups))]
  ODS_M <- states[(7*length(demogroups)+1):(8*length(demogroups))]
  
  # parameters --------------------------------------------------------------
  a<-params[1:length(demogroups)] # NWD to OD
  b<-params[(length(demogroups)+1):(2*length(demogroups))] #OD to NWD
  c <-params[(2*length(demogroups)+1):(3*length(demogroups))]#NWD to NWDS
  d <-params[(3*length(demogroups)+1):(4*length(demogroups))] #NWDS to ODS
  e <-params[(4*length(demogroups)+1):(5*length(demogroups))]#ODS to NWDS
  knwd <-params[(5*length(demogroups)+1):(6*length(demogroups))] # NWD mortality
  kod <-params[(6*length(demogroups)+1):(7*length(demogroups))] #OD mortality
  knwds <-params[(7*length(demogroups)+1):(8*length(demogroups))]#NWDS mortality
  kods <-params[(8*length(demogroups)+1):(9*length(demogroups))] #ODS mortality

  # equations ---------------------------------------------------------------
  dNWD <-  - a*NWD + b*OD - knwd*NWD - c*NWD 
  dOD <-  - kod*OD   + a*NWD - b*OD - c*OD
  dNWDS <- c*NWD - knwds*NWDS - d*NWDS + e*ODS
  dODS <-   - kods*ODS - e*ODS + d*NWDS + c*OD
  #mortality trackgin equations
  dODS_M <- kods*ODS 
  dNWDS_M <- knwds*NWDS 
  dOD_M <- kod*OD 
  dNWD_M <- knwd*NWD 

  output = c(dNWD, dOD, dNWDS, dODS, dNWD_M, dOD_M, dNWDS_M, dODS_M) 
  list(output)
}