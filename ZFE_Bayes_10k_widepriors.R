#################ZFE model###############
########### Author: Ioana Chelcea##########
############################
### Concentrations     = nmol/µl  
### Amounts            = nmol   
### Volumes            = µl       # change to µl?  
### Time               = h          
### Flow rates(k)      = µl/h     #Change to µl/h


##############
#Chorion (diffusion)
#2 compartment ZFE
#Diffusion from chorion into yolk and embryo and between each other
#Perffusion between yolk and embryo
#Diffusion between water and yolk and embryo

###############################

library(deSolve)
library(reshape2)
library(ggplot2)
library(tidyverse)
library(extrafont)
library(deBInfer)
library(coda)

ZFE.model <- function(times, initial_v,parameters) {
  
  with(as.list(c(initial_v, parameters)),{
    if(times<t_hatch){
      
      Pyw  = 10^(K_yw)
      Pew  = pew*K_ow     
      Pye  = pye*K_ow       
      Pcw  = pcw*K_ow  
      Ppw  = 10^(slope*K_ow-intercept) #Model by Nynke Kramer. Parametrize if measured C is missing
      
      ##T-dependent processes 
      kd_y         = KT(T_exp, TR, TA, kd_y_ref) #T-dependent yolk consumption rate (ml/h) for exponential model
      kg_e         = KT(T_exp, TR, TA, kg_e_ref) #T-dependent embryo growth rate (ml /h)
      ksg_e        = KT(T_exp, TR_ksg, TA, ksg_e_ref)
      epiboly_rate = KT(T_exp, TR_epiboly, TA, epiboly_rate_ref) #T-dependent rate of epiboly (as fraction covered of yolk surface area)
      
      #Processes starting at time points
      circulation = organ_start(t_circulation, times)
      metabolism  = organ_start(t_metabolism, times)
      
      ### Volumes  ###
      dVy       = (Vy*-kd_y)      #exponential decay yolk consumption 
      dVe       = kg_e            #linear growth of embryonic body
      V_ZFE     = Vy+Ve           #V of whole ZFE in ml
      dSA_e     = ifelse(epiboly <=1,0,ksg_e)
      d_epiboly = epiboly_rate
      BW        = Ve*0.9979        #density of 0.9979 from van wijk et al 2018
      
      
      Vchorion  = Vchorion_0-V_ZFE
      SA_Yolk   = pi^(1/3)*(6*Vy)^(2/3) #Based on Volume, assuming sphere
      SA_E      = ifelse(epiboly <=1,(pi^(1/3)*(6*V_ZFE)^(2/3))*epiboly,SA_e)#Assume that blastomer engulfs yolk
      Fy        = fr_yolk*(KT(T_exp, TR_Fcard, TA, F_card_ref) * (BW/BW_Fcard_ref)^(-0.1))#cardiac output adjusted for temperature and BW of the simulation
      
      ### Clearance  ###
      
      Cl_sc     = Cl_liv*Ve         #Scaled liver clearance
      Vmax      = Cl_sc*Km
      dA_cleared_acum =(Vmax*(Ae/Ve))/(Km+(Ae/Ve))*metabolism
      
      ### Compartments  ### 
      
      dA_Chorion = (Kp_Chorion*SA_Chorion*(Aw/Vw-(Achorion/Vchorion)/Pcw)
                    +Kp_E*SA_E*((Ae/Ve)/Pew-Achorion/Vchorion)
                    +Kp_Yolk*SA_Yolk*ifelse(epiboly <=1,1-epiboly,0)*((Ay/Vy)/Pyw-Achorion/Vchorion))
      
      dA_Yolk    = (Kp_Yolk*SA_Yolk*ifelse(epiboly <=1, epiboly,1)*(Ae/Ve-(Ay/Vy)/Pye)
                    +Kp_Yolk*(SA_Yolk*ifelse(epiboly <=1,1-epiboly,0))*(Achorion/Vchorion-(Ay/Vy)/Pyw)
                    +Fy*(Ae/Ve-(Ay/Vy)/Pye))
      
      dA_Embryo  = (Kp_E*SA_E*(Achorion/Vchorion-(Ae/Ve)/Pew)
                    +Kp_Yolk*SA_Yolk*ifelse(epiboly <=1, epiboly,1)*((Ay/Vy)/Pye-Ae/Ve)
                    +Fy*((Ay/Vy)/Pye-Ae/Ve)
                    -(Vmax*(Ae/Ve))/(Km+(Ae/Ve))*metabolism)
      
      
      dA_Water  = (Kp_Chorion*SA_Chorion*((Achorion/Vchorion)/Pcw-Aw/Vw)
                   +SA_plastic*(Ap/Ppw-(Aw)))
      
      dA_Plastic = SA_plastic*((Aw)-Ap/Ppw)
      
      ### Volumes, masses and concentrations  ###    
      
      A_ZFE    = Ae+Ay
      C_ZFE    = A_ZFE/(V_ZFE)*1E3           #C in whole ZFE (nmol/ml)
      C_embryo = Ae/Ve*1E3                   #C in embryo body (nmol/ml)
      C_yolk   = Ay/Vy*1E3                   #C in yolk (nmol/ml) 
      C_water  = Aw/Vw*1E3                   #C in water (nmol/ml)
      A_ZFE_tot  = A_ZFE+Achorion
      Mass_balance   = (Aw_0)-(Aw+A_ZFE+Achorion+Ap+Aclacum)#Mass balance: Amount of compound at start of experiment-(sum of amounts in all modeled compartments)
      
    } 
    else {  
      
      ###################################################
      ############### MODEL AFTER HATCHING ##########
      ###################################################
      
      Pyw  = 10^(K_yw)
      Pew  = pew*K_ow     
      Pye  = pye*K_ow       
      Pcw  = pcw*K_ow
      Ppw  = 10^(slope*K_ow-intercept) #Model by Nynke Kramer. Parametrize if measured C is missing
      
      #Processes starting at time points
      circulation = organ_start(t_circulation, times)
      metabolism  = organ_start(t_metabolism, times)
      
      ##T-dependent processes 
      kd_y   = KT(T_exp, TR, TA, kd_y_ref) #T-dependent yolk consumption rate (ml/h) for exponential model
      kg_e   = KT(T_exp, TR, TA, kg_e_ref) #T-dependent embryo growth rate (ml /h)
      ksg_e  = KT(T_exp, TR_ksg, TA, ksg_e_ref)
      
      ### Volumes  ###
      
      dVy       = (Vy*-kd_y)      #exponential decay yolk consumption 
      dVe       = kg_e            #linear growth of embryonic body
      dSA_e     = ksg_e 
      d_epiboly = 0
      SA_Yolk   = pi^(1/3)*(6*Vy)^(2/3)
      
      BW  = Ve*0.9979
      Fy  = fr_yolk*(KT(T_exp, TR_Fcard, TA, F_card_ref) * ((BW/BW_Fcard_ref)^(-0.1)))#*BW #cardiac output adjusted for temperature and BW of the simulation (ml/d)
      
      ### Clearance  ###
      
      Cl_sc     = Cl_liv*Ve         #Scaled liver clearance
      Vmax      = Cl_sc*Km
      dA_cleared_acum =(Vmax*(Ae/Ve))/(Km+(Ae/Ve))*metabolism
      
      ### Compartments  ### 
      dA_Chorion= 0
      
      dA_Yolk   = (Fy*((Ae/Ve)-(Ay/Vy)/Pye) #Amount in yolk (flows in and out from water and embryo)
                   +Kp_Yolk*SA_Yolk*(Ae/Ve-(Ay/Vy)/Pye))
      
      dA_Embryo = (SA_e*Kp_E*((Aw/Vw)-((Ae/Ve)/Pew))
                   +Kp_Yolk*SA_Yolk*((Ay/Vy)/Pye-Ae/Ve)
                   +Fy*((Ay/Vy)/Pye - (Ae/Ve))   #Amount in Embryo (flows in and out from yolk and water)
                   -(Vmax*(Ae/Ve))/(Km+(Ae/Ve))*metabolism)  #Clearance only with a fraction of the amount that enters the embryo (Assumption: everything is evenely distributed in the embryo at early life stage--> same fraction as the volume fraction)
      
      dA_Water  = (SA_e*Kp_E*((Ae/Ve)/Pew-(Aw/Vw)) #Amount in water (flows in and out from yolk and embryo)                                 
                   +SA_plastic*(Ap/Ppw-(Aw)))#Amount in water 
      
      
      dA_Plastic = SA_plastic*((Aw)-Ap/Ppw)
      
      ### Volumes, masses and concentrations  ###    
      
      V_ZFE    = Ve+Vy                       #V of whole ZFE in ml
      A_ZFE    = Ay+Ae                       #Amount in whole ZFE in µmol
      C_ZFE    = A_ZFE/(V_ZFE)*1E3           #C in whole ZFE (nmol/ml)
      C_embryo = Ae/Ve*1E3                   #C in embryo body (nmol/ml)
      C_yolk   = Ay/Vy*1E3                   #C in yolk (nmol/ml)
      C_water  = Aw/Vw*1E3                   #C in water (nmol/ml)
      
      Mass_balance   = (Aw_0)-(Aw+A_ZFE+Achorion+Ap+Aclacum)#Mass balance: Amount of compound at start of experiment-(sum of amounts in all modeled compartments)
      
      Vchorion  = 0
      epiboly   = 0
      SA_E = 0
      A_ZFE_tot = A_ZFE
    }
    
    list(c(dVy,dVe,dSA_e,d_epiboly,dA_Chorion, dA_Yolk, dA_Embryo, dA_cleared_acum, dA_Water, dA_Plastic), "Mass_balance" = Mass_balance,"Aw_0"=Aw_0, "Azfe"= A_ZFE,"A_ZFE_tot" = A_ZFE_tot,
         "C_ZFE"=C_ZFE,"Ce"= C_embryo,"Cy"= C_yolk,"C_water"=C_water, "V_ZFE"= V_ZFE, "V_Chorion"= Vchorion, "SA_e_hatch"= SA_E, "SA_yolk"= SA_Yolk, "Fy"= Fy, "kd_y"= kd_y)
    
  })}
##################################################
############### 2. Functions ####################
##################################################


#Function for switching on organs at specific timepoint
organ_start<- function(t_0,times){
  organ = ifelse(t_0<=times,1,0)
  return(organ)
}

#Temperature adjustment function(Arrhenius temperature function)
KT <- function(T, TR, TA, KR){ KR*exp((TA/TR) - (TA/T)) }


############# LIKELIHOOD MODEL ################

ZFE_obs_model <- function(data, sim.data, samp) {
  
  #Following Johnson et al. (2013) we employ a small correction ec 
  #that is needed because the DE solution can equal zero, whereas 
  #the parameter epsilon of the Poisson likelihood must be strictly positive.
  w.obs <- which(sim.data[,"time"] %in% data$time)
  ## The simulated data set at only the correct tims
  simdat <- sim.data[w.obs,]
  
  epsilon <-1e-6
  ## log.Predition 
  #log.yhat      <-log(simdat[ ,"A_ZFE_tot"]+epsilon)
  
  ## log.Observed data
 # log.y         <-log(data$A_ZFE)
  
  # The method of Maximum likelihood
 # sig2            <- 0.05 ## Model error (residuals); mostly between 0.3 and 0.5 (corresponding to coefficients of variation of about 30-50%); Bois et al. (1998, 2000)
  #llik_ZFE  <- -2*sum ((dnorm (data$A_ZFE,
 #                                    mean =log.yhat,
 #                                    sd = sqrt(sig2),
 #                                    log = TRUE)))
  #llik_ZFE <- -2*sum(dnorm(data$A_ZFE, mean = log(simdat[ ,"A_ZFE_tot"]+ epsilon), sd = samp[["sdlog.ZFE"]], log = TRUE))
  llik_ZFE <- sum(dlnorm(data$A_ZFE, meanlog = log(simdat[ ,"A_ZFE_tot"]+ epsilon), sdlog = samp[["sdlog.ZFE"]], log = TRUE))
  
  llik <-llik_ZFE
  #print(c(simdat[ ,"A_ZFE_tot"]))
  return(llik)
}




##################################################
############### 3. PARAMETERS ####################
##################################################



### Experiment speciffic parameters ####


Vw_well    = 2000 #Volume of water in the whole well (µl)
embryo_nr  = 10
Cw_0       = 5.8/1000 #nmol/µl (21 µM measured in water during experiment)
well_r     = 15.7/2   #mm radius of well
SA_plastic_exp = 0#(2*(pi*well_r^3+Vw_well)/(well_r)-pi*well_r^2)/embryo_nr #mm^2 how much plastic is the water in contact with (cylinder without top circle)

Vw         = debinfer_par(name = "Vw", var.type = "de", fixed = TRUE,
                          value = Vw_well/embryo_nr) #Volume of water in µl 
Aw_0       = debinfer_par(name = "Aw_0", var.type = "de", fixed = TRUE, value = Cw_0*Vw_well/embryo_nr)

SA_plastic = debinfer_par(name = "SA_plastic", var.type = "de", fixed = TRUE,
                          value = SA_plastic_exp)
T_exp      = debinfer_par( name = "T_exp", var.type = "de", fixed = TRUE,
                           value = 299.15, prior = "unif", hypers = list(min = 298.15, max = 300.15),
                           prop.var = c(0.0001),samp.type = "rw-ref")   # Experimental T in Kelvin

##################################################################
############### Physiological parameters #####################
#############################################################

####Times for various events####
t_hatch       = debinfer_par( name = "t_hatch", var.type = "de", fixed = TRUE, 
                              value = 60 ,prior = "norm", hypers = list(mean = 60, sd = 6),
                              prop.var = c(0.00001), samp.type = "rw") #Hacthing time
t_circulation = debinfer_par( name = "t_circulation", var.type = "de", fixed = TRUE, 
                              value = 36, prior = "norm", hypers = list(mean = 36, sd = 6 ),
                              prop.var = c(0.00001), samp.type = "rw") #Time at which perfusion starts
t_metabolism  = debinfer_par( name = "t_metabolism", var.type = "de", fixed = FALSE,
                              value = 71.76, prior = "unif", hypers = list(min = 6, max = 120),
                              prop.var = c(0.0001), samp.type = "rw-ref") #Time at which UGTs start being active


### Reference parameters for T dependency ####
### Parameter list ###
kd_y_ref <- debinfer_par(name = "kd_y_ref", var.type = "de", fixed = TRUE,
                         value = 0.0144, prior = "norm", hypers = list( mean = 0.0144, sd = 0.001),
                         prop.var = c(0.0000001), samp.type = "rw") #yolk consumption rate (µl/h) for exponential model. Fitted on data by Simeon et al 2020 and Halbach et al 2020

kg_e_ref  = debinfer_par(name = "kg_e_ref", var.type = "de", fixed = TRUE,
                         value = 0.00243,prior = "norm", hypers = list( mean = 0.00243, sd = 0.0017),
                         prop.var = c(0.0000001), samp.type = "rw")#embryo growth rate (µl /h) #based on mean slope between linear fit of data from Simeon et al 2020 (slope = 3E-6) and data from Halbach et al 2020 (slope = 2E-6)

TR      <-  debinfer_par(name = "TR", var.type = "de", fixed = TRUE,value =as.numeric(299.15)) 
# Reference T in Kelvin for growth and yolk consumption (Simeon et al 2020 and Halbach et al 2020)

TR_ksg    = debinfer_par(name = "TR_ksg", var.type = "de", fixed = TRUE,
                         value = as.numeric(301.65))#, prior = "unif", hypers = list(min = 26 + 273.15, max = 31 + 273.15),
#prop.var = c(1), samp.type = "rw")  #Ref temperature for SA growth based on Hagedorn et al 1998 (T = 28.5), Kimmel et al 1995 (T= 28.5) and  Guo et al.2017 (T not given)

ksg_e_ref = debinfer_par(name = "ksg_e_ref", var.type = "de", fixed = TRUE,
                         value = 0.0183, prior = "unif", hypers = list(min = 0.017, max =0.024),
                         prop.var= c(0.000001), samp.type = "rw")  #Surface area growth rate at 10 hpf to 3-5 dpf, linearlly fit based on data by Hagedorn et al 1998 (T = 28.5), Kimmel et al 1995 (T= 28.5) and  Guo et al.2017 (T not given)

BW_Fcard_ref = debinfer_par(name = "BW_Fcard_ref", var.type = "de", fixed = TRUE,
                            value = 0.097, prior = "norm", hypers = list(mean = 0.097, sd=	0.02),
                            prop.var = c(0.000001), samp.type = "rw")  #BW of an embryo at 48 hpf at T =  28.25 without Chorion or yolk

TR_Fcard     = debinfer_par(name = "TR_Fcard", var.type = "de", fixed = TRUE,
                            value = as.numeric(301.4))#, prior = "unif", hypers = list(min= 300.15, max = 302.15),
#prop.var = c(0.5), samp.type = "rw") #Reference T for cardiac output (Kelvin), Average of Jacob et al 2002 and Malone et al 2007

F_card_ref   = debinfer_par(name = "F_card_ref", var.type = "de", fixed = TRUE,
                            value = 2.024, prior = "norm", hypers = list(mean = 2.024, sd = 0.5),
                            prop.var= c(0.000001), samp.type = "rw") #µl/h Average of 48 hpf from Jacob et al 2002 and Malone et al 2007

fr_yolk      = debinfer_par(name = "fr_yolk", var.type = "de", fixed = FALSE,
                            value = 0.88, prior = "unif", hypers = list(min = 0.01, max = 1 ),
                            prop.var = c(0.0001), samp.type = "rw-ref")#Cardiac output fraction going to yolk

TR_epiboly   = debinfer_par(name = "TR_epiboly", var.type = "de", fixed = TRUE,
                            value = 301.65)#, prior = "unif", hypers = list(min=28 + 273.15, max= 29 + 273.15),
#prop.var = c(0.01), samp.type = "rw")  # Kimmel et al 1995

epiboly_rate_ref = debinfer_par(name = "epiboly_rate_ref", var.type = "de", fixed = TRUE,
                                value = 0.1255, prior = "unif", hypers = list(min = 0.1, max = 0.15),
                                prop.var = c(0.000001), samp.type = "rw-ref") # (in fraction of yolk are covered) fitted from Kimmel et al 1995

TA       <- debinfer_par(name = "TA", var.type = "de", fixed = TRUE,
                         value = 6930, prior = "unif", hypers = list(min = 5544, max = 8316),
                         prop.var = c(0.000001), samp.type = "rw") # Arrhenius T in Kelvin (Billat et al.)



######## Surface areas ##########

SA_Chorion = debinfer_par(name = "SA_Chorion", var.type = "de", fixed = TRUE,
                          value = 4.85, prior = "norm", hypers = list(mean = 4.85, sd = 0.8),
                          prop.var =c(0.000001), samp.type = "rw")    #Surface area in mm^2 based on average from Kimmel et al 1995 and Brox et al 2014


Vchorion_0 = debinfer_par(name = "Vchorion_0", var.type = "de", fixed = TRUE,
                          value = 0.98, prior = "norm", hypers = list(mean = 0.98, sd= 0.18),
                          prop.var = c(0.000001), samp.type = "rw")    #(µl) Volume of periviteline space (excluding embryo and yolk) at t0 from images by Kimmel et al assuming perfect sphere

######## Permeability constants ##########

Kp_Chorion  = debinfer_par(name = "Kp_Chorion", var.type = "de", fixed = FALSE,
                           value = 0.8, prior = "unif", hypers = list(min = 0.001,max = 10),
                           prop.var = c(0.0001), samp.type ="rw-ref")

Kp_Yolk     = debinfer_par(name = "Kp_Yolk", var.type = "de", fixed = FALSE,
                           value = 0.51, prior = "unif", hypers = list(min = 0.001, max = 10),
                           prop.var = c(0.0001), samp.type = "rw-ref")
Kp_E        = debinfer_par(name = "Kp_E", var.type = "de", fixed = FALSE,
                           value = 0.55, prior = "unif", hypers = list(min = 0.001, max = 10),
                           prop.var = c(0.0001),samp.type = "rw-ref")


#### Clearance ###

Cl_liv = debinfer_par(name = "Cl_liv", var.type = "de", fixed = TRUE,
                      value = 0.725, prior = "unif", hypers = list(min = 0, max = 1000),
                      prop.var= c(0.000001), samp.type = "rw")#µl/h/µl fish, adult value transformed from Chelcea et al 2022.

Km = debinfer_par(name = "Km", var.type = "de", fixed = FALSE,
                  value = 0.0039, prior = "unif", hypers = list(min = 0.0001, max = 0.01),
                  prop.var = c(0.0001), samp.type = "rw-ref") 

##################################################################
############### Chemical parameters #####################
#############################################################

K_ow  = debinfer_par(name = "K_ow", var.type = "de", fixed = TRUE,
                     value = 4.34, prior = "unif", hypers = list(min= 4.08, max = 5),
                     prop.var = c(0.000001), samp.type = "rw-ref") #log octanol-water partitioning

K_yw  = debinfer_par(name = "K_yw", var.type = "de", fixed = TRUE,
                     value = 4.01, prior = "norm", hypers = list(min= 1E-3, max = 10),
                     prop.var = c(0.000001), samp.type = "rw") #log octanol-water partitioning

### Partition coefficients (unitless) ###

#BPA Yolk-water partitioning measured by Ulrich et al. 2020 used for yolk-chorion Partitioning
pew  = debinfer_par(name = "pew", var.type = "de", fixed = FALSE,
                    value = 5.68, prior = "unif", hypers = list(min = 0.1, max = 100),
                    prop.var = c(0.00001), samp.type = "rw")    #embryo water partitionining
pye  = debinfer_par(name = "pye", var.type = "de", fixed = FALSE,
                    value = 1.34, prior = "unif", hypers = list(min = 0.1, max = 100),
                    prop.var =  c(0.00001), samp.type = "rw") #yolk embryo partitioning Base it on how much more fat content is in yolk compared to embryo?
pcw  = debinfer_par(name = "pcw", var.type = "de", fixed = FALSE,
                    value = 0.314, prior = "unif", hypers = list(min = 0.1, max = 100),
                    prop.var = c(0.00001), samp.type = "rw") #Chorion water (close to 0?)
slope = debinfer_par(name = "slope", var.type = "de", fixed = TRUE,
                     value = 0.97 , prior = "unif", hypers = list(min = 0.1, max = 100),
                     prop.var =c(0.00001), samp.type = "rw")    # for calculating plastic-water partition coefficient (Kramer 2010)
intercept = debinfer_par(name = "intercept", var.type = "de", fixed = TRUE,
                         value = 6.94, prior = "unif", hypers = list(min = 6.79, max = 100),
                         prop.var =c(0.00001), samp.type = "rw")



sdlog.ZFE <-debinfer_par(name = "sdlog.ZFE", var.type ="obs", fixed = FALSE, value = 1.5, 
                         prior = "norm", hypers = list(mean= 1.5, sd = 1), 
                         prop.var= c(0.0001), samp.type = "rw")



############INITs################

Ve       = debinfer_par(name = "Ve", var.type = "init", fixed = TRUE,
                        value = 6.73E-3, prior = "unif", hypers= list(min = 4.6E-3, max = 1.06E-2),
                        prop.var = c(0.000001), samp.type = "rw") #Volume of embryo at t0 (µl) #based on hemisphere from images (Kimmel et al 1995) at 0.17 hpf

Vy       = debinfer_par(name = "Vy", var.type = "init", fixed = TRUE,
                        value = 2.7E-01, prior = "unif", hypers = list(min= 0.16, max = 0.27),
                        prop.var = c(0.000001), samp.type = "rw")  #Volume yolk at t0-t5 in µl (Simeon et al 2020, Kimmel et al 1995,Brox et al 2014) assuming perfect sphere

SA_e    = debinfer_par(name = "SA_e", var.type = "init", fixed = TRUE,
                       value = 1.71, prior = "unif", hypers = list(min = 1.37, max = 2.052),
                       prop.var = c(0.000001), samp.type = "rw")   #Surface area of spherical blastoderm at 10 hpf (Average from Kimmel et al and Hagedorn et al)
epiboly  = debinfer_par(name = "epiboly", var.type = "init", fixed = TRUE, value = 0.00001)
Achorion = debinfer_par(name = "Achorion", var.type = "init", fixed = TRUE, value = 0)
Ay       = debinfer_par(name = "Ay", var.type = "init", fixed = TRUE,value = 0)
Ae       = debinfer_par(name = "Ae", var.type = "init", fixed = TRUE, value = 0)
Aclacum  = debinfer_par(name = "Aclacum", var.type = "init", fixed = TRUE,value =0)
Aw       = debinfer_par(name = "Aw", var.type = "init", fixed = TRUE, value = Cw_0*Vw_well/embryo_nr)
Ap       = debinfer_par(name = "Ap", var.type = "init", fixed = TRUE, value = 0) #Amount of compound in water in µmol/ml  (10 µM)



### Parameter list ###
mcmc.pars <- setup_debinfer(# Experiment speciffic parameters ####
                            
                            Vw,
                            Aw_0,
                            SA_plastic,
                            T_exp,
                            ####Times for various events####
                            t_hatch,
                            t_circulation,
                            t_metabolism,
                            ### Reference parameters for T dependency ####
                            kd_y_ref,
                            kg_e_ref, 
                            TR,
                            TR_ksg,
                            ksg_e_ref, 
                            BW_Fcard_ref,
                            TR_Fcard,
                            F_card_ref,
                            fr_yolk,
                            TR_epiboly,
                            epiboly_rate_ref,
                            TA,
                            ######## Volumes ##########
                            Vchorion_0,
                            ######## Surface areas ##########
                            SA_Chorion,
                            ######## Permeability constants ##########
                            Kp_Chorion,
                            Kp_Yolk,
                            Kp_E,
                            #### Clearance ###
                            Cl_liv,
                            Km,
                            ### Partition coefficients ###
                            K_ow,
                            K_yw,
                            pew,
                            pye,
                            pcw,
                            slope,
                            intercept,
                            ######Observations
                            sdlog.ZFE,
                            #### Initial conditions##
                            Vy,
                            Ve,
                            SA_e,
                            epiboly,
                            Achorion,
                            Ay,
                            Ae,
                            Aclacum,
                            Aw,
                            Ap)



##################################################
############### 3. RUN MODELS ####################
##################################################


#######EXPERIMENTAL Data#####
Exp_time = c(12, 24, 72, 96, 120)#hpf
Exp_data = c(0.06,
             0.06,
             0.11,
             0.10,
             0.07) #mean BPZ (nmol/embryo)
Exp_water = c((5.8/1000)*2000/10, #21µM /1000= nmol/µl *2000µl = 42.8 nmol/well
              (5.5/1000)*2000/10,
              (4.9/1000)*2000/10,
              (4.5/1000)*2000/10,
              (3.5/1000)*2000/10)
Exp_ZFE  = data.frame(Exp_time, Exp_data, Exp_water)
colnames(Exp_ZFE)=c("time", "A_ZFE", "A_Water")


iter = 10000

ZFE_mcmc <- de_mcmc(N = iter, data = Exp_ZFE, de.model = ZFE.model,
                    obs.model = ZFE_obs_model, all.params = mcmc.pars,
                    Tmax = max(Exp_ZFE$time), data.times = Exp_ZFE$time, cnt = 1000,sizestep = 1,
                    plot = FALSE, solver = "ode", verbose.mcmc = TRUE, method = "lsoda")




#plot(Exp_ZFE_plot, xlab = "Time (h)", ylab = "Amount in ZFE", xlim = c(0,168), ylim=c(0,0.2))
#Alternatively we can plot an ensemble of posterior trajectories using
#plot(post_traj, plot.type = "ensemble", col = "#FF000040")

## MCMC convergence assessment
## uncomment the following section of R code to run an additional two chains and calculate the potential scale reduction factor to assess convergence
## deBInfer currently does not provide automated facilities for running multiple chains, we need to  repeat the inference call to run additional chains



#####################PLOT


dev.off()
par(mar = c(1, 1, 1, 1))
post_prior_densplot(ZFE_mcmc, burnin = 10)
post_traj <- post_sim(ZFE_mcmc, n = 900, times = seq(0,168,by = 0.1),
                      output = "all", prob = 0.95)

post_prior_densplot(ZFE_mcmc, param = "Kp_Chorion", show.obs = FALSE, main = "Kp Chorion")
post_prior_densplot(ZFE_mcmc, param = "Kp_Yolk", show.obs = FALSE, main = "Kp Yolk")
post_prior_densplot(ZFE_mcmc, param = "Kp_E", show.obs = FALSE, main = "Kp Embryo")
post_prior_densplot(ZFE_mcmc, param = "Km", show.obs = FALSE, main = "Km")
post_prior_densplot(ZFE_mcmc, param = "fr_yolk", show.obs = FALSE, main = "Fraction of Fcard going to yolk")
post_prior_densplot(ZFE_mcmc, param = "pew", show.obs = FALSE, main = "Partition embryo-water (not log Kow adjusted)")
post_prior_densplot(ZFE_mcmc, param = "pye", show.obs = FALSE, main = "Partition yolk-embryo (not log Kow adjusted)")
post_prior_densplot(ZFE_mcmc, param = "pcw", show.obs = FALSE, main = "Partition chorion-water (not log Kow adjusted)")
post_prior_densplot(ZFE_mcmc, param = "t_metabolism", show.obs = FALSE, main = "Metabolism start")

dev.off()
par(mfrow = c(1,1))
plot(post_traj, plot.type = "medianHDI", auto.layout = FALSE,lty = c(1,2),
     col = c("red","grey"))

#post_prior_densplot(ZFE_mcmc, param = "Kp_Chorion", show.obs = FALSE, main = "t_metabolism")
#post_prior_densplot(ZFE_mcmc, param = "t_hatch", show.obs = FALSE, main = "time hatch")
summary(ZFE_mcmc)
Exp_ZFE_plot  = data.frame(Exp_time, Exp_data)
#show the trajectories of all three state variables in a single plot,
plot(Exp_ZFE_plot, xlab = "Time (h)", ylab = "Amount in ZFE", xlim = c(0,168), ylim=c(0,0.2))
for(i in seq_along(post_traj$sims)) {
  DATA1 <- as.data.frame(post_traj$sims[i])
  lines(DATA1[,7] ~ DATA1[,1], col = "yellow")#first one is number of var
  lines(DATA1[,8] ~ DATA1[,1],col = "gray")
  lines(DATA1[,14] ~ DATA1[,1],col = "blue")
}

