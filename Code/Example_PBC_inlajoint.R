#Prepare pbc2 data

library(JMbayes2)
library(INLA)
library(INLAjoint)



data(pbc2) # dataset
names(pbc2)
# extract some variable of interest without missing values
Longi <- na.omit(pbc2[, c("id", "years", "status","drug","age", 
                          "sex","year","serBilir","SGOT", "albumin", "edema",
                          "platelets", "alkaline","spiders", "ascites")])

Surv <- Longi[c(which(diff(as.numeric(Longi[,which(colnames(Longi)=="id")]))==1),
                length(Longi[,which(colnames(Longi)=="id")])),-c(7:10, 12:16)]



Surv$death <- ifelse(Surv$status=="dead",1,0) # competing event 1
Surv$trans <- ifelse(Surv$status=="transplanted",1,0) # competing event 2


##Serum bilirubin death event
{#Usual joint model with current value association
M1 <- joint(formSurv = inla.surv(years, death) ~ drug,
            formLong = serBilir ~ (1 + year)*drug +
              (1 + year|id), family = "lognormal",
            dataLong = Longi, dataSurv=Surv, id = "id", timeVar = "year", 
            assoc = "CV", basRisk = "rw2", NbasRisk=25, 
            control=list(int.strategy="eb"))
summary(M1)		

#Joint model with nonlinear current value association
M2 <- joint(formSurv = inla.surv(years, death) ~ drug,
            formLong = serBilir ~ (1 + year)*drug +
              (1 + year|id), family = "lognormal",
            dataLong = Longi, dataSurv=Surv, id = "id", timeVar = "year", 
            assoc = "NL_CV", basRisk = "rw2", NbasRisk=25, 
            control=list(int.strategy="eb"))
summary(M2)		

library(ggplot2)
#Plot the association 
plot(M2, NLeffectonly=TRUE)$NL_Association + 
  geom_hline(yintercept = M1$summary.hyperpar$mean[6], linetype = "dashed", color = "blue") +
  geom_ribbon(aes(ymin = M1$summary.hyperpar$`0.025quant`[6], ymax = M1$summary.hyperpar$`0.975quant`[6]), fill = "blue", alpha = 0.25)
}

##Serum bilirubin transplant event
{#Usual joint model with current value association
  M1 <- joint(formSurv = inla.surv(years, trans) ~ drug,
              formLong = serBilir ~ (1 + year)*drug +
                (1 + year|id), family = "lognormal",
              dataLong = Longi, dataSurv=Surv, id = "id", timeVar = "year", 
              assoc = "CV", basRisk = "rw2", NbasRisk=25, 
              control=list(int.strategy="eb"))
  summary(M1)		
  
  #Joint model with nonlinear current value association
  M2 <- joint(formSurv = inla.surv(years, trans) ~ drug,
              formLong = serBilir ~ (1 + year)*drug +
                (1 + year|id), family = "lognormal",
              dataLong = Longi, dataSurv=Surv, id = "id", timeVar = "year", 
              assoc = "NL_CV", basRisk = "rw2", NbasRisk=25, 
              control=list(int.strategy="eb"))
  summary(M2)		
  
  library(ggplot2)
  #Plot the association 
  plot(M2, NLeffectonly=TRUE)$NL_Association + 
    geom_hline(yintercept = M1$summary.hyperpar$mean[6], linetype = "dashed", color = "blue") +
    geom_ribbon(aes(ymin = M1$summary.hyperpar$`0.025quant`[6], ymax = M1$summary.hyperpar$`0.975quant`[6]), fill = "blue", alpha = 0.25)
}



## SGOT Death event # interesting
{#Usual joint model with current value association
  M1 <- joint(formSurv = inla.surv(years, death) ~ drug,
              formLong = SGOT ~ (1 + year)*drug +
                (1 + year|id), family = "lognormal",
              dataLong = Longi, dataSurv=Surv, id = "id", timeVar = "year", 
              assoc = "CV", basRisk = "rw2", NbasRisk=25, 
              control=list(int.strategy="eb"))
  summary(M1)		
  
  #Joint model with nonlinear current value association
  M2 <- joint(formSurv = inla.surv(years, death) ~ drug,
              formLong = SGOT ~ (1 + year)*drug +
                (1 + year|id), family = "lognormal",
              dataLong = Longi, dataSurv=Surv, id = "id", timeVar = "year", 
              assoc = "NL_CV", basRisk = "rw2", NbasRisk=25, 
              control=list(int.strategy="eb"))
  summary(M2)		
  
  library(ggplot2)
  #Plot the association 
  plot(M2, NLeffectonly=TRUE)$NL_Association + 
    geom_hline(yintercept = M1$summary.hyperpar$mean[6], linetype = "dashed", color = "blue") +
    geom_ribbon(aes(ymin = M1$summary.hyperpar$`0.025quant`[6], ymax = M1$summary.hyperpar$`0.975quant`[6]), fill = "blue", alpha = 0.25)
}

## SGOT transplant event
{#Usual joint model with current value association
  M1 <- joint(formSurv = inla.surv(years, trans) ~ drug,
              formLong = SGOT ~ (1 + year)*drug +
                (1 + year|id), family = "lognormal",
              dataLong = Longi, dataSurv=Surv, id = "id", timeVar = "year", 
              assoc = "CV", basRisk = "rw2", NbasRisk=25, 
              control=list(int.strategy="eb"))
  summary(M1)		
  
  #Joint model with nonlinear current value association
  M2 <- joint(formSurv = inla.surv(years, trans) ~ drug,
              formLong = SGOT ~ (1 + year)*drug +
                (1 + year|id), family = "lognormal",
              dataLong = Longi, dataSurv=Surv, id = "id", timeVar = "year", 
              assoc = "NL_CV", basRisk = "rw2", NbasRisk=25, 
              control=list(int.strategy="eb"))
  summary(M2)		
  
  library(ggplot2)
  #Plot the association 
  plot(M2, NLeffectonly=TRUE)$NL_Association + 
    geom_hline(yintercept = M1$summary.hyperpar$mean[6], linetype = "dashed", color = "blue") +
    geom_ribbon(aes(ymin = M1$summary.hyperpar$`0.025quant`[6], ymax = M1$summary.hyperpar$`0.975quant`[6]), fill = "blue", alpha = 0.25)
}



## albumin Death event # interesting
{#Usual joint model with current value association
  M1 <- joint(formSurv = inla.surv(years, death) ~ drug,
              formLong = albumin ~ (1 + year)*drug +
                (1 + year|id), family = "lognormal",
              dataLong = Longi, dataSurv=Surv, id = "id", timeVar = "year", 
              assoc = "CV", basRisk = "rw2", NbasRisk=25, 
              control=list(int.strategy="eb"))
  summary(M1)		
  
  #Joint model with nonlinear current value association
  M2 <- joint(formSurv = inla.surv(years, death) ~ drug,
              formLong = albumin ~ (1 + year)*drug +
                (1 + year|id), family = "lognormal",
              dataLong = Longi, dataSurv=Surv, id = "id", timeVar = "year", 
              assoc = "NL_CV", basRisk = "rw2", NbasRisk=25, 
              control=list(int.strategy="eb"))
  summary(M2)		
  
  library(ggplot2)
  #Plot the association 
  plot(M2, NLeffectonly=TRUE)$NL_Association + 
    geom_hline(yintercept = M1$summary.hyperpar$mean[6], linetype = "dashed", color = "blue") +
    geom_ribbon(aes(ymin = M1$summary.hyperpar$`0.025quant`[6], ymax = M1$summary.hyperpar$`0.975quant`[6]), fill = "blue", alpha = 0.25)
}

## albumin transplant event
{#Usual joint model with current value association
  M1 <- joint(formSurv = inla.surv(years, trans) ~ drug,
              formLong = albumin ~ (1 + year)*drug +
                (1 + year|id), family = "lognormal",
              dataLong = Longi, dataSurv=Surv, id = "id", timeVar = "year", 
              assoc = "CV", basRisk = "rw2", NbasRisk=25, 
              control=list(int.strategy="eb"))
  summary(M1)		
  
  #Joint model with nonlinear current value association
  M2 <- joint(formSurv = inla.surv(years, trans) ~ drug,
              formLong = albumin ~ (1 + year)*drug +
                (1 + year|id), family = "lognormal",
              dataLong = Longi, dataSurv=Surv, id = "id", timeVar = "year", 
              assoc = "NL_CV", basRisk = "rw2", NbasRisk=25, 
              control=list(int.strategy="eb"))
  summary(M2)		
  
  library(ggplot2)
  #Plot the association 
  plot(M2, NLeffectonly=TRUE)$NL_Association + 
    geom_hline(yintercept = M1$summary.hyperpar$mean[6], linetype = "dashed", color = "blue") +
    geom_ribbon(aes(ymin = M1$summary.hyperpar$`0.025quant`[6], ymax = M1$summary.hyperpar$`0.975quant`[6]), fill = "blue", alpha = 0.25)
}



## edema
{
#Usual joint model with current value association
  M1 <- joint(formSurv = inla.surv(years, death) ~ drug,
              formLong = edema ~ (1 + year)*drug +
                (1 + year|id), family = "lognormal",
              dataLong = Longi, dataSurv=Surv, id = "id", timeVar = "year",
              assoc = "CV", basRisk = "rw2", NbasRisk=25,
              control=list(int.strategy="eb"))
  summary(M1)

  #Joint model with nonlinear current value association
  M2 <- joint(formSurv = inla.surv(years, death) ~ drug,
              formLong = edema ~ (1 + year)*drug +
                (1 + year|id), family = "lognormal",
              dataLong = Longi, dataSurv=Surv, id = "id", timeVar = "year",
              assoc = "NL_CV", basRisk = "rw2", NbasRisk=25,
              control=list(int.strategy="eb"))
  summary(M2)

  library(ggplot2)
  #Plot the association
  plot(M2, NLeffectonly=TRUE)$NL_Association +
    geom_hline(yintercept = M1$summary.hyperpar$mean[6], linetype = "dashed", color = "blue") +
    geom_ribbon(aes(ymin = M1$summary.hyperpar$`0.025quant`[6], ymax = M1$summary.hyperpar$`0.975quant`[6]), fill = "blue", alpha = 0.25)


## edema transplant event
#Usual joint model with current value association
  M1 <- joint(formSurv = inla.surv(years, trans) ~ drug,
              formLong = edema ~ (1 + year)*drug +
                (1 + year|id), family = "dgp",
              dataLong = Longi, dataSurv=Surv, id = "id", timeVar = "year",
              assoc = "CV", basRisk = "rw2", NbasRisk=25,
              control=list(int.strategy="eb"))
  summary(M1)

  #Joint model with nonlinear current value association
  M2 <- joint(formSurv = inla.surv(years, trans) ~ drug,
              formLong = edema ~ (1 + year)*drug +
                (1 + year|id), family = "nmixnb",
              dataLong = Longi, dataSurv=Surv, id = "id", timeVar = "year",
              assoc = "NL_CV", basRisk = "rw2", NbasRisk=25,
              control=list(int.strategy="eb"))
  summary(M2)

  library(ggplot2)
  #Plot the association
  plot(M2, NLeffectonly=TRUE)$NL_Association +
    geom_hline(yintercept = M1$summary.hyperpar$mean[6], linetype = "dashed", color = "blue") +
    geom_ribbon(aes(ymin = M1$summary.hyperpar$`0.025quant`[6], ymax = M1$summary.hyperpar$`0.975quant`[6]), fill = "blue", alpha = 0.25)
}



## platelets Death event # interesting
{#Usual joint model with current value association
  M1 <- joint(formSurv = inla.surv(years, death) ~ drug,
              formLong = platelets ~ (1 + year)*drug +
                (1 + year|id), family = "lognormal",
              dataLong = Longi, dataSurv=Surv, id = "id", timeVar = "year", 
              assoc = "CV", basRisk = "rw2", NbasRisk=25, 
              control=list(int.strategy="eb"))
  summary(M1)		
  
  #Joint model with nonlinear current value association
  M2 <- joint(formSurv = inla.surv(years, death) ~ drug,
              formLong = platelets ~ (1 + year)*drug +
                (1 + year|id), family = "lognormal",
              dataLong = Longi, dataSurv=Surv, id = "id", timeVar = "year", 
              assoc = "NL_CV", basRisk = "rw2", NbasRisk=25, 
              control=list(int.strategy="eb"))
  summary(M2)		
  
  library(ggplot2)
  #Plot the association 
  plot(M2, NLeffectonly=TRUE)$NL_Association + 
    geom_hline(yintercept = M1$summary.hyperpar$mean[6], linetype = "dashed", color = "blue") +
    geom_ribbon(aes(ymin = M1$summary.hyperpar$`0.025quant`[6], ymax = M1$summary.hyperpar$`0.975quant`[6]), fill = "blue", alpha = 0.25)
}

## platelets transplant event
{#Usual joint model with current value association
  M1 <- joint(formSurv = inla.surv(years, trans) ~ drug,
              formLong = platelets ~ (1 + year)*drug +
                (1 + year|id), family = "lognormal",
              dataLong = Longi, dataSurv=Surv, id = "id", timeVar = "year", 
              assoc = "CV", basRisk = "rw2", NbasRisk=25, 
              control=list(int.strategy="eb"))
  summary(M1)		
  
  #Joint model with nonlinear current value association
  M2 <- joint(formSurv = inla.surv(years, trans) ~ drug,
              formLong = platelets ~ (1 + year)*drug +
                (1 + year|id), family = "lognormal",
              dataLong = Longi, dataSurv=Surv, id = "id", timeVar = "year", 
              assoc = "NL_CV", basRisk = "rw2", NbasRisk=25, 
              control=list(int.strategy="eb"))
  summary(M2)		
  
  library(ggplot2)
  #Plot the association 
  plot(M2, NLeffectonly=TRUE)$NL_Association + 
    geom_hline(yintercept = M1$summary.hyperpar$mean[6], linetype = "dashed", color = "blue") +
    geom_ribbon(aes(ymin = M1$summary.hyperpar$`0.025quant`[6], ymax = M1$summary.hyperpar$`0.975quant`[6]), fill = "blue", alpha = 0.25)
}



## alkaline Death event # kind of not really
{#Usual joint model with current value association
  M1 <- joint(formSurv = inla.surv(years, death) ~ drug,
              formLong = alkaline ~ (1 + year)*drug +
                (1 + year|id), family = "lognormal",
              dataLong = Longi, dataSurv=Surv, id = "id", timeVar = "year", 
              assoc = "CV", basRisk = "rw2", NbasRisk=25, 
              control=list(int.strategy="eb"))
  summary(M1)		
  
  #Joint model with nonlinear current value association
  M2 <- joint(formSurv = inla.surv(years, death) ~ drug,
              formLong = alkaline ~ (1 + year)*drug +
                (1 + year|id), family = "lognormal",
              dataLong = Longi, dataSurv=Surv, id = "id", timeVar = "year", 
              assoc = "NL_CV", basRisk = "rw2", NbasRisk=25, 
              control=list(int.strategy="eb"))
  summary(M2)		
  
  library(ggplot2)
  #Plot the association 
  plot(M2, NLeffectonly=TRUE)$NL_Association + 
    geom_hline(yintercept = M1$summary.hyperpar$mean[6], linetype = "dashed", color = "blue") +
    geom_ribbon(aes(ymin = M1$summary.hyperpar$`0.025quant`[6], ymax = M1$summary.hyperpar$`0.975quant`[6]), fill = "blue", alpha = 0.25)
}

## alkaline transplant event
{#Usual joint model with current value association
  M1 <- joint(formSurv = inla.surv(years, trans) ~ drug,
              formLong = alkaline ~ (1 + year)*drug +
                (1 + year|id), family = "lognormal",
              dataLong = Longi, dataSurv=Surv, id = "id", timeVar = "year", 
              assoc = "CV", basRisk = "rw2", NbasRisk=25, 
              control=list(int.strategy="eb"))
  summary(M1)		
  
  #Joint model with nonlinear current value association
  M2 <- joint(formSurv = inla.surv(years, trans) ~ drug,
              formLong = alkaline ~ (1 + year)*drug +
                (1 + year|id), family = "lognormal",
              dataLong = Longi, dataSurv=Surv, id = "id", timeVar = "year", 
              assoc = "NL_CV", basRisk = "rw2", NbasRisk=25, 
              control=list(int.strategy="eb"))
  summary(M2)		
  
  library(ggplot2)
  #Plot the association 
  plot(M2, NLeffectonly=TRUE)$NL_Association + 
    geom_hline(yintercept = M1$summary.hyperpar$mean[6], linetype = "dashed", color = "blue") +
    geom_ribbon(aes(ymin = M1$summary.hyperpar$`0.025quant`[6], ymax = M1$summary.hyperpar$`0.975quant`[6]), fill = "blue", alpha = 0.25)
}



## ascites Death event
{#Usual joint model with current value association
  M1 <- joint(formSurv = inla.surv(years, death) ~ drug,
              formLong = ascites ~ (1 + year)*drug +
                (1 + year|id), family = "binomial",
              dataLong = Longi, dataSurv=Surv, id = "id", timeVar = "year", 
              assoc = "CV", basRisk = "rw2", NbasRisk=25, 
              control=list(int.strategy="eb"))
  summary(M1)		
  
  #Joint model with nonlinear current value association
  M2 <- joint(formSurv = inla.surv(years, death) ~ drug,
              formLong = ascites ~ (1 + year)*drug +
                (1 + year|id), family = "binomial",
              dataLong = Longi, dataSurv=Surv, id = "id", timeVar = "year", 
              assoc = "NL_CV", basRisk = "rw2", NbasRisk=25, 
              control=list(int.strategy="eb"))
  summary(M2)		
  
  library(ggplot2)
  #Plot the association 
  plot(M2, NLeffectonly=TRUE)$NL_Association + 
    geom_hline(yintercept = M1$summary.hyperpar$mean[4], linetype = "dashed", color = "blue") +
    geom_ribbon(aes(ymin = M1$summary.hyperpar$`0.025quant`[4], ymax = M1$summary.hyperpar$`0.975quant`[4]), fill = "blue", alpha = 0.25)
}

## ascites transplant event # interesting
{#Usual joint model with current value association
  M1 <- joint(formSurv = inla.surv(years, trans) ~ drug,
              formLong = ascites ~ (1 + year)*drug +
                (1 + year|id), family = "binomial",
              dataLong = Longi, dataSurv=Surv, id = "id", timeVar = "year", 
              assoc = "CV", basRisk = "rw2", NbasRisk=25, 
              control=list(int.strategy="eb"))
  summary(M1)		
  
  #Joint model with nonlinear current value association
  M2 <- joint(formSurv = inla.surv(years, trans) ~ drug,
              formLong = ascites ~ (1 + year)*drug +
                (1 + year|id), family = "binomial",
              dataLong = Longi, dataSurv=Surv, id = "id", timeVar = "year", 
              assoc = "NL_CV", basRisk = "rw2", NbasRisk=25, 
              control=list(int.strategy="eb"))
  summary(M2)		
  
  library(ggplot2)
  #Plot the association 
  plot(M2, NLeffectonly=TRUE)$NL_Association + 
    geom_hline(yintercept = M1$summary.hyperpar$mean[4], linetype = "dashed", color = "blue") +
    geom_ribbon(aes(ymin = M1$summary.hyperpar$`0.025quant`[4], ymax = M1$summary.hyperpar$`0.975quant`[4]), fill = "blue", alpha = 0.25)
}


##Spiders Death event # interesting
{
M1 <- joint(formSurv = inla.surv(years, death) ~ drug,
            formLong = spiders ~ (1 + year)*drug +
              (1 + year|id), family = "binomial",
            dataLong = Longi, dataSurv=Surv, id = "id", timeVar = "year", 
            assoc = "CV", basRisk = "rw2", NbasRisk=25, 
            control=list(int.strategy="eb"))
summary(M1)		

#Joint model with nonlinear current value association
M2 <- joint(formSurv = inla.surv(years, death) ~ drug,
            formLong = spiders ~ (1 + year)*drug +
              (1 + year|id), family = "binomial",
            dataLong = Longi, dataSurv=Surv, id = "id", timeVar = "year", 
            assoc = "NL_CV", basRisk = "rw2", NbasRisk=25, 
            control=list(int.strategy="eb"))
summary(M2)		

#Spiders plot

plot(M2, NLeffectonly=TRUE)$NL_Association + 
  geom_hline(yintercept = M1$summary.hyperpar$mean[4], linetype = "dashed", color = "blue") +
  geom_ribbon(aes(ymin = M1$summary.hyperpar$`0.025quant`[4], ymax = M1$summary.hyperpar$`0.975quant`[4]), fill = "blue", alpha = 0.25)
}

##Spiders transplant event # interesting
{
  M1 <- joint(formSurv = inla.surv(years, trans) ~ drug,
              formLong = spiders ~ (1 + year)*drug +
                (1 + year|id), family = "binomial",
              dataLong = Longi, dataSurv=Surv, id = "id", timeVar = "year", 
              assoc = "CV", basRisk = "rw2", NbasRisk=25, 
              control=list(int.strategy="eb"))
  summary(M1)		
  
  #Joint model with nonlinear current value association
  M2 <- joint(formSurv = inla.surv(years, trans) ~ drug,
              formLong = spiders ~ (1 + year)*drug +
                (1 + year|id), family = "binomial",
              dataLong = Longi, dataSurv=Surv, id = "id", timeVar = "year", 
              assoc = "NL_CV", basRisk = "rw2", NbasRisk=25, 
              control=list(int.strategy="eb"))
  summary(M2)		
  
  #Spiders plot
  
  plot(M2, NLeffectonly=TRUE)$NL_Association + 
    geom_hline(yintercept = M1$summary.hyperpar$mean[4], linetype = "dashed", color = "blue") +
    geom_ribbon(aes(ymin = M1$summary.hyperpar$`0.025quant`[4], ymax = M1$summary.hyperpar$`0.975quant`[4]), fill = "blue", alpha = 0.25)
}



