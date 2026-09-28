#Prepare pbc2 data

library(JMbayes2)
library(INLA)
library(INLAjoint)
library(ggplot2)

# Setup

{

data(pbc2) # dataset
names(pbc2)
# extract some variable of interest without missing values
pbc2$edema <- ifelse(pbc2$edema == "No edema", "No", "Yes")
Longi <- na.omit(pbc2[, c("id", "years", "status","drug","age", 
                          "sex","year","serBilir","SGOT", "albumin", "edema",
                          "platelets", "alkaline","spiders", "ascites")])

Surv <- Longi[c(which(diff(as.numeric(Longi[,which(colnames(Longi)=="id")]))==1),
                length(Longi[,which(colnames(Longi)=="id")])),-c(7:10, 12:16)]



Surv$death <- ifelse(Surv$status=="dead",1,0) # competing event 1
Surv$trans <- ifelse(Surv$status=="transplanted",1,0) # competing event 2


}

##Serum bilirubin death event
{#Usual joint model with current value association
M1 <- joint(formSurv = inla.surv(years, death) ~ drug,
            formLong = serBilir ~ (1 + year)*drug +
              (1 + year|id), family = "lognormal",
            dataLong = Longi, dataSurv=Surv, id = "id", timeVar = "year", 
            assoc = "CV", basRisk = "rw2", NbasRisk=25, 
            control=list(int.strategy="eb"))

#Joint model with nonlinear current value association
M2 <- joint(formSurv = inla.surv(years, death) ~ drug,
            formLong = serBilir ~ (1 + year)*drug +
              (1 + year|id), family = "lognormal",
            dataLong = Longi, dataSurv=Surv, id = "id", timeVar = "year", 
            assoc = "NL_CV", basRisk = "rw2", NbasRisk=25, 
            control=list(int.strategy="eb"))


#Plot the non-linear association 
plot(M2$summary.random$uv1$mean, M2$summary.random$NL_CV_L1_S1$mean, main = "Serum bilirubin association: \n death event",
     xlab = "Linear predictor value", ylab = "Effect", type = "p", pch = 20)

# plot non-linear confidence intervals
pol <- matrix(0, nrow = length(M2$summary.random$uv1$mean), ncol = 3)
pol[,1] <- M2$summary.random$uv1$mean
pol[,2] <- M2$summary.random$NL_CV_L1_S1$`0.975quant`
pol[,3] <- M2$summary.random$NL_CV_L1_S1$`0.025quant`
colnames(pol) <- c("x", "top", "bottom")

pol <- pol[order(pol[,1]),]

polygon(x = c(pol[,1], rev(pol[,1])),
        y = c(pol[,2], rev(pol[,3])),
        col = adjustcolor("gray", 0.3), border = adjustcolor("black", 0.3))

# Reference line
abline(h = 0)

# plot linear association
abline(h = M1$summary.hyperpar$mean[6], col = "blue", lty = 4, lwd = 2) 

# linear confidence intervals
polygon(x = c(min(M2$summary.random$uv1$mean)-10, max(M2$summary.random$uv1$mean)+10,
              max(M2$summary.random$uv1$mean)+10, min(M2$summary.random$uv1$mean)-10),
        y = c(M1$summary.hyperpar$`0.975quant`[6], M1$summary.hyperpar$`0.975quant`[6],
              M1$summary.hyperpar$`0.025quant`[6], M1$summary.hyperpar$`0.025quant`[6]),
        col = adjustcolor("blue", 0.25), border = adjustcolor("blue", 0.25))





# plot survival curves

n1 <- length(Longi$id)
n2 <- length(Surv$id)
maxMeas <- 16

predl <- M2$summary.fitted.values$mean[1:n1] # Longitudinal preds for all measurements
preds <- M2$summary.fitted.values$mean[(1+n1):(n2+n1)] # Survival preds for all patient





# Plot IWRES against time (????) (time??)
plot(Longi$years, M2$residuals$deviance.residuals[1:n1])


# PLot IWRES against biomarker prediction (????)
plot(predl, M2$residuals$deviance.residuals[1:n1])

# Get Cox-Snell residuals

M2$basRisk



 # Plot with ggplot (nl effect only true)
  plot(M2, NLeffectonly=T)$NL_Association +
    geom_hline(yintercept = M1$summary.hyperpar$mean[6], linetype = "dashed", color = "blue") +
    geom_ribbon(aes(ymin = M1$summary.hyperpar$`0.025quant`[6], ymax = M1$summary.hyperpar$`0.975quant`[6]), fill = "blue", alpha = 0.25)
  
  # Plot with ggplot (nl effect only false)
  plot(M2, NLeffectonly=F)$NL_Association +
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
  

  #Plot the non-linear association 
  plot(M2$summary.random$uv1$mean, M2$summary.random$NL_CV_L1_S1$mean, main = "Serum bilirubin association: \n transplant event",
       xlab = "Linear predictor value", ylab = "Multiplication factor", type = "p", pch = 20)
  
  # plot the confidence intervals
  
  pol <- matrix(0, nrow = length(M2$summary.random$uv1$mean), ncol = 3)
  pol[,1] <- M2$summary.random$uv1$mean
  pol[,2] <- M2$summary.random$uv1$`0.975quant`
  pol[,3] <- M2$summary.random$uv1$`0.025quant`
  colnames(pol) <- c("x", "top", "bottom")
  
  pol <- pol[order(pol[,1]),]
  
  polygon(x = c(pol[,1], rev(pol[,1])),
          y = c(pol[,2], rev(pol[,3])),
          col = "gray", density = 20)
  
  # reference line
  abline(h = 0)
  
  # plot linear association
  abline(h = M1$summary.hyperpar$mean[6], col = "blue", lty = 4, lwd = 2) 
  
  # linear confidence interval
  polygon(x = c(min(M2$summary.random$uv1$mean)-10, max(M2$summary.random$uv1$mean)+10,
                max(M2$summary.random$uv1$mean)+10, min(M2$summary.random$uv1$mean)-10),
          y = c(M1$summary.hyperpar$`0.975quant`[6], M1$summary.hyperpar$`0.975quant`[6],
                M1$summary.hyperpar$`0.025quant`[6], M1$summary.hyperpar$`0.025quant`[6]),
          col = "blue", density = 20)
}



## SGOT Death event # kind of interesting #1
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
  

  #Plot the association 
  plot(M2$summary.random$uv1$mean, M2$summary.random$NL_CV_L1_S1$mean, main = "SGOT association: \n death event",
       xlab = "Linear predictor value", ylab = "multiplication factor", type = "p", pch = 20)
  
  abline(h = 0)
  abline(h = M1$summary.hyperpar$mean[6], col = "blue", lty = 4, lwd = 2) 
  
  polygon(x = c(min(M2$summary.random$uv1$mean)-10, max(M2$summary.random$uv1$mean)+10,
                max(M2$summary.random$uv1$mean)+10, min(M2$summary.random$uv1$mean)-10),
          y = c(M1$summary.hyperpar$`0.975quant`[6], M1$summary.hyperpar$`0.975quant`[6],
                M1$summary.hyperpar$`0.025quant`[6], M1$summary.hyperpar$`0.025quant`[6]),
          col = "blue", density = 20)
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
  

  #Plot the association 
  plot(M2$summary.random$uv1$mean, M2$summary.random$NL_CV_L1_S1$mean, main = "SGOT association: \n transplant event",
       xlab = "Linear predictor value", ylab = "multiplication factor", type = "p", pch = 20)
  
  abline(h = 0)
  abline(h = M1$summary.hyperpar$mean[6], col = "blue", lty = 4, lwd = 2) 
  
  polygon(x = c(min(M2$summary.random$uv1$mean)-10, max(M2$summary.random$uv1$mean)+10,
                max(M2$summary.random$uv1$mean)+10, min(M2$summary.random$uv1$mean)-10),
          y = c(M1$summary.hyperpar$`0.975quant`[6], M1$summary.hyperpar$`0.975quant`[6],
                M1$summary.hyperpar$`0.025quant`[6], M1$summary.hyperpar$`0.025quant`[6]),
          col = "blue", density = 20)
}



## albumin Death event 
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
  

  #Plot the association 
  plot(M2$summary.random$uv1$mean, M2$summary.random$NL_CV_L1_S1$mean, main = "Albumin association: \n death event",
       xlab = "Linear predictor value", ylab = "multiplication factor", type = "p", pch = 20)
  
  abline(h = 0)
  abline(h = M1$summary.hyperpar$mean[6], col = "blue", lty = 4, lwd = 2) 
  
  polygon(x = c(min(M2$summary.random$uv1$mean)-10, max(M2$summary.random$uv1$mean)+10,
                max(M2$summary.random$uv1$mean)+10, min(M2$summary.random$uv1$mean)-10),
          y = c(M1$summary.hyperpar$`0.975quant`[6], M1$summary.hyperpar$`0.975quant`[6],
                M1$summary.hyperpar$`0.025quant`[6], M1$summary.hyperpar$`0.025quant`[6]),
          col = "blue", density = 20)
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
  

  #Plot the association 
  plot(M2$summary.random$uv1$mean, M2$summary.random$NL_CV_L1_S1$mean, main = "Albumin association: \n transplant event",
       xlab = "Linear predictor value", ylab = "multiplication factor", type = "p", pch = 20)
  
  abline(h = 0)
  abline(h = M1$summary.hyperpar$mean[6], col = "blue", lty = 4, lwd = 2) 
  
  polygon(x = c(min(M2$summary.random$uv1$mean)-10, max(M2$summary.random$uv1$mean)+10,
                max(M2$summary.random$uv1$mean)+10, min(M2$summary.random$uv1$mean)-10),
          y = c(M1$summary.hyperpar$`0.975quant`[6], M1$summary.hyperpar$`0.975quant`[6],
                M1$summary.hyperpar$`0.025quant`[6], M1$summary.hyperpar$`0.025quant`[6]),
          col = "blue", density = 20)
}



## edema Death event # Interesting !!!!
{
#Usual joint model with current value association
  M1 <- joint(formSurv = inla.surv(years, death) ~ drug,
              formLong = edema ~ (1 + year)*drug +
                (1 + year|id), family = "binomial",
              dataLong = Longi, dataSurv=Surv, id = "id", timeVar = "year",
              assoc = "CV", basRisk = "rw2", NbasRisk=25,
              control=list(int.strategy="eb"))
  summary(M1)

  #Joint model with nonlinear current value association
  M2 <- joint(formSurv = inla.surv(years, death) ~ drug,
              formLong = edema ~ (1 + year)*drug +
                (1 + year|id), family = "binomial",
              dataLong = Longi, dataSurv=Surv, id = "id", timeVar = "year",
              assoc = "NL_CV", basRisk = "rw2", NbasRisk=25,
              control=list(int.strategy="eb"))
  summary(M2)


  #Plot the association
  
  plot(M2$summary.random$uv1$mean, M2$summary.random$NL_CV_L1_S1$mean, main = "Edema association: \n death event",
       xlab = "Linear predictor value", ylab = "multiplication factor", type = "p", pch = 20)
  
  abline(h = 0)
  abline(h = M1$summary.hyperpar$mean[4], col = "blue", lty = 4, lwd = 2) 
  
  polygon(x = c(min(M2$summary.random$uv1$mean)-10, max(M2$summary.random$uv1$mean)+10,
                max(M2$summary.random$uv1$mean)+10, min(M2$summary.random$uv1$mean)-10),
          y = c(M1$summary.hyperpar$`0.975quant`[4], M1$summary.hyperpar$`0.975quant`[4],
                M1$summary.hyperpar$`0.025quant`[4], M1$summary.hyperpar$`0.025quant`[4]),
          col = "blue", density = 20)
  
  # plot(M2, NLeffectonly=TRUE)$NL_Association + 
  # 	geom_hline(yintercept = M1$summary.hyperpar$mean[4], linetype = "dashed", color = "blue") +
  # 	geom_ribbon(aes(ymin = M1$summary.hyperpar$`0.025quant`[4], ymax = M1$summary.hyperpar$`0.975quant`[4]), fill = "blue", alpha = 0.25)
  
}

## edema transplant event
{
#Usual joint model with current value association
  M1 <- joint(formSurv = inla.surv(years, trans) ~ drug,
              formLong = edema ~ (1 + year)*drug +
                (1 + year|id), family = "binomial",
              dataLong = Longi, dataSurv=Surv, id = "id", timeVar = "year",
              assoc = "CV", basRisk = "rw2", NbasRisk=25,
              control=list(int.strategy="eb"))
  summary(M1)

  #Joint model with nonlinear current value association
  M2 <- joint(formSurv = inla.surv(years, trans) ~ drug,
              formLong = edema ~ (1 + year)*drug +
                (1 + year|id), family = "binomial",
              dataLong = Longi, dataSurv=Surv, id = "id", timeVar = "year",
              assoc = "NL_CV", basRisk = "rw2", NbasRisk=25,
              control=list(int.strategy="eb"))
  summary(M2)


  #Plot the association
  plot(M2$summary.random$uv1$mean, M2$summary.random$NL_CV_L1_S1$mean, main = "Edema association: \n transplant event",
       xlab = "Linear predictor value", ylab = "multiplication factor", type = "p", pch = 20)
  
  abline(h = 0)
  abline(h = M1$summary.hyperpar$mean[4], col = "blue", lty = 4, lwd = 2) 
  
  polygon(x = c(min(M2$summary.random$uv1$mean)-10, max(M2$summary.random$uv1$mean)+10,
                max(M2$summary.random$uv1$mean)+10, min(M2$summary.random$uv1$mean)-10),
          y = c(M1$summary.hyperpar$`0.975quant`[4], M1$summary.hyperpar$`0.975quant`[4],
                M1$summary.hyperpar$`0.025quant`[4], M1$summary.hyperpar$`0.025quant`[4]),
          col = "blue", density = 20)
  }



## platelets Death event # VERY interesting # 2

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

  
  #Plot the association 
  plot(M2$summary.random$uv1$mean, M2$summary.random$NL_CV_L1_S1$mean, main = "Platelets association: \n death event",
       xlab = "Linear predictor value", ylab = "multiplication factor", type = "p", pch = 20)
  
  abline(h = 0)
  abline(h = M1$summary.hyperpar$mean[6], col = "blue", lty = 4, lwd = 2) 
  
  polygon(x = c(min(M2$summary.random$uv1$mean)-10, max(M2$summary.random$uv1$mean)+10,
                max(M2$summary.random$uv1$mean)+10, min(M2$summary.random$uv1$mean)-10),
          y = c(M1$summary.hyperpar$`0.975quant`[6], M1$summary.hyperpar$`0.975quant`[6],
                M1$summary.hyperpar$`0.025quant`[6], M1$summary.hyperpar$`0.025quant`[6]),
          col = "blue", density = 20)
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
  

  #Plot the association 
  plot(M2$summary.random$uv1$mean, M2$summary.random$NL_CV_L1_S1$mean, main = "Platelets association: \n transplant event",
       xlab = "Linear predictor value", ylab = "multiplication factor", type = "p", pch = 20)
  
  abline(h = 0)
  abline(h = M1$summary.hyperpar$mean[6], col = "blue", lty = 4, lwd = 2) 
  
  polygon(x = c(min(M2$summary.random$uv1$mean)-10, max(M2$summary.random$uv1$mean)+10,
                max(M2$summary.random$uv1$mean)+10, min(M2$summary.random$uv1$mean)-10),
          y = c(M1$summary.hyperpar$`0.975quant`[6], M1$summary.hyperpar$`0.975quant`[6],
                M1$summary.hyperpar$`0.025quant`[6], M1$summary.hyperpar$`0.025quant`[6]),
          col = "blue", density = 20)
}



## alkaline Death event # kind of not really # 3
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
  
  #Plot the association 
  plot(M2$summary.random$uv1$mean, M2$summary.random$NL_CV_L1_S1$mean, main = "alkaline association: \n death event",
       xlab = "Linear predictor value", ylab = "multiplication factor", type = "p", pch = 20)
  
  abline(h = 0)
  abline(h = M1$summary.hyperpar$mean[6], col = "blue", lty = 4, lwd = 2) 
  
  polygon(x = c(min(M2$summary.random$uv1$mean)-10, max(M2$summary.random$uv1$mean)+10,
                max(M2$summary.random$uv1$mean)+10, min(M2$summary.random$uv1$mean)-10),
          y = c(M1$summary.hyperpar$`0.975quant`[6], M1$summary.hyperpar$`0.975quant`[6],
                M1$summary.hyperpar$`0.025quant`[6], M1$summary.hyperpar$`0.025quant`[6]),
          col = "blue", density = 20)
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
  
  #Plot the association 
  plot(M2$summary.random$uv1$mean, M2$summary.random$NL_CV_L1_S1$mean, main = "Alkaline association: \n transplant event",
       xlab = "Linear predictor value", ylab = "multiplication factor", type = "p", pch = 20)
  
  abline(h = 0)
  abline(h = M1$summary.hyperpar$mean[6], col = "blue", lty = 4, lwd = 2) 
  
  polygon(x = c(min(M2$summary.random$uv1$mean)-10, max(M2$summary.random$uv1$mean)+10,
                max(M2$summary.random$uv1$mean)+10, min(M2$summary.random$uv1$mean)-10),
          y = c(M1$summary.hyperpar$`0.975quant`[6], M1$summary.hyperpar$`0.975quant`[6],
                M1$summary.hyperpar$`0.025quant`[6], M1$summary.hyperpar$`0.025quant`[6]),
          col = "blue", density = 20)
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
  
  #Plot the association 
  plot(M2$summary.random$uv1$mean, M2$summary.random$NL_CV_L1_S1$mean, main = "Ascites association: \n death event",
       xlab = "Linear predictor value", ylab = "multiplication factor", type = "p", pch = 20)
  
  abline(h = 0)
  abline(h = M1$summary.hyperpar$mean[4], col = "blue", lty = 4, lwd = 2) 
  
  polygon(x = c(min(M2$summary.random$uv1$mean)-10, max(M2$summary.random$uv1$mean)+10,
                max(M2$summary.random$uv1$mean)+10, min(M2$summary.random$uv1$mean)-10),
          y = c(M1$summary.hyperpar$`0.975quant`[4], M1$summary.hyperpar$`0.975quant`[4],
                M1$summary.hyperpar$`0.025quant`[4], M1$summary.hyperpar$`0.025quant`[4]),
          col = "blue", density = 20)
}

## ascites transplant event # very interesting  # 4
														# Stupid warning ????
														# ?????
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
  
  #Plot the association 
  plot(M2$summary.random$uv1$mean, M2$summary.random$NL_CV_L1_S1$mean, main = "Ascites association: \n transplant event",
       xlab = "Linear predictor value", ylab = "multiplication factor", type = "p", pch = 20)
  
  abline(h = 0)
  abline(h = M1$summary.hyperpar$mean[4], col = "blue", lty = 4, lwd = 2) 
  
  polygon(x = c(min(M2$summary.random$uv1$mean)-10, max(M2$summary.random$uv1$mean)+10,
                max(M2$summary.random$uv1$mean)+10, min(M2$summary.random$uv1$mean)-10),
          y = c(M1$summary.hyperpar$`0.975quant`[4], M1$summary.hyperpar$`0.975quant`[4],
                M1$summary.hyperpar$`0.025quant`[4], M1$summary.hyperpar$`0.025quant`[4]),
          col = "blue", density = 20)
}


##Spiders Death event # interesting # 5
											#????
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

plot(M2$summary.random$uv1$mean, M2$summary.random$NL_CV_L1_S1$mean, main = "Spiders association: \n death event",
     xlab = "Linear predictor value", ylab = "multiplication factor", type = "p", pch = 20)

abline(h = 0)
abline(h = M1$summary.hyperpar$mean[4], col = "blue", lty = 4, lwd = 2) 

polygon(x = c(min(M2$summary.random$uv1$mean)-10, max(M2$summary.random$uv1$mean)+10,
              max(M2$summary.random$uv1$mean)+10, min(M2$summary.random$uv1$mean)-10),
        y = c(M1$summary.hyperpar$`0.975quant`[4], M1$summary.hyperpar$`0.975quant`[4],
              M1$summary.hyperpar$`0.025quant`[4], M1$summary.hyperpar$`0.025quant`[4]),
        col = "blue", density = 20)
}

##Spiders transplant event #  pretty interesting # 6
													 # STUPID WARNING ?????
													# ?????????
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
  
  plot(M2$summary.random$uv1$mean, M2$summary.random$NL_CV_L1_S1$mean, main = "Spiders association: \n transplant event",
       xlab = "Linear predictor value", ylab = "multiplication factor", type = "p", pch = 20)
  
  abline(h = 0)
  abline(h = M1$summary.hyperpar$mean[4], col = "blue", lty = 4, lwd = 2) 
  
  polygon(x = c(min(M2$summary.random$uv1$mean)-10, max(M2$summary.random$uv1$mean)+10,
                max(M2$summary.random$uv1$mean)+10, min(M2$summary.random$uv1$mean)-10),
          y = c(M1$summary.hyperpar$`0.975quant`[4], M1$summary.hyperpar$`0.975quant`[4],
                M1$summary.hyperpar$`0.025quant`[4], M1$summary.hyperpar$`0.025quant`[4]),
          col = "blue", density = 20)
}



