# TB screening in ambulatory PLHIV study
# amb_plhiv_screeing_code.R
# University of Cape Town

# 2. Diagnostic test accuracy (indirect comparisons), including yield calculation and sensitivity analyses
     #Datasets: 
     #4_sens_spec_precalc_new.csv 
          #main analysis with culture as reference standard:
     #4_sens_spec_precalc_new_cult_xpert.csv 
          #sensitivity analysis with culture or Xpert as reference standard
     #4_sens_spec_precalc_new_xpert_conf_test.csv
          #sensitivity analysis with only Xpert as reference standard
     #7_head_to_head_crp_algorithm_precalc.csv
          #sensitivity analysis comparing W4SS followed by Xpert with CRP (>=10mg) followed by Xpert

# 2. Diagnostic test accuracy (indirect comparisons), including yield calculation and sensitivity analyses -------------------------------------------------------------------

    library(reshape2)
    library(tidyverse)
    library(dplyr)
    library(mada)
    library(meta)
    library(lme4)
    library(msm)
    library(altmeta)
    library(lmtest)
    #read in dataset depending on analysis
    
    #main analysis with culture as reference standard ####################################################################################
    merged  = read.csv("4_sens_spec_precalc_new.csv", stringsAsFactors = FALSE)
    
    #sensitivity analysis with culture or Xpert as reference standard #####################################################################
    #merged  = read.csv("4_sens_spec_precalc_new_cult_xpert.csv", stringsAsFactors = FALSE)
    
    #sensitivity analysis with only Xpert as reference standard ###########################################################################
    #merged  = read.csv("4_sens_spec_precalc_new_xpert_conf_test.csv", stringsAsFactors = FALSE)
    
    #sensitivity analysis comparing W4SS followed by Xpert with CRP (>=10mg) followed by Xpert ############################################
    #merged  = read.csv("7_head_to_head_crp_algorithm_precalc.csv", stringsAsFactors = FALSE)
    
    #subset group and test
    #e.g
    #X <- subset(merged, merged$test=="who" & merged$ID == "all")
    X <- subset(merged, merged$test==vv & merged$ID == yy)
    
#calculation for >=2 studies (bivariate GLMM and model assuming no correlation between sensitivity and specificity to simplify the model when non-convergence )
      N <- length(X$TP)
      number_of_people <- sum(X$tot)
      
      ## In order to specify the generalized linear model, first, we need to set up the data 
      ## Set up the data
      ## Generate 5 new variables of type long. We need these before we can reshape the data.
      # ÃÂ	n1 is number diseased
      # ÃÂ	n0 is number without disease
      # ÃÂ	true1 is number of true positives
      # ÃÂ	true0 is the number of true negatives
      # ÃÂ	study is the unique identifier for each study. _n will generate a sequence of numbers. 
      
      X$n1 <- X$TP+X$FN
      X$n0 <- X$FP+X$TN
      X$true1 <- X$TP
      X$true0 <- X$TN 
      X$study <- 1:N
      
      ## Reshape the data from wide to long format ###
      Y = reshape(X, direction = "long", varying = list( c("n1" , "n0") , c( "true1","true0" ) ) ,
                  timevar = "sens" , times = c(1,0) , v.names = c("n","true") ) 
      
      ## Sort data by study to cluster the 2 records per study together ###
      Y = Y[order(Y$id),]  
      Y$spec<- 1-Y$sens
      
      ## Perform meta-analysis ## 
      MA_Y = glmer( formula = cbind(  true , n - true ) ~ 0 + sens + spec + (0+sens + spec|study),
                    data = Y , family = binomial , nAGQ=1 , verbose=0, control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 2e5)))
      
      #In the case of non-convergence, assume no correlation between measures of sensitivity and specificity to simplify the model using this approach:
      #MA_Y = glmer( formula = cbind( true , n - true ) ~ 0 + sens + spec + (0+sens|study) + (0+spec|study),
      #              data = Y , family = binomial , nAGQ = 1 , verbose = 0 , control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 2e5)))
      
      ma_Y = summary(MA_Y)
      
      lsens = ma_Y$coeff[1,1]
      lspec = ma_Y$coeff[2,1]
      
      se.lsens = ma_Y$coeff[1,2]
      se.lspec = ma_Y$coeff[2,2] 
      
      ## Then we can manually create 95% confidence intervals for logit sens and spec
      logit_Sens = c(lsens, lsens-qnorm(0.975)*se.lsens, lsens+qnorm(0.975)*se.lsens ) 
      logit_Spec = c(lspec, lspec-qnorm(0.975)*se.lspec, lspec+qnorm(0.975)*se.lspec ) 
      Sens = plogis( logit_Sens ) 
      Spec = plogis( logit_Spec ) 

      Sens
      Spec

#calculation if only 1 study - individual sensitivities and specificities with 95% CI constructed using the wilson method
      
      mrfit <- reitsma(X)
      mrfit.1 <- summary(mrfit, correction = 0, correction.control = "all", digits=2)
      mrfit.1$coefficients[3,] #sensitivity
      1-mrfit.1$coefficients[4,] #specificity

#difference from W4SS p-value using bivariate meta-regression with test type as covariate      
      
        #subset group and test to compare with w4ss
        #e.g 
        #yy <- "all"
        #vv <- "crp_10"
      
        X <- subset(merged, merged$ID == yy)
        X <- subset(X, X$test=="who" | X$test==vv)
      
        N <- length(X$TP)
        number_of_people <- sum(X$tot)
        
        ## In order to specify the generalized linear model, first, we need to set up the data 
        ## Set up the data
        ## Generate 5 new variables of type long. We need these before we can reshape the data.
        # ÃÂ	n1 is number diseased
        # ÃÂ	n0 is number without disease
        # ÃÂ	true1 is number of true positives
        # ÃÂ	true0 is the number of true negatives
        # ÃÂ	study is the unique identifier for each study. _n will generate a sequence of numbers. 
        
        X$n1 <- X$TP+X$FN
        X$n0 <- X$FP+X$TN
        X$true1 <- X$TP
        X$true0 <- X$TN 
        X$study <- 1:N
        
        ## Reshape the data from wide to long format ###
        Y = reshape(X, direction = "long", varying = list( c("n1" , "n0") , c( "true1","true0" ) ) ,
                    timevar = "sens" , times = c(1,0) , v.names = c("n","true") ) 
        
        
        ## Sort data by study to cluster the 2 records per study together ###
        
        Y = Y[order(Y$id),]  
        Y$spec<- 1-Y$sens
        
        
        ### Fit the model without the covariate ###
        
        #(A = glmer( formula = cbind(  true , n - true ) ~ 0 + sens + spec + (0+sens + spec|study),
        #            data = Y , family = binomial  ) ) 
        
        ### Add covariate terms to the model for both logit sensitivity and logit specificity. 
        ### This model assumes equal variances for both tests. ***
        
        Y$test <- as.factor(Y$test)                            # Convert character column to factor
        
        Y$T1  <- 2 - as.numeric(Y$test) 
        Y$T2 <- 1 - Y$T1 
        
        Y$seT1  <- (Y$T1)*(Y$sens) 
        Y$seT2 <- (Y$T2)*(Y$sens) 
        
        Y$spT1  <- (Y$T1)*(Y$spec) 
        Y$spT2 <- (Y$T2)*(Y$spec) 
        
        (B = glmer( formula = cbind(  true , n - true ) ~ 0 + seT1 + seT2 + spT1 + spT2 + 
                      (0+sens + spec|study), data = Y , family = binomial, nAGQ=1 , verbose=0, control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 2e5)))) 
        
        ## Is there a statistically significant difference in sensitivity between T1 and T2?  
        
        (C = glmer( formula = cbind(  true , n - true ) ~ 0 + sens + spT1 + spT2 + 
                      (0+sens + spec|study), data = Y , family = binomial,nAGQ=1 , verbose=0, control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 2e5)))) 
        
        sens_diff <- lrtest(B,C)
        sens_diff
        ## Is there a statistically significant difference in specificity between T1 and T2?  
        
        (D = glmer( formula = cbind(  true , n - true ) ~ 0 + seT1 + seT2 + spec + 
                      (0+sens + spec|study), data = Y , family = binomial,nAGQ=1 , verbose=0, control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 2e5)))) 
        
        lrtest(B,D)
        spec_diff <- lrtest(B,D)    
        spec_diff
        
        ##when non-convergence, assume no correlation between sensitivity and specificity to simplify the model as below
        ##
        ##
        #(B = glmer( formula = cbind(  true , n - true ) ~ 0 + seT1 + seT2 + spT1 + spT2 + (0+sens|study) + (0+spec|study)
        #            , data = Y , family = binomial, nAGQ=1 , verbose=0, control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 2e5)))) 
        
        ## Is there a statistically significant difference in sensitivity between T1 and T2?  
        
        #(C = glmer( formula = cbind(  true , n - true ) ~ 0 + sens + spT1 + spT2 + (0+sens|study) + (0+spec|study)
        #            , data = Y , family = binomial,nAGQ=1 , verbose=0, control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 2e5)))) 
        
        #sens_diff <- lrtest(B,C)
        
        ## Is there a statistically significant difference in specificity between T1 and T2?  
        
        #(D = glmer( formula = cbind(  true , n - true ) ~ 0 + seT1 + seT2 + spec + (0+sens|study) + (0+spec|study)
        ##            , data = Y , family = binomial,nAGQ=1 , verbose=0, control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 2e5)))) 
        
        #lrtest(B,D)
        #spec_diff <- lrtest(B,D)
        
        
#Yield calculation####################################################################################
    
    #first make datframe with sensitivity and specificity of test, then calculate other values
    #e.g for W4SS
    #yy <- 0.82 
    #vv <- 0.42
    yield_calc = data.frame(sensitivity = yy, specificity = vv)
    
  

    #values are per 1000#  
    #Calculate number screened positive, TP, FP, PPV, TN, FN, NPV#
    #let y be the prevalence e.g y <- 0.1#
    y <- 0.1
    
    yield_calc$numberscreen_pos <- (yield_calc$sensitivity*(y)*1000) + ((1-y)*1000*(1-(yield_calc$specificity))) #if applicable to calculate
    yield_calc$numberscreen_pos <- round(yield_calc$numberscreen_pos, digits=0)
    
    yield_calc$true_positive <- ((y)*1000*(yield_calc$sensitivity))
    yield_calc$true_positive <- round(yield_calc$true_positive, digits=0)
    
    yield_calc$false_positive <- ((1-y)*1000*(1-(yield_calc$specificity)))
    yield_calc$false_positive <- round(yield_calc$false_positive, digits=0)
    
    yield_calc$ppv <- 100*(yield_calc$sensitivity*(y))/(((1-yield_calc$specificity)*(1-y)) +(yield_calc$sensitivity*y))
    yield_calc$ppv <- round(yield_calc$ppv,digits=1)
    
    yield_calc$true_negative <- (((1-y)*1000)-yield_calc$false_positive)
    yield_calc$true_negative <- round(yield_calc$true_negative, digits=0)
    
    yield_calc$false_negative <- ((y*1000)-yield_calc$true_positive)
    yield_calc$false_negative <- round(yield_calc$false_negative, digits=0)
    
    yield_calc$npv <- 100*(yield_calc$specificity*(1-y))/(((1-yield_calc$sensitivity)*y) + (yield_calc$specificity*(1-y)))
    yield_calc$npv <- round(yield_calc$npv,digits=1)
    
    yield_calc$NNS <- 1000/yield_calc$true_positive #if applicable to calculate
    yield_calc$NNS <- round(yield_calc$NNS,digits=0)
    
