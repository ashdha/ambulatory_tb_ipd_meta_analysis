# TB screening in ambulatory PLHIV study
# amb_plhiv_screeing_code.R
# University of Cape Town


# 3. Diagnostic test accuracy (direct comparisons)
     #Datasets: 
     #7_head_to_head_meta_a_new.csv

# 3. Diagnostic test accuracy (direct comparisons) -------------------------------------------------------------------
    
    library(reshape2)
    library(tidyverse)
    library(dplyr)
    library(mada)
    library(meta)
    library(lme4)
    library(msm)
    library(altmeta)
    library(lmtest)
    
    
    merged  = read.csv("7_head_to_head_precalc.csv", stringsAsFactors = FALSE)
    merged %>%
      dplyr::select(ID, author, TP, FP, FN, TN, tot, test, who_status) -> merged
    
    
  #difference from W4SS using bivariate meta-regression with test type as covariate to calculate accuracy estimates and p-value      

  #calculation for >=2 studies (bivariate GLMM)
    
    #e.g.
    #yy <- "all"
    #vv <- "crp_10"
    X <- subset(merged, merged$ID == yy)
    X <- subset(X, X$test==vv)
    
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
    N <- length(X$TP)
    X$study <- 1:N
    
    ## Reshape the data from wide to long format ###
    Y = reshape(X, direction = "long", varying = list( c("n1" , "n0") , c( "true1","true0" ) ) ,
                timevar = "sens" , times = c(1,0) , v.names = c("n","true") ) 
    
    
    ## Sort data by study to cluster the 2 records per study together ###
    
    Y = Y[order(Y$id),]  
    Y$spec<- 1-Y$sens
    
    ### Add covariate terms to the model for both logit sensitivity and logit specificity. 
    ### This model assumes equal variances for both tests. ***
    
    Y$who_status <- as.factor(Y$who_status)                            # Convert character column to factor
    
    Y$T1  <- 2 - as.numeric(Y$who_status) #let T1 be test to compare & T2 who screen
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
    
    #sensitivity and specificity estimates
    cB = summary(B)
    
    lsensT1 = cB$coeff[1,1]
    lspecT1 = cB$coeff[3,1]
    lsensT2 = cB$coeff[2,1]
    lspecT2 = cB$coeff[4,1]
    
    se.lsensT1 = cB$coeff[1,2]
    se.lspecT1 = cB$coeff[3,2]
    se.lsensT2 = cB$coeff[2,2]
    se.lspecT2 = cB$coeff[4,2] 
    
    ## Then we can manually create 95% confidence intervals for logit sens and spec
    logit_SensT1 = c(lsensT1, lsensT1-qnorm(0.975)*se.lsensT1, lsensT1+qnorm(0.975)*se.lsensT1 ) 
    logit_SensT2 = c(lsensT2, lsensT2-qnorm(0.975)*se.lsensT2, lsensT2+qnorm(0.975)*se.lsensT2 ) 
    
    logit_SpecT1 = c(lspecT1, lspecT1-qnorm(0.975)*se.lspecT1, lspecT1+qnorm(0.975)*se.lspecT1 ) 
    logit_SpecT2 = c(lspecT2, lspecT2-qnorm(0.975)*se.lspecT2, lspecT2+qnorm(0.975)*se.lspecT2 ) 
    
    SensT1 = plogis( logit_SensT1 )
    SensT2 = plogis( logit_SensT2 )
  
    SpecT1 = plogis( logit_SpecT1 ) 
    SpecT2 = plogis( logit_SpecT2 ) 
    
    SensT1
    SpecT1
    
    SensT2
    SpecT2
    
#calculation for when only 1 study - individual sensitivities and specificities with 95% CI constructed using the wilson method in mada package
    
    #e.g.
    #yy <- "art"
    #vv <- "crp_10"
    X <- subset(merged, merged$ID == yy)
    X <- subset(X, X$test==vv)
    
    #p values#
    d <- reitsma(X, formula = cbind(tsens, tfpr) ~ who_status)
    e <- summary(d)
    print(e)
    
    #fit model for test 1#
    subgroup_test1 <- subset(X, X$who_status=="no")
    mrfit <- reitsma(subgroup_test1)
    mrfit.T1 <- summary(mrfit, correction = 0, correction.control = "all", digits=2)
    mrfit.T1$coefficients[3,] #sensitivity
    1-mrfit.T1$coefficients[4,] #specificity
    
    #fit model for test 2 (who)#
    subgroup_test2 <- subset(X, X$who_status=="yes")
    mrfit1 <- reitsma(subgroup_test2)
    mrfit.T2 <- summary(mrfit1, correction = 0, correction.control = "all", digits=2)
    mrfit.T2$coefficients[3,] #sensitivity
    1-mrfit.T2$coefficients[4,] #specificity
    
    
