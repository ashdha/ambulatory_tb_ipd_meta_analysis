# TB screening in ambulatory PLHIV study
# amb_plhiv_screeing_code.R
# University of Cape Town

# 5. Additional diagnostic accuracy estimates (Table S8)
     #Datasets: 
     #4_sens_spec_precalc_new.csv


# 5. Additional diagnostic accuracy estimates (Table S8) -------------------------------------------------------------------
    
    library(tidyverse)
    library(dplyr)
    library(mada)
    library(meta)
    library(lme4)
    library(msm)
    library(altmeta)
    library(lmtest)
    library(reshape2)
    library(Hmisc)
    
    #yy is  subgroup and vv is test
    #e.g 
    #yy <- "all"
    #vv <- "who"
    
    merged  = read.csv("4_sens_spec_precalc_new.csv", stringsAsFactors = FALSE)
    ind_test <- subset(merged, merged$test==vv & merged$ID == yy)
    
#calculate univariate DOR and perform eggers test, minimum 3 studies##############################
    x <- metabin(tp, n1, fp, n2, data=ind_test, sm="OR", 
                 studlab=paste(author), addincr=FALSE, allstudies=TRUE)
    summary(x)
    e <- metabias(x, method.bias = "linreg", k.min =3) #minimum 3 studies
    print(e)
    
#perform trim-and-fill analysis,  minimum 3 studies##############################################
    if (x$k >2){tf1 <- trimfill(x, comb.fixed=TRUE)}
    if (x$k <=2){tf1 <- list(TE.random = NA, lower.random = NA, upper.random = NA)}
    summary(tf1)
    
#positive screen meta-anlaysis####################################################################
    prop.calc <- metaprop(pos, tot, data=ind_test, sm = "PLOGIT", method = "GLMM", studlab=paste(author))
    summary(prop.calc)
    
#metaregression####################################################################################
    #ref_std or prev_std
    metareg(x, ref_std)
    metareg(x, prev_std)
    

#alternative PPV/NPV calculation methods (trivariate GLMM)#########################################
    #adapted from cochrane tutorial
    
    #yy is  subgroup and vv is test
    #e.g 
    yy <- "all"
    vv <- "who"
    
    #read in dataset
    ind_test <- subset(merged, merged$test==vv & merged$ID == yy)
    
#if >=2 studies available#
    
    X <- ind_test
    ## In order to specify the generalized linear model, first, we need to set up the data 
    
    ## Set up the data
    ## Generate 5 new variables of type long before reshaping the data.
    # n2 is total number
    # n1 is number positive test
    # n0 is number with negative test
    # true2 is the number of test positives
    # true1 is number of true positives
    # true0 is the number of true negatives
    # study is the unique identifier for each study. _n will generate a sequence of numbers. 
    
    X$n2 <- X$tot
    X$n1 <- X$tp+X$fp
    X$n0 <- X$fn+X$tn
    X$true2 <- X$tp+X$fp
    X$true1 <- X$tp
    X$true0 <- X$tn 
    X$study <- 1:length(X$author)
    
    X %>%
      dplyr::select(-ID, -test) -> X
    
    Y = reshape(X, direction = "long", varying = list( c("n1" , "n0", "n2") , c( "true1","true0", "true2" ) ) ,
                timevar = "ppv" , times = c(1,0,2) , v.names = c("n","true") )
    
    ## Sort data by study to cluster the 2 records per study together ###
    Y = Y[order(Y$id),]  
    Y$npv<- 1-Y$ppv
    Y$prev <- Y$ppv*Y$npv
    Y$prev [ Y$prev==-2] <- 1
    Y$npv [ Y$npv==-1] <- 0
    Y$ppv [ Y$ppv==2] <- 0
    
    ## Perform meta-analysis ## 
    (MA_Y = glmer( formula = cbind(true , n - true ) ~ 0 + ppv + npv + prev + (0+ppv + npv + prev|study),
                   data = Y , family = binomial, control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 2e5)))) 
    
    
    ma_Y = summary(MA_Y)
    ma_Y$coefficients
    #confint(MA_Y, method=c("boot"), boot.type=c("basic"))
    ## Therefore, to extract the coefficients 
    ma_Y$coeff 
    (lppv = ma_Y$coeff[1,1])
    (lnpv = ma_Y$coeff[2,1]) 
    (lprev = ma_Y$coeff[3,1]) 
    
    se.lppv = ma_Y$coeff[1,2]
    se.lnpv = ma_Y$coeff[2,2] 
    se.lprev = ma_Y$coeff[3,2] 
    
    ## Then we can manually create 95% confidence intervals for logit ppv and npv
    ppv = c(lppv, lppv-qnorm(0.975)*se.lppv, lppv+qnorm(0.975)*se.lppv ) 
    npv = c(lnpv, lnpv-qnorm(0.975)*se.lnpv, lnpv+qnorm(0.975)*se.lnpv ) 
    Prev = c(lprev, lprev-qnorm(0.975)*se.lprev, lprev+qnorm(0.975)*se.lprev ) 
    
    ## R has a built in logit and inv.logit function (use qlogis and plogis )
    #make dataframe with n, number studies, ppv, npv
    
    no_st <- as.data.frame (ma_Y$ngrps)
    tot <- as.data.frame (sum(X$tot))
    
    ppv.test <- 100*t(plogis(ppv ))
    ppv.test <- as.data.frame(ppv.test)
    ppv.test %>%  
      dplyr::rename("ppv" = V1, "ppv_lb" = V2, "ppv_ub" = V3)  -> ppv.test
    #ppv.test <- round(ppv.test,digits=1)
    
    npv.test <-  100*t(plogis( npv ))
    npv.test <- as.data.frame(npv.test)
    npv.test %>%  
      dplyr::rename("npv" = V1, "npv_lb" = V2, "npv_ub" = V3)  -> npv.test
    #npv.test <- round(npv.test,digits=1)
    
    prev.test <-  100*t(plogis( Prev ))
    prev.test <- as.data.frame(prev.test)
    prev.test %>%  
      dplyr::rename("prev" = V1, "prev_lb" = V2, "prev_ub" = V3)  -> prev.test
    #prev.test <- round(prev.test,digits=1))

    trivariate.results <- cbind(ppv.test, npv.test, prev.test)
    trivariate.results

#results if only 1 study available   
    X <- ind_test
    
    X$n2 <- X$tot # test prev
    X$n1 <- X$tp+X$fp # ppv
    X$n0 <- X$fn+X$tn # npv
    X$true2 <- X$tp+X$fp # test prev
    X$true1 <- X$tp # ppv
    X$true0 <- X$tn # npv
    X$study <- 1:length(X$author)
    
    
    ppv <- X$true1/X$n1
    npv <- X$true0/X$n0
    prev <- X$true2/X$n2
    
    #calculate CIs using wilson method
    ppv.test <- 100*binconf(X$true1, X$n1)
    npv.test <- 100*binconf(X$true0, X$n0)
    prev.test <- 100*binconf(X$true2, X$n2)


    trivariate.results <- data.frame(ppv = ppv.test[1], ppv_lb = ppv.test[2], ppv_ub = ppv.test[3], 
                                     npv = npv.test[1], npv_lb = npv.test[2], npv_ub = npv.test[3],
                                     prev = prev.test[1], prev_lb = prev.test[2], prev_ub = prev.test[3])
    
    trivariate.results
    
    
    
    