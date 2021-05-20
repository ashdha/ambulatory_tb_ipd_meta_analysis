# TB screening in ambulatory PLHIV study
# amb_plhiv_screeing_code.R
# University of Cape Town

# 4. Forest plots, ROC curves, and Funnel plots
     #Datasets: 
     #4_sens_spec_precalc_new.csv
    
# 4. Forest plots, ROC curves, and Funnel plots -------------------------------------------------------------------
    
    library(tidyverse)
    library(mada)
    library(meta)
    library(lme4)
    library(altmeta)

#Forest plot####################################################################################

    subgroup  = read.csv("4_sens_spec_precalc_new.csv", stringsAsFactors = FALSE)

    #remove those rows of 2x2 tables with >=3 cells with zero events, as not compatible with mada#
    subgroup1 <- subset(subgroup, subgroup$TP ==0 & subgroup$TN ==0 & subgroup$FP ==0 & subgroup$FN ==0)
    subgroup2 <- subset(subgroup, subgroup$TP ==0 & subgroup$TN ==0 & subgroup$FP ==0)
    subgroup3 <- subset(subgroup, subgroup$TP ==0 & subgroup$TN ==0 & subgroup$FN ==0)
    subgroup4 <- subset(subgroup, subgroup$TP ==0 & subgroup$FP ==0 & subgroup$FN ==0)
    subgroup5 <- subset(subgroup, subgroup$TN ==0 & subgroup$FP ==0 & subgroup$FN ==0)
    subgroup <- setdiff(subgroup, subgroup1)
    subgroup <- setdiff(subgroup, subgroup2)
    subgroup <- setdiff(subgroup, subgroup3)
    subgroup <- setdiff(subgroup, subgroup4)
    subgroup <- setdiff(subgroup, subgroup5)
    rm("subgroup1", "subgroup2", "subgroup3", "subgroup4", "subgroup5")
    
    #subset relevant subgroup and test
    #e.g a  <- subset(subgroup, ID == "all" & test=="who")
    a <- subset(z, z$test==v)
    
    #forest plot for sensitivity and specificity
    mada::forest(madad(a), type = "sens", main="Sensitivity", cex=0.8, snames = a$author)
    mada::forest(madad(a), type = "spec", main="Specificity", cex=0.8, snames = a$author)
    
#ROC currve######################################################################################

    merged  = read.csv("4_sens_spec_precalc_new.csv", stringsAsFactors = FALSE)
    
    #subset relevant subgroup and test
    #yy is  subgroup and vv is test for subset 
    #e.g 
    #yy <- "all"
    #vv <- "who"
    
    X <- subset(merged, merged$test==vv & merged$ID == yy)
    
    a <- meta.dt(tp, fp, fn, tn, X, method = "biv.glmm") #bivariate meta-analysis
      
    j.main <- paste0(yy, " - ROC curve for \n", vv)
    par(pty="s") #plot square region
      
    #plot roc curve
    plot(a, cex.studies = 0.8, pch.studies = 2, col.studies = "black",
           roc=TRUE, lwd.roc=2, col.overall = "forest green", pch.overall = 16, 
           cex.overall = 1.2, confid=FALSE, col.confid="red", main=j.main, cex.main=0.8)
    
#Funnel plot####################################################################################
    
    #read in dataframe
    subgroup  = read.csv("4_sens_spec_precalc_new.csv")
    
    #subset relevant subgroup and test
    #y is  subgroup and v is test
    #e.g 
    #y <- "all"
    #v <- "who"
    t <- subset(subgroup, subgroup$ID == "all" & subgroup$test== "who")
    
    ms1 <- metabin(tp, n1, n2-tn, n2, data=t, sm="OR", 
                       studlab=paste(author), addincr=FALSE, allstudies=TRUE)
    funnel(ms1, main = j)

