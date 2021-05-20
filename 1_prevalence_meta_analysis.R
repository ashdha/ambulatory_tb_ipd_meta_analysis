# TB screening in ambulatory PLHIV study
# amb_plhiv_screeing_code.R
# University of Cape Town

# 1. Prevalence meta-analysis
     #Datasets: 
     #table_prev_groups.csv

# 1. Prevalence meta-analysis -------------------------------------------------------------------
        
    library(meta)
    library(tidyverse)
    #read in dataset
    subgroup <- read.csv("table_prev_groups.csv")

    #prevelance analysis for culture as definition of TB###########################################

    #subset relevant subgroup,
    #e.g 
    #x  <- subset(subgroup, id == "art")
    x  <- subset(subgroup, id == y)
    
    #perform meta-analysis
    x <- metaprop(r1, n1, data=x, sm = "PLOGIT", method = "GLMM", studlab=paste(author))
    summary(x)
    #perform eggers test, minimun 3 studies
    e <- metabias(x, method.bias = "linreg", k.min =3)
    print(e)
    #subgroup analysis, 
    #e.g 
    #v <- "cd4_high_200"
    #w <- "cd4_less_200"
    v  <- subset(subgroup, id == v)
    w  <- subset(subgroup, id == w)
    x <- rbind(v, w )
    
    y <- metaprop(r1, n1, data=x, sm = "PFT", byvar=id, studlab=paste(author))
    summary(y)
    
    #prevelance analysis for culture or xpert as definition of TB##################################
    
    #subset relevant subgroup,
    #e.g 
    #x  <- subset(subgroup, id == "all")
    x  <- subset(subgroup, id == y)
    
    #perform meta-analysis
    x <- metaprop(r2, n2, data=x, sm = "PLOGIT", method = "GLMM", studlab=paste(author))
    summary(x)
    #perform eggers test, minimun 3 studies
    e <- metabias(x, method.bias = "linreg", k.min =3)
    print(e)
    #subgroup analysis, 
    #e.g 
    v <- "cd4_high_200"
    w <- "cd4_less_200"
    v  <- subset(subgroup, id == v)
    w  <- subset(subgroup, id == w)
    x <- rbind(v, w )
    
    y <- metaprop(r1, n1, data=x, sm = "PFT", byvar=id, studlab=paste(author))
    summary(y)
    
