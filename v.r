
Break = "\n<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>\n"

notice = "   Welcome to 'EDP 380C.12 Survey of Multivariate Methods'.
   Programs developed by Reza Norouzian, Copyright (C) 2019-present"

message(Break, notice, Break)


##----------------------------------------- Data Generating Function -------------------------------------------------->
gen_mv_data <- function(n, mu1, mu2, mu3 = NULL, mu4 = NULL, cov1, cov2 = cov1, cov3 = NULL, cov4 = NULL, k = 2) {
  if(k == 2) {
  group1 <- as.data.frame(mvrnorm(n, mu1, cov1))
  group1$k <- 1
  group2 <- as.data.frame(mvrnorm(n, mu2, cov2))
  group2$k <- 2
  
   data <- rbind(group1, group2)
   
  } else if(k == 3) {
    group1 <- as.data.frame(mvrnorm(n, mu1, cov1))
    group1$k <- 1
    group2 <- as.data.frame(mvrnorm(n, mu2, cov2))
    group2$k <- 2
    group3 <- as.data.frame(mvrnorm(n, mu3, cov3))
    group3$k <- 3
    
    data <- rbind(group1, group2, group3)
    
  } else if(k == 4) {
    group1 <- as.data.frame(mvrnorm(n, mu1, cov1))
    group1$k <- 1
    group2 <- as.data.frame(mvrnorm(n, mu2, cov2))
    group2$k <- 2
    group3 <- as.data.frame(mvrnorm(n, mu3, cov3))
    group3$k <- 3
    group4 <- as.data.frame(mvrnorm(n, mu4, cov4))
    group4$k <- 4
    
    data <- rbind(group1, group2, group3, group4)
    
    
} else NULL
}  
## Example to generate 2 DV's for 2 groups
# ex1 <- gen_mv_data(n = 5, mu1 = c(1,2), mu2 = c(2,3), cov1 = matrix(c(3,1,1,3),2,2), cov2 = matrix(c(5,2,2,5),2,2))



## Example to generate 3 DV's for 4 groups with 10 subjects each
# ex2 <- gen_mv_data(n = 10, mu1 = c(1,4,1), mu2 = c(8,9,1), mu3 = c(4,8,1), mu4 = c(1,4,1),
#                    cov1 = matrix(c(3,1,1,1,3,1,1,1,3),3,3), cov2 = matrix(c(3,1,1,1,3,1,1,1,3),3,3), cov3 = matrix(c(3,1,1,1,3,1,1,1,3),3,3),
#                    cov4 = matrix(c(3,1,1,1,3,1,1,1,3),3,3), k = 4)


#=======================================================================

G_matrix <- function(fit, digits = 8){

if(inherits(fit, c("lmerMod", "lmerModLmerTest", "lme4"))){
  
  vc <- VarCorr(fit)
  out <- as.matrix(Matrix::bdiag(vc))
  if(is.null(unlist(dimnames(out)))) {
    nm <- unlist(lapply(vc, function(x) attributes(x)$dimnames[1]))
    dimnames(out) <- list(nm, nm)
  }
  round(out, digits)
  
} else if(inherits(fit, "lme")) { round(getVarCov(fit), digits) }

}
                        
#=======================================================================                        
                                             
needzzsf <- c('car','psych','reshape','tidyverse','lme4','nlme','MASS','CCA','matrixcalc', 'mvoutlier', 'vegan', 'haven',
          'parallel','rela','gplots','ICSNP','mvtnorm','mvnormtest','normtest', 'micompr', 'heplots', 'HSAUR', 'bbmle',
          'normwhn.test','nortest','biotools','effects','ez','yacca')

                        
not.have23 <- needzzsf[!(needzzsf %in% installed.packages()[,"Package"])]
if(length(not.have23)) install.packages(not.have23)


suppressWarnings(
suppressMessages({ 
  
  for(i in needzzsf){
    library(i, character.only = TRUE)
  }
}))
