
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
  
  vc <- VarCorr(fit)
  
  if(inherits(fit, c("lmerMod", "lmerModLmerTest", "lme4"))){
    
    out <- as.matrix(Matrix::bdiag(vc))
    if(is.null(unlist(dimnames(out)))) {
      nm <- unlist(lapply(vc, function(x) attributes(x)$dimnames[1]))
      dimnames(out) <- list(nm, nm)
    }
    round(out, digits)
    
  } else if(inherits(fit, "lme") & length(fit$group) < 2) { round(getVarCov(fit), digits) } else { vc }
  
}
                        
#=======================================================================
                        
do_context <- function(data, context_vars, group_id){
  
  all_names <- names(data)
  
  id <- grep("id|group|grp", all_names, value = TRUE, ignore.case = TRUE)
  
  if(!all(group_id %in% all_names)) { 
    
    stop(paste(toString(dQuote(group_id)), "not found for 'group_id' in the 'data'.", if(length(id)>0) 
      paste("\nPossibly you meant to use", toString(dQuote(id)), "as 'group_id', no?")), call. = FALSE) 
    
    }
  
  ok <- context_vars %in% all_names
  
  if(!all(ok)) message(paste("\n", toString(dQuote(context_vars[!ok])), "not found in the 'data' thus ignored."))
  
  context_vars <- context_vars[ok] 
  
  dum_vars <- all_names[sapply(data, function(i)is.character(i)|is.factor(i))]
  
  dum_names <- context_vars[context_vars %in% dum_vars]
  
  is_dum <- length(dum_names) > 0
  
  num_names <- context_vars[!(context_vars %in% dum_vars)]
  
  is_num <- length(num_names) > 0
  
  
  if(is_num){
    data <- data %>%
      group_by(across(all_of(group_id))) %>%                          
      mutate(across(all_of(num_names), list(wthn = ~ . - mean(.), btw = ~ mean(.)))) %>% 
      as.data.frame()
  }
  
  
  if(is_dum){
    data <- data %>%
      dummy_cols(select_columns = dum_names) %>% 
      group_by(across(all_of(group_id))) %>%                          
      mutate(across(starts_with(paste0(dum_names, "_")), list(wthn = ~ . - mean(.), btw = ~ mean(.)))) %>% 
      as.data.frame()
  }
  
  return(data)
}

#=================================================================================================================================  
                               
                               
hetro_var <- function(fit){
  
    if(!inherits(fit, c("lme", "gls"))) stop("Only 'lme()' & 'gls()' models are accepted.", call. = FALSE)
  
  coef(fit$modelStruct$varStruct, uncons = FALSE, allCoef = TRUE)
  
}
     
#=================================================================================================================================

                               
rho_lme <- function(fit) {
  
  if(!inherits(fit, c("lme", "gls"))) stop("Only 'lme()' & 'gls()' models are accepted.", call. = FALSE)
  
  coef(fit$modelStruct$corStruct, uncons = FALSE, allCoef = TRUE)

}                               
 
                               
#=================================================================================================================================
                               
cov_str <- function(fit, cov = TRUE, time_var = "time", hlm = TRUE){
  
  rho <- rho_lme(fit)
  hetro <- hetro_var(fit)
  sig <- sigma(fit)
  dat <- getData(fit)
  if(!(time_var %in% names(dat))) stop("Your 'time_var' doesn't exist in your data.", call. = FALSE)
  time_vals <- unique(dat[[time_var]])
  
  if(is.null(rho) & is.null(hetro)) return(id_cor(fit, cov = cov, time_var = time_var, hlm = hlm))
  if(is.null(rho) & !is.null(hetro)) {
    
    res <- id_cor(fit, cov = cov, time_var = time_var, hlm = hlm)
    diag(res) <- sum_ranef_var(fit)+(sig*hetro)^2
    return(res)  
    
  }
  corm <- corMatrix(fit$modelStruct$corStruct)[[1]]
  
  res <- corm*sig^2+if(hlm)sum_ranef_var(fit) else 0 *if(is.null(hetro)) 1 else t(t(hetro))%*%t(hetro)
  if(!is.null(hetro)) diag(res) <- sum_ranef_var(fit)+(sigma(fit)*hetro_var(fit))^2
    
  if(!cov) res <- cov2cor(res)
  rownames(res) <- colnames(res) <- paste0(time_var,time_vals)
  return(res)
  
  }                               
                               
#=================================================================================================================================  
                               
id_cor <- function(fit, cov = TRUE, time_var = "time", hlm = TRUE){
  
  sig <- sigma(fit)^2
  vrs <- as.numeric(VarCorr(fit)[,"Variance"])
  time_vals <- unique(getData(fit)[[time_var]])
  steps <- length(time_vals)
  x <- diag(steps)
  res <- sig*x
  
  diag(res) <- sum_ranef_var(fit, resid = TRUE) 
  res[col(res)!=row(res)] <- if(hlm) sum_ranef_var(fit) else 0
  
  rownames(res) <- colnames(res) <- paste0(time_var,time_vals)
  if(!cov) res <- cov2cor(res)
  return(res)
} 
                               
#=================================================================================================================================
           
sum_ranef_var <- function(fit, resid = FALSE){
vrs <- as.numeric(VarCorr(fit)[,"Variance"])
sum(if(resid) vrs else rev(vrs)[-1], na.rm = TRUE)
}                               
                                                       

    
#===================================================================================================================================
    
    
# Some hack to turn off unnneeded tick mark on the 3rd and 4th axes of plot effects
                                    
plot.efflist <- function (x, selection, rows, cols, graphics = TRUE, 
                          lattice, ...) 
{
  lattice <- if (missing(lattice)) 
    list()
  else lattice
  if (!missing(selection)) {
    if (is.character(selection)) 
      selection <- gsub(" ", "", selection)
    return(plot(x[[selection]], lattice = lattice, ...))
  }
  effects <- gsub(":", "*", names(x))
  neffects <- length(x)
  mfrow <- mfrow(neffects)
  if (missing(rows) || missing(cols)) {
    rows <- mfrow[1]
    cols <- mfrow[2]
  }
  for (i in 1:rows) {
    for (j in 1:cols) {
      if ((i - 1) * cols + j > neffects) 
        break
      more <- !((i - 1) * cols + j == neffects)
      lattice[["array"]] <- list(row = i, col = j, 
                                 nrow = rows, ncol = cols, more = more)
      pp <- plot(x[[(i - 1) * cols + j]], lattice = lattice, 
                 ...)
      # hack to turn off opposite side tick marks
      pp$x.scales$tck=c(1,0)
      pp$y.scales$tck=c(1,0)
      print(pp)
    }
  }
}
environment(plot.efflist) <- asNamespace("effects")    
    
    
#========================================================================                        
                                             
needzzsf <- c('car','psych','reshape','tidyverse','lme4','nlme','MASS','CCA','matrixcalc', 'mvoutlier', 'vegan', 'haven', 'fastDummies', "emmeans", 'sjPlot', 'lmerTest', 'reghelper',
          'parallel','rela','gplots','ICSNP','mvtnorm','mvnormtest','normtest', 'micompr', 'heplots', 'HSAUR', 'bbmle', 'jtools', 'stargazer', 'interactions',
          'normwhn.test','nortest','biotools','effects','ez','yacca')

                        
not.have23 <- needzzsf[!(needzzsf %in% installed.packages()[,"Package"])]
if(length(not.have23)) install.packages(not.have23)


suppressWarnings(
suppressMessages({ 
  
  for(i in needzzsf){
    library(i, character.only = TRUE)
  }
}))
    
source('https://raw.githubusercontent.com/rnorouzian/A/main/a.r')    
