
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
                          
G_pca <- function(fit) {
  
  if(!inherits(fit, c("lmerMod", "lmerModLmerTest", "lme4"))) stop("Non-lme4 model detected.", call. = FALSE)
  obj <- summary(lme4::rePCA(fit))
  model <- lme4::VarCorr(fit)
  if(length(obj) == length(model)) {
    obj <- Map(function(x, z) {
      colnames(x$importance) <- paste(z, colnames(model[[z]]), sep = '_')
      x
    }, obj, names(obj))
  }
  else if(length(obj) == 1) {
    colnames(obj[[1]]$importance) <- unlist(mapply(paste, names(model), sapply(model, colnames), MoreArgs = list(sep = '_')))
  }
  return(obj)
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
                               
cov_str_gls <- function(fit, cov = TRUE){
  corm <- corMatrix(fit$modelStruct$corStruct)[[5]]
  varests <- hetro_var(fit=fit)
  covm <- corm*fit$sigma^2*if(!is.null(varests))t(t(varests))%*%t(varests) else 1
  return(covm)
}


#=================================================================================================================================                               
                               
cov_str <- function(fit, cov = TRUE, time_var = "time", hlm = TRUE){
  
  if(inherits(fit, "gls")) return(cov_str_gls(fit=fit, cov=cov))
  
  rho <- rho_lme(fit)
  hetro <- hetro_var(fit)
  sig <- sigma(fit)
  dat <- getData(fit)
  if(!(time_var %in% names(dat))) stop("'time_var' doesn't exist in the data.", call. = FALSE)
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
                                                       
#==================================================================================================================================

multilogit <- function (...){
  
  X <- list(...)
  K <- length(X)
  X <- as.data.frame(X)
  N <- nrow(X)
  if(N == 1){
    f <- exp(X[1, ])
    below <- sum(f)
    as.numeric(f/below)
  } else {
    f <- lapply(1:N, function(i) exp(X[i, ]))
    below <- sapply(1:N, function(i) sum(f[[i]]))
    p <- sapply(1:N, function(i) unlist(f[[i]])/below[i])
    p <- t(as.matrix(p))
    colnames(p) <- NULL
    p
  }
}
                
                
#====================================================================================================================
             
                
inv.multilogit <- function(x, lambda = 1, diff = TRUE, log = FALSE){
  
  x <- round(x)
  
  if(length(x) == 1){ x <- 0:x  ;
  message("Note: ", length(x), " categories were assumed.")
  }
  
  if(diff){ 
    x <- x - min(x)
    f <- exp(lambda * x)
  }
  if(!log){
    output <- f/sum(f)
  } else {
    output <- log(f) - log(sum(f))
  }
  output
}
         
                   
#===========================================================================================================================
                              
                              
logit <- function(x){ 
    
 return(stats::qlogis(x)) 
    
}


#===========================================================================================================================
                              
                              
inv.logit <- function(x, percent = FALSE, digits = 4){
  
  p <- stats::plogis(x)
  return(if(percent) 
  noquote(paste0(round(p*1e2, digits = digits), "%")) else p)
  
}            
         
#===================================================================================================================================
    
    
# Some hack to turn off unnneeded tick mark on the 3rd and 4th axes of plot effects
                                    
plot.efflist <- function (x, selection, rows, cols, graphics = TRUE, 
                          lattice, rug = FALSE, multiline = TRUE, ...) 
{
  lattice <- if (missing(lattice)) 
    list()
  else lattice
  if (!missing(selection)) {
    if (is.character(selection)) 
      selection <- gsub(" ", "", selection)
    pp <- plot(x[[selection]], lattice = lattice, rug = rug, multiline=multiline, ...)
    pp$x.scales$tck=c(1,0)
    pp$y.scales$tck=c(1,0)
    return(pp)
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
      pp <- plot(x[[(i - 1) * cols + j]], lattice = lattice, rug = rug, multiline = multiline,
                 ...)
      # hack to turn off opposite side tick marks
      pp$x.scales$tck=c(1,0)
      pp$y.scales$tck=c(1,0)
      print(pp)
    }
  }
}
environment(plot.efflist) <- asNamespace("effects")
    
#===========================================================================================================================


odiag <- function(x) x[(n <- nrow(x))^2-(1:n)*(n-1)]

#===========================================================================================================================

lo_ave_up <- function(data = NULL, vars, vals = NULL, digits = 9){
  
  if(is.null(vals)){
    sapply(vars, function(x) 
      round(setNames(mean(data[[x]]) + c(-1, 0, 1)*sd(data[[x]]), 
               paste0(x, c('-1SD', '.Mean', '+1SD'))), digits), simplify = FALSE) 
  } else {
    
    setNames(lapply(vars, function(i) vals), vars)
  }
}             

#============================================================================================================================
         
         
inv <- function (X, tol = sqrt(.Machine$double.eps)) 
{
  if (length(dim(X)) > 2L || !(is.numeric(X) || is.complex(X))) 
    stop("'X' must be a numeric or complex matrix")
  if (!is.matrix(X)) 
    X <- as.matrix(X)
  Xsvd <- svd(X)
  if (is.complex(X)) 
    Xsvd$u <- Conj(Xsvd$u)
  Positive <- Xsvd$d > max(tol * Xsvd$d[1L], 0)
  if (all(Positive)) 
    Xsvd$v %*% (1/Xsvd$d * t(Xsvd$u))
  else if (!any(Positive)) 
    array(0, dim(X)[2L:1L])
  else Xsvd$v[, Positive, drop = FALSE] %*% ((1/Xsvd$d[Positive]) * 
       t(Xsvd$u[, Positive, drop = FALSE]))
}   
                    
                    
#========================================================================                        
                    
what_estimate <- function(formula, data){
  
mm <- model.matrix(as.formula(formula), data = data)

noquote(colnames(mm))
}
     
                    
#========================================================================                        
                    
                    
form_length <- function(...){
  
  sapply(list(...),
         function(x) nchar(as.character(x)[[3]]))
  
}                    
 
         
#========================================================================         
         
         
edf <- function(emm_fit) min(as.data.frame(emm_fit)$df, na.rm = TRUE)

#========================================================================
         
sigma_effsize <- function(fit) { 
 
  vc <- VarCorr(fit)
   
  if(inherits(fit, "lme")) sqrt(sum(as.numeric(vc[,"Variance"]), na.rm = TRUE)) else 
    sqrt(sum(as.numeric(c(attr(vc[[1]], "stddev"), attr(vc, "sc")))^2, na.rm = TRUE))
}         

         
#=================================================================================================================================  

lmectl <- function(maxIter = 200, msMaxIter = 1000, niterEM = 500,
                       msMaxEval = 4000)
  {
  lmeControl(maxIter = maxIter, msMaxIter = msMaxIter, niterEM = niterEM,
                         msMaxEval = msMaxEval)
}         
  
         
         
#==== Developing rePCA & isSingular methods for lme models ====================================================================================================


Tlist_lme <- function(fit) rev(pdMatrix(fit$modelStruct$reStruct, factor = TRUE))

theta_lme <- function(fit) sapply(Tlist_lme(fit), function(i) i[lower.tri(i, diag = TRUE)])

lowerbd <- function(x){
     dd <- diag(0, nrow=nrow(x))
     dd[lower.tri(dd)] <- -Inf
     dd[lower.tri(dd, diag=TRUE)]
   }

lwr_lme <- function(fit) sapply(Tlist_lme(fit), lowerbd)
                                  
         
isSingular_lme <- function(fit, tol = 1e-04){ 
  
  lwr <- lwr_lme(fit)
  theta <- theta_lme(fit)
  any(theta[lwr == 0] < tol)
}



rePCA_lme <- function(x){

chfs <- Tlist_lme(x)
nms <- names(chfs)
unms <- unique(nms)
names(unms) <- unms

svals <- function(m) { # this is applied each of the RE matrices
  vv <- svd(m, nv = 0L)
  names(vv) <- c("sdev", "rotation")
  vv$center <- FALSE
  vv$scale <- FALSE
  class(vv) <- "prcomp"
  vv
}

structure(lapply(unms, function(m)
  svals(Matrix::bdiag(chfs[which(nms ==
        m)]))), class = "prcomplist")
}


#====================================================================================================              
         
rm_cor_test_dups <- function(data) {
  data %>%
    filter(var1 != var2) %>%
    mutate(col1 = pmin(var1, var2),
           col2 = pmax(var1, var2)) %>%
    distinct(col1, col2, .keep_all = TRUE) %>%
    select(-col1, -col2)
}         

                 
#======================================================================================================
                 
                 
cor2G2 <- function(sd.int = 6, sd.slope = .01, rho = .3){
    
    cormat <-  matrix(c(sd.int, rho, rho, sd.slope), 2, 2)
    res <- data.frame(lme4::sdcor2cov(cormat), row.names = c("V1", "V2") ) 
    colnames(res) <- rownames(res)
    as.matrix(res)
  }
  
#======================================================================================================                 
                 
cor_DV <- function(rho = .7){  
  
  mu <- rep(0,2)
  Sigma <- cor2G2(sd.int = 2, sd.slope = .275, rho = rho)
  
  rawvars <- as.data.frame(mvrnorm(n=1000, mu=mu, Sigma=Sigma, empirical = TRUE))
  
  p <- ggplot(rawvars)+aes(V1, V2) + geom_point()+stat_ellipse(color=2, size = 2) +
    theme(axis.text.x=element_blank(),
          axis.text.y=element_blank())
  print(p)
  round(cor(rawvars), 6)
}

#cor_DV(.9)                 
   
#========================================================================                        
                 
z2 <- function(n, mu0, y_bar, SIGMA){
  
  ddf = nrow(y_bar)
  
  z_2 = n*t(matrix(y_bar-mu0))%*%solve(SIGMA)%*%matrix(y_bar-mu0)
  
  crit_value = qchisq(.05, df = ddf, lower.tail = FALSE)
  
  p_value = pchisq(z_2, df = ddf, lower.tail = FALSE)
  val <- round(c(z_2, crit_value),2)
  
  ci <- max(qchisq(.99, ddf), val)
  curve(dchisq(x, ddf), -.01, ci, n=1e3, panel.last = abline(v=c(z_2, crit_value), col=1:2, lty = 3:2),
        ylab = "Density", xlab = paste0("X^2 ", "(df=", ddf,")"),lwd = 2)
  pu <- mean(par("usr")[3:4])
  
  
  text(c(z_2, crit_value), pu, c(paste("X^2=",val[1]),paste("Crit.X2=",val[2])), pos=4, srt = 90, xpd = NA)
  data.frame(z2 = z_2, crit.value = crit_value, p.value = p_value)
  }  
  
#========================================================================
                 

m <- function(N, p, k) N - 1 - (p + k)/2
s <- function(N, p, k) sqrt( (p^2*(k-1)^2 -4) / (p^2 + (k-1)^2 -5) )

lambda2F <- function(lambda, N, p, k) 
  { 
  
val =  ((1-sqrt(lambda))/ (sqrt(lambda)))*((m(N,p,k)*s(N,p,k)-(p*(k-1)/2) +1 )/(p*(k-1)))
df1 = p*(k-1)  
df2 = (m(N,p,k)*s(N,p,k)-(p*(k-1))/2)+1

c(F_value=val, df1=df1, df2=df2)
}                 
                 
 
#========================================================================
                 
lambda2x2 <- function(lambda, N, p, k){
  
  X2 = -((N-1)-(.5*(p+k)))*log(lambda)
  df = p*(k-1)
  c(X2=X2, df=df)
}  
                 
                 
#========================================================================                 
                 
trim <- function(X){
  X <- setNames(X, trimws(names(X)))
  y <- sapply(names(X), function(x) is.character(as.vector(X[[x]])))
  X[y] <- lapply(X[y], trimws)
  return(X)
}                 
              
#=========================================================================
              
s_eta <- function(p, k) sqrt( ((p^2)*(k-1)^2 -4 ) / ((p^2)+(k-1)^2 -5 ) )              
                            
              
#=========================================================================              
              
box_m2 <- function(data, group){
    dat <- dplyr::select(data, -{{group}})
    rstatix::box_m(dat, data %>% pull({{group}}))
}
 
              
#========================================================================              
              
F2T2 <- function(F_value, n1, n2, p) F_value / (( n1+n2 - p-1 )/(p*(n1+n2-2)))

F2D2 <- function(F_value, n1, n2, p) ((n1 + n2)/(n1*n2))*F2T2(F_value, n1, n2, p)
           
#========================================================================          
              
#lme2spss <- lmeControl(sigma=1e-5, opt ="optim", returnObject=TRUE)     
              
              
#========================================================================              
              
plot.prof <- function(fit){
  
  if(!inherits(fit, c("lmerMod", "lmerModLmerTest", "lme4", "glmmTMB", "glmerMod"))) stop("Model not supported.", call. = FALSE)
  
  pp <- profile(fit, signames = FALSE)
  
  dd <- as.data.frame(pp)
  
  if(".zeta" %in% names(dd)) names(dd)[which(names(dd) == ".zeta")] <- "value"
  if(inherits(fit, "glmmTMB")) dd$value <- sqrt(dd$value)
  
  ggplot2::ggplot(dd,aes(.focal, value)) +  geom_hline(yintercept = 0, colour = 8, linetype = 2) +
    geom_line(colour = 2) + geom_point(colour = 2) +
    facet_wrap(~.par, scale = "free_x") + xlab("Parameter Value") +
    ylab("Zeta (Normal)")
}
              
#========================================================================   
              
exam.efa <- function(x, factors, data = NULL, covmat = NULL, n.obs = NA,
                     subset = NULL, na.action = "na.omit", start = NULL,
                     scores = c("none", "regression", "Bartlett"),
                     rotation = "varimax", control = NULL, cutoff = .5, digits = 6, plot = TRUE, file.name = NULL, ...){
  
  
  cc <- match.call(expand.dots = FALSE)
  cc[[1]] <- quote(factanal)
  fit <- eval.parent(cc)
  fit$call <- match.call(expand.dots = FALSE)
  
  
  res <- round(transform(subset(as.data.frame.table(fit[[2]]), Freq >= cutoff),
                         Var2 = match(Var2, unique(Var2)), Var1 = as.numeric(Var1)), digits = digits)
  
  names(res) <- c("Item", "Factor", "Loading")
  
  rownames(res) <- NULL
  
  file.name <- trimws(file.name)  
  
  if(length(file.name) != 0){
    
    nm <- paste0(file.name, ".csv")
    ur <- try(write.csv(res, nm, row.names = FALSE), silent = TRUE)
    if(inherits(ur, "try-error")) stop(paste0("\nClose the Excel file '", nm, "' and try again OR pick another file name."), call. = FALSE)
    message(paste0("\nNote: Check folder '", basename(getwd()),"' for the Excel file '", nm, "'.\n"))
  } 
  
  f <- res$Factor
  i <- res$Item
  u <- unique(f)
    
  if(plot){
    
    graphics.off()
    org.par <- par(no.readonly = TRUE)
    on.exit(par(org.par))
    par(mar = c(3.8, 1.1, 1.5, 4.1), mgp = c(1.7, .5, 0))
    
    y <- unlist(lapply(u, function(j) seq_len(sum(f == j))))
    
    plot(f, y, las = 1, pch = 22, cex = 1.2, xlim = c(-.1, max(f)+.1), axes = FALSE, xlab = NA, main = NA, font.lab = 2, ylab = NA)
    
    at <- mean(u)
    
    mtext("FACTORS", 1, line = 2.5, font = 2, at = at)
    
    text(f, 0, f, pos = 1, xpd = NA, font = 2, cex = 1.75)
    
    rect(u-.5, 0, u+.5, tapply(y, f, FUN = max)+.5, col = adjustcolor(1:8, .25), xpd = NA, border = NA)
    
    dup <- duplicated(i) | duplicated(i, fromLast = TRUE)
    
    text(f, y, i, pos = 4, cex = .7, xpd = NA, font = 2, col = ifelse(dup, 2, 1))
    
    legend(at, par('usr')[4], legend = "ITEMS", pch = 22, horiz = TRUE, bty = "n", text.font = 2, xpd = NA, pt.cex = 1.4, yjust = .5)
  }
  return(res)
}                                
              
#========================================================================              

# 'sjPlot', 'sjstats'    
    
needzzsf <- c('car','psych','tidyverse','lme4','nlme','MASS','matrixcalc', 'haven', 'lmerTest', 'reghelper',
          'parallel','rela','bbmle', 'jtools','interactions','GPArotation', 'rstatix',
          'effects', 'effectsize', "paran", 'Hotelling','ICC')

                        
not.have23 <- needzzsf[!(needzzsf %in% installed.packages()[,"Package"])]
if(length(not.have23)) install.packages(not.have23)


suppressWarnings(
suppressMessages({ 
  
  for(i in needzzsf){
    library(i, character.only = TRUE)
  }
}))
  
options(dplyr.summarise.inform = FALSE)
              
source('https://raw.githubusercontent.com/rnorouzian/A/main/a.r') 
         
         
formals(lmerControl)$check.nobs.vs.nRE <- "ignore"       
formals(glmerControl)$check.nobs.vs.nRE <- "ignore"
formals(glmerControl)$check.nobs.vs.nlev <- "ignore"
formals(lmerControl)$check.nobs.vs.nlev <- "ignore"           
