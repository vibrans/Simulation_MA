library(SimDesign)
setwd("E:/knobivan/Documents/MA/R/Simulation/parallel/Simulation_MA_4/Daten/")
############################ data generation ##################################################
R <- 1000
Design <- expand.grid(
  N = c(100, 500, 1000),
  
  # reliabilities: 0.9, 0.8, 0.75, 0.625, 0.5, 1/3 (if var(xi) = 1 and lambda11 = 1)
  errvar = c(1/9, 0.25, 1/3, 0.6, 1, 2),
  
  # effect sizes when sd(Y|X=0) = sqrt(1.43): 
  # 0, 0.25, 0.50, 0.8, 1
  effect = c(0, 0.3, 0.6, 0.95, 1.2)      
)


# Design <- expand.grid(
#   N = c(100),
#   
#   # reliabilities: 0.9, 0.8, 0.75, 0.625, 0.5, 1/3 (if var(xi) = 1)
#   errvar = c(0.6, 1, 2),
#   
#   # effect sizes when sd(Y|X=0) = sqrt(1.43): 
#   # 0, 0.25, 0.50, 0.8, 1
#   effect = c(0.6, 0.95)      
# )

Generate <- function(condition, fixed_objects=NULL){
  # true score and indicator variables
  data = MASS::mvrnorm(n=condition$N, mu=c(0, 0), Sigma=matrix(c(1, 0.5, 0.5, 1), nrow=2), empirical=TRUE)
  xi1 = data[, 1]
  xi2 = data[, 2]
  
  y11 <- 0 + 1*xi1 + rnorm(condition$N, 0, sqrt(condition$errvar))
  y21 <- 0.5 + 0.9*xi1 + rnorm(condition$N, 0, sqrt(condition$errvar))
  y31 <- 0.4 + 0.7*xi1 + rnorm(condition$N, 0, sqrt(condition$errvar))
  
  y12 <- 0 + 1*xi2 + rnorm(condition$N, 0, sqrt(condition$errvar))
  y22 <- 0.2 + 0.7*xi2 + rnorm(condition$N, 0, sqrt(condition$errvar))
  y32 <- 0.3 + 0.8*xi2 + rnorm(condition$N, 0, sqrt(condition$errvar))
  
  # assignment to treatment
  x <- rbinom(condition$N, 1, pnorm(0.4*xi1 + 0.6*xi2))
  
  # outcome
  y <- 0.5*xi1 + 0.7*xi2 + condition$effect*x + rnorm(condition$N, 0, 0.8)
  
  id <- 1:condition$N
  dat <- data.frame(id, y, x, y11, y21, y31, y12, y22, y32)
}

####################################################################################
Analyse <- function(condition, dat, fixed_objects=NULL){
  # JM
  r1 <- computeRaykovJM(dat, pred="regression")
  names(r1) <- c("JMreg_Xconv", "glmJMreg_Xconv", 
                 "m1JMreg_Xconv", "m1JMreg_ave", "m1JMreg_se",
                 "m2JMreg_Xconv", "m2JMreg_ave", "m2JMreg_se")
  r2 <- computeLatentPSonegroup(dat)
  r3 <- computeLatentPS(dat)
  r4 <- computeEffectLiteR(dat)
  sd0 <- sd(dat[dat$x==0, "y"])
  c(r1, r2, r3, r4, sd0=sd0, glass=condition$effect/sd0)
}

# input: Xconv ... names of convergence-indicator-columns for one estimationmodel
# task: in later computations rows with nonconverged models should be ignored
# as well as rows with implausible values
# output: vector of rownumbers of innocent trials (all models converged)
# bug: if none of the trials are nonconverged return value is numeric(0) causes further problems
rowsXconv <- function(convvars) {
  if (length(convvars)==1){
    rX <- as.numeric(rownames(results[!(results[, convvars]==1), , drop=F]), do.NULL=F)
  }else{
    ex <- quote(results[, convvars[1]]==1)
    for (i in 2:length(convvars)){
      ex <- rlang::expr(!!ex | results[, convvars[!!eval(i)]]==1)
    }
    rX <- as.numeric(rownames(results[!(eval(ex)), , drop=F], do.NULL=F))
  }
  rX
}


Summarise <- function(condition, results, fixed_objects=NULL){
  ## character vector for names of average effects
  # must be in correct order!
  labels_ave <- c("m1JMreg_ave", "m2JMreg_ave",
                  "m3_ave", "m4_ave", "m5_ave")
  
  labels_se <- c("m1JMreg_se", "m2JMreg_se",
                 "m3_se", "m4_se", "m5_se")
  
  
  ## set rownames of resultsmatrix 
  rownames(results) <- 1:R
  
  ## set environment of external functions to local environment
  # so that they can find resultsmatrix
  environment(rowsXconv) <- environment()
  
  ## number of nonconverged models
  SX <- sapply(X=c(
    # JM
    "JMreg_Xconv", "glmJMreg_Xconv", "m1JMreg_Xconv", "m2JMreg_Xconv",
    
    # the rest
    "m3_Xconv",
    "mProbit_Xconv", "m4_Xconv",
    "m5_Xconv"),
    FUN=function(x) {
      S <- sum(results[, x])
      names(S) <- paste0("S_", x)
      return(S)
    },
    USE.NAMES=F)
  
  ## find rownumbers of converged trials
  # conv ... list of vectors with rownumbers of all converged trials modelspecific
  conv <- sapply(X=list(c("JMreg_Xconv", "glmJMreg_Xconv", "m1JMreg_Xconv"),
                        c("JMreg_Xconv", "glmJMreg_Xconv", "m2JMreg_Xconv"),
                        #---
                        c("m3_Xconv"),
                        c("mProbit_Xconv", "m4_Xconv"),
                        c("m5_Xconv")),
                 FUN=function(x) rowsXconv(x))
  if (!is.list(conv)){
    conv <- lapply(seq_len(ncol(conv)), function(i) conv[,i])  # make list if not already
  }
  
  ## find rownumbers of converged trials and nonoutlier trials modelspecific
  rno <- mapply(FUN=function(x, y){
    wh <- robustbase::adjboxStats(results[x, y], coef=12)
    u <- as.numeric(rownames(results[results[, y] <= wh$fence[[1]], , drop=F]))
    o <- as.numeric(rownames(results[results[, y] >= wh$fence[[2]], , drop=F]))
    d <- setdiff(x, c(u, o))
    return(d)
  },
  x=conv, y=labels_ave)
  if (!is.list(rno)){
    rno <- lapply(seq_len(ncol(rno)), function(i) rno[,i])  # make list if not already
  }
  
  
  # number of outliers (with medcouple) without considering the nonconverged trials
  SO <- mapply(FUN=function(x, y, z){
    so <- length(x) - length(y)
    names(so) <- paste0("SO_", gsub("_ave", "", z, perl=TRUE))
    return(so)
  },
  x=conv,
  y=rno,
  z=labels_ave
  )
  
  ## means of estimated average effects
  M <- mapply(FUN=function(x, y){
    m <- mean(results[y, x])
    names(m) <- paste0("M_", gsub("_ave", "", x, perl=TRUE))
    return(m)
  },
  x=labels_ave, y=rno,
  USE.NAMES=FALSE)
  
  ## SE with corrected N
  SE <- mapply(FUN=function(x, y){
    se <- sd(results[x, y])
    names(se) <- paste0("SE_", gsub("_ave", "", y, perl=TRUE))
    return(se)
  },
  x=rno, y=labels_ave,
  USE.NAMES=FALSE)
  
  ## biases
  B <- sapply(X=labels_ave,
              FUN=function(x) {
                z <- paste0("M_", gsub("_ave", "", x, perl=TRUE))
                b <- M[[z]] - condition$effect
                names(b) <- paste0("B_", gsub("_ave", "", x, perl=TRUE))
                return(b)
              }, USE.NAMES=F)
  
  ## relative biases
  RB <- mapply(FUN=function(x, y){
    rb <- x/condition$effect
    names(rb) <- paste0("RB_", gsub("B_", "", y, perl=TRUE))
    return(rb)
  },
  x=B, y=names(B),
  USE.NAMES=FALSE)
  
  ## mean square errors of average effects
  MSE <- mapply(FUN=function(x, y){
    mse <- sum((results[y, x] - condition$effect)^2)/(length(results[y, x]))
    names(mse) <- paste0("MSE_", gsub("_ave", "", x, perl=TRUE))
    return(mse)
  }, x=labels_ave, y=rno,
  USE.NAMES=FALSE)
  
  ## relative efficiency
  # reference estimator is the effectLiteR estimator
  RE <- sapply(X=c("MSE_m1JMreg", "MSE_m2JMreg",
                   "MSE_m3", "MSE_m4", "MSE_m5"),
               FUN=function(x) {
                 re <- MSE[[x]]/MSE[["MSE_m1JMreg"]]
                 names(re) <- paste0("RE_", gsub("MSE_", "", x, perl=TRUE))
                 return(re)
               }, USE.NAMES=F)
  
  # empirical type-I error rate and False Coverage Rate with empirical SE
  # when nulleffect variable takes on empirical type-I-error rate
  # when there actually is an effect variable takes on the value of the Power
  t <- qt(0.975, condition$N-1)
  AlPow <- mapply(
    FUN = function(x, y, z) {
      I <- abs((results[y, x]) / results[y , z]) > t
      alpow <- sum(I)/length(y)
      names(alpow) <- paste0("AlPow_", gsub("_ave", "", x, perl=TRUE))
      return(alpow)
    },
    x = labels_ave,
    y = rno,
    z = labels_se,
    USE.NAMES = FALSE)
  
  
  ## coverage with 95% CI = 1-FCR
  COV <- mapply(
    FUN = function(x, y, z) {
      ci <- t*results[y, z]
      I <- results[y, x] - ci < condition$effect &
        condition$effect < results[y, x] + ci
      cov <- sum(I)/length(y)
      names(cov) <- paste0("COV_", gsub("_ave", "", x, perl=TRUE))
      return(cov)
    },
    x = labels_ave,
    y = rno,
    z = labels_se,
    USE.NAMES = FALSE)
  
  
  ## bias of standard error
  BSE <- mapply(FUN=function(x, y, z){
    m <- mean(results[y, x])
    bse <- m - z
    names(bse) <- paste0("BSE_", gsub("_se", "", x, perl=TRUE))
    return(bse)
  },
  x=labels_se, y=rno, z=SE,
  USE.NAMES=FALSE)
  
  
  ## results
  ret <- c(SX, SO, M, SE, B, RB, MSE, RE, AlPow, COV, BSE)
  return(ret)
}


# mfs ... factor score model
# m1 ... Propensity Score Modell mit Factor Scores
# m2 ... doubly robust Propensity Score Modell mit Factor Scores
# m1SO; m1SM; m1JO; m1JM
# m3 ... latentes Propensity Score Einschritt-Eingruppen-Modell
# mProbit ... sem for computing coefficients of probit regression
# m4 ... latentes Propensity Score Mehrschritt-Mehrgruppen-Modell
# m5 ... EffectLiteR

results <- runSimulation(design=Design, replications=R, 
                         generate=Generate, analyse=Analyse, summarise=Summarise, 
                         parallel=T, verbose=T, 
                         filename="Simulation_MA_4",
                         packages=c("MASS", "lavaan", "EffectLiteR", "robustbase"),
                         save_seeds=T, save_results=T, save=T)

####################################### Raykov ###########################################

##---------------- Raykov joint multigroup models -----------------------------##

computeRaykovJM <- function(d, pred){
  ### catching errors and warnings
  # flips
  fs_Xconv <- 0
  glm_Xconv <- 0
  m1_Xconv <- 0
  m2_Xconv <- 0
  
  # handlers
  h_fs <- function(w){
    if(any( grepl( "model has NOT converged", w) )) { # NONconvergence of fs model
      fs_Xconv <<- 1
    }
    invokeRestart( "muffleWarning" )
  }
  
  h_glm <- function(w){
    if(any( grepl( "did not converge", w) )) {
      glm_Xconv <<- 1
    }
    invokeRestart("muffleWarning")
  }
  
  #---------
  h_m1 <- function(w){
    if(any( grepl( "model has NOT converged", w) )) {
      m1_Xconv <<- 1
    }
    invokeRestart("muffleWarning")
  }
  
  #---------
  h_m2 <- function(w){
    if(any( grepl( "model has NOT converged", w) )) {
      m2_Xconv <<- 1
    }
    invokeRestart("muffleWarning")
  }
  
  
  #---------
  
  # step 1: estimate factor scores with joint multigroup model
  mfs <- '
  xi1 =~ c(1,1)*y11 + c(la21,la21)*y21 + c(la31,la31)*y31
  xi2 =~ c(1,1)*y12 + c(la22,la22)*y22 + c(la32,la32)*y32
  
  y11 ~ c(nu11,nu11)*1
  y21 ~ c(nu21,nu21)*1
  y31 ~ c(nu31,nu31)*1
  xi1 ~ c(m0xi1, m1xi1)*1 + c(0, NA)*1
  
  y11 ~~ y11
  y21 ~~ y21
  y31 ~~ y31
  xi1 ~~ c(NA,NA)*xi1
  
  y12 ~ c(nu12,nu12)*1
  y22 ~ c(nu22,nu22)*1
  y32 ~ c(nu32,nu32)*1
  xi2 ~ c(m0xi2, m1xi2)*1 + c(0, NA)*1
  
  y12 ~~ y12
  y22 ~~ y22
  y32 ~~ y32
  xi2 ~~ c(NA,NA)*xi2
  
  xi1 ~~ c(NA,NA)*xi2
  '
  
  fit_fs <- withCallingHandlers(sem(mfs, data=d, group="x", group.label=c("0","1"),
                                    group.equal=c("intercepts","loadings")), warning=h_fs)
  
  # compute factor scores and merge with df
  fs <- withCallingHandlers(data.frame(do.call("rbind", lavPredict(fit_fs, method=pred))),
                            warning=function(w) invokeRestart("muffleWarning"))
  fs$id <- unlist(lavInspect(fit_fs, "case.idx"))
  colnames(fs) <- c("xi1_fs", "xi2_fs", "id")
  d <- merge(d, fs)
  
  
  # step 3: estimate probit coefficients
  fit_MPS <- withCallingHandlers(glm(x ~ xi1_fs + xi2_fs, family=binomial(link='probit'),
                                     data=d), warning=h_glm)
  # step 4: compute MPS
  d$MPS <- pnorm(coef(fit_MPS)[1] + coef(fit_MPS)[2]*d$xi1_fs + coef(fit_MPS)[3]*d$xi2_fs)
  
  # ----------------- #
  
  #### modified Raykov #####
  # modified procedure: for better comparability regarding our proposal
  # - probit coefficients estimated with lavaan
  # - coefficients of probit used because data simulated with probit
  # - probit transformed MPS as only covariate in regression
  d$probitMPS <- coef(fit_MPS)[1] + coef(fit_MPS)[2]*d$xi1_fs + coef(fit_MPS)[3]*d$xi2_fs
  
  m1 <- '
  y ~ c(a01,a11)*probitMPS
  y ~ c(a00,a10)*1
  probitMPS ~ c(probitMPS_0,probitMPS_1)*1
  
  group % c(gw0,gw1)*w
  N := exp(gw0) + exp(gw1)
  relfreq0 := exp(gw0)/N
  relfreq1 := exp(gw1)/N
  
  probitMPS := probitMPS_0*relfreq0 + probitMPS_1*relfreq1
  
  g10 := a10 - a00
  g11 := a11 - a01
  ave := g10 + g11*probitMPS'
  
  fit_m1 <- withCallingHandlers(sem(model=m1, data=d, group="x", group.label=c("0","1")),
                                warning=h_m1)
  
  m1_ave <- parameterEstimates(fit_m1)[parameterEstimates(fit_m1)$lhs=="ave", "est"]
  m1_se <- parameterEstimates(fit_m1)[parameterEstimates(fit_m1)$lhs=="ave", "se"]
  
  # ----------------- #
  #### original Raykov #####
  # (Raykov computes his regression with MPS,
  # and estimated True Score Variables and manifest nonfallible variables)
  # use lavaan instead of lm, because of lost se of ave when ave computed from
  # lm-with-interaction estimated coefficients
  # also (regarding comparability) smaller SE than in lavaan
  m2 <- '
  y ~ c(a01,a11)*MPS + c(a02,a12)*xi1_fs + c(a03,a13)*xi2_fs
  y ~ c(a00,a10)*1
  MPS ~ c(mMPS_0,mMPS_1)*1
  xi1_fs ~ c(m0xi1_fs,m1xi1_fs)*1
  xi2_fs ~ c(m0xi2_fs,m1xi2_fs)*1
  
  group % c(gw0,gw1)*w
  N := exp(gw0) + exp(gw1)
  relfreq0 := exp(gw0)/N
  relfreq1 := exp(gw1)/N
  
  mMPS := mMPS_0*relfreq0 + mMPS_1*relfreq1
  mxi1_fs := m0xi1_fs*relfreq0 + m1xi1_fs*relfreq1
  mxi2_fs := m0xi2_fs*relfreq0 + m1xi2_fs*relfreq1
  
  g10 := a10 - a00
  g11 := a11 - a01
  g12 := a12 - a02
  g13 := a13 - a03
  ave := g10 + g11*mMPS + g12*mxi1_fs + g13*mxi2_fs'
  
  fit_m2 <- withCallingHandlers(sem(model=m2, data=d, group="x", group.label=c("0","1")), 
                                warning=h_m2)
  
  m2_ave <- parameterEstimates(fit_m2)[parameterEstimates(fit_m2)$lhs=="ave", "est"]
  m2_se <- parameterEstimates(fit_m2)[parameterEstimates(fit_m2)$lhs=="ave", "se"]
  
  # -------- #
  
  return(c(fs_Xconv,
           glm_Xconv,
           m1_Xconv, m1_ave, m1_se,
           m2_Xconv, m2_ave, m2_se))
}

#################### latent PS approach as onegroup-onestep model #################################
# reason: no treatment x covariate-IA in simulated data
# one step until now only possible if no multigroup model is estimated
# BUT: onegroup models have problems with dependent parameters; let's have a look
computeLatentPSonegroup <- function(d){
  
  ### catching errors and warnings
  # flips
  m3_Xconv <- 0
  
  # handlers
  h_m3 <- function(w){
    if(any( grepl( "model has NOT converged", w) )) { # NONconvergence of fs model
      m3_Xconv <<- 1
    }
    invokeRestart( "muffleWarning" )
  }
  
  m3 <- '
  # für die PSs
  xi1 =~ 1*y11 + y21 + y31
  xi2 =~ 1*y12 + y22 + y32
  
  y11 ~ 0*1
  xi1 ~ NA*1
  
  y12 ~ 0*1
  xi2 ~ NA*1
  x ~ c1*xi1 + c2*xi2
  #x ~~ NA*xi1   # besserer Modelfit als ohne diese Zeile, aber anderer Schätzer für x = AVE
  #x ~~ NA*xi2
  x ~ c0*1
  
  # effect regression
  probit <~ c1*xi1 + c2*xi2
  probit ~ c0*1
  
  xi1 ~~ 0*probit
  xi2 ~~ 0*probit
  y ~ ave*x + probit
  '
  fit_m3 <- withCallingHandlers(sem(m3, data=d), warning=h_m3)
  
  m3_ave <- parameterEstimates(fit_m3)[parameterEstimates(fit_m3)$label=="ave", "est"]
  m3_se <- parameterEstimates(fit_m3)[parameterEstimates(fit_m3)$label=="ave", "se"]
  
  return(c(m3_Xconv=m3_Xconv, m3_ave=m3_ave, m3_se=m3_se))
}


################################ latent PS approach ##########################################
computeLatentPS <- function(d){
  ### catching errors and warnings
  # flips
  mProbit_Xconv <- 0
  m4_Xconv <- 0
  
  # handlers
  h_mProbit <- function(w){
    if(any( grepl( "model has NOT converged", w) )) { # NONconvergence of fs model
      mProbit_Xconv <<- 1
    }
    invokeRestart( "muffleWarning" )
  }
  
  h_m4 <- function(w){
    if(any( grepl( "model has NOT converged", w) )) { # NONconvergence of fs model
      m4_Xconv <<- 1
    }
    invokeRestart( "muffleWarning" )
  }
  
  mProbit <- '
  xi1 =~ 1*y11 + y21 + y31
  xi2 =~ 1*y12 + y22 + y32
  
  y11 ~ 0*1
  xi1 ~ NA*1
  
  y21 ~ 0*1
  xi2 ~ NA*1
  
  x ~ xi1 + xi2
  
  #x ~~ NA*xi1   ##########
  #x ~~ NA*xi2 
  '
  fit_Probit <- withCallingHandlers(sem(mProbit, data=d, ordered="x", parameterization="theta"),
                                    warning=h_mProbit)
  
  m4 <- paste0('
               xi1 =~ c(1,1)*y11 + c(la21,la21)*y21 + c(la31,la31)*y31
               xi1 ~ c(mxi10, mxi11)*1
               y11 ~ c(0,0)*1
               y21 ~ c(nu21,nu21)*1
               y31 ~ c(nu31,nu31)*1
               
               xi2 =~ c(1,1)*y12 + c(la22,la22)*y22 + c(la32,la32)*y32
               xi2 ~ c(mxi20, mxi21)*1
               y12 ~ c(0,0)*1
               y22 ~ c(nu22,nu22)*1
               y32 ~ c(nu32,nu32)*1
               
               probit <~ ', coef(fit_Probit)["x~xi1"], '*xi1 + ', coef(fit_Probit)["x~xi2"], '*xi2
               probit ~ ', coef(fit_Probit)["x|t1"]*(-1), '*1
               xi1 ~~ 0*probit
               xi2 ~~ 0*probit
               
               y ~ c(a01,a11)*probit ## with interaction
               y ~ c(a00,a10)*1
               
               
               group % c(gw0,gw1)*w
               N := exp(gw0) + exp(gw1)
               relfreq0 := exp(gw0)/N
               relfreq1 := exp(gw1)/N
               
               mprobit0 := ', coef(fit_Probit)["x|t1"]*(-1), ' + ', coef(fit_Probit)["x~xi1"], '*mxi10 + ', coef(fit_Probit)["x~xi2"], '*mxi20
               mprobit1 := ', coef(fit_Probit)["x|t1"]*(-1), ' + ', coef(fit_Probit)["x~xi1"], '*mxi11 + ', coef(fit_Probit)["x~xi2"], '*mxi21
               mprobit := mprobit0*relfreq0 + mprobit1*relfreq1
               
               g10 := a10 - a00
               g11 := a11 - a01
               
               ave := g10 + g11*mprobit  ## average effect
               ')
  
  fit_m4 <- withCallingHandlers(sem(m4, data=d, group="x", group.label=c("0","1")),
                                warning=h_m4)
  
  # results
  m4_ave <- parameterEstimates(fit_m4)[parameterEstimates(fit_m4)$lhs=="ave", "est"]
  m4_se <- parameterEstimates(fit_m4)[parameterEstimates(fit_m4)$lhs=="ave", "se"]
  
  return(c(mProbit_Xconv=mProbit_Xconv, 
           m4_Xconv=m4_Xconv, m4_ave=m4_ave, m4_se=m4_se))
}


############################## EffectLiteR approach ##########################################

computeEffectLiteR <- function(d){
  ### catching errors and warnings
  # flip
  m5_Xconv <- 0
  
  # handler
  h <- function(w){
    if(any( grepl( "model has NOT converged", w) )) {
      m5_Xconv <<- 1
    }
    # w$message <- strtrim(w$message, 20)
    # warning(w$message)
    invokeRestart("muffleWarning")
  }
  
  mm <- '
  xi1 =~ c(1,1)*y11 + c(la21,la21)*y21 + c(la31,la31)*y31
  xi2 =~ c(1,1)*y12 + c(la22,la22)*y22 + c(la32,la32)*y32
  
  y11 ~ c(0,0)*1
  y21 ~ c(nu21,nu21)*1
  y31 ~ c(nu31,nu31)*1
  
  y12 ~ c(0,0)*1
  y22 ~ c(nu22,nu22)*1
  y32 ~ c(nu32,nu32)*1
  '
  
  fit_m5 <- withCallingHandlers(effectLite(y="y", x="x", z=c("xi1","xi2"), measurement=mm, data=d),
                                warning=h)
  
  m5_ave <- fit_m5@results@Egx[1, "Estimate"]
  m5_se <- fit_m5@results@Egx[1, "SE"]
  
  return(c(m5_Xconv=m5_Xconv, m5_ave=m5_ave, m5_se=m5_se))
}


