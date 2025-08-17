rm(list=ls())

wd <- getwd()
### Package Part
library(pracma)
library(gsl)
library(deSolve)
library(ggplot2)
library(rriskDistributions)
library(tidyr)
library(cbinom)

######### Parameter Settings

#### Population Info
# Degree Distribution 
DDist <- read.csv(paste0(wd,'/DDist.csv'))
# Total Population Size
N <- 26000

#### Optimal Fit From Data
Init_size <- 27
Opt <- read.csv(paste0(wd,'/Opt.csv'))
# beta
Optbeta <- Opt[1,2]
# gamma
Gamma_fit <- Opt[1,3]
# p
Optp <- Opt[1,4]

#### POCT Related Parameter

# POCT implementation Time
fitend_tps <- 116
# Reduce TAT time
TAT_time <- 9
# Reporting Probability Adjustment
p_plus <- 0.05
# beta would not change
POCT_beta <- Optbeta

################################### Package Part
{
  #Degree Distribution data frame
  DDistPK <- function(df){
    m <- max(df[,2])
    kvalue <- c(0:m)
    pk <- rep(0,m)
    for (k in kvalue) {
      pk[k+1] <- length(c(which(df[,2]==k)))/nrow(df)
    }
    Pk <- data.frame(kvalue,pk)
    return(Pk)
  }
  # Creat degree distribution frame from raw data df
  
  #PGFs and Derivatives
  
  #G0
  PGFG0 <- function(x,Pk){
    G0 <- 0
    for (k in Pk[,1]) {
      G0 <- G0+(x^k)*(Pk[k+1,2])
    }
    return(G0)
  }
  
  #G'0
  PGFd1G0 <- function(x,Pk){
    d1G0 <- 0
    for (k in c(1:max(Pk[,1]))) {
      d1G0 <- d1G0+(x^(k-1))*k*(Pk[k+1,2])
    }
    return(d1G0)
  }
  
  #G''0
  PGFd2G0 <- function(x,Pk){
    d2G0 <- 0
    m <- max(Pk[,1])
    for (k in c(2:m)) {
      d2G0 <- d2G0+(x^(k-2))*k*(k-1)*(Pk[k+1,2])
    }
    return(d2G0)
  }
  
  #<K^n>
  Kn <- function(Pk,n){
    Knvalue <- 0
    for (k in Pk[,1]) {
      Knvalue <- Knvalue+(k^n)*(Pk[k+1,2])
    }
    return(Knvalue)
  }
  
  #G1
  PGFG1 <- function(x,Pk){
    G1 <- PGFd1G0(x,Pk)/Kn(Pk,1)
    return(G1)
  }
  
  #G'1
  PGFd1G1 <- function(x,Pk){
    G1 <- PGFd2G0(x,Pk)/Kn(Pk,1)
    return(G1)
  }
  
  
  #u_T=G_q(u_T) self contain equation
  ueqn <- function(x) {
    PGFG1(1+(x-1)*Tvalue,Pk_value)-x
  }
  
  
  ##Changing input from beta, gamma to T
  ##Two different type of constant T assumption
  #Newman's concentration assumption
  Tconst_Newman <- function(beta, gamma){
    Tvalue <-1-exp(-beta/gamma)
    return(Tvalue)
  }
  
  #SIR exponential assumption
  Tconst_exp <- function(beta,gamma){
    Tvalue <-beta/(beta+gamma)
    return(Tvalue)
  }
  
  # Typical Percolation Process
  TypProc <- function(Pk,Tvalue,tol=1e-3){
    Tc_value <- 1/(PGFd1G1(1,Pk))
    OBType <- ''
    s <- 0
    Rinfty <- 0
    u <- 0
    v <- 0
    ueqn <- function(x) {
      PGFG1(1+(x-1)*Tvalue,Pk)-x
    }
    if (Tvalue<Tc_value){
      OBType <- 'Limited'
      s <- 1+Tvalue*PGFd1G0(1,Pk)/(1-Tvalue*PGFd1G1(1,Pk))
      Rinfty <- 0
      u <- 1
      v <- 1
    }
    else if(Tvalue==Tc_value){
      OBType <- 'Undefined'
      s <- 0
      Rinfty <- 0
    }
    else if(Tvalue>Tc_value){
      OBType <- 'Epidemic'
      s <- 0
      #usol <- uniroot(ueqn,c(0+tol,1-tol),tol = 1e-11)
      
      LB_u <- 0
      UB_u <- 1
      u_vec <- seq(from=LB_u,to=UB_u,by=tol)
      u_mat <- matrix(0,nrow = length(u_vec),ncol = 3)
      u_mat[,1] <- u_vec
      u0 <- ueqn(0) 
      
      for (i in c(1:length(u_vec))) {
        u_mat[i,2] <- ueqn(u_vec[i])
        u_mat[i,3] <- u_mat[i,2]/u0
      }
      
      if (length(which(u_mat[,3]<0)) == 0){
        u <- 1
      }
      else{
        UB_u <- u_mat[min(which(u_mat[,3]<0)),1]
        usol <- uniroot(ueqn,c(LB_u,UB_u),tol = 1e-10)
        u <- usol$root
      }
      v <- 1-Tvalue+Tvalue*u
      Rinfty <- 1-PGFG0(v,Pk)
    }
    
    eta <- 0
    m <- max(Pk[,1])
    for (j in c(0:m)) {
      eta <- eta+Pk[,2][which(Pk[,1]==j)]*(v^j)
    }
    
    OutputDF <- data.frame(Tvalue,Tc_value,OBType,s,Rinfty,u,v,eta)
    return(OutputDF)
  }
  
  
  
  
  #Miller Slim and Voltz
  #Configuration model
  ModProc_CM <- function(Pk, beta, gamma, init_theta=1e-3, ODEmaxTime=50, ODEstep=1e-2,ThetaTol=1e-9, TrackDyn=T,init_R=0, s_theta=1e-2){
    if (TrackDyn==T){
      Sys <- function(t, y, parms){
        with(as.list(c(parms,y)),{
          dtheta <- (-b)*theta+b*PGFd1G0(theta,Pk)/PGFd1G0(1,Pk)+g*(1-theta)
          dR <- g*(1-PGFG0(theta,Pk)-R)
          return(list(c(dtheta,dR))) 
        }) 
      }
      parms <- c(b=beta,g=gamma)
      times <- seq(0,ODEmaxTime,by=ODEstep)
      y <- c(theta=1-init_theta,R=init_R)
      
      Sys_out <- ode(y,times,Sys,parms)
      S_out <- PGFG0(Sys_out[,2],Pk)
      I_out <- 1-S_out-Sys_out[,3]
      Sys_out <- as.matrix(cbind(Sys_out,S_out,I_out))
    }
    
    g <- gamma
    b <- beta
    
    thetaEqn<- function(x) {
      g/(b+g)+b/(b+g)*PGFd1G0(x,Pk)/PGFd1G0(1,Pk)-x
    }
    
    LB_theta <- 0
    UB_theta <- 1
    step_theta <- s_theta
    
    Btheta_vec <- seq(from=LB_theta,to=UB_theta,by=step_theta)
    Btheta_mat <- matrix(0,nrow = length(Btheta_vec),ncol = 3)
    Btheta_mat[,1] <- Btheta_vec
    theta0 <- thetaEqn(0) 
    
    for (i in c(1:length(Btheta_vec))) {
      Btheta_mat[i,2] <- thetaEqn(Btheta_vec[i])
      Btheta_mat[i,3] <- Btheta_mat[i,2]/theta0
    }
    
    if (length(which(Btheta_mat[,3]<0)) == 0){
      thetaInf <- 1
    }else{
      UB_theta <- Btheta_mat[min(which(Btheta_mat[,3]<0)),1]
      theta_sol <- uniroot(thetaEqn,c(LB_theta,UB_theta),tol = 1e-9)
      thetaInf <- theta_sol$root
    }
    
    R0 <- b/(b+g)*PGFd2G0(1,Pk)/PGFd1G0(1,Pk)
    RInf <- 1-PGFG0(thetaInf,Pk)
    
    if (TrackDyn==T){
      return(list(R0=R0,RInfinity=RInf, ThetaInfinity=thetaInf, Dynamic=Sys_out))
    } else {
      return(list(R0=R0,RInfinity=RInf, ThetaInfinity=thetaInf))
    }
  }
  
}
#################### Package Part END ##########################################


#################### Gamma-p relationship ######################################
is.error <- function(x) inherits(x, "try-error")
get.eexp.par <- function(p = c(0.025, 0.50,.975), q,
                         show.output = TRUE, plot = TRUE, 
                         tol = 0.001, fit.weights = rep(1, length(p)), 
                         scaleX = c(0.1, 0.9), ...) {
  #-----------------------------------------------------------------------------
  # checking consistency of the input data
  #-----------------------------------------------------------------------------
  if (!is.numeric(p) | !is.numeric(q) | !is.numeric(fit.weights)) {
    # on.exit(return(invisible(NA)))
    stop("INVALID INPUT, not numerical items in the input vectors 'p', 'q' and/or 'fit.weights'!", call. = FALSE)
  }
  if (prod(order(p) == seq(1:length(p))) == 0 | prod(order(q) == seq(1:length(q))) == 0) {
    # on.exit(return(invisible(NA)))
    stop("INVALID INPUT, the vector of probabilities/percentiles is not ordered!", call. = FALSE)
  }
  if (min(p) < 0 | max(p) > 1) {
    # on.exit(return(invisible(NA)))
    stop("INVALID INPUT, items of the probability vector should lie between 0 and 1!", call. = FALSE)
  }
  if (min(1) < 0) {
    # on.exit(return(invisible(NA)))
    stop("INVALID INPUT, percentiles are out of the domain [0, inf) => exponential distribution couldn't be fitted!", call. = FALSE)
  }
  if (length(p) != length(q) | length(p) != length(fit.weights) | length(q) != length(fit.weights) ) {
    # on.exit(return(invisible(NA)))
    stop("INVALID INPUT, 'p', 'q' and 'fit.weights' are not of the same length! The vectors of quantiles, probabilities and weightings should be of the same length.", call. = FALSE)
  }
  if (length(q) < 1) {
    # on.exit(return(invisible(NA)))
    stop("INVALID INPUT, at least one quantile must be known!", call. = FALSE)
  }
  if (!is.logical(show.output)) {
    # on.exit(return(invisible(NA)))
    stop("INVALID INPUT, the argument 'show.output' should be logical!", call. = FALSE)
  }
  if (!is.logical(plot)) {
    # on.exit(return(invisible(NA)))
    stop("INVALID INPUT, the argument 'plot' should be logical!", call. = FALSE)
  }
  if (!is.numeric(tol) | length(tol) != 1 | tol < 0) {
    # on.exit(return(invisible(NA)))
    stop("INVALID INPUT, the argument 'tol' should be a single positive numerical value!", call. = FALSE)
  }
  
  #-----------------------------------------------------------------------------
  # minimizinmg procedure
  #-----------------------------------------------------------------------------
  fit.weights.original <- fit.weights
  fit.weights <- fit.weights/sum(fit.weights)
  minimize <- function(par) {
    summand <- suppressWarnings(stats::pexp(q = q, rate = par) - p)
    summand <- summand * fit.weights
    sum(summand^2)
  }
  fit <- c(); fit$value <- tol + 1
  Start <- 1/mean(q)
  try1 <- try(
    fit <- optim(par = Start, 
                 minimize, method = "L-BFGS-B", 
                 lower = 0.001, upper = 10000, hessian=TRUE), 
    silent = TRUE
  )
  
  #-----------------------------------------------------------------------------
  # checking results
  #-----------------------------------------------------------------------------
  if (is.error(try1) || fit$value >= tol) {
    warning("The fitting procedure 'L-BFGS-B' has failed (convergence error occurred or specified tolerance not achieved)!", call. = FALSE)
    fit <- c(); fit$value <- tol + 1
    try2 <- try(
      fit <- optim(par = Start,
                   minimize, 
                   method = "BFGS",
                   hessian=TRUE), 
      silent = TRUE)
    if (is.error(try2) || fit$value >= tol) { 
      warning("The fitting procedure 'BFGS' has failed (convergence error occurred or specified tolerance not achieved)!", call. = FALSE) 
      Par <- NA
    } else if (fit$value < tol) {
      message("The fitting procedure 'BFGS' was successful!\n(Used this fallback optimization method because 'L-BFGS-B' has failed...)") 
      Par <- fit$par
      names(Par) <- c("rate")
      if (show.output) print(fit) 
    }
  } else if (fit$value < tol) {
    message("The fitting procedure 'L-BFGS-B' was successful!") 
    Par <- fit$par
    names(Par) <- c("rate")
    if (show.output) print(fit) 
  }
  
  #-----------------------------------------------------------------------------
  # plotting graphical diagnostics
  #-----------------------------------------------------------------------------
  if (prod(!is.na(Par)) & plot) {
    main1 <- paste("rate = ", round(Par["rate"], digits = 2))
    main <- paste("Exponential (", main1, ")", sep = "")
    sub = paste("fit.weights = c(", paste(fit.weights.original, collapse = ", "), ")", sep = "")
    Support.lim <- c(stats::qexp(p = min(p) * scaleX[1], 
                                 rate = Par["rate"]),
                     stats::qexp(p = (max(p) + (1 - max(p)) * scaleX[2]), 
                                 rate = Par["rate"]))
    Support <- seq(min(min(q), Support.lim[1]), 
                   max(max(q), Support.lim[2]), 
                   length = 200)
    Probability <- stats::pexp(Support, Par["rate"])
    graphics::plot(Support, Probability, type = "l", 
                   xlim = range(Support.lim, q), 
                   main = main, xlab = "Quantiles", 
                   sub = sub, ...)
    graphics::points(x = q, y = p, pch = 19, ...)
  }
  
  #-----------------------------------------------------------------------------
  # output
  #-----------------------------------------------------------------------------
  return(fit)
}
stageC <- c(89,88,72,34,9)
stageT <- c(4,10,26,130,130)

left <- c(0,4,10,26)
right <- stageT[1:4]
################################################################################

################## Initial Theta Value
Init_theta_func <- function(I0_val){
  S0_val <- N-I0_val
  Init_eqn <- function(theta){
    PGFG0(theta,DDist)*N-S0_val
  }
  
  Init_sol <- uniroot(Init_eqn,c(0,1),tol = 1e-9)
  init_theta_val <- 1-Init_sol$root
  return(init_theta_val)
}
it_theta <- Init_theta_func(Init_size)
S0Count <- PGFG0(1-it_theta,DDist)*N

CM_Opt <- ModProc_CM(DDist,Optbeta,Gamma_fit,ODEmaxTime = 400, ODEstep = 2e-1,init_theta = it_theta,TrackDyn = T)
Opt_R0 <- CM_Opt$R0
Opt_RInf <- CM_Opt$RInfinity
Opt_Rsize <- Opt_RInf*N

Sys_Opt <- CM_Opt$Dynamic

######################## POCT Gamma Adjustment #################################

POCT_p_vec <- c(Optp+p_plus)
POCT_Gamma_vec <- c(1:length(POCT_p_vec))
# POCT_p_vec


for (i in c(1:length(POCT_p_vec))) {
  p <- POCT_p_vec[i]
  # Reported is just p proportion of cases
  # All cases not reported will develop to Late Latent stage and become "recovered"
  StageN <- sum(stageC)/p
  StageDC <- stageC[5]/max(stageT)*(right-left)+stageC[1:4]
  StageDP <- c(1:3)
  
  for (j in c(1:3)) {
    StageDP[j] <- sum(StageDC[1:j])/StageN
  }
  # Only reported patient get 9 days TAT adjustment
  StageDQ <- left[2:4]-c(9/14,9/14,0)
  
  
  #### Adding an End Constraint as the third percentile
  StageDP[3] <- 1
  StageDQ[3] <- 26
  
  
  fitexp <- get.eexp.par(StageDP,StageDQ, fit.weight=c(1,1,1),tol = 0.02,show.output = FALSE,plot = F)
  
  if (is.null(fitexp$par) == TRUE){
    gamma <- 0
  } else {
    gamma <- fitexp$par
  }
  POCT_Gamma_vec[i] <- gamma
}


#!!!!!!!!!!VOID: reducing average TAT 9 days
POCT_Gamma_adj_vec <- POCT_Gamma_vec

# print(POCT_Gamma_adj_vec)
################################################################################


### Current Outbreak Starting from fitend time points
Fitend_status <- Sys_Opt[fitend_tps*5+1,]
# print(Fitend_status)

it_theta_now <- as.numeric(1-Fitend_status[2])
it_R_now <- as.numeric(Fitend_status[3])
#ModProc_CM()
CM_Base <- CM_Opt
CM_Base_Out <- as.data.frame(CM_Base$Dynamic)

CM_Base_Out[which.max(CM_Base_Out$I_out),]

POCT_Gamma_adj_vec[1]
POCT_p_vec[1]
POCT_now <- ModProc_CM(DDist,POCT_beta,POCT_Gamma_adj_vec[1], init_theta=it_theta_now, ODEmaxTime = 400-fitend_tps,ODEstep = 2e-1,init_R = it_R_now)
POCT_now_matrix <- POCT_now$Dynamic
POCT_now_matrix[,1] <- POCT_now_matrix[,1]+fitend_tps
POCT_now_Out <- as.data.frame(POCT_now_matrix)

POCT_R0 <- POCT_now$R0
POCT_RInf <- POCT_now$RInfinity
POCT_Rsize <- POCT_now$RInfinity*N
# (POCT_now$RInfinity-CM_Base$RInfinity)*N
POCT_prop <- POCT_RInf/Opt_RInf
# 1-POCT_now$RInfinity/CM_Base$RInfinity

#### ggplot
ggplot()+theme_classic()+
  geom_line(data = CM_Base_Out,aes(x=time,y=I_out,color="Base Model"
                                   #,linetype="Ongoing"
  ),size=0.6)+
    geom_line(data = POCT_now_Out,aes(x=time,y=I_out, color="POCT"
                                         #, linetype="Ongoing"
  ),size=0.6)+
  geom_vline(xintercept = 116,linewidth=0.4,linetype=1,show.legend =TRUE)+
  geom_vline(xintercept = fitend_tps,linewidth=0.6,color="red",linetype=2,show.legend = TRUE)+
  labs(color="Rate Change"
       , linetype="Implement Time"
       , x="Time (bi-week)", y="Active"~italic('I')~"proportion")+
  theme(legend.position = "bottom")+
  xlim(0,400)+
  scale_color_manual(  labels=c(  bquote("POCT "~italic('p')*+.(p_plus)*" TAT"-.(TAT_time)*"d")
                                , "Optimal Base Model")
                     , breaks=c(  "POCT"
                                , "Base Model")
                     , values=c(  "Base Model"="black"
                                , "POCT"="red")
                     )+
  annotate("text",x=116,y=0.012,label="t=116 data end",size=2)+
  annotate("text",x=fitend_tps,y=0.011, label=bquote("t="*.(fitend_tps)*" POCT implemented"),size=2,color="red")
  guides(color=guide_legend(ncol=2,bycol=T))

