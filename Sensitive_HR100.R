rm(list=ls())

wd <- getwd()
### Package Part
library(pracma)
library(gsl)
library(deSolve)
library(igraph)
library(foreach)
library(parallel)
library(iterators)
library(doParallel)
library(ggplot2)
library(rriskDistributions)
library(tidyr)
library(cbinom)
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
#################### Package Part END ###########################################

################################ Distribution part####################################################
{
  
  ########### Degree Distribution
  # input
  
  Degree_frame <- read.csv(paste0(wd,'/Degree Seq 19-23 Init.csv'))
  n <- length(Degree_frame[,1])
  Deg_Seq <- sort(Degree_frame[,1])
  Init_Seq <- sort(Degree_frame[,2])
  #Degree_frame
  
  # data histogram for plot capped at max degree=200
  deg <- c(1:200)
  prob <- c(1:200)
  Deg_data <- data.frame(deg, prob)
  Init_data <- data.frame(deg, prob)
  
  for (i in c(1:200)) {
    Deg_data[i,2] <- length(which(Deg_Seq==i))/length(Deg_Seq)
    Init_data[i,2]  <- length(which(Init_Seq==i))/length(Init_Seq)
  }
  #Deg_data
  #Init_data
  
  # fit the power law at x_min=1
  # Using the plfit :
  # Power laws, Pareto distributions and Zipf's law, M. E. J. Newman, Contemporary Physics, 46, 323-351, 2005.
  # Aaron Clauset, Cosma R .Shalizi and Mark E.J. Newman: Power-law distributions in empirical data. SIAM Review 51(4):661-703, 2009.
  
  # Fit_mle <- fit_power_law(Deg_Seq, xmin = 1, start = 2, force.continuous = FALSE, implementation = c("R.mle"))
  Fit_plfit <- fit_power_law(Deg_Seq, xmin = 1, implementation = 'plfit')
  Init_plfit <- fit_power_law(Init_Seq, xmin = 1, implementation = 'plfit')
  
  # Alpha values
  #alpha_mle <- Fit_mle@coef
  alpha_plfit <- Fit_plfit$alpha
  alpha_Init <- Init_plfit$alpha
  #alpha_mle
  alpha_plfit
  alpha_Init
  
  
  # fitted distribution pmf curves
  Seq_init <- data.frame(deg,prob)
  Seq_plfit <- data.frame(deg,prob)
  
  for (i in c(1:200)) {
    Seq_init[i,2] <- i^(-alpha_Init)
    Seq_plfit[i,2] <- i^(-alpha_plfit)
  }
  
  # Normalize scaling factor
  scal_init <- sum(Seq_init$prob)
  Seq_init$Nprob <- Seq_init$prob/scal_init
  sum(Seq_init$Nprob)
  
  scal_plfit <- sum(Seq_plfit$prob)
  Seq_plfit$Nprob <- Seq_plfit$prob/scal_plfit
  sum(Seq_plfit$Nprob)
  
  # Compare two method for data with outlier 
  
  # plot(Seq_init$deg,Seq_init$Nprob,xlim = c(1,20),type = "b",pch=0, col='black', xlab = 'Degree', ylab = 'Pmf')
  # points(Init_data$prob, type = 'p', pch=1, col='black')
  # points(Seq_plfit$deg, Seq_plfit$Nprob, type = "b",pch=3, lty=2, col='red')
  # points(Deg_data$prob, type = 'p', pch=2, col='red')
  # legend('topright',legend=c('Init_plfit','Init_data','Adj_plfit', 'Adj_data'), col=c('black','black','red', 'red'), pch = c(0,1,3,2))
  # 
  # alpha_Init
  # alpha_plfit
  
  # length(Deg_Seq)
  ## Using plfit method as degree distribution
  kvalue <- append(Seq_plfit$deg,0,0)
  Pk <- append(Seq_plfit$Nprob,0,0)
  DDist <- data.frame(kvalue,Pk)
}
################################ Distribution part END ##############

################################ Continuous Binomial
{
  ### Recoding Continuous Binomial Density by estimation from CDF
  Cum_cb <- function(N,x,p){
    if (x<=0) {
      return(0)
    }
    else if (x>N+1){
      return(1)
    }
    else{
      out <- 1-pbeta(p,x,N+1-x)
      return(out)
    }
  }
  
  Den_cb_est <- function(N,x,p){
    out<-Cum_cb(N,x+1,p)-Cum_cb(N,x,p)
    return(out)
  }
  
  ### Compare with cbinom
  x_seq<-seq(0,11,0.01)
  
  dcb_out<-dcbinom(x_seq,10,0.4)
  
  est_out<-c(1:length(x_seq))
  
  for (i in c(1:length(x_seq))) {
    x <- x_seq[i]
    est_out[i] <- Den_cb_est(10,x,0.4)
  }
}
################################ Continuous Binomial END

################################ Initial Condition
{
Casedata <- read.csv(paste0(wd,"/Time Test Adj.csv"))
tps <- length(Casedata[,1])
R_data <- Casedata$Test.Adj
# Fitting End date week 230, 2023/05/30, point 116
# Afterward POCT is implemented and might change the dynamic
fitend_tps <- 116
R_datafit <- R_data[1:fitend_tps]

## Temporary Population: 23K by Dr. Cater Estmation
N <- 26000

timeseq <- seq(1,(tps+1)*10,10)
fitend_timeseq <- seq(1,(fitend_tps+1)*10,10)
# Initial Condition Solver based on I0=1-S0-(R0=0)
Init_theta_func <- function(I0_val){
  S0_val <- N-I0_val
  Init_eqn <- function(theta){
    PGFG0(theta,DDist)*N-S0_val
  }
  
  Init_sol <- uniroot(Init_eqn,c(0,1),tol = 1e-9)
  init_theta_val <- 1-Init_sol$root
  return(init_theta_val)
}
it_theta <- Init_theta_func(27)
S0Count <- PGFG0(1-it_theta,DDist)*N
S0Count
}

alpha_Init
# Modified get.exp.par function from rriskDistributions
# https://github.com/BfRstats/rriskDistributions/blob/master/R/rriskDistributions_functions.R
# Fitting Lst square of erros
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

############## Fitted Result HPC
#Test500_Init14_3QT <-read.csv(paste0(wd,"/3QT_500_Init14_combined.csv"), header=T)
Test100_Init27_3QT <-read.csv(paste0(wd,"/3QT_HighRes100_Init27_combined.csv"), header=T)
#Test500_Init14_3QT_MASIR <-read.csv(paste0(wd,"/3QT_500_Init14_MASIR_combined.csv"), header=T)
Test100_Init27_3QT_MASIR <-read.csv(paste0(wd,"/3QT_HighRes100_Init27_MASIR_combined.csv"), header=T)
Test100_Init117_3QT_MASIR <-read.csv(paste0(wd,"/3QT_HighRes100_Init117_MASIR_combined.csv"), header=T)


#readmat<-Test25_Init10
readmat <- Test100_Init27_3QT
#readmat
length(readmat[,1])
readmat[1,]

beta_value <- readmat[,4]
gamma_value <- readmat[,3]
p_value <- readmat[,2]
lik_value <- readmat[,5]
log_diff <- -log(max(lik_value)-lik_value+0.001)


df <- data.frame(beta_value,gamma_value,p_value,lik_value,log_diff)
#dfs <- subset(df,df$beta_value<0.15&df$beta_value>0.05)
max(lik_value)
which.max(lik_value)
Opt <- df[which.max(lik_value),]

ggplot(df, aes(x = beta_value, y = p_value, fill = log_diff))+
  geom_tile() + theme_bw()+
  scale_fill_gradient2(low = "white",mid="orange", midpoint = 1, high = "red")+
  labs(fill="Adjusted Loglik",x=""~beta~"(bi-week\u207b\u00b9)", y="p")

Opt
# ggplot(df, aes(x = beta_value, y = p_value, fill = lik_value))+
#   geom_tile(color = "black") +
#   scale_fill_gradient(low = "white", high = "blue")

# ggplot(df, aes(x = beta_value, y = p_value, fill = gamma_value))+
#   geom_tile(color = "black") +
#   scale_fill_gradient(low = "white", high = "blue")


### Optimal Parameter Case
Optbeta <- Opt[1,1]
# beta = 0.0037 bi-weekly
Gamma_fit <- Opt[1,2]
# gamma=0.005329582 bi-weekly
Optp <- Opt[1,3]
# p = 0.14
Optp
Gamma_fit
(1/Gamma_fit)/26*365
Opt
# Checking Error
CM_Opt<- ModProc_CM(DDist,Optbeta,Gamma_fit,ODEmaxTime = 500, ODEstep = 1e-1,init_theta = it_theta,TrackDyn = T)
CM_Opt$R0
CM_Opt$RInfinity
CM_Opt$RInfinity*N

Sys_Opt <- CM_Opt$Dynamic
Opt_times_fit <- Sys_Opt[fitend_timeseq,]
Opt_times <- Sys_Opt[timeseq,]

Rcount_fit <- Opt_times_fit[,3]*N
Opt_times_fit <- cbind(Opt_times_fit,Rcount_fit)
Rcount <- Opt_times[,3]*N
Opt_times <- cbind(Opt_times,Rcount)
inc_count <- rep(0,tps)
inc_count_fit <- rep(0,fitend_tps)
for (i in c(1:fitend_tps)){
  inc_count_fit[i] <- Opt_times_fit[i+1,6]-Opt_times_fit[i,6]
}
for (i in c(1:tps)){
  inc_count[i] <- Opt_times[i+1,6]-Opt_times[i,6]
}

# Incident Compare
inc_exp <- inc_count*Optp

CR_data<-R_data
for (i in c(1:length(R_data))){
  CR_data[i] <- sum(R_data[1:i])
}


t_data <- seq(1,tps)
df_inc <- as.data.frame(cbind(t_data,inc_count,inc_exp))
df_data <- as.data.frame(cbind(t_data,R_data,CR_data))

df_datafit <- subset(df_data,df_data$t_data <= fitend_tps)
df_notfit <- subset(df_data,df_data$t_data > fitend_tps)

# ggplot()+
#   geom_line(data = df_inc,aes(x=t_data, y=inc_count,color="Network Model", linetype="Prediction"))+
#   geom_line(data = df_inc,aes(x=t_data, y=inc_exp,color="Network Model", linetype="Expectation"))+
#   geom_vline(aes(xintercept=fitend_tps),linewidth=0.08,linetype="dotted")+
#   geom_point(data=df_datafit,aes(x=t_data,y=R_data,color="Fitted data"), alpha=0.5)+
#   geom_point(data=df_notfit,aes(x=t_data,y=R_data,color="Unfitted data"), alpha=0.5)+
#   scale_linetype_manual(values=c("Expectation"="dashed","Prediction"="solid"))+
#   scale_color_manual(breaks=c("Network Model","Fitted data","Unfitted data"),values=c("Network Model"="red","Fitted data"="blue","Unfitted data"="green"))+
#   labs(x="Time(bi-week)",y="Casecount",color="Color",linetype="Line Type")

ggplot()+theme_bw()+
  geom_line(data = df_inc,aes(x=t_data, y=inc_exp,color="Network model reporting expectation"))+
  geom_vline(aes(xintercept=fitend_tps),linewidth=0.8,linetype="dotted")+
  geom_point(data=df_datafit,aes(x=t_data,y=R_data,color="Fitted data"), alpha=0.5)+
  geom_point(data=df_notfit,aes(x=t_data,y=R_data,color="Unfitted data"), alpha=0.5)+
  scale_linetype_manual(values=c("Expectation"="dashed","Prediction"="solid"))+
  scale_color_manual(breaks=c("Network model reporting expectation","Fitted data","Unfitted data"),values=c("Network model reporting expectation"="red","Fitted data"="blue","Unfitted data"="green"))+
  labs(x="Time (bi-week)",y="Casecount",color="Color",linetype="Line Type")





# plot(seq(1,tps), inc_count, type='l', col="red", ylim = c(0,20),xlab ="t", ylab='Incident R')
# points(seq(1,tps), inc_count*Optp, type='l', col="green")
# abline(v=116)
# points(seq(1,fitend_tps),R_datafit)
# points(seq(fitend_tps+1,tps), R_data[seq(fitend_tps+1,tps)], col='blue')
# legend('topleft',legend=c('Prediction Curve', 'Expected Reporting Curve','fitted data','unfitted data'), col=c('red','green','black','blue'), lty=c(1,1,0,0), pch = c(NA,NA,1,1),cex=0.5)

# Cumulative Compare
#R_data

t_Sys <- Sys_Opt[,1]
CR_Sys <- Sys_Opt[,3]*N
CR_Exp <- Sys_Opt[,3]*N*Optp
R_Sys <- Sys_Opt[,3]
I_Sys <- Sys_Opt[,5]
df_Sys <- as.data.frame(cbind(t_Sys,CR_Sys,CR_Exp,R_Sys,I_Sys))

# Compare fitting
# ggplot()+
#   geom_line(data = df_Sys,aes(x=t_Sys, y=CR_Sys,color="Network Model", linetype="Prediction"))+
#   geom_line(data = df_Sys,aes(x=t_Sys, y=CR_Exp,color="Network Model", linetype="Expectation"))+
#   geom_vline(aes(xintercept=fitend_tps),linewidth=0.08,linetype="dotted")+
#   geom_point(data=df_datafit,aes(x=t_data,y=CR_data,color="Fitted data"), alpha=0.25)+
#   geom_point(data=df_notfit,aes(x=t_data,y=CR_data,color="Unfitted data"), alpha=0.25)+
#   scale_linetype_manual(values=c("Expectation"="dashed","Prediction"="solid"))+
#   scale_color_manual(breaks=c("Network Model","Fitted data","Unfitted data"),values=c("Network Model"="red","Fitted data"="blue","Unfitted data"="green"))+
#   labs(x="Time(bi-week)",y="Cumulative R",color="Color",linetype="Line Type")+
#   xlim(0,tps)+
#   ylim(0,1000)

ggplot()+theme_bw()+
  geom_line(data = df_Sys,aes(x=t_Sys, y=CR_Exp,color="Network model reporting expectation"))+
  geom_vline(aes(xintercept=fitend_tps),linewidth=0.8,linetype="dotted")+
  geom_point(data=df_datafit,aes(x=t_data,y=CR_data,color="Fitted data"), alpha=0.25)+
  geom_point(data=df_notfit,aes(x=t_data,y=CR_data,color="Unfitted data"), alpha=0.25)+
  scale_linetype_manual(values=c("Expectation"="dashed","Prediction"="solid"))+
  scale_color_manual(breaks=c("Network model reporting expectation","Fitted data","Unfitted data"),values=c("Network model reporting expectation"="red","Fitted data"="blue","Unfitted data"="green"))+
  labs(x="Time (bi-week)",y="Cumulative R",color="Color",linetype="Line Type")+
  xlim(0,tps)


# Tendency
ggplot()+theme_bw()+
  geom_line(data = df_Sys,aes(x=t_Sys, y=R_Sys,color="R"))+
  geom_line(data = df_Sys,aes(x=t_Sys, y=I_Sys,color="I"))+
  geom_vline(aes(xintercept=fitend_tps),linewidth=0.5,linetype="dotted")+
  geom_vline(aes(xintercept=tps),linewidth=0.5,linetype="dotted")+
  labs(x="Time (bi-week)",y="Proportion",color="Compartment")+
  xlim(0,350)




# plot(Sys_Opt[,1],Sys_Opt[,3]*N, col="red", ylim = c(0,400), xlab ="t (bi-week)", ylab='Cumulative R')
# abline(v=116)
# points(Sys_Opt[,1], Sys_Opt[,3]*N*Optp, col="green")
# points(seq(1,fitend_tps),CR_data[seq(1,fitend_tps)], col='black', pch=1)
# points(seq(fitend_tps+1,tps), CR_data[seq(fitend_tps+1,tps)], col='blue',pch=1)
# legend('topleft',legend=c('Prediction Curve','Expected Reporting Curve','fitted data','unfitted data'), col=c('red','green','black','blue'), pch = c(19,19,1,1),cex=0.5)

#### Fully mixed/Mass Action SIR Model
MASIR_Proc <- function(beta, gamma, init_S=1e-3, ODEmaxTime=50, ODEstep=1e-2,TrackDyn=T){
  if (TrackDyn==T){
    Sys <- function(t, y, parms){
      with(as.list(c(parms,y)),{
        dS <- (-b)*I*S
        dI <- b*I*S-g*I
        dR <- g*I
        return(list(c(dS,dI,dR)))
      })
    }
    parms <- c(b=beta,g=gamma)
    times <- seq(0,ODEmaxTime,by=ODEstep)
    y <- c(S=init_S,I=1-init_S,R=0)
    
    Sys_out <- ode(y,times,Sys,parms)
  }
  
  g <- gamma
  b <- beta
  S_0 <- init_S
  R_0 <- 0
  
  R0 <- b/g
  RInf <- Sys_out[length(Sys_out[,4]),4]
  
  if (TrackDyn==T){
    return(list(R0=R0,RInfinity=RInf, Dynamic=Sys_out))
  } else {
    return(list(R0=R0,RInfinity=RInf))
  }
}

#readmat<-Test25_Init10
readmat_MASIR <- Test100_Init27_3QT_MASIR
readmat_MBSIR <- Test100_Init117_3QT_MASIR

#readmat
beta_value_MASIR <- readmat_MASIR[,4]
gamma_value_MASIR <- readmat_MASIR[,3]
p_value_MASIR <- readmat_MASIR[,2]
lik_value_MASIR <- readmat_MASIR[,5]
log_diff_MASIR <- -log(max(lik_value_MASIR)-lik_value_MASIR+0.001)
R0_value_MASIR <- beta_value_MASIR/gamma_value_MASIR

df_MASIR <- data.frame(beta_value_MASIR,gamma_value_MASIR,p_value_MASIR,lik_value_MASIR,log_diff_MASIR,R0_value_MASIR)
#dfs_MASIR <- subset(df_MASIR,df_MASIR$beta_value_MASIR<0.15&df_MASIR$beta_value_MASIR>0.05)
max(lik_value_MASIR)
which.max(lik_value_MASIR)
Opt_MASIR <- df_MASIR[which.max(lik_value_MASIR),]

ggplot(df_MASIR, aes(x = beta_value_MASIR, y = p_value_MASIR, fill = log_diff_MASIR))+
  geom_tile() + theme_bw()+
  scale_fill_gradient2(low = "white",mid="orange", midpoint = 3, high = "red")+
  labs(fill="Adjusted Loglik",x=""~beta~"(bi-week\u207b\u00b9)", y="p")

# Init117
beta_value_MBSIR <- readmat_MBSIR[,4]
gamma_value_MBSIR <- readmat_MBSIR[,3]
p_value_MBSIR <- readmat_MBSIR[,2]
lik_value_MBSIR <- readmat_MBSIR[,5]
log_diff_MBSIR <- -log(max(lik_value_MBSIR)-lik_value_MBSIR+0.001)
R0_value_MBSIR <- beta_value_MBSIR/gamma_value_MBSIR

df_MBSIR <- data.frame(beta_value_MBSIR,gamma_value_MBSIR,p_value_MBSIR,lik_value_MBSIR,log_diff_MBSIR,R0_value_MBSIR)
#dfs_MBSIR <- subset(df_MBSIR,df_MBSIR$beta_value_MBSIR<0.15&df_MBSIR$beta_value_MBSIR>0.05)
max(lik_value_MBSIR)
which.max(lik_value_MBSIR)
Opt_MBSIR <- df_MBSIR[which.max(lik_value_MBSIR),]

ggplot(df_MBSIR, aes(x = beta_value_MBSIR, y = p_value_MBSIR, fill = log_diff_MBSIR))+
  geom_tile() + theme_bw()+
  scale_fill_gradient2(low = "white",mid="orange", midpoint = 3, high = "red")+
  labs(fill="Adjusted Loglik",x=""~beta~"(bi-week\u207b\u00b9)", y="p")



# ggplot(df, aes(x = beta_value, y = p_value, fill = lik_value))+
#   geom_tile(color = "black") +
#   scale_fill_gradient(low = "white", high = "blue")
# 
# ggplot(df, aes(x = beta_value, y = p_value, fill = gamma_value))+
#   geom_tile(color = "black") +
#   scale_fill_gradient(low = "white", high = "blue")




Opt_MASIR

### Optimal Parameter Case
MAOptbeta <- Opt_MASIR[1,1]
# beta = 0.0037 bi-weekly
MAGamma_fit <- Opt_MASIR[1,2]
# gamma=0.005329582 bi-weekly
MAOptp <- Opt_MASIR[1,3]
# p = 0.14

MAGamma_fit

# Checking Error
I0_MA <-27
MA_Opt<- MASIR_Proc(MAOptbeta, MAGamma_fit, init_S = (N-I0_MA)/N, ODEmaxTime=500, ODEstep=1e-1,TrackDyn = T)
MA_Opt$R0
MA_Opt$RInfinity
MA_Opt$RInfinity*N

MASys_Opt <- MA_Opt$Dynamic
MAOpt_times_fit <- MASys_Opt[fitend_timeseq,]
MAOpt_times <- MASys_Opt[timeseq,]

MAOpt_times_fit
MARcount_fit <- MAOpt_times_fit[,4]*N
MAOpt_times_fit <- cbind(MAOpt_times_fit,MARcount_fit)
MARcount <- MAOpt_times[,4]*N
MAOpt_times <- cbind(MAOpt_times,MARcount)


MAinc_count <- rep(0,tps)
MAinc_count_fit <- rep(0,fitend_tps)
for (i in c(1:fitend_tps)){
  MAinc_count_fit[i] <- MAOpt_times_fit[i+1,5]-MAOpt_times_fit[i,5]
}
for (i in c(1:tps)){
  MAinc_count[i] <- MAOpt_times[i+1,5]-MAOpt_times[i,5]
}

LogLik <- c(1:fitend_tps)
MAinc_count_fit
for (k in c(1:fitend_tps)){
  LogLik[k] <- log(Den_cb_est(MAinc_count_fit[k],R_datafit[k],MAOptp))
}
sum(LogLik)

# Init117
Opt_MBSIR

### Optimal Parameter Case
MBOptbeta <- Opt_MBSIR[1,1]
# beta = 0.0037 bi-weekly
MBGamma_fit <- Opt_MBSIR[1,2]
# gamma=0.005329582 bi-weekly
MBOptp <- Opt_MBSIR[1,3]
# p = 0.14

MBGamma_fit

# Checking Error
I0_MB <- 117
MB_Opt<- MASIR_Proc(MBOptbeta, MBGamma_fit, init_S = (N-I0_MB)/N, ODEmaxTime=500, ODEstep=1e-1,TrackDyn = T)
MB_Opt$R0
MB_Opt$RInfinity
MB_Opt$RInfinity*N

MBSys_Opt <- MB_Opt$Dynamic
MBOpt_times_fit <- MBSys_Opt[fitend_timeseq,]
MBOpt_times <- MBSys_Opt[timeseq,]

MBOpt_times_fit
MBRcount_fit <- MBOpt_times_fit[,4]*N
MBOpt_times_fit <- cbind(MBOpt_times_fit,MBRcount_fit)
MBRcount <- MBOpt_times[,4]*N
MBOpt_times <- cbind(MBOpt_times,MBRcount)


MBinc_count <- rep(0,tps)
MBinc_count_fit <- rep(0,fitend_tps)
for (i in c(1:fitend_tps)){
  MBinc_count_fit[i] <- MBOpt_times_fit[i+1,5]-MBOpt_times_fit[i,5]
}
for (i in c(1:tps)){
  MBinc_count[i] <- MBOpt_times[i+1,5]-MBOpt_times[i,5]
}

LogLik <- c(1:fitend_tps)
MBinc_count_fit
for (k in c(1:fitend_tps)){
  LogLik[k] <- log(Den_cb_est(MBinc_count_fit[k],R_datafit[k],MBOptp))
}
sum(LogLik)
Opt_MBSIR

# Incident Compare
MAinc_exp <- MAinc_count*MAOptp
MBinc_exp <- MBinc_count*MBOptp

df_MAinc <- as.data.frame(cbind(t_data,MAinc_count,MAinc_exp))
df_MBinc <- as.data.frame(cbind(t_data,MBinc_count,MBinc_exp))

# ggplot()+
#   geom_line(data = df_inc,aes(x=t_data, y=inc_count,color="Network Model", linetype="Prediction"))+
#   geom_line(data = df_inc,aes(x=t_data, y=inc_exp,color="Network Model", linetype="Expectation"))+
#   geom_line(data = df_MAinc,aes(x=t_data, y=MAinc_count,color="MA I0=23", linetype="Prediction"), alpha=0.75)+
#   geom_line(data = df_MAinc,aes(x=t_data, y=MAinc_exp,color="MA I0=23", linetype="Expectation"))+
#   geom_line(data = df_MBinc,aes(x=t_data, y=MBinc_count,color="MA I0=100", linetype="Prediction"))+
#   geom_line(data = df_MBinc,aes(x=t_data, y=MBinc_exp,color="MA I0=100", linetype="Expectation"))+
#   geom_vline(aes(xintercept=fitend_tps),linewidth=0.08,linetype="dotted")+
#   geom_point(data=df_datafit,aes(x=t_data,y=R_data,color="Fitted data"), alpha=0.15)+
#   geom_point(data=df_notfit,aes(x=t_data,y=R_data,color="Unfitted data"), alpha=0.15)+
#   scale_linetype_manual(values=c("Expectation"="dashed","Prediction"="solid"))+
#   scale_color_manual(labels=c("Network Model",bquote("MA "~I[0]*"=23"), bquote("MA "~I[0]*"=100"),"Fitted data","Unfitted data"),breaks=c("Network Model","MA I0=23", "MA I0=100","Fitted data","Unfitted data"),values=c("Network Model"="red", "MA I0=23"="orange","MA I0=100"="purple","Fitted data"="blue","Unfitted data"="green"))+
#   labs(x="Time(bi-week)",y="Casecount",color="Color",linetype="Line Type")

ggplot()+theme_bw()+
  #geom_line(data = df_inc,aes(x=t_data, y=inc_count,color="Network model reported expectation", linetype="Prediction"))+
  geom_line(data = df_inc,aes(x=t_data, y=inc_exp,color="Network model reported expectation"))+
  #geom_line(data = df_MAinc,aes(x=t_data, y=MAinc_count,color="MA I0=23", linetype="Prediction"), alpha=0.75)+
  geom_line(data = df_MAinc,aes(x=t_data, y=MAinc_exp,color="MA I0=27"))+
  #geom_line(data = df_MBinc,aes(x=t_data, y=MBinc_count,color="MA I0=100", linetype="Prediction"))+
  geom_line(data = df_MBinc,aes(x=t_data, y=MBinc_exp,color="MA I0=117"))+
  geom_vline(aes(xintercept=fitend_tps),linewidth=0.5,linetype="dotted")+
  geom_point(data=df_datafit,aes(x=t_data,y=R_data,color="Fitted data"), alpha=0.4)+
  geom_point(data=df_notfit,aes(x=t_data,y=R_data,color="Unfitted data"), alpha=0.4)+
  scale_linetype_manual(values=c("Expectation"="dashed","Prediction"="solid"))+
  scale_color_manual(labels=c("Network model reported expectation",bquote("MA "~I[0]*"=27 model reported expectation") , bquote("MA "~I[0]*"=117 model reported expectation"),"Fitted data","Unfitted data"),breaks=c("Network model reported expectation","MA I0=27", "MA I0=117","Fitted data","Unfitted data"),values=c("Network model reported expectation"="red", "MA I0=27"="orange","MA I0=117"="purple","Fitted data"="blue","Unfitted data"="green"))+
  labs(x="Time (bi-week)",y="Casecount",color="Color",linetype="Line Type")


# Cumulative Compare
#R_data
MASys_Opt[1,]
t_MASys <- MASys_Opt[,1]
t_MBSys <- MBSys_Opt[,1]
CR_MASys <- MASys_Opt[,4]*N
CR_MBSys <- MBSys_Opt[,4]*N
CR_MAExp <- MASys_Opt[,4]*N*MAOptp
CR_MBExp <- MBSys_Opt[,4]*N*MBOptp

R_MASys <- MASys_Opt[,4]
I_MASys <- MASys_Opt[,3]
R_MBSys <- MBSys_Opt[,4]
I_MBSys <- MBSys_Opt[,3]

df_MASys <- as.data.frame(cbind(t_MASys,CR_MASys,CR_MAExp,R_MASys,I_MASys))
df_MBSys <- as.data.frame(cbind(t_MBSys,CR_MBSys,CR_MBExp,R_MBSys,I_MBSys))

# Compare fitting
# ggplot()+
#   geom_line(data = df_Sys,aes(x=t_Sys, y=CR_Sys,color="Network Model", linetype="Prediction"))+
#   geom_line(data = df_Sys,aes(x=t_Sys, y=CR_Exp,color="Network Model", linetype="Expectation"))+
#   geom_line(data = df_MASys,aes(x=t_MASys, y=CR_MASys,color="MA I0=23", linetype="Prediction"))+
#   geom_line(data = df_MASys,aes(x=t_MASys, y=CR_MAExp,color="MA I0=23", linetype="Expectation"))+
#   geom_line(data = df_MBSys,aes(x=t_MBSys, y=CR_MBSys,color="MA I0=100", linetype="Prediction"))+
#   geom_line(data = df_MBSys,aes(x=t_MBSys, y=CR_MBExp,color="MA I0=100", linetype="Expectation"))+
#   geom_vline(aes(xintercept=fitend_tps),linewidth=0.08,linetype="dotted")+
#   geom_point(data=df_datafit,aes(x=t_data,y=CR_data,color="Fitted data"), alpha=0.25)+
#   geom_point(data=df_notfit,aes(x=t_data,y=CR_data,color="Unfitted data"), alpha=0.25)+
#   scale_linetype_manual(values=c("Expectation"="dashed","Prediction"="solid"))+
#   scale_color_manual(labels=c("Network Model",bquote("MA "~I[0]*"=23"), bquote("MA "~I[0]*"=100"),"Fitted data","Unfitted data"),breaks=c("Network Model","MA I0=23", "MA I0=100","Fitted data","Unfitted data"),values=c("Network Model"="red", "MA I0=23"="orange","MA I0=100"="purple","Fitted data"="blue","Unfitted data"="green"))+
#   labs(x="Time(bi-week)",y="Cumulative R",color="Color",linetype="Line Type")+
#   xlim(0,tps)+
#   ylim(0,800)

ggplot()+theme_bw()+
  #geom_line(data = df_Sys,aes(x=t_Sys, y=CR_Sys,color="Network Model", linetype="Prediction"))+
  geom_line(data = df_Sys,aes(x=t_Sys, y=CR_Exp,color="Network model reported expectation"))+
  #geom_line(data = df_MASys,aes(x=t_MASys, y=CR_MASys,color="MA I0=23", linetype="Prediction"))+
  geom_line(data = df_MASys,aes(x=t_MASys, y=CR_MAExp,color="MA I0=27"))+
  #geom_line(data = df_MBSys,aes(x=t_MBSys, y=CR_MBSys,color="MA I0=100", linetype="Prediction"))+
  geom_line(data = df_MBSys,aes(x=t_MBSys, y=CR_MBExp,color="MA I0=117"))+
  geom_vline(aes(xintercept=fitend_tps),linewidth=0.5,linetype="dotted")+
  geom_point(data=df_datafit,aes(x=t_data,y=CR_data,color="Fitted data"), alpha=0.25)+
  geom_point(data=df_notfit,aes(x=t_data,y=CR_data,color="Unfitted data"), alpha=0.25)+
  scale_linetype_manual(values=c("Expectation"="dashed","Prediction"="solid"))+
  scale_color_manual(labels=c("Network model reported expectation",bquote("MA "~I[0]*"=27 model reported expectation") , bquote("MA "~I[0]*"=117 model reported expectation"),"Fitted data","Unfitted data"),breaks=c("Network model reported expectation","MA I0=27", "MA I0=117","Fitted data","Unfitted data"),values=c("Network model reported expectation"="red", "MA I0=27"="orange","MA I0=117"="purple","Fitted data"="blue","Unfitted data"="green"))+
  labs(x="Time (bi-week)",y="Casecount",color="Color",linetype="Line Type")+
  xlim(0,tps)+
  ylim(0,400)


# R Tendency
ggplot()+theme_bw()+
  geom_line(data = df_Sys,aes(x=t_Sys, y=R_Sys,color="Network Model"))+
  geom_line(data = df_MASys,aes(x=t_MASys, y=R_MASys,color="MA I0=27"))+
  geom_line(data = df_MBSys,aes(x=t_MBSys, y=R_MBSys,color="MA I0=117"))+
  geom_vline(aes(xintercept=fitend_tps),linewidth=0.5,linetype="dotted")+
  geom_vline(aes(xintercept=tps),linewidth=0.5,linetype="dotted")+
  scale_color_manual(labels=c("Network Model",bquote("MA "~I[0]*"=27"), bquote("MA "~I[0]*"=117")),breaks=c("Network Model","MA I0=27", "MA I0=117","Fitted data","Unfitted data"),values=c("Network Model"="red", "MA I0=27"="orange","MA I0=117"="purple"))+
  labs(x="Time (bi-week)",y="R Proportion",color="Color")+
  xlim(0,500)

# I Tendency
ggplot()+theme_bw()+
  geom_line(data = df_Sys,aes(x=t_Sys, y=I_Sys,color="Network Model"))+
  geom_line(data = df_MASys,aes(x=t_MASys, y=I_MASys,color="MA I0=27"))+
  geom_line(data = df_MBSys,aes(x=t_MBSys, y=I_MBSys,color="MA I0=117"))+
  geom_vline(aes(xintercept=fitend_tps),linewidth=0.5,linetype="dotted")+
  geom_vline(aes(xintercept=tps),linewidth=0.5,linetype="dotted")+
  scale_color_manual(labels=c("Network Model",bquote("MA "~I[0]*"=27"), bquote("MA "~I[0]*"=117")),breaks=c("Network Model","MA I0=27", "MA I0=117","Fitted data","Unfitted data"),values=c("Network Model"="red", "MA I0=27"="orange","MA I0=117"="purple"))+
  labs(x="Time (bi-week)",y="I Proportion",color="Color")+
  xlim(0,500)


####################### Confidence Interval
# delta in parameters
d <- 1e-5

# MSV
Opt

# beta derivative
MSV_Opt_plusBeta <- ModProc_CM(DDist,Optbeta+d/2,Gamma_fit,ODEmaxTime = 500, ODEstep = 1e-1,init_theta = it_theta,TrackDyn = T)
MSV_R_Opt_plusBeta <- MSV_Opt_plusBeta$Dynamic[fitend_timeseq,3]*N
MSV_Opt_minusBeta <- ModProc_CM(DDist,Optbeta-d/2,Gamma_fit,ODEmaxTime = 500, ODEstep = 1e-1,init_theta = it_theta,TrackDyn = T)
MSV_R_Opt_minusBeta <- MSV_Opt_minusBeta$Dynamic[fitend_timeseq,3]*N

# P derivative
################### Fitting gamma based on p
Optp
MSV_p_vec <- c(Optp-d/2,Optp,Optp+d/2)
MSV_Gamma_vec <- c(1:length(MSV_p_vec))
MSV_p_vec
for (i in c(1:length(MSV_p_vec))) {
  p <- MSV_p_vec[i]
  # Reported is just p proportion of cases
  # All cases not reported will develop to Late Latent stage and become "recovered"
  StageN <- sum(stageC)/p
  StageDC <- stageC[5]/max(stageT)*(right-left)+stageC[1:4]
  StageDP <- c(1:3)
  
  for (j in c(1:3)) {
    StageDP[j] <- sum(StageDC[1:j])/StageN
  }
  StageDQ <- left[2:4]
  
  #### Adding an End Constraint as the third percentile
  StageDP[3] <- 1
  StageDQ[3] <- 26
  
  
  fitexp <- get.eexp.par(StageDP,StageDQ, fit.weight=c(1,1,1),tol = 0.02,show.output = FALSE,plot = F)
  
  if (is.null(fitexp$par) == TRUE){
    gamma <- 0
  } else {
    gamma <- fitexp$par
  }
  MSV_Gamma_vec[i] <- gamma
}

MSV_Opt_plusP <- ModProc_CM(DDist,Optbeta,MSV_Gamma_vec[3],ODEmaxTime = 500, ODEstep = 1e-1,init_theta = it_theta,TrackDyn = T)
MSV_R_Opt_plusP <- MSV_Opt_plusP$Dynamic[fitend_timeseq,3]*N
MSV_Opt_minusP <- ModProc_CM(DDist,Optbeta,MSV_Gamma_vec[1],ODEmaxTime = 500, ODEstep = 1e-1,init_theta = it_theta,TrackDyn = T)
MSV_R_Opt_minusP <- MSV_Opt_minusP$Dynamic[fitend_timeseq,3]*N


# Jacobian Jf
MSV_loglik_Opt_plusBeta <- rep(0,fitend_tps)
MSV_loglik_Opt_minusBeta <- rep(0,fitend_tps)
MSV_loglik_Opt_plusP <- rep(0,fitend_tps)
MSV_loglik_Opt_minusP <- rep(0,fitend_tps)

for (i in c(1:fitend_tps)){
  MSV_loglik_Opt_plusBeta[i] <- Den_cb_est(MSV_R_Opt_plusBeta[i+1]-MSV_R_Opt_plusBeta[i], R_datafit[i], MSV_p_vec[2])
  MSV_loglik_Opt_minusBeta[i] <- Den_cb_est(MSV_R_Opt_minusBeta[i+1]-MSV_R_Opt_minusBeta[i], R_datafit[i], MSV_p_vec[2])
  MSV_loglik_Opt_plusP[i] <- Den_cb_est(MSV_R_Opt_plusP[i+1]-MSV_R_Opt_plusP[i], R_datafit[i], MSV_p_vec[3])
  MSV_loglik_Opt_minusP[i] <- Den_cb_est(MSV_R_Opt_minusP[i+1]-MSV_R_Opt_minusP[i], R_datafit[i], MSV_p_vec[1])
}

MSV_Jf_Betarow <- (MSV_loglik_Opt_plusBeta-MSV_loglik_Opt_minusBeta)/d
MSV_Jf_Prow <- (MSV_loglik_Opt_plusP-MSV_loglik_Opt_minusP)/d

MSV_Jf <- cbind(MSV_Jf_Betarow,MSV_Jf_Prow)

MSV_GTGmat <- t(MSV_Jf) %*% MSV_Jf
#MSV_Qmat <- 1/fitend_tps*MSV_Hmat

MSV_V0mat <- inv(MSV_GTGmat)
MSV_V0mat
# Confidence Interval of Single parameter
MSV_V_Beta <- MSV_V0mat[1,1]
MSV_V_P <- MSV_V0mat[2,2]
MSV_Beta_UpStat <- qnorm(0.975, 0, sqrt(MSV_V_Beta))
MSV_Beta_LowStat <- qnorm(0.025, 0, sqrt(MSV_V_Beta))
MSV_P_UpStat <- qnorm(0.975, 0, sqrt(MSV_V_P))
MSV_P_LowStat <- qnorm(0.025, 0, sqrt(MSV_V_P))

# Beta
HRn <-100
MSV_Beta_CI_seq <- seq(0,0.00228,1e-6)
MSV_Beta_teststat <- sqrt(fitend_tps)*(Optbeta-MSV_Beta_CI_seq)
MSV_Beta_LB <- min(MSV_Beta_CI_seq[which(MSV_Beta_LowStat < MSV_Beta_teststat & MSV_Beta_UpStat > MSV_Beta_teststat)])
MSV_Beta_UB <- max(MSV_Beta_CI_seq[which(MSV_Beta_LowStat < MSV_Beta_teststat & MSV_Beta_UpStat > MSV_Beta_teststat)])

MSV_Beta_LB
MSV_Beta_UB
Optbeta

# P
MSV_P_CI_seq <- seq(0,1,1e-4)
MSV_P_teststat <- sqrt(fitend_tps)*(Optp-MSV_P_CI_seq)
MSV_P_LB <- min(MSV_P_CI_seq[which(MSV_P_LowStat < MSV_P_teststat & MSV_P_UpStat > MSV_P_teststat)])
MSV_P_UB <- max(MSV_P_CI_seq[which(MSV_P_LowStat < MSV_P_teststat & MSV_P_UpStat > MSV_P_teststat)])

MSV_P_LB
MSV_P_UB
Optp

# MASIR
Opt_MASIR
MAGamma_fit
MAOptp
# beta derivative
MASIR_Opt_plusBeta <- MASIR_Proc(MAOptbeta+d/2,MAGamma_fit,init_S = (N-I0_MA)/N, ODEmaxTime = fitend_tps, ODEstep = 1e-1,TrackDyn = T)
MASIR_R_Opt_plusBeta <- MASIR_Opt_plusBeta$Dynamic[fitend_timeseq,4]*N
MASIR_Opt_minusBeta <- MASIR_Proc(MAOptbeta-d/2,MAGamma_fit,init_S = (N-I0_MA)/N, ODEmaxTime = fitend_tps, ODEstep = 1e-1,TrackDyn = T)
MASIR_R_Opt_minusBeta <- MASIR_Opt_minusBeta$Dynamic[fitend_timeseq,4]*N

# P derivative
################### Fitting gamma based on p
MAOptp
MASIR_p_vec <- c(MAOptp-d/2,MAOptp,MAOptp+d/2)
MASIR_Gamma_vec <- c(1:length(MASIR_p_vec))
MASIR_p_vec
for (i in c(1:length(MASIR_p_vec))) {
  p <- MASIR_p_vec[i]
  # Reported is just p proportion of cases
  # All cases not reported will develop to Late Latent stage and become "recovered"
  StageN <- sum(stageC)/p
  StageDC <- stageC[5]/max(stageT)*(right-left)+stageC[1:4]
  StageDP <- c(1:3)
  
  for (j in c(1:3)) {
    StageDP[j] <- sum(StageDC[1:j])/StageN
  }
  StageDQ <- left[2:4]
  
  #### Adding an End Constraint as the third percentile
  StageDP[3] <- 1
  StageDQ[3] <- 26
  
  
  fitexp <- get.eexp.par(StageDP,StageDQ, fit.weight=c(1,1,1),tol = 0.02,show.output = FALSE,plot = F)
  
  if (is.null(fitexp$par) == TRUE){
    gamma <- 0
  } else {
    gamma <- fitexp$par
  }
  MASIR_Gamma_vec[i] <- gamma
}

MASIR_Opt_plusP <- MASIR_Proc(MAOptbeta,MASIR_Gamma_vec[3],init_S = (N-I0_MA)/N, ODEmaxTime = fitend_tps, ODEstep = 1e-1,TrackDyn = T)
MASIR_R_Opt_plusP <- MASIR_Opt_plusP$Dynamic[fitend_timeseq,4]*N
MASIR_Opt_minusP <- MASIR_Proc(MAOptbeta,MASIR_Gamma_vec[1],init_S = (N-I0_MA)/N, ODEmaxTime = fitend_tps, ODEstep = 1e-1,TrackDyn = T)
MASIR_R_Opt_minusP <- MASIR_Opt_minusP$Dynamic[fitend_timeseq,4]*N

# Jacobian Jf
MASIR_loglik_Opt_plusBeta <- rep(0,fitend_tps)
MASIR_loglik_Opt_minusBeta <- rep(0,fitend_tps)
MASIR_loglik_Opt_plusP <- rep(0,fitend_tps)
MASIR_loglik_Opt_minusP <- rep(0,fitend_tps)

for (i in c(1:fitend_tps)){
  MASIR_loglik_Opt_plusBeta[i] <- Den_cb_est(MASIR_R_Opt_plusBeta[i+1]-MASIR_R_Opt_plusBeta[i], R_datafit[i], MASIR_p_vec[2])
  MASIR_loglik_Opt_minusBeta[i] <- Den_cb_est(MASIR_R_Opt_minusBeta[i+1]-MASIR_R_Opt_minusBeta[i], R_datafit[i], MASIR_p_vec[2])
  MASIR_loglik_Opt_plusP[i] <- Den_cb_est(MASIR_R_Opt_plusP[i+1]-MASIR_R_Opt_plusP[i], R_datafit[i], MASIR_p_vec[3])
  MASIR_loglik_Opt_minusP[i] <- Den_cb_est(MASIR_R_Opt_minusP[i+1]-MASIR_R_Opt_minusP[i], R_datafit[i], MASIR_p_vec[1])
}

MASIR_Jf_Betarow <- (MASIR_loglik_Opt_plusBeta-MASIR_loglik_Opt_minusBeta)/d
MASIR_Jf_Prow <- (MASIR_loglik_Opt_plusP-MASIR_loglik_Opt_minusP)/d

MASIR_Jf <- cbind(MASIR_Jf_Betarow,MASIR_Jf_Prow)

MASIR_GTGmat <- t(MASIR_Jf) %*% MASIR_Jf  
#MASIR_Qmat <- 1/fitend_tps*MASIR_Hmat

MASIR_V0mat <- inv(MASIR_GTGmat)

# Confidence Interval of Single parameter
MASIR_V_Beta <- MASIR_V0mat[1,1]
MASIR_V_P <- MASIR_V0mat[2,2]
MASIR_Beta_UpStat <- qnorm(0.975, 0, sqrt(MASIR_V_Beta))
MASIR_Beta_LowStat <- qnorm(0.025, 0, sqrt(MASIR_V_Beta))
MASIR_P_UpStat <- qnorm(0.975, 0, sqrt(MASIR_V_P))
MASIR_P_LowStat <- qnorm(0.025, 0, sqrt(MASIR_V_P))

# Beta
MASIR_Beta_CI_seq <- seq(0.03,0.14,1e-5)
MASIR_Beta_teststat <- sqrt(fitend_tps)*(MAOptbeta-MASIR_Beta_CI_seq)
MASIR_Beta_LB <- min(MASIR_Beta_CI_seq[which(MASIR_Beta_LowStat < MASIR_Beta_teststat & MASIR_Beta_UpStat > MASIR_Beta_teststat)])
MASIR_Beta_UB <- max(MASIR_Beta_CI_seq[which(MASIR_Beta_LowStat < MASIR_Beta_teststat & MASIR_Beta_UpStat > MASIR_Beta_teststat)])

MASIR_Beta_LB
MASIR_Beta_UB
MAOptbeta

# P
MASIR_P_CI_seq <- seq(0,0.5,1e-4)
MASIR_P_teststat <- sqrt(fitend_tps)*(MAOptp-MASIR_P_CI_seq)
MASIR_P_LB <- min(MASIR_P_CI_seq[which(MASIR_P_LowStat < MASIR_P_teststat & MASIR_P_UpStat > MASIR_P_teststat)])
MASIR_P_UB <- max(MASIR_P_CI_seq[which(MASIR_P_LowStat < MASIR_P_teststat & MASIR_P_UpStat > MASIR_P_teststat)])

MASIR_P_LB
MASIR_P_UB
Optp
MAOptp

# MASIR Init100
Opt_MBSIR
MBGamma_fit
# beta derivative
MBSIR_Opt_plusBeta <- MASIR_Proc(MBOptbeta+d/2,MBGamma_fit,init_S = (N-I0_MB)/N, ODEmaxTime = fitend_tps, ODEstep = 1e-1,TrackDyn = T)
MBSIR_R_Opt_plusBeta <- MBSIR_Opt_plusBeta$Dynamic[fitend_timeseq,4]*N
MBSIR_Opt_minusBeta <- MASIR_Proc(MBOptbeta-d/2,MBGamma_fit,init_S = (N-I0_MB)/N, ODEmaxTime = fitend_tps, ODEstep = 1e-1,TrackDyn = T)
MBSIR_R_Opt_minusBeta <- MBSIR_Opt_minusBeta$Dynamic[fitend_timeseq,4]*N

# P derivative
################### Fitting gamma based on p
MBOptp
MBSIR_p_vec <- c(MBOptp-d/2,MBOptp,MBOptp+d/2)
MBSIR_Gamma_vec <- c(1:length(MBSIR_p_vec))
MBSIR_p_vec
for (i in c(1:length(MBSIR_p_vec))) {
  p <- MBSIR_p_vec[i]
  # Reported is just p proportion of cases
  # All cases not reported will develop to Late Latent stage and become "recovered"
  StageN <- sum(stageC)/p
  StageDC <- stageC[5]/max(stageT)*(right-left)+stageC[1:4]
  StageDP <- c(1:3)
  
  for (j in c(1:3)) {
    StageDP[j] <- sum(StageDC[1:j])/StageN
  }
  StageDQ <- left[2:4]
  
  #### Adding an End Constraint as the third percentile
  StageDP[3] <- 1
  StageDQ[3] <- 26
  
  
  fitexp <- get.eexp.par(StageDP,StageDQ, fit.weight=c(1,1,1),tol = 0.03,show.output = FALSE,plot = F)
  
  if (is.null(fitexp$par) == TRUE){
    gamma <- 0
  } else {
    gamma <- fitexp$par
  }
  MBSIR_Gamma_vec[i] <- gamma
}

MBSIR_Opt_plusP <- MASIR_Proc(MBOptbeta,MBSIR_Gamma_vec[3],init_S = (N-I0_MB)/N, ODEmaxTime = fitend_tps, ODEstep = 1e-1,TrackDyn = T)
MBSIR_R_Opt_plusP <- MBSIR_Opt_plusP$Dynamic[fitend_timeseq,4]*N
MBSIR_Opt_minusP <- MASIR_Proc(MBOptbeta,MBSIR_Gamma_vec[1],init_S = (N-I0_MB)/N, ODEmaxTime = fitend_tps, ODEstep = 1e-1,TrackDyn = T)
MBSIR_R_Opt_minusP <- MBSIR_Opt_minusP$Dynamic[fitend_timeseq,4]*N

# Jacobian Jf
MBSIR_loglik_Opt_plusBeta <- rep(0,fitend_tps)
MBSIR_loglik_Opt_minusBeta <- rep(0,fitend_tps)
MBSIR_loglik_Opt_plusP <- rep(0,fitend_tps)
MBSIR_loglik_Opt_minusP <- rep(0,fitend_tps)

for (i in c(1:fitend_tps)){
  MBSIR_loglik_Opt_plusBeta[i] <- Den_cb_est(MBSIR_R_Opt_plusBeta[i+1]-MBSIR_R_Opt_plusBeta[i], R_datafit[i], MBSIR_p_vec[2])
  MBSIR_loglik_Opt_minusBeta[i] <- Den_cb_est(MBSIR_R_Opt_minusBeta[i+1]-MBSIR_R_Opt_minusBeta[i], R_datafit[i], MBSIR_p_vec[2])
  MBSIR_loglik_Opt_plusP[i] <- Den_cb_est(MBSIR_R_Opt_plusP[i+1]-MBSIR_R_Opt_plusP[i], R_datafit[i], MBSIR_p_vec[3])
  MBSIR_loglik_Opt_minusP[i] <- Den_cb_est(MBSIR_R_Opt_minusP[i+1]-MBSIR_R_Opt_minusP[i], R_datafit[i], MBSIR_p_vec[1])
}

MBSIR_Jf_Betarow <- (MBSIR_loglik_Opt_plusBeta-MBSIR_loglik_Opt_minusBeta)/d
MBSIR_Jf_Prow <- (MBSIR_loglik_Opt_plusP-MBSIR_loglik_Opt_minusP)/d

MBSIR_Jf <- cbind(MBSIR_Jf_Betarow,MBSIR_Jf_Prow)

MBSIR_GTGmat <- t(MBSIR_Jf) %*% MBSIR_Jf  
#MBSIR_Qmat <- 1/fitend_tps*MBSIR_Hmat

MBSIR_V0mat <- inv(MBSIR_GTGmat)

# Confidence Interval of Single parameter
MBSIR_V_Beta <- MBSIR_V0mat[1,1]
MBSIR_V_P <- MBSIR_V0mat[2,2]
MBSIR_Beta_UpStat <- qnorm(0.975, 0, sqrt(MBSIR_V_Beta))
MBSIR_Beta_LowStat <- qnorm(0.025, 0, sqrt(MBSIR_V_Beta))
MBSIR_P_UpStat <- qnorm(0.975, 0, sqrt(MBSIR_V_P))
MBSIR_P_LowStat <- qnorm(0.025, 0, sqrt(MBSIR_V_P))

# Beta
MBSIR_Beta_CI_seq <- seq(0.03,0.14,1e-5)
MBSIR_Beta_teststat <- sqrt(fitend_tps)*(MBOptbeta-MBSIR_Beta_CI_seq)
MBSIR_Beta_LB <- min(MBSIR_Beta_CI_seq[which(MBSIR_Beta_LowStat < MBSIR_Beta_teststat & MBSIR_Beta_UpStat > MBSIR_Beta_teststat)])
MBSIR_Beta_UB <- max(MBSIR_Beta_CI_seq[which(MBSIR_Beta_LowStat < MBSIR_Beta_teststat & MBSIR_Beta_UpStat > MBSIR_Beta_teststat)])

MBSIR_Beta_LB
MBSIR_Beta_UB
MBOptbeta

# P
MBSIR_P_CI_seq <- seq(0,0.5,1e-4)
MBSIR_P_teststat <- sqrt(fitend_tps)*(MBOptp-MBSIR_P_CI_seq)
MBSIR_P_LB <- min(MBSIR_P_CI_seq[which(MBSIR_P_LowStat < MBSIR_P_teststat & MBSIR_P_UpStat > MBSIR_P_teststat)])
MBSIR_P_UB <- max(MBSIR_P_CI_seq[which(MBSIR_P_LowStat < MBSIR_P_teststat & MBSIR_P_UpStat > MBSIR_P_teststat)])

MBSIR_P_LB
MBSIR_P_UB
Optp
MBOptp


######################### Sensitive Analysis of MSV

# disturbance

# 1% each side to alpha
# use single variable confidence interval for beta and p
# Only do that for I0=117 case
# Bar figure

##### Beta
Optbeta

MSV_Beta_LB
MSV_Beta_UB

CM_Opt$RInfinity
CM_Opt$RInfinity*N


##### MSV plus beta
# MSV_Opt_plusBeta
MSV_Sen_plusBeta <- ModProc_CM(DDist,MSV_Beta_UB,Gamma_fit,ODEmaxTime = 500, ODEstep = 1e-1,init_theta = it_theta,TrackDyn = T)
df_Sys_plusB <- as.data.frame(MSV_Sen_plusBeta$Dynamic)

MSV_Sen_plusBeta$RInfinity
MSV_Sen_plusBeta$RInfinity*N

##### MSV minus beta
#MSV_Opt_minusBeta
MSV_Sen_minusBeta <- ModProc_CM(DDist,MSV_Beta_LB,Gamma_fit,ODEmaxTime = 500, ODEstep = 1e-1,init_theta = it_theta,TrackDyn = T)
df_Sys_minusB <- as.data.frame(MSV_Sen_minusBeta$Dynamic)
MSV_Sen_minusBeta$RInfinity
MSV_Sen_minusBeta$RInfinity*N


# R Tendency
ggplot()+theme_bw()+
  geom_line(data = df_Sys,aes(x=t_Sys, y=R_Sys,color="Opt"))+
  geom_line(data = df_Sys_plusB,aes(x=time, y=R,color="plus"))+
  geom_line(data = df_Sys_minusB,aes(x=time, y=R,color="minus"))+
  scale_color_manual(labels=c("Optimal Parameters",bquote(~beta*"="~hat(beta)+d), bquote(~beta*"="~hat(beta)-d)),breaks=c("Opt","plus", "minus"),values=c("Opt"="black", "plus"="red","minus"="blue"))+
  labs(x="Time (bi-week)",y="R Proportion",color="Color")+
  xlim(0,500)
######### Both large and small grass
######### Geom_Shade

# I Tendency
ggplot()+theme_bw()+
  geom_line(data = df_Sys,aes(x=t_Sys, y=I_Sys,color="Opt"))+
  geom_line(data = df_Sys_plusB,aes(x=time, y=I_out,color="plus"))+
  geom_line(data = df_Sys_minusB,aes(x=time, y=I_out,color="minus"))+
  scale_color_manual(labels=c("Optimal Parameters",bquote(~beta*"="~hat(beta)+d), bquote(~beta*"="~hat(beta)-d)),breaks=c("Opt","plus", "minus"),values=c("Opt"="black", "plus"="red","minus"="blue"))+
  labs(x="Time (bi-week)",y="R Proportion",color="Color")+
  xlim(0,500)
#+ylim(0.06,0.07)


# p-gamma

##### p
Optp
MSV_P_UB
MSV_P_LB

##### gamma
Opt$gamma_value
MSV_p_sen_vec <- c(MSV_P_LB,Optp,MSV_P_UB)
MSV_Gamma_sen_vec <- c(1:length(MSV_p_sen_vec))
MSV_p_sen_vec
for (i in c(1:length(MSV_p_sen_vec))) {
  p <- MSV_p_sen_vec[i]
  # Reported is just p proportion of cases
  # All cases not reported will develop to Late Latent stage and become "recovered"
  StageN <- sum(stageC)/p
  StageDC <- stageC[5]/max(stageT)*(right-left)+stageC[1:4]
  StageDP <- c(1:3)
  
  for (j in c(1:3)) {
    StageDP[j] <- sum(StageDC[1:j])/StageN
  }
  StageDQ <- left[2:4]
  
  #### Adding an End Constraint as the third percentile
  StageDP[3] <- 1
  StageDQ[3] <- 26
  
    fitexp <- get.eexp.par(StageDP,StageDQ, fit.weight=c(1,1,1),tol = 0.02,show.output = FALSE,plot = F)
  
  if (is.null(fitexp$par) == TRUE){
    gamma <- 0
  } else {
    gamma <- fitexp$par
  }
  MSV_Gamma_sen_vec[i] <- gamma
}


##### MSV plus p
# MSV_Opt_plusp
MSV_Sen_plusP <- ModProc_CM(DDist,Optbeta,MSV_Gamma_sen_vec[3],ODEmaxTime = 500, ODEstep = 1e-1,init_theta = it_theta,TrackDyn = T)
df_Sys_plusP <- as.data.frame(MSV_Sen_plusP$Dynamic)

MSV_Sen_plusP$RInfinity
MSV_Sen_plusP$RInfinity*N

##### MSV minus p
#MSV_Opt_minusp
MSV_Sen_minusP <- ModProc_CM(DDist,Optbeta,MSV_Gamma_sen_vec[1],ODEmaxTime = 500, ODEstep = 1e-1,init_theta = it_theta,TrackDyn = T)
df_Sys_minusP <- as.data.frame(MSV_Sen_minusP$Dynamic)
MSV_Sen_minusP$RInfinity
MSV_Sen_minusP$RInfinity*N

# R Tendency
ggplot()+theme_bw()+
  geom_line(data = df_Sys,aes(x=t_Sys, y=R_Sys,color="Opt"))+
  geom_line(data = df_Sys_plusP,aes(x=time, y=R,color="plus"))+
  geom_line(data = df_Sys_minusP,aes(x=time, y=R,color="minus"))+
  scale_color_manual(labels=c("Optimal Parameters",bquote(~p*"="~hat(p)+d), bquote(~p*"="~hat(p)-d)),breaks=c("Opt","plus", "minus"),values=c("Opt"="black", "plus"="red","minus"="blue"))+
  labs(x="Time (bi-week)",y="R Proportion",color="Color")+
  xlim(0,500)
#+ylim(0.06,0.07)


# I Tendency
ggplot()+theme_bw()+
  geom_line(data = df_Sys,aes(x=t_Sys, y=I_Sys,color="Opt"))+
  geom_line(data = df_Sys_plusP,aes(x=time, y=I_out,color="plus"))+
  geom_line(data = df_Sys_minusP,aes(x=time, y=I_out,color="minus"))+
  scale_color_manual(labels=c("Optimal Parameters",bquote(~p*"="~hat(p)+d), bquote(~p*"="~hat(p)-d)),breaks=c("Opt","plus", "minus"),values=c("Opt"="black", "plus"="red","minus"="blue"))+
  labs(x="Time (bi-week)",y="R Proportion",color="Color")+
  xlim(0,500)
#+ylim(0.06,0.07)



# alpha
alpha_plfit
alpha_plus <- alpha_plfit*(1+0.02)
alpha_minus <- alpha_plfit*(1-0.02)
alpha_plus
alpha_minus

alpha_plfit


# fitted distribution pmf curves
Seq_plfit_plus <- data.frame(deg,prob)
Seq_plfit_minus <- data.frame(deg,prob)

for (i in c(1:200)) {
  Seq_plfit_plus[i,2] <- i^(-alpha_plus)
  Seq_plfit_minus[i,2] <- i^(-alpha_minus)
}

# Normalize scaling factor
scal_plus <- sum(Seq_plfit_plus$prob)
Seq_plfit_plus$Nprob <- Seq_plfit_plus$prob/scal_plus
sum(Seq_plfit_plus$Nprob)

scal_minus <- sum(Seq_plfit_minus$prob)
Seq_plfit_minus$Nprob <- Seq_plfit_minus$prob/scal_minus
sum(Seq_plfit_minus$Nprob)

# length(Deg_Seq)
## Using plfit method as degree distribution
Pk_plus <- append(Seq_plfit_plus$Nprob,0,0)
Pk_minus <- append(Seq_plfit_minus$Nprob,0,0)

Pk <- Pk_plus
DDist_plus <- data.frame(kvalue,Pk)

Pk <- Pk_minus
DDist_minus <- data.frame(kvalue,Pk)

MSV_Sen_plusAlpha <- ModProc_CM(DDist_plus,Optbeta,Gamma_fit,ODEmaxTime = 500, ODEstep = 1e-1,init_theta = it_theta,TrackDyn = T)
MSV_Sen_minusAlpha <- ModProc_CM(DDist_minus,Optbeta,Gamma_fit,ODEmaxTime = 500, ODEstep = 1e-1,init_theta = it_theta,TrackDyn = T)


##### MSV plus Alpha
# MSV_Opt_plusAlpha
df_Sys_plusA <- as.data.frame(MSV_Sen_plusAlpha$Dynamic)
MSV_Sen_plusAlpha$RInfinity
MSV_Sen_plusAlpha$RInfinity*N

##### MSV minus beta
#MSV_Opt_minusAlpha
df_Sys_minusA <- as.data.frame(MSV_Sen_minusAlpha$Dynamic)
MSV_Sen_minusAlpha$RInfinity
MSV_Sen_minusAlpha$RInfinity*N

# R Tendency
ggplot()+theme_bw()+
  geom_line(data = df_Sys,aes(x=t_Sys, y=R_Sys,color="Opt"))+
  geom_line(data = df_Sys_plusA,aes(x=time, y=R,color="plus"))+
  geom_line(data = df_Sys_minusA,aes(x=time, y=R,color="minus"))+
  scale_color_manual(labels=c("Optimal Parameters",bquote(~alpha*"="~hat(alpha)+d), bquote(~alpha*"="~hat(alpha)-d)),breaks=c("Opt","plus", "minus"),values=c("Opt"="black", "plus"="red","minus"="blue"))+
  labs(x="Time (bi-week)",y="R Proportion",color="Color")+
  xlim(0,500)

# I Tendency
ggplot()+theme_bw()+
  geom_line(data = df_Sys,aes(x=t_Sys, y=I_Sys,color="Opt"))+
  geom_line(data = df_Sys_plusA,aes(x=time, y=I_out,color="plus"))+
  geom_line(data = df_Sys_minusA,aes(x=time, y=I_out,color="minus"))+
  scale_color_manual(labels=c("Optimal Parameters",bquote(~alpha*"="~hat(alpha)+d), bquote(~alpha*"="~hat(alpha)-d)),breaks=c("Opt","plus", "minus"),values=c("Opt"="black", "plus"="red","minus"="blue"))+
  labs(x="Time (bi-week)",y="R Proportion",color="Color")+
  xlim(0,500)


##########MBSIR I0=117

##### Beta
MBOptbeta
MBSIR_Beta_LB
MBSIR_Beta_UB

MB_Opt$RInfinity
MB_Opt$RInfinity*N


##### MSV plus beta
# MSV_Opt_plusBeta
MBSIR_Sen_plusBeta <- MASIR_Proc(MBSIR_Beta_UB,MBGamma_fit,init_S = (N-I0_MB)/N, ODEmaxTime = 500, ODEstep = 1e-1,TrackDyn = T)
df_MBSIR_plusB <- as.data.frame(MBSIR_Sen_plusBeta$Dynamic)

MBSIR_Sen_plusBeta$RInfinity
MBSIR_Sen_plusBeta$RInfinity*N

##### MSV minus beta
#MSV_Opt_minusBeta
MBSIR_Sen_minusBeta <- MASIR_Proc(MBSIR_Beta_LB,MBGamma_fit,init_S = (N-I0_MB)/N, ODEmaxTime = 500, ODEstep = 1e-1,TrackDyn = T)
df_MBSIR_minusB <- as.data.frame(MBSIR_Sen_minusBeta$Dynamic)

MBSIR_Sen_minusBeta$RInfinity
MBSIR_Sen_minusBeta$RInfinity*N

# R Tendency
ggplot()+theme_bw()+
  geom_line(data = MBSys_Opt,aes(x=time, y=R,color="Opt"))+
  geom_line(data = df_MBSIR_plusB,aes(x=time, y=R,color="plus"))+
  geom_line(data = df_MBSIR_minusB,aes(x=time, y=R,color="minus"))+
  scale_color_manual(labels=c("Optimal Parameters",bquote(~beta*"="~hat(beta)+d), bquote(~beta*"="~hat(beta)-d)),breaks=c("Opt","plus", "minus"),values=c("Opt"="black", "plus"="red","minus"="blue"))+
  labs(x="Time (bi-week)",y="R Proportion",color="Color")+
  xlim(0,500)
######### Both large and small grass
######### Geom_Shade

# I Tendency
ggplot()+theme_bw()+
  geom_line(data = MBSys_Opt,aes(x=time, y=I,color="Opt"))+
  geom_line(data = df_MBSIR_plusB,aes(x=time, y=I,color="plus"))+
  geom_line(data = df_MBSIR_minusB,aes(x=time, y=I,color="minus"))+
  scale_color_manual(labels=c("Optimal Parameters",bquote(~beta*"="~hat(beta)+d), bquote(~beta*"="~hat(beta)-d)),breaks=c("Opt","plus", "minus"),values=c("Opt"="black", "plus"="red","minus"="blue"))+
  labs(x="Time (bi-week)",y="R Proportion",color="Color")+
  xlim(0,500)
#+ylim(0.06,0.07)


# p-gamma

##### p
MBOptp
MBSIR_P_UB
MBSIR_P_LB

##### gamma
Opt_MBSIR$gamma_value_MBSIR
MBSIR_p_sen_vec <- c(MBSIR_P_LB,MBOptp,MBSIR_P_UB)
MBSIR_Gamma_sen_vec <- c(1:length(MBSIR_p_sen_vec))
MBSIR_p_sen_vec

for (i in c(1:length(MBSIR_p_sen_vec))) {
  p <- MBSIR_p_sen_vec[i]
  # Reported is just p proportion of cases
  # All cases not reported will develop to Late Latent stage and become "recovered"
  StageN <- sum(stageC)/p
  StageDC <- stageC[5]/max(stageT)*(right-left)+stageC[1:4]
  StageDP <- c(1:3)
  
  for (j in c(1:3)) {
    StageDP[j] <- sum(StageDC[1:j])/StageN
  }
  StageDQ <- left[2:4]
  
  #### Adding an End Constraint as the third percentile
  StageDP[3] <- 1
  StageDQ[3] <- 26
  
  fitexp <- get.eexp.par(StageDP,StageDQ, fit.weight=c(1,1,1),tol = 0.03,show.output = FALSE,plot = F)
  
  if (is.null(fitexp$par) == TRUE){
    gamma <- 0
  } else {
    gamma <- fitexp$par
  }
  MBSIR_Gamma_sen_vec[i] <- gamma
}

MBSIR_Gamma_sen_vec

##### MBSIR plus p
MBSIR_Sen_plusP <- MASIR_Proc(MBOptbeta,MBSIR_Gamma_sen_vec[3],init_S = (N-I0_MB)/N, ODEmaxTime = 500, ODEstep = 1e-1,TrackDyn = T)
df_MBSIR_plusP <- as.data.frame(MBSIR_Sen_plusP$Dynamic)

MBSIR_Sen_plusP$RInfinity
MBSIR_Sen_plusP$RInfinity*N

##### MBSIR minus p
MBSIR_Sen_minusP <- MASIR_Proc(MBOptbeta,MBSIR_Gamma_sen_vec[1],init_S = (N-I0_MB)/N, ODEmaxTime = 500, ODEstep = 1e-1,TrackDyn = T)
df_MBSIR_minusP <- as.data.frame(MBSIR_Sen_minusP$Dynamic)

MBSIR_Sen_minusP$RInfinity
MBSIR_Sen_minusP$RInfinity*N

# R Tendency
ggplot()+theme_bw()+
  geom_line(data = MBSys_Opt,aes(x=time, y=R,color="Opt"))+
  geom_line(data = df_MBSIR_plusP,aes(x=time, y=R,color="plus"))+
  geom_line(data = df_MBSIR_minusP,aes(x=time, y=R,color="minus"))+
  scale_color_manual(labels=c("Optimal Parameters",bquote(~p*"="~hat(p)+d), bquote(~p*"="~hat(p)-d)),breaks=c("Opt","plus", "minus"),values=c("Opt"="black", "plus"="red","minus"="blue"))+
  labs(x="Time (bi-week)",y="R Proportion",color="Color")+
  xlim(0,500)
######### Both large and small grass
######### Geom_Shade

# I Tendency
ggplot()+theme_bw()+
  geom_line(data = MBSys_Opt,aes(x=time, y=I,color="Opt"))+
  geom_line(data = df_MBSIR_plusP,aes(x=time, y=I,color="plus"))+
  geom_line(data = df_MBSIR_minusP,aes(x=time, y=I,color="minus"))+
  scale_color_manual(labels=c("Optimal Parameters",bquote(~p*"="~hat(p)+d), bquote(~p*"="~hat(p)-d)),breaks=c("Opt","plus", "minus"),values=c("Opt"="black", "plus"="red","minus"="blue"))+
  labs(x="Time (bi-week)",y="R Proportion",color="Color")+
  xlim(0,500)
#+ylim(0.06,0.07)


##########MASIR I0=27

##### Beta
MAOptbeta
MASIR_Beta_LB
MASIR_Beta_UB

MA_Opt$RInfinity
MA_Opt$RInfinity*N


##### MSV plus beta
# MSV_Opt_plusBeta
MASIR_Sen_plusBeta <- MASIR_Proc(MASIR_Beta_UB,MAGamma_fit,init_S = (N-I0_MA)/N, ODEmaxTime = 500, ODEstep = 1e-1,TrackDyn = T)
df_MASIR_plusB <- as.data.frame(MASIR_Sen_plusBeta$Dynamic)

MASIR_Sen_plusBeta$RInfinity
MASIR_Sen_plusBeta$RInfinity*N

##### MSV minus beta
#MSV_Opt_minusBeta
MASIR_Sen_minusBeta <- MASIR_Proc(MASIR_Beta_LB,MAGamma_fit,init_S = (N-I0_MA)/N, ODEmaxTime = 500, ODEstep = 1e-1,TrackDyn = T)
df_MASIR_minusB <- as.data.frame(MASIR_Sen_minusBeta$Dynamic)

MASIR_Sen_minusBeta$RInfinity
MASIR_Sen_minusBeta$RInfinity*N

# p-gamma

##### p
MAOptp
MASIR_P_UB
MASIR_P_LB

##### gamma
Opt_MASIR$gamma_value_MASIR
MASIR_p_sen_vec <- c(MASIR_P_LB,MAOptp,MASIR_P_UB)
MASIR_Gamma_sen_vec <- c(1:length(MASIR_p_sen_vec))
MASIR_p_sen_vec

for (i in c(1:length(MASIR_p_sen_vec))) {
  p <- MASIR_p_sen_vec[i]
  # Reported is just p proportion of cases
  # All cases not reported will develop to Late Latent stage and become "recovered"
  StageN <- sum(stageC)/p
  StageDC <- stageC[5]/max(stageT)*(right-left)+stageC[1:4]
  StageDP <- c(1:3)
  
  for (j in c(1:3)) {
    StageDP[j] <- sum(StageDC[1:j])/StageN
  }
  StageDQ <- left[2:4]
  
  #### Adding an End Constraint as the third percentile
  StageDP[3] <- 1
  StageDQ[3] <- 26
  
  fitexp <- get.eexp.par(StageDP,StageDQ, fit.weight=c(1,1,1),tol = 0.03,show.output = FALSE,plot = F)
  
  if (is.null(fitexp$par) == TRUE){
    gamma <- 0
  } else {
    gamma <- fitexp$par
  }
  MASIR_Gamma_sen_vec[i] <- gamma
}

MASIR_Gamma_sen_vec

##### MASIR plus p
MASIR_Sen_plusP <- MASIR_Proc(MAOptbeta,MASIR_Gamma_sen_vec[3],init_S = (N-I0_MA)/N, ODEmaxTime = 500, ODEstep = 1e-1,TrackDyn = T)
df_MASIR_plusP <- as.data.frame(MASIR_Sen_plusP$Dynamic)

MASIR_Sen_plusP$RInfinity
MASIR_Sen_plusP$RInfinity*N

##### MASIR minus p
MASIR_Sen_minusP <- MASIR_Proc(MAOptbeta,MASIR_Gamma_sen_vec[1],init_S = (N-I0_MA)/N, ODEmaxTime = 500, ODEstep = 1e-1,TrackDyn = T)
df_MASIR_minusP <- as.data.frame(MASIR_Sen_minusP$Dynamic)

MASIR_Sen_minusP$RInfinity
MASIR_Sen_minusP$RInfinity*N



library(forcats)
library(ggbreak)

MSV_base <- round(CM_Opt$RInfinity*N,0)
MSV_plusB <- round(MSV_Sen_plusBeta$RInfinity*N,0)
MSV_minusB <- round(MSV_Sen_minusBeta$RInfinity*N,0)
MSV_plusP <- round(MSV_Sen_plusP$RInfinity*N,0)
MSV_minusP <- round(MSV_Sen_minusP$RInfinity*N,0)
MSV_plusA <- round(MSV_Sen_plusAlpha$RInfinity*N,0)
MSV_minusA <- round(MSV_Sen_minusAlpha$RInfinity*N,0)
MBSIR_base <- round(MB_Opt$RInfinity*N,0)
MBSIR_plusB <- round(MBSIR_Sen_plusBeta$RInfinity*N,0)
MBSIR_minusB <- round(MBSIR_Sen_minusBeta$RInfinity*N,0)
MBSIR_plusP <- round(MBSIR_Sen_plusP$RInfinity*N,0)
MBSIR_minusP <- round(MBSIR_Sen_minusP$RInfinity*N,0)
MASIR_base <- round(MA_Opt$RInfinity*N,0)
MASIR_plusB <- round(MASIR_Sen_plusBeta$RInfinity*N,0)
MASIR_minusB <- round(MASIR_Sen_minusBeta$RInfinity*N,0)
MASIR_plusP <- round(MASIR_Sen_plusP$RInfinity*N,0)
MASIR_minusP <- round(MASIR_Sen_minusP$RInfinity*N,0)



Variation <- c("MSV beta","MSV p","MSV alpha","MASIR I0=27 beta","MASIR I0=27 p","MASIR I0=117 beta","MASIR I0=117 p")
y <-    c(         1,  2,  3,  4,  5, 6, 7)
tmin <- c(MSV_minusB, MSV_plusP,  MSV_plusA, MASIR_minusB,  MASIR_plusP, MBSIR_minusB,  MBSIR_plusP) 
tmax <- c( MSV_plusB,MSV_minusP, MSV_minusA,  MASIR_plusB, MASIR_minusP,  MBSIR_plusB, MBSIR_minusP)
tpoint <-c( MSV_base,  MSV_base,   MSV_base,   MASIR_base,   MASIR_base,   MBSIR_base,   MBSIR_base)
Model <- c(    "MSV",     "MSV",      "MSV",      "MASIR",      "MASIR",      "MBSIR",      "MBSIR")

df <- as.data.frame(cbind(y,tmin,tmax,tpoint))
df$stage <- Variation
df$Model <- Model
df

ggplot(data=df,aes(y=fct_inorder(stage)))+
  theme_classic()+
  #geom_point(aes(x=tpoint,color=Model),size=1.8,shape=2)+
  geom_text(aes(x=tmax,label=tmax,color=Model),size=2, hjust=-0.4,vjust=-0.9, show.legend = FALSE)+
  geom_text(aes(x=tmin,label=tmin,color=Model),size=2, hjust=1.2,vjust=-0.9, show.legend = FALSE)+
  geom_errorbarh(aes(xmin=tmin,xmax=tmax,color=Model),height=0.4)+
  labs(x="Final Infected Size",y="Parameter Being Varied",color="Model")+
  scale_x_break(c(1900,14000),scales=1.5)+
  expand_limits(x = c(1550,19000))+
  geom_vline(xintercept=MSV_base,linewidth=0.5,linetype="dashed")+
  geom_vline(xintercept=MBSIR_base,linewidth=0.5,linetype="dashed",color="red")+
  geom_vline(xintercept=MASIR_base,linewidth=0.5,linetype="dashed",color="blue")+
  scale_x_continuous(breaks=c(1400,1500,1600,1718,1800,1900,14000, 15171, 16000,17000,17501,18000,19000), sec.axis = sec_axis(~ . /N, name = "Final Infected Proportion"))+
  scale_y_discrete(labels=c(bquote(~beta),bquote(~p), bquote(~alpha),bquote(~beta),bquote(~p),bquote(~beta),bquote(~p)))+
  scale_color_manual(labels=c("Network Model", bquote("MA"~I[0]*"= 27"), bquote("MA"~I[0]*"= 117")),breaks=c("MSV", "MASIR", "MBSIR"),values=c("MASIR"="blue", "MBSIR"="red", "MSV"="black"))

MSV_base
MASIR_base
MBSIR_base
