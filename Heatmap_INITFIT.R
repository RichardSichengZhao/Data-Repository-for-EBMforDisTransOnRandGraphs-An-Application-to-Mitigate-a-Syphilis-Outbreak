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
  ModProc_CM <- function(Pk, beta, gamma, init_theta=1e-3, ODEmaxTime=50, ODEstep=1e-2,ThetaTol=1e-9, TrackDyn=T){
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
      y <- c(theta=1-init_theta,R=0)
      
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
    step_theta <- 1e-2
    
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
  # alpha_plfit
  # alpha_Init
  
  
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

}

############## Fitted Result HPC
# Test25_Init10 <- read.csv(paste0(wd,"/Test25_Init10_combined.csv"), header=T)
#Test100_Init10 <- read.csv(paste0(wd,"/Test100_Init10_combined.csv"), header=T)
#Test100_Init10_2QT <- read.csv(paste0(wd,"/2QT_Test100_Init10_combined.csv"), header=T)
#Test100_Init5_2QT <- read.csv(paste0(wd,"/2QT_Test100_Init5_combined.csv"), header=T)
#Test100_Init5 <- read.csv(paste0(wd,"/Test100_Init5_combined.csv"), header=T)
#Test100_Init10_3QT <-read.csv(paste0(wd,"/3QT_Test100_Init10_combined.csv"), header=T)
#Test100_Init10_3QT_MASIR <-read.csv(paste0(wd,"/3QT_Test100_Init10_MASIR_combined.csv"), header=T)
Test100_InitFIT_3QT <-read.csv(paste0(wd,"/3QT_Test100_InitFIT_combined.csv"), header=T)
Test100_InitFIT_3QT_MASIR <-read.csv(paste0(wd,"/3QT_Test100_InitFIT_MASIR_combined.csv"), header=T)

#readmat<-Test25_Init10
readmat <- Test100_InitFIT_3QT

readmat[1,]
length(readmat[,1])

beta_value <- readmat[,4]
gamma_value <- readmat[,3]
p_value <- readmat[,2]
I0_value <- readmat[,5]
lik_value <- readmat[,6]
# log_diff <- -log(max(lik_value)-lik_value+0.1)

df <- data.frame(beta_value,gamma_value,p_value,I0_value,lik_value)

ListInit <- list(0)
ListMaxInd <- list(0)
MaxLik_seq <- c(0)
for (i in c(1:max(I0_value))) {
  dfi <- subset(df,df$I0_value==i)
  max_ind <- which.max(dfi$lik_value)
  lik_vec <- dfi$lik_value
  ListMaxInd[[i]] <- dfi[max_ind,]
  MaxLik_seq[i] <- dfi[max_ind,5]
  log_diff <- -log(max(lik_vec)-lik_vec+0.1)
  ListInit[[i]] <- cbind(dfi,log_diff)
}
ListMaxInd
MaxLik_seq

Init_seq <- c(1:max(I0_value))
df_Init <- cbind(Init_seq,MaxLik_seq)
df_Init
ggplot(df_Init, aes(x=Init_seq,y=MaxLik_seq))+theme_bw()+
  geom_point()+
  annotate("point",x=which.max(MaxLik_seq),y=max(MaxLik_seq),colour="red",size=3)+
  ylim(-200,-175)+
  xlim(10,50)+
  labs(x="Initial Active Infected Casecount "~I[0], y="Maxium Log-likelihood")


# plot(c(1:max(I0_value)),MaxLik_seq)
# points(14, MaxLik_seq[14],col="red",pty=2)
#dfs <- subset(df,df$beta_value<0.15&df$beta_value>0.05)


max(lik_value)
which.max(lik_value)
Opt <- df[which.max(lik_value),]
Opt
# Optimal initial condition I0=14

ggplot(ListInit[[27]], aes(x = beta_value, y = p_value, fill = log_diff))+
  geom_tile() +
  scale_fill_gradient2(low = "white",mid="orange", midpoint = -2, high = "red")

Opt
# ggplot(df, aes(x = beta_value, y = p_value, fill = lik_value))+
#   geom_tile(color = "black") +
#   scale_fill_gradient(low = "white", high = "blue")

# ggplot(df, aes(x = beta_value, y = p_value, fill = gamma_value))+
#   geom_tile(color = "black") +
#   scale_fill_gradient(low = "white", high = "blue")




Opt
### Optimal Parameter Case
Optbeta <- Opt[1,1]
# beta = 0.0037 bi-weekly
Gamma_fit <- Opt[1,2]
# gamma=0.005329582 bi-weekly
Optp <- Opt[1,3]
# p = 0.14
it_theta <- Init_theta_func(Opt[1,4])
# S0Count <- PGFG0(1-it_theta,DDist)*N
# S0Count

Gamma_fit
(1/Gamma_fit)/26

# Checking Error
CM_Opt<- ModProc_CM(DDist,Optbeta,Gamma_fit,ODEmaxTime = tps, ODEstep = 1e-1,init_theta = it_theta,TrackDyn = T)
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
plot(seq(1,tps), inc_count, type='l', col="red", ylim = c(0,20),xlab ="t", ylab='Incident R')
points(seq(1,tps), inc_count*Optp, type='l', col="green")
abline(v=116)
points(seq(1,fitend_tps),R_datafit)
points(seq(fitend_tps+1,tps), R_data[seq(fitend_tps+1,tps)], col='blue')
legend('topleft',legend=c('Prediction Curve', 'Expected Reporting Curve','fitted data','unfitted data'), col=c('red','green','black','blue'), lty=c(1,1,0,0), pch = c(NA,NA,1,1),cex=0.5)

# Cumulative Compare
#R_data
CR_data<-R_data
for (i in c(1:length(R_data))){
  CR_data[i] <- sum(R_data[1:i])
}

plot(Sys_Opt[,1],Sys_Opt[,3]*N, col="red", ylim = c(0,400), xlab ="t (bi-week)", ylab='Cumulative R')
abline(v=116)
points(Sys_Opt[,1], Sys_Opt[,3]*N*Optp, col="green")
points(seq(1,fitend_tps),CR_data[seq(1,fitend_tps)], col='black', pch=1)
points(seq(fitend_tps+1,tps), CR_data[seq(fitend_tps+1,tps)], col='blue',pch=1)
legend('topleft',legend=c('Prediction Curve','Expected Reporting Curve','fitted data','unfitted data'), col=c('red','green','black','blue'), pch = c(19,19,1,1),cex=0.5)

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



# Test100_InitFIT_3QT_MASIR125 <-read.csv(paste0(wd,"/3QT_Test100_InitFIT_MASIR_combined_1-125.csv"), header=T)
#readmat_MASIR <- Test100_InitFIT_3QT_MASIR125
readmat_MASIR <- Test100_InitFIT_3QT_MASIR
readmat_MASIR[1,]

beta_value_MASIR <- readmat_MASIR[,4]
gamma_value_MASIR <- readmat_MASIR[,3]
p_value_MASIR <- readmat_MASIR[,2]
I0_value_MASIR <- readmat_MASIR[,5]
lik_value_MASIR <- readmat_MASIR[,6]
#log_diff_MASIR <- -log(max(lik_value_MASIR)-lik_value_MASIR+0.01)

df_MASIR <- data.frame(beta_value_MASIR,gamma_value_MASIR,p_value_MASIR,I0_value_MASIR,lik_value_MASIR)
#dfs_MASIR <- subset(df_MASIR,df_MASIR$beta_value_MASIR<0.15&df_MASIR$beta_value_MASIR>0.05)

# ggplot(dfs_MASIR, aes(x = beta_value_MASIR, y = p_value_MASIR, fill = log_diff_MASIR))+
#   geom_tile() +
#   scale_fill_gradient2(low = "white",mid="orange", midpoint = -1, high = "red")

max(lik_value_MASIR)
which.max(lik_value_MASIR)
Opt_MASIR <- df_MASIR[which.max(lik_value_MASIR),]

ListInit_MASIR <- list(0)
ListMaxInd_MASIR <- list(0)
MaxLik_seq_MASIR <- c(0)
I0_vec_MASIR <- unique(I0_value_MASIR)

for (i in c(1:length(I0_vec_MASIR))) {
  dfi <- subset(df_MASIR,df_MASIR$I0_value_MASIR==I0_vec_MASIR[i])
  max_ind <- which.max(dfi$lik_value_MASIR)
  lik_vec <- dfi$lik_value_MASIR
  ListMaxInd_MASIR[[i]] <- dfi[max_ind,]
  MaxLik_seq_MASIR[i] <- dfi[max_ind,5]
  log_diff_MASIR <- -log(max(lik_vec)-lik_vec+0.1)
  ListInit_MASIR[[i]] <- cbind(dfi,log_diff_MASIR)
}
#ListMaxInd_MASIR[[14]]
#ListMaxInd_MASIR[[105]]


MaxLik_seq_MASIR
which.max(MaxLik_seq_MASIR)

Init_seq_MASIR <- I0_vec_MASIR
df_Init_MASIR <- cbind(Init_seq_MASIR,MaxLik_seq_MASIR)
df_Init_MASIR
ggplot(df_Init_MASIR, aes(x=Init_seq_MASIR,y=MaxLik_seq_MASIR))+
  theme_bw()+
  geom_point()+
  annotate("point",x=I0_vec_MASIR[which.max(MaxLik_seq_MASIR)],y=max(MaxLik_seq_MASIR),colour="red",size=3)+
  annotate("point",x=27,y=ListMaxInd_MASIR[[27]]$lik_value_MASIR,colour="green",size=3)+
  ylim(-200,-175)+
  labs(x="Initial Active Infected Casecount", y="Maxium Log-likelihood")

ggplot(ListInit_MASIR[[27]], aes(x = beta_value_MASIR, y = p_value_MASIR, fill = log_diff_MASIR))+
  geom_tile() +
  scale_fill_gradient2(low = "white",mid="orange", midpoint = -2, high = "red")
ggplot(ListInit_MASIR[[117]], aes(x = beta_value_MASIR, y = p_value_MASIR, fill = log_diff_MASIR))+
  geom_tile() +
  scale_fill_gradient2(low = "white",mid="orange", midpoint = -2, high = "red")
#ListInit_MASIR


Opt_MASIR
ListMaxInd_MASIR[[117]]
ListMaxInd_MASIR[[27]]

#Opt_MASIR <- ListMaxInd_MASIR[[117]]
Opt_MASIR <- ListMaxInd_MASIR[[27]]

Opt_MASIR
### Optimal Parameter Case
MAOptbeta <- Opt_MASIR[1,1]
# beta = 0.0037 bi-weekly
MAGamma_fit <- Opt_MASIR[1,2]
# gamma=0.005329582 bi-weekly
MAOptp <- Opt_MASIR[1,3]
# p = 0.14

MAGamma_fit
(1/MAGamma_fit)/26
MAOptbeta

# Checking Error
MA_Opt<- MASIR_Proc(MAOptbeta, MAGamma_fit, init_S = (N-Opt_MASIR[1,4])/N, ODEmaxTime=1000, ODEstep=1e-1,TrackDyn = T)
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
Opt_MASIR
# Incident Compare
plot(seq(1,tps), inc_count, type='l', col="red", ylim = c(0,50),xlab ="t", ylab='Incident R')
lines(seq(1,tps),inc_count*Optp, type = 'l',lty=2, col="red")
lines(seq(1,tps), MAinc_count, type = 'l', col='blue')
lines(seq(1,tps), MAinc_count*MAOptp, type = 'l',lty=2, col='blue')
abline(v=116)
points(seq(1,fitend_tps),R_datafit)
points(seq(fitend_tps+1,tps), R_data[seq(fitend_tps+1,tps)], col='green')
legend('topleft',legend=c('MSV Prediction Curve', 'MA Prediction Curve','fitted data','unfitted data'), col=c('red','blue','black','green'), lty=c(1,1,0,0), pch = c(NA,NA,1,1),cex=0.5)
Sys_Opt
# Cumulative Compare
plot(Sys_Opt[,1],Sys_Opt[,3]*N, col="red", type = 'l',ylim = c(0,500), xlab ="t (bi-week)", ylab='Cumulative R')
lines(MASys_Opt[,1], MASys_Opt[,4]*N, type = 'l', col='blue')
lines(Sys_Opt[,1],Sys_Opt[,3]*N*Optp, type = 'l',lty=2, col='red')
lines(MASys_Opt[,1], MASys_Opt[,4]*N*MAOptp, type = 'l',lty=2, col='blue')
abline(v=116)
points(seq(1,fitend_tps),CR_data[seq(1,fitend_tps)], col='black', pch=1)
points(seq(fitend_tps+1,tps), CR_data[seq(fitend_tps+1,tps)], col='green',pch=1)
legend('topleft',legend=c('MSV Prediction Curve', 'MA Prediction Curve','fitted data','unfitted data'), col=c('red','blue','black','green'), lty=c(1,1,0,0), pch = c(NA,NA,1,1),cex=0.5)
MAOptp
Optp
