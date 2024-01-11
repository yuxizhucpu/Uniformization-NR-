
getallinten_second_map <- function(params,dt){
  
  #linear algebra X*Beta, get each intensity
  paramatrix=matrix(params,1,6,byrow = T)
  instantinten <- dt[,1] %*% paramatrix
  
  #for 3 * 3 model, 9 columns are needed.
  q11=-(instantinten[,1]+instantinten[,2])
  q12=instantinten[,1]
  q13=instantinten[,2]
  
  q21=instantinten[,3]
  q22=-instantinten[,3]-instantinten[,4]
  q23=instantinten[,4]
  
  q31=instantinten[,5]
  q32= instantinten[,6]
  q33=-instantinten[,5]-instantinten[,6]
  
  instantinten=cbind(instantinten,q11,q22,q33)
  
  #get max(Q)
  lambda_max=apply(abs(instantinten),1,max)
  
  #get jump probability matrix
  p11=1+(q11/lambda_max)
  p12=q12/lambda_max
  p13=q13/lambda_max
  
  p21=q21/lambda_max
  p22=1+(q22/lambda_max)
  p23=q23/lambda_max
  
  p31=q31/lambda_max
  p32=q32/lambda_max
  p33=1+q33/lambda_max
  
  dt1 <- dt%>% mutate(keyformap = cur_group_id(.,time , trans))
  
  jump_prob_row=cbind(p11,p12,p13,p21,p22,p23,p31,p32,p33,q12,q13,q23,lambda_max,dt[,-1],dt1$keyformap)
  return(jump_prob_row)
  
}

#Format data into the form needed by Uniformization
getallinten_second <- function(params,dt){
  
  #linear algebra X*Beta, get each intensity
  paramatrix=matrix(params,1,6,byrow = T)
  instantinten <- dt[,1] %*% paramatrix
  
  #for 3 * 3 model, 9 columns are needed.
  q11=-(instantinten[,1]+instantinten[,2])
  q12=instantinten[,1]
  q13=instantinten[,2]
  
  q21=instantinten[,3]
  q22=-instantinten[,3]-instantinten[,4]
  q23=instantinten[,4]
  
  q31=instantinten[,5]
  q32= instantinten[,6]
  q33=-instantinten[,5]-instantinten[,6]
  
  instantinten=cbind(instantinten,q11,q22,q33)
  
  #get max(Q)
  lambda_max=apply(abs(instantinten),1,max)
  
  #get jump probability matrix
  p11=1+(q11/lambda_max)
  p12=q12/lambda_max
  p13=q13/lambda_max
  
  p21=q21/lambda_max
  p22=1+(q22/lambda_max)
  p23=q23/lambda_max
  
  p31=q31/lambda_max
  p32=q32/lambda_max
  p33=1+q33/lambda_max
  
  jump_prob_row=cbind(p11,p12,p13,p21,p22,p23,p31,p32,p33,q12,q13,q23,lambda_max,dt[,-1])
  return(jump_prob_row)
  
}

#Format data into the form needed by Ward, Canonical decomposition
getallinten_ward_second <- function(params,dt){
  
  #linear algebra X*Beta, get each intensity
  paramatrix=matrix(params,1,6,byrow = T)
  instantinten <- dt[,1] %*% paramatrix
  
  #for 3 * 3 model, 9 columns are needed.
  q11=-(instantinten[,1]+instantinten[,2])
  q12=instantinten[,1]
  q13=instantinten[,2]
  
  q21=instantinten[,3]
  q22=-instantinten[,3]-instantinten[,4]
  q23=instantinten[,4]
  
  q31=instantinten[,5]
  q32= instantinten[,6]
  q33=-instantinten[,5]-instantinten[,6]
  
  instantinten=cbind(instantinten,q11,q22)
  
  jump_prob_row=cbind(q11,q12,q13,q21,q22,q23,q31,q32,q33,dt[,-1])
  return(jump_prob_row)
  
}

#Get intensity matrix
Qmatrix <- function(a,b,c,d,e,f){
  Q=matrix(c(-(a+b),a,b,
             c,-(c+d),d,
             e,f,-(e+f)),nrow=3,ncol=3,byrow = T)
  return(Q)
}
#Transfer data from long to wide
msm_long_to_wide <- function(df){
  df1=df[1:(dim(df)[1]-1),c(1,2,3)]
  colnames(df1) <- c("subject","start","from")
  df2=df[2:(dim(df)[1]),c(2,3)]
  colnames(df2) <- c("stop","to")
  df_final=cbind(df1,df2)
  return(df_final)
}
#Function for simulation study
generate_dt <- function(x){
  store_all_data=as.list(NULL)
  
  jj=1
  while (jj<=x) {
    
    sim.df <- data.frame(subject = rep(1:500, rep(7,500)), time = rep(sort(sample(0:120,7)), 500))
    qmatrix <- rbind(c(-0.03, 0.02, 0.01 ),
                     c(0.01, -0.06, 0.05 ),
                     c(0.01, 0.03, -0.04))
    sim.df1=simmulti.msm(sim.df, qmatrix)
    
    sim.df1.wide=ddply(sim.df1,.(subject),msm_long_to_wide)
    sim.df1.wide$time=sim.df1.wide$stop-sim.df1.wide$start
    sim.df1.wide$trans=1
    sim.df1.wide$trans[sim.df1.wide$from==1&sim.df1.wide$to==1]=1
    sim.df1.wide$trans[sim.df1.wide$from==1&sim.df1.wide$to==2]=2
    sim.df1.wide$trans[sim.df1.wide$from==1&sim.df1.wide$to==3]=3
    sim.df1.wide$trans[sim.df1.wide$from==2&sim.df1.wide$to==1]=4
    sim.df1.wide$trans[sim.df1.wide$from==2&sim.df1.wide$to==2]=5
    sim.df1.wide$trans[sim.df1.wide$from==2&sim.df1.wide$to==3]=6
    sim.df1.wide$trans[sim.df1.wide$from==3&sim.df1.wide$to==1]=7
    sim.df1.wide$trans[sim.df1.wide$from==3&sim.df1.wide$to==2]=8
    sim.df1.wide$trans[sim.df1.wide$from==3&sim.df1.wide$to==3]=9
    
    dataset1=cbind(1,sim.df1.wide[,c(6,7)])
    
    params1=c(0.02,0.01,0.01,0.05,0.01,0.03)
    
    jump_matrix=getallinten_second(params1,dataset1)
    
    test=as.matrix(jump_matrix )
    
    es_res=likelihood_part2(as.matrix(jump_matrix ))
    
    socrevec=es_res[1:6]
    
    m_upper=es_res[7:27]
    l <- length(m_upper)
    n <- length(socrevec)
    # Reconstruct 
    m2 <- matrix(NA, n, n)
    m2[lower.tri(m2, diag = TRUE)] <- m_upper
    m2=t(m2)
    m2[lower.tri(m2)] <- t(m2)[lower.tri(m2)]  # If symmetric, fill also upper half
    infmatrix=m2
    
    if (sum(eigen(-m2)$value<=0)==0){
      store_all_data[[jj]]=sim.df1
      jj=jj+1
    }
    #print(jj)
  }
  
  return(store_all_data)
  
  
}

#Likelihood function based on Ward approximation
glike_ward=function(y,dataset2){
  
  jump_matrix=getallinten_ward_second(y,dataset2)
  test=as.matrix(jump_matrix)
  
  return(likelihood_ward_no_map(test))
}

#Likelihood function based on Canonical decomposition
glike_decomp=function(y,dataset2){
  jump_matrix=getallinten_ward_second(y,dataset2)
  test=as.matrix(jump_matrix )
  return(likelihood_decompno_map(test))
}

#Likelihood function based on Uniformization
glike=function(y,dataset2){
  
  jump_matrix=getallinten_second(y,dataset2)
  test=as.matrix(jump_matrix )
  
  return(likelihoodno_map(test))
}

#Optimization based on Ward approximation (N-M)
like_each_ward=function(trial){
  
  sim.df1=store_all_data1[[trial]]
  sim.df1.wide=ddply(sim.df1,.(subject),msm_long_to_wide)
  sim.df1.wide$time=sim.df1.wide$stop-sim.df1.wide$start
  sim.df1.wide$trans=1
  sim.df1.wide$trans[sim.df1.wide$from==1&sim.df1.wide$to==1]=1
  sim.df1.wide$trans[sim.df1.wide$from==1&sim.df1.wide$to==2]=2
  sim.df1.wide$trans[sim.df1.wide$from==1&sim.df1.wide$to==3]=3
  sim.df1.wide$trans[sim.df1.wide$from==2&sim.df1.wide$to==1]=4
  sim.df1.wide$trans[sim.df1.wide$from==2&sim.df1.wide$to==2]=5
  sim.df1.wide$trans[sim.df1.wide$from==2&sim.df1.wide$to==3]=6
  sim.df1.wide$trans[sim.df1.wide$from==3&sim.df1.wide$to==1]=7
  sim.df1.wide$trans[sim.df1.wide$from==3&sim.df1.wide$to==2]=8
  sim.df1.wide$trans[sim.df1.wide$from==3&sim.df1.wide$to==3]=9
  dataset1=cbind(1,sim.df1.wide[,c(6,7)])
  
  theta0=c(0.02,0.01,0.01,0.05,0.01,0.03)
  Q= Qmatrix(theta0[1],theta0[2],theta0[3],theta0[4],theta0[5],theta0[6])
  iniQ=crudeinits.msm(state~time, subject = subject, data = sim.df1, qmatrix = Q)
  theta0=round(c(iniQ[1,2],iniQ[1,3],iniQ[2,1],iniQ[2,3],iniQ[3,1],iniQ[3,2]),5)
  
  rr4=optim(theta0,glike_ward,dataset2=dataset1,method = "Nelder-Mead",hessian = T,control=list(fnscale=6000,maxit=1000))  
  return(c(rr4$par,diag(solve(rr4$hessian))) )
}

#Optimization based on Canonical decomposition (N-M)
like_each_decomp=function(trial){
  
  sim.df1=store_all_data1[[trial]]
  sim.df1.wide=ddply(sim.df1,.(subject),msm_long_to_wide)
  sim.df1.wide$time=sim.df1.wide$stop-sim.df1.wide$start
  sim.df1.wide$trans=1
  sim.df1.wide$trans[sim.df1.wide$from==1&sim.df1.wide$to==1]=1
  sim.df1.wide$trans[sim.df1.wide$from==1&sim.df1.wide$to==2]=2
  sim.df1.wide$trans[sim.df1.wide$from==1&sim.df1.wide$to==3]=3
  sim.df1.wide$trans[sim.df1.wide$from==2&sim.df1.wide$to==1]=4
  sim.df1.wide$trans[sim.df1.wide$from==2&sim.df1.wide$to==2]=5
  sim.df1.wide$trans[sim.df1.wide$from==2&sim.df1.wide$to==3]=6
  sim.df1.wide$trans[sim.df1.wide$from==3&sim.df1.wide$to==1]=7
  sim.df1.wide$trans[sim.df1.wide$from==3&sim.df1.wide$to==2]=8
  sim.df1.wide$trans[sim.df1.wide$from==3&sim.df1.wide$to==3]=9
  dataset1=cbind(1,sim.df1.wide[,c(6,7)])
  
  theta0=c(0.02,0.01,0.01,0.05,0.01,0.03)
  Q= Qmatrix(theta0[1],theta0[2],theta0[3],theta0[4],theta0[5],theta0[6])
  iniQ=crudeinits.msm(state~time, subject = subject, data = sim.df1, qmatrix = Q)
  theta0=round(c(iniQ[1,2],iniQ[1,3],iniQ[2,1],iniQ[2,3],iniQ[3,1],iniQ[3,2]),5)
  
  rr4=optim(theta0,glike_decomp,dataset2=dataset1,method = "Nelder-Mead",hessian = T,control=list(fnscale=6000,maxit=1000))  
  return(c(rr4$par,diag(solve(rr4$hessian))) )
}

#Optimization based on Uniformization (N-M)
like_each=function(trial){
  
  sim.df1=store_all_data1[[trial]]
  sim.df1.wide=ddply(sim.df1,.(subject),msm_long_to_wide)
  sim.df1.wide$time=sim.df1.wide$stop-sim.df1.wide$start
  sim.df1.wide$trans=1
  sim.df1.wide$trans[sim.df1.wide$from==1&sim.df1.wide$to==1]=1
  sim.df1.wide$trans[sim.df1.wide$from==1&sim.df1.wide$to==2]=2
  sim.df1.wide$trans[sim.df1.wide$from==1&sim.df1.wide$to==3]=3
  sim.df1.wide$trans[sim.df1.wide$from==2&sim.df1.wide$to==1]=4
  sim.df1.wide$trans[sim.df1.wide$from==2&sim.df1.wide$to==2]=5
  sim.df1.wide$trans[sim.df1.wide$from==2&sim.df1.wide$to==3]=6
  sim.df1.wide$trans[sim.df1.wide$from==3&sim.df1.wide$to==1]=7
  sim.df1.wide$trans[sim.df1.wide$from==3&sim.df1.wide$to==2]=8
  sim.df1.wide$trans[sim.df1.wide$from==3&sim.df1.wide$to==3]=9
  dataset1=cbind(1,sim.df1.wide[,c(6,7)])
  
  theta0=c(0.02,0.01,0.01,0.05,0.01,0.03)
  Q= Qmatrix(theta0[1],theta0[2],theta0[3],theta0[4],theta0[5],theta0[6])
  iniQ=crudeinits.msm(state~time, subject = subject, data = sim.df1, qmatrix = Q)
  theta0=round(c(iniQ[1,2],iniQ[1,3],iniQ[2,1],iniQ[2,3],iniQ[3,1],iniQ[3,2]),5)
  
  
  rr4=optim(theta0,glike,dataset2=dataset1,method = "Nelder-Mead",hessian = T,control=list(fnscale=6000,maxit=1000))  
  return(c(rr4$par,diag(solve(rr4$hessian))) )
}

#Optimization based on Ward approximation (BFGS)
like_each_ward1=function(trial){
  
  sim.df1=store_all_data1[[trial]]
  sim.df1.wide=ddply(sim.df1,.(subject),msm_long_to_wide)
  sim.df1.wide$time=sim.df1.wide$stop-sim.df1.wide$start
  sim.df1.wide$trans=1
  sim.df1.wide$trans[sim.df1.wide$from==1&sim.df1.wide$to==1]=1
  sim.df1.wide$trans[sim.df1.wide$from==1&sim.df1.wide$to==2]=2
  sim.df1.wide$trans[sim.df1.wide$from==1&sim.df1.wide$to==3]=3
  sim.df1.wide$trans[sim.df1.wide$from==2&sim.df1.wide$to==1]=4
  sim.df1.wide$trans[sim.df1.wide$from==2&sim.df1.wide$to==2]=5
  sim.df1.wide$trans[sim.df1.wide$from==2&sim.df1.wide$to==3]=6
  sim.df1.wide$trans[sim.df1.wide$from==3&sim.df1.wide$to==1]=7
  sim.df1.wide$trans[sim.df1.wide$from==3&sim.df1.wide$to==2]=8
  sim.df1.wide$trans[sim.df1.wide$from==3&sim.df1.wide$to==3]=9
  dataset1=cbind(1,sim.df1.wide[,c(6,7)])
  
  theta0=c(0.02,0.01,0.01,0.05,0.01,0.03)
  Q= Qmatrix(theta0[1],theta0[2],theta0[3],theta0[4],theta0[5],theta0[6])
  iniQ=crudeinits.msm(state~time, subject = subject, data = sim.df1, qmatrix = Q)
  theta0=round(c(iniQ[1,2],iniQ[1,3],iniQ[2,1],iniQ[2,3],iniQ[3,1],iniQ[3,2]),5)
  
  rr4=optim(theta0,glike_ward,dataset2=dataset1,method = "BFGS",hessian = T,control=list(fnscale=6000,maxit=1000))  
  return(c(rr4$par,diag(solve(rr4$hessian))) )
}

#Optimization based on Canonical decomposition (BFGS)
like_each_decomp1=function(trial){
  
  sim.df1=store_all_data1[[trial]]
  sim.df1.wide=ddply(sim.df1,.(subject),msm_long_to_wide)
  sim.df1.wide$time=sim.df1.wide$stop-sim.df1.wide$start
  sim.df1.wide$trans=1
  sim.df1.wide$trans[sim.df1.wide$from==1&sim.df1.wide$to==1]=1
  sim.df1.wide$trans[sim.df1.wide$from==1&sim.df1.wide$to==2]=2
  sim.df1.wide$trans[sim.df1.wide$from==1&sim.df1.wide$to==3]=3
  sim.df1.wide$trans[sim.df1.wide$from==2&sim.df1.wide$to==1]=4
  sim.df1.wide$trans[sim.df1.wide$from==2&sim.df1.wide$to==2]=5
  sim.df1.wide$trans[sim.df1.wide$from==2&sim.df1.wide$to==3]=6
  sim.df1.wide$trans[sim.df1.wide$from==3&sim.df1.wide$to==1]=7
  sim.df1.wide$trans[sim.df1.wide$from==3&sim.df1.wide$to==2]=8
  sim.df1.wide$trans[sim.df1.wide$from==3&sim.df1.wide$to==3]=9
  dataset1=cbind(1,sim.df1.wide[,c(6,7)])
  
  theta0=c(0.02,0.01,0.01,0.05,0.01,0.03)
  Q= Qmatrix(theta0[1],theta0[2],theta0[3],theta0[4],theta0[5],theta0[6])
  iniQ=crudeinits.msm(state~time, subject = subject, data = sim.df1, qmatrix = Q)
  theta0=round(c(iniQ[1,2],iniQ[1,3],iniQ[2,1],iniQ[2,3],iniQ[3,1],iniQ[3,2]),5)
  
  rr4=optim(theta0,glike_decomp,dataset2=dataset1,method = "BFGS",hessian = T,control=list(fnscale=6000,maxit=1000))  
  return(c(rr4$par,diag(solve(rr4$hessian))) )
}

#Optimization based on Canonical Uniformization (BFGS)
like_each1=function(trial){
  
  sim.df1=store_all_data1[[trial]]
  sim.df1.wide=ddply(sim.df1,.(subject),msm_long_to_wide)
  sim.df1.wide$time=sim.df1.wide$stop-sim.df1.wide$start
  sim.df1.wide$trans=1
  sim.df1.wide$trans[sim.df1.wide$from==1&sim.df1.wide$to==1]=1
  sim.df1.wide$trans[sim.df1.wide$from==1&sim.df1.wide$to==2]=2
  sim.df1.wide$trans[sim.df1.wide$from==1&sim.df1.wide$to==3]=3
  sim.df1.wide$trans[sim.df1.wide$from==2&sim.df1.wide$to==1]=4
  sim.df1.wide$trans[sim.df1.wide$from==2&sim.df1.wide$to==2]=5
  sim.df1.wide$trans[sim.df1.wide$from==2&sim.df1.wide$to==3]=6
  sim.df1.wide$trans[sim.df1.wide$from==3&sim.df1.wide$to==1]=7
  sim.df1.wide$trans[sim.df1.wide$from==3&sim.df1.wide$to==2]=8
  sim.df1.wide$trans[sim.df1.wide$from==3&sim.df1.wide$to==3]=9
  dataset1=cbind(1,sim.df1.wide[,c(6,7)])
  
  theta0=c(0.02,0.01,0.01,0.05,0.01,0.03)
  Q= Qmatrix(theta0[1],theta0[2],theta0[3],theta0[4],theta0[5],theta0[6])
  iniQ=crudeinits.msm(state~time, subject = subject, data = sim.df1, qmatrix = Q)
  theta0=round(c(iniQ[1,2],iniQ[1,3],iniQ[2,1],iniQ[2,3],iniQ[3,1],iniQ[3,2]),5)
  
  
  rr4=optim(theta0,glike,dataset2=dataset1,method = "BFGS",hessian = T,control=list(fnscale=6000,maxit=1000))  
  return(c(rr4$par,diag(solve(rr4$hessian))) )
}

#Optimization based on Uniformization TBNR
like_each_nr <- function(trial){
  
  sim.df1=store_all_data1[[trial]]
  sim.df1.wide=ddply(sim.df1,.(subject),msm_long_to_wide)
  sim.df1.wide$time=sim.df1.wide$stop-sim.df1.wide$start
  sim.df1.wide$trans=1
  sim.df1.wide$trans[sim.df1.wide$from==1&sim.df1.wide$to==1]=1
  sim.df1.wide$trans[sim.df1.wide$from==1&sim.df1.wide$to==2]=2
  sim.df1.wide$trans[sim.df1.wide$from==1&sim.df1.wide$to==3]=3
  sim.df1.wide$trans[sim.df1.wide$from==2&sim.df1.wide$to==1]=4
  sim.df1.wide$trans[sim.df1.wide$from==2&sim.df1.wide$to==2]=5
  sim.df1.wide$trans[sim.df1.wide$from==2&sim.df1.wide$to==3]=6
  sim.df1.wide$trans[sim.df1.wide$from==3&sim.df1.wide$to==1]=7
  sim.df1.wide$trans[sim.df1.wide$from==3&sim.df1.wide$to==2]=8
  sim.df1.wide$trans[sim.df1.wide$from==3&sim.df1.wide$to==3]=9
  dataset1=cbind(1,sim.df1.wide[,c(6,7)])
  
  theta0=c(0.02,0.01,0.01,0.05,0.01,0.03)
  Q= Qmatrix(theta0[1],theta0[2],theta0[3],theta0[4],theta0[5],theta0[6])
  iniQ=crudeinits.msm(state~time, subject = subject, data = sim.df1, qmatrix = Q)
  theta0=round(c(iniQ[1,2],iniQ[1,3],iniQ[2,1],iniQ[2,3],iniQ[3,1],iniQ[3,2]),5)
  
  params1=theta0
  jump_matrix=getallinten_second(params1,dataset1)
  test=as.matrix(jump_matrix )
  
  
  seeyixa=likelihood_part2(test)
  
  
  socrevec=seeyixa[1:6]
  
  m_upper=seeyixa[7:27]
  l <- length(m_upper)
  n <- length(socrevec)
  # Reconstruct 
  m2 <- matrix(NA, n, n)
  m2[lower.tri(m2, diag = TRUE)] <- m_upper
  m2=t(m2)
  m2[lower.tri(m2)] <- t(m2)[lower.tri(m2)]  # If symmetric, fill also upper half
  infmatrix=m2
  
  
  params2=solve(infmatrix,infmatrix %*%params1-socrevec)
  itra=1
  while  (sqrt(sum((params2-params1)^2))>0.0001) {
    params1=params2
    jump_matrix=getallinten_second(params1,dataset1)
    
    test=as.matrix(jump_matrix )
    seeyixa=likelihood_part2(test)
    
    socrevec=seeyixa[1:6]
    
    m_upper=seeyixa[7:27]
    
    l <- length(m_upper)
    n <- length(socrevec)
    # Reconstruct 
    m2 <- matrix(NA, n, n)
    m2[lower.tri(m2, diag = TRUE)] <- m_upper
    m2=t(m2)
    m2[lower.tri(m2)] <- t(m2)[lower.tri(m2)]  # If symmetric, fill also upper half
    infmatrix=m2
    
    params2=solve(infmatrix,infmatrix %*%params1-socrevec)
    
    itra=itra+1
    print(itra)
    
  }
  
  
  return(c(params2,diag(solve(infmatrix))))
  
  
}

#Optimization based on Canonical decomposition Fisher
like_each_fisher <- function(trial){
  skip_to_next <- FALSE
  withTimeout(  
  tryCatch(
    
    {
      sim.df1=store_all_data1[[trial]]
      sim.df1.wide=ddply(sim.df1,.(subject),msm_long_to_wide)
      sim.df1.wide$time=sim.df1.wide$stop-sim.df1.wide$start
      sim.df1.wide$trans=1
      sim.df1.wide$trans[sim.df1.wide$from==1&sim.df1.wide$to==1]=1
      sim.df1.wide$trans[sim.df1.wide$from==1&sim.df1.wide$to==2]=2
      sim.df1.wide$trans[sim.df1.wide$from==1&sim.df1.wide$to==3]=3
      sim.df1.wide$trans[sim.df1.wide$from==2&sim.df1.wide$to==1]=4
      sim.df1.wide$trans[sim.df1.wide$from==2&sim.df1.wide$to==2]=5
      sim.df1.wide$trans[sim.df1.wide$from==2&sim.df1.wide$to==3]=6
      sim.df1.wide$trans[sim.df1.wide$from==3&sim.df1.wide$to==1]=7
      sim.df1.wide$trans[sim.df1.wide$from==3&sim.df1.wide$to==2]=8
      sim.df1.wide$trans[sim.df1.wide$from==3&sim.df1.wide$to==3]=9
      dataset1=cbind(1,sim.df1.wide[,c(6,7)])
      
      theta0=c(0.02,0.01,0.01,0.05,0.01,0.03)
      Q= Qmatrix(theta0[1],theta0[2],theta0[3],theta0[4],theta0[5],theta0[6])
      iniQ=crudeinits.msm(state~time, subject = subject, data = sim.df1, qmatrix = Q)
      theta0=round(c(iniQ[1,2],iniQ[1,3],iniQ[2,1],iniQ[2,3],iniQ[3,1],iniQ[3,2]),5)
      
      Q= Qmatrix(theta0[1],theta0[2],theta0[3],theta0[4],theta0[5],theta0[6])
      model <- msm(state~time, subject = subject, data = sim.df1, qmatrix = Q,opt.method="fisher",analyticp = FALSE,control=list(fnscale=6000,maxit=1000))
      p          <- model$estimates
      p.se       <- sqrt(diag(model$covmat))
      estimationt=model$estimates.t
      estimatestion.se=(exp( p+1.96*p.se)-model$estimates.t)/1.96
      
      return(c(estimationt ,  estimatestion.se))
      
    }
    , error = function(e) { skip_to_next <<- TRUE; print(e$message)}), timeout = 60)

}

#figures
benchmark_plot <- function(x){
  #results 
  log=TRUE
  object=x
  y_max=1.05 * max(object$time)
  y_min <- 0
  object$ntime <- convert_to_unit(object$time, "ms")
  #object=object[  !(object$ntime >=10000 & object$expr=="fishe" ),]
  
  object <- object %>% filter(!(expr %in% c("generate")))
  dt=data.frame("expr"=object$expr, "ntime"=object$ntime)
  dt=dt[dt$expr %in% c("Ward", "diagnoal", "uniformization1", "uniformization2","fisher" ),]
  dt$expr=as.factor(dt$expr)
  dt$ntime=log(dt$ntime)
  plot1=dt %>% mutate(expr = factor(expr, levels=c("uniformization2", "uniformization1","fisher","diagnoal","Ward"))) %>%
    ggplot2::ggplot(ggplot2::aes_string(x="expr", y="ntime",fill="expr")) + geom_boxplot(# custom boxes alpha=0.8,
      
      # Notch?
      notch=TRUE,
      notchwidth = 0.8,
      
      # custom outliers
      outlier.colour="#EB984E",
      outlier.fill="#EB984E",
      outlier.size=3)+  scale_fill_manual(values=c("#A569BD","#FFEE58","#5DADE2","#73C6B6","#85929E"))+
    ggplot2::coord_flip()+
    ggplot2::theme_bw() + theme(panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
  
  return(plot1)
}

