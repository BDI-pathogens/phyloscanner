# transform parameters to R
transform_parameters = function(alpha, rho, mu, baseline, beta, xi,all=T){
  tmp_alpha = log(alpha)
  tmp_rho = log(rho)
  tmp_rho = aperm(rho,c(1,3,2))
  dim(tmp_rho) =  c(dim(tmp_rho)[1], (dim(tmp_rho)[2]* dim(tmp_rho)[3]))
  tmp_beta = aperm(beta,c(1,3,2))
  dim(tmp_beta) =  c(dim(tmp_beta)[1], (dim(tmp_beta)[2]* dim(tmp_beta)[3]))
  if(all==T){
    tmp_xi = log(xi/(1-xi))
  }else{
    tmp_xi = xi
  }
  return(list(hyperpars=as.matrix(cbind(tmp_alpha, tmp_rho,mu,baseline)),
              beta=as.matrix(tmp_beta),
              xi=as.matrix(tmp_xi)))
}

# changes to proposal density by variable transformation
transform_parameters_density_trans = function(trans_pars, samples_dim){
  ngp = samples_dim$alpha
  xi = 1/ (1 + exp(-trans_pars[[3]]))
  ans = -apply(trans_pars[[1]][,1:(3*ngp)], 1, sum) - apply(log(xi), 1, sum)  -apply(log(1-xi), 1, sum)
  return(ans)
}

# transform parameters back
transform_back_parameters = function(trans_pars, samples_dim){
  trans_pars[[1]][,1:(samples_dim$alpha + samples_dim$rho)]= exp(trans_pars[[1]][,1:(samples_dim$alpha + samples_dim$rho)])
  tmp = exp(trans_pars[[3]])
  trans_pars[[3]] = tmp/(1+tmp)
  return(trans_pars)
}

# fit proposal density
fit_proposal_density = function(trans_pars){
  trans_pars_fit = mvnorm.mle(do.call(cbind, trans_pars))
  tmp = try(chol(trans_pars_fit$sigma),silent=TRUE)
  if(is(tmp, 'try-error')){
    trans_pars_fit$sigma = trans_pars_fit$sigma + 10^(-10) * diag(nrow(trans_pars_fit$sigma))
  }
  return(trans_pars_fit)
}

# generate the transformed parameters from the proposal
proposal_density_mvnorm_rand = function(proposal_density_pars, samples_dim, N2){
  pars = rmvnorm(N2,proposal_density_pars$mu, proposal_density_pars$sigma)
  ans = list()
  dim1 = samples_dim$alpha + samples_dim$rho + samples_dim$mu + samples_dim$baseline
  dim2 = samples_dim$beta
  dim3 = samples_dim$xi
  ans[[1]] = pars[,1:dim1]
  ans[[2]] = pars[,(dim1 + 1):(dim1 + dim2)]
  ans[[3]] = pars[,(dim1 + dim2 + 1):length(proposal_density_pars$mu)]
  return(ans)
}

# kernel 
kernel_se_spd <- function(alpha,rho,slambda){
  spd = 2 * pi * alpha^2 * prod(rho) * exp(-0.5 * (slambda^2 %*% matrix(rho,ncol = 1)^2))
  return(spd)
}

# rep functions
rep.row<-function(x,n){
  matrix(rep(x,each=n),nrow=n)
}
rep.col<-function(x,n){
  matrix(rep(x,each=n), ncol=n, byrow=TRUE)
}

# log poisson density with a log rate
dpois_lograte =function(x,llbd){
  ans = x * llbd - exp(llbd) - lgamma(x+1)
  return(ans)
}

# log unnormalised posterior density
log_posterior_density= function(trans_pars, data.fit, samples_dim, index, verbose){
  # transform parameters back
  samples = transform_back_parameters(trans_pars, samples_dim)
  # #gp
  ngp = samples_dim$alpha
  # #intercept
  nintercept = samples_dim$mu
  for (k in 1:ngp){
    spd = sapply(1:nrow(samples[[1]]), 
                 function(x){kernel_se_spd(samples[[1]][x,k], samples[[1]][x,(ngp+2*k-1):(ngp+2*k)], data.fit$slambda)})
    assign(paste0('spd',k),spd)
    f = data.fit$Xgp %*% (sqrt(get(paste0('spd',k))) * t(samples[[2]][,(900*(k-1)+1):(900*k)]))
    assign(paste0('f',k),f)
  }
  plrate = array(dim = c(nrow(samples[[1]]),dim(data.fit$y)))
  if(index=='1'){
    tmp_invgamma =  dinvgamma(samples[[1]][,3],2.93854,8.30173, log = T) + dinvgamma(samples[[1]][,4],2.60877,7.73394, log = T)+
      dinvgamma(samples[[1]][,5],3.48681,9.21604, log = T) + dinvgamma(samples[[1]][,6],2.60877,7.73394, log = T) 
    for (i in 1:ncol(data.fit$y)) {
      if(i %in% data.fit$id_mf){
        plrate[,,i] = t(f1 + rep.row(samples[[1]][,ncol(samples[[1]])], nrow(f2)) + rep.row(samples[[1]][,(3*ngp+1)], nrow(f2)))  + log(samples[[3]][,data.fit$xi_id_src[,i]]) + log(samples[[3]][,data.fit$xi_id_rec[,i]])
      }else{
        plrate[,,i] = t(f2 + rep.row(samples[[1]][,ncol(samples[[1]])], nrow(f2)))  + log(samples[[3]][,data.fit$xi_id_src[,i]]) + log(samples[[3]][,data.fit$xi_id_rec[,i]])
      }
    }
  }else if(index=='2'){
    tmp = dinvgamma(samples[[1]][,c(ngp + 1,ngp + 3)],2.93854,8.30173, log = T) 
    dim(tmp) = c(nrow(samples[[1]]), 2)
    tmp_invgamma = apply(tmp,1,sum)
    tmp = dinvgamma(samples[[1]][,c(ngp + 2,ngp + 4)],2.60877,7.73394, log = T) 
    dim(tmp) = c(nrow(samples[[1]]), 2)
    tmp_invgamma = tmp_invgamma + apply(tmp,1,sum)
    tmp = dinvgamma(samples[[1]][,c(ngp + 5,ngp + 7)],3.48681,9.21604, log = T) 
    dim(tmp) = c(nrow(samples[[1]]), 2)
    tmp_invgamma = tmp_invgamma + apply(tmp,1,sum)
    tmp = dinvgamma(samples[[1]][,c(ngp + 6,ngp + 8)],2.60877,7.73394, log = T) 
    dim(tmp) = c(nrow(samples[[1]]), 2)
    tmp_invgamma = tmp_invgamma + apply(tmp,1,sum)
    for (i in 1:ncol(data.fit$y)){
      if(i %in% data.fit$id_mf_h){
        plrate[,,i] = t(f1 + rep.row(samples[[1]][,ncol(samples[[1]])], nrow(f2)) + rep.row(samples[[1]][,(3*ngp+1)], nrow(f2)))  + log(samples[[3]][,data.fit$xi_id_src[,i]]) + log(samples[[3]][,data.fit$xi_id_rec[,i]])
      }else if(i %in% data.fit$id_mf_l){
        plrate[,,i] = t(f2 + rep.row(samples[[1]][,ncol(samples[[1]])], nrow(f2))+ rep.row(samples[[1]][,(3*ngp+1)], nrow(f2))) + log(samples[[3]][,data.fit$xi_id_src[,i]]) + log(samples[[3]][,data.fit$xi_id_rec[,i]])
      }else if(i %in% data.fit$id_fm_h){
        plrate[,,i] = t(f3 + rep.row(samples[[1]][,ncol(samples[[1]])], nrow(f2)))  + log(samples[[3]][,data.fit$xi_id_src[,i]]) + log(samples[[3]][,data.fit$xi_id_rec[,i]])
      }else if(i %in% data.fit$id_fm_l){
        plrate[,,i] = t(f4 + rep.row(samples[[1]][,ncol(samples[[1]])], nrow(f2)))  + log(samples[[3]][,data.fit$xi_id_src[,i]]) + log(samples[[3]][,data.fit$xi_id_rec[,i]])
      }
    }
  }else if(index=='3'){
    tmp_invgamma =  dinvgamma(samples[[1]][,3],2.93854,8.30173, log = T) + dinvgamma(samples[[1]][,4],2.60877,7.73394, log = T)+
      dinvgamma(samples[[1]][,5],3.48681,9.21604, log = T) + dinvgamma(samples[[1]][,6],2.60877,7.73394, log = T) 
    for (i in 1:ncol(data.fit$y)) {
      if(i %in% data.fit$id_mf){
        plrate[,,i] = t(f1 + rep.row(samples[[1]][,ncol(samples[[1]])], nrow(f2)) + rep.row(samples[[1]][,(3*ngp+1)], nrow(f2)))  + log(samples[[3]][,data.fit$xi_id_src[,i]]) + log(samples[[3]][,data.fit$xi_id_rec[,i]])
      }else{
        plrate[,,i] = t(f2 + rep.row(samples[[1]][,ncol(samples[[1]])], nrow(f2)))  + log(samples[[3]][,data.fit$xi_id_src[,i]]) + log(samples[[3]][,data.fit$xi_id_rec[,i]])
      }
      
      if(i %in% data.fit$id_hh){
        plrate[,,i] = plrate[,,i] + rep.col(samples[[1]][,(3*ngp+2)],ncol(plrate[,,i]))
      }else if(i %in% data.fit$id_hl){
        plrate[,,i] = plrate[,,i] + rep.col(samples[[1]][,(3*ngp+3)],ncol(plrate[,,i]))
      }else if(i %in% data.fit$id_lh){
        plrate[,,i] = plrate[,,i] + rep.col(samples[[1]][,(3*ngp+4)],ncol(plrate[,,i]))
      }else if(i %in% data.fit$id_eh){
        plrate[,,i] = plrate[,,i] + rep.col(samples[[1]][,(3*ngp+5)],ncol(plrate[,,i]))
      }else if(i %in% data.fit$id_el){
        plrate[,,i] = plrate[,,i] + rep.col(samples[[1]][,(3*ngp+6)],ncol(plrate[,,i]))
      }
    }
  }else if(index=='4'){
    tmp = dinvgamma(samples[[1]][,c(ngp + 1,ngp + 3)],2.93854,8.30173, log = T) 
    dim(tmp) = c(nrow(samples[[1]]), 2)
    tmp_invgamma = apply(tmp,1,sum)
    tmp = dinvgamma(samples[[1]][,c(ngp + 2,ngp + 4)],2.60877,7.73394, log = T) 
    dim(tmp) = c(nrow(samples[[1]]), 2)
    tmp_invgamma = tmp_invgamma + apply(tmp,1,sum)
    tmp = dinvgamma(samples[[1]][,c(ngp + 5,ngp + 7)],3.48681,9.21604, log = T) 
    dim(tmp) = c(nrow(samples[[1]]), 2)
    tmp_invgamma = tmp_invgamma + apply(tmp,1,sum)
    tmp = dinvgamma(samples[[1]][,c(ngp + 6,ngp + 8)],2.60877,7.73394, log = T) 
    dim(tmp) = c(nrow(samples[[1]]), 2)
    tmp_invgamma = tmp_invgamma + apply(tmp,1,sum)
    for (i in 1:ncol(data.fit$y)){
      if(i %in% data.fit$id_mf_h){
        plrate[,,i] = t(f1 + rep.row(samples[[1]][,ncol(samples[[1]])], nrow(f2)) + rep.row(samples[[1]][,(3*ngp+1)], nrow(f2)))  + log(samples[[3]][,data.fit$xi_id_src[,i]]) + log(samples[[3]][,data.fit$xi_id_rec[,i]])
      }else if(i %in% data.fit$id_mf_l){
        plrate[,,i] = t(f2 + rep.row(samples[[1]][,ncol(samples[[1]])], nrow(f2)) + rep.row(samples[[1]][,(3*ngp+1)], nrow(f2)))  + log(samples[[3]][,data.fit$xi_id_src[,i]]) + log(samples[[3]][,data.fit$xi_id_rec[,i]])
      }else if(i %in% data.fit$id_fm_h){
        plrate[,,i] = t(f3 + rep.row(samples[[1]][,ncol(samples[[1]])], nrow(f2)))  + log(samples[[3]][,data.fit$xi_id_src[,i]]) + log(samples[[3]][,data.fit$xi_id_rec[,i]])
      }else if(i %in% data.fit$id_fm_l){
        plrate[,,i] = t(f4 + rep.row(samples[[1]][,ncol(samples[[1]])], nrow(f2)))  + log(samples[[3]][,data.fit$xi_id_src[,i]]) + log(samples[[3]][,data.fit$xi_id_rec[,i]])
      }
      if(i %in% data.fit$id_hh){
        plrate[,,i] = plrate[,,i] + rep.col(samples[[1]][,(3*ngp+2)],ncol(plrate[,,i]))
      }else if(i %in% data.fit$id_hl){
        plrate[,,i] = plrate[,,i] + rep.col(samples[[1]][,(3*ngp+3)],ncol(plrate[,,i]))
      }else if(i %in% data.fit$id_lh){
        plrate[,,i] = plrate[,,i] + rep.col(samples[[1]][,(3*ngp+4)],ncol(plrate[,,i]))
      }else if(i %in% data.fit$id_eh){
        plrate[,,i] = plrate[,,i] + rep.col(samples[[1]][,(3*ngp+5)],ncol(plrate[,,i]))
      }else if(i %in% data.fit$id_el){
        plrate[,,i] = plrate[,,i] + rep.col(samples[[1]][,(3*ngp+6)],ncol(plrate[,,i]))
      }
    }
  }else if(index=='5'){
    tmp = dinvgamma(samples[[1]][,c(ngp + 1,ngp + 3)],2.93854,8.30173, log = T) 
    dim(tmp) = c(nrow(samples[[1]]), 2)
    tmp_invgamma = apply(tmp,1,sum)
    tmp = dinvgamma(samples[[1]][,c(ngp + 2,ngp + 4)],2.60877,7.73394, log = T) 
    dim(tmp) = c(nrow(samples[[1]]), 2)
    tmp_invgamma = tmp_invgamma + apply(tmp,1,sum)
    tmp = dinvgamma(samples[[1]][,c(ngp + 5,ngp + 7)],3.48681,9.21604, log = T) 
    dim(tmp) = c(nrow(samples[[1]]), 2)
    tmp_invgamma = tmp_invgamma + apply(tmp,1,sum)
    tmp = dinvgamma(samples[[1]][,c(ngp + 6,ngp + 8)],2.60877,7.73394, log = T) 
    dim(tmp) = c(nrow(samples[[1]]), 2)
    tmp_invgamma = tmp_invgamma + apply(tmp,1,sum)
    tmp = dinvgamma(samples[[1]][,c(ngp + 9,ngp + 11,ngp + 13)],2.93854,8.30173, log = T) 
    dim(tmp) = c(nrow(samples[[1]]), 3)
    tmp_invgamma = tmp_invgamma + apply(tmp,1,sum)
    tmp = dinvgamma(samples[[1]][,c(ngp + 10,ngp + 12,ngp + 14)],2.60877,7.73394, log = T) 
    dim(tmp) = c(nrow(samples[[1]]), 3)
    tmp_invgamma = tmp_invgamma + apply(tmp,1,sum)
    for (i in 1:ncol(data.fit$y)){
      if(i %in% data.fit$id_mf_h){
        plrate[,,i] = t(f1 + rep.row(samples[[1]][,ncol(samples[[1]])], nrow(f2)) + rep.row(samples[[1]][,(3*ngp+1)], nrow(f2)))  + log(samples[[3]][,data.fit$xi_id_src[,i]]) + log(samples[[3]][,data.fit$xi_id_rec[,i]])
      }else if(i %in% data.fit$id_mf_l){
        plrate[,,i] = t(f2 + rep.row(samples[[1]][,ncol(samples[[1]])], nrow(f2))+ rep.row(samples[[1]][,(3*ngp+1)], nrow(f2)))   + log(samples[[3]][,data.fit$xi_id_src[,i]]) + log(samples[[3]][,data.fit$xi_id_rec[,i]])
      }else if(i %in% data.fit$id_fm_h){
        plrate[,,i] = t(f3 + rep.row(samples[[1]][,ncol(samples[[1]])], nrow(f2)))  + log(samples[[3]][,data.fit$xi_id_src[,i]]) + log(samples[[3]][,data.fit$xi_id_rec[,i]])
      }else if(i %in% data.fit$id_fm_l){
        plrate[,,i] = t(f4 + rep.row(samples[[1]][,ncol(samples[[1]])], nrow(f2)))  + log(samples[[3]][,data.fit$xi_id_src[,i]]) + log(samples[[3]][,data.fit$xi_id_rec[,i]])
      }
      if(i %in% data.fit$id_se){
        plrate[,,i] = plrate[,,i] + t(f5)
      }else if(i %in% data.fit$id_sl){
        plrate[,,i] = plrate[,,i] + t(f6)
      }else if(i %in% data.fit$id_sh){
        plrate[,,i] = plrate[,,i] + t(f7)
      }
    }
  }else if(index=='6'){
    tmp = dinvgamma(samples[[1]][,c(ngp + 1,ngp + 3)],2.93854,8.30173, log = T) 
    dim(tmp) = c(nrow(samples[[1]]), 2)
    tmp_invgamma = apply(tmp,1,sum)
    tmp = dinvgamma(samples[[1]][,c(ngp + 2,ngp + 4)],2.60877,7.73394, log = T) 
    dim(tmp) = c(nrow(samples[[1]]), 2)
    tmp_invgamma = tmp_invgamma + apply(tmp,1,sum)
    tmp = dinvgamma(samples[[1]][,c(ngp + 5,ngp + 7)],3.48681,9.21604, log = T) 
    dim(tmp) = c(nrow(samples[[1]]), 2)
    tmp_invgamma = tmp_invgamma + apply(tmp,1,sum)
    tmp = dinvgamma(samples[[1]][,c(ngp + 6,ngp + 8)],2.60877,7.73394, log = T) 
    dim(tmp) = c(nrow(samples[[1]]), 2)
    tmp_invgamma = tmp_invgamma + apply(tmp,1,sum)
    tmp = dinvgamma(samples[[1]][,c(ngp + 9,ngp + 11,ngp + 13)],2.93854,8.30173, log = T) 
    dim(tmp) = c(nrow(samples[[1]]), 3)
    tmp_invgamma = tmp_invgamma + apply(tmp,1,sum)
    tmp = dinvgamma(samples[[1]][,c(ngp + 10,ngp + 12,ngp + 14)],2.60877,7.73394, log = T) 
    dim(tmp) = c(nrow(samples[[1]]), 3)
    tmp_invgamma = tmp_invgamma + apply(tmp,1,sum)
    for (i in 1:ncol(data.fit$y)){
      if(i %in% data.fit$id_mf_h){
        plrate[,,i] = t(f1 + rep.row(samples[[1]][,ncol(samples[[1]])], nrow(f2)) + rep.row(samples[[1]][,(3*ngp+1)], nrow(f2)))  + log(samples[[3]][,data.fit$xi_id_src[,i]]) + log(samples[[3]][,data.fit$xi_id_rec[,i]])
      }else if(i %in% data.fit$id_mf_l){
        plrate[,,i] = t(f2 + rep.row(samples[[1]][,ncol(samples[[1]])], nrow(f2)) + rep.row(samples[[1]][,(3*ngp+1)], nrow(f2)))  + log(samples[[3]][,data.fit$xi_id_src[,i]]) + log(samples[[3]][,data.fit$xi_id_rec[,i]])
      }else if(i %in% data.fit$id_fm_h){
        plrate[,,i] = t(f3 + rep.row(samples[[1]][,ncol(samples[[1]])], nrow(f2)))  + log(samples[[3]][,data.fit$xi_id_src[,i]]) + log(samples[[3]][,data.fit$xi_id_rec[,i]])
      }else if(i %in% data.fit$id_fm_l){
        plrate[,,i] = t(f4 + rep.row(samples[[1]][,ncol(samples[[1]])], nrow(f2)))  + log(samples[[3]][,data.fit$xi_id_src[,i]]) + log(samples[[3]][,data.fit$xi_id_rec[,i]])
      }
      if(i %in% data.fit$id_se){
        plrate[,,i] = plrate[,,i] + t(f5)
      }else if(i %in% data.fit$id_sl){
        plrate[,,i] = plrate[,,i] + t(f6)
      }else if(i %in% data.fit$id_sh){
        plrate[,,i] = plrate[,,i] + t(f7)
      }
      if(i %in% data.fit$id_hh){
        plrate[,,i] = plrate[,,i] + rep.col(samples[[1]][,(3*ngp+2)],ncol(plrate[,,i]))
      }else if(i %in% data.fit$id_hl){
        plrate[,,i] = plrate[,,i] + rep.col(samples[[1]][,(3*ngp+3)],ncol(plrate[,,i]))
      }else if(i %in% data.fit$id_lh){
        plrate[,,i] = plrate[,,i] + rep.col(samples[[1]][,(3*ngp+4)],ncol(plrate[,,i]))
      }else if(i %in% data.fit$id_eh){
        plrate[,,i] = plrate[,,i] + rep.col(samples[[1]][,(3*ngp+5)],ncol(plrate[,,i]))
      }else if(i %in% data.fit$id_el){
        plrate[,,i] = plrate[,,i] + rep.col(samples[[1]][,(3*ngp+6)],ncol(plrate[,,i]))
      }
    }
  }else if(index=='7'){
    tmp_invgamma =  dinvgamma(samples[[1]][,2],2.93854,8.30173, log = T) + dinvgamma(samples[[1]][,3],2.60877,7.73394, log = T)
    for (i in 1:ncol(data.fit$y)) {
      plrate[,,i] =t(f1 + rep.row(samples[[1]][,ncol(samples[[1]])], nrow(f1)))  + 
        log(samples[[3]][,data.fit$xi_id_src[,i]]) + log(samples[[3]][,data.fit$xi_id_rec[,i]])
    }
  }
  tmp_pois = sapply(1:dim(plrate)[1], function(x)sum(dpois_lograte(data.fit$y, plrate[x,,])))
  if(verbose){
    cat('likelihood : ', paste0(median(tmp_pois), ' [ ', quantile(tmp_pois,probs = 0.025),' - ',quantile(tmp_pois,probs = 0.975) ,' ] '), '\n')
  }
  tmp_hn = dhnorm(samples[[1]][,1:ngp],10, log = T)
  dim(tmp_hn) = c(nrow(samples[[1]]),ngp)
  tmp_hn = apply(tmp_hn, 1, sum)
  if(verbose){
    cat('prior - alpha : ', paste0(median(tmp_hn), ' [ ', quantile(tmp_hn,probs = 0.025),' - ',quantile(tmp_hn,probs = 0.975) ,' ] '), '\n')
  }
  tmp_n = dnorm(samples[[1]][,(3*ngp+1):ncol(samples[[1]])], 0, 10, log = T) 
  dim(tmp_n) = c(nrow(samples[[1]]),ncol(samples[[1]])-3*ngp)
  tmp_n = apply(tmp_n, 1, sum)
  if(verbose){
    cat('prior - intercept : ',  paste0(median(tmp_n), ' [ ', quantile(tmp_n,probs = 0.025),' - ',quantile(tmp_n,probs = 0.975) ,' ] '), '\n')
  }
  tmp_beta = dbeta(samples[[3]], rep.row(data.fit$shape[,1], nrow(samples[[3]])), 
                   rep.row(data.fit$shape[,2], nrow(samples[[3]])), log=T)
  ans = tmp_pois + 
    apply(dnorm(samples[[2]], 0, 1, log = T), 1, sum) + 
    tmp_hn + tmp_n + tmp_invgamma +
    apply(tmp_beta, 1, sum) 
  
  return(ans)
}

# log density of proposal
log_proposal_density = function(trans_pars, trans_pars_fit, samples_dim){
  x = do.call(cbind,trans_pars)
  ans = mvtnorm::dmvnorm(x,trans_pars_fit$mu,trans_pars_fit$sigma,log=T) + transform_parameters_density_trans(trans_pars, samples_dim)
  return(ans)
}

# set up
setwd('/rds/general/user/xx4515/ephemeral/age35_gp_2f_tmp/')
files = list.files('/rds/general/user/xx4515/ephemeral/age35_gp_2f_tmp/')
args <- commandArgs(trailingOnly = TRUE)
print(args)
files = files[grepl(paste0('gpa', args[1]), files)]
randomnumber <- round(runif(1) * 100000000)
set.seed(randomnumber)
files <- sample(files)
ptm <- Sys.time()
verbose <- T
# load packages
library(rstan)
library(Rfast)
library(MASS)
library(extraDistr)
library(data.table)
library(Brobdingnag)
library(coda)
library(mvtnorm)

# load a set of samples
fit_list = list()
for (i in 1:4) {
  cat('load ', files[i], '\n')
  load(files[i])
  fit_list[[i]] = fit
}
fit = sflist2stanfit(fit_list)

if(args[1]=='7'){
  samples = extract(fit, pars = c('beta', 'baseline','alpha', 'rho', 'xi'))
  samples$mu <- matrix(nrow = length(samples$baseline), ncol=0)
}else{
  samples = extract(fit, pars = c('beta', 'mu', 'baseline','alpha', 'rho', 'xi'))
}
# print(lapply(samples, dim))
samples_dim = lapply(samples, function(x) {
  y = dim(x)[2:3]
  y[is.na(y)] = 1
  return(prod(y))
})


# load input
load('input.rda')

# determine N1
tmp = lapply(samples, function(x)
  matrix(x, nrow = length(samples$baseline)))
tmp = do.call(cbind, tmp)
N1 = median(apply(tmp, 2, function(x)
  effectiveSize(mcmc(x))))
cat('N1 = ', N1)

# setup
N2 = nrow(samples$beta)
s1 = N1 / (N1 + N2)
s2 = N2 / (N1 + N2)
tol = 10 ^ (-10)

# initial
nfile = floor(length(files) / 4)
log_marginal_likelihood = matrix(nrow = 1, ncol = nfile)
log_marginal_likelihood[1] = -2e3

# proposal
if (args[1] == '7') {
  dim(samples$rho) <- c(dim(samples$rho), 1)
  dim(samples$beta) <- c(dim(samples$beta), 1)
}
trans_pars = transform_parameters(
  samples$alpha,
  samples$rho,
  samples$mu,
  samples$baseline,
  samples$beta,
  samples$xi,
  all = T
)
trans_pars_fit_mvn = fit_proposal_density(trans_pars)


for (i in 2:nfile) {
  cat('process ', i, 'th out of ', nfile, ' files \n')
  # proposed samples
  proposed_trans_pars = matrix(rnorm(N2 * length(trans_pars_fit_mvn$mu), 0, 1), ncol =
                                 length(trans_pars_fit_mvn$mu))
  
  # posterior samples
  fit_list = list()
  for (k in 1:4) {
    j = ((4 * (i - 1) + 1):(i * 4))[k]
    cat('load ', files[j], '\n')
    load(files[j])
    fit_list[[k]] = fit
  }
  fit = sflist2stanfit(fit_list)
  if(args[1]=='7'){
    samples = extract(fit, pars = c('beta', 'baseline','alpha', 'rho', 'xi'))
    samples$mu <- matrix(nrow = length(samples$baseline), ncol=0)
  }else{
    samples = extract(fit, pars = c('beta', 'mu', 'baseline','alpha', 'rho', 'xi'))
  }
  if (args[1] == '7') {
    dim(samples$rho) <- c(dim(samples$rho), 1)
    dim(samples$beta) <- c(dim(samples$beta), 1)
  }
  trans_pars <-
    transform_parameters(
      samples$alpha,
      samples$rho,
      samples$mu,
      samples$baseline,
      samples$beta,
      samples$xi,
      all = T
    )
  
  # evaluate q(thetastar)
  dpost_samples = log_posterior_density(trans_pars, data.fit, samples_dim, args[1], verbose)
  
  
  # posterior mean and covariance
  dim1 = samples_dim$alpha + samples_dim$rho + samples_dim$mu + samples_dim$baseline
  dim2 = samples_dim$beta
  dim3 = samples_dim$xi
  tmp_mu = list()
  tmp_mu[[1]] = trans_pars_fit_mvn$mu[1:dim1]
  tmp_mu[[2]] = trans_pars_fit_mvn$mu[(dim1 + 1):(dim2 + dim1)]
  tmp_mu[[3]] = trans_pars_fit_mvn$mu[(dim1 + dim2 + 1):length(trans_pars_fit_mvn$mu)]
  
  # evaluate the two mu term
  tmp_trans_pars = list()
  for (kk in 1:3) {
    tmp_trans_pars[[kk]] = 2 * rep.row(tmp_mu[[kk]], nrow(trans_pars[[kk]])) - trans_pars[[kk]]
  }
  dpost_samples2 = log_posterior_density(tmp_trans_pars, data.fit, samples_dim, args[1], verbose)
  R = t(chol(trans_pars_fit_mvn$sigma))
  Rinv = forwardsolve(R, diag(nrow(R)))
  Rlogdet = determinant(R, logarithm = T)
  
  # evaluate the minus transform
  tmp = R %*% t(proposed_trans_pars)
  tmp_prop_trans_pars = rep.col(trans_pars_fit_mvn$mu, N2) - tmp
  tmp_prop_trans_pars_l = list()
  tmp_prop_trans_pars_l[[1]] = t(tmp_prop_trans_pars[1:dim1, ])
  tmp_prop_trans_pars_l[[2]] = t(tmp_prop_trans_pars[(1 + dim1):(dim1 + dim2), ])
  tmp_prop_trans_pars_l[[3]] = t(tmp_prop_trans_pars[(dim1 + dim2 + 1):nrow(tmp_prop_trans_pars), ])
  dpost_proposed_samples1 = log_posterior_density(tmp_prop_trans_pars_l, data.fit, samples_dim, args[1], verbose)
  
  # evaluate the plus transform
  tmp_prop_trans_pars = rep.col(trans_pars_fit_mvn$mu, N2) + tmp
  tmp_prop_trans_pars_u = list()
  tmp_prop_trans_pars_u[[1]] = t(tmp_prop_trans_pars[1:dim1, ])
  tmp_prop_trans_pars_u[[2]] = t(tmp_prop_trans_pars[(1 + dim1):(dim1 + dim2), ])
  tmp_prop_trans_pars_u[[3]] = t(tmp_prop_trans_pars[(dim1 + dim2 + 1):nrow(tmp_prop_trans_pars), ])
  dpost_proposed_samples2 = log_posterior_density(tmp_prop_trans_pars_u, data.fit, samples_dim, args[1], verbose)
  
  dprop_proposed_samples = dmvnorm(proposed_trans_pars, rep(0, ncol(proposed_trans_pars)), diag(ncol(proposed_trans_pars)), log =T) +
    transform_parameters_density_trans(trans_pars, samples_dim)
  
  tmp_samples <- do.call(cbind, trans_pars)
  tmp_samples <-
    t(Rinv %*% t(tmp_samples - rep.row(
      trans_pars_fit_mvn$mu, nrow(tmp_samples)
    )))
  dprop_samples1 = dmvnorm(tmp_samples, rep(0, ncol(tmp_samples)), diag(ncol(tmp_samples)), log =
                             T) +
    transform_parameters_density_trans(trans_pars, samples_dim)
  if (verbose) {
    cat(
      'posterior density of samples',
      paste0(
        median(dpost_samples),
        ' [ ',
        quantile(dpost_samples, probs = 0.025),
        ' - ',
        quantile(dpost_samples, probs = 0.975) ,
        ' ] '
      ),
      '\n'
    )
    cat(
      'posterior density of 2* mu - samples',
      paste0(
        median(dpost_samples2),
        ' [ ',
        quantile(dpost_samples2, probs = 0.025),
        ' - ',
        quantile(dpost_samples2, probs = 0.975) ,
        ' ] '
      ),
      '\n'
    )
    cat(
      'posterior density of proposed samples (minus)',
      paste0(
        median(dpost_proposed_samples1),
        ' [ ',
        quantile(dpost_proposed_samples1, probs = 0.025),
        ' - ',
        quantile(dpost_proposed_samples1, probs = 0.975) ,
        ' ] '
      ),
      '\n'
    )
    cat(
      'posterior density of proposed samples (plus)',
      paste0(
        median(dpost_proposed_samples2),
        ' [ ',
        quantile(dpost_proposed_samples2, probs = 0.025),
        ' - ',
        quantile(dpost_proposed_samples2, probs = 0.975) ,
        ' ] '
      ),
      '\n'
    )
    cat(
      'proposal density of normalised samples',
      paste0(
        median(dprop_samples1),
        ' [ ',
        quantile(dprop_samples1, probs = 0.025),
        ' - ',
        quantile(dprop_samples1, probs = 0.975) ,
        ' ] '
      ),
      '\n'
    )
    cat(
      'proposal density of proposed samples',
      paste0(
        median(dprop_proposed_samples),
        ' [ ',
        quantile(dprop_proposed_samples, probs = 0.025),
        ' - ',
        quantile(dprop_proposed_samples, probs = 0.975) ,
        ' ] '
      ),
      '\n'
    )
  }
  
  ll1 = Rlogdet$modulus - log(2) + dpost_samples + dpost_samples2 - dprop_samples1
  ll2 = Rlogdet$modulus - log(2) + dpost_proposed_samples1 + dpost_proposed_samples2 - dprop_proposed_samples
  if (verbose) {
    cat('ll1',
        paste0(
          median(ll1),
          ' [ ',
          quantile(ll1, probs = 0.025),
          ' - ',
          quantile(ll1, probs = 0.975) ,
          ' ] '
        ),
        '\n')
    cat('ll2',
        paste0(
          median(ll2),
          ' [ ',
          quantile(ll2, probs = 0.025),
          ' - ',
          quantile(ll2, probs = 0.975) ,
          ' ] '
        ),
        '\n')
  }
  
  if (i == 2) {
    # star
    ll1star = median(ll1[ll1 > quantile(ll1, 0.025) &
                           ll1 < quantile(ll1, 0.975)])
    ll2star = median(ll2[ll2 > quantile(ll2, 0.025) &
                           ll2 < quantile(ll2, 0.975)])
    if (verbose) {
      cat('lstar values are  ', ll1star, ' and ', ll2star, ' \n')
    }
  }
  
  tmp1 = as.brob(ll2 - ll2star)
  tmp2 = as.brob(log_marginal_likelihood[i - 1] - ll2star)
  tmp3 = as.brob(ll1 - ll1star)
  tmp4 = as.brob(log_marginal_likelihood[i - 1] - ll1star)
  tmp5 = as.brob(-ll1star)
  tmp1 = exp(tmp1)
  tmp2 = exp(tmp2)
  tmp3 = exp(tmp3)
  tmp4 = exp(tmp4)
  tmp5 = exp(tmp5)
  tmp6 = sum(tmp1 / (s1 * tmp1 + s2 * tmp2)) / N2
  tmp7 = sum(tmp5 / (s1 * tmp3 + s2 * tmp4)) / N1
  log_marginal_likelihood[i] = log(tmp6) - log(tmp7)
  if (verbose) {
    cat(
      'the log marginal likelihood is ',
      (log(tmp6)),
      ' - ',
      log(tmp7),
      ' = ',
      log_marginal_likelihood[i],
      ' \n'
    )
  }
  gc()
}

time <- Sys.time() - ptm
save(log_marginal_likelihood,
     file = paste0('bs_warp', args[1], '_results_', randomnumber, '.rda'))

