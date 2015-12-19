# Authors: Maes Benjamin , Prévot Clément
# Date  : 2015/12/19
# 
# Calibration d'un modele a l'aide d'une methode d'optimisation
#

# Approximation lineaire des donnees
# data : data$V1 contient les diametres pupillaires pour les 25 premiers instants
# t    : temps (nombre reel) ou désire la valeur du diametre pupillaire 
#
# Valeur finale : fonction du diametre pupillaire 
linear_approx <- function(data) {
  f <- function(t) {
    i <- floor(t) + 1  # indice dans data (qui commence a l'indice 1)
    if (i <= 1)
      y <- data[1,]
    else
      if (i>=25)
        y <- data[25,]
      else
        # interpolation lineaire entre data[i] et data[i+1] 
        y <- data[i,] + (t - (i-1)) * (data[i + 1,] - data[i,])
  
    return(y)
  }
  
  return(f)
}

# Troisieme modele qui combine les donnees des taches simples
# dm    : diametre pupillaire lors de la digit memory task
# ws    : diametre pupillaire lors de la word seach task
# theta : parametres du modele
#
# Valeur finale : 
#     vecteur avec les diametres pupillaires issus du modele pour les instants de 0 a 24
model3 <- function(f_dm, f_ws, theta) {
  alpha1 <- theta[1]
  omega1 <- theta[2]
  T <- theta[3]
  alpha2 <- theta[4]
  omega2 <- (24 - omega1 * T)/(24 - T)
  k1 <- theta[5]
  k2 <- theta[6]
  
  vy <- c()
  for(t in 0:24) {
    if (t < T){
      y <- k1 * (alpha1 * f_dm(t) + (1 - alpha1) * f_ws(omega1 * t))
    }else{
      y <- k2 * (alpha2 * f_dm(t) + (1 - alpha2) * f_ws(omega1 * T + omega2 * (t - T)))
    }
    vy <- c(vy, y)
  }
  
  return(vy)
}

mean_square_error <- function(y1, y2) {
  return(mean((y1 - y2)^2))  
}

# Calcul de la qualite du modele 3: 
#    erreur quadratique moyenne par rapport aux donne de la tache dual
# dm   : digit memory dataframe 
# ws   : word search dataframe 
# dual : double task dataframe 
#
# Valeur finale : fonction qui calcule la mse en fonction de theta
mse_model3 <- function(dm, ws, dual) {
  f_dm <- linear_approx(dm)
  f_ws <- linear_approx(ws)

  mse <- function(theta) {
    # calcul de la sortie du modele
    y <- model3(f_dm, f_ws, theta)
    # erreur quadratique
    return(mean_square_error(y, dual$V1))
  }
  
  return(mse)
}


#==================================================================#
#        (1 + 1)-ES - Algorithme et analyse des performances       #
#==================================================================#

# Algorithme 
onePlusOneES <- function(f,x0,params){
  m <- x0
  fm <- f(m)
  nbEval <- 0
  
  sigma<-params$sigma
  nbEvalMax <- params$nbEvalMax
  
  fdyn<-c(fm)
  
  while (nbEval < nbEvalMax){ 
    x_prime <- m + sigma* rnorm(6)
    fx <- f(x_prime)
    if(fx < fm){
      m <- x_prime
      fm <- fx
    }
    fdyn <-c(fdyn,fm)
    nbEval <- nbEval +1
  } 
  
  return( list(x=m,fx=fm,fdyn=fdyn))
}

# Analyse des performances
analyseOnePlusOneES <-function(datas){
  
  # Mean square error function for the model3
  mse_5bfd <- mse_model3(datas$dm_5, datas$ws_bfd, datas$dual_5bfd)
  mse_5bff <- mse_model3(datas$dm_5, datas$ws_bff, datas$dual_5bff)
  mse_5hfd <- mse_model3(datas$dm_5, datas$ws_hfd, datas$dual_5hfd)
  mse_5hff <- mse_model3(datas$dm_5, datas$ws_hff, datas$dual_5hff)
  
  mse_9bfd <- mse_model3(datas$dm_9, datas$ws_bfd, datas$dual_9bfd)
  mse_9bff <- mse_model3(datas$dm_9, datas$ws_bff, datas$dual_9bff)
  mse_9hfd <- mse_model3(datas$dm_9, datas$ws_hfd, datas$dual_9hfd)
  mse_9hff <- mse_model3(datas$dm_9, datas$ws_hff, datas$dual_9hff)
  
  params <- list(nbEvalMax = 2000 , sigma = 1)
  
  res_5bfd <- onePlusOneES(f =mse_5bfd,x0 =  datas$dual_5bfd$V1, params = params)
  res_5bff <- onePlusOneES(f =mse_5bff,x0 =  datas$dual_5bff$V1, params = params)
  res_5hfd <- onePlusOneES(f =mse_5hfd,x0 =  datas$dual_5hfd$V1, params = params)
  res_5hff <- onePlusOneES(f =mse_5hff,x0 =  datas$dual_5hff$V1, params = params)
  
  plot(res_5hfd$fdyn, type="o",pch=30,main = "Evolution de l'erreur moyenne quadratique", col = "red", log = "yx",xlab = "nb Eval", ylab = "erreur moyenne quadratique")
  lines(res_5hff$fdyn, col = "blue")
  lines(res_5bfd$fdyn, col = "green")
  lines(res_5bff$fdyn, col = "orange")
  dev.print(jpeg,filename="result/OnePlusOne_result/mean_square_error_5.jpg",quality=100,units="px",width=2500,res=300)
  
  
  res_9bfd <- onePlusOneES(f =mse_9bfd,x0 =  datas$dual_9bfd$V1, params = params)
  res_9bff <- onePlusOneES(f =mse_9bff,x0 =  datas$dual_9bff$V1, params = params)
  res_9hfd <- onePlusOneES(f =mse_9hfd,x0 =  datas$dual_9hfd$V1, params = params)
  res_9hff <- onePlusOneES(f =mse_9hff,x0 =  datas$dual_9hff$V1, params = params)  
  
  plot(res_9hfd$fdyn, type="o",pch=30,main = "Evolution de l'erreur moyenne quadratique", col = "red", log = "yx",xlab = "nb Eval", ylab = "erreur moyenne quadratique")
  lines(res_9hff$fdyn, col = "blue")
  lines(res_9bfd$fdyn, col = "green")
  lines(res_9bff$fdyn, col = "orange")
  dev.print(jpeg,filename="result/OnePlusOne_result/mean_square_error_9.jpg",quality=100,units="px",width=2500,res=300)
}

#==================================================================#
# (1 + 1)-ES - Fifth rule - Algorithme et analyse des performances #
#==================================================================#

# Algorithme 
onePlusOneFifthES <- function(){
}

# Analyse des performances
analyseOnePlusOneFifthES <-function(datas){
  
}

#==================================================================#
#       (μ/μ, λ)-ES - Algorithme et analyse des performances       #
#==================================================================#

# Algorithme 
mumMuLambdaES <- function(){
}

# Analyse des performances
analyseMuMuLambdaES <-function(datas){
  
}

drawGraphics <- function(datas){
  plot(datas$dm_5$V1, type="o",main = "encod5-bf_d", col = "red", ylab = "diameter", xlab = "time")
  lines(datas$ws_bfd$V1, type = "o", col = "green")
  lines(datas$dual_5bfd$V1, type="o", col = "blue")
  dev.print(jpeg,filename="result/encod5-bf_d.jpg",quality=100,units="px",width=2500,res=300)
  
  plot(datas$dm_5$V1, type="o",main = "encod5-bf_f", col = "red", ylab = "diameter", xlab = "time")
  lines(datas$ws_bff$V1, type = "o", col = "green")
  lines(datas$dual_5bff$V1, type="o", col = "blue")
  dev.print(jpeg,filename="result/encod5-bf_f.jpg",quality=100,units="px",width=2500,res=300)
  
  plot(datas$dm_5$V1, type="o",main = "encod5-hf_d", col = "red", ylab = "diameter", xlab = "time")
  lines(datas$ws_hfd$V1, type = "o", col = "green")
  lines(datas$dual_5hfd$V1, type="o", col = "blue")
  dev.print(jpeg,filename="result/encod5-hf_d.jpg",quality=100,units="px",width=2500,res=300)
  
  plot(datas$dm_5$V1, type="o",main = "encod5-hf_f", col = "red", ylab = "diameter", xlab = "time")
  lines(datas$ws_hff$V1, type = "o", col = "green")
  lines(datas$dual_5hff$V1, type="o", col = "blue")
  dev.print(jpeg,filename="result/encod5-hf_f.jpg",quality=100,units="px",width=2500,res=300)
  
  plot(datas$dm_9$V1, type="o",main = "encod9-bf_d", col = "red", ylab = "diameter", xlab = "time")
  lines(datas$ws_bfd$V1, type = "o", col = "green")
  lines(datas$dual_9bfd$V1, type="o", col = "blue")
  dev.print(jpeg,filename="result/encod9-bf_d.jpg",quality=100,units="px",width=2500,res=300)
  
  plot(datas$dm_9$V1, type="o",main = "encod9-bf_f", col = "red", ylab = "diameter", xlab = "time")
  lines(datas$ws_bff$V1, type = "o", col = "green")
  lines(datas$dual_9bff$V1, type="o", col = "blue")
  dev.print(jpeg,filename="result/encod9-bf_f.jpg",quality=100,units="px",width=2500,res=300)
  
  plot(datas$dm_9$V1, type="o",pch="46",main = "encod9-hf_d", col = "red", ylab = "diameter", xlab = "time")
  lines(datas$ws_hfd$V1, type = "o", col = "green")
  lines(datas$dual_9hfd$V1, type="o", col = "blue")
  dev.print(jpeg,filename="result/encod9-hf_d.jpg",quality=100,units="px",width=2500,res=300)
  
  plot(datas$dm_9$V1, type="o",main = "encod9-hf_f", col = "red", ylab = "diameter", xlab = "time")
  lines(datas$ws_hff$V1, type = "o", col = "green")
  lines(datas$dual_9hff$V1, type="o", col = "blue")
  dev.print(jpeg,filename="result/encod9-hf_f.jpg",quality=100,units="px",width=2500,res=300)
}

#==================================================================#
#                             MAIN                                 #
#==================================================================#
main <- function() {
  # chargement des donnees
  dm_5 <- read.table("data/digitMem_encod5.dat") 
  dm_9 <- read.table("data/digitMem_encod9.dat") 
  
  ws_hfd <- read.table("data/wordSearch_hf_d.dat")
  ws_hff <- read.table("data/wordSearch_hf_f.dat")
  ws_bfd <- read.table("data/wordSearch_bf_d.dat")
  ws_bff <- read.table("data/wordSearch_bf_f.dat")
  
  dual_5bfd <- read.table("data/dualTask_encod5-bf_d.dat")
  dual_5bff <- read.table("data/dualTask_encod5-bf_f.dat")
  dual_5hfd <- read.table("data/dualTask_encod5-hf_d.dat")
  dual_5hff <- read.table("data/dualTask_encod5-hf_f.dat")
  
  dual_9bfd <- read.table("data/dualTask_encod9-bf_d.dat")
  dual_9bff <- read.table("data/dualTask_encod9-bf_f.dat")
  dual_9hfd <- read.table("data/dualTask_encod9-hf_d.dat")
  dual_9hff <- read.table("data/dualTask_encod9-hf_f.dat")
  #====================================================================================}#
  
  # les differents données
  datas <- list(dm_5=dm_5, dm_9=dm_9, ws_bff=ws_bff, ws_bfd=ws_bfd, ws_hff=ws_hff, ws_hfd=ws_hfd, dual_5hff=dual_5hff, dual_5hfd=dual_5hfd, dual_5bff=dual_5bff, dual_5bfd=dual_5bfd, dual_9hff=dual_9hff, dual_9hfd=dual_9hfd, dual_9bff=dual_9bff, dual_9bfd=dual_9bfd)
  
  # traçage des graphique basique  -  diametre de la pupille en fonction du temps #
  drawGraphics(datas)
  
  analyseOnePlusOneES(datas)
  
  #analyseOnePlusOneFifthES(datas)
  
  #analyseMuMuLambda(datas)
  
}

main()