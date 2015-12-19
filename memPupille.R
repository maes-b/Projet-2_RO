# Author: Sebastien Verel
# Date  : 2016/12/14
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
    if (t < T)
      y <- k1 * (alpha1 * f_dm(t) + (1 - alpha1) * f_ws(omega1 * t))
    else
      y <- k2 * (alpha2 * f_dm(t) + (1 - alpha2) * f_ws(omega1 * T + omega2 * (t - T)))
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

main <- function() {
  # chargement des donnees
  dm <- read.table("../data/digitMem_encod5.dat")
  ws <- read.table("../data/wordSearch_bf_f.dat")
  dual <- read.table("../data/dualTask_encod5-bf_f.dat")

  # Tracer des courbes experimentales
  plot(dm$V1, type="o", col = "red", ylab = "diameter", xlab = "time")
  lines(ws$V1, type = "o", col = "green")
  lines(dual$V1, type="o", col = "blue")  
  
  # test du modele 3 avec des parametres trouves "manuellement"
  f_dm <- linear_approx(dm)
  f_ws <- linear_approx(ws)  
  y <- model3(f_dm, f_ws, c(0.28, 2.1, 7, 2.1, 1, 1))
  lines(y, col="black")
  
  # Mean square error function for the model3
  mse <- mse_model3(dm, ws, dual)
  
  mse(c(0.28, 2.1, 7, 2.1, 1, 1))
}


