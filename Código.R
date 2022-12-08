Pob <- function(b, sd, tipo, prob, N){
  # Semilla para fijar los números aleatorios
  set.seed(16)
  
  # Librerías que utiliza la función
  library(dplyr)
  
  # Parámetros del modelo
  parametros <- cbind(b,sd,tipo,prob)
  
  # Matriz que guardará el resultado de las covariables
  covariables <- matrix(0,N,length(b))
    
  # Generando los valores de cada sujeto en la población
  for (i in 1:length(b)) {
    covariables[,i] <- if(parametros[i,3] == 1){rnorm(n = N, mean = parametros[i,1], sd = parametros[i,2])}else{(rbinom(n = N, size = 1, prob = parametros[i,4])) * (sd = parametros[i,4]) }
  }
  
  # Obteniendo el valor esparado de cada sujeto
  mu <- exp(as.matrix(apply(X = covariables, MARGIN = 1, FUN = sum)))
  
  # Obteniendo la información de conteo para cada sujeto
  y <- rpois(n = N, lambda = mu)

  #Construyendo la base de datos 
  base <-  as.data.frame(cbind(y, covariables))
  
  # Dicotomizando la variable dependiente
  base$y_dic <- if_else(condition = base$y < 1, true = 0, false = 1)
  
  return(base)
  
}

base <- Pob(b  = c(-0.55608, 0.37697, 0.10030, 0.05444), # Vector que contiene los valores con el promedio de cada una de las covariables (variables numericas)
            sd = c(0.16204, 0.15148, 0.02275, 0.01349),  # Vector que contiene los valores con la desviación estandar para cada variable (variables numericas)
            tipo = c(1,1,1,0),                           # Vector que indica con uno (1) y cero (0) si la variables es numerica (1) o si es dicotoma (0) 
            prob = c(0,0,0,0.25),                        # Vector que contiene las probabilidades para las variables dicotomicas
            N = 10000)                                   # Tamaño de la población a simular





m1 <- glm(formula = numerolesiones ~ educmad + sexo + estcivil + estrato, data = datos, family="poisson")
summary(m1)
# m2 <- glm(formula = numerolesiones_2 ~ educmad + sexo + estcivil + estrato, data = datos, family="poisson")
cov.m1 <- vcovHC(m1, type="HC0")
std.err <- sqrt(diag(cov.m1))
# stargazer(m1,m2, type = "text")

datos$fitted_1 <- predict(m1,type=c("response"), dispersion = std.err)
datos$fitted_2 <- predict(m1,type=c("response"))




# Ajuntando el modelo logistico












