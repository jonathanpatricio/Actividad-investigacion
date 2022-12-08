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


comp_modelos <- function(base, n, repeticiones){
  # Semilla para fijar los resultados
  set.seed(16)
  
  # Capturando la hora de inicio del la función
  Inicio <- DescTools::Now()
  
  # Librerías que utiliza la función
  library(sandwich)
  library(DescTools)
  library(doParallel)
  library(parallel)
  library(cutpointr)
  
  # Matrices que guardarán los resultados de las hipótesis planteadas para 2 y 3 tratamientos
  R2_Poisson   <- matrix(0,length(n),repeticiones)
  R2_Logistico <- matrix(0,length(n),repeticiones)
  
  Brier_Poisson   <- matrix(0,length(n),repeticiones)
  Brier_Logistico <- matrix(0,length(n),repeticiones)
  
  AUC_Poisson   <- matrix(0,length(n),repeticiones)
  AUC_Logistico <- matrix(0,length(n),repeticiones)
  
  
  
  # Bucle
  for (i in 1:length(n)){for (j in 1:repeticiones) {

      # Obteniendo la muestra
      sample <- base[sample(x = 1:nrow(base), size = n[[i]], replace = FALSE),]
      
      # Ajustando los modelos
      Poisson   <- glm(formula = y_dic ~ ., data = sample, family = "poisson")
      std.err   <- sqrt(diag(vcovHC(Poisson, type="HC0"))) # Errores estándar robustos del modelo de Poisson
      Logistico <- glm(formula = y_dic ~ ., data = sample, family = binomial(link = "logit"))
      
      
      
      #Evaluando el rendimiento global de los modelos
      # R2 Nagelkerke obtenidos en ambos modelos 
      R2_Poisson   [i,j] <- PseudoR2(Poisson,   which = "Nagelkerke")
      R2_Logistico [i,j] <- PseudoR2(Logistico, which = "Nagelkerke")
      
      # Brier score
      Brier_Poisson   [i,j] <- 0
      Brier_Logistico [i,j] <- 0
      
      # Evaluando la discriminación de los modelos
      sample$Poisson <- predict(Poisson, type=c("response"))
      sample$Logistico <- predict(Logistico, type=c("response"))
      
      cut_point_Poisson <- cutpointr(sample,                     # base de datos con la que trabajo
                                     Poisson,                    #variable que guarda las predicciones del modelo
                                     y_dic,                      # variable dependiente convertida como numérica
                                     direction = ">=", 
                                     pos_class = 0,              # valores en mi variable observada un caso positivo
                                     neg_class = 1,              # valores en mi variable observada un caso negativo
                                     method = maximize_metric,   # Metodo para maximizar (se)
                                     metric = youden)            # youden = sensitivity + specificity - 1
      
      cut_point_Logistico <- cutpointr(sample,                   # base de datos con la que trabajo
                                       Logistico,                #variable que guarda las predicciones del modelo
                                       y_dic,                    # variable dependiente convertida como numérica
                                       direction = ">=", 
                                       pos_class = 0,            # valores en mi variable observada un caso positivo
                                       neg_class = 1,            # valores en mi variable observada un caso negativo
                                       method = maximize_metric, # Metodo para maximizar (se)
                                       metric = youden)          # youden = sensitivity + specificity - 1
        
      # AUC Roc curve
      AUC_Poisson   [i,j] <- cut_point_Poisson$AUC
      AUC_Logistico [i,j] <- cut_point_Logistico$AUC
    
      print(c(n[[i]],j))
      
    }
    
  }
  
  #Base de datos para 2 tratamientos
  R2_Poisson    <- as.matrix(apply(X = R2_Poisson,      MARGIN = 1, FUN = mean))
  R2_Logistico  <- as.matrix(apply(X = R2_Logistico,    MARGIN = 1, FUN = mean))
  
  AUC_Poisson   <- as.matrix(apply(X = AUC_Poisson,     MARGIN = 1, FUN = mean))
  AUC_Logistico <- as.matrix(apply(X = AUC_Logistico,   MARGIN = 1, FUN = mean))

  
  Base <- as.data.frame(cbind(n = n))
  
  Base <- mutate(Base,R2_Poisson)
  Base <- mutate(Base,R2_Logistico)
  Base <- mutate(Base,AUC_Poisson)
  Base <- mutate(Base,AUC_Logistico)

  # Capturando la hora de término del la función
  Fin <- DescTools::Now()
  
  # Calculando la duración
  Duración <- Fin - Inicio; print(Duración)
  
  # Retornando los resultados
  return(Resultados = Base)
  
}

Comparaciones <- comp_modelos(base = base, n = c(1000), repeticiones = 100)
