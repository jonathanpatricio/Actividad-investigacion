Pob <- function(b, sd, tipo, prob, N, dic){
  # Semilla para fijar los números aleatorios
  set.seed(16)
  
  # Librerías que utiliza la función
  library(dplyr)
  
  # Parámetros del modelo
  parametros <- cbind(b,sd,tipo,prob)
  
  # Matriz que guardará el resultado de las covariables
  covariables <- matrix(0, N, length(b))
  
  # Generando los valores de cada sujeto en la población
  for (i in 1:length(b)) {
    covariables[,i] <- if(parametros[i,3] == 1){rnorm(n = N, mean = parametros[i,1], sd = parametros[i,2])}else
    {(rbinom(n = N, size = 1, prob = parametros[i,4])) * (parametros[i,1]) }
  }
  
  # Obteniendo el valor esperado para cada sujeto
  mu <- exp(as.matrix(apply(X = covariables, MARGIN = 1, FUN = sum)))
  
  # Obteniendo la información de conteo para cada sujeto
  y <- rpois(n = N, lambda = mu)
  
  # Construyendo la base de datos 
  base <-  as.data.frame(cbind(covariables, y))
  
  # Dicotomizando la variable dependiente
  base$y_dic <- if_else(condition = base$y < dic, true = 0, false = 1)
  base_completa <- base
  base <- base[,-c(length(b)+1)]
  return(list(base_completa = base_completa, base = base))
}

poblacion <- Pob(b  = c(0.987, 0.097, 0.021, 0.360, 0.173, -0.010, -0.002, 0.231, 0.089, 0.004, 0.009, 0.012), 
                 sd = c(0.023, 0.038, 0.039, 0.033, 0.029,  0.001,  0.046, 0.046, 0.040, 0.001, 0.001, 0.001), 
                 tipo = c(1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1 ), 
                 prob = c(0, 0.208, 0.208, 0.259, 0.613, 0, 0, 0, 0, 0, 0, 0), 
                 dic = 1,
                 N  = 5670)

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
  library(cutpointr)
  library(caret)
  library(rms)
  library(generalhoslem)
  
  # Matrices que guardarán los resultados de las comparaciones
  R2_Poisson   <- matrix(0,length(n),repeticiones)
  R2_Logistico <- matrix(0,length(n),repeticiones)
  
  Brier_Poisson   <- matrix(0,length(n),repeticiones)
  Brier_Logistico <- matrix(0,length(n),repeticiones)
  
  AUC_Poisson   <- matrix(0,length(n),repeticiones)
  AUC_Logistico <- matrix(0,length(n),repeticiones)
  
  Se_Poisson   <- matrix(0,length(n),repeticiones)
  Se_Logistico <- matrix(0,length(n),repeticiones)
  
  Es_Poisson   <- matrix(0,length(n),repeticiones)
  Es_Logistico <- matrix(0,length(n),repeticiones)
  
  VPP_Poisson   <- matrix(0,length(n),repeticiones)
  VPP_Logistico <- matrix(0,length(n),repeticiones)
  
  VPN_Poisson   <- matrix(0,length(n),repeticiones)
  VPN_Logistico <- matrix(0,length(n),repeticiones)
  
  HL_test_poisson  <- matrix(0,length(n),repeticiones)
  HL_test_logistic <- matrix(0,length(n),repeticiones)
  
  
  # Bucle
  for (i in 1:length(n)){for (j in 1:repeticiones) {

      # Obteniendo la muestra
      sample <- base[sample(x = 1:nrow(base), size = n[[i]], replace = FALSE),]
      
      # Ajustando los modelos
      Poisson   <- glm(formula = y_dic ~ ., data = sample, family = "poisson")
      std.err   <- sqrt(diag(vcovHC(Poisson, type="HC0"))) # Errores estándar robustos del modelo de Poisson
      Logistico <- glm(formula = y_dic ~ ., data = sample, family = binomial(link = "logit"))
      
      # Evaluando la discriminación de los modelos
      sample$Poisson   <- predict(Poisson,   type = c("response")) # Obteniendo las predicciones para el modelo Poisson
      sample$Logistico <- predict(Logistico, type = c("response")) # Obteniendo las predicciones para el modelo Logístico 
      
      cut_point_Poisson <- cutpointr(sample,                     # base de datos con la que se construyó el modelo 
                                     Poisson,                    # variable que guarda las predicciones del modelo
                                     y_dic,                      # variable dependiente convertida como numérica
                                     direction = ">=", 
                                     pos_class = 1,              # valores en mi variable observada un caso positivo
                                     neg_class = 0,              # valores en mi variable observada un caso negativo
                                     method = maximize_metric,   # Método para maximizar (se)
                                     metric = youden)            # youden = sensitivity + specificity - 1
      
      cut_point_Logistico <- cutpointr(sample,                   # base de datos con la que se construyó el modelo 
                                       Logistico,                # variable que guarda las predicciones del modelo
                                       y_dic,                    # variable dependiente convertida como numérica
                                       direction = ">=", 
                                       pos_class = 1,            # valores en mi variable observada un caso positivo
                                       neg_class = 0,            # valores en mi variable observada un caso negativo
                                       method = maximize_metric, # Metodo para maximizar (se)
                                       metric = youden)          # youden = sensitivity + specificity - 1
        
      
      # AUC Roc curve
      AUC_Poisson   [i,j] <- cut_point_Poisson$AUC
      AUC_Logistico [i,j] <- cut_point_Logistico$AUC
      
      # Resultados de la matriz de confusión 
      sample$Poisson_dic   <- as.numeric(if_else(condition = sample$Poisson   >= cut_point_Poisson$optimal_cutpoint, true = 1, false = 0 ))
      sample$Logistico_dic <- as.numeric(if_else(condition = sample$Logistico >= cut_point_Logistico$optimal_cutpoint, true = 1, false = 0 ))
      
      Conf_matriz_poisson  <-  confusionMatrix(as.factor(sample$Poisson_dic),   as.factor(sample$y_dic), positive = "1")
      Conf_matriz_logistic <-  confusionMatrix(as.factor(sample$Logistico_dic), as.factor(sample$y_dic), positive = "1")
      
      Se_Poisson    [i,j] <- as.matrix(Conf_matriz_poisson$byClass)[1,1]
      Se_Logistico  [i,j] <- as.matrix(Conf_matriz_logistic$byClass)[1,1]
      
      Es_Poisson    [i,j] <- as.matrix(Conf_matriz_poisson$byClass)[2,1]
      Es_Logistico  [i,j] <- as.matrix(Conf_matriz_logistic$byClass)[2,1]
      
      VPP_Poisson   [i,j] <- as.matrix(Conf_matriz_poisson$byClass)[3,1]
      VPP_Logistico [i,j] <- as.matrix(Conf_matriz_logistic$byClass)[3,1]
      
      VPN_Poisson   [i,j] <- as.matrix(Conf_matriz_poisson$byClass)[4,1]
      VPN_Logistico [i,j] <- as.matrix(Conf_matriz_logistic$byClass)[4,1]
      
      #Evaluando el rendimiento global de los modelos
      
      # R2 Nagelkerke obtenidos en ambos modelos 
      R2_Poisson   [i,j] <- PseudoR2(Poisson,   which = "Nagelkerke")
      R2_Logistico [i,j] <- PseudoR2(Logistico, which = "Nagelkerke")
      
      # Brier score
      Brier_Poisson   [i,j] <- BrierScore(resp = sample$y_dic, pred = sample$Poisson, scaled = TRUE)
      Brier_Logistico [i,j] <- BrierScore(resp = sample$y_dic, pred = sample$Logistico, scaled = TRUE)
      
      # Evaluando la calibración en ambos modelos
      HL_test_poisson  [i,j] <- if((hoslem.test(sample$y_dic, fitted(Poisson),   g = 10)$p.value) < 0.05 ) 1 else 0 
      HL_test_logistic [i,j] <- if((hoslem.test(sample$y_dic, fitted(Logistico),   g = 10)$p.value) < 0.05 ) 1 else 0
    
      
      # Capturando el avance del proceso
      print(c(n[[i]],j))
      
    }
    
  }
  
  # Obteniendo los valores promedios de las comparaciones
  R2_Poisson    <- as.matrix(apply(X = R2_Poisson,   MARGIN = 1, FUN = mean))
  R2_Logistico  <- as.matrix(apply(X = R2_Logistico, MARGIN = 1, FUN = mean))
  
  Brier_Poisson    <- as.matrix(apply(X = Brier_Poisson,   MARGIN = 1, FUN = mean))
  Brier_Logistico  <- as.matrix(apply(X = Brier_Logistico, MARGIN = 1, FUN = mean))
  
  AUC_Poisson   <- as.matrix(apply(X = AUC_Poisson,  MARGIN = 1, FUN = mean))
  AUC_Logistico <- as.matrix(apply(X = AUC_Logistico,MARGIN = 1, FUN = mean))
  
  Se_Poisson    <- as.matrix(apply(X = Se_Poisson,   MARGIN = 1, FUN = mean))
  Se_Logistico  <- as.matrix(apply(X = Se_Logistico, MARGIN = 1, FUN = mean))
  
  Es_Poisson    <- as.matrix(apply(X = Es_Poisson,   MARGIN = 1, FUN = mean))
  Es_Logistico  <- as.matrix(apply(X = Es_Logistico, MARGIN = 1, FUN = mean))
  
  VPP_Poisson   <- as.matrix(apply(X = VPP_Poisson,  MARGIN = 1, FUN = mean))
  VPP_Logistico <- as.matrix(apply(X = VPP_Logistico,MARGIN = 1, FUN = mean))
  
  VPN_Poisson   <- as.matrix(apply(X = VPN_Poisson,  MARGIN = 1, FUN = mean))
  VPN_Logistico <- as.matrix(apply(X = VPN_Logistico,MARGIN = 1, FUN = mean))
  
  HL_test_poisson   <- as.matrix(apply(X = HL_test_poisson,  MARGIN = 1, FUN = mean))
  HL_test_logistic  <- as.matrix(apply(X = HL_test_logistic, MARGIN = 1, FUN = mean))

  # Base de datos con los resultados
  Base <- as.data.frame(cbind(n = n))
  Base <- mutate(Base,R2_Poisson)
  Base <- mutate(Base,R2_Logistico)
  
  Base <- mutate(Base,Brier_Poisson)
  Base <- mutate(Base,Brier_Logistico)
  
  Base <- mutate(Base,AUC_Poisson)
  Base <- mutate(Base,AUC_Logistico)
  
  Base <- mutate(Base,Se_Poisson)
  Base <- mutate(Base,Se_Logistico)
  
  Base <- mutate(Base,Es_Poisson)
  Base <- mutate(Base,Es_Logistico)
  
  Base <- mutate(Base,VPP_Poisson)
  Base <- mutate(Base,VPP_Logistico)
  
  Base <- mutate(Base,VPN_Poisson)
  Base <- mutate(Base,VPN_Logistico)
  
  Base <- mutate(Base,HL_test_poisson)
  Base <- mutate(Base,HL_test_logistic)
  
  # Capturando la hora de término del la función
  Fin <- DescTools::Now()
  
  # Calculando la duración
  Duración <- Fin - Inicio; print(Duración)
  
  # Retornando los resultados
  return(Resultados = Base)
}

Comparaciones <- comp_modelos(base = base, n = c(1000, 1500), repeticiones = 100)
