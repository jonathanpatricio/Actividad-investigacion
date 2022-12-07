#Librerías a utilizar
library(readxl)
library(sandwich)
library(dplyr)
library(stargazer)

#Leyendo la base de datos
datos <- read_excel("examen.xls")

datos$sexo <- factor(datos$sexo, labels = c("Masculino", "Femenino"))
datos$educmad <- factor(datos$educmad, labels = c("12 años o menos","Más de 12 años"))
datos$raza <- factor(datos$raza, labels = c("Negra", "Blanca"))
datos$estrato <- factor(datos$estrato, labels = c("Estratos 1 y 2","Estratos 3 y 4","Estratos 5 y 6"))
datos$estcivil <- factor(datos$estcivil, labels = c("Soltero","Casado/Unión Libre","Separado/Divorciado"))
datos$lesion <- factor(datos$lesion, labels = c("No", "Si"))
datos$numerolesiones_2 <- if_else(condition = datos$numerolesiones < 1, true = 0, false = 1)

# Ajustando el modelo de poisson

m1 <- glm(formula = numerolesiones ~ educmad + sexo + estcivil + estrato, data = datos, family="poisson")
summary(m1)
# m2 <- glm(formula = numerolesiones_2 ~ educmad + sexo + estcivil + estrato, data = datos, family="poisson")
cov.m1 <- vcovHC(m1, type="HC0")
std.err <- sqrt(diag(cov.m1))
# stargazer(m1,m2, type = "text")

datos$fitted_1 <- predict(m1,type=c("response"), dispersion = std.err)
datos$fitted_2 <- predict(m1,type=c("response"))




# Ajuntando el modelo logistico












