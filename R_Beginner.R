#    Name: Jonathan Medina
#    Purpose: Curso Master en RStudio: de principiante a experto

#    To: My dear students
#    This program is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 2 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

# PIDIENDO AYUDA EN RSTUDIO #

help.start()
?print
for (i in 1:8){
print(i)
}
1:8
length(1:8)

help()
help('+')

apropos("class")
val = 1:8
class(val)

example("read.table")
example("read.csv")

??regresion # == help.search("regresion")
??print # == help.search("print")

??regression # == help.search()
help.search("regression") # == ??

library(help = "stats")

vignette()
vignette("withr")

# OBJETOS TIPOS DE DATOS Y OPERACIONES BASICAS #

3+7
3-7

variable_mia <- (180 / 6) - 15
180/6
30-15

variable_mia

(180 / 6) - 15 -> variable_mia

z <- c(561, 1025, 1729, 2485, 2830)
print(z) # impresion explicita
z # impresion implicita

1:100
seq(from = 1, to = 100, by = 1) # == 1:100

c(1.1, 2.2, 5.5, 4.8) - c(1, 1, 1, 1) # 1i + 1j + 1k + 1l

c(1.1, 2.2, 5.5, 4.8) - c(1, 3) # 1.1 - 1, 2.2 - 3, 5.5 - 1, 4.8 - 3

ls()
# character == cadena de caracteres

"Hola Mundo!" # Cadena de caracteres o character

colores <- c("rojo", "verde", "azul", "azul", "rojo")
colores 

complejo <- 2 + 1i # Numero complejo tiene parte real y parte imaginaria

variable_mia == 15 # TRUE 
variable_mia == 16 # FALSE

length(complejo)
2 # 2 + 0i = 2 + 0 = 2
length(2)
2 + 1i # Longitud se convierte en 1
length(complejo) <- 3
complejo
complejo[1]
complejo[2]
complejo[2] <- 3 + 4i
complejo
complejo[3] <- 5 + 7i # cplx (complex) == complejo

2^1024 # Grande, tiende al infinito
-2^1024 # Grande, pero es negativo

# Todo numero dividido para cero es igual a infinito
# 0/0 ?
0/0 # NaN, Non Adding value

Sys.Date()
fecha_curso <- Sys.Date()
fecha_curso
class(fecha_curso) # No es character

valores <- vector(mode = "numeric", length = 10)
class(valores)

vector(mode = "character", length = 10)

vector(mode = "logical", length = 10)

list(0, "Hola", TRUE)

c(T, 19, 1 + 3i) # T == TRUE, 19 = 19.0, T == 1 y F == 0, T = 1 + 0i; 19 + 0i

valores <- vector(mode = "numeric", length = 5)
valores

as.logical(valores)

valores <- as.logical(valores)
valores

m <- matrix(data = 1:12, nrow = 4, ncol = 3) # 4*3 = 12
m
1:12

n <- matrix(data = 1:12, nrow = 4, ncol = 3, byrow = TRUE) # byrow = T, TRUE
n

?matrix

colores
class(colores)
factor(colores)

data.frame(llave = z, color = colores) # NA 

# SUBCONJUNTOS DE DATOS #

mi_vector <- 11:30
11:30 # seq(from = 11, to = 30, by = 1)
length(11:30)
mi_vector # Forma implicita
print(mi_vector)

mi_vector[3] # Posicion
rev(mi_vector)[3]
mi_vector[-3]
mi_vector
mi_vector[-c(3:5)]

mi_vector[1]
mi_vector[1:5]
mi_vector[c(4, 6, 13)] # c = combine, concatenate; ascendente y en orden

mi_vector[c(6, 13, 4)]
mi_vector[[3]] # == mi_vector[3]
mi_vector[[3]] == mi_vector[3] # Resultado = 13

mi_vector[rep(c(TRUE, FALSE), 10)]

mi_vector[c(FALSE, FALSE, TRUE)] # Aquellos valores cuyas posiciones sean multiplos de 3

mi_vector > 20
mi_vector[mi_vector > 20]

mi_arreglo <- array(seq(from = 1, to = 18, by = 1), dim = c(3, 3, 2)) # seq()
mi_arreglo

mi_arreglo[1, 3, 2] # Primer valor corresponde a la fila; segundo valor corresponde a la columna; tercer valor corresponde al arreglo

mi_arreglo[1:2, 1:2, 1]

mi_matriz <- matrix(data = 1:9, ncol = 3, nrow = 3) # data = seq(from = 1, to = 9, by = 1)
mi_matriz

mi_matriz[1, ]
mi_matriz[ ,1]

mi_matriz[2:3, ] # Primer espacio corresponde a las filas, segundo espacio corresponde a las columnas

mi_matriz[c(1, 3), ]

carro <- list(color = "rojo", nllantas = 4, marca = "Renault", ncil = 4)

carro$color # "rojo"
carro[c("ncil", "nllantas")]
carro$marca
carro[["marca"]]

carro[["mar", exact = FALSE]]

camioneta <- list(color = "azul", nllantas = 4, marca = "BMW", ncil = 6)

cochera <- list(carro, camioneta)
cochera

cochera[[c(2, 1)]]

# mph = miles per hour; km per hour 
data("cars")
View(cars)
cars$speed
length(cars$speed)

cars$dist>100
cars$dist[cars$dist>100]
cars$speed[cars$dist>100]

# IMPORTACION DE DATOS #
# 1ra Forma
Base_Datos <- read.table(file = "C:\\Users\\USUARIO\\OneDrive - Escuela Superior Politécnica del Litoral\\Archivos adjuntos\\R\\Libros\\1er libro leido\\RBook\\Vegetation2.txt",
                         header = TRUE) # , sep = ","
names(Base_Datos)
str(Base_Datos)

# 2da Forma
setwd("C:\\Users\\USUARIO\\OneDrive - Escuela Superior Politécnica del Litoral\\Archivos adjuntos\\R\\Libros\\1er libro leido\\RBook\\")
Datos2 <- read.table(file = "Vegetation2.txt", header = TRUE)

m <- mean(Datos2$R)
m
m1 <- mean(Datos2$R[Datos2$Transect == 1])
m1
Datos2$R[Datos2$Transect == 1]
m2 <- mean(Datos2$R[Datos2$Transect == 2])
m2
m3 <- mean(Datos2$R[Datos2$Transect == 3])
m3
unique(Datos2$Transect)
m4 <- mean(Datos2$R[Datos2$Transect == 4])
m4
m5 <- mean(Datos2$R[Datos2$Transect == 5])
m5
m6 <- mean(Datos2$R[Datos2$Transect == 6])
m6
m7 <- mean(Datos2$R[Datos2$Transect == 7])
m7
m8 <- mean(Datos2$R[Datos2$Transect == 8])
m8
c(m1, m2, m3, m4, m5, m6, m7, m8)
# La funcion que se denomina tapply
tapply(Datos2$R, Datos2$Transect, mean)
tapply(X = Datos2$R, INDEX = Datos2$Transect, FUN = mean)

Me <- tapply(X = Base_Datos$R, INDEX = Base_Datos$Transect, FUN = mean)
Me
Sd <- tapply(X = Base_Datos$R, INDEX = Base_Datos$Transect, FUN = sd)
Sd
Le <- tapply(X = Base_Datos$R, INDEX = Base_Datos$Transect, FUN = length)
Le
cbind(Me, Sd, Le)

# Lapply y Sapply
# Sapply
sapply(Base_Datos[ , 5:9], FUN = mean)
lapply(Base_Datos[, 5:9], FUN = mean)

# Artificio de datos para visualizar un data frame con sapply
sapply(cbind(Base_Datos$R, Base_Datos$ROCK, Base_Datos$LITTER, Base_Datos$ML, Base_Datos$BARESOIL), FUN = mean)
sapply(data.frame(cbind(Base_Datos$R, Base_Datos$ROCK, Base_Datos$LITTER, Base_Datos$ML, Base_Datos$BARESOIL)), FUN = mean)

# Exportacion de Datos 1ra Parte 
setwd("C:/Users/USUARIO/OneDrive - Escuela Superior Politécnica del Litoral/Archivos adjuntos/R/Libros/1er libro leido/RBook/")
Squid <- read.table(file = "squidGSI.txt", header = TRUE, sep = "")
?read.table

names(Squid)
colnames(Squid) # colnames() == names()

str(Squid)

# Manera erronea de importar
Squid <- read.table(file = "squidGSI.txt", header = TRUE, sep = ",")

mean(GSI, data = Squid) # media muestral o promedio, m = 1, 2, 3, 4; media = (1+2+3+4)/n = 4
mean(Squid$GSI)
boxplot(GSI ~ factor(Location), data = Squid)
boxplot(GSI ~ Location, data = Squid)

Squid$Sex
Squid$Month
max(Squid$Month)
Sel <- Squid$Sex == 1
unique(Squid$Sex)
Sel
SquidF <- Squid[Sel, ]

SquidM <- Squid[Squid$Sex == 2, ]
SquidM2 <- Squid[Squid$Sex == 2, c("Month", "Location", "Sex")]
SquidM2 <- Squid[Squid$Sex == 2, c(3, 4, 5)]

SquidF.OR.1 <- Squid[Squid$Sex == 1 & Squid$Location == 1, ] # & == y
unique(Squid$Location)

SquidF <- Squid[Squid$Sex == 1, ] # Squid[Sel, ], Sel <- Squid$Sex == 1
SquidF1 <- SquidF[SquidF$Location == 1, ]

Squid$Month[1:50]
unique(Squid$Month)
order(unique(Squid$Month))
order(Squid$Month[1:50], decreasing = TRUE)
?order

Ord1 <- order(Squid$Month)
Squid[Ord1, ]
View(Squid[Ord1, ])

Squid$Location <- factor(Squid$Location)
class(Squid$Location)
unique(Squid$Location)
unique(Squid$Sex)
Squid$fSex <- factor(Squid$Sex)
unique(Squid$fSex)
Squid$fSex <- factor(Squid$Sex, labels = c("M", "F"))
unique(Squid$fSex)

Squid$fLocation <- factor(Squid$Location,
                          levels = c(2, 3, 1, 4))
unique(Squid$fLocation)

## Exportacion
Sq1 <- read.table(file = "squid1.txt", header = TRUE)
Sq2 <- read.table(file = "squid2.txt", header = TRUE)
SquidMerged <- merge(Sq1, Sq2, by = "Sample")
?merge
SquidMerged <- merge(Sq1, Sq2, by = "Sample", all = TRUE) # all = FALSE

SquidM <- Squid[Squid$Sex == 1, ]
write.table(SquidM, 
            file = "MaleSquid.txt",
            sep = " ", quote = FALSE, append = FALSE, na = "NA")

write.table(SquidM, 
            file = "MaleSquid.txt",
            sep = " ", quote = TRUE, append = FALSE, na = "NA")

write.table(SquidM, 
            file = "MaleSquid.txt",
            sep = ",", quote = TRUE, append = FALSE, na = "NA")

# La funcion Plot
setwd("C:/Users/Toshiba/OneDrive - Escuela Superior Politécnica del Litoral/Archivos adjuntos/R/Libros/1er libro leido/RBook/")
Veg <- read.table(file = "Vegetation2.txt", header = TRUE)
plot(Veg$BARESOIL, Veg$R)
?plot
plot(Veg$R, Veg$BARESOIL)
plot(x = Veg$BARESOIL, y = Veg$R) # Grafica de dispersion
plot(y = Veg$R, x = Veg$BARESOIL)
plot(R ~ BARESOIL, data = Veg)
plot(x = Veg$BARESOIL, y = Veg$R, xlab = "Exposed Soil", ylab = "Species richness", main = "Sactter plot",
     xlim = c(0, 45), ylim = c(4, 19))
xlim <- c(min(Veg$BARESOIL), max(Veg$BARESOIL), na.rm = TRUE) # na.rm = remove nA's
min(Veg$BARESOIL)
max(Veg$BARESOIL)
NAVA <- vector()
n <- 1 # 2
dim(Veg)[2]
for (i in 1:dim(Veg)[2]){
  NAVA[n] <- sum(is.na(Veg[,i]))
  n = n + 1
}
NAVA
plot(x = Veg$BARESOIL, y = Veg$R, 
     xlab = "Exposed soil",
     ylab = "Species richness",
     main = "Scatter plot",
     xlim = c(0, 45), ylim = c(4, 19), pch = 18) # pch puede tomar 25 valores distintos. 1:25

Veg$Transect
unique(Veg$Transect)
plot(x = Veg$BARESOIL, y = Veg$R, xlab = "Exposed soil", ylab = "Species richness", 
     main = "Scatter plot", xlim = c(0, 45), ylim = c(4, 19), pch = Veg$Transect)

# Use of Vector for pch
unique(Veg$Time)
Veg$Time2 <- Veg$Time
unique(Veg$Time2)
Veg$Time2[Veg$Time <= 1974] <- 1
Veg$Time2[Veg$Time > 1974] <- 16
Veg$Time2

# 
plot(x = Veg$BARESOIL, y = Veg$R,
     xlab = "Exposed soil",
     ylab = "Species richness", main = "Scatter plot",
     xlim = c(0, 45), ylim = c(4, 19),
     pch = Veg$Time2)

# Loops y functions
setwd("C:/Users/Toshiba/OneDrive - Escuela Superior Politécnica del Litoral/Archivos adjuntos/R/Libros/1er libro leido/RBook/")
Owls <- read.table(file = "Owls.txt", header = TRUE)
names(Owls)
str(Owls)
