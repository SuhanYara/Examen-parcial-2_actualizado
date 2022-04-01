###########
# PARTE 1. 
###########

# La tabla alojada en el objeto resultados.RDS (anexo a esta entrada) tiene identificadores de  genes y  
# contrastes de medidas de expresión así como p-values para un experimento de microarreglos en donde se 
# quiere comparar una condición A vs WT. El contraste realizado es A-WT. Escribe un programa en R que:
  
# Lea la tabla
library(limma) #Abrimos la librería limma para leer los datos
library(Biobase)
load("01_Raw Data/resultados.RDS") #Cargamos así, porque no hay un título que conocemos
fit2

# Hacemos revisión de datos por valor p y logFC
topTable(fit2, sort.by = "p")
topTable(fit2, number=20, sort.by = "logFC")


# 1. Haga un boxplot de expresión ####
##############################################################

# Primero para ello debemos guardarlo en un nuevo objeto/tabla para que la función boxplot lo reconozca
Expresiones<-topTable(fit2, number = nrow(fit2), sort.by = "logFC")
Expresiones
boxplot(Expresiones$logFC, main = "Expresión de genes", xlab = "Genes", ylab = "logFC")

# Comentario ####
# Es notorio que las expresiones se dieron entre valor de 3 a -3, 6 veces arriba o abajo de la base. 


# 2. Grafique un PCA ####
##############################################################
fit2
Componentes<- princomp(fit2$design)
plot(Componentes$loadings, main="Componentes principales de la expresionen logFC")

# Comentario ####
# Aquí debido a los datos, considero que el análisis de componentes, en mi gráfica no es representativa,
# si vamos cambiando el contenido de fit2, el design es el único que separa KO de WT, los demás se van hacia un sólo grupo de datos
# o marca varios, para cada gen. Entonces en este creo que es lo más cercano a diferenciar los dos tratamientos usados.



# 3. Encuentre las sondas que se sobre expresan y sub expresan ####
##############################################################
# De la tabla que creamos 
Expresiones
# Podemos extraer los genes regulados, seleccionando concidicones de interés
# Primero vamos a agregar una columna a la tabla para especificar si Up o down

Expresiones$diffexpressed <- "NO" # Con esto decimos que no
# Y condicionamos los valores de expresion (logFC) y los de el valor P para definir Up o Down
Expresiones$diffexpressed[Expresiones$logFC > 1 & Expresiones$P.Value < 0.05] <- "UP"
Expresiones$diffexpressed[Expresiones$logFC < -1 & Expresiones$P.Value < 0.05] <- "DOWN"

# Ahora solicitamos las listas
which(Expresiones$diffexpressed == "DOWN")
which(Expresiones$diffexpressed == "UP")
# Aunque sólo da el número o posición en la lista de genes

# Guardamos en un vector la lista
which(Expresiones$diffexpressed == "DOWN")-> DownRegulated
Expresiones[DownRegulated,] # solicitamos los nombres, que podríamos guardar nuevamente en otra lista.

which(Expresiones$diffexpressed == "UP")-> UpRegulated
Expresiones[UpRegulated,]

# 4. Cuente cuántas sondas se sobre expresan y cuántas se sub expresan ####
##############################################################
# Contamos de los vectores creados
length(DownRegulated)
length(UpRegulated)

# Con la función  podemos ver genes Upregulated o Downrgulated, ya que identifica genes diferencialmente expresados según su valor p
Up_downR<- decideTests(fit2, p.value=0.05, lfc =1 )
summary(Up_downR) # Aquí podemos ver los resultados
# Notese que aquí, el número de genes es mayor, a lo que nosotros obtuvimos
# Se deba a las variable o indicador de cambio que insertemos incluso el método

# 5. Genera una figura de volcán manualmente, que incluya todas las sondas ####
##############################################################
# A partir de los datos que usamos arriba, usamos ggplot, para hacer un gráfico volcan

library(ggplot2)
# Creamos un objeto de la gráfica que haremos para después insertale loca cambios requeridos
Volcan<-ggplot(data=Expresiones, aes(x=logFC, y=-log10(P.Value))) + geom_point() + theme_minimal() 

mycolors <- c("green", "red", "black") # Insertamos los colores que deseamos para cada regulación
names(mycolors) <- c("DOWN", "UP", "NO") # Colocamos un vector con los nombres que etiquetaremos para cada tipo de regulación

# Añadimos otro objeto para que los valores del color se relacionen con los de la regulación 
Coloritos <- Volcan + scale_colour_manual(values = mycolors)

# Creamos una nueva columna "etiquetas, que contendrá el nombre de los genes expresados diferencialmente 
# (NA en caso de que no lo estén)

Expresiones$etiquetas <- NA
Expresiones$etiquetas[Expresiones$diffexpressed != "NO"] <- Expresiones$ID.Symbol[Expresiones$diffexpressed != "NO"]
# Aquí, tomamos los nombres de los genes que están bajos las condiciones que ya habíamos determinado para agregarlos en
# la nueva columna que creamos, y que después sólo esos aparezcan en la gráfica. 

ggplot(data=Expresiones, aes(x=logFC, y=-log10(P.Value), col=diffexpressed, label=etiquetas)) + geom_point() + theme_minimal() + geom_text()

# Podemos agregar líneas de límite de logFC y p-valor
ggplot(data=Expresiones, aes(x=logFC, y=-log10(P.Value), col=diffexpressed, label=etiquetas,)) + geom_point() + theme_minimal() + geom_text() + geom_vline(xintercept=c(-1, 1), col="red") + geom_hline(yintercept=-log10(0.05), col="red")


# 6. Realiza un análisis de GO para las tres ontologías para los genes diferencialmente expresados ####
##############################################################
# Guardamos las listas y vamos a hacer un análisis GO
Downlist<- Expresiones[DownRegulated,] 
Uplist<- Expresiones[UpRegulated,]

# Extraemos los nombres Y GUARDAMOS EL OBJETO, y despuñe slo enviamos a un txt a usar...
NamesDown<-Downlist$ID.Symbol
NamesUp<-Uplist$ID.Symbol

write.table(NamesDown, "03_Output/NamesDown.txt", fileEncoding = "UTF-8", quote = FALSE)
write.table(NamesUp, "03_Output/NamesUp.txt", fileEncoding = "UTF-8", quote = FALSE)

# Nos vamos a http://www.pantherdb.org/ y caragmos indepdniente las listas de Up y Down, y descargamos para Función Biológica (BF)
# MFunción molecular (MF), y componente celular (CC), los guardamos en RAw_Data

# Genera gráficas y/o tablas de las categorías sobre o sub expresadas ####
BFUp<-read.table(file ="01_Raw Data/BF_Up.txt", sep = "\t")
MFUp<-read.table(file ="01_Raw Data/MF_Up.txt", sep = "\t")
CCUp<-read.table(file ="01_Raw Data/CC_Up.txt", sep = "\t")
BFdw<-read.table(file ="01_Raw Data/BF_Dw.txt", sep = "\t")
MFdw<-read.table(file ="01_Raw Data/MF_Dw.txt", sep = "\t")
CCdw<-read.table(file ="01_Raw Data/CC_Dw.txt", sep = "\t")

#POdemos explorar cada una:
BFUp # y Extraer de cada una datos como el porcentaje de genes o número de genes cuya función ubicamos
MFUp
CCUp
BFdw
BFdw
BFdw
#3 POnemos nombre a las columnas

colnames(BFUp)<-c("","Función biológica",  "No.genes", "Porcentaje")
colnames(MFUp)<-c("","Función biológica",  "No.genes", "Porcentaje")
colnames(CCUp)<-c("","Función biológica",  "No.genes", "Porcentaje")
colnames(BFdw)<-c("","Función biológica",  "No.genes", "Porcentaje")
colnames(BFdw)<-c("","Función biológica",  "No.genes", "Porcentaje")
colnames(BFdw)<-c("","Función biológica",  "No.genes", "Porcentaje")

# Hacemos histogramas, aquí sólo dejaré uno. Como ejemplo, con las tablas ya podemos darnos una idea 
hist(BFUp$No.genes, main = "Histograma de genes sobre expresados y su función biológica", xlab = "No. genes", ylab = "Frecuencia")



# A partir de este resultado elabora una hipótesis de lo que trata el experimento####
# Muchos genes está asociados a respuestas a estimulos, desarrollo y metabolismo, creo que 
# podría ser un experimento de genes que controlan el desarrollo


