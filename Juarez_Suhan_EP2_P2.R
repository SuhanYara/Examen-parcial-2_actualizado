###########
# PARTE 2. 
###########

# Elabora un programa completo de la práctica. Siempre puedes usar las versiones más ligeras de los archivos que se utilizan en la práctica. 
# En particular resuelve las preguntas contenidas en la práctica. Discute tus resultados.
# http://datos.langebio.cinvestav.mx/~cei/cursos/R/sesion08adv2.html

# Problema: El propósito de estos ejercicios es obtener la anotación del genoma con el que estamos trabajando, y utilizar esta anotación para 
# catalogar nuestro mapeo de la secuenciación masiva. El resultado de esta práctica será obtener una figura mostrando la frecuencia con la que 
# hemos secuenciado fragmentos de distintos tipos de elementos genómicos (regiones codificantes, miRNAs, rRNA, snoRNA, etc).

# Cargamos lo paquetes

library(IRanges)
library(rtracklayer)
library(GenomicRanges)
library(Rsamtools)
library(GenomeInfoDb)
# Primero generamos un objeto base Iranges
x <- IRanges(start=c(11,35,40), end=c(20,50,63)) # Colocamos el inicio y final
x

# Existen diferentes funciones para extraer información de ese objeto 
start(x) # Sólo marca un vector con los valores de inicio
end(x)   # o sólo los finales
width(x) # Ancho de cada rango
range(x) # Rango total de todos
coverage(x) # La cobertura, su suma en cada posición
reduce(x) # Une los rangos que están encimados

# Ahora, para ver la relación de esto 
exons <- reduce(x) # Creamos un objeto que se llame exones
reads <- IRanges(start=c(1,21,30,50,80), width=20) # Otro lecturas 
reads

# Verificamos las coincidencias de esos valores de exones con los de la base de datos en lecturas
countOverlaps(exons, reads)

# Ahora si cargamos datos experimentales ####
load(file = "01_Raw Data/human.Rdata")


# Si queremos revisar los datos podemos usar diferentes funciones como:
seqnames(human) # Nombres de las secuncias
ranges(human) # Rangos
strand(human) # Tipo de cadena
mcols(human)  # Metadatos globales
# De este último acceder de distintas formas
table(mcols(human)$gene_biotype)

# Y quedarnos con algunas...
mcols(human) <- mcols(human)[,c("source","gene_id","gene_name","gene_biotype")]

mcols(human)
human

# Ejercicios ####
#1. ¿Cómo le harían para quedarse exclusivamente con las anotaciones de "miRNA"?
Cadenas<-strand(human)

#2. ¿ y solamente aquellas anotaciones de la cadena "-"?
cadenasNeg<-which(Cadenas == "-")
cadenasNeg # Aquí me quedaría con el número de gen con cadenas negativas 


#2. Anotación de secuencias mapeadas ####
# Vamos a utilizar el paquete Rsamtools cuando tengamos que importar datos de archivos SAM/BAM. 
library(Rsamtools)
# Para no importar absolutamente toda la información del archivo BAM, vamos a seleccionar cuáles datos queremos. 
what <- c("rname", "strand", "pos", "qwidth")
param <- ScanBamParam(what=what)
bam <- scanBam(file ="01_Raw Data/human_mapped_small.bam", param=param)

# Finalmente vamos a querer guardarlos en un objeto GRanges, por lo que necesitamos el nombre del cromosoma 
# (o secuencia de referencia), la cadena, la posición y la longitud de la secuencia que fue mapeada.
class(bam)
lapply(bam, names)

# Usamos la información anterior para construir un archivo Granges

mapGR <- GRanges(seqnames <- bam[[1]]$rname, ranges <- IRanges(start=bam[[1]]$pos, width=bam[[1]]$qwidth), strand   <- bam[[1]]$strand)
mapGR

# Contanis kas secuencias alineadoas con las de las anotaciones que yta teniamos antes
mcols(human)$counts<-countOverlaps(human, mapGR)
mcols(human)
# Sale un Warning, se refiere a que no coinciden absolutamente todos los nombres de los cromosomas 


# Agregamos valores de cada gen
# Por tipo de gen
typeCounts <- aggregate(mcols(human)$counts, by=list("biotype"=mcols(human)$gene_biotype), sum)
typeCounts

# Por gen individual
geneCounts = aggregate(mcols(human)$counts, by=list("id"=mcols(human)$gene_name), sum)
head(geneCounts)

# Las cuentas por gen las podríamos usar directamente en un proceso de análisis de expresión diferencial, por ejemplo:
minCount <- 40000
typeCountsHigh <- typeCounts[typeCounts$x > minCount,]
typeCountsHigh <- typeCountsHigh[order(typeCountsHigh$x),]
typeCountsHigh <- rbind(data.frame("biotype"="other",
                                  "x"=sum(typeCounts$x[typeCounts$x <= minCount])),
                       typeCountsHigh)

## Graficamos en uno de pastel
pie(typeCountsHigh$x, labels=typeCountsHigh$biotype, col=rev(rainbow(nrow(typeCountsHigh))),
    main="Number of aligned reads per biotype")
