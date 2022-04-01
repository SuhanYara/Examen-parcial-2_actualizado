###########
# PARTE 3. 
###########

# Genera un programa en R que modele lo siguiente: Se tienen tres genes, A,B, y C.
# El gen B está regulado negativamente por A, el A, negativamente por C y el C, negativamente por B.
library(BoolNet)
Genes<-loadNetwork("01_Raw Data/Genes")
Genes

# 1. ¿Cuántos atractores tiene esta red? ####
##################################################
plotNetworkWiring(Genes) # Dibujamos la gráfica

atractores<-getAttractors(Genes)
atractores

# Comentario ####
# Tiene dos atractores: uno de 2 estados y el segundo de 6 estados



# 2. ¿Cual es el estado más probable? ####
##################################################
# Comentario ####
# El segundo atractor, es el que tiene 75% de probabilidad, donde inicioa 100 ABC rescptivamente



#3. ¿Hay atractores cíclicos?
##################################################
plotStateGraph(atractores) #3 Con este podemo ver gráficamente
# Comentario ####
# Si, el segundo atractor es ciclico, como se nota en la fígura, contiene 6 estados.


# 4. Dibuja los atractores ####
plotAttractors(atractores) 
# Comentario ####
# Con al función anterior dibujamos los atractores, en cada estado. 

# 5. Dibuja todos los atractores juntos ####
# Comentario ####
#### Se usó la función plotplotStateGraph, en la parte 3.

# No hubiera esperado tantos estados para 3 genes en esas condiciones, para relaciones simple
# Pero don esto se denota la importancia de la generación de redes

