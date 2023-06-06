#Código usado para el proyecto

# Descargamos las bibliotecas necesarias para el funcionamiento del programa.
if(!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("msa")
install.packages ("seqinr")
install.packages ("ape")
install.packages("phangorn", dep=TRUE)

#Activamos las bibliotecas
library(phangorn)
library(msa)
library(seqinr)
library(ape)


#Extraemos la información del archivo con múltiples secuencias en formato FASTA
secuencia <- readAAStringSet("D:/Proyecto Final Genomica/secuencias.dna")
secuencia

#Alineamos las secuencias con ClustalW,
alineClustalW <- msa(secuencia, "ClustalW")
alineClustalW

#imprimimos el alineamiento
print(alineClustalW, show="complete")

#Generamos una matriz de consenso 
conMat <- consensusMatrix(alineClustalW)
md <- dim(conMat)
md

#Observamos la puntuación de conservación de alineamiento múltiple
data(BLOSUM62)
blos <- msaConservationScore(alineClustalW, BLOSUM62)
blos

#Convertimos a un objeto utilizable por seqinr
alineamientoCW <- msaConvert(alineClustalW, type="seqinr::alignment")

#Generamos la matriz de distancias
d <- dist.alignment(alineamientoCW, "identity")

#Creamos  distintos árboles con la matriz de distancias
clusW_NJ <- nj(d)
clusW_UPGMA <- upgma(d)
clusW_WPGMA <- wpgma(d)

#Graficamos los resultados
plot(clusW_UPGMA, main="Árbol UPGMA")
plot(clusW_UPGMA, main="Árbol WPGMA")
par(mar = c(2,2,2,2))
plot(clusW_NJ, main="Árbol con NJ")
plot(clusW_NJ)

#referencias
toBibtex(citation("msa"))



