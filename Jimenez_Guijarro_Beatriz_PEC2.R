## ----setup, include=FALSE-------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, message = FALSE, warning = FALSE,
                      cache = TRUE, comment = NA)


## ----echo=TRUE, eval=FALSE------------------------------------------------------------------
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# install.packages("dplyr")
# install.packages("ggplot2")
# install.packages("ggrepel")
# install.packages("gplots")
# BiocManager::install("Biobase")
# BiocManager::install("oligo")
# BiocManager::install("arrayQualityMetrics")
# BiocManager::install("genefilter")
# BiocManager::install("pd.mouse430.2")
# BiocManager::install("mouse4302.db")
# BiocManager::install("limma")
# BiocManager::install("AnnotationDbi")
# BiocManager::install("org.Mm.eg.db")
# BiocManager::install("clusterProfiler")
# BiocManager::install("enrichplot")


## -------------------------------------------------------------------------------------------
# Eliminar las muestras tomadas a las 2 horas de "allTargets.txt"

# Cargar el archivo "allTargets.txt"
allTargets <- read.table("allTargets.txt", header = TRUE, sep = " ", stringsAsFactors = FALSE)

# Filtrar las muestras eliminando aquellas que tienen "hour 2" en la columna 'time'
filteredTargets <- allTargets[allTargets$time != "hour 2", ]

# Mostrar la cantidad de muestras eliminadas y las restantes
print("Muestras después de eliminar las que están a las 2 horas:")
print(paste("Muestras originales:", nrow(allTargets)))
print(paste("Muestras restantes:", nrow(filteredTargets)))

# Mostrar la nueva tabla filtrada "filteredTargets"
filteredTargets


## -------------------------------------------------------------------------------------------
# Añadir los nombres de los archivos .CEL a "filteredTargets"

# Obtener la lista de archivos .CEL en la carpeta "GSE38531"
files_cel <-list.files("GSE38531") 

# Extraer los primeros 9 caracteres de los nombres de archivos .CEL
sampleNames <- substr(files_cel, 1, 9) 

# Filtrar solo los archivos cuyos nombres coinciden con los valores en "sample"
matching_files <- files_cel[sampleNames %in% filteredTargets$sample]
matching_sampleNames <- sampleNames[sampleNames %in% filteredTargets$sample]

# Verificar coincidencia entre los nombres de archivos .CEL y las filas de "filteredTargets"
# y unir los nombres de los archivos .CEL a "filteredTargets" creando un nuevo dataframe "targets"
if (length(matching_files) > 0) {
  
  # Reordenar los nombres y archivos .CEL para que coincidan con el orden de "filteredTargets"
  order <- match(filteredTargets$sample, matching_sampleNames)
  # Unir los nombres completos de los archivos al dataframe "targets" en el orden original
  targets <- cbind(
    filenames = matching_files[order],
    filteredTargets)
  
} else {
  # Mensaje de error si no coinciden
  warning("Los nombres de los archivos no coinciden con las muestras del data frame.")
}

# Mostrar el resultado final de "targets"
head(targets)


## -------------------------------------------------------------------------------------------
# Aplicar "selectSamples" para obtener 24 muestras de "targets"

# Semilla con mi DNI (sin letra)
seed <- 49218850

# Llamada al archivo "selectSamples.R" y a su función interna "filter_microarray"
source("selectSamples.R")
finalTargets <- filter_microarray(targets, seed)

# Ordenar las muestras de "finalTargets" en el orden original de "allTargets.txt"
finalTargets <- finalTargets[order(as.numeric(sub(".*\\.", "", rownames(finalTargets)))), ]

# Mostrar las dimensiones finales de "finalTargets" y el numero de muestras por
# combinación "infection" (infección) - "agent" (tratamiento)
print(paste("Muestras finales:", nrow(finalTargets)))
print("Muestras para cada combinación infection (infección) - agent (tratamiento):")
table(finalTargets$infection, finalTargets$agent)


## -------------------------------------------------------------------------------------------
library(dplyr)

# Crear las columnas "shortName" y "group" en el archivo "finalTargets"
finalTargets <- finalTargets %>%
  mutate(
    # Crear "shortName" extrayendo el nombre entre el último guion bajo y ".CEL"
    shortName = sub(".*_(.*?)_Mouse430\\+2\\.CEL", "\\1", filenames),
    
    # Crear "group", de tipo factor, según las combinaciones de "infection" y "agent"
    group = factor(case_when(
      infection == "uninfected" & agent == "untreated" ~ "Ctl",
      infection == "uninfected" & agent == "linezolid" ~ "Ctl-Lin",
      infection == "uninfected" & agent == "vancomycin" ~ "Ctl-Van",
      infection == "S. aureus USA300" & agent == "untreated" ~ "24h",
      infection == "S. aureus USA300" & agent == "linezolid" ~ "24h-Lin",
      infection == "S. aureus USA300" & agent == "vancomycin" ~ "24h-Van",
      TRUE ~ NA_character_ # Valores predeterminados como NA
    ), levels = c("Ctl", "Ctl-Lin", "Ctl-Van", "24h", "24h-Lin", "24h-Van"))
  )

# Mostrar el archivo "finalTargets" modificado y con las muestras seleccionadas
finalTargets



## ----echo=TRUE------------------------------------------------------------------------------
library(Biobase)

# Crear objeto de tipo "AnnotatedDataFrame" de "finalTargets"
dataTargets <- AnnotatedDataFrame(finalTargets)


## -------------------------------------------------------------------------------------------
library(oligo)

# Obtener las rutas de los archivos .CEL correspondientes a las muestras seleccionadas
celFiles <- paste0("GSE38531/", pData(dataTargets)$filenames)

# Crear el objeto ExpressionSet con los datos en bruto de los archivos .CEL y los datos
# de "dataTargets"
rawData <- read.celfiles(celFiles, phenoData = dataTargets)



## -------------------------------------------------------------------------------------------
# Modificar nombres de las muestras por su "shortName" del objeto "dataTargets"
rownames(pData(rawData)) <- dataTargets$shortName
colnames(rawData) <-rownames(pData(rawData))

# Mostrar la información sobre los datos en bruto "rawData"
show(rawData)


## -------------------------------------------------------------------------------------------
# Asignar colores a los grupos
colors <- c(rep("red",4), rep("yellow",4), rep("blue",4),
            rep("green",4), rep("orange",4), rep("purple",4))


## -------------------------------------------------------------------------------------------
# Boxplots de los datos en bruto de "rawData"
boxplot(rawData, col = colors, names = dataTargets$shortName,
        main="Distribución de intensidad de los datos en bruto",
        las=2, cex.axis=0.8)


## -------------------------------------------------------------------------------------------
# Diagrama de densidad de la señal de las muestras en "rawData"
hist(rawData, col = colors, lty =1:nrow(pData(rawData)),
     main ="Densidad de la señal de las muestras en bruto")

# Leyenda para el diagrama de densidad
legend(x ="topright",legend = dataTargets$shortName, col = colors,
       lty =1:nrow(pData(rawData)), cex =0.5)


## -------------------------------------------------------------------------------------------
# Análisis de las componentes principales de "rawData"
pc <- prcomp(t(exprs(rawData)), scale. = FALSE)
summary(pc)


## -------------------------------------------------------------------------------------------
library(ggplot2)
library(ggrepel)

# Crear función para realizar el gráfico de las componentes principales (PC)
plotPCA <- function (datos, labels, factor, title, scale,colores, size = 1.5, glineas = 0.25) {
  data <- prcomp(t(datos),scale=scale)
  #Ajustes del gráfico
  dataDf <- data.frame(data$x)
  Group <- factor
  loads <- round(data$sdev^2/sum(data$sdev^2)*100,1)
  # Gráfico principal
  p1 <- ggplot(dataDf,aes(x=PC1, y=PC2)) +
    theme_classic() +
    geom_hline(yintercept = 0, color = "gray70") +
    geom_vline(xintercept = 0, color = "gray70") +
    geom_point(aes(color = Group), alpha = 0.55, size = 3) +
    coord_cartesian(xlim = c(min(data$x[,1])-5,max(data$x[,1])+5)) +
    scale_fill_discrete(name = "Group")
  # Evitar superposición de etiquetas
  p1 + geom_text_repel(aes(y = PC2 + 0.25, label = labels),segment.size = 0.25, size = size) +
    labs(x = c(paste("PC1",loads[1],"%")),y=c(paste("PC2",loads[2],"%"))) +
    ggtitle(paste("Principal Component Analysis for: ",title,sep=" ")) +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_color_manual(values=colores)
}


## -------------------------------------------------------------------------------------------
# Crear el gráfico de componentes principales de "rawData" a partir de la función
plotPCA(exprs(rawData), labels = dataTargets$shortName, factor = dataTargets$group,
        title="rawData", scale = FALSE, size = 3,
        colores = c("red", "yellow", "blue", "green", "orange", "purple"))


## -------------------------------------------------------------------------------------------
# Cálculo de las distancias y composición del cluster
dist_clust <- hclust(dist(t(exprs(rawData))), method = "average")

# Cluster jerárquico de los datos en bruto "rawData"
plot(dist_clust, labels = dataTargets$shortName, hang = -1, xlab = "Distancia euclídea",
     main = "Clúster jerárquico de los datos en bruto",)


## -------------------------------------------------------------------------------------------
# Control de calidad de los datos en bruto
library(arrayQualityMetrics)
arrayQualityMetrics(rawData, outdir = "arrayQualityMetrics_PEC2",
                    force = TRUE, intgroup = "infection")


## -------------------------------------------------------------------------------------------
# Heatmap basado en diatancias de los datos en bruto
man_dist <- dist(t(exprs(rawData)))
heatmap(as.matrix(man_dist), col = heat.colors(16))


## -------------------------------------------------------------------------------------------
# Normalización de los datos en bruto "rawData" mediante el paquete "oligo"
library(oligo)
esetData_rma <- oligo::rma(rawData)
esetData_rma


## -------------------------------------------------------------------------------------------
# Boxplots de los datos normalizados "esetData_rma"
boxplot(esetData_rma, col = colors, names = dataTargets$shortName,
        main="Distribución de intensidad de los datos normalizados",
        las=2, cex.axis=0.8)


## -------------------------------------------------------------------------------------------
# Diagrama de densidad de la señal de las muestras normalizadas en "esetData_rma"
hist(esetData_rma, col = colors, lty =1:nrow(pData(esetData_rma)),
     main ="Densidad de la señal de las muestras normalizadas")

# Leyenda para el diagrama de densidad
legend(x ="topright",legend = dataTargets$shortName, col = colors,
       lty =1:nrow(pData(esetData_rma)), cex =0.5)


## -------------------------------------------------------------------------------------------
# Análisis de las componentes principales de los datos normalizados "esetData_rma"
pc <- prcomp(t(exprs(esetData_rma)), scale. = FALSE)
summary(pc)


## -------------------------------------------------------------------------------------------
# Crear el gráfico de componentes principales de "esetData_rma" a partir de la función
plotPCA(exprs(esetData_rma), labels = dataTargets$shortName, factor = dataTargets$group,
        title="esetData_rma", scale = FALSE, size = 3,
        colores = c("red", "yellow", "blue", "green", "orange", "purple"))


## -------------------------------------------------------------------------------------------
# Cálculo de las distancias y composición del cluster
dist_clust1 <- hclust(dist(t(exprs(esetData_rma))), method = "average")

# Cluster jerárquico de los datos en bruto "rawData"
plot(dist_clust1, labels = dataTargets$shortName, hang = -1, xlab = "Distancia euclídea",
     main = "Clúster jerárquico de los datos normalizados",)


## -------------------------------------------------------------------------------------------
# Control de calidad de los datos normalizados
library(arrayQualityMetrics)
arrayQualityMetrics(esetData_rma, outdir = "arrayQualityMetrics_Norm_PEC2",
                    force = TRUE, intgroup = "infection")


## -------------------------------------------------------------------------------------------
library(genefilter)
library(mouse4302.db)

# Modificación de las anotaciones de "esetData_rma" por la base de datos "mouse4302.db"
annotation(esetData_rma) <- "mouse4302.db"

# Filtrado de los datos normalizados "esetData_rma"

filtered <- nsFilter(esetData_rma,
                     require.entrez = TRUE, # Requerir que las sondas tengan ID Entrez
                     remove.dupEntrez = TRUE, # Eliminar duplicados basados en ID Entrez
                     var.filter = TRUE,  # Activa el filtrado por variabilidad
                     var.func = IQR, # Función estadística de filtrado por función
                     var.cutoff = 0.90,  # Elige el 10% superior (umbral: 90% menos variable)
                     # var.cutoff es un cuantil de todos los valores de var.func
                     filterByQuantile = TRUE, 
                     feature.exclude = "^AFFX" # Excluir sondas de control Affymetrix
)


## -------------------------------------------------------------------------------------------
# Objetos devueltos por la función "nsFilter"
names(filtered)


## -------------------------------------------------------------------------------------------
# Clase del objeto "eset" devuelto por la función "nsFilter"
class(filtered$eset)

# Información del objeto "eset"
filtered$eset


## -------------------------------------------------------------------------------------------
# Contenido del informe de filtrado devuelto por la función "nsFilter"
print(filtered$filter.log)


## ----echo=TRUE------------------------------------------------------------------------------
# Almacenar genes filtrados en "esetData_filtered"
esetData_filtered <-filtered$eset


## -------------------------------------------------------------------------------------------
# Guardar los datos de "esetData_rma" y "esetData_filtered"
save(esetData_rma, esetData_filtered, file="./results/normalized.filtered.Data.Rda")

# Guardar los valores de expresión de "esetData_rma" en un archivo .csv
write.csv(exprs(esetData_rma), file="./results/normalized.Data.csv")

# Guardar los valores de expresión de "esetData_filtered" en un archivo .csv
write.csv(exprs(esetData_filtered), file="./results/filtered.Data.csv")


## -------------------------------------------------------------------------------------------
# Matriz de diseño de "esetData_filtered"
design_matrix <- model.matrix(~0+group, pData(esetData_filtered))
colnames(design_matrix) <- c("Ctl", "Ctl.Lin", "Ctl.Van", "inf24h", "inf24h.Lin", "inf24h.Van")
print(design_matrix)


## -------------------------------------------------------------------------------------------
library(limma)

# Matriz de contrastes
cont_matrix <- makeContrasts (inf24h_vs_Ctl = inf24h - Ctl,
                              inf24h.Lin_vs_Ctl.Lin = inf24h.Lin - Ctl.Lin,
                              inf24h.Van_vs_Ctl.Van = inf24h.Van - Ctl.Van,
                              levels = design_matrix)

print(cont_matrix)


## -------------------------------------------------------------------------------------------
library(limma)

# Ajustar el modelo lineal
fit <- lmFit(esetData_filtered, design_matrix)

# Aplicar los contrastes
fit_contrasts <- contrasts.fit(fit, cont_matrix)

# Realizar la estimación bayesiana para calcular valores p y estadísticos
fit_contrasts <- eBayes(fit_contrasts)

# Comprobar la clase del objeto final
class(fit_contrasts)


## -------------------------------------------------------------------------------------------
# Tabla para el contraste Infectados vs no infectados sin tratamiento
topTab_inf24hvsCtl <- topTable (fit_contrasts, number=nrow(fit_contrasts),
                                coef="inf24h_vs_Ctl", adjust="fdr")

# Cabecera de la tabla
head(topTab_inf24hvsCtl)


## -------------------------------------------------------------------------------------------
# Tabla para el contraste Infectados vs no infectados tratados con LINEZOLID
topTab_inf24h.LinvsCtl.Lin <- topTable (fit_contrasts, number=nrow(fit_contrasts),
                                coef="inf24h.Lin_vs_Ctl.Lin", adjust="fdr")

# Cabecera de la tabla
head(topTab_inf24h.LinvsCtl.Lin)


## -------------------------------------------------------------------------------------------
# Tabla para el contraste Infectados vs no infectados tratados con VANCOMICINA
topTab_inf24h.VanvsCtl.Van <- topTable (fit_contrasts, number=nrow(fit_contrasts),
                                coef="inf24h.Van_vs_Ctl.Van", adjust="fdr")

# Cabecera de la tabla
head(topTab_inf24h.VanvsCtl.Van)


## -------------------------------------------------------------------------------------------
library(mouse4302.db)
library(AnnotationDbi)

# Extraer los símbolos genéticos asociados a las sondas en el modelo "fit_contrasts"
geneSymbols <- AnnotationDbi::select(mouse4302.db, rownames(fit_contrasts), c("SYMBOL"))

# Obtener un vector de los símbolos genéticos para etiquetar los puntos en el volcanoplot
SYMBOLS <- geneSymbols$SYMBOL

# volcanoplot para la comparación Infectados vs no infectados sin tratamiento
volcanoplot(fit_contrasts, coef=1, highlight=6, names=SYMBOLS, 
            main=paste("Genes diferencialmente expresados",
            colnames(cont_matrix)[1], sep="\n"))

# Añadir líneas verticales en el gráfico en los valores de logFC -1 y 1
abline(v=c(-1,1))


## -------------------------------------------------------------------------------------------
library(mouse4302.db)

# volcanoplot para la comparación Infectados vs no infectados tratados con LINEZOLID
volcanoplot(fit_contrasts, coef=1, highlight=6, names=SYMBOLS, 
            main=paste("Genes diferencialmente expresados",
            colnames(cont_matrix)[2], sep="\n"))

# Añadir líneas verticales en el gráfico en los valores de logFC -1 y 1
abline(v=c(-1,1))


## -------------------------------------------------------------------------------------------
# volcanoplot para la comparación Infectados vs no infectados tratados con VANCOMICINA
volcanoplot(fit_contrasts, coef=1, highlight=6, names=SYMBOLS, 
            main=paste("Genes diferencialmente expresados",
            colnames(cont_matrix)[3], sep="\n"))

# Añadir líneas verticales en el gráfico en los valores de logFC -1 y 1
abline(v=c(-1,1))


## -------------------------------------------------------------------------------------------
library(limma)
# Comparaciones múltiples del modelo "fit_contrasts" con "decideTests"
res <- decideTests(fit_contrasts, method = "separate",
                   adjust.method = "fdr", p.value = 0.1, lfc = 1)


## -------------------------------------------------------------------------------------------
# Resumen del análisis realizado con "decideTests"
sum_res_rows <- apply(abs(res),1,sum)
res_selected <- res[sum_res_rows!=0,]
print(summary(res_selected))


## -------------------------------------------------------------------------------------------
# Diagrama de Venn para las comparaciones múltiples
vennDiagram (res_selected[,1:3], cex=0.9, main ="Genes en común entre las tres comparaciones")


## -------------------------------------------------------------------------------------------
library(AnnotationDbi)
# Extraer los datos de expresión de las sondas de "esetData_filtered"
# seleccionando únicamente los genes que están presentes en "res_selected"
heatmap_data <- exprs(esetData_filtered)[rownames(exprs(esetData_filtered)) %in% rownames(res_selected),]

# Extraer los símbolos genéticos asociados a los genes seleccionados para el heatmap
geneSymbols1 <- AnnotationDbi::select(mouse4302.db, rownames(heatmap_data), c("SYMBOL"))

# Obtener un vector de los símbolos genéticos para etiquetar los puntos en el heatmap
SYMBOLS1 <- geneSymbols1$SYMBOL

# Renombrar las filas del heatmap para usar los símbolos genéticos en lugar de los IDs de genes
rownames(heatmap_data) <- SYMBOLS1


## -------------------------------------------------------------------------------------------
library(gplots)

# Heatmap de los datos seleccionados como diferencialmente expresados "res_selected"
heatmap.2(heatmap_data, Rowv = TRUE, Colv = TRUE, dendrogram = "both",
          main = "Genes diferencialmente expresados \n p-valor < 0,1, logFC >=1",
          scale = "row", col = rainbow(100), sepcolor = "white", sepwidth = c(0.05,0.05),
          cexRow = 0.5, cexCol = 0.9, key = TRUE, keysize = 1.5, density.info = "histogram",
          ColSideColors = colors, tracecol = NULL, srtCol = 30)


## -------------------------------------------------------------------------------------------
library(mouse4302.db)

# Tipos de anotaciones de la base de datos "mouse4302.db"
keytypes(mouse4302.db)


## -------------------------------------------------------------------------------------------
library(AnnotationDbi)

# Obtenemos las anotaciones de la "topTable" de Infectados vs no infectados sin tratamiento
geneAnots_inf24hvsCtl <- AnnotationDbi::select(mouse4302.db, rownames(topTab_inf24hvsCtl),
                                   keytype = "PROBEID",
                                   c("SYMBOL", "ENTREZID", "ENSEMBL", "GENENAME"))

head(geneAnots_inf24hvsCtl)

# Guardar los resultados completos en un archivo .csv
write.csv(geneAnots_inf24hvsCtl, file="./results/topAnnotated_inf24hvsCtl.csv")


## -------------------------------------------------------------------------------------------
library(AnnotationDbi)

# Obtenemos las anotaciones de la "topTable" de Infectados vs no infectados tratados con LINEZOLID
geneAnots_inf24h.LinvsCtl.Lin <- AnnotationDbi::select(mouse4302.db,
                                            rownames(topTab_inf24h.LinvsCtl.Lin),
                                            keytype = "PROBEID",
                                            c("SYMBOL", "ENTREZID", "ENSEMBL", "GENENAME"))

head(geneAnots_inf24h.LinvsCtl.Lin)

# Guardar los resultados completos en un archivo .csv
write.csv(geneAnots_inf24h.LinvsCtl.Lin, file="./results/topAnnotated_inf24h.LinvsCtl.Lin.csv")


## -------------------------------------------------------------------------------------------
library(AnnotationDbi)

# Obtenemos las anotaciones de la "topTable" de Infectados vs no infectados tratados con LINEZOLID
geneAnots_inf24h.VanvsCtl.Van <- AnnotationDbi::select(mouse4302.db,
                                            rownames(topTab_inf24h.VanvsCtl.Van),
                                            keytype = "PROBEID",
                                            c("SYMBOL", "ENTREZID", "ENSEMBL", "GENENAME"))

head(geneAnots_inf24h.VanvsCtl.Van)

# Guardar los resultados completos en un archivo .csv
write.csv(geneAnots_inf24h.VanvsCtl.Van, file="./results/topAnnotated_inf24h.VanvsCtl.Van.csv")


## -------------------------------------------------------------------------------------------
library(mouse4302.db)
library(AnnotationDbi)

# Crear una lista con las listas de genes diferencialmente expresados ("topTable")
# para cada comparación
listOfTables <- list(inf24hvsCtl = topTab_inf24hvsCtl,
                     inf24h.LinvsCtl.Lin  = topTab_inf24h.LinvsCtl.Lin,
                     inf24h.VanvsCtl.Van = topTab_inf24h.VanvsCtl.Van)

# Inicializar una lista vacía para almacenar los genes seleccionados en cada comparación
listOfSelected <- list()

# Iterar a través de las "topTable"
for (i in 1:length(listOfTables)){
  
  # Seleccionar la "topTable" correspondiente a la comparación actual
  topTab <- listOfTables[[i]]
  
  # Obtener los nombres de los genes presentes en la "topTable"
  probesUniverse <- rownames(topTab)
  
  # Mapear los genes a sus correspondientes identificadores ENTREZ utilizando "mouse4302.db"
  entrezUniverse<- AnnotationDbi::select(mouse4302.db, probesUniverse, "ENTREZID")
  
  # Extraer los identificadores ENTREZ del resultado
  entrezUniverse <- entrezUniverse$ENTREZID
  
  # Eliminar identificadores ENTREZ duplicados
  entrezUniverse <- entrezUniverse[!duplicated(entrezUniverse)]
  
  # Identificar genes significativamente expresados con p-valor < 0.1 y logFC > 1
  whichGenes<- topTab["adj.P.Val"]<0.1 & topTab["logFC"] > 1
  
  # Filtrar los identificadores ENTREZ de los genes significativos
  topGenes <- entrezUniverse[whichGenes]
  
  # Eliminar duplicados en la lista de genes seleccionados
  topGenes <- topGenes[!duplicated(topGenes)]
  
  # Almacenar los genes seleccionados en la lista de resultados
  listOfSelected[[i]] <- topGenes
  
  # Nombrar cada elemento de la lista con el nombre de la comparación
  names(listOfSelected)[i] <- names(listOfTables)[i]
}

# Imprimir la cantidad de genes seleccionados para cada comparación
sapply(listOfSelected, length)


## -------------------------------------------------------------------------------------------
library(org.Mm.eg.db)
library(clusterProfiler)
library(enrichplot)

# Iterar sobre la lista de genes seleccionados para cada comparación
for (i in 1:length(listOfSelected)){
  
  # Realizar análisis de enriquecimiento GO (Gene Ontology) para cada comparación
  enrich_go <- enrichGO(
    gene = as.integer(listOfSelected[[i]]), #Lista de genes seleccionados (conversión a enteros) 
    universe = entrezUniverse, # Genes "universales" (todos los de la "topTable")
    keyType = "ENTREZID", # Identificador de genes ENTREZ
    OrgDb = org.Mm.eg.db, # Base de datos de anotación para el organismo (Mus musculus)
    ont = "BP", # Ontología de interés (Biological Process)
    pAdjustMethod = "BH", # Método de ajuste de p-valor 
    qvalueCutoff = 0.25, # Valor q máximo permitido
    readable = TRUE) # Convertir IDs de genes a nombres legibles
  
  # Separador para cada comparación
  cat("##################################")
  cat("\nComparison: ", names(listOfSelected)[i],"\n")
  
  # Extraer los resultados del análisis de enriquecimiento en otro data frame duplicado
  enrich_results <- as.data.frame(enrich_go)

  # Eliminar la columna "geneID" del data frame de resultados duplicado
  enrich_results$geneID <- NULL
  
  # Mostrar las primeras filas del resultado sin la columna "geneID"
  print(head(enrich_results))
  
  # Si hay resultados significativos
  if (length(rownames(enrich_go@result)) != 0) {
    
    # Guardar los resultados de cada comparación en formato .csv en la carpeta "results"
    write.csv(as.data.frame(enrich_go), 
              file =paste0("./results/","enrichGO.Results_", names(listOfSelected)[i], ".csv"), 
              row.names = FALSE)
    
    # Crear un archivo .pdf para cada comparación para almacenar gráficos de enriquecimiento
    pdf(file = paste0("./results/","EnrichmentGO_all_plots_", names(listOfSelected)[i], ".pdf"))
      
      # Generar y guardar un gráfico de barras de los términos enriquecidos
      print(barplot(enrich_go, showCategory = 15, font.size = 4,
          title = paste0("Gene Ontology Enrichment Analysis (ORA) for ",
                         names(listOfSelected)[i],". Barplot")))
      
      # Generar y guardar un gráfico de red de genes relacionados con términos enriquecidos
      print(cnetplot(enrich_go, categorySize = "geneNum", schowCategory = 15,
                     vertex.label.cex = 0.75,
                     title = paste0("Gene Ontology Enrichment Analysis (ORA) for ",
                                    names(listOfSelected)[i],". Gene Network")))
      
      # Generar y guardar un gráfico de puntos para la visualización de términos enriquecidos
      print(dotplot(enrich_go, showCategory = 15,
                    title = paste0("Gene Ontology Enrichment Analysis (ORA) for ",
                                   names(listOfSelected)[i],". Dotplot")))
      
      # Generar y guardar un gráfico jerárquico de términos GO enriquecidos
      print(goplot(enrich_go, showCategory = 15,
                   title = paste0("Gene Ontology Enrichment Analysis (ORA) for ",
                                   names(listOfSelected)[i],". Go Hierarchy")))
      
      # Generar y guardar un mapa de enriquecimiento
      print(emapplot(pairwise_termsim(enrich_go), cex_label_category=0.5,
                     title = paste0("Gene Ontology Enrichment Analysis (ORA) for ",
                                   names(listOfSelected)[i],". Enrichment Map")))

    # Cerrar el archivo pdf
    dev.off()
  }
  
  # Mostrar un gráfico de red de genes relacionados con términos enriquecidos para cada comparación
  print(cnetplot(enrich_go, categorySize = "geneNum", schowCategory = 15,
                 vertex.label.cex = 0.75,
                 title = paste0("Gene Ontology Enrichment Analysis (ORA) for ",
                                    names(listOfSelected)[i],". Gene Network")))
}



## ----eval=FALSE-----------------------------------------------------------------------------
# knitr::purl("Jimenez_Guijarro_Beatriz_PEC2.Rmd", output = "Jimenez_Guijarro_Beatriz_PEC2.R")


## ----codeScript, file="Jimenez_Guijarro_Beatriz_PEC2.R", echo=TRUE, eval=FALSE--------------
# 

