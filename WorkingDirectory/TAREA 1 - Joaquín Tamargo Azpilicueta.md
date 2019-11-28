TAREA 1. MICROARRAYS.
===============================================

* Alumno/a: Joaquín Tamargo Azpilicueta
* E-mail: joatamazp@alum.us.es
* Biología Molecular de Sistemas, Grado en Bioquímica (Universidad de Sevilla)
------

## Tabla de Contenidos
1. [Introducción](#introduccion)
2. [Análisis de calidad](#calidad)
3. [Preprocesamiento de los datos](#preprocesamiento)
4. [Análisis de expresión génica diferencial](#expresados)
5. [Enriquecimiento de términos de ontología de genes](#discusion)
6. [Conclusiones](#conclusion)
7. [Bibliografía](#bibliografia)

## Información complementaria
1. [README](./readme.r)
2. [Script completo](./script.r)

--------

# Análisis de expresión a genoma completo de los efectos de la sacarosa en  plántulas de *Arabidopsis thaliana* sometidas a condiciones de anoxia

## 1. Introducción 
<a class="anchor" id="introduccion"></a>

### 1.1. Antecedentes

Las plantas son sometidas constantemente a condiciones edáficas y climáticas no óptimas, que afectan a la homeostasis celular y que, en último térmico, impiden su crecimiento y su integridad$^1$. Para poder sobrevivir, las plantas deben responder y adaptarse de forma eficaz a las circunstancias. Esta adaptación puede venir dada por diferenciaciones específicas o por la modulación del comportamiento fisiológico en respuesta a una señal de estrés. Se han identificado cambios en el modelo de expresión de genes en respuesta a distintas condiciones de estrés (e.g. frío y estrés osmótico)$^2$, pero particularmente las condiciones anaerobias dan lugar a cambios rápidos en la expresión de genes relacionados con la respuesta adaptativa a bajos niveles de oxígeno. Estos efectos parecen ser atenuados cuando se les suministra azúcar exógeno (Figura 1). 

Para poder comprender cómo tenía lugar a nivel transcriptómico esta atenuación, en la literatura anterior al estudio que se validará a continuación (Loreti, Poggi, Novi, Alpi & Perata, 2005) se puede encontrar que se usaron microarrays, pero las sondas utilizadas se limitaban al estudio de la expresión diferencial de no más de 3500 genes relacionados con el metabolismo, factores de traducción y de transcripción, que no daban cuenta de los efectos a nivel del transcriptoma global. 

Por ello, con el fin de conocer con más profundidad el rol de los azúcares en la adaptación a la anaerobiosis, en este artículo se hace uso del chip de microarrays de Affymetrix ATH1. Con él, se pudo determinar el efecto a nivel de expresión del genoma global de la anaerobiosis y del suministro de sacarosa exógena a plantas en esta circunstancia. Más aún, con experimentos auxiliares ajenos a la transcriptómica se validaron los resultados. La disponibilidad de azúcares exógenos aumenta la tolerancia a la anoxia, y este efecto no puede conseguirse con otros azúcares como glucosa (Figura 1).

 <img src="/WorkingDirectory/Data/imagen1.jpg" alt="imagen" width="400"/>

<span style="font-size:0.85em">**Figura 1.** Efectos de la sacarosa exógena en la supervivencia de las plántulas. Estas fueron germinadas durante 4 días con o sin sacarosa según el caso, tras lo que se transfirieron a una cámara de anoxia durante 1, 2 o 3 días. Figura tomada de Loreti et al. (2005).</span>

### 1.2. Objetivo del estudio

A continuación, se realiza una validación de los datos de los experimentos de transcriptómica de este estudio con el fin de comprobar si la adición de sacarosa atenúa los efectos producidos por hipoxia.

### 1.3. Diseño experimental

Se usan homogeneizados de plántulas de *Arabidopsis thaliana* ecotipo Columbia *glabra* para evitar la variabilidad biológica debida a las diferenciaciones inherentes a los tejidos de la planta. Se incubaron durante 2 días en oscuridad a 4 ºC  y se transfirieron a una cámara a 23ºC durante 4 días antes de los tratamientos. En este punto, en el caso que procediera, se añadió sacarosa 90 mM, monitorizando la cantidad de azúcar con el fin de que este persistiera a lo largo de todo el experimento. Acto seguido, a las plántulas que correspondiera se las incubó en una cámara cerrada de anerobiosis durante 6 horas, que es el tiempo más apropiado porque representa la respuesta molecular  a la anaerobiosis prolongada antes de que se de la muerte inducida por anoxia. Para el estudio por microarrays se usan 3 condiciones, con 2 réplicas efectuadas independientemente. En cada una de ellas se hicieron 4 cultivos de plántulas que se juntaron para el protocolo de extracción de RNA. De esta forma, al final del proceso, para los estudios de transcriptómica existen dos réplicas para cada una de estos tres experimentos:

- Plántulas en aerobiosis sin sacarosa
- Plántulas en anaerobiosis sin sacarosa
- Plántulas en anaerobiosis con sacarosa
- <font color=#811700 size=1.9%> *No se realizó ningún control para verificar los efectos de la sacarosa en condiciones de aerobiosis.*</font>

### 1.4. Flujo de trabajo

En suma, en este trabajo se validarán los datos de transcriptómica del artículo de Loreti et al. (2005), almacenados en el [accesion number GSE2133](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE2133 "GSE2133 Accesion Viewer") de la base de datos Gene Expression Omnibus [(GEO)](https://www.ncbi.nlm.nih.gov/geo/). Se analizará la calidad de los datos obtenidos, se normalizarán los datos para que sean comparables y se estimarán los niveles de expresión. A partir de ellos, se establecerán contrastes para poder, finalmente, conocer qué genes (y en la medida de lo posible, la función de las proteínas que codifiquen) son sobreexpresados o reprimidos en cada condición. Todos estos procesos se resumen en el siguiente flujo de trabajo:

![jobflow.png](jobflow.png)

De esta manera, en primer lugar se cargarán los datos descargados desde la plataforma Gene Expression Omnibus (GEO) del NCBI, accediendo al [accesion number GSE2133](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE2133 "GSE2133 Accesion Viewer"). *Se descargaron para el análisis los datos desde [aquí](https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE2133&format=file)*. Se leerán los datos brutos con *ReadAffy* del paquete [*affy*](https://bioconductor.org/packages/release/bioc/html/affy.html) Con la función *cdfName* se muestra el código que identifica el diseño de la placa: [ATH1-121501](http://bioconductor.org/packages/release/data/annotation/html/ath1121501.db.html) y evaluando la variable microarray.raw.data generada tras usar la función *ReadAffy* se muestra la información de la placa, como el número de genes (22810).


```R
suppressMessages(library(affy))
suppressMessages(library(simpleaffy))
suppressMessages(library(affyPLM))
suppressMessages(library(limma))
suppressMessages(library(annaffy))
suppressMessages(library(gplots))
suppressMessages(library(RColorBrewer))
suppressMessages(library(treemap))

microarray.raw.data <- ReadAffy(verbose=TRUE)

microarray.raw.data

cdfName(microarray.raw.data)

```

    Warning message:
    “package ‘AnnotationDbi’ was built under R version 3.6.1”
    Warning message:
    “package ‘IRanges’ was built under R version 3.6.1”
    Warning message:
    “package ‘S4Vectors’ was built under R version 3.6.1”


    1 reading /Users/ParaisoNoFiscal/Biologia Molecular de Sistemas - Jupyter/TAREA 1/GSM38613.CEL.gz ...instantiating an AffyBatch (intensity a 506944x6 matrix)...done.
    Reading in : /Users/ParaisoNoFiscal/Biologia Molecular de Sistemas - Jupyter/TAREA 1/GSM38613.CEL.gz
    Reading in : /Users/ParaisoNoFiscal/Biologia Molecular de Sistemas - Jupyter/TAREA 1/GSM38614.CEL.gz
    Reading in : /Users/ParaisoNoFiscal/Biologia Molecular de Sistemas - Jupyter/TAREA 1/GSM38615.CEL.gz
    Reading in : /Users/ParaisoNoFiscal/Biologia Molecular de Sistemas - Jupyter/TAREA 1/GSM38616.CEL.gz
    Reading in : /Users/ParaisoNoFiscal/Biologia Molecular de Sistemas - Jupyter/TAREA 1/GSM38617.CEL.gz
    Reading in : /Users/ParaisoNoFiscal/Biologia Molecular de Sistemas - Jupyter/TAREA 1/GSM38618.CEL.gz


    Warning message:
    “replacing previous import ‘AnnotationDbi::tail’ by ‘utils::tail’ when loading ‘ath1121501cdf’”
    Warning message:
    “replacing previous import ‘AnnotationDbi::head’ by ‘utils::head’ when loading ‘ath1121501cdf’”
    
    



    AffyBatch object
    size of arrays=712x712 features (20 kb)
    cdf=ATH1-121501 (22810 affyids)
    number of samples=6
    number of genes=22810
    annotation=ath1121501
    notes=



'ATH1-121501'


## 2. Análisis de calidad
<a class="anchor" id="calidad"></a>

En todos los experimentos existen dos fuentes principales de variabilidad: la variabilidad biológica, que es de nuestro interés de estudio, y la variabilidad experimental que viene dada por limitaciones en los factores técnicos. Este ruido debe ser corregido en la medida de lo posible con un preprocesamiento de los datos brutos. Para conocer en qué extensión la variabilidad experimental está afectando a nuestros datos, procedemos a elaborar un control de calidad de los datos. 

En primer lugar, se procede a una visualización del estado general de los microarrays. En términos generales, se observa un buen estado (Figura 2), mas allá de unas marcas en el centro inherentes al diseño de los *chips* (Figura 3). Es destacable la diferencia de color entre las dos placas obtenidas con el experimento en anaerobiosis con sacarosa, lo que podría deberse a procesos de degradación en la placa de la izquierda.


```R
par(mfrow = c(3,2))
par(mar = c(2, 7, 2, 2))
image(microarray.raw.data[,1],col=rainbow(100), main = "", ylab = "- Sucrose 
Aerobiosis", cex.lab=1.55)
image(microarray.raw.data[,2],col=rainbow(100), main = "")
image(microarray.raw.data[,3],col=rainbow(100), main = "", ylab = "- Sucrose
Anaerobiosis", cex.lab=1.55)
image(microarray.raw.data[,4],col=rainbow(100), main = "")
image(microarray.raw.data[,5],col=rainbow(100), main = "", ylab = "+ 90 mM Sucrose
Anaerobiosis", cex.lab=1.55)
image(microarray.raw.data[,6],col=rainbow(100), main = "")
```


![png](output_5_0.png)


<span style="font-size:0.85em">**Figura 2.** </span> Aspecto general de los chips usados para los experimentos de microarrays.


```R
par(mfrow = c(1,1))
image(microarray.raw.data[,1], col = rainbow(100), main = "- Sucrose     Aerobiosis", cex.main = 1.5)
```


![png](output_7_0.png)


<span style="font-size:0.85em">**Figura 3.** </span> Detalle del microarray correspondiente al fichero GSM38613, en condiciones de aerobiosis sin sacarosa. Se aprecia un buen estado de la placa más allá del prominente daño del centro propio del diseño mismo de la placa. Los distintos colores indican diferente intensidad de fluorescencia.

Para poder conocer más profundamente el estado de las muestras y poder identificar procesos de degradación del RNA, se procede a un control de calidad implementado dentro del paquete [simppleaffy](https://bioconductor.org/packages/release/bioc/html/simpleaffy.html)de Bioconductor de los datos crudos:


```R
plot(qc(microarray.raw.data))
```


![png](output_9_0.png)


<span style="font-size:0.85em"> **Figura 4. Resumen del control de calidad.** Se muestra, en la primera columna, el nombre del fichero de cada muestra. En la segunda columna, se representa el porcentaje de detección (sondas para las que se detecta fluorescencia) y la fluorescencia del fondo. En la tercera columna, se representan algunos valores relativos a posibles procesos de degradación del RNA. Los círculos y triángulos representan la proporción de fluorescencia entre los entremos 3' y 5' de los genes para *gapdh* y *actin*, respectivamente. Se espera que el extremo 3' esté sobreexpresado con respecto al extremo 5' por los propios procesos de síntesis de cDNA y degradación de RNA. No obstante, con procesos de degradación masiva se aprecia que el extremo 3' está muy altamente sobreexpresado con respecto al 5'. En <span style="color:blue">azul</span> aparecen las medidas que superan los estándares de calidad y en <span style="color:red">rojo</span> los que no.</span>

En este control de calidad, se aprecia un porcentaje de detección homogéneo y una fluorescencia de fondo con valores muy diferentes entre las distintas muestras. Los valores para este último parámetro abarcan un gran rango en las distintas muestras y de hecho aparecen como potencialmente problemáticos, pero pueden ser fácilmente normalizados y corregidos por el algoritmo descrito en el apartado 3 (preprocesamiento de los datos). Los valores para el factor de proporcionalidad (la relación de fluorescencia de la celda con respecto a la de menor fluorescencia) son adecuados. 

No obstante, aparecen valores potencialmente problemáticos para las sondas de control de degradación en la muestra GSM38617 (muestra en anaerobiosis y con sacarosa) y para GSM38613. Los avisos de este último pueden  omitirse porque los valores de degradación para la sonda de la actina (el otro control) son correctos. Por otra parte, el RNA de la muestra del fichero GSM38617 (Anaerobiosis +Suc) parece haber sufrido un grave proceso de degradación y debe ser descartada para continuar con el estudio. Este proceso degradativo es visible en las imágenes de las placas (Figura 2), estando ausente en la otra réplica. No se mencionan estos potenciales defectos en la calidad de las muestras en el artículo. La degradación de ambos controles implica procesos degradativos profundos, y consecuentemente no se tendrán en cuenta los datos respectivos a esta réplica en la extracción de genes diferencialmente expresados ([apartado 4](#expresados)).

A continuación, con el fin de ilustrar que las intensidades de fluorescencia para las distintas muestras es normalmente muy distinta y que por tanto se requiere un procesamiento de estos datos antes de trabajar con ellos en los contrastes de hipótesis, se representa la intensidad de las distintas muestras mediante un diagrama de cajas y bigotes y un histograma, descriptores globales de la distribución de las intensidades.


```R
boxplot(microarray.raw.data,col=rainbow(8),las=2,ylab="Luminescence", ylim = c(5, 16), cex.lab = 1.4, xaxt = "n")
    text(x=1.5,y=16,"Aerobiosis",cex = 1.4, col = "dark blue")
    text(x=1.5,y=15,"- Sucrose",cex = 1.4)
    text(x=4.5,y=16,"Anaerobiosis",cex = 1.4, col = "dark blue")
    text(x=3.5,y=15,"- Sucrose",cex = 1.4)
    text(x=5.5,y=15,"+ Sucrose",cex = 1.4)
    lines(x=c(0.7,2.3), y=c(15.5,15.5), type = "l", lwd = 2)
    lines(x=c(2.7,6.3), y=c(15.5,15.5), type = "l", lwd = 2)
    lines(x=c(0.7,2.3), y=c(14.5,14.5), type = "l", lwd = 2)
    lines(x=c(2.7,4.3), y=c(14.5,14.5), type = "l", lwd = 2)
    lines(x=c(4.7,6.3), y=c(14.5,14.5), type = "l", lwd = 2)

hist(microarray.raw.data,col=rainbow(8), cex.lab = 1.4, ylab = "Density")
legend(x=11.5, y=0.5, legend = AffyRNAdeg(microarray.raw.data)$sample.names, col=rainbow(8), pch = 19)
```


![png](output_11_0.png)



![png](output_11_1.png)


<span style="font-size:0.85em">**Figura 5.** Diagrama de cajas y bigotes (arriba) e histograma (abajo) de los datos crudos de microarrays. Las muestras han sido nombradas como GSM386 y un identificador del 13 al 18; las muestras 13 y 14 corresponden a las muestras en aerobiosis sin sacarosa; 15, 16, 17 y 18 en anaerobiosis, siendo la 15 y 16 sin sacarosa y la 17 y 18 con sacarosa.</span>

En ambos gráficos se aprecia una gran diferencia entre las intensidades registradas. Más aún, los datos referentes a la muestra GSM38617 dan intensidades aparentemente menores, siendo esto congruente con la imagen de la placa (Figura 2) y el análisis de calidad (Figura 4). Deben, por tanto, normalizarse todos los datos de las intensidades antes de pasar a hacer contrastes de hipótesis. 

## 3. Procesamiento de los datos
<a class="anchor" id="preprocesamiento"></a>

### 3.1 Preprocesamiento de los datos

Para que todas las muestras sean comparables, se procesan los datos brutos mediante un *Robust Multiarray Average* (RMA), un algoritmo que realiza una corrección o sustracción de la fluorescencia de fondo, normalización de la fluorescencia de cada sonda y una estimación de los niveles de expresión de los genes representados por las distintas sondas del microarray utilizado en *log2*. A continuación, se representan con un *boxplot* y un *histograma* estos datos con el fin de conocer si son comparables.


```R
microarray.processed.data <- rma(microarray.raw.data)
boxplot(microarray.processed.data,col=rainbow(8),las=2,ylab="Luminescence", ylim = c(0, 16), cex.lab = 1.4, xaxt = "n")
    text(x=1.5,y=16,"Aerobiosis",cex = 1.4, col = "dark blue")
    text(x=1.5,y=15,"- Sucrose",cex = 1.4)
    text(x=4.5,y=16,"Anaerobiosis",cex = 1.4, col = "dark blue")
    text(x=3.5,y=15,"- Sucrose",cex = 1.4)
    text(x=5.5,y=15,"+ Sucrose",cex = 1.4)
    lines(x=c(0.7,2.3), y=c(15.5,15.5), type = "l", lwd = 2)
    lines(x=c(2.7,6.3), y=c(15.5,15.5), type = "l", lwd = 2)
    lines(x=c(0.7,2.3), y=c(14.5,14.5), type = "l", lwd = 2)
    lines(x=c(2.7,4.3), y=c(14.5,14.5), type = "l", lwd = 2)
    lines(x=c(4.7,6.3), y=c(14.5,14.5), type = "l", lwd = 2)

hist(microarray.processed.data,col=rainbow(8), cex.lab = 1.4, ylab = "Density", xlab = "Intensity")
legend(x=10, y=0.25, legend = AffyRNAdeg(microarray.raw.data)$sample.names, col=rainbow(8), pch = 19)
```

    Background correcting
    Normalizing
    Calculating Expression



![png](output_13_1.png)



![png](output_13_2.png)


<span style="font-size:0.85em">**Figura 6.** Diagrama de cajas y bigotes (arriba) e histograma (abajo) de los datos ya preprocesados de microarrays. Las muestras han sido nombradas como GSM386 y un identificador del 13 al 18; las muestras 13 y 14 corresponden a las muestras en aerobiosis sin sacarosa; 15, 16, 17 y 18 en anaerobiosis, siendo la 15 y 16 sin sacarosa y la 17 y 18 con sacarosa.</span>

Con el preprocesmiento de los datos se aprecia homogeneidad en las intensidades de todos las condiciones y réplicas, y ya pueden utilizarse para los contrastes y conocer la expresión diferencial de genes. Antes, no obstante, deben estimarse los niveles de expresión para cada uno de ellos.

### 3.2 Estimación de los niveles de expresión

El RMA, además de la normalización de los datos de fluorescencia, permite obtener una estimación de los niveles de expresión (en *log2*), pudiendo extraerse como una matriz usando la función *exprs* del paquete [affy](https://bioconductor.org/packages/release/bioc/html/affy.html) de Bioconductor. Tendrá tantas columnas como ficheros de datos (6 columnas) y tantas filas como sondas (22810). Por defecto, el nombre de las columnas es el nombre del fichero .CEL, por lo que se renombrarán para que sean más informativos con la condición a la que corresponden y el número de la réplica.


```R
expression.level <- exprs(microarray.processed.data)
sampleID <- c("Aer_-Suc1","Aer_-Suc2","Ana_-Suc1","Ana_-Suc2","Ana_+Suc1","Ana_+Suc2")
colnames(expression.level) <- sampleID
head(expression.level)
```


<table>
<caption>A matrix: 6 × 6 of type dbl</caption>
<thead>
	<tr><th></th><th scope=col>Aer_-Suc1</th><th scope=col>Aer_-Suc2</th><th scope=col>Ana_-Suc1</th><th scope=col>Ana_-Suc2</th><th scope=col>Ana_+Suc1</th><th scope=col>Ana_+Suc2</th></tr>
</thead>
<tbody>
	<tr><th scope=row>244901_at</th><td>4.020316</td><td>5.672301</td><td>5.038545</td><td>6.296530</td><td>6.535084</td><td>5.967054</td></tr>
	<tr><th scope=row>244902_at</th><td>4.021045</td><td>5.446902</td><td>5.074191</td><td>6.593599</td><td>5.861396</td><td>5.988343</td></tr>
	<tr><th scope=row>244903_at</th><td>4.580053</td><td>5.391182</td><td>6.116024</td><td>6.748804</td><td>7.250701</td><td>5.972113</td></tr>
	<tr><th scope=row>244904_at</th><td>4.822901</td><td>4.984814</td><td>5.147224</td><td>5.233352</td><td>6.233128</td><td>5.009946</td></tr>
	<tr><th scope=row>244905_at</th><td>4.189475</td><td>3.976261</td><td>4.759150</td><td>3.979296</td><td>5.280234</td><td>4.222614</td></tr>
	<tr><th scope=row>244906_at</th><td>5.394555</td><td>6.789372</td><td>6.285982</td><td>7.034621</td><td>7.284569</td><td>6.436245</td></tr>
</tbody>
</table>



En primer lugar, se calculan los valores medios de expresión para cada condición y se almacenan estos datos en una matriz. Debe matizarse que los potenciales procesos degradativos señalados en la Figura 4 pueden llevar a falsas conclusiones. Por ello, y a falta de réplicas adicionales, se asumirá que la media corresponde únicamente a la segunda réplica de la condición en anaerobiosis con sacarosa. Las conclusiones extraidas de estos resultados deben extraerse con extrema cautela, a la espera de poder replicar el experimento en un futuro.


```R
aer.no.suc <- (expression.level[,"Aer_-Suc1"]+expression.level[,"Aer_-Suc2"])/2
ana.no.suc <- (expression.level[,"Ana_-Suc1"]+expression.level[,"Ana_-Suc2"])/2
ana.wi.suc <- expression.level[,"Ana_+Suc2"]

mean.expression<-matrix(data = c(aer.no.suc,ana.no.suc,ana.wi.suc), nrow = nrow(expression.level), ncol = 3)
colnames(mean.expression)<-c("aer.no.suc", "ana.no.suc", "ana.wi.suc")
rownames(mean.expression)<-rownames(expression.level)
head(mean.expression)
```


<table>
<caption>A matrix: 6 × 3 of type dbl</caption>
<thead>
	<tr><th></th><th scope=col>aer.no.suc</th><th scope=col>ana.no.suc</th><th scope=col>ana.wi.suc</th></tr>
</thead>
<tbody>
	<tr><th scope=row>244901_at</th><td>4.846309</td><td>5.667538</td><td>5.967054</td></tr>
	<tr><th scope=row>244902_at</th><td>4.733974</td><td>5.833895</td><td>5.988343</td></tr>
	<tr><th scope=row>244903_at</th><td>4.985618</td><td>6.432414</td><td>5.972113</td></tr>
	<tr><th scope=row>244904_at</th><td>4.903858</td><td>5.190288</td><td>5.009946</td></tr>
	<tr><th scope=row>244905_at</th><td>4.082868</td><td>4.369223</td><td>4.222614</td></tr>
	<tr><th scope=row>244906_at</th><td>6.091963</td><td>6.660301</td><td>6.436245</td></tr>
</tbody>
</table>



Es interesante obtener una previsualización con el fin de determinar si el efecto global de los tratamientos tienen un efecto inhibidor o activador sobre el transcriptoma. Para realizar una previsualización comparativa, se usan gráficos de dispersión. 

En este caso, nos limitaremos a comparar el único control (aerobiosis, sin sacarosa) con los dos tratamientos y los dos tratamientos entre sí. 

1. Al comparar la situación en aerobiosis y en anaerobiosis sin sacarosa, se comprueban los efectos que tiene la falta de oxígeno en las plántulas.
 
2. Al comparar la situación en anaerobiosis con y sin sacarosa, se comprueban los efectos de la adición de sacarosa en plantas sometidas a anaerobiosis.
 
3. Al comparar la situación en aerobiosis sin sacarosa y anaerobiosis con sacarosa, se puede estimar, si existe, el efecto atenuador de la adición de sacarosa en plantas sometidas a estrés por anaerobiosis.

Acerca de este último punto, conocer el efecto de la sacarosa en condiciones de normoxia es necesario para la medida en la que esta afecta a estas plantas de la misma manera que en el artículo se comprueba si hay efectos diferenciales entre las plantas crecidas en hipoxia cuando se les suministra sacarosa exógena. 

La última comparación (aerobiosis -sacarosa / anaerobiosis +sacarosa) no tiene sentido por sí sola, puesto que no podríamos conocer hasta qué punto cada una de las dos condiciones diferentes ha afectado al transcriptoma. Podemos evaluarla solo si tenemos los efectos de la anaerobiosis, obtenidos con el primer contraste. Curiosamente, en la Figura 1 se muestra el crecimiento de las plantas con sacarosa y sin sacarosa tanto en condiciones aerobias como en anoxia, pero no hay información acerca del transcriptoma de las plantas en condiciones aerobias con sacarosa. No se mencionan en el artículo razones para no haber realizado este experimento.


```R
par(mfrow=c(1,3))
plot(main = "Anaerobiosis effects", x = aer.no.suc, y = ana.no.suc,xlab="Aerobiosis - Sucrose",ylab="Anaerobiosis - Sucrose",pch=19,cex=0.5, cex.lab = 1.5)
plot(main = "Sucrose effects 
on anaerobiosis", x = ana.no.suc, y = ana.wi.suc,xlab="Anaerobiosis - Sucrose",ylab="Anaerobiosis + Sucrose",pch=19,cex=0.5, cex.lab = 1.5)
plot(main = "Effects of sucrose on
anaerobic conditions", x = aer.no.suc, y = ana.wi.suc,xlab="Aerobiosis - Sucrose",ylab="Anaerobiosis + Sucrose",pch=19,cex=0.5, cex.lab = 1.5)

```


![png](output_19_0.png)


<span style="font-size:0.85em">**Figura 7.** Diagrama de dispersión de los datos de expresión media, donde cada punto representa un transcrito. Puntos que estén sobre la diagonal no se expresan de forma diferencial. Aquellos puntos que queden por encima son sobreexpresados en la segunda condición respecto a la primera. Contrariamente, aquellos puntos que quedan por debajo de la diagonal estarán reprimidos con respecto a la primera condición. </span>

En condiciones de anaerobiosis, no se aprecia un efecto muy sustancial cuando se añade sacarosa (Figura 7, centro) si se compara con los efectos que tiene la anaerobiosis (Figura 7, izquierda). No es fácilmente distinguible el efecto que tiene la sacarosa sobre las condiciones anaerobias (Figura 7, derecha), aunque es apreciable un aparente descenso de la sobreexpresión de algunos genes y una mayor represión de otros. Para poder determinar este efecto, y conocer en profundidad qué genes están afectados, se debe hacer una selección de los genes expresados de forma diferencial (DEGs).

## 4. Análisis de expresión génica diferencial
<a class="anchor" id="expresados"></a>

[LIMMA](https://bioconductor.org/packages/release/bioc/html/limma.html) (*LInear Models for Microarray Analysis*) es un paquete de Bioconductor para el análisis de análisis de datos de expresión génica de microarrays, concretamente para la evaluación de la expresión diferencial. Para poder proceder a e la selección de genes diferencialmente expresados, primero se genera una matriz que represente el diseño experimental, asignando a cada condición experimental un número entero:

1. Aerobiosis Sin sacarosa
2. Anaerobiosis Sin sacarosa
3. Anaerobiosis Con sacarosa (90 mM)

Se asignará uno de estos 3 números a cada una de las muestras con su etiqueta correspondiente, usando la función *model.matrix* contenida en el paquete LIMMA. A continuación, atendiendo al diseño experimental, con la función *lmFit*, se ajusta la estimación de los niveles de expresión de cada gen a un modelo lineal. Esto permitirá obtener básicamente la media de las réplicas en cada condición.

*En el análisis de calidad (apartado 2), se señaló que la muestra GSM38617 puede haber sufrido graves procesos de degradación. Por ello, solo se tendrán en cuenta los datos de la segunda réplica (GSM38618), con el fin de poder continuar con el estudio de la expresión diferencial.*


```R
expression.level <- cbind(expression.level, expression.level[,6])
expression.level <- expression.level[,-5]
sampleID2 <- c("Aer_-Suc1","Aer_-Suc2","Ana_-Suc1","Ana_-Suc2","Ana_+Suc2","Ana_+Suc2_Bis")
colnames(expression.level) <- sampleID2

head(expression.level)

experimental.design <- model.matrix(~ -1+factor(c(1,1,2,2,3,3)))
colnames(experimental.design) <- c("Aer_no_Suc","Ana_no_Suc","Ana_wi_Suc")
linear.fit <- lmFit(expression.level, experimental.design)
```


<table>
<caption>A matrix: 6 × 6 of type dbl</caption>
<thead>
	<tr><th></th><th scope=col>Aer_-Suc1</th><th scope=col>Aer_-Suc2</th><th scope=col>Ana_-Suc1</th><th scope=col>Ana_-Suc2</th><th scope=col>Ana_+Suc2</th><th scope=col>Ana_+Suc2_Bis</th></tr>
</thead>
<tbody>
	<tr><th scope=row>244901_at</th><td>4.020316</td><td>5.672301</td><td>5.038545</td><td>6.296530</td><td>5.967054</td><td>5.967054</td></tr>
	<tr><th scope=row>244902_at</th><td>4.021045</td><td>5.446902</td><td>5.074191</td><td>6.593599</td><td>5.988343</td><td>5.988343</td></tr>
	<tr><th scope=row>244903_at</th><td>4.580053</td><td>5.391182</td><td>6.116024</td><td>6.748804</td><td>5.972113</td><td>5.972113</td></tr>
	<tr><th scope=row>244904_at</th><td>4.822901</td><td>4.984814</td><td>5.147224</td><td>5.233352</td><td>5.009946</td><td>5.009946</td></tr>
	<tr><th scope=row>244905_at</th><td>4.189475</td><td>3.976261</td><td>4.759150</td><td>3.979296</td><td>4.222614</td><td>4.222614</td></tr>
	<tr><th scope=row>244906_at</th><td>5.394555</td><td>6.789372</td><td>6.285982</td><td>7.034621</td><td>6.436245</td><td>6.436245</td></tr>
</tbody>
</table>



A continuación, para determinar los genes que se expresan diferencialmente, se comparan (como se mencionó en la sección 3.2) el único control (aerobiosis, sin sacarosa) con los dos tratamientos y los dos tratamientos entre sí. Esto es:

* Aerobiosis y anaerobiosis sin sacarosa.
* Anaerobiosis con y sin sacarosa, se comprueban los efectos de la adición de sacarosa en plantas sometidas a anaerobiosis.
* Al comparar la situación en aerobiosis sin sacarosa y anaerobiosis con sacarosa, se puede estimar, si existe, el efecto atenuador de la adición de sacarosa en plantas sometidas a estrés por anaerobiosis.

Para especificar los contrastes se usa *makeContrasts*, que construye la matriz de contrastes. Esta corresponde a los contrastes del conjunto de datos que se especifican.


```R
contrast.matrix <- makeContrasts(Ana_no_Suc-Aer_no_Suc, Ana_wi_Suc-Ana_no_Suc, Ana_wi_Suc-Aer_no_Suc, levels=c("Aer_no_Suc","Ana_no_Suc", "Ana_wi_Suc"))
```

En el presente artículo usaron para los contrastes de hipótesis el criterio de *fold-change*, sin especificar qué umbral escogieron ya que se atuvieron a los criterios de determinación de la sobreexpresión o represión implantados en el software usado (Microarray Analysis Suite 5.0). El método del fold-change es adecuado al tratarse de organismos modelo y al disponer de un limitado número de réplicas para cada experimento independiente (más en nuestro caso, donde una de las réplicas debe ser descartada al haber sufrido graves procesos de degradación). Es por ello que a continuación se reproduce este mismo método. No obstante, a todos los efectos debe tenerse en cuenta que este método no tiene controles de falsos positivos, algo que fácilmente podría ser subsanado con un mayor número de réplicas y el cálculo de estadísticos apropiados.

Se calcula el fold-change y los p-valores de cada gen con las funciones *constrasts.fit* y *eBayes*. A continuación, con *topTable* se obtiene una tabla con información sobre la expresión diferencial de los genes en los contrastes asignados anteriomente a un número. Para poder conocer qué genes se expresan diferencialmente, se extrae el fold-change y el identificador de cada sonda.

Los genes expresados diferencialmente pueden estar sobreexpresados, teniendo entonces un fold-change mayor que el umbral o reprimidos, teniendo un fold-change menor que el umbral. Se escogerá un umbral de 2, puesto que:

* los genes que superen este umbral (fc > log$_2$(2) => fc > 1) estarían expresados más del dos veces con respecto al fenotipo sin tratar. 

* De forma similar, los genes que sean menores que el umbral (fc < -log$_2$(2) => fc < -1) estarían expresados menos de la mitad con respecto al fenotipo sin tratar.

La función DEGs.calling recibe, a modo de comprobación, la comparación que se está evaluando (*condition.id*), los resultados del cálculo usando las funciones constrasts.fit y eBayes del fold-change y los p-valores correspondientes para cada gen en cada uno de los constrastes especificados, almacenados estos resultados en *contrast.results*. Se deben especificar el número de genes de la placa de microarrays (*gene.number* = 22810 en el caso de ATH1), así como las anotaciones específicas de esta (*chip* = ath1121501.db). Es necesario establecer un umbral (*log2threshold*). Adicionalmente, es necesario especificar el coeficiente de los contrastes del modelo experimental (*coef*). Además, se pueden añadir parámetros para modificar el gráfico. Finalmente, la función *DEGs.calling* da lugar al número de genes activados y reprimidos para cada comparación, un gráfico de *scatterplot* que recoge estos datos y almacena en txt y html los genes activados o reprimidos.

Estos genes que han sido seleccionados según su fold-change se representarán en una tabla, usando la función *aafTableAnn*, del paquete [annaffy](https://bioconductor.org/packages/release/bioc/html/annaffy.html) y se guardarán en formato HTML con el nombre indicado en *output.activated.txt/html* (es similar con los reprimidos). La función *aafTableAnn* recibirá un vector que almacena el nombre de los genes activados o reprimidos y los datos de las anotaciones del chip de microarrays que se estén utilizando. En este caso, se usa el [Affymetrix Arabidopsis ATH1 Genome Array](http://bioconductor.org/packages/release/data/annotation/html/ath1121501.db.html), cuyos datos se almacenan en la base de datos del chip ath1121501.


```R
contrast.linear.fit <- contrasts.fit(linear.fit, contrast.matrix)
contrast.results <- eBayes(contrast.linear.fit)

DEGs.calling <- function(condition.id, contrast.results, gene.number, coef, log2threshold, title.plot,
                        control.mean.expression, treat.mean.expression, xlab, ylab, pch = 19, cex = 0.5,
                        cex.lab = 1.5, col = "grey", chip, output.activated.txt, output.activated.html,
                        output.repressed.txt, output.repressed.html)
{
    dif.exp.info <- topTable(contrast.results, number = gene.number, coef, sort.by="logFC")
    genes.ids    <- rownames(dif.exp.info)
    fold.change  <- dif.exp.info[["logFC"]]
    names(x = fold.change) <- genes.ids
    
    ## Activated genes
    
    activated.genes  <- genes.ids[fold.change > log2threshold]
    length_activated <- length(activated.genes)
    
    ## Repressed genes
    
    repressed.genes  <- genes.ids[fold.change < -log2threshold]
    length_repressed <- length(repressed.genes)
    
    ## Result
    
    result <- list(activated.genes = activated.genes, repressed.genes = repressed.genes,
                length_activated = length_activated, length_repressed = length_repressed)
    
    ## Plot
    
    plot(main = title.plot, x = control.mean.expression, y = treat.mean.expression,xlab = xlab, 
         ylab = ylab, pch=19,cex=0.5, cex.lab = 1.5, col="grey")
    
    points(x = control.mean.expression[activated.genes],
       y = treat.mean.expression[activated.genes],pch=19,cex=0.5,col="red")
    points(control.mean.expression[repressed.genes],
       treat.mean.expression[repressed.genes],pch=19,cex=0.5,col="blue")
    
    ## Storage of DEGs data in .txt and .html format
    
        # Activated genes
        
        activated.genes.table <- aafTableAnn(probeids = activated.genes, chip = chip, aaf.handler())
    
        saveText(activated.genes.table, file = output.activated.txt)
        saveHTML(activated.genes.table, file = output.activated.html)
    
        # Repressed genes
        
        repressed.genes.table <- aafTableAnn(probeids = repressed.genes, chip = chip, aaf.handler())
    
        saveText(repressed.genes.table, file = output.repressed.txt)
        saveHTML(repressed.genes.table, file = output.repressed.html)
    
    return(list(length_activated = length_activated, length_repressed = length_repressed))
    
}

## PRIMER CONTRASTE: AEROBIOSIS -SACAROSA / ANAEROBIOSIS -SACAROSA

DEGs.calling(condition.id = "Aer.no.Ana.no", contrast.results = contrast.results,
            gene.number = 22810, coef = 1, log2threshold = 1, title.plot = "Anaerobiosis effects",
            control.mean.expression = aer.no.suc, treat.mean.expression = ana.no.suc, 
            xlab="Aerobiosis -Sucrose",ylab="Anaerobiosis -Sucrose", pch = 19, cex = 0.5,
            cex.lab = 1.5, col = "grey", chip = "ath1121501.db", output.activated.txt = "activated_genes_Aer_no_Ana_no.txt",
            output.activated.html = "activated_genes_Aer_no_Ana_no.html", 
            output.repressed.txt = "repressed_genes_Aer_no_Ana_no.txt", 
            output.repressed.html = "repressed_genes_Aer_no_Ana_no.html")

text(aer.no.suc["252746_at"]+0.15,ana.no.suc["252746_at"]+0.3,"SUS1", col="black", cex=0.7)
text(aer.no.suc["245998_at"]+0.3,ana.no.suc["245998_at"]-0.3,"SUGT", col="black", cex=0.7)
text(aer.no.suc["253416_at"]+0.3,ana.no.suc["253416_at"]-0.3,"PDC1", col="black", cex=0.7)
text(aer.no.suc["264953_at"]+0.3,ana.no.suc["264953_at"]-0.3,"ADH", col="black", cex=0.7)
text(aer.no.suc["248138_at"]+0.3,ana.no.suc["248138_at"]-0.3,"PDC2", col="black", cex=0.7)
text(aer.no.suc["260847_s_at"]+0.3,ana.no.suc["260847_s_at"]-0.3,"AAT", col="black", cex=0.7)


## SEGUNDO CONTRASTE: ANAEROBIOSIS -SACAROSA / ANAEROBIOSIS +SACAROSA

DEGs.calling(condition.id = "Ana.no.Ana.wi", contrast.results = contrast.results,
            gene.number = 22810, coef = 2, log2threshold = 1, title.plot = "Sucrose effects on anaerobiosis",
            control.mean.expression = ana.no.suc, treat.mean.expression = ana.wi.suc, 
            xlab="Anaerobiosis -Sucrose",ylab="Anaerobiosis +Sucrose", chip = "ath1121501.db",
            output.activated.txt = "activated_genes_Ana_no_Ana_wi.txt",
            output.activated.html = "activated_genes_Ana_no_Ana_wi.html", 
            output.repressed.txt = "repressed_genes_Ana_no_Ana_wi.txt", 
            output.repressed.html = "repressed_genes_Ana_no_Ana_wi.html")

text(ana.no.suc["252746_at"]+0.15,ana.wi.suc["252746_at"]+0.3,"SUS1", col="black", cex=0.7)
text(ana.no.suc["245998_at"]+0.3,ana.wi.suc["245998_at"]-0.3,"SUGT", col="black", cex=0.7)
text(ana.no.suc["253416_at"]+0.3,ana.wi.suc["253416_at"]-0.3,"PDC1", col="black", cex=0.7)
text(ana.no.suc["264953_at"]+0.3,ana.wi.suc["264953_at"]-0.3,"ADH", col="black", cex=0.7)
text(ana.no.suc["248138_at"]+0.3,ana.wi.suc["248138_at"]-0.3,"PDC2", col="black", cex=0.7)


## TERCER CONTRASTE: AEROBIOSIS -SACAROSA / ANAEROBIOSIS +SACAROSA

DEGs.calling(condition.id = "Aer.no.Ana.wi", contrast.results = contrast.results,
            gene.number = 22810, coef = 3, log2threshold = 1, title.plot = "Effects of sucrose on
anaerobic conditions", control.mean.expression = aer.no.suc, treat.mean.expression = ana.wi.suc, 
            xlab="Aerobiosis -Sucrose",ylab="Anaerobiosis +Sucrose", chip = "ath1121501.db", 
            output.activated.txt = "activated_genes_Aer_no_Ana_wi.txt",
            output.activated.html = "activated_genes_Aer_no_Ana_wi.html", 
            output.repressed.txt = "repressed_genes_Aer_no_Ana_wi.txt", 
            output.repressed.html = "repressed_genes_Aer_no_Ana_wi.html")

text(aer.no.suc["252746_at"]+0.15,ana.wi.suc["252746_at"]+0.3,"SUS1", col="black", cex=0.7)
text(aer.no.suc["245998_at"]+0.3,ana.wi.suc["245998_at"]-0.3,"SUGT", col="black", cex=0.7)
text(aer.no.suc["253416_at"]+0.3,ana.wi.suc["253416_at"]-0.3,"PDC1", col="black", cex=0.7)
text(aer.no.suc["264953_at"]+0.3,ana.wi.suc["264953_at"]-0.3,"ADH", col="black", cex=0.7)
text(aer.no.suc["248138_at"]+0.3,ana.wi.suc["248138_at"]-0.3,"PDC2", col="black", cex=0.7)


```

    Loading required package: ath1121501.db
    
    Loading required package: org.At.tair.db
    
    
    
    
    
    Warning message in chkPkgs(chip):
    “The ath1121501.db package does not appear to contain annotation data.”
    Warning message in result_fetch(res@ptr, n = n):
    “SQL statements must be issued with dbExecute() or dbSendStatement() instead of dbGetQuery() or dbSendQuery().”
    Warning message in result_fetch(res@ptr, n = n):
    “SQL statements must be issued with dbExecute() or dbSendStatement() instead of dbGetQuery() or dbSendQuery().”



<dl>
	<dt>$length_activated</dt>
		<dd>424</dd>
	<dt>$length_repressed</dt>
		<dd>839</dd>
</dl>




![png](output_25_2.png)



<dl>
	<dt>$length_activated</dt>
		<dd>250</dd>
	<dt>$length_repressed</dt>
		<dd>269</dd>
</dl>




![png](output_25_4.png)



<dl>
	<dt>$length_activated</dt>
		<dd>520</dd>
	<dt>$length_repressed</dt>
		<dd>890</dd>
</dl>




![png](output_25_6.png)


<span style="font-size:0.85em">**Figura 7.** Diagrama de dispersión de los datos de expresión media, donde cada punto representa un transcrito. Puntos que estén sobre la diagonal no se expresan de forma diferencial. Se han seleccionado aquellos puntos que queden por encima y que superan el umbral 2 del fold-change (<span style="color:red">rojo</span>), es decir, que son sobreexpresados en la segunda condición respecto a la primera. Contrariamente, aquellos puntos que quedan por debajo de la diagonal y son menores que el umbral -2 del fold-change (<span style="color:blue">azul</span>) estarán reprimidos con respecto a la primera condición. </span>

* **Aerobiosis y anaerobiosis sin sacarosa.**

Las condiciones anaerobias tienen un efecto global represor (se reprime la expresión a menos de la mitad en 839 genes), aunque algunos de los genes que están sobreexpresados en esta situación (424) aparezcan mucho más representados.

* **Anaerobiosis con y sin sacarosa, se comprueban los efectos de la adición de sacarosa en plantas sometidas a anaerobiosis.**

La sacarosa ejerce sobre el transcriptoma de plantas anaerobias un efecto moderado, activando y reprimiendo aproximadamente a unos 250 genes.


```R
par(mfrow=c(1,2))

## PRIMER CONTRASTE: AEROBIOSIS -SACAROSA / ANAEROBIOSIS -SACAROSA

DEGs.calling(condition.id = "Aer.no.Ana.no", contrast.results = contrast.results,
            gene.number = 22810, coef = 1, log2threshold = 2, title.plot = "Anaerobiosis effects",
            control.mean.expression = aer.no.suc, treat.mean.expression = ana.no.suc, 
            xlab="Aerobiosis -Sucrose",ylab="Anaerobiosis -Sucrose", pch = 19, cex = 0.5,
            cex.lab = 1.5, col = "grey", chip = "ath1121501.db", 
            output.activated.txt = "activated_genes_Aer_no_Ana_no_th2.txt",
            output.activated.html = "activated_genes_Aer_no_Ana_no_th2.html", 
            output.repressed.txt = "repressed_genes_Aer_no_Ana_no_th2.txt", 
            output.repressed.html = "repressed_genes_Aer_no_Ana_no_th2.html")

## TERCER CONTRASTE: AEROBIOSIS -SACAROSA / ANAEROBIOSIS +SACAROSA

DEGs.calling(condition.id = "Aer.no.Ana.wi", contrast.results = contrast.results,
            gene.number = 22810, coef = 3, log2threshold = 2, title.plot = "Effects of sucrose on
anaerobic conditions", control.mean.expression = aer.no.suc, treat.mean.expression = ana.wi.suc, 
            xlab="Aerobiosis -Sucrose",ylab="Anaerobiosis +Sucrose", chip = "ath1121501.db", 
            output.activated.txt = "activated_genes_Aer_no_Ana_wi_th2.txt",
            output.activated.html = "activated_genes_Aer_no_Ana_wi_th2.html", 
            output.repressed.txt = "repressed_genes_Aer_no_Ana_wi_th2.txt", 
            output.repressed.html = "repressed_genes_Aer_no_Ana_wi_th2.html")
```


<dl>
	<dt>$length_activated</dt>
		<dd>100</dd>
	<dt>$length_repressed</dt>
		<dd>123</dd>
</dl>




<dl>
	<dt>$length_activated</dt>
		<dd>101</dd>
	<dt>$length_repressed</dt>
		<dd>129</dd>
</dl>




![png](output_27_2.png)


<span style="font-size:0.85em">**Figura 8.** Diagrama de dispersión de los datos de expresión media, donde cada punto representa un transcrito. Puntos que estén sobre la diagonal no se expresan de forma diferencial. Se han seleccionado aquellos puntos que queden por encima y que superan el umbral 4 del fold-change (<span style="color:red">rojo</span>), es decir, que son sobreexpresados en la segunda condición respecto a la primera. Contrariamente, aquellos puntos que quedan por debajo de la diagonal y son menores que el umbral -4 del fold-change (<span style="color:blue">azul</span>) estarán reprimidos con respecto a la primera condición. </span>


* **Al comparar la situación en aerobiosis sin sacarosa y anaerobiosis con sacarosa, se puede estimar, si existe, el efecto atenuador de la adición de sacarosa en plantas sometidas a estrés por anaerobiosis.**

Cuando se añade sacarosa a plantas anaerobias y se compara su transcriptoma con el de plantas aerobias, se aprecia que la sacarosa aparentemente aumenta la sobreexpresión de algunos genes (520, en comparación con los 424 activados en anaerobiosis únicamente) y la represión de otros (890, en comparación con los 839 activados en anaerobiosis) (Figura 7).

Si se aumenta el umbral del fold-change a 2 (los genes para considerarse sobreexpresados deben estar expresados 4 veces más que en el fenotipo normal y para ser reprimidos deben estarlo menos de 4 veces que el fenotipo normal), se aprecia una diferencia muy pequeña que podría deberse a efectos espúreos y no a la adición de sacarosa (Figura 8). De esta manera, la adición de sacarosa no tiene un efecto atenuador global, sino que afectará específicamente a distintos genes.

---

Se han generado los documentos txt y html que almacenan cada una de las sondas, con su correspondiente símbolo si existiera, una breve descripción del gen, el cromosoma donde se encuentra y su localización en cada uno de ellos, el nombre del gen, una referencia bibliográfica, su función y la ruta donde participa. Pueden verse o descargarse directamente desde aquí:

### Genes activados

 .| Primer contraste | Segundo contraste | Tercer contraste
-- | -- | -- | --
HTML | [Ver](activated_genes_Aer_no_Ana_no.html) | [Ver](activated_genes_Ana_no_Ana_wi.html) | [Ver](activated_genes_Aer_no_Ana_wi.html)
Text | [Descargar](activated_genes_Aer_no_Ana_no.txt) | [Descargar](activated_genes_Ana_no_Ana_wi.txt) | [Descargar](activated_genes_Aer_no_Ana_wi.txt)

### Genes reprimidos

  .| Primer contraste | Segundo contraste | Tercer contraste
-- | -- | -- | --
HTML | [Ver](repressed_genes_Aer_no_Ana_no.html) | [Ver](repressed_genes_Ana_no_Ana_wi.html) | [Ver](repressed_genes_Aer_no_Ana_wi.html)
Text | [Descargar](repressed_genes_Aer_no_Ana_no.txt) | [Descargar](repressed_genes_Ana_no_Ana_wi.txt) | [Descargar](repressed_genes_Aer_no_Ana_wi.txt)

---

A continuación, podría ser útil representar en un diagramma de Venn aquellos genes diferencialmente expresados. Idealmente, podrían representarse los genes sobreexpresados (o reprimidos) por la presencia de sacarosa en condiciones aerobias y anaerobias. En la intersección de los dos círculos correspondientes a estos genes aparecerían genes que se sobreexpresan en ambos casos, independientemente de las condiciones de oxígeno; en el resto de ambos círculos quedarían los genes diferencialmente expresados únicamente en cada una de las condiciones por la presencia de sacarosa. No obstante, el diseño experimental de los estudios de transcriptómica de este artículo no incluía un control positivo para la adición de azúcar exógeno, de tal manera que la elaboración de este diagrama no es posible. 

## Análisis del mapa de calor

La elaboración de un mapa de calor permite la representación de los distintos genes diferencialmente expresados (DEGs) y su distinto nivel de expresión en cada una de las condiciones. Para ello, el primer paso que se lleva a cabo es el almacenamiento de los genes diferencialmente expresados en todas las condiciones en un vector donde no haya repeticiones de valores. Se extraen los niveles de expresión de los DEGs en cada una de las condiciones y se normaliza con la función *scale*. Finalmente, con *heatmap.2* del paquete [gplots](https://cran.r-project.org/web/packages/gplots/index.html), de [CRAN](https://cran.r-project.org/) se elabora el mapa de calor. Los dendrogramas a la izquierda del mapa de calor muestra cómo las distintas condiciones y los genes están agrupados independientemente. Los patrones en el mapa de calor, que muestra los valores para la expresión de cada condición y gen, indicarían una posible correlación. 


```R
DEGs_act_rep <- function(contrast.results, gene.number, coef, log2threshold)
{
 activated <- list()
 repressed <- list()
 j <- 1
    for (i in 1:length(coef))
    {
    dif.exp.info <- topTable(contrast.results, number = gene.number, coef = i, sort.by="logFC")
    genes.ids    <- rownames(dif.exp.info)
    fold.change  <- dif.exp.info[["logFC"]]
    names(x = fold.change) <- genes.ids
    
    ## Activated genes
    
    activated.genes  <- genes.ids[fold.change > log2threshold]
    length_activated <- length(activated.genes)
    
    ## Repressed genes
    
    repressed.genes  <- genes.ids[fold.change < -log2threshold]
    length_repressed <- length(repressed.genes)
    
    activated[[j]] <- activated.genes
    repressed[[j]] <- repressed.genes
        
    j <- j+1
    }
    
    list(activated = activated, repressed = repressed)
}

result <- DEGs_act_rep(contrast.results = contrast.results, gene.number = 22810, coef = c(1,2,3), log2threshold = 1)

activated <- result$activated
repressed <- result$repressed

complete.DEGs <- unique(c(activated[[1]], repressed[[1]], activated[[2]], repressed[[2]], activated[[3]], repressed[[3]]))
colnames(mean.expression) <- c("AER -Suc", "ANA -Suc", "ANA +Suc")
DEG.expression <- mean.expression[complete.DEGs,c("AER -Suc", "ANA -Suc", "ANA +Suc")]
normalized.DEG.expression <- t(scale(t(DEG.expression)))

heatmap.2(normalized.DEG.expression,Colv=FALSE,dendrogram="row",
          labRow=c(""),density.info="none",trace="none",
          col=heat.colors(100)[100:1],margins = c(8,8),cexCol=1.2)
```


![png](output_29_0.png)


<span style="font-size:0.85em">**Figura 9.** Mapa de calor de los genes expresados diferencialmente en cada una de las tres condiciones: aerobiosis sin sacarosa (AER -Suc), anaerobiosis sin (ANA -Suc) y con sacarosa (ANA +Suc). Los dendrogramas agrupan agrupaciones independientes de los datos de expresión génica en cada condición. Se muestra una leyenda con los colores que designan valores entre -1 y +1 siendo -1 el estado más reprimido y +1 el más sobreexpresado. </span>

En el mapa de calor (Figura 9) se aprecia que un gran número de genes expresados en aerobiosis están agrupados independientemente, indicando una posible sobreexpresión en relación a los primeros. Esto es fácilmente distinguible atendiendo a las ramificaciones del dendrograma, que se divide al inicio en tres ramas. El clado superior se relaciona con genes sobreexpresados en aerobiosis con respecto a la anaerobiosis. Este se subdivide en diferentes grupos. Algunos de ellos, ubicados en la parte superior del mapa de calor, tienen una expresión mayor en condiciones de anaerobiosis con sacarosa que sin este azúcar. No obstante, también existen ejemplos (en la parte media-inferior del primer clado) de una mayor expresión en condiciones aerobias que en anaerobiosis sin sacarosa, y este último mayor que cuando se añade sacarosa. Finalmente, hay un grupo de genes fuertemente reprimidos que son más activados en condiciones de anaerobiosis sin sacarosa que en anerobiosis.

En cuanto al clado central, se aprecian niveles de sobreexpresión de genes en la condición con sacarosa que en las otras dos. En este punto, es necesario reseñar que la ausencia de un control donde se monitorice la expresión diferencial de genes en condiciones aerobias con sacarosa impide que se pueda determinar con certeza si la sacarosa ha ejercido su efecto en respuesta a las condiciones anaerobias. Los niveles de expresión de estos genes sobreexpresados en anaerobiosis con sacarosa son, por otro lado, parecidos en las condiciones aerobia y anaerobia sin sacarosa, lo que podría indicar que la sacarosa ha estimulado la síntesis de genes relacionados con el metabolismo de la sacarosa, independientemente de la cantidad de oxígeno a las que se le someta a la planta.

Quizás, uno de los patrones más interesantes para nuestro estudio sean los datos de expresión correspondientes al último clado. En él, se muestran genes que están muy reprimidos en condiciones aerobias, pero que se expresan en mayor o menor medida en condiciones anaerobias. Si atendemos a la parte más inferior de este último clado, se aprecia sobreexpresión de genes en la condición anaerobia con sacarosa que están más atenuados. Hay una gran diversidad de explicaciones para este fenómeno, pero desde hace años se conoce la existencia de formación de supercomplejos macromoleculares en la membrana interna mitocondrial para mejorar la eficiencia del transporte electrónico mitocondrial bajo condiciones de hipoxia$^3$ (Figura 10). Estas agrupaciones han recibido el nombre de *respirasomas* y siguen siendo un tema de interés científico en campos como la medicina puesto que tras él pueden subyacer algunas claves para entender algunas patologías humanas como la isquemia$^4$. De la misma forma, bajo condiciones de hipoxia, es necesario optimizar el oxígeno usado en la respiración. Este metabolismo puede verse aumentado si se aporta de forma exógena un azúcar porque provocará un aumento de la tasa glucolítica y por tanto también de la oxidación por la ruta aerobia, de tal manera que podría hipotetizarse que este grupo de genes participarían en procesos metabólicos de este tipo.

<img src="gr3.jpg" alt="ya queda menos" width="400"/>

<span style="font-size:0.85em">**Figura 10.** Moléculas de citocromo *c* viajando entre los complejos respiratorios. A) Modelo fluido del CIII dimérico (verde) y CIV dimérico (morado) embebidos en su membrana mitocondrial. B) Modelo sólido de CIII dimérico y CIV monomérico. P y D indican loso sitios de unión del citocromo proximal y distal, respectivamente, en cada complejo. "Pool" denota la población de citocromos en el espacio intermembrana. Las imágenes fueron creadas con el software UCFS Chimera usando las estructuras recogidas en la PDB para CIII, CIV (5XTH) y citocromo *c* (2N9I). Figura y texto adaptado de Pérez-Mejías, G., Guerra-Castellano, A., Díaz-Quintana, A., De la Rosa, M., & Díaz-Moreno, I. (2019). Cytochrome c: Surfing Off of the Mitochondrial Membrane on the Tops of Complexes III and IV. *Computational And Structural Biotechnology Journal*, 17, 654-660. doi: 10.1016/j.csbj.2019.05.002 </span>

También es relevante en este último clado los genes de la parte superior, que aparecen muy representados en la condición de anaerobiosis sin azúcar pero con los mismos niveles de expresión en aerobiosis y cuando se aplica sacarosa a estas plantas sometidas a anaerobiosis. Estos genes serían los genes de interés de estudio, puesto que su sobreexpresión en plantas sometidas a anaerobiosis se ve reducida a niveles basales cuando se suministra sacarosa exógena.

## 5. Anotación funcional y discusión
<a class="anchor" id="discusion"></a>

Ya habiendo sido revisados superficialmente los datos de expresión, es necesario conocer profundamente qué genes particularmente son afectados en cada caso. Para ello, haremos uso de herramientas de ontología de genes. La ontología es el estudio de un objeto y sus relaciones dentro de un dominio del conocimiento. En este caso, la ontología de genes (GO, Gene Ontology) describe nuestro conocimiento del dominio biológico desde tres perspectivas:

* **Procesos biológicos.** Describen eventos que ocurren de principio a fin tales como la transducción de señales.
* **Componentes celulares.** Describe las localizaciones celulares en las que el producto génico realiza una función, bien en compartimentos celulares o en complejos macromoleculares estables en los que participe.
* **Funciones moleculares.** Describen actividades que ocurren a nivel molecular tales como actividades enzimáticas (e.g. kinasas).

"El análisis de enriquecimiento de conjunto de genes (GSEA, Gene Set Enrichment Analysis) es un método computacional que determina si un conjunto de genes muestra diferencias significativas, concordantes entre dos estados biológicos"$^5$. Haciendo un test exacto de Fisher puede estudiarse si existe enriquecimiento para cada GO de interés, pudiéndose determinar si en un conjunto de genes dados hay términos de GO que aparecen con significancia estadística con respecto al conjunto de genes que representan el universo total.

De forma parecida, se usará [GOrilla](http://cbl-gorilla.cs.technion.ac.il/#ref) para analizar el enriquecimiento de conjunto de genes los contrastes que se han indicado.

### Primer contraste: genes diferencialmente expresados en anaerobiosis

Para ello, se introducen en [GOrilla](http://cbl-gorilla.cs.technion.ac.il/#ref) manualmente la lista de genes activados recogidos en [activated_genes_Aer_no_Ana_no.txt](activated_genes_Aer_no_Ana_no.html) ([descargable](activated_genes_Aer_no_Ana_no.txt)) y los reprimidos [repressed_genes_Aer_no_Ana_no.txt](repressed_genes_Aer_no_Ana_no.html) ([descargable](repressed_genes_Aer_no_Ana_no.txt)).

#### *Genes activados*

<img src="activated_genes_Aer_no_Ana_no.png" alt="Activated Aer- Ana-">

En el diagrama, se representa con colores el p-valor, que abarca desde <10$^-$$^9$ (rojo) hasta >10$^-$$^3$ (blanco), aunque se haya restringido la búsqueda a p-valores < 10$^-$$^4$. En él, puede apreciarse que la anaerobiosis estimula la **respiración anaerobia**, la respuesta al **ácido abscísico** (hormona liberada en respuesta al estrés de la planta), y la **respuesta celular a la hipoxia**.

Con la herramienta [ReViGO](http://revigo.irb.hr/) (Reduce and Visualize GO) se elaboró un diagrama que agrupa los genes activados activados o reprimidos en función de su función biológica. En este caso, se modificó ligeramente el script original de ReViGO para obtener un archivo jpeg haciendo uso del script siguiente: [REVIGO_aer_no_ana_no_oe.r](./REVIGO_aer_no_ana_no_oe.r). Con el fin de simplificar, para los siguientes contrastes no se detalla el script, aunque se facilite el acceso a él desde un hipervínculo y se presente una figura con el resultado.

<img src="revigo_aer_no_ana_no_oe.jpeg" alt="Smiley face" height="960" width="720">


```R
revigo.names <- c("term_ID","description","freqInDbPercent","abslog10pvalue","uniqueness","dispensability","representative");
revigo.data <- rbind(
c("GO:0007623","circadian rhythm",0.539,5.7932,0.821,0.000,"circadian rhythm"),
c("GO:0009061","anaerobic respiration",0.047,5.2815,0.815,0.000,"anaerobic respiration"),
c("GO:0009987","cellular process",58.826,7.3468,0.971,0.000,"cellular process"),
c("GO:0048511","rhythmic process",0.600,6.6091,0.930,0.000,"rhythmic process"),
c("GO:0050896","response to stimulus",25.943,29.2993,0.948,0.000,"response to stimulus"),
c("GO:0070482","response to oxygen levels",0.285,55.9355,0.560,0.000,"response to oxygen levels"),
c("GO:0009628","response to abiotic stimulus",7.946,43.2916,0.572,0.344,"response to oxygen levels"),
c("GO:0071456","cellular response to hypoxia",0.117,56.8633,0.416,0.370,"response to oxygen levels"),
c("GO:0033993","response to lipid",3.267,4.1649,0.516,0.298,"response to oxygen levels"),
c("GO:0006950","response to stress",14.156,24.6861,0.546,0.497,"response to oxygen levels"),
c("GO:1901700","response to oxygen-containing compound",6.504,5.3372,0.497,0.644,"response to oxygen levels"),
c("GO:0070887","cellular response to chemical stimulus",5.287,49.4034,0.486,0.581,"response to oxygen levels"),
c("GO:0097305","response to alcohol",2.469,4.0110,0.524,0.689,"response to oxygen levels"),
c("GO:0009651","response to salt stress",2.270,4.0424,0.477,0.617,"response to oxygen levels"),
c("GO:0033554","cellular response to stress",3.772,36.6308,0.528,0.219,"response to oxygen levels"),
c("GO:0009408","response to heat",0.859,7.8570,0.515,0.439,"response to oxygen levels"),
c("GO:0051716","cellular response to stimulus",12.637,35.8539,0.542,0.434,"response to oxygen levels"),
c("GO:0042221","response to chemical",12.434,34.6655,0.552,0.481,"response to oxygen levels"),
c("GO:0009266","response to temperature stimulus",2.287,5.5834,0.486,0.546,"response to oxygen levels"),
c("GO:0009737","response to abscisic acid",2.443,4.0110,0.524,0.629,"response to oxygen levels"),
c("GO:0080167","response to karrikin",0.531,6.9172,0.541,0.466,"response to oxygen levels"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$abslog10pvalue <- as.numeric( as.character(stuff$abslog10pvalue) );
stuff$freqInDbPercent <- as.numeric( as.character(stuff$freqInDbPercent) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );
jpeg( filename = "revigo_aer_no_ana_no_oe.jpeg", width = 960, height =  720, 
     units = "px", pointsize = 12, quality = 75)

treemap(
    dtf = stuff,
    index = c("representative","description"),
    vSize = "abslog10pvalue",
    type = "categorical",
    vColor = "representative",
    title = "Overexpressed genes on Aerobiosis -Suc / Anaerobiosis -Suc contrast",
    fontsize.title = 25,
    fontsize.labels = 23,
    inflate.labels = FALSE,      
    lowerbound.cex.labels = 0,   
    bg.labels = "#CCCCCCAA",     
    position.legend = "none"
)

dev.off()
```


<strong>pdf:</strong> 2


Se aprecian prácticamente los mismos resultados que con el diagrama de GOrilla, estando muy representada la respuesta a condiciones de bajo oxígeno, respiración anaerobia y de respuesta a estímulos abióticos. Es importante remarcar que la expresión de genes en respuesta a hipoxia no es exclusiva de plantas: en animales, por ejemplo, el descubrimiento de los procesos que están involucrados en la adaptación celular a bajas concentraciones de oxígeno valió el premio Nobel del 2019 de Fisiología y Medicina a los investigadores William Kaelin Jr., Gregg L. Semenza y Peter Ratcliffe. No es de extrañar que estos sistemas sean ubicuos en muchos organismos puesto que de la disponibilidad de oxígeno es un factor que su propia supervivencia.

Son, por tanto, muy evidentes los profundos cambios en el modelo de expresión génica que sufren las plantas sometidas a anoxia.

<img src="GO_1.png" alt="absiakabamoya">

Dentro del GO:0071453, donde se encuentra la función de respuesta a niveles de oxígeno, se encuentran algunas proteínas de choque térmico (heat shock proteins, HSP), algunas de las cuales actúan como chaperonas ayudando al plegamiento correcto de las proteínas en condiciones de anoxia. Las proteínas de esta familia se encuentran sobrerrepresentadas en condiciones de anoxia. Algunas HSP son inducidas por hipoxia en animales, y un choque térmico moderado puede aumentar la resistencia a la hipoxia, como lo demuestran los autores tratando algunas de las plántulas durante 1.5 h antes de las 6 horas de anoxia incubándolas a 38 ºC. Este tratamiento aumentó la tolerancia a la ausencia de oxígeno, y la adición de sacarosa aumentó la producción de HSP (Figura 11).

<img src="hsp.png" alt="" height = "480" width = "360">

<span style="font-size:0.85em">**Figura 11.** Efectos del tratamiento previo sobre la tolerancia a la anoxia en plántulas de *Arabidopsis*. Las semillas se trataron durante 1.5 h a 38 ºC antes del tratamiento de anoxia. Figura tomada de Loreti, E. et al. (2005).</span>


#### *Genes reprimidos*

<img src="repressed_genes_Aer_no_Ana_no.png" alt="Activated Aer- Ana-">

<img src="revigo_aer_no_ana_no_re.jpeg" alt="" height="960" width="720">


```R

```

En condiciones de hipoxia, se activan genes de la ruta de señalización de la auxina (principalmente las proteínas pertenecientes a la familia de proteínas parecidas a SAUR de respuesta a auxina), fitohormonas reguladoras del crecimiento de los tejidos vegetales. En este contexto, es coherente que la planta detenga su crecimiento en condiciones de estrés. Estos genes fueron seleccionados por los autores del artículo para tratar de comprender si la inhibición de estos genes era menor cuando se añadía sacarosa a las plántulas (Figura 12), como se detalla posteriormente.

### Segundo contraste: efectos del azúcar exógeno sobre plantas sometidas a anaerobiosis

Se pueden descargar los archivos correspondientes a los [genes activados](activated_genes_Ana_no_Ana_wi.html) [aquí](activated_genes_Ana_no_Ana_wi.txt) y los [genes reprimidos](repressed_genes_Ana_no_Ana_wi.html) [aquí](repressed_genes_Ana_no_Ana_wi.txt). 

#### *Genes activados*

<img src="revigo_ana_no_ana_wi_oe.jpeg" alt="" height="960" width="720">

Un bajo número (250) de genes están sobreexpresados cuando se añade sacarosa a las plántulas. El efecto a nivel génico que este azúcar tiene está de alguna manera relacionado con la respuesta a estímulos abióticos, aumentando su respuesta al estrés.

#### *Genes reprimidos*

<img src="repressed_genes_Ana_no_Ana_wi.png" alt="Activated Ana- Ana+" height="600" width="400">

<img src="revigo_ana_no_ana_wi_re.jpeg" alt="" height="960" width="720">


Por otro lado, en cuanto a la represión de genes, se ve cómo se afectan directamente los genes de respuesta a niveles de oxígeno que se activaban cuando se sometía a las plántulas a anaerobiosis. De esta manera, se observa que hay una inhibición muy clara de los genes de respuesta a hipoxia, disminuyendo así mismo su respuesta a los estímulos.

En el artículo, se asegura que la adición de sacarosa "reduce eficazmente los efectos negativos de la anoxia" en el grupo de genes involucrados en la fisiología de la auxina (Figura 12). Sin embargo, no se encuentran GOs reprimidas significativamente relacionadas con estos procesos.

<img src="GO_2.png" alt="absiakabamoya">

Por otro lado, se activan genes del GO:0050896 como un gran número de peroxidasas (AT1G49570, AT5G64120, AT5G64110, etc.), factores de transcripción, etc. Particularmente, AT3G58450 corresponde a una proteína de la familia de proteínas universales del estrés. Estos resultados apuntan a que la adición de sacarosa exógena aumenta la resistencia de la planta a la anoxia al reprimir genes relacionados con la respuesta a niveles bajos de oxígeno y al estrés.

### Tercer contraste: diferencias de expresión génica entre plantas sometidas a aerobiosis sin azúcar exógeno y sometidas a anaerobiosis con azúcar exógeno.

Se pueden descargar los archivos correspondientes a los [genes activados](activated_genes_Aer_no_Ana_wi.html) [aquí](activated_genes_Aer_no_Ana_wi.txt) y los [genes reprimidos](repressed_genes_Aer_no_Ana_wi.html) [aquí](repressed_genes_Aer_no_Ana_wi.txt). 

#### *Genes activados*

<img src="activated_genes_Aer_no_Ana_wi.png" alt="">

<img src="revigo_aer_no_ana_wi_oe.jpeg" alt="" height="960" width="720">

<img src="GO_3.png" alt="">

La respuesta celular desencadenada por las condiciones de anoxia debe ser algo más leve que la que tienen lugar cuando no se suministra a las plantas sacarosa (basada esta suposición en el p-valor del test usado). De hecho, tiene lugar la activación de los mismos GO, pero cuando se añade sacarosa los p-valores disminuyen drásticamente, dando a entender que su sobreexpresión es más moderada.

#### *Genes reprimidos*

<img src="repressed_genes_Aer_no_Ana_wi.png" alt="">

<img src="revigo_aer_no_ana_wi_re.jpeg" alt="" height="960" width="720">

La adición de sacarosa a plántulas sometidas a hipoxia inhibe al expresión de genes relacionados con las vías de señalización activadas por auxina y a la respuesta a estímulos, como ocurría cuando no se añadía sacarosa. No obstante, la represión es menor cuando se añade sacarosa, como se explica en el artículo (Figura 12).

<img src="sucrose_auxin.png" alt="" height="480" width="360">

<span style="font-size:0.85em">**Figura 12.** Genes significativamente reprimidos por la anoxia y anotados como genes codificantes de proteínas relacionadas con la fisiología de la auxina fueron seleccionados, y su cambio de expresión (*fold change*) bajo anoxia fue representado contra el cambio de expresión en condiciones de anoxia con sacarosa. Figura tomada de Loreti, E. et al. (2005).</span>


## Discusión

En el artículo se revisan los aspectos metabólicos de las plantas sometidas a anaerobiosis: metabolismo de lípidos (la), de sacarosa, el metabolismo fermentativo y de algunos aminoácidos como alanina, la cadena transportadora de electrones mitocondrial, e incluso algunos aspectos como la fisiología de la auxina o la expresión de proteínas de shock térmico. A continuación se presentan los *scatterplot* de la figura 7 sobre la que se han añadido los nombres de:

Affymetrix code | AGI CODE | Description | Abbreviation
-- | -- | -- | --
260869_at |	At1g43800 | stearoyl acyl carrier protein desaturase, putative similar to stearoyl acyl carrier protein desaturase | SACPD
252746_at |	At3g43190	|  sucrose synthase -like protein SUCROSE SYNTHASE (SUCROSE-UDP GLUCOSYLTRANSFERASE) | SUS1
245998_at |	At5g20830 | sucrose-UDP glucosyltransferase | SUGT
264953_at |	At1g77120 |	alcohol dehydrogenase | ADH
253416_at |	At4g33070 |	pyruvate decarboxylase-1 | PDC1
248138_at |	At5g54960 |	pyruvate decarboxylase  | PDC2
260847_s_at	| At1g17290 |	alanine aminotransferase, putative similar to alanine aminotransferase | AAT
244951_s_at	| ATMG00180	|	cytochrome c biogenesis orf452 | Cc
257337_at |	AtMg00060	|	NADH dehydrogenase subunit 5 (nad5) | NAD5
258879_at |	At3g03270 |		universal stress protein (USP) family protein / early nodulin ENOD18 family protein | USP
252515_at |	At3g46230 |	heat shock protein 17  | HSP17
258930_at | At3g10040 |	unknown protein predicted by genscan | OxySt1?
264968_at |	At1g67360 |	stress related protein, putative similar to stress related protein GI:5802955 | OxySt2?
263150_at |	At1g54050 |	heat-shock protein, putative similar to heat-shock protein GI:472939 | HSPOx



```R
## PRIMER CONTRASTE: AEROBIOSIS -SACAROSA / ANAEROBIOSIS -SACAROSA

DEGs.calling(condition.id = "Aer.no.Ana.no", contrast.results = contrast.results,
            gene.number = 22810, coef = 1, log2threshold = 1, title.plot = "Anaerobiosis effects",
            control.mean.expression = aer.no.suc, treat.mean.expression = ana.no.suc, 
            xlab="Aerobiosis -Sucrose",ylab="Anaerobiosis -Sucrose", pch = 19, cex = 0.5,
            cex.lab = 1.5, col = "grey", chip = "ath1121501.db", output.activated.txt = "activated_genes_Aer_no_Ana_no.txt",
            output.activated.html = "activated_genes_Aer_no_Ana_no.html", 
            output.repressed.txt = "repressed_genes_Aer_no_Ana_no.txt", 
            output.repressed.html = "repressed_genes_Aer_no_Ana_no.html")

text(aer.no.suc["260869_at"]+0.15,ana.no.suc["260869_at"]+0.3,"SACPD", col="black", cex=0.7)
text(aer.no.suc["252746_at"]-0.15,ana.no.suc["252746_at"]+0.3,"SUS1", col="black", cex=0.7)
text(aer.no.suc["245998_at"]+0.15,ana.no.suc["245998_at"]+0.3,"SUGT", col="black", cex=0.7)
text(aer.no.suc["264953_at"]+0.15,ana.no.suc["264953_at"]+0.3,"ADH", col="black", cex=0.7)
text(aer.no.suc["253416_at"]+0.15,ana.no.suc["253416_at"]+0.3,"PDC1", col="black", cex=0.7)
text(aer.no.suc["248138_at"]+0.15,ana.no.suc["248138_at"]+0.3,"PDC2", col="black", cex=0.7)
text(aer.no.suc["260847_s_at"]+0.15,ana.no.suc["260847_s_at"]+0.3,"AAT", col="black", cex=0.7)
text(aer.no.suc["244951_s_at"]+0.15,ana.no.suc["244951_s_at"]+0.3,"Cc", col="black", cex=0.7)
text(aer.no.suc["257337_at"]+0.15,ana.no.suc["257337_at"]+0.3,"NAD5", col="black", cex=0.7)
text(aer.no.suc["258879_at"]+0.15,ana.no.suc["258879_at"]+0.3,"USP", col="black", cex=0.7)
text(aer.no.suc["252515_at"]+0.15,ana.no.suc["252515_at"]-0.1,"HSP17", col="black", cex=0.7)
text(aer.no.suc["258930_at"]+0.15,ana.no.suc["258930_at"]+0.3,"OxySt1", col="chartreuse4", cex=0.7)
text(aer.no.suc["264968_at"]+0.15,ana.no.suc["264968_at"]+0.3,"OxySt2", col="chartreuse4", cex=0.7)
text(aer.no.suc["263150_at"]+0.15,ana.no.suc["263150_at"]+0.3,"HSPOx", col="chartreuse4", cex=0.7)

## SEGUNDO CONTRASTE: ANAEROBIOSIS -SACAROSA / ANAEROBIOSIS +SACAROSA

DEGs.calling(condition.id = "Ana.no.Ana.wi", contrast.results = contrast.results,
            gene.number = 22810, coef = 2, log2threshold = 1, title.plot = "Sucrose effects on anaerobiosis",
            control.mean.expression = ana.no.suc, treat.mean.expression = ana.wi.suc, 
            xlab="Anaerobiosis -Sucrose",ylab="Anaerobiosis +Sucrose", chip = "ath1121501.db",
            output.activated.txt = "activated_genes_Ana_no_Ana_wi.txt",
            output.activated.html = "activated_genes_Ana_no_Ana_wi.html", 
            output.repressed.txt = "repressed_genes_Ana_no_Ana_wi.txt", 
            output.repressed.html = "repressed_genes_Ana_no_Ana_wi.html")

text(ana.no.suc["260869_at"]+0.15,ana.wi.suc["260869_at"]+0.3,"SACPD", col="black", cex=0.7)
text(ana.no.suc["252746_at"]+0.15,ana.wi.suc["252746_at"]+0.3,"SUS1", col="black", cex=0.7)
text(ana.no.suc["245998_at"]+0.15,ana.wi.suc["245998_at"]+0.3,"SUGT", col="black", cex=0.7)
text(ana.no.suc["264953_at"]+0.15,ana.wi.suc["264953_at"]+0.3,"ADH", col="black", cex=0.7)
text(ana.no.suc["253416_at"]+0.15,ana.wi.suc["253416_at"]+0.3,"PDC1", col="black", cex=0.7)
text(ana.no.suc["248138_at"]+0.15,ana.wi.suc["248138_at"]+0.3,"PDC2", col="black", cex=0.7)
text(ana.no.suc["260847_s_at"]+0.15,ana.wi.suc["260847_s_at"]+0.3,"AAT", col="black", cex=0.7)
text(ana.no.suc["244951_s_at"]+0.15,ana.wi.suc["244951_s_at"]+0.3,"Cc", col="black", cex=0.7)
text(ana.no.suc["257337_at"]+0.15,ana.wi.suc["257337_at"]+0.3,"NAD5", col="black", cex=0.7)
text(ana.no.suc["258879_at"]+0.15,ana.wi.suc["258879_at"]+0.3,"USP", col="black", cex=0.7)
text(ana.no.suc["252515_at"]+0.15,ana.wi.suc["252515_at"]+0.3,"HSP17", col="black", cex=0.7)
text(ana.no.suc["258930_at"]+0.15,ana.wi.suc["258930_at"]+0.3,"OxySt1", col="chartreuse4", cex=0.7)
text(ana.no.suc["264968_at"]+0.15,ana.wi.suc["264968_at"]+0.3,"OxySt2", col="chartreuse4", cex=0.7)
text(ana.no.suc["263150_at"]+0.15,ana.wi.suc["263150_at"]+0.3,"HSPOx", col="chartreuse4", cex=0.7)


## TERCER CONTRASTE: AEROBIOSIS -SACAROSA / ANAEROBIOSIS +SACAROSA

DEGs.calling(condition.id = "Aer.no.Ana.wi", contrast.results = contrast.results,
            gene.number = 22810, coef = 3, log2threshold = 1, title.plot = "Effects of sucrose on
anaerobic conditions", control.mean.expression = aer.no.suc, treat.mean.expression = ana.wi.suc, 
            xlab="Aerobiosis -Sucrose",ylab="Anaerobiosis +Sucrose", chip = "ath1121501.db", 
            output.activated.txt = "activated_genes_Aer_no_Ana_wi.txt",
            output.activated.html = "activated_genes_Aer_no_Ana_wi.html", 
            output.repressed.txt = "repressed_genes_Aer_no_Ana_wi.txt", 
            output.repressed.html = "repressed_genes_Aer_no_Ana_wi.html")

text(aer.no.suc["260869_at"]+0.15,ana.wi.suc["260869_at"]+0.3,"SACPD", col="black", cex=0.7)
text(aer.no.suc["252746_at"]+0.15,ana.wi.suc["252746_at"]+0.3,"SUS1", col="black", cex=0.7)
text(aer.no.suc["245998_at"]+0.15,ana.wi.suc["245998_at"]+0.3,"SUGT", col="black", cex=0.7)
text(aer.no.suc["264953_at"]+0.15,ana.wi.suc["264953_at"]+0.3,"ADH", col="black", cex=0.7)
text(aer.no.suc["253416_at"]+0.15,ana.wi.suc["253416_at"]+0.3,"PDC1", col="black", cex=0.7)
text(aer.no.suc["248138_at"]+0.15,ana.wi.suc["248138_at"]+0.3,"PDC2", col="black", cex=0.7)
text(aer.no.suc["260847_s_at"]+0.15,ana.wi.suc["260847_s_at"]+0.3,"AAT", col="black", cex=0.7)
text(aer.no.suc["244951_s_at"]+0.15,ana.wi.suc["244951_s_at"]+0.3,"Cc", col="black", cex=0.7)
text(aer.no.suc["257337_at"]+0.15,ana.wi.suc["257337_at"]+0.3,"NAD5", col="black", cex=0.7)
text(aer.no.suc["258879_at"]+0.15,ana.wi.suc["258879_at"]+0.3,"USP", col="black", cex=0.7)
text(aer.no.suc["252515_at"]+0.15,ana.wi.suc["252515_at"]+0.3,"HSP17", col="black", cex=0.7)
text(aer.no.suc["258930_at"]+0.15,ana.wi.suc["258930_at"]+0.3,"OxySt1", col="chartreuse4", cex=0.7)
text(aer.no.suc["264968_at"]+0.15,ana.wi.suc["264968_at"]+0.3,"OxySt2", col="chartreuse4", cex=0.7)
text(aer.no.suc["263150_at"]+0.15,ana.wi.suc["263150_at"]+0.3,"HSPOx", col="chartreuse4", cex=0.7)

```


<dl>
	<dt>$length_activated</dt>
		<dd>424</dd>
	<dt>$length_repressed</dt>
		<dd>839</dd>
</dl>




![png](output_42_1.png)



<dl>
	<dt>$length_activated</dt>
		<dd>250</dd>
	<dt>$length_repressed</dt>
		<dd>269</dd>
</dl>




![png](output_42_3.png)



<dl>
	<dt>$length_activated</dt>
		<dd>520</dd>
	<dt>$length_repressed</dt>
		<dd>890</dd>
</dl>




![png](output_42_5.png)


<span style="font-size:0.85em">**Figura 8.** Diagrama de dispersión de los datos de expresión media, donde cada punto representa un transcrito. Puntos que estén sobre la diagonal no se expresan de forma diferencial. Se han seleccionado aquellos puntos que queden por encima y que superan el umbral 2 del fold-change (<span style="color:red">rojo</span>), es decir, que son sobreexpresados en la segunda condición respecto a la primera. Contrariamente, aquellos puntos que quedan por debajo de la diagonal y son menores que el umbral -2 del fold-change (<span style="color:blue">azul</span>) estarán reprimidos con respecto a la primera condición. Los nombres de las proteínas OxySt1, OxySt2 y HSPOx están escritas en color <span style="color:#009548">verde</span> porque aún se desconoce con certeza su función, aunque se tenga conocimiento de que intervienen en el estrés por anoxia. </span>

La adición de sacarosa mitiga en algunos genes la expresión inducida por anoxia, pero no así en otras proteínas. La adición de sacarosa provoca la represión de algunos genes que intervienen en la respuesta a oxígeno como [At3g10040](https://www.uniprot.org/uniprot/Q8RWY5/publications), cuya función sigue sin ser clara más allá de que es un factor de transcripción o genes que participan en el metabolismo fermentativo (ADH). Los mecanismos de la modulación de la transcripción por azúcar exógeno aún están en investigación. 

## 6. Conclusión
<a class="anchor" id="conclusion"></a>

Los estudios del transcriptoma a genoma completo suponen un gran avance en biología molecular en tanto que permiten conocer los efectos de un tratamiento a nivel global, si bien la técnica de microarrays se limita a sondas conocidas y otras limitaciones que técnicas más novedosas como la secuenciación masiva de RNA (RNAseq) han superado con creces. Ninguna de estas técnicas tienen en cuenta la posible regulación de la eficiencia de la traducción de los transcritos, pero aportan información de expresión de transcritos que frecuentemente se correlacionan con los datos de expresión de proteínas. 

En este contexto, la adición de sacarosa exógena a plántulas de *Arabidopsis* tiene un efecto ligeramente atenuador de genes relacionados con la respuesta a condiciones de oxígeno bajo como pueden ser algunas proteínas de la familia de HSP. No obstante, este efecto atenuador es muy leve y apenas afecta a 250 genes comparado con el fuerte efecto que suponen las condiciones de anoxia en estas plantas. Muchos de estos genes deben aún ser estudiados con profundidad para poder comprender qué efecto tiene su expresión en la célula.

Por otro lado, el tratamiento con un choque térmico de 1.5 h a 38 ºC hace que las plántulas desarrollen una mayor tolerancia a las condiciones de hipoxia. Cuando se les añade sacarosa exógena, esta resistencia aumenta quizás porque este azúcar favorece la transcripción o la acumulación de transcritos de proteínas de choque térmico, que actúan como chaperonas en condiciones de anoxia. 

## 7. Bibliografía
<a class="anchor" id="bibliografia"></a>

Se validaron los resultados del artículo:

**Loreti, E., Poggi, A., Novi, G., Alpi, A., & Perata, P. (2005). A Genome-Wide Analysis of the Effects of Sucrose on Gene Expression in Arabidopsis Seedlings under Anoxia. *Plant Physiology*, 137(3), 1130-1138. doi: 10.1104/pp.104.057299**

Adicionalmente, se consultaron:

1. Mickelbart, M., Hasegawa, P., & Bailey-Serres, J. (2015). Genetic mechanisms of abiotic stress tolerance that translate to crop yield stability. *Nature Reviews Genetics*, 16(4), 237-251. doi: 10.1038/nrg3901

2. Kumar, D., Hazra, S., Datta, R., & Chattopadhyay, S. (2016). Transcriptome analysis of Arabidopsis mutants suggests a crosstalk between ABA, ethylene and GSH against combined cold and osmotic stress. *Scientific Reports*, 6(1). doi: 10.1038/srep36867

3. Pérez-Mejías, G., Guerra-Castellano, A., Díaz-Quintana, A., De la Rosa, M., & Díaz-Moreno, I. (2019). Cytochrome c: Surfing Off of the Mitochondrial Membrane on the Tops of Complexes III and IV. *Computational And Structural Biotechnology Journal*, 17, 654-660. doi: 10.1016/j.csbj.2019.05.002

4. Guerra-Castellano, A., Díaz-Quintana, A., Pérez-Mejías, G., Elena-Real, C., González-Arzola, K., & García-Mauriño, S. et al. (2018). Oxidative stress is tightly regulated by cytochrome c phosphorylation and respirasome factors in mitochondria. *Proceedings Of The National Academy Of Sciences*, 115(31), 7955-7960. doi: 10.1073/pnas.1806833115

5. Subramanian, A., Tamayo, P., Mootha, V., Mukherjee, S., Ebert, B., & Gillette, M. et al. (2005). Gene set enrichment analysis: A knowledge-based approach for interpreting genome-wide expression profiles. Proceedings Of The National Academy Of Sciences, 102(43), 15545-15550. doi: 10.1073/pnas.0506580102

**Transparencias del curso de biología molecular de sistemas, Francisco J. Romero Campero e Ignacio Pérez Hurtado de Mendoza, CC BY-NC-ND 3.0**

Para el análisis bioinformático, se usaron las siguientes herramientas y paquetes:

6. Eran Eden*, Roy Navon*, Israel Steinfeld, Doron Lipson and Zohar Yakhini. "GOrilla: A Tool For Discovery And Visualization of Enriched GO Terms in Ranked Gene Lists", [*BMC Bioinformatics* 2009, 10:48](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-10-48). 
7. Supek, F., Bošnjak, M., Škunca, N., & Šmuc, T. (2011). REVIGO Summarizes and Visualizes Long Lists of Gene Ontology Terms. *Plos ONE*, 6(7), e21800. doi: 10.1371/journal.pone.0021800
8. Gautier, L., Cope, L., Bolstad, B. M., and Irizarry, R. A. 2004. affy---analysis of Affymetrix GeneChip data at the probe level. Bioinformatics 20, 3 (Feb. 2004), 307-315.
9. Bolstad, BM (2004) Low Level Analysis of High-density Oligonucleotide Array Data: Background, Normalization and Summarization. Dissertation. University of California, Berkeley.
10.  Bolstad BM, Collin F, Brettschneider J, Simpson K, Cope L, Irizarry RA, and Speed TP. (2005) Quality Assessment of Affymetrix GeneChip Data in Bioinformatics and Computational Biology Solutions Using R and Bioconductor. Gentleman R, Carey V, Huber W, Irizarry R, and Dudoit S. (Eds.), Springer, New York.
11. Brettschneider J, Collin F, Bolstad BM, and Speed TP. (2007) Quality assessment for short oligonucleotide arrays. Technometrics. In press.
12. Crispin J Miller (2019). simpleaffy: Very simple high level analysis of Affymetrix data. http://www.bioconductor.org, http://bioinformatics.picr.man.ac.uk/simpleaffy/.
13. Ritchie, M.E., Phipson, B., Wu, D., Hu, Y., Law, C.W., Shi, W., and Smyth, G.K. (2015). limma powers differential expression analyses for RNA-sequencing and microarray studies. Nucleic Acids Research 43(7), e47.
14. Colin A. Smith (2019). annaffy: Annotation tools for Affymetrix biological metadata. R package version 1.56.0.
15. Carlson M (2016). ath1121501.db: Affymetrix Arabidopsis ATH1 Genome Array annotation data (chip ath1121501). R package version 3.2.3.

**Se recogen los paquetes para su sencilla instalación en el script [README](./readme.r) y el script completo en [script.r](script.r).**
