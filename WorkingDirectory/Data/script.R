################################################################
#                                                              #
#                            TAREA 1                           #
#                     MICROARRAYS - SCRIPT                     #
#                                                              #
#    Biología Molecular de Sistemas - Grado en Bioquímica      #    
#     Joaquín Tamargo Azplicueta - joatamazp@alum.us.es        #
#                                                              #
################################################################

## En este script solo se recogen las funciones usadas para la
## elaboración de esta tarea. Se describen con detalle en el 
## texto las imágenes y gráficos obtenidos.

## Se cargan los datos del chip desde GSE2133, accesible desde
## (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE2133)

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

##################### ANÁLISIS DE CALIDAD ######################

## Evaluación del estado general de los chips. Se observa buen
## estado general, aunque el color tan distinto de los chips
## correspondientes a la condición de anaerobiosis con sacarosa
## pueda deberse a procesos de degradación. 

image(microarray.raw.data[,1],col=rainbow(100), main = "", ylab = "- Sucrose 
Aerobiosis", cex.lab=1.55)
image(microarray.raw.data[,2],col=rainbow(100), main = "")
image(microarray.raw.data[,3],col=rainbow(100), main = "", ylab = "- Sucrose
Anaerobiosis", cex.lab=1.55)
image(microarray.raw.data[,4],col=rainbow(100), main = "")
image(microarray.raw.data[,5],col=rainbow(100), main = "", ylab = "+ 90 mM Sucrose
Anaerobiosis", cex.lab=1.55)
image(microarray.raw.data[,6],col=rainbow(100), main = "")

## Elaboración de boxplot e histograma de las intensidades en 
## cada placa. Se aprecia gran dispersión así que es necesaria
## su normalización.

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

## Control de calidad. 

plot(qc(microarray.raw.data))

## Valores potencialmente problemáticos para GSM38617 (muestra 
## en anaerobiosis y con sacarosa) por posibles procesos de 
## degradación del RNA.

############### PREPROCESAMIENTO DE LOS DATOS ##################

## Robust Multiarray Average (RMA) para sustraer la fluorescencia
## de fondo, normalizar la fluorescencia de cada sonda y estimar
## los niveles de expresión de los genes representados por las
## distintas sondas en log2:

microarray.processed.data <- rma(microarray.raw.data)

##    Sustracción de la fluorescencia de fondo y normalización
##    de la fluorescencia para cada sonda

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

##    Estimación de los niveles de expresión (en log2)

expression.level <- exprs(microarray.processed.data)
sampleID <- c("Aer_-Suc1","Aer_-Suc2","Ana_-Suc1","Ana_-Suc2","Ana_+Suc1","Ana_+Suc2")
colnames(expression.level) <- sampleID
head(expression.level)

## Para la previsualización de los datos en un scatterplot,
## se calculan valores medios de expresión para condición
## (en este caso se asume que los procesos de degradación han
## tenido lugar y no se considerarán los datos para la primera
## réplica de la condición en anaerobiosis con sacarosa).  
## Los datos se almacenan en una matriz y se realizan los gráficos.

aer.no.suc <- (expression.level[,"Aer_-Suc1"]+expression.level[,"Aer_-Suc2"])/2
ana.no.suc <- (expression.level[,"Ana_-Suc1"]+expression.level[,"Ana_-Suc2"])/2
ana.wi.suc <- expression.level[,"Ana_+Suc2"]

mean.expression<-matrix(data = c(aer.no.suc,ana.no.suc,ana.wi.suc), nrow = nrow(expression.level), ncol = 3)
colnames(mean.expression)<-c("aer.no.suc", "ana.no.suc", "ana.wi.suc")
rownames(mean.expression)<-rownames(expression.level)
head(mean.expression)

## PRIMER CONTRASTE: AEROBIOSIS - ANAEROBIOSIS (sin sacarosa)
## El objetivo de este contraste es conocer los genes reprimidos o
## sobreexpresados en condiciones de anaerobiosis.

plot(main = "Anaerobiosis effects", x = aer.no.suc, y = ana.no.suc,xlab="Aerobiosis - Sucrose",ylab="Anaerobiosis - Sucrose",pch=19,cex=0.5, cex.lab = 1.5)

## SEGUNDO CONTRASTE: ANAEROBIOSIS (con sacarosa) -
##                    ANAEROBIOSIS (sin sacarosa)
## El objetivo de este contraste es conocer los genes reprimidos o
## sobreexpresados cuando se suministra sacarosa a plantas 
## sometidas a anaerobiosis.

plot(main = "Sucrose effects on anaerobiosis", x = ana.no.suc, y = ana.wi.suc,xlab="Anaerobiosis - Sucrose",ylab="Anaerobiosis + Sucrose",pch=19,cex=0.5, cex.lab = 1.5)

## TERCER CONTRASTE: AEROBIOSIS - ANAEROBIOSIS (con sacarosa)
## El objetivo de este contraste es conocer si la adición de 
## sacarosa restaura o alivia el estrés producido por la
## anaerobiosis.

plot(main = "Effects of sucrose on
anaerobic conditions", x = aer.no.suc, y = ana.wi.suc,xlab="Aerobiosis - Sucrose",ylab="Anaerobiosis + Sucrose",pch=19,cex=0.5, cex.lab = 1.5)

######## SELECCIÓN DE GENES DIFERENCIALMENTE EXPRESADOS ########

## Elaboración de la matriz de los niveles de expresión.
## Es probable que GSM38617 haya sufrido graves procesos de 
## degradación por lo que solo se tendrán en cuenta los datos 
## de la segunda réplica en las mismas condiciones (GSM38618).

expression.level <- cbind(expression.level, expression.level[,6])
expression.level <- expression.level[,-5]
sampleID2 <- c("Aer_-Suc1","Aer_-Suc2","Ana_-Suc1","Ana_-Suc2","Ana_+Suc2","Ana_+Suc2_Bis")
colnames(expression.level) <- sampleID2

head(expression.level)

experimental.design <- model.matrix(~ -1+factor(c(1,1,2,2,3,3)))
colnames(experimental.design) <- c("Aer_no_Suc","Ana_no_Suc","Ana_wi_Suc")
linear.fit <- lmFit(expression.level, experimental.design)

contrast.matrix <- makeContrasts(Ana_no_Suc-Aer_no_Suc,
                                 Ana_wi_Suc-Ana_no_Suc, 
                                 Ana_wi_Suc-Aer_no_Suc,
                                 levels=c("Aer_no_Suc","Ana_no_Suc", "Ana_wi_Suc"))

contrast.linear.fit <- contrasts.fit(linear.fit, contrast.matrix)
contrast.results <- eBayes(contrast.linear.fit)




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

## PRIMER CONTRASTE: AEROBIOSIS -SACAROSA / ANAEROBIOSIS -SACAROSA

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


## SEGUNDO CONTRASTE: ANAEROBIOSIS -SACAROSA / ANAEROBIOSIS +SACAROSA

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


## TERCER CONTRASTE: AEROBIOSIS -SACAROSA / ANAEROBIOSIS +SACAROSA

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


## MAPA DE CALOR. Se describe detalladamente su utilidad y los resultados de él en el texto.

library(gplots)
library(RColorBrewer)

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

############### ANOTACIÓN FUNCIONAL Y DISCUSIÓN ####################################################

## Se detallan los procesos seguidos en el texto. Se usó GOrilla
## (http://cbl-gorilla.cs.technion.ac.il/) y ReViGO 
## (http://revigo.irb.hr/) para determinar los genes 
## diferencialmente expresados. Como ejemplo:

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

####################### BIBLIOGRAFIA ################################################################
##
## 1. Gautier, L., Cope, L., Bolstad, B. M., and Irizarry, R. A. 2004. affy---analysis of Affymetrix
##    GeneChip data at the probe level. Bioinformatics 20, 3 (Feb. 2004), 307-315.
## 2. Bolstad, BM (2004) Low Level Analysis of High-density Oligonucleotide Array Data: Background,
##    Normalization and Summarization. Dissertation. University of California, Berkeley.
## 3.  Bolstad BM, Collin F, Brettschneider J, Simpson K, Cope L, Irizarry RA, and Speed TP. (2005)
##    Quality Assessment of Affymetrix GeneChip Data in Bioinformatics and Computational Biology
##    Solutions Using R and Bioconductor. Gentleman R, Carey V, Huber W, Irizarry R, and Dudoit S.
##    (Eds.), Springer, New York.
## 4. Brettschneider J, Collin F, Bolstad BM, and Speed TP. (2007) Quality assessment for short
##    oligonucleotide arrays. Technometrics. In press.
## 5. Crispin J Miller (2019). simpleaffy: Very simple high level analysis of Affymetrix data.
##    http://www.bioconductor.org, http://bioinformatics.picr.man.ac.uk/simpleaffy/.
## 6. Ritchie, M.E., Phipson, B., Wu, D., Hu, Y., Law, C.W., Shi, W., and Smyth, G.K. (2015). limma
##    powers differential expression analyses for RNA-sequencing and microarray studies. Nucleic Acids
##    Research 43(7), e47.
## 7. Colin A. Smith (2019). annaffy: Annotation tools for Affymetrix biological metadata. R package
##    version 1.56.0.
## 8. Carlson M (2016). ath1121501.db: Affymetrix Arabidopsis ATH1 Genome Array annotation data
##    (chip ath1121501). R package version 3.2.3.
## 9. Loreti, E., Poggi, A., Novi, G., Alpi, A., & Perata, P. (2005). A Genome-Wide Analysis of 
##    the Effects of Sucrose on Gene Expression in Arabidopsis Seedlings under Anoxia. Plant 
##    Physiology, 137(3), 1130-1138. doi: 10.1104/pp.104.057299
##
####################################################################################################