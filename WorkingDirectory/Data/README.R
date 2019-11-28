################################################################
#                                                              #
#                            TAREA 1                           #
#                MICROARRAYS - PAQUETES REQUERIDOS             #
#                                                              #
#    Biología Molecular de Sistemas - Grado en Bioquímica      #    
#     Joaquín Tamargo Azplicueta - joatamazp@alum.us.es        #
#                                                              #
################################################################

## Con el fin de facilitar la revisión de nuestros resultados, se han recopilado los paquetes 
## necesarios en este anexo. Estos paquetes son: affy (1.62.0), affyPLM (1.60.0),
## simpleaffy (2.60.0), limma (3.40.6), annaffy (1.56.0) y las anotaciones del chip ATH1 
## de Affymetrix (ath1121501.db).

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("affy")
BiocManager::install("affyPLM")
BiocManager::install("simpleaffy")
BiocManager::install("limma")
BiocManager::install("annaffy")
BiocManager::install("ath1121501.db")

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
##
####################################################################################################