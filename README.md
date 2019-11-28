# Microarray Analysis in RStudio - Straightforward instructions by Affy

In this repository, I thoroughly analyze the microarray data from this article: A Genome-Wide Analysis of the Effects of Sucrose on Gene Expression in Arabidopsis Seedlings under Anoxia. *Plant Physiology*, 137(3), 1130-1138. doi: 10.1104/pp.104.057299. 

Plants are constantly exposed to non-optimal edafic and climatic conditions that affect cell homeostasis and, ultimately, impair their growth and integrity. In order to survive, plants must adapt in the most effective way to this circumstances; i.e. specific differentiation or modulation of phisiologycal behaviour in response to stress signaling. Earlier, changes in genes expression model in some stress conditions had been identified such as anaerobic conditions, that rapidly trigger an adaptative response to hipoxia. These symptoms seem to be attenuated when exogen sugar is added.

Lately, whole-transcriptomic analysis such as microarrays have made it possible to understand more profoundly which genes (and therefore, supposedly proteins) are implicated in that response. In the analysed article, Affymetrix ATH1 microarray chip is used to determine the whole-genome expression effect of hypoxia and supply of exogen sugar in this circumstance. Furthermore, I would like to underline that auxiliary experiments were carried out in the article in order to biologically prove the results provided by the transcriptomic experiment.

The whole analysis have been done in R. To make it easy for the reader to corroborate the results in my script, packages have been compiled under. These packeges are: affy (1.62.0), affyPLM (1.60.0), simpleaffy (2.60.0), limma (3.40.6), annaffy (1.56.0) and Affymetrix ATH1 microarray annotation chip (ath1121501.db).

``` {r}
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("affy")
BiocManager::install("affyPLM")
BiocManager::install("simpleaffy")
BiocManager::install("limma")
BiocManager::install("annaffy")
BiocManager::install("ath1121501.db"

```
-----

### Bibliography
1. Gautier, L., Cope, L., Bolstad, B. M., and Irizarry, R. A. 2004. affy---analysis of Affymetrix GeneChip data at the probe level. Bioinformatics 20, 3 (Feb. 2004), 307-315.
2. Bolstad, BM (2004) Low Level Analysis of High-density Oligonucleotide Array Data: Background, Normalization and Summarization. Dissertation. University of California, Berkeley.
3.  Bolstad BM, Collin F, Brettschneider J, Simpson K, Cope L, Irizarry RA, and Speed TP. (2005) Quality Assessment of Affymetrix GeneChip Data in Bioinformatics and Computational Biology Solutions Using R and Bioconductor. Gentleman R, Carey V, Huber W, Irizarry R, and Dudoit S. (Eds.), Springer, New York.
4. Brettschneider J, Collin F, Bolstad BM, and Speed TP. (2007) Quality assessment for short oligonucleotide arrays. Technometrics. In press.
5. Crispin J Miller (2019). simpleaffy: Very simple high level analysis of Affymetrix data. http://www.bioconductor.org, http://bioinformatics.picr.man.ac.uk/simpleaffy/.
6. Ritchie, M.E., Phipson, B., Wu, D., Hu, Y., Law, C.W., Shi, W., and Smyth, G.K. (2015). limma powers differential expression analyses for RNA-sequencing and microarray studies. Nucleic Acids Research 43(7), e47.
7. Colin A. Smith (2019). annaffy: Annotation tools for Affymetrix biological metadata. R package version 1.56.0.
8. Carlson M (2016). ath1121501.db: Affymetrix Arabidopsis ATH1 Genome Array annotation data (chip ath1121501). R package version 3.2.3.
