# Microarray Analysis in RStudio - Straightforward instructions by Affy

##### Joaquín Tamargo Azpilicueta (joatamazp@alum.us.es), Universidad de Sevilla. 28th November, 2019

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

Reviewed article:

**Loreti, E., Poggi, A., Novi, G., Alpi, A., & Perata, P. (2005). A Genome-Wide Analysis of the Effects of Sucrose on Gene Expression in Arabidopsis Seedlings under Anoxia. *Plant Physiology*, 137(3), 1130-1138. doi: 10.1104/pp.104.057299**

Aditionally:

1. Mickelbart, M., Hasegawa, P., & Bailey-Serres, J. (2015). Genetic mechanisms of abiotic stress tolerance that translate to crop yield stability. *Nature Reviews Genetics*, 16(4), 237-251. doi: 10.1038/nrg3901

2. Kumar, D., Hazra, S., Datta, R., & Chattopadhyay, S. (2016). Transcriptome analysis of Arabidopsis mutants suggests a crosstalk between ABA, ethylene and GSH against combined cold and osmotic stress. *Scientific Reports*, 6(1). doi: 10.1038/srep36867

3. Pérez-Mejías, G., Guerra-Castellano, A., Díaz-Quintana, A., De la Rosa, M., & Díaz-Moreno, I. (2019). Cytochrome c: Surfing Off of the Mitochondrial Membrane on the Tops of Complexes III and IV. *Computational And Structural Biotechnology Journal*, 17, 654-660. doi: 10.1016/j.csbj.2019.05.002

4. Guerra-Castellano, A., Díaz-Quintana, A., Pérez-Mejías, G., Elena-Real, C., González-Arzola, K., & García-Mauriño, S. et al. (2018). Oxidative stress is tightly regulated by cytochrome c phosphorylation and respirasome factors in mitochondria. *Proceedings Of The National Academy Of Sciences*, 115(31), 7955-7960. doi: 10.1073/pnas.1806833115

5. Subramanian, A., Tamayo, P., Mootha, V., Mukherjee, S., Ebert, B., & Gillette, M. et al. (2005). Gene set enrichment analysis: A knowledge-based approach for interpreting genome-wide expression profiles. Proceedings Of The National Academy Of Sciences, 102(43), 15545-15550. doi: 10.1073/pnas.0506580102

**Molecular Systems Biology [Biología Molecular de Sistemas, 3rd year Biochemistry Degree, Universidad de Sevilla] slides, made by Francisco J. Romero Campero & Ignacio Pérez Hurtado de Mendoza, CC BY-NC-ND 3.0**

For the bioinformatic analysis, these tools and packages have been used:

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

**There are two R scripts in Data folder that may be useful: `README.r` contains the download instructions for the packages used and `script.r`contains the R raw script.**
