################## Applied Bioinformatics - Rcode ##################
#### Nikita Benner, Maximilian Opitz, Zoé Zick ####

###### Preperation before starting the whole Rscript part ######

# Setting up a project-specific library to avoid any problems

Sys.unsetenv("R_LIBS_USER")
dir.create("RLibrary")
.libPaths()
.libPaths(paste(getwd(), "RLibrary", sep="/"))
setRepositories()

# installing necessary packages 

#install.packages('BiocManager')
#install.packages('tidyverse')
#install.packages("tximport")
#install.packages('ensembldb')
#install.packages('rhdf5')
#install.packages('datapasta')
#install.packages("edgeR")
#install.packages('matrixStats')
#install.packages('cowplot')
#install.packages('DT')
#install.packages('plotly')
#install.packages('gt')
#install.packages("reshape2")
#install.packages('heatmaply')
#install.packages("gplots")
#install.packages("GSEABase")
#install.packages("Biobase")
#install.packages("GSVA")
#install.packages("gprofiler2")
#install.packages("clusterProfiler")
#install.packages("msigdbr")
#install.packages("enrichplot")
#install.packages("qusage")
#install.packages("patchwork")

# loading all the installed packages 

library(rhdf5) # Used in Part 1 / used for functions contained to handle hd5 files of kallisto output
library(tidyverse) # Used in all parts / Extensive package containing alot of different useful things 
library(tximport) # Used in Part 1 / package for getting Kallisto results into R
library(ensembldb) # Used in Part 1 / package to handle ensembl annotated data
library(edgeR) # Used in Part 2 and 4 / only used for the DGEList object and for normalization methods
library(matrixStats) # Used in Part 2 / used for easy calculation of stats on rows or columns of a data matrix
library(cowplot) # Used in Part 2 / package used to combine multiple plots into one 
library(DT) # Used in Part 3 - 5 / to create interactive tables
library(plotly) # Used in Part 3 and 4 / to make interactive plots
library(gt) # Used in Part 3 and 4/ ggplot, but for tables
library(limma) # Used in Part 4 and 5 / venerable package for differential gene expression using linear modeling
library(reshape2) # Used in Part 4
library(heatmaply) # Used in Part 4
library(gplots) # Used in Part 5 / to create heatmaps
#library(GSEABase) #USed for Part 5 / functions and methods for Gene Set Enrichment Analysis
#library(Biobase) #Used for Part 5 / base functions for bioconductor; required by GSEABase
#library(GSVA) # Used for Part 5 / Gene Set Variation Analysis, a non-parametric and unsupervised method for estimating variation of gene set enrichment across samples.
library(gprofiler2) #Used for Part 5 / tools for accessing the GO enrichment results using g:Profiler web resources
library(clusterProfiler) # USed for Part 5 / provides a suite of tools for functional enrichment analysis
#library(msigdbr) # Used for Part 5 / access to msigdb collections directly within R
#library(enrichplot) # Used for Part 5 / great for making the standard GSEA enrichment plots
#library(qusage) # Used for Part 5 / Quantitative Set Analysis for Gene Expression


###### Part 1: annotation of the data ######

# read in the studydesign file created for our study 
targets <- read_tsv("studydesign_canna.txt")

# reading in the kallisto abundance files created using the file.path function 
path <- file.path(targets$sra_accession, "abundance.tsv") 

# check if the path if correct as savety measure 
all(file.exists(path)) 

# getting the annotations fitting to cannabis sativa female using BiomaRt
library(biomaRt)

# choosing the plant ensemble as we are dealing with plants - accessing the plant genomes
listMarts(host="https://plants.ensembl.org") 

# defining the plant mart as biomart
myMart <- useMart(biomart="plants_mart", host="https://plants.ensembl.org")

# take the ensemble annotations for cannabis sativa female
csfem.anno <- useMart(biomart="plants_mart", dataset = "csfemale_eg_gene", host="https://plants.ensembl.org")
csfem.attributes <- listAttributes(csfem.anno)

# looking at all the available attributes withing the cannabis sativa female annotation
Tx.csfem <- getBM(attributes=c("ensembl_transcript_id", "ensembl_gene_id", "description"), mart = csfem.anno)

# turining this into a table 
Tx.csfem <- as_tibble(Tx.csfem)

# renaming the two columns fitting to the studydesigne file
Tx.csfem <- dplyr::rename(Tx.csfem, target_id = ensembl_transcript_id, gene_name = ensembl_gene_id)

# and reorder the columns 
Tx.csfem <- dplyr::select(Tx.csfem, "target_id", "gene_name")

# finally import the kallisto transcipt counts into R by copying the adundance files into the wd
Txi_gene <- tximport(path, 
                     type = "kallisto", 
                     tx2gene = Tx.csfem, #Mapping from transcript IDs to gene IDs 
                     txOut = FALSE, #How does the result change if this =FALSE vs =TRUE?
                     countsFromAbundance = "lengthScaledTPM",
                     ignoreTxVersion = TRUE)

Txi_trans <- tximport(path, 
                      type = "kallisto", 
                      tx2gene = Tx.csfem, #Mapping from transcript IDs to gene IDs 
                      txOut = TRUE, #How does the result change if this =FALSE vs =TRUE?
                      countsFromAbundance = "lengthScaledTPM",
                      ignoreTxVersion = TRUE)


###### Part 2: normalization, filtering and visualization of the data ######

# examining the data we created in Part 1

myTPM_gene <- Txi_gene$abundance

mycountsg <- Txi_gene$counts

colSums(myTPM_gene)

colSums(mycountsg)


# re-save the labels for the samples from the "targets"

targets
sampleLables <- targets$sample

# create the summary statistics for the transcript and genes
# add those to the data matrix and transform this data matrix

myTPM_gene.stats <-  transform(myTPM_gene, 
                               SD=rowSds(myTPM_gene), 
                               AVG=rowMeans(myTPM_gene),
                               MED=rowMedians(myTPM_gene))

head(myTPM_gene.stats)


# make differential gene expression list from the counts created above
DEGlist_gene <- DGEList(mycountsg)

cpm <- cpm(DEGlist_gene)

# create the log counts per million and put those into a data frame
log2.cpm_g <- cpm(DEGlist_gene, log=TRUE)


log2.cpm_g.df <- as_tibble(log2.cpm_g, rownames = "geneID")

# add the sampleLables to the dataframes
colnames(log2.cpm_g.df) <- c("geneID", sampleLables)


# create pivot table 
log2.cpm_g.df.pivot <- pivot_longer(log2.cpm_g.df, # dataframe to be pivoted
                                 cols = FL21:TR43, # column names to be stored as a SINGLE variable
                                 names_to = "samples", # name of that new variable (column)
                                 values_to = "expression") # name of new variable (column) storing all the values (data)


# now plot the pivot data tables
gene_p1 <- ggplot(log2.cpm_g.df.pivot) +
  aes(x=samples, y=expression, fill=samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median", 
               geom = "point", 
               shape = 95, 
               size = 10, 
               color = "black", 
               show.legend = FALSE) +
  labs(y="log2 expression", x = "sample",
       title="Log2 Counts per Million (CPM)",
       subtitle="unfiltered, non-normalized",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()

# now the filtering starts 
table(rowSums(DEGlist_gene$counts==0)==12)

keepers <- rowSums(cpm>1)>=6

# now filter according to the said condition above
DEGlist_gene.filtered <- DEGlist_gene[keepers,]
dim(DEGlist_gene.filtered) #1907   12

log2.cpm_g.filtered <- cpm(DEGlist_gene.filtered, log=TRUE)
log2.cpm_g.filtered.df <- as_tibble(log2.cpm_g.filtered, rownames = "geneID")
colnames(log2.cpm_g.filtered.df) <- c("geneID", sampleLables)

# now pivot the filtered data 
# pivot this FILTERED data, just as you did earlier
log2.cpm_g.filtered.df.pivot <- pivot_longer(log2.cpm_g.filtered.df, # dataframe to be pivoted
                                           cols = FL21:TR43, # column names to be stored as a SINGLE variable
                                           names_to = "samples", # name of that new variable (column)
                                           values_to = "expression") # name of new variable (column) storing all the values (data)

# and plot all of this 

gene_p2 <- ggplot(log2.cpm_g.filtered.df.pivot) +
  aes(x=samples, y=expression, fill=samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median", 
               geom = "point", 
               shape = 95, 
               size = 10, 
               color = "black", 
               show.legend = FALSE) +
  labs(y="log2 expression", x = "sample",
       title="Log2 Counts per Million (CPM)",
       subtitle="filtered, non-normalized",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()

# and now we normalize the data whoohoo
DEGlist_gene.filtered.norm <- calcNormFactors(DEGlist_gene.filtered, method = "TMM")

# use cpm to get counts per milion of the normalized data 
log2.cpm_g.filtered.norm <- cpm(DEGlist_gene.filtered.norm, log=TRUE)
log2.cpm_g.filtered.norm.df <- as_tibble(log2.cpm_g.filtered.norm, rownames = "geneID")
colnames(log2.cpm_g.filtered.norm.df) <- c("geneID", sampleLables)

# pivot the normalized data as done before
log2.cpm_g.filtered.norm.df.pivot <- pivot_longer(log2.cpm_g.filtered.norm.df,
                                                cols = FL21:TR43,
                                                names_to = "samples", 
                                                values_to = "expression") 


gene_p3 <- ggplot(log2.cpm_g.filtered.norm.df.pivot) + 
  aes(x=samples, y=expression, fill=samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median", 
               geom = "point", 
               shape = 95, 
               size = 10, 
               color = "black", 
               show.legend = FALSE) +
  labs(y="log2 expression", x = "sample",
       title="Log2 Counts per Million (CPM)",
       subtitle="filtered, TMM normalized",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()


# finally plot all the plots created together in one using the cowplot package 
plot_grid(gene_p1, gene_p2, gene_p3, labels = c('A', 'B', 'C'), label_size = 12)



###### Part 3: data exploration and data wrangling ######

# get the interesting variables from the study design file 
targets
organ <- targets$organ
organ <- factor(organ)

cond <- targets$condition
cond <- factor(cond)


group <- factor(targets$group)
# Hierarchical Clustering will now begin
# need to transform the log2.cpm_g.filtered.norm data frame into matrices

distance <- dist(t(log2.cpm_g.filtered.norm), method = "euclidean") 
clusters <- hclust(distance, method = "complete") 

plot(clusters, labels=sampleLables, main = "Hierarchical cluster of Cannabis Sativa Female Data")
#  weird plot -> Clustering doesnt make much sense so far -> requires interpretation

# And now PCA aka Principle Component Analysis

# calculate the PCA results 
pca.res.g <- prcomp(t(log2.cpm_g.filtered.norm), scale.=F, retx=T)

# captures the eigenvalues of the PCA results 
pc.var.g <- pca.res.g$sdev^2 

#calculate the percentage variance explained by each PC using the eigenvalues 
pc.per.g <- round(pc.var.g/sum(pc.var.g)*100, 1) 

pca.res.g.df <- as_tibble(pca.res.g$x)


pca.plot.g <- ggplot(pca.res.g.df) +
  aes(x=PC1, y=PC2, label=sampleLables, color = organ) +
  geom_point(size=4) +
  stat_ellipse() +
  xlab(paste0("PC1 (",pc.per.g[1],"%",")")) + 
  ylab(paste0("PC2 (",pc.per.g[2],"%",")")) +
  labs(title="PCA plot - organs",
       caption=paste0("produced on ", Sys.time())) +
  coord_fixed() +
  theme_bw()

ggplotly(pca.plot.g)


pca.plot.g2 <- ggplot(pca.res.g.df) +
  aes(x=PC1, y=PC2, label=sampleLables, color = cond) +
  geom_point(size=4) +
  stat_ellipse() +
  xlab(paste0("PC1 (",pc.per.g[1],"%",")")) + 
  ylab(paste0("PC2 (",pc.per.g[2],"%",")")) +
  labs(title="PCA plot - stages",
       caption=paste0("produced on ", Sys.time())) +
  coord_fixed() +
  theme_bw()

ggplotly(pca.plot.g2)



pca.plot.g3 <- ggplot(pca.res.g.df) +
  aes(x=PC1, y=PC2, label=sampleLables, color = group) +
  geom_point(size=4) +
  #stat_ellipse() +
  xlab(paste0("PC1 (",pc.per.g[1],"%",")")) + 
  ylab(paste0("PC2 (",pc.per.g[2],"%",")")) +
  labs(title="PCA plot - organs and stages combined",
       caption=paste0("produced on ", Sys.time())) +
  coord_fixed() +
  theme_bw()

ggplotly(pca.plot.g3)
# use dplyr 'mutate' function to add new columns based on existing data which is given in the studydesign file 

mydata.g.df <- mutate(log2.cpm_g.filtered.norm.df,
                      flower.AVG = (FL21 + FL22 + FL23 + FL41 + FL42 + FL43)/6, 
                      trichomes.AVG = (TR21 + TR22 + TR23 + TR41 + TR42 + TR43)/6,
                    #now make columns comparing each of the averages above that you're interested in
                      LogFC = (trichomes.AVG - flower.AVG)) %>% #note that this is the first time you've seen the 'pipe' operator
  mutate_if(is.numeric, round, 2)

# create an interactive table based on the mutated data

datatable(mydata.g.df[,c(1,12:14)], 
          extensions = c('KeyTable', "FixedHeader"), 
          filter = 'top',
          options = list(keys = TRUE, 
                         searchHighlight = TRUE, 
                         pageLength = 10, 
                         #dom = "Blfrtip", 
                         #buttons = c("copy", "csv", "excel"),
                         lengthMenu = c("10", "25", "50", "100")))

write_tsv(mydata.g.df[,c(1,12:14)],"Mutated_data_table.txt")

###### Part 4: identify differentially expressed genes (DEGs) and differential transcript usage (DTU) ######

##### DEG #####

# setting up the design matrix 
group <- factor(targets$group)
cond <- factor(targets$condition)

design <- model.matrix(~0 + group) # design <- model.matrix(~block + treatment)
#colnames(design) <- levels(group)

# model the mean-variance relationship using VOOM from the limma package 
v.DEGList.filtered.norm.g <- voom(DEGlist_gene.filtered.norm, design, plot = FALSE)

# fit a linear model to the data 
fit_g <- lmFit(v.DEGList.filtered.norm.g, design)

# create a contrast matrix  
contrast.matrix <- makeContrasts(organs = (grouptrichome_S2+grouptrichome_S4)/2 - (groupflower_S2+groupflower_S4)/2,
                                 stages = (grouptrichome_S2+groupflower_S2)/2 - (grouptrichome_S4+groupflower_S4)/2,
                                 levels=colnames(design))

# extract the linear model fit for each 
fits_g <- contrasts.fit(fit_g, contrast.matrix)

# calculate and save the bayesian stats for each fit
ebFit_g <- eBayes(fits_g)
write.fit(ebFit_g, file="lmfit_results_gene.txt")


# create a toptable to view the DEGs and convert those of cause into tibbles 

myTopHits_g <- topTable(ebFit_g, adjust ="BH", coef=1, number=40000, sort.by="logFC")
myTopHits.df.g <- myTopHits_g %>%
  as_tibble(rownames = "geneID")

gt(myTopHits.df.g)

# finally create another interactive (or not interactive) plot using the toptable data creaated above 
vplot_g <- ggplot(myTopHits_g) +
  aes(y=-log10(adj.P.Val), x=logFC, text = paste("Symbol:", myTopHits.df.g$geneID)) +
  geom_point(size=2) +
  geom_hline(yintercept = -log10(0.05), linetype="longdash", colour="grey", size=1) +
  geom_vline(xintercept = 1, linetype="longdash", colour="#6633CC", size=1) +
  geom_vline(xintercept = -1, linetype="longdash", colour="#FF00CC", size=1) +
  #annotate("rect", xmin = 1, xmax = 12, ymin = -log10(0.01), ymax = 7.5, alpha=.2, fill="#BE684D") +
  #annotate("rect", xmin = -1, xmax = -12, ymin = -log10(0.01), ymax = 7.5, alpha=.2, fill="#2C467A") +
  labs(title="Volcano plot - Cannabis Sativa Female",
       subtitle = "Cannabis Sativa Female",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()

ggplotly(vplot_g)


# now pull out the DEGs
results_g <- decideTests(ebFit_g, method="global", adjust.method="BH", p.value=0.05, lfc=1)
head(results_g)
summary(results_g)


# now retrieve the expression data for the DEGs 
colnames(v.DEGList.filtered.norm.g$E) <- sampleLables


diffGenes_g_org <- v.DEGList.filtered.norm.g$E[results_g[,1] !=0,]
diffGenes_g_stg <- v.DEGList.filtered.norm.g$E[results_g[,2] !=0,]

# as always convert them into tibbles 
diffGenes.df.g.org <- as_tibble(diffGenes_g_org, rownames = "geneID")
diffGenes.df.g.stg <- as_tibble(diffGenes_g_stg, rownames = "geneID")

# create an interactive table for the DEGs and create a file containing them 
datatable(diffGenes.df.g.org,
          extensions = c('KeyTable', "FixedHeader"),
          caption = 'Table 1: DEGs in Cannabis Sativa Female',
          options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 10, lengthMenu = c("10", "25", "50", "100"))) %>%
  formatRound(columns=c(2:11), digits=2)


datatable(diffGenes.df.g.stg,
          extensions = c('KeyTable', "FixedHeader"),
          caption = 'Table 1: DEGs in Cannabis Sativa Female',
          options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 10, lengthMenu = c("10", "25", "50", "100"))) %>%
  formatRound(columns=c(2:11), digits=2)

write_tsv(diffGenes.df.g.org,"DiffGenes_g_org.txt") 
write_tsv(diffGenes.df.g.stg,"DiffGenes_g_stg.txt")



###### Part 5: functional enrichment ######

##### GO #####

# starting with the GO enrichment using the gProfiler2
# chose first the top genes using toptable for the GO 
myTopHits_g <- topTable(ebFit_g, adjust ="BH", coef=1, number=50, sort.by="logFC")

# using the gost function of gProfiler2 to create the GO functional enrichment using our data 
gost.res.g <- gost(rownames(myTopHits_g), organism = "csfemale", correction_method = "fdr", significant = F)

# create an interactive (or not) plot using the gost data 
gostplot(gost.res.g, interactive = T, capped = F)
gene_gost_plot <- gostplot(gost.res.g, interactive = F, capped = F)

# save the non-interactive plot and a table of it 
publish_gostplot(
  gene_gost_plot, 
  highlight_terms = c(),
  filename = NULL,
  width = NA,
  height = NA)


# publish_gosttable(
#    gost.res.g,
#    highlight_terms = NULL,
#    use_colors = TRUE,
#    show_columns = c("source", "term_name", "term_size", "intersection_size"),
#    filename = NULL,
#    ggplot=TRUE)
# 
# 
# write_tsv(gost.res.g,"Gost_table.txt")
