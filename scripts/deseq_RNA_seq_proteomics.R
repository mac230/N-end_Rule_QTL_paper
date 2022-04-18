## -----
## load packages and set up

## check for Bioconductor and install if not available
ifelse(!requireNamespace("BiocManager", quietly = TRUE),
       install.packages("BiocManager",
                        dependencies = TRUE,
                        repos = "http://cran.wustl.edu/",
                        quiet = TRUE),
       paste0("Bioconductor available"))
require("BiocManager")

## -----
## load packages or install if not available
## have to split these out by bioconductor vs. non-bioconductor
## non-bioconductor
package_installer <- function(x){
    if(!requireNamespace(x, quietly = TRUE))
        install.packages(x, dependencies = TRUE,
                         repos = "http://cran.wustl.edu/",
                         quiet = TRUE, INSTALL_opts = '--no-lock')}

required_packages <- c("tidyverse", "DESeq2", "readxl", "lattice",
                       "readr", "qvalue", "sva", "tximport")

sapply(X = required_packages, FUN = package_installer)
sapply(X = required_packages, FUN = require, character.only = TRUE)


## bioconductor
bioc_package_installer <- function(x){if(!requireNamespace(x))
                                          BiocManager::install(x, INSTALL_opts = '--no-lock')}
bioc_packages <-  c("qvalue", "sva", "DESeq2", "tximportData",
                    "BUSpaRse", "biomaRt", "org.Sc.sgd.db")
sapply(X = bioc_packages, FUN = bioc_package_installer)
sapply(X = bioc_packages, FUN = require, character.only = TRUE)
## BiocManager::install("BUSpaRse")
## source("~/emacs/N-end_Rule_QTL_paper/UBR1_RNA-seq_analysis/deseq/get_TX.R")
## source("~/emacs/N-end_Rule_QTL_paper/UBR1_RNA-seq_analysis/deseq/utils.R")
## install.packages("remotes")
## remotes::install_github("lambdamoses/BUStoolsR")


## -----
## directories
base_dir   <- "~/emacs/N-end_Rule_QTL_paper/UBR1_RNA-seq_analysis/fastp/"

## -----
## sample information
## read in a list of samples and covariates
samples <- read.csv(paste0(base_dir, "ubr1_RNA-seq_sample_info.csv"),
                    header = T, sep = ",")

## make sure "wild-type" is the reference level for our factor
samples$Group <- relevel(samples$Group, ref = "wild_type")

## make sure "Batch" is a factor
samples$Batch <- as.factor(samples$Batch)

## scale the concentration variable
samples$Scale_Conc <- scale(samples$Concentration, center = T)

## make a list of abundance files for reading in using 'paste0'
files <- paste0(base_dir, "deseq_input/", samples$Sample, "_abundance.tsv")
names(files) <- as.character(samples$Sample)
## make sure you got the paths/names correct 
sapply(X = files, FUN = file.exists)


## -----
## annotations and other info
gene_annot <- read.table(file = paste0(base_dir, "ensemblGenes_ensembl83_160307_MOD.txt"),
                         sep = "\t", stringsAsFactors = F,
                         head = T)

rownames(gene_annot) <- gene_annot[, 1]
all_names <- gene_annot[, 1] 
names(all_names) <- gene_annot[, 1]
all_names[which(all_names == "")] <- names(all_names)[which(all_names == "")]


## -----
## read in the counts; creates a sample (columns) x gene (rows) matrix 
gene_counts <- sapply(X = files, FUN = function(x) {
                          curr_data <- read.table(x, sep = "\t", head = T)
                          curr_data[, "est_counts"]
                      })

## add in transcript names (these all have the same order in ea. abundance file)
unique(read.table(files[1], sep="\t", head=TRUE, stringsAsFactors = F)$target_id ==
    read.table(files[10], sep="\t", head=TRUE, stringsAsFactors = F)$target_id)

## add row names to transcript abundances
raw_gene_names <- read.table(files[1], sep="\t", head=TRUE, stringsAsFactors = F)$target_id
final_gene_names <- sapply(X = 1:length(raw_gene_names),
                           FUN = function(n){
                               strsplit(x = raw_gene_names[n], split = "_")[[1]][1] 
})
rownames(gene_counts) <- final_gene_names

head(gene_counts)


## -----
## "transcript integrity number" (TIN)-based filtering
## read TIN's for ea. gene, convert to a matrix that contains only the TIN score
## for each sample and drops the other columns of the associated "abundance.tsv" file
tins <- sapply(X = samples$Sample,
               FUN = function(x) {
                   curr_tin <- read.table(paste0(base_dir, "deseq_input/", x,
                                                 "_pseudoalignments.tin.xls"),
                                          head = T, stringsAsFactors = F)
                   rownames(curr_tin) <- curr_tin[, "geneID"]
                   curr_tin <- curr_tin[rownames(gene_counts), ]
                   curr_tin[, "TIN"]
})

## give the TIN scores the same gene names as our counts file
rownames(tins) <- final_gene_names

nrow(tins) ; nrow(gene_counts)

mean_tin <- rowMeans(tins)


## create a vector of effective gene length for filtering 
eff_gene_length <- read.table(paste0(base_dir, "deseq_input/",
                                     samples$Sample[1], "_abundance.tsv"),
                              sep = "\t", header = T)[, "eff_length"]
names(eff_gene_length) <- final_gene_names


## how many genes will each step of our filter below remove?
nrow(tins) ; length(mean_tin[mean_tin > 1])            ## ~150
length(eff_gene_length) ; length(eff_gene_length > 0)  ## 0
length(gene_counts[rowMin(round(gene_counts)) > 0, 1]) ## ~250

genes_filtered <- gene_counts[rowMin(round(gene_counts)) > 0
                              & apply(tins, 1, min) > 0
                              & eff_gene_length > 0, ]

c_filter <- rowSums(genes_filtered) > 1000
genes_filtered <- genes_filtered[c_filter, ]


## remove any rows w/ NAs
genes_filtered <- genes_filtered[complete.cases(genes_filtered), ]

nrow(gene_counts) ; nrow(genes_filtered)
## down to 6206 from 6612

write.table(x = genes_filtered,
            file = paste0(base_dir, "/gene_counts_filtered.csv"),
            append = F, quote = F, sep = ",",
            row.names = T, col.names = T)

col_data <- data.frame(samples[, c("OD", "Batch", "RIN", "Scale_Conc",
                                   "Concentration", "Group", "S_Area"), ])

## make sure column names match the order specified in 'col_data'
data.frame(x = colnames(genes_filtered),
           y = samples$Sample,
           w = col_data$Group)

## build the deseq model
## note: "Results tables are generated using the function results, which
## extracts a results table with log2 fold changes, p values and adjusted p
## values. With no additional arguments to results, the log2 fold change and
## Wald test p value will be for the last variable in the design formula, and if
## this is a factor, the comparison will be the last level of this variable over
## the reference level (see previous note on factor levels). However, the order
## of the variables of the design do not matter so long as the user specifies
## the comparison to build a results table for, using the name or contrast
## arguments of results." from -
## https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
dds <- DESeqDataSetFromMatrix(countData = round(genes_filtered),
                              colData = col_data,
                              design = ~ OD + Batch + RIN + Group)

dds <- DESeq(dds, betaPrior = T)
res <- results(dds)
summary(res)

dat <- counts(dds, normalized = T)
idx <- rowMeans(dat) > 1
dat <- dat[idx, ]

mod  <- model.matrix(~ OD + Batch + RIN + Group,
                     colData(dds))
mod0 <- model.matrix(~ OD + Batch + RIN,
                     colData(dds))
svseq <- svaseq(dat, mod, mod0, n.sv = 2)

col_data_sva <- col_data
col_data_sva$SV1 <- svseq$sv[, 1]
col_data_sva$SV2 <- svseq$sv[, 2]


## -----
## final model that incorporates the surrogate vars 
dds_sva <- DESeqDataSetFromMatrix(countData = round(genes_filtered),
                                  colData = col_data_sva,
                                  design = ~ OD + Batch + RIN + SV1 + SV2 + Group)

dds_sva <- DESeq(dds_sva, betaPrior = T)
res_sva <- results(dds_sva) ; summary(res_sva)

vsd <- vst(dds_sva, blind = F)
plotPCA(vsd, intgroup = "Group")
plotPCA(vsd, intgroup = c("Group", "RNA_Batch"))
plotPCA(vsd, intgroup = c("Group", "Harvest_Batch"))


## -----
## work with the model output: 
## function to return short name for a vector
## of systematic gene names.  'n' below is the
## max value in a logical vector resulting from
## a test of whether a systematic name from the
## result list matches a name in the complete gene
## list - return the short name if there is one;
## if not, return the empty string. 
gene_checker <- function(x) {
    n <- which.max(x == gene_annot$geneID)
    gene_annot[n, "geneName"]
}

res_sva        <- as.data.frame(res_sva)

## add variables for plotting
res_sva$geneID <- rownames(res_sva)
res_sva$log2_fold_change <- log2(1 / 2^(res_sva$log2FoldChange))
res_sva$log10_pvalue     <- -1 * log10(res_sva$pvalue)
res_sva$log10_padj       <- -1 * log10(res_sva$padj)

res_sva$gene_name <- sapply(X = res_sva$geneID,
                            FUN = gene_checker)

## use 'res_out' for overplotting significant genes
res_out <- res_sva[res_sva$padj <= 0.1, ]
res_out$geneID <- rownames(res_out)

write.table(x = res_out,
            file = paste0(base_dir, "DE_genes_RNA_seq.csv"),
            append = F, quote = F, sep = ",",
            row.names = T, col.names = T)

write.table(x = res_out,
            file = paste0("~/emacs/N-end_Rule_QTL_paper/tables_drafts/", "RNA_seq_results.csv"),
            append = F, quote = F, sep = ",",
            row.names = T, col.names = T)


## -----
## plot RNA-seq results
text_genes <- res_out[res_out$log10_padj > 3.46, ]
text_genes[1 + nrow(text_genes), ] <- res_out["YGR184C", ]

pdf(file = paste0("~/emacs/N-end_Rule_QTL_paper/figures_drafts/",
                  gsub(pattern = " .*",
                       replacement = "",
                       x = Sys.time()),
                  "_DE_genes_RNA_volcano_final.pdf"))
print(
xyplot(log10_padj ~ log2_fold_change,
       data = res_sva,
       type = c("g", "p"),
       xlim = c(-0.8, 0.8),
       ylab = expression("Log"["10"]*" RNA "*italic("p")*" Value (BY vs. BY "*italic("UBR1")* "-469A>T)"),
       xlab = expression("Log"["2"]*" RNA Fold Change (BY vs. BY "*italic("UBR1")* "-469A>T)"),
       pch = 21,
       ## equiv. to 'bg' in base R 
       fill = gray(0.9, alpha = 0.5),
       col = gray(0.2, alpha = 0.5),
       cex = 1.2,
       scales = list(x = list(at = seq(from = -0.75, to = 0.75, by = 0.25)),
                     y = list(at = seq(from = 0, to = 25, by = 5)),
                     tck = c(1, 0)),
       key = list(corner = c(0.9975, 0.99),
                  points = list(pch = 21,
                                col = gray(0.2, alpha = 0.5),
                                fill = c(gray(0.7, alpha = 0.5),
                                         "#DDCC77"),
                                cex = 2.5),
                  text = list(labels = c(expression("Corrected "*italic("p")*" > 0.1"),
                                         expression("Corrected "*italic("p")*" < 0.1"))),
                  padding.text = 2.5,
                  between = 0.9,
                  background = gray(1.0, alpha = 0.75)),
       par.settings = list(axis.text = list(cex = 1.2),
                           par.ylab.text = list(cex = 1.25),
                           par.xlab.text = list(cex = 1.25)),
       panel = function(...){
           panel.abline(h = 0, v = 0, lwd = 3, col = gray(0.4))
           panel.abline(h = seq(from = 0, to = 25, by = 5),
                        v = seq(from = -0.75, to = 0.75, by = 0.25),
                        col = gray(0.9),
                        lwd = 1)
           panel.xyplot(...)
           panel.abline(h = 1, col = gray(0.4), lty = 2, lwd = 3)
           panel.points(x = res_out$log2_fold_change,
                        y = res_out$log10_padj,
                        pch = 21,
                        fill = "#DDCC77",
                        col = gray(0.4, alpha = 0.5),
                        lwd = 1.2,
                        cex = 1.75)
           panel.text(ifelse(text_genes$gene_name != "",
                             text_genes$gene_name,
                             text_genes$geneID),
                      x = text_genes$log2_fold_change + 0.07,
                      y = text_genes$log10_padj)
})
)
dev.off()


## -----
## load in proteomics data
res_sva <- read.table(file = paste0(base_dir, "DE_genes_RNA_seq.csv"),
                      sep = ",", stringsAsFactors = F)

prot_dir  <- "~/emacs/N-end_Rule_QTL_paper/UBR1_RNA-seq_analysis/"
prot_file <- "formatted_proteomics_data.csv"
prot_dat  <- read.csv(paste0(prot_dir, prot_file),
                      sep = ",", stringsAsFactors = F,
                      header = T)

## need to format the proteomics data for downstream analysis
## first, get systematic gene names from the protein names field:
x <- prot_dat$Gene.ID
y <- gsub(pattern = ".*sce:", replacement = "", x = x)
length(y)                ## 3783
length(y[nchar(y) < 10]) ## 3775, so 8 w/ no yeast gene name

## what are the 8?
## A5PJX3 = albumin (bovine)
## Q06891 = raf1 (yeast); R0030W (2 micron plasmid gene)
## B2RA01 = keratin (human, i.e., contamination)
## P03871 = rep1 (yeast); R0020C (2 micron plasmid gene)
## P13645 = keratin (human)
## P35908 = keratin (human, type II)
## O00109 = keratin (human, cuticular)
## Q9NSB4 = keratin (human)

prot_dat$geneID <- gsub(pattern = ".*sce:",
                        replacement = "",
                        x = prot_dat$Gene.ID)

prot_name_filter <- nchar(prot_dat$geneID) < 10
length(prot_name_filter[prot_name_filter == F])

prot_dat <- prot_dat[prot_name_filter, ]

## convert the log2 fold change
prot_all$log2_fold_change <- log2(1 / prot_all$fold_change)
prot_all$log10_padj <- prot_all$log10_padj

## take the subset of the proteomics needed and match the columns from the RNA data
prot_all <- prot_dat[, c("Gene.Symbol", "fold_change", "padj", "geneID")]
colnames(prot_all)        <- c("gene_name", "fold_change", "padj", "geneID")
prot_all$log2_fold_change <- log2(1 / prot_all$fold_change)
prot_all$log10_padj       <- -1 * log10(prot_all$padj)
rownames(prot_all)        <- prot_all$geneID


## several hundred rows have no 'padj' value,
## so need the two-step filter here
prot_filter   <- prot_all$padj <= 0.1 & !is.na(prot_all$padj)
prot_filtered <- prot_all[prot_filter, ]
nrow(prot_filtered) ## 39, as expected

write.table(x = prot_filtered,
            file = paste0(base_dir, "DE_genes_proteomcis.csv"),
            append = F, quote = F, sep = ",",
            row.names = T, col.names = T)


## -----
## plot proteomics results
text_proteins <- prot_filtered[prot_filtered$log10_padj > 1.5, ]
text_proteins[1 + nrow(text_proteins), ] <- prot_filtered["YGR184C", ]

pdf(file = paste0("~/emacs/N-end_Rule_QTL_paper/figures_drafts/",
                  gsub(pattern = " .*",
                       replacement = "",
                       x = Sys.time()),
                  "_DE_genes_protein_volcano_final.pdf"))
print(
xyplot(log10_padj_protein_protein ~ log2_fold_change_protein_protein,
       data = prot_all,
       type = c("p"),
       xlim = c(-1.25, 1.1),
       ylim = c(-0.5, 9), 
       ylab = expression("Log"["10"]*" Protein "*italic("p")*" Value (BY vs. BY "*italic("UBR1")* "-469A>T)"),
       xlab = expression("Log"["2"]*" Protein Fold Change (BY vs. BY "*italic("UBR1")* "-469A>T)"),
       pch = 21,
       ## equiv. to 'bg' in base R 
       fill = gray(0.9, alpha = 0.5),
       col = gray(0.2, alpha = 0.5),
       cex = 1.2,
       scales = list(x = list(at = seq(from = -1.2, to = 0.8, by = 0.4),
                              labels = as.character(c(-1.7, seq(from = -0.8,
                                                          to = 0.8,
                                                          by = 0.4)))),
                     y = list(at = seq(from = 0, to = 14, by = 2),
                              labels = as.character(c(seq(from = 0,
                                                          to = 6,
                                                          by = 2), 13.5))),
                     tck = c(1, 0)),
       key = list(corner = c(0.9975, 0.99),
                  points = list(pch = 21,
                                col = gray(0.2, alpha = 0.5),
                                fill = c(gray(0.7, alpha = 0.5),
                                         "#6699CC"),
                                cex = 2.5),
                  text = list(labels = c(expression("Corrected "*italic("p")*" > 0.1"),
                                         expression("Corrected "*italic("p")*" < 0.1"))),
                  padding.text = 2.5,
                  between = 0.9,
                  background = gray(1.0, alpha = 0.75)),
       par.settings = list(axis.text = list(cex = 1.2),
                           par.ylab.text = list(cex = 1.25),
                           par.xlab.text = list(cex = 1.25)),
       panel = function(...){
           panel.abline(h = 0, v = 0, lwd = 3, col = gray(0.4))
           panel.abline(v = seq(from = -1.2, to = 1.6, by = 0.4),
                        h = seq(from = 0, to = 14, by = 2), 
                        col = gray(0.9),
                        lwd = 1)
           panel.xyplot(...)
           panel.abline(h = 1, col = gray(0.4), lty = 2, lwd = 3)
           panel.points(x = prot_filtered$log2_fold_change_protein,
                        y = prot_filtered$log10_padj_protein,
                        pch = 21,
                        fill = "#6699CC",
                        col = gray(0.2, alpha = 0.5),
                        lwd = 1.1,
                        cex = 1.75)
           panel.text(ifelse(text_proteins$gene_name != "",
                             text_proteins$gene_name,
                             text_proteins$geneID),
                      x = text_proteins$log2_fold_change_protein,
                      y = text_proteins$log10_padj_protein + 0.2)
           ## off-axis vps36 plotting
           panel.points(x = -1.2, y = 8,
                        pch = 21,
                        fill = "#6699CC",
                        col = gray(0.2, alpha = 0.5),
                        lwd = 1.1,
                        cex = 1.75)
           panel.text("VPS36", x = -1.1, y = 8.3)
})
)
dev.off()


## format the RNA data the same way
str(res_sva)
RNA_all <- res_sva[, c("gene_name", "geneID",
                       "log2_fold_change", "padj",
                       "log10_padj")]
RNA_all$geneID <- rownames(res_sva)

## format these for later use in plotting
colnames(prot_all) <- paste0(colnames(prot_all), "_protein")
colnames(prot_all) <- gsub(pattern = "_protein_protein",
                           replacement = "_protein",
                           x = colnames(prot_all))
colnames(RNA_all)  <- paste0(colnames(RNA_all), "_RNA")
str(RNA_all) ; str(prot_all)
prot_all

## -----
## now, get the set of genes detected
## at both the protein and RNA levels 
both_all_only <- list()
for(i in 1:nrow(prot_all)) {
    x <- rownames(prot_all)[i]
    y <- x == rownames(RNA_all)
    if (max(y) > 0) {
        both_all_only[[i]] <- cbind(prot_all[i, ], RNA_all[which.max(y), ])
}
}

both_all_only <- do.call("rbind", both_all_only)
nrow(both_all_only) ; nrow(prot_all)

## now, get the genes that are DE at protein and RNA levels
both_de <- list()
for(i in 1:nrow(both_all_only)) {
    x <- both_all_only[i, ]
    y <- x$padj_protein <= 0.1 & x$padj_RNA <= 0.1 & !is.na(x$padj_protein)
    if (y > 0) {
        both_de[[i]] <- both_all_only[i, ]
}
}

both_de <- do.call("rbind", both_de)
## TMA10, HSP26, UBR1

## DE only at the RNA level
rna_de <- list()
for(i in 1:length(rownames(both_all_only))) {
    x <- both_all_only[i, ]
    y <- x$padj_protein > 0.1 & x$padj_RNA <= 0.1
    if (!is.na(x$padj_protein) & y > 0) {
        rna_de[[i]] <- both_all_only[i, ]
}
}

rna_de <- do.call("rbind", rna_de)

## DE only at the protein level
protein_de <- list()
for(i in 1:length(rownames(both_all_only))) {
    x <- both_all_only[i, ]
    y <- x$padj_protein <= 0.1 & x$padj_RNA > 0.1
    if (!is.na(x$padj_protein) & y > 0) {
        protein_de[[i]] <- both_all_only[i, ]
}
}

protein_de <- do.call("rbind", protein_de)

## -----
## plot protein vs. RNA abundance plot 
pdf(file = paste0("~/emacs/N-end_Rule_QTL_paper/figures_drafts/",
                  gsub(pattern = " .*",
                       replacement = "",
                       x = Sys.time()),
                  "_protein_RNA_scatter_final.pdf"))
print(
xyplot(log2_fold_change_protein ~ log2_fold_change_RNA,
       data = both_all_only,
       type = c("p"),
       pch = 21,
       ylab = expression("Log"["2"]*" Fold Change Protein (BY vs. BY "*italic("UBR1")* "-469A>T)"),
       xlab = expression("Log"["2"]*" Fold Change RNA (BY vs. BY "*italic("UBR1")* "-469A>T)"),
       fill = gray(0.9, alpha = 0.5),
       col = gray(0.2, alpha = 0.5),
       cex = 1.1,
       xlim = c(-0.45, 0.45),
       ylim = c(-1, 1),
       scales = list(x = list(at = seq(from = -1, to = 1, by = 0.2)),
                     y = list(at = seq(from = -1, to = 1, by = 0.2)),
                     tck = c(1, 0)),
       key = list(corner = c(0.0001, 0.998),
                  points = list(pch = 21,
                                col = gray(0.2, alpha = 0.5),
                                fill = c("#882255",
                                         "#DDCC77",
                                         "#6699CC",
                                         gray(0.7, alpha = 0.5)),
                                cex = 2.4),
                  text = list(labels = c(expression("Corrected "*italic("p")*" RNA and Protein < 0.1"),
                                         expression("Corrected "*italic("p")*" RNA only < 0.1"),
                                         expression("Corrected "*italic("p")*" Protein only < 0.1"),
                                         expression("Corrected "*italic("p")*" RNA and Protein > 0.1"))),
                  padding.text = 1.5,
                  between = 0.2),
                  ## background = gray(1.0, alpha = 0.75)),
       par.settings = list(axis.text = list(cex = 1.2),
                           par.ylab.text = list(cex = 1.25),
                           par.xlab.text = list(cex = 1.25)),
       panel = function(...) {
           panel.abline(h = 0, v = 0, lwd = 2, col = gray(0.4))
           panel.abline(h = 0, v = 0, lwd = 1, col = gray(0.9))
           panel.abline(h = seq(from = -1, to = 1, by = 0.2),
                        v = seq(from = -1, to = 1, by = 0.2),
                        col = gray(0.9))
           panel.xyplot(...)
           ## regression line for protein vs. mRNA 
           ## panel.abline(lm(both_all_only$log2_fold_change_protein
           ##                 ~ both_all_only$log2_fold_change_RNA),
           ##              lty = 2, col = gray(0.2))
           panel.points(x = both_de[, "log2_fold_change_RNA"],
                        y = both_de[, "log2_fold_change_protein"],
                        pch = 21, fill = "#882255",
                        col = gray(0.2, alpha = 0.5), cex = 2)
           panel.points(x = rna_de[, "log2_fold_change_RNA"],
                        y = rna_de[, "log2_fold_change_protein"],
                        pch = 21, fill = "#DDCC77",
                        col = gray(0.2, alpha = 0.5), cex = 2)                        
           panel.points(x = protein_de[, "log2_fold_change_RNA"],
                        y = protein_de[, "log2_fold_change_protein"],
                        pch = 21, fill = "#6699CC",
                        col = gray(0.2, alpha = 0.5), cex = 2)
           panel.text(both_de$gene_name_RNA,
                      x = both_de$log2_fold_change_RNA + 0.07,
                      y = both_de$log2_fold_change_protein,
                      col = "#882255")

}
)
)
dev.off()


## -----
## output supp tables
head(prot_all)

prot_out <- prot_all
colnames(prot_out) <- gsub(pattern = "_protein_protein",
                           replacement = "",
                           colnames(prot_all))

## write proteomics output to table
prot_out <- prot_out[, c(1, 4, 2, 5, 3, 6)]
ord_prot <- order(prot_out$padj, decreasing = F)
prot_out <- prot_out[ord_prot, c(2, 1, 3:6)]
head(prot_out)
write.table(x = prot_out,
            file = "~/emacs/N-end_Rule_QTL_paper/tables_drafts/all_protein_out.csv",
            append = F, quote = F, sep = ",",
            row.names = F, col.names = T)

save(prot_all, file = "~/emacs/N-end_Rule_QTL_paper/tables_drafts/all_proteomics_out.rdata")
save(res_sva, file = "~/emacs/N-end_Rule_QTL_paper/tables_drafts/all_RNA_out.rdata")
## load("~/emacs/N-end_Rule_QTL_paper/tables_drafts/all_proteomics_out.rdata")
## load("~/emacs/N-end_Rule_QTL_paper/tables_drafts/all_RNA_out.rdata")

RNA_out <- res_sva
RNA_out <- RNA_out[, c(10, 11, 1, 2, 4, 5, 6, 9)]
ord_rna <- order(RNA_out$padj, decreasing = F)
RNA_out <- RNA_out[ord_rna, ]
RNA_out$fold_change <- 2^(RNA_out$log2FoldChange)
RNA_out <- RNA_out[, c("geneID", "gene_name", "fold_change", "log2FoldChange", "padj", "log10_padj")]

write.table(x = RNA_out,
            file = "~/emacs/N-end_Rule_QTL_paper/tables_drafts/all_RNA_out.csv",
            append = F, quote = F, sep = ",",
            row.names = F, col.names = T)



## -----
## output lists for gene ontology analysis
ontology_dir <- paste0(base_dir, "gene_ontology/")

## RNA reference set
RNA_ref_set <- as.data.frame(rownames(RNA_all),
                             stringsAsFactors = F)

## all rows should have gene names:
nrow(RNA_ref_set)
length(RNA_ref_set[(nchar(RNA_ref_set[, 1]) > 5), 1])

write.table(x = RNA_ref_set,
            file = paste0(ontology_dir, "RNA_ref_set.txt"),
            append = F, quote = F, sep = "",
            row.names = F, col.names = F)

RNA_ref_set <- read.table(file = paste0(ontology_dir, "RNA_ref_set.txt"),
                          header = F,
                          quote = "",
                          stringsAsFactors = F)

## RNA test set 
RNA_test_set <- as.data.frame(rownames(res_out),
                             stringsAsFactors = F)

## all rows should have gene names:
nrow(RNA_test_set)
length(RNA_test_set[(nchar(RNA_test_set[, 1]) > 5), 1])

write.table(x = RNA_test_set,
            file = paste0(ontology_dir, "RNA_test_set.txt"),
            append = F, quote = F, sep = "",
            row.names = F, col.names = F)

RNA_test_set <- read.table(file = paste0(ontology_dir, "RNA_test_set.txt"),
                           header = F,
                           quote = "",
                           stringsAsFactors = F)

## protein reference set
prot_ref_set <- as.data.frame(rownames(prot_all),
                             stringsAsFactors = F)

## all rows should have gene names:
nrow(prot_ref_set)
length(prot_ref_set[(nchar(prot_ref_set[, 1]) > 5), 1])

write.table(x = prot_ref_set,
            file = paste0(ontology_dir, "prot_ref_set.txt"),
            append = F, quote = F, sep = "",
            row.names = F, col.names = F)

prot_ref_set <- read.table(file = paste0(ontology_dir, "prot_ref_set.txt"),
                           header = F,
                           quote = "",
                           stringsAsFactors = F)


## protein test set
prot_test_set <- as.data.frame(rownames(prot_filtered),
                             stringsAsFactors = F)

## all rows should have gene names:
nrow(prot_test_set)
length(prot_test_set[(nchar(prot_test_set[, 1]) > 5), 1])

write.table(x = prot_test_set,
            file = paste0(ontology_dir, "prot_test_set.txt"),
            append = F, quote = F, sep = "",
            row.names = F, col.names = F)

prot_test_set <- read.table(file = paste0(ontology_dir, "prot_test_set.txt"),
                            header = F,
                            quote = "",
                            stringsAsFactors = F)


## both ref set
both_ref_set <- as.data.frame(rownames(both_all_only),
                              stringsAsFactors = F)

## all rows should have gene names:
nrow(both_ref_set)
length(both_ref_set[(nchar(both_ref_set[, 1]) > 5), 1])

write.table(x = both_ref_set,
            file = paste0(ontology_dir, "both_ref_set.txt"),
            append = F, quote = F, sep = "",
            row.names = F, col.names = F)

both_ref_set <- read.table(file = paste0(ontology_dir, "both_ref_set.txt"),
                           header = F,
                           quote = "",
                           stringsAsFactors = F)

## both test set
both_test_set <- as.data.frame(rownames(both_only),
                              stringsAsFactors = F)

## all rows should have gene names:
nrow(both_test_set)
length(both_test_set[(nchar(both_test_set[, 1]) > 5), 1])

write.table(x = both_test_set,
            file = paste0(ontology_dir, "both_test_set.txt"),
            append = F, quote = F, sep = "",
            row.names = F, col.names = F)

both_test_set <- read.table(file = paste0(ontology_dir, "both_test_set.txt"),
                            header = F,
                            quote = "",
                            stringsAsFactors = F)

## systematic names don't seem to work well
## with pantherdb, so write out to entrez instead
## list of alternative geneIDs we can use: 
columns(org.Sc.sgd.db)

set_names <- c("RNA_reference_set_ensembl.txt", "protein_reference_set_ensembl.txt", 
               "RNA_test_set_ensembl.txt", "protein_test_set_ensembl.txt", 
               "both_reference_set_ensembl.txt", "both_test_set_ensembl.txt") 

set_list <- list(RNA_ref_set, prot_ref_set,
                 RNA_test_set, prot_test_set,
                 both_ref_set, both_test_set)

Map(f = function(x, y, ...) {
        ## convert systematic names to entrez via 'select': 
        out <- select(org.Sc.sgd.db,
                      keys = x$V1,
                      columns = c("ENTREZID"))
        ## keep only the entrezid:
        out <- out$SGD
        write.table(x = out, file = paste0(ontology_dir, y),
                    append = F, quote = F, sep = "",
                    row.names = F, col.names = F)
    }, x = set_list, y = set_names)
