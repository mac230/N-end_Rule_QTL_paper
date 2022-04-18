## -----
## setup for reading in results
library("tidyverse")
base_dir <- "~/emacs/N-end_Rule_QTL_paper/UBR1_RNA-seq_analysis/fastp"
ontology_dir <- paste0(base_dir, "/gene_ontology/GO_results/")


## -----
## read in each of the 3 results files
## 'biological process'
## the tables output by panther have strange
## formatting, so have to reformat them a 
## fair bit to get them into a readable format
biop <- read.table(file = dir(ontology_dir, full.names = T)[1],
                   sep = "\t", skip = 6, header = T, stringsAsFactors = F)

col_splits <- unlist(lapply(X = strsplit(x = colnames(biop)[3:ncol(biop)],
                                         split = "RNA_test_set_ensembl.txt.."),
                            FUN = function(x) {
                                x[2]
                            }))

colnames(biop) <- c("GO_Term", "ref_list", col_splits)

go_term <- unlist(lapply(X = strsplit(x = biop$GO_Term,
                                      split = " \\(.*"),
                         FUN = function(x) {
                             x[1]
                         }))

biop$go_term <- go_term
biop$ontology <- rep("bio_process", nrow(biop))


## 'cellular component'
cell_comp <- read.table(file = dir(ontology_dir, full.names = T)[2],
                        sep = "\t", skip = 6, header = T, stringsAsFactors = F)

col_splits <- unlist(lapply(X = strsplit(x = colnames(cell_comp)[3:ncol(cell_comp)],
                                         split = "RNA_test_set_ensembl.txt.."),
                            FUN = function(x) {
                                x[2]
                            }))

colnames(cell_comp) <- c("GO_Term", "ref_list", col_splits)

go_term <- unlist(lapply(X = strsplit(x = cell_comp$GO_Term,
                                      split = " \\(.*"),
                         FUN = function(x) {
                             x[1]
                         }))

cell_comp$go_term <- go_term
cell_comp$ontology <- rep("cellular_component", nrow(cell_comp))


## 'reactome pathway'
reactome <- read.table(file = dir(ontology_dir, full.names = T)[3],
                        sep = "\t", skip = 6, header = T, stringsAsFactors = F)

col_splits <- unlist(lapply(X = strsplit(x = colnames(reactome)[3:ncol(reactome)],
                                         split = "RNA_test_set_ensembl.txt.."),
                            FUN = function(x) {
                                x[2]
                            }))

colnames(reactome) <- c("GO_Term", "ref_list", col_splits)

go_term <- unlist(lapply(X = strsplit(x = reactome$GO_Term,
                                      split = " \\(.*"),
                         FUN = function(x) {
                             x[1]
                         }))

reactome$go_term <- go_term
reactome$ontology <- rep("reactome", nrow(reactome))


## -----
## merge into a single dataframe
all_terms <- list(biop, cell_comp, reactome)
all_terms <- do.call("rbind", all_terms)


## pick specific terms from the full list for plotting 
plot_rows <- c(7, 5, 4, 2, 1, 21, 19, 17)
plot_rows <- rev(plot_rows)
plot_terms <- all_terms[plot_rows, ]
 
plot_terms$go_fac <- factor(plot_terms$go_term,
                            levels = plot_terms$go_term,
                            labels = str_to_title(plot_terms$go_term))

plot_terms$fold.Enrichment.

go_cols <- ifelse(plot_terms$ontology == "bio_process",
                  gray(0.95),
                  "#D6EDE9")

better_y_lab <- c("Hsp90 Chaperone Cycle for\nSteroid Hormone Receptors", 
                  "Hsf1-Dependent\nTransactivation", 
                  "Cellular Response\nto Heat Stress", 
                  "Chaperone Cofactor-\nDependent Protein Refolding", 
                  "De Novo Posttranslational\nProtein Folding", 
                  "Chaperone-Mediated\nProtein Folding", 
                  "Cellular Response to\nMisfolded or Unfolded Protein",
                  "Ion Transmembrane\nTransport")

pdf(file = paste0("~/emacs/N-end_Rule_QTL_paper/figures_drafts/",
                  gsub(pattern = " .*",
                       replacement = "",
                       x = Sys.time()),
                  "_gene_ontology_plot_final.pdf"), width = 7)
print(
barchart(go_fac ~ fold.Enrichment.,
         data = plot_terms,
         horizontal = T,
         xlim = c(-2, 32.5),
         axes = T,
         stack = F,
         scales = list(tck = c(1, 0),
                       x = list(at = seq(from = 0, to = 30, by = 5),
                                cex = 1.25),
                       y = list(labels = rep("",
                                             length(levels(plot_terms$go_fac))),
                                cex = 1)),
         key = list(corner = c(0.99, 0.99),
                    rectangles = list(size = 5,
                                      border = T,
                                      col = c(gray(0.95),
                                              "#D6EDE9"),
                                      lwd = 10),
                    background = "white",
                    between = 0.5,
                    text = list(labels = c("GO Biological Process",
                                           "Reactome Pathway"),
                                cex = 1.2)),
         par.settings = list(par.xlab.text = list(cex = 1.25),
                             par.ylab.text = list(cex = 1.1),
                             ## allow text to be placed 
                             ## outside the plotting area 
                             clip = list(panel = F),
                             ## padding for text 
                             layout.widths = list(left.padding = 30)), 
         box.ratio = 2.5,
         col = go_cols,
         xlab = "Fold Enrichment",
         panel = function(...) {
             panel.abline(v = seq(from = 5, to = 30, by = 5),
                          lty = 1, col = gray(0.9))
             panel.barchart(...)
             ## panel.abline(h = 3.5, lty = 2, col = gray(0.4), lwd = 2)
             panel.text(paste0(rev(plot_terms[, 3]), " / ",
                               rev(plot_terms[, 4])),
                        x = rep(1, 8),
                        y = seq(from = 8, to = 1, by = -1))
             panel.text(paste(better_y_lab), x = -2.5,
                        y = seq(from = 1, to = 8, by = 1),
                        adj = 1,
                        cex = 1.2) 
         })
)

dev.off()
