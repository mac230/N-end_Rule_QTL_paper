## -----
## cdf plot of QTLs
s
## -----
## setup
library("lattice")
library("latticeExtra")
library("RColorBrewer")
library("grid")
library("VariantAnnotation")


## -----
## read in peaks
base_dir <- "~/data/illumina/"
table_dir  <- "2021.10.30_all_UPS_rdata/peaks/merged_delta_AF_peak_tables/"
table_file <- "averaged_merged_all_peaks_table.csv"
avg_table <- read.csv(file = paste0(base_dir,
                                    table_dir,
                                    table_file),
                      header = T, sep = ",",
                      quote = "")
str(avg_table) ; nrow(avg_table)


## -----
## sort peaks by pathway
path_test <- function(x) {
    ifelse(x == "Ala_TFT" | x == "Cys_TFT" | x == "Gly_TFT" |
           x == "Met_TFT" | x == "Pro_TFT" | x == "Ser_TFT" |
           x == "Thr_TFT" | x == "Val",
           "Ac/N-end Pathway", "Arg/N-end Pathway")
}


## create Ac/ vs. Arg/N-end pathway variable
avg_table$path <- path_test(avg_table$reporter)
avg_table$path <- as.factor(avg_table$path)
nrow(avg_table[avg_table$path == "Arg/N-end Pathway", ])
nrow(avg_table[avg_table$path == "Ac/N-end Pathway", ])

pdf(file = paste0("~/emacs/N-end_Rule_QTL_paper/figures_drafts/",
                  gsub(pattern = " .*",
                       replacement = "",
                       x = Sys.time()),
                  "_delta_AF_ecdf.pdf"))

## "~/emacs/N-end_Rule_QTL_paper/figures_drafts/QTL_ecdf_by_pathway.pdf",
##     width = 10, height = 5)

print(
ecdfplot(~ delta_AF,
         groups = path,
         data = avg_table,
         col = c("#6699CC", "#AA4499"),
         lwd = 1.5,
         xlim = c(-0.75, 0.75),
         ylim = c(-0.1, 1.05),
         xlab = expression(paste("Effect Size (", Delta, "RM Allele Frequency)")),
         ylab = "Cumulative Density",
         key = list(## corner = c(0.005, 0.935),
                    lines = list(size = 3.5,
                                 lwd = 3.5,
                                 col = c("#AA4499", "#6699CC")),
                    text = list(labels = c("Arg/N-end Pathway",
                                           "Ac/N-end Pathway"),
                                cex = 1.7),
                    padding.text = 2.5),
         scales = list(x = list(at = seq(from = -0.7, to = 0.7, by = 0.20),
                                labels = sprintf("%.1f", seq(from = -0.7, to = 0.7, by = 0.20)),
                                cex = 1.7),
                       y = list(at = seq(from = 0, to = 1, by = 0.2),
                                cex = 1.7),
                       tck = c(1, 0)),
         par.settings = list(par.ylab.text = list(cex = 1.7),
                             par.xlab.text = list(cex = 1.7)),
         panel = function(...) {
             panel.abline(h = 0.5, v = 0,
                          lwd = 1,
                          col = gray(0.5),
                          lty = 2)
             panel.ecdfplot(...)
             panel.abline(h = c(0, 1), col = gray(0.7))
             panel.text(c("RM allele = Lower UPS Activity",
                          "RM allele = Higher UPS Activity"),
                        x = c(-0.375, 0.375),
                        y = c(-0.05, -0.05),
                        cex = 1.2)
             ## panel.arrows(x1 = -0.74, x2 = -0.5,
             ##              y1 = -0.05, y2 = -0.05,
             ##              length = 0.08,
             ##              col = gray(0.4))
             ## panel.arrows(x1 = -0.01, x2 = -0.165,
             ##              y1 = -0.05, y2 = -0.05,
             ##              length = 0.08,
             ##              col = gray(0.4))
             ## panel.arrows(x1 = 0.74, x2 = 0.5,
             ##              y1 = -0.05, y2 = -0.05,
             ##              length = 0.08,
             ##              col = gray(0.4))
             ## panel.arrows(x1 = 0.01, x2 = 0.165,
             ##              y1 = -0.05, y2 = -0.05,
             ##              length = 0.08,
             ##              col = gray(0.4))
}
)
)
dev.off()


pdf(file = paste0("~/emacs/N-end_Rule_QTL_paper/figures_drafts/",
                  gsub(pattern = " .*",
                       replacement = "",
                       x = Sys.time()),
                  "_LOD_ecdf.pdf"))
print(
ecdfplot(~ LOD,
         groups = path,
         data = avg_table,
         col = c("#6699CC", "#AA4499"),
         lwd = 1.5,
         xlim = c(-5, 255),
         ylim = c(-0.05, 1.05),
         xlab = "LOD",
         ylab = "Cumulative Density",
         key = list(## corner = c(0.005, 0.999),
                    lines = list(size = 3.5,
                                 lwd = 3.5,
                                 col = c("#AA4499", "#6699CC")),
                    text = list(labels = c("Arg/N-end Pathway",
                                           "Ac/N-end Pathway"),
                                cex = 1.7),
                    padding.text = 2.5),
         scales = list(x = list(at = seq(from = 0, to = 250, by = 25),
                                labels = as.character(seq(from = 0, to = 250, by = 25)),
                                cex = 1.7),
                       y = list(at = seq(from = 0, to = 1, by = 0.2),
                                cex = 1.7),
                       tck = c(1, 0)),
         par.settings = list(par.ylab.text = list(cex = 1.7),
                             par.xlab.text = list(cex = 1.7)),
         panel = function(...) {
             ## panel.abline(h = 0.5, v = 0,
             ##              lwd = 1,
             ##              col = gray(0.7),
             ##              lty = 2)
             ## panel.abline(h = c(0, 1), col = gray(0.2))
             panel.ecdfplot(...)
             panel.abline(h = c(0, 1), col = gray(0.7))             
             ## panel.text(c("RM allele = Lower UPS Activity",
             ##              "RM allele = Higher UPS Activity"),
             ##            x = c(-0.25, 0.25),
             ##            y = c(-0.05, -0.05),
             ##            cex = 1.2)
             ## panel.arrows(x1 = -0.74, x2 = -0.5,
             ##              y1 = -0.05, y2 = -0.05,
             ##              length = 0.08,
             ##              col = gray(0.4))
             ## panel.arrows(x1 = -0.01, x2 = -0.165,
             ##              y1 = -0.05, y2 = -0.05,
             ##              length = 0.08,
             ##              col = gray(0.4))
             ## panel.arrows(x1 = 0.74, x2 = 0.5,
             ##              y1 = -0.05, y2 = -0.05,
             ##              length = 0.08,
             ##              col = gray(0.4))
             ## panel.arrows(x1 = 0.01, x2 = 0.165,
             ##              y1 = -0.05, y2 = -0.05,
             ##              length = 0.08,
             ##              col = gray(0.4))
}
))
dev.off()
