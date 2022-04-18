## -----
## code to make a heatmap of QTL peaks across multiple reporters 

## -----
## <<heatmap_setup>>
## ifelse to determine which system we're on
## if not on home machine, assume msi
system   <- Sys.info()["nodename"]
base_dir <- "~/data/illumina/"
## determine if we're on msi; set wd appropriately 
base_dir <- ifelse(system != "mahlon-linux",
                   "/home/albertf/mahlon/data/illumina/",
                   base_dir)

#############
## USER INPUT
#############
## set the specific directory you'll work 
## in and name the comparison table
## TRAILING SLASH AT END OF DIR
## below, your project, e.g.,
## "2021.08.17_FPFA002_TDH3pr_Arg_N-end_TFT_sorts/"
proj           <- "2021.10.30_all_UPS_rdata/peaks/merged_delta_AF_peak_tables/"
proj_dir       <- paste0(base_dir, proj)
peaks_table    <- read.table(file = paste0(proj_dir, "all_merged_peaks.csv"),
                             sep = ",", header = T)
#################
## END USER INPUT
#################


## -----
## <<avg_peak_table>>
## now make a table of the average LOD, delta_AF, and peak 
## positions that we'll use to construct the heatamp 
avg_table <- data.frame(reporter = peaks_table$reporter,
                        chr = peaks_table$chr,
                        LOD = (peaks_table$rep_1_LOD + peaks_table$rep_2_LOD) / 2,
                        delta_AF = (peaks_table$rep_1_delta_AF + peaks_table$rep_2_delta_AF) / 2,
                        left_Index = (peaks_table$rep_1_left_Index + peaks_table$rep_2_left_Index) / 2,
                        max_Index = (peaks_table$rep_1_max_Index + peaks_table$rep_2_max_Index) / 2 ,
                        right_Index = (peaks_table$rep_1_right_Index + peaks_table$rep_2_right_Index) / 2)

## write out the averaged table
write.table(x = avg_table,
            file = paste0(proj_dir, "averaged_merged_all_peaks_table.csv"),
            append = F, quote = F, sep = ",",
            row.names = F, col.names = T)

avg_table <- read.table(file = paste0(proj_dir, "averaged_merged_all_peaks_table.csv"),
                        header = T, sep = ",")

## -----
## <<common_functions_code_
## needed tables, files, etc....
## SNPs is a giant table w/ SNP positions 
SNPs           <- read.table("~/QTL_scripts/SNPs_Maggie_170809_BY_positions.txt",
                             stringsAsFactors=FALSE,
                             head=FALSE)

## as of 8/31/17, the SNPs seem not to be fully
## filtered and are out of sorting order
## I think the next section duplicates this code
## w/o the loop.  Frank's comment suggests he
## thinks the 'for' loop isn't working, but I
## think it is.  Does Maggie's code fix it?
## I guess not, given the dates.  
for (thisChr in unique(SNPs[,1])){
    SNPs[SNPs[,1] == thisChr, 2] <- sort(SNPs[SNPs[,1] == thisChr, 2])
}
SNPs <- rbind(SNPs[SNPs[, 1] == "chrI", ],
              SNPs[SNPs[, 1] == "chrII", ],
              SNPs[SNPs[, 1] == "chrIII", ],
              SNPs[SNPs[, 1] == "chrIV", ],
              SNPs[SNPs[, 1] == "chrV", ],
              SNPs[SNPs[, 1] == "chrVI", ],
              SNPs[SNPs[, 1] == "chrVII", ],
              SNPs[SNPs[, 1] == "chrVIII", ],
              SNPs[SNPs[, 1] == "chrIX", ],
              SNPs[SNPs[, 1] == "chrX", ],
              SNPs[SNPs[, 1] == "chrXI", ],
              SNPs[SNPs[, 1] == "chrXII", ],
              SNPs[SNPs[, 1] == "chrXIII", ],
              SNPs[SNPs[, 1] == "chrXIV", ],
              SNPs[SNPs[, 1] == "chrXV", ],
              SNPs[SNPs[, 1] == "chrXVI", ])

library("lattice")
library("latticeExtra")
library("RColorBrewer")
library("grid")

## check for Bioconductor and install if not available
ifelse(!requireNamespace("BiocManager", quietly = TRUE),
       install.packages("BiocManager",
                        dependencies = TRUE,
                        repos = "http://cran.wustl.edu/",
                        quiet = TRUE),
       paste0("Bioconductor available"))
require("BiocManager")

bioc_package_installer <- function(x){if(!requireNamespace(x))
                                          BiocManager::install(x,
                                                               INSTALL_opts = '--no-lock')}
bioc_package_installer("VariantAnnotation")

library("VariantAnnotation")
source("~/QTL_scripts/gTest.R")
source("~/QTL_scripts/x_qtl_seq_functions_170831.R")
source("~/QTL_scripts/mp_JB_170901.R")
source("~/QTL_scripts/peaksFromVector.R")
("~/QTL_scripts/chr_cutoffs.Rdata")
load("~/QTL_scripts/chr_labels.Rdata")

## data frame with all yeast genes, plus
## chr, pos., strand, and names 
geneInfo <- read.table("~/QTL_scripts/ensemblGenes_ensembl83_160307_MOD.txt",
                       stringsAsFactors=FALSE,
                       sep="\t",
                       header=TRUE)
## rownames become systemtatic names 
rownames(geneInfo) <- geneInfo[,"geneID"]
## "geneName" is the common name, e.g., 'HOG1'
## for some (many?) rows of 'geneInfo', there
## is no 'geneName', so it's just an empty string
## e.g., head(allNames)
allNames <- geneInfo[, "geneName"]

names(allNames) <- geneInfo[,1]
allNames[which(allNames == "")] <- names(allNames)[which(allNames == "")]
allNamesInv <- names(allNames)
names(allNamesInv) <- allNames

sepBetweenChr <- 0
trimFromEnd   <- 15e3
obsMin        <- 10
LoessSpan     <- 0.1
## same as in Albert 2014
AF_thres       <- 0.09653124
## multipool LOD threshold, our usual value 
multi_thres    <- 4.5


## -----
## <<heatmap_setup>>
## build a dummy dataframe to drop QTLs into
## ~120 bins if you use 100 kb windows for the map

## read in chr lengths
## need this for building the heatmap
chr_lengths <- read.table("~/QTL_scripts/sacCer3ChromLengths.txt", header = F)

## don't need mitochondrial chromosome 
chr_lengths <- chr_lengths[1:16, ]

## add the gcoords and numeric (not roman) chromosomes:
chr_lengths$chr <- 1:16
chr_lengths$gcoords <- getGcoords(chr = 1:16,
                                  pos = chr_lengths$chr,
                                  spacing = sepBetweenChr)
chr_lengths <- chr_lengths[, c(3, 1, 2, 4)]
names(chr_lengths) <- c("chr", "chr_r", "length", "gcoords")

## <<build_heatmap_dataframe>>
k <- 1
reporters <- levels(avg_table$reporter)

final <- list()
for (k in 1:length(reporters)) {    

## make a list of 16 containing the chr and bin as a dataframe
## have to round up the chr lengths so the bins are all the same
## size (1e5 bp)

## we need to make an extra bin for ea. chromosome, since the
## lengths don't end in 1e5 multiples.  'round_up' allows you
## to do this using a multiple of your choice.
## So, chrI = 230218 bp and round_up(230218, 1e5) = 3e5
round_up <- function(from, to) {
    ceiling(from / to) * to
    }

round_down <- function(from, to) {
    floor(from / to) * to
    }    

## assign each bin to the correct chromosome
chr_indices <- unlist(sapply(X = 1:16, FUN = function(x) {
                                 rep(x,
                                     times = length(seq(from = 1e5,
                                                        to = round_up(chr_lengths$length[x],
                                                                      to = 1e5),
                                                        by = 1e5)))
                      }))

## create chromosome bins - this will serve as a conditioning factor
bins <- seq(from = 1e5,
            to = length(chr_indices) * 1e5,
            by = 1e5)

## when we assign peaks to bins, we need to add the length of
## all preceding chromosomes plus the peak position.  we'll use
## the bin lengths, not the chromosome lengths for this purpose.
## so, a QTL on chrII at base 150 = 3e5 (length of chr I bins) + 150
bin_sums <- c(0, unlist(sapply(X = 2:16, FUN = function(x) {
                              max(bins[chr_indices == x - 1])
                              })))

## assign the peaks for ea. reporter to a dataframe 
out <- avg_table[avg_table$reporter == reporters[k], ]

## assign each peak a bin for the heatmap 
for (b in 1:nrow(out)) {
    out$bin[b] <- round_down(from = out$max_Index[b] + bin_sums[out$chr[b]],
                             to = 1e5)
}

## dummy dataframe that contains all bins    
chr_bins <- data.frame(reporter = rep(reporters[k], length(bins)),
                       chr = chr_indices,
                       LOD = rep(0, length(bins)),
                       delta_AF = rep(0, length(bins)),
                       left_Index = rep(0, length(bins)),
                       max_Index = rep(0, length(bins)),
                       right_Index = rep(0, length(bins)),
                       bin = bins)

## the bins above aren't exact, so assign them to one
## of the bins we defined using 'which.min'
## assign peaks to bins via 'which.min'
for (r in 1:nrow(out)) {
    ind <- out$chr[r]
    chrs <- chr_bins$bin[chr_bins$chr == ind & chr_bins$bin > out$bin[r]]
    out$bin[r] <- chrs[which.min(abs(out$bin[r] - chrs))]
}

## if there's a peak, replace the current row with the peak row
## if not, just leave everything at 0
for (n in 1:nrow(chr_bins)) {
    test <- min(abs(chr_bins[n, "bin"] - out$bin))
    ind  <- which.min(abs(chr_bins[n, "bin"] - out$bin))
    chr_bins[n, ] <- if(test == 0) out[ind, ] else chr_bins[n, ]
}

## assign the dataframe w/ each peak and peak bin position to a list
## that we'll collapse after we've run through all the reporters 
final[[k]] <- chr_bins
}

## collapse output to a single frame
test <- do.call("rbind", final)

## remove the "_" from the factor level
levels(test$reporter) <- gsub(pattern = "_", replacement = " ",
                              x = levels(test$reporter))

## square root transform LOD scores for better visualization
test$LOD_mult <- test$delta_AF / abs(test$delta_AF)

## make the LOD negative if the high - low BY allele frequency < 1
## 2021.01.20 - this doesn't really work and won't be in the final
## heatmap 
for (i in 1:length(test$LOD_mult)) {
    test$LOD_mult[i] <- if(is.nan(test$LOD_mult[i]) == T) 0 else test$LOD_mult[i]
    }

test$sq_LOD <- sqrt(test$LOD) * test$LOD_mult

## <<alphabetical_order>>
## put the amino acids in alphabetical order,
## but also split by pathway (Arg/N-end vs. Ac/N-end)
## have to reverse it due to the way lattice plots
## Arg/N-end: arg, asn, asp, gln, glu, his, ile, leu, lys, phe, trp, tyr
## Ac/N-end: ala, cys, gly, met, pro, ser, thr, val
aa_order  <- c(2, 3, 4, 6, 7, 9, 10, 11, 12, 14, 18, 19, 1, 5, 8, 13, 15, 16, 17, 20)
## aa_order[order(aa_order)] == 1:20
test$aa_factor <- factor(test$reporter,
                    levels = levels(test$reporter)[rev(aa_order)])
## levels(test$aa_factor)

## <<type_I_II_order>>
## still alphabetical, but w/ Arg/N-end
## split by type I vs II ubr1 rec site
## type I  = arg asn asp gln glu his lys 
## type II = ile leu phe trp tyr 
## Ac/N-end: ala, cys, gly, met, pro, ser, thr, val
type_order <- c(2, 3, 4, 6, 7, 9, 12, 10, 11, 14, 18, 19, 1, 5, 8, 13, 15, 16, 17, 20)
## type_order[order(type_order)] == 1:20
## levels(test$reporter)[type_order]
test$type_factor <- factor(test$reporter,
                           levels = levels(test$reporter)[rev(type_order)])
## levels(test$type_factor)

## <<dynamic_range_order>>
## timer dynamic range order (high to low):
## Arg/N-end: asn, trp, asp, phe, tyr, lys, arg, his, gln, ile, leu, glu
## Ac/N-end: cys, met, pro, val, ser, thr, ala, gly
deg_order <- c(3, 18, 4, 14, 19, 12, 2, 9, 6, 10, 11, 7, 5, 13, 15, 20, 16, 17, 1, 8)
test$deg_factor <- factor(test$reporter,
                          levels = levels(test$reporter)[rev(deg_order)])
## deg_order[order(deg_order)] == 1:20
## levels(test$reporter)[deg_order]

## position where the labels for the x axis (chr) go
chr_labels <- sapply(X = 1:16, FUN = function(x) {
                         mean(test$bin[test$chr == x])
                     })

## position where the dividing lines on the x axis go
chr_cutoffs <- sapply(X = 1:15, FUN = function(x) {
                          max(test$bin[test$chr == x]) + 5e4
})


sq_LOD_ats  <- c(seq(from = -16, to = -2, length.out = 4),
                 seq(from = 2, to = 16, length.out = 4))
sq_LOD_cols <- c(brewer.pal(10, "RdBu")[1:4],
                 "white", 
                 brewer.pal(10, "RdBu")[7:10])


## this object demarcates the bins on the colorkey
## 0.075 is the smallest delta_AF value for the QTLs
delta_AF_ats  <- c(seq(from = -0.7, to = -0.075, length.out = 4),
                   seq(from = 0.075, to = 0.7, length.out = 4))
## "RdBu" is colorlbind friendly:
## display.brewer.all(colorblindFriendly = TRUE)
## "PiYG" should work for the TFT plots
## look at a palette
## display.brewer.pal(10, "RdBu")
delta_AF_cols <- c(brewer.pal(10, "RdBu")[1:4],
                   "white", 
                   brewer.pal(10, "RdBu")[7:10])
{
pdf(file = "~/Desktop/2022.01.20_delta_AF_heatmap_final.pdf",
    height = 5.5, width = 11)
print( 
levelplot(delta_AF ~ bin * type_factor,
          data = test,
          xlab = "Chromosome",
          col.regions = delta_AF_cols,
          at = delta_AF_ats,
          colorkey = list(col = delta_AF_cols,
                          at = 0:7,
                          labels = as.character(c(-0.70, -0.50, -0.25, -0.075,
                                                  0.075, 0.25, 0.50, 0.70)),
                          ## reduce the height of the colorkey a bit:
                          height = 0.5,
                          title = "delta AF"),
          scales = list(x = list(at = chr_labels,
                                 tck = c(1, 0),
                                 alternating = 1,
                                 labels = as.roman(1:16)),
                        y = list(alternating = 1,
                                 tck = c(1, 0))),
          par.settings = list(layout.widths = list(left.padding = 3,
                                                   right.padding = 3.5),
                              par.ylab.text = list(cex = 1.25,
                                                   col = "white")),
          panel = function(...){
              panel.levelplot(...)
              panel.abline(v = chr_cutoffs,
                           lty = 1, col = gray(0.4))
              panel.abline(h = c(8.5, 13.5),
                           lty = 1, col = gray(0.4))
          })
)
## reporter text 
grid.text(x = rep(0.25, 3),
          y = c(4.35, 2.97, 1.6),
          label = c("Type I Arg/N-end\nReporters",
                    "Type II Arg/N-end\nReporters",
                    "Ac/N-end\nReporters"),
          default.units = "in",
          rot = 90,
          gp = gpar(cex = 0.85),
          vp = "plot_01.toplevel.vp")
## heatmap key text
grid.text(x = 10.85,
          y = 2.92,
          label = expression(paste(Delta, "Allele Frequency (High vs. Low Pool)")),
          default.units = "in",
          rot = 90,
          gp = gpar(cex = 0.85),
          vp = "plot_01.toplevel.vp")
## type I Arg/N-end reporter vertical line
grid.lines(x = c(0.5, 0.5),
           y = c(3.65, 4.95),
           default.units = "in",
           gp = gpar(col = gray(0.4), lwd = 1.5))
## type II Arg/N-end reporter vertical line
grid.lines(x = c(0.5, 0.5),
           y = c(2.55, 3.45),
           default.units = "in",
           gp = gpar(col = gray(0.4), lwd = 1.5))
## Ac/N-end reporter vertical line
grid.lines(x = c(0.5, 0.5),
           y = c(0.81, 2.35),
           default.units = "in",
           gp = gpar(col = gray(0.4), lwd = 1.5))
dev.off()
}

## find the viewport for the plot you're working on
grid.ls(viewport = T, grobs = F)


## this object demarcates the bins on the colorkey
## 0.075 is the smallest delta_AF value for the QTLs
delta_AF_ats  <- c(seq(from = -0.7, to = -0.075, length.out = 4),
                   seq(from = 0.075, to = 0.7, length.out = 4))
## "RdBu" is colorlbind friendly:
## display.brewer.all(colorblindFriendly = TRUE)
## "PiYG" should work for the TFT plots
## look at a palette
## display.brewer.pal(10, "RdBu")
delta_AF_cols <- c(brewer.pal(10, "RdBu")[1:4],
                   "white", 
                   brewer.pal(10, "RdBu")[7:10])
{
pdf(file = "~/Desktop/2022.01.20_delta_AF_heatmap_final.pdf",
    height = 5.5, width = 11)
print( 
levelplot(delta_AF ~ bin * type_factor,
          data = test,
          xlab = "chromosome",
          col.regions = delta_AF_cols,
          at = delta_AF_ats,
          colorkey = list(col = delta_AF_cols,
                          at = 0:7,
                          labels = as.character(c(-0.70, -0.50, -0.25, -0.075,
                                                  0.075, 0.25, 0.50, 0.70)),
                          ## reduce the height of the colorkey a bit:
                          height = 0.5,
                          title = "delta AF"),
          scales = list(x = list(at = chr_labels,
                                 tck = c(1, 0),
                                 alternating = 1,
                                 labels = as.roman(1:16)),
                        y = list(alternating = 1,
                                 tck = c(1, 0))),
          par.settings = list(layout.widths = list(left.padding = 3,
                                                   right.padding = 3.5),
                              par.ylab.text = list(cex = 1.25,
                                                   col = "white")),
          panel = function(...){
              panel.levelplot(...)
              panel.abline(v = chr_cutoffs,
                           lty = 1, col = gray(0.4))
              panel.abline(h = c(8.5, 13.5),
                           lty = 1, col = gray(0.4))
          })
)
## reporter text 
grid.text(x = rep(0.25, 3),
          y = c(4.35, 2.97, 1.6),
          label = c("Type I Arg/N-end\nReporters",
                    "Type II Arg/N-end\nReporters",
                    "Ac/N-end\nReporters"),
          default.units = "in",
          rot = 90,
          gp = gpar(cex = 0.85),
          vp = "plot_01.toplevel.vp")
## heatmap key text
grid.text(x = 10.85,
          y = 2.92,
          label = expression(paste(Delta, "Allele Frequency (High vs. Low Pool)")),
          default.units = "in",
          rot = 90,
          gp = gpar(cex = 0.85),
          vp = "plot_01.toplevel.vp")
## type I Arg/N-end reporter vertical line
grid.lines(x = c(0.5, 0.5),
           y = c(3.65, 4.95),
           default.units = "in",
           gp = gpar(col = gray(0.4), lwd = 1.5))
## type II Arg/N-end reporter vertical line
grid.lines(x = c(0.5, 0.5),
           y = c(2.55, 3.45),
           default.units = "in",
           gp = gpar(col = gray(0.4), lwd = 1.5))
## Ac/N-end reporter vertical line
grid.lines(x = c(0.5, 0.5),
           y = c(0.81, 2.35),
           default.units = "in",
           gp = gpar(col = gray(0.4), lwd = 1.5))
dev.off()
}
