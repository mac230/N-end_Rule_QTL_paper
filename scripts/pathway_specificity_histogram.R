## -----
## code to analyze how many reporters per QTL we see


## -----
## setup
library("lattice")
library("latticeExtra")
library("RColorBrewer")
library("grid")
library("VariantAnnotation")
source("~/QTL_scripts/gTest.R")
source("~/QTL_scripts/x_qtl_seq_functions_170831.R")
source("~/QTL_scripts/mp_JB_170901.R")
source("~/QTL_scripts/peaksFromVector.R")
load("~/QTL_scripts/chr_labels.Rdata")
hist_colors <- c(gray(0.5), gray(0.9))

## -----
## ifelse to determine which system we're on
## if not on home machine, assume msi
system   <- Sys.info()["nodename"]
base_dir <- "~/data/illumina/"
## determine if we're on msi; set wd appropriately
base_dir <- ifelse(system != "mahlon-linux",
                   "/home/albertf/mahlon/data/illumina/",
                   base_dir)
table_dir  <- "2021.10.30_all_UPS_rdata/peaks/merged_delta_AF_peak_tables/"
table_file <- "averaged_merged_all_peaks_table.csv"
avg_table <- read.csv(file = paste0(base_dir,
                                    table_dir,
                                    table_file),
                      header = T, sep = ",",
                      quote = "")

## -----
## calculate the median number of QTLs per reporter
reporter_list <- list()
for (i in 1:length(unique(avg_table$reporter))) {
    reporter <- levels(avg_table$reporter)[i]
    out <- avg_table[avg_table$reporter == reporter, ]
    reporter_list[[i]] <- nrow(out)
}
reporter_list <- do.call("rbind", reporter_list)
median(reporter_list)
## median n. QTLs per reporter = 7

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
                                  spacing = 0)
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
peaks <- do.call("rbind", final)

## remove any bins that don't contain a peak
peaks_final <- peaks[peaks$LOD > 0, ]

## simplify the reporter variable w/ 'gsub'
peaks_final$aa <- gsub(pattern = "_TFT",
                 replacement = "",
                 x = peaks_final$reporter)
peaks_final$aa <- as.factor(peaks_final$aa)

## add in the pathways
path_test <- function(x) {
    ifelse(x == "Ala" | x == "Cys" | x == "Gly" |
           x == "Met" | x == "Pro" | x == "Ser" |
           x == "Thr" | x == "Val",
           "Ac/N-end Pathway", "Arg/N-end Pathway")
}
peaks_final$pathway <- path_test(peaks_final$aa)
peaks_final$pathway <- as.factor(peaks_final$pathway)

## make sure the output matches the delta_AF heatmap:
mark <- unique(peaks_final$bin)
ord <- order(mark)
mark <- mark[ord]
mark_out <- vector()
mark_chr <- vector()
for (i in 1:length(mark)) {
    mark_out[i] <- nrow(peaks_final[peaks_final$bin == mark[i], ])
    mark_chr[i] <- unique(peaks_final$chr[peaks_final$bin == mark[i]])
}
## rbind(mark_chr, mark_out)
##          [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12] [,13]
## mark_chr    1    2    2    2    4    4    4    4    4     4     5     7     7
## mark_out    7    1    3    2    3    1    1    1    6     1     7     1     7
##          [,14] [,15] [,16] [,17] [,18] [,19] [,20] [,21] [,22] [,23] [,24]
## mark_chr     7     7     7     8     8     9     9    10    10    10    10
## mark_out     2     9     8     1     2     2     6    10     1     1     4
##          [,25] [,26] [,27] [,28] [,29] [,30] [,31] [,32] [,33] [,34] [,35]
## mark_chr    11    11    11    11    12    12    12    12    13    13    13
## mark_out     1     2     2     1     1     1    14     1     3     3     3
##          [,36] [,37] [,38] [,39] [,40] [,41] [,42] [,43]
## mark_chr    14    15    15    15    15    16    16    16
## mark_out    10     8     2     1     3     1     4     1

## now, individual pathways
arg_peaks_final <- peaks_final[peaks_final$pathway == "Arg/N-end Pathway", ]
arg_mark <- unique(arg_peaks_final$bin)
arg_ord <- order(arg_mark)
arg_mark <- arg_mark[arg_ord]
arg_mark_out <- vector()
arg_mark_chr <- vector()
for (i in 1:length(arg_mark)) {
    arg_mark_out[i] <- nrow(arg_peaks_final[arg_peaks_final$bin == arg_mark[i], ])
    arg_mark_chr[i] <- unique(arg_peaks_final$chr[arg_peaks_final$bin == arg_mark[i]])
}
## rbind(arg_mark_chr, arg_mark_out)
##              [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12]
## arg_mark_chr    2    2    4    4    4    4    4    7    7     7     8     8
## arg_mark_out    1    3    3    1    1    2    1    1    4     8     1     2
##              [,13] [,14] [,15] [,16] [,17] [,18] [,19] [,20] [,21] [,22] [,23]
## arg_mark_chr     9    10    10    10    10    11    11    12    12    13    13
## arg_mark_out     2     5     1     1     2     1     2     6     1     2     3
##              [,24] [,25] [,26]
## arg_mark_chr    13    14    15
## arg_mark_out     3     7     7


ac_peaks_final <- peaks_final[peaks_final$pathway == "Ac/N-end Pathway", ]
ac_mark <- unique(ac_peaks_final$bin)
ac_ord <- order(ac_mark)
ac_mark <- ac_mark[ac_ord]
ac_mark_out <- vector()
ac_mark_chr <- vector()
for (i in 1:length(ac_mark)) {
    ac_mark_out[i] <- nrow(ac_peaks_final[ac_peaks_final$bin == ac_mark[i], ])
    ac_mark_chr[i] <- unique(ac_peaks_final$chr[ac_peaks_final$bin == ac_mark[i]])
}

## rbind(ac_mark_chr, ac_mark_out)
##             [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12]
## ac_mark_chr    1    2    4    4    5    7    7    7    9    10    10    11
## ac_mark_out    7    2    1    4    7    7    2    5    6     5     2     2
##             [,13] [,14] [,15] [,16] [,17] [,18] [,19] [,20] [,21] [,22] [,23]
## ac_mark_chr    11    12    12    12    13    14    15    15    15    15    16
## ac_mark_out     1     1     1     8     1     3     1     2     1     3     1
##             [,24] [,25]
## ac_mark_chr    16    16
## ac_mark_out     4     1

all_peaks_out <- data.frame(peaks = mark_out, chr = mark_chr)
all_peaks_out$path <- rep(x = "All", times = nrow(all_peaks_out))

arg_peaks_out <- data.frame(peaks = arg_mark_out, chr = arg_mark_chr)
arg_peaks_out$path <- rep(x = "Arg/N-end Pathway", times = nrow(arg_peaks_out))

ac_peaks_out <- data.frame(peaks = ac_mark_out, chr = ac_mark_chr)
ac_peaks_out$path <- rep(x = "Ac/N-end Pathway", times = nrow(ac_peaks_out))

per_q <- list()
per_q[[1]] <- all_peaks_out
per_q[[2]] <- ac_peaks_out
per_q[[3]] <- arg_peaks_out
per_q_final <- do.call("rbind", per_q)
per_q_final$path_factor <- factor(per_q_final$path,
                                  levels = c("All",
                                             "Ac/N-end Pathway",
                                             "Arg/N-end Pathway"),
                                  labels = c("All",
                                             "Ac/N-end Pathway",
                                             "Arg/N-end Pathway"))

per_q_breaks <- seq(from = 0, to = 15, by = 1)
per_q_colors <- c(gray(0.5), gray(0.7), gray(0.9))
per_q_hist <- histogram(~ peaks | path_factor,
                      type = "count",
                      data = per_q_final,
                      xlim = c(0, 15),
                      ylim = c(0, 18),
                      layout = c(3, 1),
                      xlab = "Reporters per QTL",
                      scales = list(x = list(at = seq(from = 0.5, to = 14.5, by = 1),
                                             labels = as.character(1:15)),
                                    y = list(at = seq(from = 0, to = 20, by = 1)),
                                    alternating = F,
                                    tck = c(1, 0)),
                      breaks = per_q_breaks,
                      par.strip.text = list(cex = 1.7),
                      par.settings = list(strip.background = list(col = gray(0.9)),
                                          axis.text = list(cex = 1.5),
                                          par.ylab.text = list(cex = 1.8),
                                          par.xlab.text = list(cex = 1.8)
                                          ),
                      panel = function(x, ..., col) {
                          panel.histogram(x, ..., col = per_q_colors[packet.number()])
                      }
                      )
print(per_q_hist)

pdf(file = paste0("/home/mahlon/emacs/N-end_Rule_QTL_paper/figures_drafts/",
                  "supplementary_figure_004_reporters_per_QTL.pdf"),
    width = 15)
print(per_q_hist)
dev.off()


## -----
## compare RM enriched vs. RM depleted QTLs
pos_neg <- list()
for (i in 1:length(levels(avg_table$reporter))) {
    reporter <- levels(avg_table$reporter)[i]
    dat <- avg_table[avg_table$reporter == reporter, ]
    pos <- nrow(dat[dat$delta_AF > 0, ])
    neg <- nrow(dat[dat$delta_AF < 0, ])
    frac <- pos / (pos + neg)
    pos_neg[[i]] <- frac
}
pos_neg_out <- do.call("rbind", pos_neg)
pos_neg_out <- data.frame(frac = pos_neg_out,
                          reporter = levels(avg_table$reporter))

## for 15 reporters, RM allele associated w/ higher UPS activity
nrow(pos_neg_out[pos_neg_out$frac > 0.5, ])

## 89 of 149 peaks RM allele associated w/ higher UPS activity
nrow(avg_table[avg_table$delta_AF > 0, ]) / 149

avg_table[avg_table$path == "Ac/N-end Pathway" & avg_table$chr == 1, ]
avg_table[avg_table$path == "Ac/N-end Pathway" &
          avg_table$chr == 7 &
          avg_table$max_Index < 250000, ]
