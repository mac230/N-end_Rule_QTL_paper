## -----
## DOA10


## -----
## setup
library("lattice")
dirs <- c("~/data/flow/2022.04.06_DOA10_Thr_1043_CRISPR_swap_flow",
          "~/data/flow/2021.07.20_DOA10_Gly_1047_CRISPR_swap_flow")

## -----
## read in allele effects and create dataframe
frame_dirs <- paste0(dirs, "/dataframes")

effect_list <- list()

for(x in 1:length(frame_dirs)){

    ## load in allele effects for ea. reporter 
    load(dir(frame_dirs[x],
             pattern = ".*effects.R",
             full.names = T),
         verbose = T)
    
    ## assign to the list
    effect_list[[x]] <- allele_effects

    ## remove object before going on to next reporter 
    rm(allele_effects)

}

effects <- do.call("rbind", effect_list)


## -----
## subset and order
effects_sub <- effects[effects$allele != "BY_full", ]

effects_sub$allele_drop <- droplevels(effects_sub$allele)

effects_sub$allele_ordered <- factor(effects_sub$allele,
                                     levels = levels(effects_sub$allele_drop)[c(4, 1, 2, 3)])

effects_final <- effects_sub[order(effects_sub$allele_ordered,
                                   decreasing = F), ]

effects_final <- effects_final[ , -7]

effects_final


## -----
## plot
pdf(file = "~/emacs/N-end_Rule_QTL_paper/figures_drafts/fine-mapping_summary_DOA10_barchart.pdf",
    height = 4.5, width = 7)
barchart(1:nrow(effects_final) ~ effects_final$effect,
         origin = 0,
         ylab = "",
         xlab = "UPS Activity Relative to BY Allele Median (SD units)",
         xlim = c(-0.6, 2.2),
         box.ratio = 2,
         scales = list(tck = c(1, 0),
                       y = list(at = seq(from = 1.5, to = 7.5, by = 2),
                                labels = c(expression(italic("DOA10")*" Y1186F"),
                                           expression(italic("DOA10")*" K1012N"),
                                           expression(italic("DOA10")*" Q410E"),
                                           expression(italic("DOA10")*" RM")),
                                cex = 1.25),
                       x = list(cex = 1.25,
                                at = seq(from = -4.5, to = 4.5, by = 0.5))),
         par.settings = list(par.xlab.text = list(cex = 1.25),
                             par.ylab.text = list(cex = 1.25),
                             clip = list(panel = F),
                             layout.widths = list(right.padding = 4.5),
                             layout.heights = list(top.padding = 5)),
         panel = function(...) {
             panel.abline(h = seq(from = 2.5, to = 7.5, by = 2),
                          lty = 1, lwd = 1, col = gray(0))
             panel.abline(v = 0, lwd = 1.5, col = gray(0.4))
             panel.barchart(...,
                            col = c("plum1"),
                            reference = F)
             panel.segments(x0 = effects_final$effect - effects_final$sem,
                            x1 = effects_final$effect + effects_final$sem,
                            y0 = 1:nrow(effects_final),
                            y1 = 1:nrow(effects_final),
                            lty = 1, col = gray(0.2))
             panel.text(x = 2.35,
                        y = 1:8,
                        labels = c("Thr", "Gly"),
                        cex = 1.25)
             ## TFT
             panel.text(x = 2.35,
                        y = 9,
                        labels = "TFT",
                        fontface = "bold",
                        cex = 1.25)
             ## legend
             ## pvals
             panel.text(x = -0.35,
                        y = 9,
                        labels = expression(italic("p")*"-value"),
                        cex = 1.25)
             panel.text(x = c(-0.36, -0.36, -0.36, -0.39, -0.36, -0.36),
                        y = c(8, 7, 6, 4, 2, 1),
                        labels = c("8.8e-13", "2.6e-16",
                                   "4.8e-10",
                                   "1.8e-7",
                                   "2.6e-16", "8.7e-13"),
                        fontface = "plain",
                        cex = 1.1)
         })
## code end 
dev.off()


## -----
## UBC6



## -----
## setup
library("lattice")
dirs <- c("~/data/flow/2021.07.16_UBC6_Ala_1049_CRISPR_swap_flow",
          "~/data/flow/2021.07.15_UBC6_Thr_1043_CRISPR_swap_flow")

## -----
## read in allele effects and create dataframe
frame_dirs <- paste0(dirs, "/dataframes")

effect_list <- list()

for(x in 1:length(frame_dirs)){

    ## load in allele effects for ea. reporter 
    load(dir(frame_dirs[x],
             pattern = ".*effects.R",
             full.names = T),
         verbose = T)
    
    ## assign to the list
    effect_list[[x]] <- allele_effects

    ## remove object before going on to next reporter 
    rm(allele_effects)

}

effects <- do.call("rbind", effect_list)


## -----
## subset and order
effects_sub <- effects[effects$allele != "BY_full", ]

effects_sub$allele_drop <- droplevels(effects_sub$allele)

effects_sub$allele_ordered <- factor(effects_sub$allele,
                                     levels = levels(effects_sub$allele_drop)[c(2, 3, 1, 4)])

effects_final <- effects_sub[order(effects_sub$allele_ordered,
                                   decreasing = T), ]

effects_final <- effects_final[ , -7]

effects_final


## -----
## plot
pdf(file = "~/emacs/N-end_Rule_QTL_paper/figures_drafts/fine-mapping_summary_UBC6_barchart.pdf",
    height = 4.5, width = 7)
barchart(1:nrow(effects_final) ~ effects_final$effect,
         origin = 0,
         ylab = "",
         xlab = "UPS Activity Relative to BY Allele Median (SD units)",
         xlim = c(-0.26, 1.2),
         box.ratio = 2,
         scales = list(tck = c(1, 0),
                       y = list(at = seq(from = 1.5, to = 7.5, by = 2),
                                labels = c(expression(italic("UBC6")*" RM term"),
                                           expression(italic("UBC6")*" D229G"),
                                           expression(italic("UBC6")*" RM pr"),
                                           expression(italic("UBC6")*" RM")),
                                cex = 1.25),
                       x = list(cex = 1.25,
                                at = seq(from = 0, to = 1, by = 0.25),
                                labels = c("0", "0.25", "0.50", "0.75", "1.0"))),
         par.settings = list(par.xlab.text = list(cex = 1.25),
                             par.ylab.text = list(cex = 1.25),
                             clip = list(panel = F),
                             layout.widths = list(right.padding = 4.5),
                             layout.heights = list(top.padding = 5)),
         panel = function(...) {
             panel.abline(h = seq(from = 2.5, to = 7.5, by = 2),
                          lty = 1, lwd = 1, col = gray(0))
             panel.abline(v = 0, lwd = 1.5, col = gray(0.4))
             panel.barchart(...,
                            col = c("plum1"),
                            reference = F)
             panel.segments(x0 = effects_final$effect - effects_final$sem,
                            x1 = effects_final$effect + effects_final$sem,
                            y0 = 1:nrow(effects_final),
                            y1 = 1:nrow(effects_final),
                            lty = 1, col = gray(0.2))
             panel.text(x = 1.27,
                        y = 1:8,
                        labels = c("Ala", "Thr"),
                        cex = 1.25)
             ## TFT
             panel.text(x = 1.27,
                        y = 9,
                        labels = "TFT",
                        fontface = "bold",
                        cex = 1.25)
             ## legend
             ## pvals
             panel.text(x = -0.125,
                        y = 9,
                        labels = expression(italic("p")*"-value"),
                        cex = 1.25)
             panel.text(x = c(-0.15, -0.16, -0.16, -0.16, -0.175),
                        y = c(8, 7, 4, 3, 1),
                        labels = c("3.0e-18", "3.4e-8",
                                   "6.1e-9", "1.6e-9",
                                   "0.011"),
                        fontface = "plain",
                        cex = 1.1)
         })
## code end 
dev.off()


## -----
## NTA1


## -----
## setup
library("lattice")
dirs <- c("~/data/flow/2021.07.07_NTA1_Asp_0659_CRISPR_swap_flow",
          "~/data/flow/2021.07.08_NTA1_Asn_0657_CRISPR_swap_flow")

## -----
## read in allele effects and create dataframe
frame_dirs <- paste0(dirs, "/dataframes")

effect_list <- list()

for(x in 1:length(frame_dirs)){

    ## load in allele effects for ea. reporter 
    load(dir(frame_dirs[x],
             pattern = ".*effects.R",
             full.names = T),
         verbose = T)
    
    ## assign to the list
    effect_list[[x]] <- allele_effects

    ## remove object before going on to next reporter 
    rm(allele_effects)

}

effects <- do.call("rbind", effect_list)


## -----
## subset and order
effects_sub <- effects[effects$allele != "BY_full", ]

effects_sub$allele_drop <- droplevels(effects_sub$allele)

effects_sub$allele_ordered <- factor(effects_sub$allele,
                                     levels = levels(effects_sub$allele_drop)[c(3, 4, 1, 2)])

effects_final <- effects_sub[order(effects_sub$allele_ordered,
                                   decreasing = T), ]

effects_final <- effects_final[ , -7]

effects_final <- effects_final[effects_final$reporter != "Asp_TFT", ]

effects_final


## -----
## plot
pdf(file = "~/emacs/N-end_Rule_QTL_paper/figures_drafts/fine-mapping_summary_NTA1_barchart.pdf",
    height = 3.5, width = 7)
barchart(1:nrow(effects_final) ~ effects_final$effect,
         origin = 0,
         ylab = "",
         xlab = "UPS Activity Relative to BY Allele Median (SD units)",
         xlim = c(-1, 1.8),
         box.ratio = 2,
         scales = list(tck = c(1, 0),
                       y = list(at = 1:4,
                                labels = c(expression(italic("NTA1")*" E129G"),
                                           expression(italic("NTA1")*" D111E"),
                                           expression(atop(NA, atop(italic("NTA1")*" RM",
                                                      "promoter"))),
                                           expression(italic("NTA1")*" RM")),
                                cex = 1.25),
                       x = list(cex = 1.25,
                                at = seq(from = -4, to = 4, by = 0.5))),
         par.settings = list(par.xlab.text = list(cex = 1.25),
                             par.ylab.text = list(cex = 1.25),
                             clip = list(panel = F),
                             layout.widths = list(right.padding = 4.5),
                             layout.heights = list(top.padding = 5)),
         panel = function(...) {
             panel.abline(h = seq(from = 1.5, to = 3.5, by = 1),
                          lty = 1, lwd = 1, col = gray(0))
             panel.abline(v = 0, lwd = 1.5, col = gray(0.4))
             panel.barchart(...,
                            col = c("tan1"),
                            reference = F)
             panel.segments(x0 = effects_final$effect - effects_final$sem,
                            x1 = effects_final$effect + effects_final$sem,
                            y0 = 1:nrow(effects_final),
                            y1 = 1:nrow(effects_final),
                            lty = 1, col = gray(0.2))
             ## panel.text(x = 1.95,
             ##            y = 1:4,
             ##            labels = c("Asn"),
             ##            cex = 1.25)
             ## TFT
             panel.text(x = 0.5,
                        y = 4.84,
                        labels = "Asn TFT",
                        fontface = "bold",
                        cex = 1.25)
             ## legend
             ## pvals
             panel.text(x = -0.75,
                        y = 4.87,
                        labels = expression(italic("p")*"-value"),
                        cex = 1.25)
             panel.text(x = c(-0.78, -0.81, -0.78),
                        y = c(4, 2, 1),
                        labels = c("3.8e-11", "7.6e-3", "2.1e-10"),
                        fontface = "plain",
                        cex = 1.1)
         })
## code end 
dev.off()
