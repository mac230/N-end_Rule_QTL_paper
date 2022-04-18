## -----
## setup
library("lattice")
dirs <- c("~/data/flow/2021.02.12_UBR1_CRISPR_swap_0657_TFT_flow",
          "~/data/flow/2021.02.06_UBR1_CRISPR_swap_0659_TFT_flow",
          "~/data/flow/2021.02.11_UBR1_CRISPR_swap_0665_TFT_flow",
          "~/data/flow/2021.02.10_UBR1_CRISPR_swap_0663_TFT_flow")


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
                                     levels = levels(effects_sub$allele_drop)[c(1, 3, 2, 4)])

effects_final <- effects_sub[order(effects_sub$allele_ordered,
                                   decreasing = T), ]

effects_final <- effects_final[ , -7]

effects_final


## -----
## plot
pdf(file = "~/emacs/N-end_Rule_QTL_paper/figures_drafts/fine-mapping_summary_full_gene_barchart.pdf")
barchart(1:nrow(effects_final) ~ effects_final$effect,
         origin = 0,
         ylab = "",
         xlab = "UPS Activity Relative to BY Allele Median (SD units)",
         xlim = c(-2.4, 2.2),
         box.ratio = 2,
         scales = list(tck = c(1, 0),
                       y = list(at = seq(from = 2.5, to = 14.5, by = 4),
                                labels = c(expression(italic("UBR1")*" RM term"),
                                           expression(italic("UBR1")*" RM ORF"),
                                           expression(italic("UBR1")*" RM pr"),
                                           expression(italic("UBR1")*" RM")),
                                cex = 1.25),
                       x = list(cex = 1.25,
                                at = seq(from = -4.5, to = 4.5, by = 0.5))),
         par.settings = list(par.xlab.text = list(cex = 1.25),
                             par.ylab.text = list(cex = 1.25),
                             clip = list(panel = F),
                             layout.widths = list(right.padding = 4.5),
                             layout.heights = list(top.padding = 5)),
         panel = function(...) {
             panel.abline(h = seq(from = 4.5, to = 12.5, by = 4),
                          lty = 1, lwd = 1, col = gray(0))
             panel.abline(v = 0, lwd = 1.5, col = gray(0.4))
             panel.barchart(...,
                            col = c("tan1", "tan1", "tomato", "tomato"),
                            reference = F)
             panel.segments(x0 = effects_final$effect - effects_final$sem,
                            x1 = effects_final$effect + effects_final$sem,
                            y0 = 1:nrow(effects_final),
                            y1 = 1:nrow(effects_final),
                            lty = 1, col = gray(0.2))
             panel.text(x = 2.45,
                        y = 1:16,
                        labels = c("Asn", "Asp", "Phe", "Trp"),
                        cex = 1.25)
             ## TFT
             panel.text(x = 2.45,
                        y = 17,
                        labels = "TFT",
                        fontface = "bold",
                        cex = 1.25)
             ## legend
             panel.rect(xleft = c(-2.3, 0.1),
                        xright = c(-2, 0.4),
                        ytop = c(18.5, 18.5),
                        ybottom = c(17.8, 17.8),
                        col = c("tan1", "tomato"))
             panel.text(x = c(-1.15, 1.27),
                        y = c(18.15, 18.15),
                        labels = c("Type I Arg/N-end",
                                   "Type II Arg/N-end"),
                        fontface = "plain",
                        cex = 1.25)
             ## pvals
             panel.text(x = -2,
                        y = 17,
                        labels = expression(italic("p")*"-value"), 
                        cex = 1.25)
             panel.text(x = c(-2.05, -2.05, -2.1, -2.1, -2.05, -2.1, -2.05, -2.1, -2.15),
                        y = c(16:11, 8, 7, 5),
                        labels = c("1.8e-19", "5.9e-14",  "4.6e-4", "7.7e-5",
                                   "5.8e-16", "9.4e-4",
                                   "9.9e-12", "5.0e-4", "0.028"),
                        fontface = "plain",
                        cex = 1.1)
         })
## code end 
dev.off()

## -----
## single variant fine-mapping 


## -----
## setup
library("lattice")
dirs <- c("~/data/flow/2021.03.19_UBR1_single_variant_0665_TFT_freezer_fresh",
          "~/data/flow/2021.03.19_UBR1_single_variant_0663_TFT_freezer_fresh")


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
effects_sub <- effects_sub[effects_sub$allele != "RM_full", ]

effects_sub$allele_drop <- droplevels(effects_sub$allele)

effects_sub$allele_ordered <- factor(effects_sub$allele,
                                     levels = levels(effects_sub$allele_drop)[c(3, 1, 2)])

effects_final <- effects_sub[order(effects_sub$allele_ordered,
                                   decreasing = T), ]

effects_final <- effects_final[ , -7]

## reporter ordering 
## effects_final <- effects_final[c(2, 1, 4, 3, 6, 5, 8, 7), ]

## -----
## plot
pdf(file = "~/emacs/N-end_Rule_QTL_paper/figures_drafts/fine-mapping_summary_single_variant_barchart_final.pdf",
    height = 4.5, width = 7)
barchart(1:nrow(effects_final) ~ effects_final$effect,
         origin = 0,
         ylab = "",
         xlab = "UPS Activity Relative to BY Allele Median (SD units)",
         xlim = c(-1.6, 1.6),
         box.ratio = 2,
         scales = list(tck = c(1, 0),
                       y = list(at = seq(from = 1.5, to = 5.5, by = 2),
                                labels = c(expression(italic("UBR1")*" RM -197A>C"),
                                           expression(italic("UBR1")*" RM -469A>T"),
                                           expression(italic("UBR1")*" RM pr")),
                                cex = 1.25),
                       x = list(cex = 1.25,
                                at = seq(from = -4.5, to = 4.5, by = 0.5))),
         par.settings = list(par.xlab.text = list(cex = 1.25),
                             par.ylab.text = list(cex = 1.25),
                             clip = list(panel = F),
                             layout.widths = list(right.padding = 4.5),
                             layout.heights = list(top.padding = 5)),
         panel = function(...) {
             panel.abline(h = seq(from = 2.5, to = 4.5, by = 2),
                          lty = 1, lwd = 1, col = gray(0))
             panel.abline(v = 0, lwd = 1.5, col = gray(0.4))
             panel.barchart(...,
                            col = "tomato",
                            reference = F)
             panel.segments(x0 = effects_final$effect - effects_final$sem,
                            x1 = effects_final$effect + effects_final$sem,
                            y0 = 1:nrow(effects_final),
                            y1 = 1:nrow(effects_final),
                            lty = 1, col = gray(0.2))
             panel.text(x = 1.8,
                        y = 1:6,
                        labels = c("Phe", "Trp"),
                        cex = 1.25)
             ## TFT
             panel.text(x = 1.8,
                        y = 6.87,
                        labels = "TFT",
                        fontface = "bold",
                        cex = 1.25)
             ## ## legend
             ## panel.rect(xleft = c(-2.2, 0.1),
             ##            xright = c(-1.9, 0.4),
             ##            ytop = c(18.5, 18.5),
             ##            ybottom = c(17.8, 17.8),
             ##            col = c("tan1", "tomato"))
             ## panel.text(x = c(-1.05, 1.27),
             ##            y = c(18.15, 18.15),
             ##            labels = c("Type I Arg/N-end",
             ##                       "Type II Arg/N-end"),
             ##            fontface = "plain",
             ##            cex = 1.25)
             ## pvals
             panel.text(x = -1.27,
                        y = 6.87,
                        labels = expression(italic("p")*"-value"),
                        cex = 1.25)
             panel.text(x = c(rep(-1.3, 4)),
                        y = c(6:3),
                        labels = c("5.0e-19", "3.6e-11",
                                   "3.4e-19", "5.8e-16"),
                        fontface = "plain",
                        cex = 1.1)
         })
## code end 
dev.off()
