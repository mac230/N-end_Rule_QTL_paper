## -----
##############
## USER INPUT:
##############

## the only things that should need to be changed for this script to work
## no trailing '/' at the end!
base_dir <- "~/data/flow/2021.01.20_sort_figure"
setwd(base_dir)
needed_dirs <- c("/fcs", "/results", "/tables", "/scripts")
dir_maker <- function(x){
    ifelse(!dir.exists(paths = paste0("./", x)),
           dir.create(path = paste0("./", x)),
           paste0("dir ", paste0(getwd(), x), " exists_"))
}
sapply(X = needed_dirs, FUN = dir_maker)
work_dir       <- paste0(base_dir, "/fcs")
results_dir    <- paste0(base_dir, "/results")
tables_dir     <- paste0(base_dir, "/tables")

## -----
## <<Required_Packages>>
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
packages <- c("colorspace", "lattice", "ggvis", "dygraphs", "DescTools", "viridis")
sapply(X = packages, FUN = package_installer)
sapply(X = packages, FUN = require, character.only = TRUE)


## -----
## bioconductor
bioc_package_installer <- function(x){if(!requireNamespace(x))
                                          BiocManager::install(x, INSTALL_opts = '--no-lock')}
bioc_packages <-  c("flowCore", "flowViz", "flowUtils", "flowStats", "flowFP", "geneplotter", "ggcyto")
sapply(X = bioc_packages, FUN = bioc_package_installer)
sapply(X = bioc_packages, FUN = require, character.only = TRUE)


## -----
## required for merging flowsets into a single flowframe
source(file = "https://raw.githubusercontent.com/mac230/flow_scripts/master/set2frame.R")


## -----
## read in the data
dat <- read.flowSet(path = work_dir,
                    min.limit = 1,
                    alter.names = T)

##-----
## <<TFT_Transformation>>
## use the transform function to get the TFT/PSV parameters we want
## start by converting 0's in fluors to 1's via truncate transform
trunc.trans   <- truncateTransform("Convert 0's to 1's.", a = 1)
trunc.fluors  <- function(x){
    transform(x,
              `GFP.A` = trunc.trans(`GFP.A`),
              `mCherry.A` = trunc.trans(`mCherry.A`))}
dat <- fsApply(x = dat, FUN = trunc.fluors)

PSV.TFT.transform <- function(x){
    transform(x,
              `log_GFP` = log10(`GFP.A`),
              `log_RFP` = log10(`mCherry.A`),
              `TFT_ratio` = log(`mCherry.A`/`GFP.A`, base = 2),
              `PSV_ratio` = log(`GFP.A`/`mCherry.A`, base = 2),
              ## 'no log' TFT ratio
              `nl_TFT_ratio` = (`mCherry.A`/`GFP.A`)
              )}

dat <- fsApply(x = dat, FUN = PSV.TFT.transform)

trunc.gfp <- truncateTransform("make set GFP positive", a = 2.5)
trunc.gfp.f <-  function(x){
    transform(x,
              `log_GFP` = trunc.gfp(`log_GFP`))}

dat <- fsApply(x = dat, FUN = trunc.gfp.f)

## -----
## <<FSC_Gate>>
## a function to gate the cells to include only haploids.
## we identify these as a sharp peak in the lower end of
## the fsc density plot.  I take 10% above and below the
## max density value
fsc.gate.generator <- function(x){
    fsc.dens  <- density(exprs(x[, 1]))
    ## return the index of the maximum y value of the density estimate
    fsc.max   <- fsc.dens[[1]][which.max(fsc.dens[[2]])]
    fsc.upper <- (fsc.max * 0.25) + fsc.max
    fsc.lower <- fsc.max - (fsc.max * 0.25)
    fsc.gate  <- c(fsc.lower, fsc.upper)
}

curv.split <- function(x){
    split(x, f = rectangleGate("FSC.A" = fsc.gate.generator(x)),
          population = "defaultRectangleGate+",
          flowSet = T, codeflowSet = T)}
curv.set <- fsApply(x = dat, FUN = curv.split)


tft.gate.generator <- function(x){
    tft.dens  <- density(exprs(x[, 15]))
    ## return the index of the maximum y value of the density estimate
    tft.max   <- tft.dens[[1]][which.max(tft.dens[[2]])]
    tft.upper <- (tft.max * 0.025) + tft.max
    tft.lower <- tft.max - (tft.max * 0.025)
    tft.gate  <- c(tft.lower, tft.upper)
}

tft.gate.generator(curv.set[[1]][[1]])

dat_f <- as.data.frame(exprs(dat[[1]]))
high_FSC <- quantile(x = dat_f$FSC.A, probs = 0.70)
low_FSC <- quantile(x = dat_f$FSC.A, probs = 0.30)

dat_g <- dat_f[dat_f$FSC.A > low_FSC & dat_f$FSC.A < high_FSC, ]

## 'high_TFT' = high UPS activity, so low TFT ratio
low_col <- "#ec6ab7ff"
high_col <-  "#5AC168ff"

low_t_col <- "#ec6ab722"
high_t_col <-  "#5AC16811"

dat_g <- dat_g[dat_g$log_GFP > 2.6, ]

high_TFT_val <- quantile(x = dat_g$TFT_ratio, probs = 0.02)
low_TFT_val  <- quantile(x = dat_g$TFT_ratio, probs = 0.98)


{
pdf(file = "~/emacs/N-end_Rule_QTL_paper/figures_drafts/TFT_sort_density_plot_gate_shading.pdf")
print(
densityplot(dat_g$TFT_ratio,
            plot.points = F,
            xlab = expression("log"["2"]*" RFP / GFP"),
            key = list(lines = list(col = c(low_col,
                                            high_col),
                                    lty = 2,
                                    lwd = 4),
                       text = list(labels = c("2% Low UPS activity gate",
                                              "2% High UPS activity gate"),
                                   cex = 2)),
            scales = list(tck = c(1, 0),
                          x = list(cex = 2,
                                   at = -4:1),
                          y = list(cex = 2)),
            xlim = c(-5, 2),
            col = gray(0.1),
            par.settings = list(par.xlab.text = list(cex = 2),
                                par.ylab.text = list(cex = 2)),
            lwd = 2,
            panel = function(...) {
                panel.densityplot(...)
                panel.abline(v = c(high_TFT_val,
                                   low_TFT_val),
                             col = c(high_col,
                                     low_col),
                             lwd = 4, lty = 2)
                panel.abline(h = 0,
                             col = gray(0.9),
                             lwd = 2)
})
)
## comment this out to get rid of the shading
grid.rect(x = c(0.27, 0.850),
          y = c(0.505, 0.505),
          width = c(0.22, 0.21),
          height = c(0.707, 0.707),
          default.units = "npc",
          draw = T,
          gp = gpar(col = c(low_t_col, high_t_col),
                    fill = c(low_t_col, high_t_col))
)
dev.off()
}

high_TFT <- dat_g[dat_g$TFT_ratio <= high_TFT_val, ]
low_TFT <- dat_g[dat_g$TFT_ratio >= low_TFT_val, ]


## -----
## plot
## 2021.10.06 - update to no longer plot neg. ctrl.
{
pdf("~/emacs/N-end_Rule_QTL_paper/figures_drafts/TFT_sort_xy_plot.pdf")
print(
xyplot(log_RFP ~ log_GFP,
       data = dat_g,
       type = c("p", "g"),
       xlim = c(1.5, 4.5),
       ylim = c(1.5, 4.5),
       pch = 19,
       col = gray(0.3, alpha = 0.2),
       cex = 0.25,
       xlab = expression("log"["10"]*" GFP"),
       ylab = expression("log"["10"]*" RFP"),
       key = list(lines = list(col = c(low_col,
                                       high_col),
                               lwd = 10,
                               size = 2.5),
                  text = list(labels = c("2% Low UPS Activity Gate",
                                         "2% High UPS Activity Gate"),
                                         cex = 2),
                                         corner = c(0.02, 0.98)),
       scales = list(tck = c(1, 0),
                     alternating = F,
                     x = list(cex = 2),
                     y = list(cex = 2)),
       par.settings = list(par.ylab.text = list(cex = 2),
                           par.xlab.text = list(cex = 2)),
       panel = function(...) {
           panel.xyplot(...)
           panel.points(x = high_TFT$log_GFP, y = high_TFT$log_RFP,
                        pch = 19, col = high_col, cex = 0.25)
           panel.points(x = low_TFT$log_GFP, y = low_TFT$log_RFP,
                        pch = 19, col = low_col, cex = 0.25)
}
))
grid.text(label = c("20,000 cells",
                    "each collected",
                    "from high and low",
                    "UPS activity gates"),
          x = c(0.31, 0.335, 0.37, 0.375),
          y = c(0.80, 0.75, 0.70, 0.65),
          rot = 0,
          default.units = "npc",
          gp = gpar(col = "black",
                    cex = 2,
                    just = "left"))
dev.off()
}
