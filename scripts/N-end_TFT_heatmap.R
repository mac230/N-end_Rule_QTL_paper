## -----
## load all the required packages
source("~/emacs/R/functions/load_flow_packages.R")

#############
## USER INPUT
#############
reporter_names <- c("Gln_TFT", "Asn_TFT", "Glu_TFT", "Asp_TFT",
                    "Ile_TFT", "His_TFT", "Tyr_TFT", "Trp_TFT",
                    "Leu_TFT", "Phe_TFT", "Lys_TFT", "Arg_TFT",
                    "Thr_TFT", "Met_TFT", "Pro_TFT", "Val_TFT",
                    "Gly_TFT", "Ser_TFT", "Ala_TFT", "Cys_TFT")
base_dir       <- "~/data/flow/N-end_all_reporters/"
work_dir       <- paste0(base_dir, "fcs/")
frame_dir      <- paste0(base_dir, "dataframes/")
gated_dir      <- paste0(frame_dir, "gated/")
ungated_dir    <- paste0(frame_dir, "ungated/")
#################
## END USER INPUT
#################


## -----
## we'll loop over all the reporter names to
## [1] read in all the fcs files for a given reporter
## [2] gate each file based on fsc
## [3] write each collection of files to a data table
## [4] add variables for downstream processing (strain/reporter names, e.g.)
for (k in 1:length(reporter_names)) {

## -----
## [x]
## name the strains based on reporter, then set
## regex for getting flowsets of the different strains
## generally, should name fcs files as follows:
## strain    - by, rm, rpn4, ...
## reporter  - Ala TFT, Arg TFT, ....
## replicate - 001, 002, etc... per strain
## e.g., "RM_Arg_TFT_001.fcs"

## this is my current approach to reading in files.  the idea is to
## put all the various strains I've used in the past here and filter this
## complete set to those strains present in the actual data I'm analyzing.  I
## pre-filter using 'grepl' because 'read.flowset' throws an error if any term
## you supply it doesn't match.
no_reporter  <- paste0(".*untagged.*", reporter_names[k], ".*fcs") 
by_strain    <- paste0("BY.*", reporter_names[k], ".*fcs") 
rm_strain    <- paste0("RM.*", reporter_names[k], ".*fcs") 
rpn4_strain  <- paste0("rpn4.*", reporter_names[k], ".*fcs") 
ubr1_strain  <- paste0("ubr1.*", reporter_names[k], ".*fcs") 
doa10_strain <- paste0("doa10.*", reporter_names[k], ".*fcs") 
pop_1_strain <- paste0("SFA.*pop_001.*", reporter_names[k], ".*fcs") 
pop_5_strain <- paste0("SFA.*pop_005.*", reporter_names[k], ".*fcs") 
pop_6_strain <- paste0("SFA.*pop_006.*", reporter_names[k], ".*fcs") 

all_strains <- c(no_reporter, by_strain, rm_strain,
                 rpn4_strain, ubr1_strain, doa10_strain,
                 pop_1_strain, pop_5_strain, pop_6_strain)

## 'dir' lists the contents of a directory - test whether each strain regex
## matches any files in the list produced by 'dir'.  if a match occurs, the
## value returned by 'max' will be 1, else 0.  take only the strains that
## returned a match in the 2nd step below w/ 'all.strains <- ...'
true_strains <- sapply(all_strains, function(strain){
                           as.logical(max(grepl(pattern = strain,
                                                x = dir(path = work_dir,
                                                        pattern = ".*.fcs",
                                                        include.dirs = F,
                                                        ignore.case = T,
                                                        recursive = F,
                                                        no.. = T))))
                       })

all_strains <- as.list(all_strains[true_strains])

all_set     <- lapply(all_strains, function(strain){
                          read.flowSet(files = NULL,
                                       path = work_dir,
                                       pattern = strain,
                                       alter.names = T,
                                       min.limit = 1)
                      })
## str(all_set[[1]]@phenoData@data$name)

name_list <- strsplit(x = names(true_strains)[true_strains == T],
                      split = "\\.\\*")

names(all_set) <- unlist(lapply(X = name_list, FUN = function(x) {
                                    x[1] }))

## -----
## <<Color_Setup>>
## linking colors to strain names in R
## I think I should be able to make something
## akin to an lisp association list where
## there is a strain name and associated color
col_untagged <- c(color = gray(0.7),   name = "no reporter")
## 2021.01.27 - new colors for BY/RM
## col_by       <- c(color = "#7A9BCCFF", name = ".*BY.*")
col_by       <- c(color = "#2166ACFF", name = ".*BY.*")
## col_rm       <- c(color = "#CC7AAAFF", name = ".*RM.*")
col_rm       <- c(color = "#BF3232FF", name = ".*RM.*")
## col_rpn4     <- c(color = "#CCAB7AFF", name = ".*rpn4.*")
col_rpn4     <- c(color = "#CCA14AFF", name = ".*rpn4.*")
## col_ubr1     <- c(color = "#88CCBBFF", name = ".*ubr1.*")
col_ubr1     <- c(color = "#78BCABFF", name = ".*ubr1.*")
## col_doa10    <- c(color = "#A3CC7AFF", name = ".*doa10.*")
col_doa10    <- c(color = "#83BC6AFF", name = ".*doa10.*")
col_pop_1    <- c(color = gray(0.7),   name = ".*population.*1.*")
col_pop_5    <- c(color = "#AA1111FF", name = ".*population.*5.*")
col_pop_6    <- c(color = gray(0),     name = ".*population.*6.*")

cols_list    <- list(col_untagged, col_by, col_rm, 
                     col_rpn4, col_ubr1, col_doa10,
                     col_pop_1, col_pop_5, col_pop_6)

col_out <- sapply(X = cols_list, FUN = function(x){
                      grepl(pattern = x["name"],
                            x = name_list)
                  })

col_out <- as.logical(unlist(sapply(1:ncol(col_out),
                                    FUN = function(x){
                                        max(col_out[, x])
                                    })))

all_cols <- unlist(sapply(X = cols_list[col_out],
                          FUN = function(x){identity(x["color"])}))

names(all_cols) <- names(all_set)

## output a dummmy plot to assess strain/color mapping
## setwd(results.dir)
## pdf(file = "color_mapping.pdf", height = 7, width = 7, bg = "transparent")
barplot(rep(4, length(name_list)), col = all_cols, ylim = c(0, 7))
box()
legend(x = "topleft", legend = names(all_set),
       lty = 1, lwd = 7.5, col = all_cols,
       bg = "white")
legend(x = "topright", y = NA,
       legend = unlist(lapply(X = cols_list, FUN = function(x){identity(x)["name"]})),
       col = unlist(lapply(X = cols_list, FUN = function(x){identity(x)["color"]})),
       lty = 1, lwd = 7.5,  bg = "white")
## dev.off()


## a function to gate the cells to include only haploids.
## we identify these as a sharp peak in the lower end of
## the fsc density plot.  I take 10% above and below the
## max density value
fsc_gate_generator <- function(fl_frame){
    fsc_dens  <- density(exprs(fl_frame[, 1]))
    ## return the index of the maximum y value of the density estimate
    fsc_max   <- fsc_dens[[1]][which.max(fsc_dens[[2]])]
    fsc_upper <- (fsc_max * 0.10) + fsc_max
    fsc_lower <- fsc_max - (fsc_max * 0.10)
    fsc_gate  <- c(fsc_lower, fsc_upper)
}

fsc_split <- function(x){
    split(x, f = rectangleGate("FSC.A" = fsc_gate_generator(x)),
          population = "defaultRectangleGate+",
          flowSet = T, codeflowSet = T)}

## gate all samples on FSC
fsc_set <- lapply(all_set, fsApply, fsc_split)

## convert gated samples from flowsets to flowframes 
fsc_frame <- lapply(fsc_set, function(set) {
                        lapply(set, set2Frame)
                    })

## pull dataframes from flowframes for gated and
## ungated sets of samples, then rename 
fsc_gate_exprs <- vector(mode = "list", length = length(all_set))
no_gate_exprs  <- fsc_gate_exprs

for(j in 1:length(fsc_frame)) {
    for (i in 1:length(fsc_frame[[j]])) {
        fsc_gate_exprs[[j]][[i]] <- as.data.frame(exprs(fsc_frame[[j]][[i]]))
        fsc_gate_exprs[[j]][[i]]$strain <- as.factor(names(fsc_frame[j]))
        fsc_gate_exprs[[j]][[i]]$replicate <- as.factor(i)
    }
    fsc_gate_exprs[[j]] <- do.call("rbind", fsc_gate_exprs[[j]])
} 

## bind into a single dataframe
fsc_gate_exprs <- do.call("rbind", fsc_gate_exprs)

## add transformation parameters
fsc_gate_exprs$log_GFP   <- log10(fsc_gate_exprs$eGFP.A)
fsc_gate_exprs$log_RFP   <- log10(fsc_gate_exprs$mCherry.A)
fsc_gate_exprs$TFT_ratio <- log2(fsc_gate_exprs$mCherry.A / fsc_gate_exprs$eGFP.A)
fsc_gate_exprs$reporter  <- as.factor(rep(x = gsub(pattern = "_",
                                                  replacement = " ",
                                                  x = reporter_names[k]),
                                         times = nrow(fsc_gate_exprs)))

## nested loops for the ungated data
for (j in 1:length(all_set)){
    for (i in 1:length(all_set[[j]])) {
        no_gate_exprs[[j]][[i]] <- as.data.frame(exprs(all_set[[j]][[i]]))
        no_gate_exprs[[j]][[i]]$strain <- as.factor(names(all_set[j]))
        no_gate_exprs[[j]][[i]]$replicate <- as.factor(i)
    }
    no_gate_exprs[[j]] <- do.call("rbind", no_gate_exprs[[j]])
}

## bind into a single dataframe
no_gate_exprs <- do.call("rbind", no_gate_exprs)

## add transformation parameters and reporter var
no_gate_exprs$log_GFP   <- log10(no_gate_exprs$eGFP.A)
no_gate_exprs$log_RFP   <- log10(no_gate_exprs$mCherry.A)
no_gate_exprs$TFT_ratio <- log2(no_gate_exprs$mCherry.A / no_gate_exprs$eGFP.A)
no_gate_exprs$reporter  <- as.factor(rep(x = gsub(pattern = "_",
                                                  replacement = " ",
                                                  x = reporter_names[k]),
                                         times = nrow(no_gate_exprs)))

## write the ungated data to the appropriate dir
write.table(x = no_gate_exprs,
            file = paste0(ungated_dir,
                          reporter_names[k],
                          "_all_ungated.csv"),
            append = F, sep = ",",
            quote = F, row.names = F)

## write the gated data to the appropriate dir
write.table(x = fsc_gate_exprs,
            file = paste0(gated_dir,
                          reporter_names[k],
                          "_all_gated.csv"),
            append = F, sep = ",",
            quote = F, row.names = F)

    }

## now, read each reporter's dataframe in and
## combine into a single dataframe
## generate a list of files in a directory using
## the 'dir' command, e.g.:
## dir(gated_dir)
## dir(ungated_dir)

## '_u' = ungated
out_u <- vector(mode = "list", length = length(dir(ungated_dir)))
for (o in 1:length(dir(gated_dir))) {
    out_u[[o]] <- read.table(file = paste0(ungated_dir, dir(ungated_dir)[o]),
                             header = T, sep = ",")
       }

## '_g' = gated
out_g <- vector(mode = "list", length = length(dir(gated_dir)))
for (o in 1:length(dir(gated_dir))) {
    out_g[[o]] <- read.table(file = paste0(gated_dir, dir(gated_dir)[o]),
                             header = T, sep = ",")
       }

## for the ungated set, it'll be
## 1e5 cells * 5 strains/reporter * 20 reporters = 1e7 rows
## fsc gating reduces 1e5 to ~2e4, so ~2e6 rows
## nrow(out_all) = 2284942
out_all <- do.call("rbind", out_g)
str(out_all)
levels(out_all$reporter)

## split by type I/II Arg/N-end 
aa_order <- c(2, 3, 4, 6, 7, 9, 12, 10, 11, 14, 18, 19, 1, 5, 8, 13, 15, 16, 17, 20)

## get the strain factor in the desired order 
out_all$strain <- factor(out_all$strain,
                         levels = levels(out_all$strain)[c(1, 3, 4, 5, 2)])

strain_paste <- expand.grid(unique(out_all$replicate),
                          levels(out_all$strain))

strain_paste <- paste0(strain_paste$Var2, "_", strain_paste$Var1)

out_all$strain_rep <- factor(paste0(out_all$strain, "_", out_all$replicate),
                             levels = strain_paste)

rep_cols <- unlist(lapply(X = 1:length(all_cols), FUN = function(x) {
                       rep(all_cols[x],
                           times = sum(grepl(pattern = names(all_cols[x]), 
                                             x = levels(out_all$strain_rep))))
                   }))

## need to order levels of 'strain_rep' like 'strain'
## 2021.01.20 - split by type I/II Arg/N-end 
aa_order <- c(2, 3, 4, 6, 7, 9, 12, 10, 11, 14, 18, 19, 1, 5, 8, 13, 15, 16, 17, 20)
out_all$o_reporter <- factor(out_all$reporter,
                             levels = levels(out_all$reporter)[aa_order])
levels(out_all$o_reporter)
params <- colnames(out_all)[unlist(lapply(X = out_all, FUN = is.numeric))]
params[10] <- "log2 TFT Ratio"

## <<session_restore_point>>
## from this point, I've:
## 1. loaded all N-end TFT flowsets into R
## 2. assigned colors
## 3. gated the cells on the basis of FSC
## 4. build dataframes from the flowsets
## 5. *ordered the reporters by type and alphabetical*
## save.session(file = "~/Desktop/2021.01.29_all_flow_sets_loaded")
source("~/emacs/R/functions/load_flow_packages.R")
restore.session(file = "~/Desktop/2021.01.29_all_flow_sets_loaded")
## ls()

## <<density_plot_final>>
## have to set up parameters for a custom legend first...
## drop everything for building a custom legend
## into a list.  loop over the list by position
## ('top' or 'bottom') and strain to make 2
## legends in the panels 
top_legend_params <- list()
top_legend_params$names <- c("BY", "RM",
                         expression(paste("BY rpn4", Delta)), 
                         expression(paste("BY ubr1", Delta)), 
                         expression(paste("BY doa10", Delta)))
top_legend_params$line_x1     <- rep(0.8, 5)
top_legend_params$line_x2     <- rep(1.15, 5)    
top_legend_params$x_positions <- rep(x = 1.2, times = 5)
top_legend_params$y_positions <- rev(seq(from = 12.7, to = 13.45, length.out = 5))
top_legend_params$color       <- all_cols
bot_legend_params <- top_legend_params
top_legend_params$y_positions <- rev(seq(from = 4.75, to = 5.5, length.out = 5))
legend_params <- list(top_legend_params, bot_legend_params)

{
pdf(file = "~/Desktop/2022.01.29_density_full_final.pdf", height = 14.5, width = 14)
print(

densityplot(~ TFT_ratio | reporter,
            groups = strain_rep,
            data = out_all,
            xlim = c(-7, 1.5),
            ## set alternating = F for one side, same side labeling
            Scales = list(alternating = 3),
            grid = T,
            plot.points = F,
            lwd = 2,
            main = list(label = "Arg/N-end Reporters"),
            sub = list(label = "Ac/N-end Reporters"),
            between = list(x = c(0, 0, 0),
                           y = c(0, 0, 3)),
            as.table = T,
            ylab = ",,",
            xlab = gsub(pattern = "_",
                        replacement = " ",
                        params[10]),
            index.cond = list(aa_order),
            par.settings = list(strip.background = list(col = gray(0.9)),
                                clip = list(panel = FALSE),
                                par.main.text = list(font = 2,
                                                     cex = 1.25,
                                                     just = "center", 
                                                     x = grid::unit(7, "in")),
                                par.sub.text = list(font = 2,
                                                    just = "center",
                                                    cex = 1.25,
                                                    x = grid::unit(7.05, "in"),
                                                    y = grid::unit(5.95, "in")),
                                axis.text = list(cex = 1),
                                par.ylab.text = list(cex = 1.25,
                                                     col = "white"),
                                par.xlab.text = list(cex = 1.25)),
            ## legend = list(inside = list(fun = grid.legend,
            ##                             args = list(labels = c("BY", "RM",
            ##                                   expression(paste("BY rpn4", Delta)), 
            ##                                   expression(paste("BY ubr1", Delta)), 
            ##                                   expression(paste("BY doa10", Delta))),
            ##                                   do.lines = T,
            ##                                   nrow = 5,
            ##                                   draw = T,
            ##                                   hgap = 1,
            ##                                   vgap = 0.25,
            ##                                   gp = gpar(col = all_cols,
            ##                                             lwd = 5,
            ##                                             cex = 1,
            ##                                             lineend = "butt",
            ##                                             npc = 50
            ##                                                   )))),
            ## key = list(text = list(c("BY", "RM",
            ##                           expression(paste("BY rpn4", Delta)),
            ##                           expression(paste("BY ubr1", Delta)),
            ##                           expression(paste("BY doa10", Delta)))),
            ##             lines = list(col = all_cols,
            ##                          lwd = 5),
            ##            corner = c(0, 1),
            ##            y = 0.98),
            panel = function(x, y, q, subscripts, ...) {
                panel.grid(h = -1, v = -1)
                panel.densityplot(x,
                                  plot.points = F,
                                  groups = out_all$strain_rep,
                                  subscripts = subscripts,
                                  as.table = T,
                                  lty = 1,
                                  col = rep_cols,
                                  lwd = 0.75)
                panel.densityplot(x,
                                  plot.points = F,
                                  groups = out_all$strain,
                                  subscripts = subscripts,
                                  as.table = t,
                                  lty = 1,
                                  lwd = 1.5,
                                  col = all_cols,
                                  ylim = c(0, 2))
            })
)


## <<plot_legend>>
## this gets placed outside the 'print' call
## write out a legend using 'grid.text' and 'grid.lines'
## because I separate the 2 N-end pathways, doubling
## the legend makes for a better visual layout.
## turns out there's not a straightforward way to
## double or position legends built using the lattice
## 'key' argument or grid.legend.  the code below works,
## but it's a bit hacky...
for (k in 1:length(legend_params)) {
    for (l in 1:length(legend_params[[1]])) {
        grid.text(label = legend_params[[k]]$names[l],
                  x = legend_params[[k]]$x_positions[l],
                  y = legend_params[[k]]$y_positions[l],
                  default.units = "in",
                  just = "left",
                  gp = gpar(col = "black", cex = 1))
        grid.lines(x = c(legend_params[[k]]$line_x1[l],
                         legend_params[[k]]$line_x2[l]),
                   y = c(legend_params[[k]]$y_positions[l],
                         legend_params[[k]]$y_positions[l]),
                   default.units = "in",
                   gp = gpar(lwd = 4, lineend = "butt",
                             col = legend_params[[k]]$color[l]))
    }
}

## replace improper y label positioning w/ custom text
## and double it since we've split the plot into 2 panels
grid.text(label = "Density",
          x = c(0.15, 0.15),
          y = c(7.5, 2),
          rot = 90,
          default.units = "in",
          gp = gpar(col = "black", cex = 1.25)
          )

dev.off()
}

## <<example_density_plot_trp_met_final>>
## grab a set of reporters that illustrate
## the deletion phenotypes.  for now, we'll
## use Met and Trp TFTs for this purpose.
out_tw <- out_all[out_all$reporter == "Met TFT" | out_all$reporter == "Trp TFT", ]

{
pdf(file = "~/Desktop/2022.01.29_example_density_2_panel_new_colors.pdf", height = 8, width = 5)
print(
densityplot(~ TFT_ratio | reporter,
            groups = strain_rep,
            data = out_tw,
            xlim = c(-8, 1.5),
            ## set alternating = F for one side, same side labeling
            scales = list(alternating = 1,
                          tck = c(1, 0)),
            grid = T,
            plot.points = F,
            lwd = 2,
            as.table = F,
            ylab = ",",
            xlab = "log2 TFT ratio",
            par.settings = list(strip.background = list(col = gray(0.9)),
                                clip = list(panel = FALSE),
                                par.main.text = list(font = 2,
                                                     cex = 1.25,
                                                     just = "center", 
                                                     x = grid::unit(7, "in")),
                                par.sub.text = list(font = 2,
                                                    just = "center",
                                                    cex = 1.25,
                                                    x = grid::unit(7.05, "in"),
                                                    y = grid::unit(5.95, "in")),
                                axis.text = list(cex = 1),
                                par.ylab.text = list(cex = 1.25,
                                                     col = "white"),
                                par.xlab.text = list(cex = 1.25)),
            panel = function(x, y, q, subscripts, ...) {
                panel.grid(h = -1, v = -1)
                panel.densityplot(x,
                                  plot.points = F,
                                  groups = out_tw$strain_rep,
                                  subscripts = subscripts,
                                  as.table = T,
                                  lty = 1,
                                  col = rep_cols,
                                  lwd = 0.75)
                panel.densityplot(x,
                                  plot.points = F,
                                  groups = out_tw$strain,
                                  subscripts = subscripts,
                                  as.table = t,
                                  lty = 1,
                                  lwd = 2,
                                  col = all_cols,
                                  ylim = c(0, 2))
            })
)
## this gets placed outside the 'print' call
## write out a legend using 'grid.text' and 'grid.lines'
## because I separate the 2 N-end pathways, doubling
## the legend makes for a better visual layout.
## turns out there's not a straightforward way to
## double or position legends built using the lattice
## 'key' argument or grid.legend.  the code below works,
## but it's a bit hacky...
## legend parameters
## drop everything for building a custom legend
## into a list.  loop over the list by position
## ('top' or 'bottom') and strain to make 2
## legends in the panels 
ex_top_legend_params <- list()
ex_top_legend_params$names <- c("BY", "RM",
                                expression(paste("BY rpn4", Delta)), 
                                expression(paste("BY ubr1", Delta)), 
                                expression(paste("BY doa10", Delta)))
## starting point of the lines in 'x'
ex_top_legend_params$line_x1     <- rep(0.9, 5)
## endind position of the lines in 'x'
ex_top_legend_params$line_x2     <- rep(1.25, 5)    
## position of the text in x - rep this 5 times for 5 strains
ex_top_legend_params$x_positions <- rep(x = 1.35, times = 5)
## position of the text in x - rep this 5 times for 5 strains
ex_top_legend_params$y_positions <- rev(seq(from = 6.55, to = 7.25, length.out = 5))
ex_top_legend_params$color       <- all_cols
ex_bot_legend_params <- ex_top_legend_params
ex_top_legend_params$y_positions <- rev(seq(from = 3.1, to = 3.85, length.out = 5))
ex_legend_params <- list(ex_top_legend_params, ex_bot_legend_params)
for (k in 1:length(ex_legend_params)) {
    for (l in 1:length(ex_legend_params[[1]])) {
        grid.text(label = ex_legend_params[[k]]$names[l],
                  x = ex_legend_params[[k]]$x_positions[l],
                  y = ex_legend_params[[k]]$y_positions[l],
                  default.units = "in",
                  just = "left",
                  gp = gpar(col = "black", cex = 1))
        grid.lines(x = c(ex_legend_params[[k]]$line_x1[l],
                         ex_legend_params[[k]]$line_x2[l]),
                   y = c(ex_legend_params[[k]]$y_positions[l],
                         ex_legend_params[[k]]$y_positions[l]),
                   default.units = "in",
                   gp = gpar(lwd = 4, lineend = "butt",
                             col = ex_legend_params[[k]]$color[l]))
    }
}
## replace improper y label positioning w/ custom text
## and double it since we've split the plot into 2 panels
grid.text(label = "Density",
          x = c(0.2, 0.2),
          y = c(2.3, 5.7),
          rot = 90,
          gp = gpar(cex = 1.25),
          default.units = "in")
dev.off()
}


## <<strip_plots>>
## these plots are built by extracting the mean/median
## of each biological replicate of each strain.  Thus,
## we reduce 10,000 observations of a replicate to a
## single value.  w/ 8 biological replicates per strain,
## we can make a nice stripplot of strain * reporter
## for the different parameters.  We'll also use this
## to make a levelplot/heatmap

## 'aggregate' creates a new dataframe from x by applying FUN to
## all unique combinations of the factors supplied to the 'by'
## argument - in this case, grab the mean of numeric data and
## keep everything else a factor 
out_agg <- aggregate.data.frame(x = out_all,
                                by = list(out_all$strain,
                                          out_all$reporter,
                                          out_all$replicate),
                            FUN = function(x) {
                                ifelse(is.numeric(x), mean(x), as.factor(x))
                            })

## 'aggregate' seems to strip the levels from factors, so add
## these back using the values present in the original dataframe
out_agg$strain <- factor(out_agg$strain,
                     levels = unique(out_agg$strain),
                     labels = levels(out_all$strain))

out_agg$reporter <- factor(out_agg$reporter,
                     levels = unique(out_agg$reporter),
                     labels = levels(out_all$reporter))

out_agg$strain_rep <- factor(out_agg$strain_rep,
                             levels = unique(out_agg$strain_rep),
                             labels = levels(out_all$strain_rep))

## need a color vector that maps to strains that's the
## length of levels(strain) * levels(replicate) (usually 40)
out_agg_cols <- vector()
for (i in 1:nrow(out_agg)) {
    out_agg_cols[i] <- all_cols[out_agg$strain[i] == names(all_cols)]
}

## convert to Z scores 
i <- 1

out_z <- list()
reporter <- levels(out_agg$reporter)
for (i in 1:length(reporter)) {
    
    z_mean <- mean(out_agg$TFT_ratio[out_agg$reporter == reporter[i]])
    z_sd   <- sd(out_agg$TFT_ratio[out_agg$reporter == reporter[i]])
    z_out  <- function(x) { ((x - z_mean) / (z_sd))  }
    out_subset <- out_agg[out_agg$reporter == reporter[i], ]
    
    for (j in 1:nrow(out_subset)) {
        out_subset$zTFT[j] <- z_out(out_subset$TFT_ratio[j])
        out_z[[i]] <- out_subset
}}

out_z <- do.call("rbind", out_z)
range(out_z$zTFT)

new_aa_order <- c(2, 3, 4, 6, 7, 9, 10, 11, 12, 14, 18, 19, 1, 5, 8, 13, 15, 16, 17, 20)
out_z$r_order <- factor(out_z$reporter,
                        levels = levels(out_z$reporter)[rev(new_aa_order)])

out_z$t_order <- factor(out_z$reporter,
                        levels = levels(out_z$reporter)[rev(aa_order)])

## aggregate into a smaller data frame for the actual heatmap
out_m <- aggregate.data.frame(x = out_z,
                               by = list(out_z$strain,
                                         out_z$reporter),
                               FUN = function(x) {
                                   ifelse(is.numeric(x), mean(x), as.factor(x))
                               })

out_m$strain <- factor(out_m$strain,
                        levels = unique(out_m$strain),
                        labels = levels(out_all$strain))

out_m$reporter <- factor(out_m$reporter,
                          levels = unique(out_m$reporter),
                          labels = levels(out_all$reporter))

out_m$r_order <- factor(out_m$reporter,
                        levels = levels(out_m$reporter)[rev(new_aa_order)])

out_m$t_order <- factor(out_m$reporter,
                        levels = levels(out_m$reporter)[rev(aa_order)])


hm_cols <- rev(brewer.pal(n = 11, name = "PiYG"))

{
pdf(file = "~/Desktop/2022.01.20_TFT_heat_final.pdf", height = 14, width = 6)
print(
levelplot(zTFT ~ strain * t_order,
          strip = T,
          xlab = "normalized TFT score",
          ylab = "Reporter",
          data = out_m,
          pretty = T,
          col.regions = hm_cols,
          ylab.right = "Z score",
          scales = list(alternating = F,
                        x = list(labels = c("BY", "RM",
                                            expression(paste("BY rpn4", Delta)), 
                                            expression(paste("BY ubr1", Delta)), 
                                            expression(paste("BY doa10", Delta))),
                                 rot = 45),
                        ## scales only on left/bottom
                        tck = c(1, 0)),
          par.settings = list(strip.background = list(col = gray(0.9)),
                              clip = list(panel = FALSE),
                              par.sub.text = list(font = 2,
                                                  just = "center",
                                                  cex = 1.25,
                                                  x = grid::unit(7.05, "in"),
                                                  y = grid::unit(6.25, "in")),
                              layout.heights = list(xlab.key.padding = 1),
                              layout.widths = list(left.padding = 5),
                              axis.text = list(cex = 1.25),
                              par.ylab.text = list(cex = 1.25,
                                                   col = "white"),
                              par.xlab.text = list(cex = 1.25)), 
##          col.regions = gray(level = 29:0/29),
          at = seq(from = -2, to = 2, length.out = 10),
          colorkey = list(at = seq(from = -2, to = 2, length.out = 11),
                          cex = 1.25,
                          title = "normalized TFT score",
                          space = "bottom",
                          columns = 1,
                          row = 10),
          border = gray(0.3),
          as.table = T,
          ## index.cond = list(aa_order),
          border.lwd = 2)
)
## Ac/N-end reporter lines
grid.lines(x = c(0.65, 0.65),
           y = c(2.4, 6.42),
           default.units = "in",
           gp = gpar(lwd = 2, col = gray(0.4)))
## Ac/N-end reporter text
grid.text(x = 0.3, y = 4.41,
          label = "Ac/N-end\nReporters",
          default.units = "in",
          rot = 90,
          gp = gpar(cex = 1.25))
## type I Arg/N-end reporter lines
grid.lines(x = c(0.65, 0.65),
           y = c(9.83, 13.25),
           default.units = "in",
           gp = gpar(lwd = 2, col = gray(0.4)))
## type I Arg/N-end reporter text
grid.text(x = 0.3, y = 11.5,
          label = "Type I Arg/N-end\nReporters",
          default.units = "in",
          rot = 90,
          gp = gpar(cex = 1.25))
## type II Arg/N-end reporter lines
grid.lines(x = c(0.65, 0.65),
           y = c(6.95, 9.25),
           default.units = "in",
           gp = gpar(lwd = 2, col = gray(0.4)))
## type II Arg/N-end reporter text
grid.text(x = 0.3, y = 8.1,
          label = "Type II Arg/N-end\nReporters",
          default.units = "in",
          rot = 90,
          gp = gpar(cex = 1.25))
dev.off()
}


## <<Z_score_stripplot_final>>
{
pdf(file = "~/Desktop/2022.01.20_z_strip_full_final.pdf", height = 14.5, width = 14)
print(
stripplot(zTFT ~ strain | reporter,
          data = out_z,
          col = gray(0),
          fill = out_agg_cols,
          pch = 21,
          cex = 1.1,
          scales = list(alternating = 3,
                        x = list(labels = c("BY", "RM",
                         expression(paste("BY rpn4", Delta)), 
                         expression(paste("BY ubr1", Delta)), 
                         expression(paste("BY doa10", Delta))),
                         rot = 45)),
          layout = c(4, 5),
          grid = T,
          main = list(label = "Arg/N-end Reporters"),
          sub = list(label = "Ac/N-end Reporters"),
          between = list(x = c(0, 0, 0),
                         y = c(0, 0, 3)),
          as.table = T,
          par.settings = list(strip.background = list(col = gray(0.9)),
                              clip = list(panel = FALSE),
                              par.main.text = list(font = 2,
                                                   cex = 1.25,
                                                   just = "center", 
                                                   x = grid::unit(7, "in"),
                                                   y = grid::unit(13, "in")),
                              par.sub.text = list(font = 2,
                                                  just = "center",
                                                  cex = 1.25,
                                                  x = grid::unit(7.05, "in"),
                                                  y = grid::unit(6.25, "in")),
                              axis.text = list(cex = 1),
                              par.ylab.text = list(cex = 1.25,
                                                   col = "white"),
                              par.xlab.text = list(cex = 1.25)),
          jitter.data = T,
          ## factor for jittering
          factor = 1.5,
          index.cond = list(aa_order),
          horizontal = F)
)
## main title, not sure why this doesn't show up otherwise
## probably some par setting re: going outside the grid....
grid.text(label = "Arg/N-end Reporters",
          x = 7.0,
          y = 14.2,
          default.units = "in",
          gp = gpar(cex = 1.25, font = 2))
## y axis labels for the strip
grid.text(label = "normalized TFT score",
          x = c(0.15, 0.15),
          y = c(10.17, 3.9),
          rot = 90,
          default.units = "in",
          gp = gpar(col = "black", cex = 1.25)
          )
dev.off()
}


## <<horizontal_N-end_heatmap>>
## 1. reverse levels on strain
out_m$strain_rev <- factor(out_m$strain,
                           levels = levels(out_m$strain)[5:1])
## 2. reverse levels on reporter
out_m$t_order_rev <- factor(out_m$t_order,
                           levels = levels(out_m$t_order)[20:1])

{
pdf(file = "~/Desktop/2022.01.20_N-end_horizontal_heatmap.pdf", height = 5, width = 14)
print(
levelplot(zTFT ~ t_order_rev * strain_rev,
          data = out_m,
          strip = T,
          ylab = "Strain",
          pretty = T,
          scales = list(alternating = F,
                        y = list(labels = c(expression(paste("BY doa10", Delta)),
                                            expression(paste("BY ubr1", Delta)),
                                            expression(paste("BY rpn4", Delta)),
                                            "RM", "BY")),
                        x = list(labels = gsub(pattern = " TFT", replacement = "",
                                               x = levels(out_m$t_order_rev)),
                                 rot = 45),
                        ## scales only on left/bottom
                        tck = c(1, 0)),
          par.settings = list(strip.background = list(col = gray(0.9)),
                              clip = list(panel = FALSE),
                              par.sub.text = list(font = 2,
                                                  just = "center",
                                                  cex = 1.5,
                                                  x = grid::unit(7.05, "in"),
                                                  y = grid::unit(6.25, "in")),
                              layout.heights = list(xlab.key.padding = 1,
                                                    bottom.padding = 3),
                              layout.widths = list(left.padding = 0,
                                                   right.padding = 10),
                              axis.text = list(cex = 1.5),
                              par.ylab.text = list(cex = 1.5,
                                                   col = "white"),
                              par.xlab.text = list(cex = 1.5,
                                                   col = "white")), 
          ## at = 1 less than number of ats in colorkey
          col.regions = brewer.pal(n = 11, name = "PiYG"),
          at = seq(from = 0, to = 3.6, length.out = 10),
          colorkey = list(at = 0:9,
                          labels = as.character(format(seq(from = 0, to = 3.6, length.out = 10),
                                                       nsmall = 1)),
                          cex = 1.5,
                          title = "Normalized UPS\nPathway Activity",
                          space = "right",
                          columns = 1,
                          row = 9),
          border = gray(0.3),
          as.table = T,
          border.lwd = 2)
)
## 3. add rotated text to colorkey
grid.text(label = "Normalized UPS\nPathway Activity",
          x = 13.6,
          y = 2.98,
          default.units = "in",
          rot = 90,
          gp = gpar(cex = 1.5, font = 1))
## 4. add pathway lines 
## Type I Arg/N-end text
grid.text(label = "Type I Arg/N-end\nReporters",
          x = 3.6125,
          y = 0.35,
          default.units = "in",
          gp = gpar(cex = 1.5, font = 1))
## Type II Arg/N-end text
grid.text(label = "Type II Arg/N-end\nReporters",
          x = 6.725,
          y = 0.35,
          default.units = "in",
          gp = gpar(cex = 1.5, font = 1))
## Ac/N-end text
grid.text(label = "Ac/N-end\nReporters",
          x = 10.1,
          y = 0.35,
          default.units = "in",
          gp = gpar(cex = 1.5, font = 1))
## Type I Arg/N-end lines
grid.lines(x = c(2.05, 5.175),
          y = 0.75,
          default.units = "in",
          gp = gpar(col = gray(0), lwd = 2))
## Type II Arg/N-end lines
grid.lines(x = c(5.7, 7.75),
          y = 0.75,
          default.units = "in",
          gp = gpar(col = gray(0), lwd = 2))
## Ac/N-end lines
grid.lines(x = c(8.28, 11.92),
          y = 0.75,
          default.units = "in",
          gp = gpar(col = gray(0), lwd = 2))
dev.off()
}

## <<example_stripplots_final>>
## use met and trp for this purpose
## these will go into figure 1, but
## are just for illustration purposes...
out_x <- out_z[out_z$reporter == "Met TFT" | out_z$reporter == "Trp TFT", ]

{
pdf(file = "~/Desktop/2022.01.20_example_strip_2_panel_final.pdf", height = 8, width = 5)
print(
stripplot(zTFT ~ strain | reporter,
          data = out_x,
          col = gray(0),
          fill = out_agg_cols,
          pch = 21,
          cex = 1.25,
          scales = list(alternating = 1,
                        tck = c(1, 0),
                        x = list(labels = c("BY", "RM",
                         expression(paste("BY rpn4", Delta)), 
                         expression(paste("BY ubr1", Delta)), 
                         expression(paste("BY doa10", Delta))),
                         rot = 45)),
          layout = c(1, 2),
          grid = T,
          as.table = F,
          jitter = T,
          factor = 1.7,
          ylab = "Z-score",
          xlab = "Strain",
          par.settings = list(strip.background = list(col = gray(0.9)),
                              clip = list(panel = FALSE),
                              par.main.text = list(font = 2,
                                                   cex = 1.25,
                                                   just = "center", 
                                                   x = grid::unit(7, "in"),
                                                   y = grid::unit(13, "in")),
                              par.sub.text = list(font = 2,
                                                  just = "center",
                                                  cex = 1.25,
                                                  x = grid::unit(7.05, "in"),
                                                  y = grid::unit(6.25, "in")),
                              axis.text = list(cex = 1.25),
                              par.ylab.text = list(cex = 1.25,
                                                   col = "white"),
                              par.xlab.text = list(cex = 1.25)),
          horizontal = F)
)
grid.text(label = "normalized TFT score",
          x = c(0.2, 0.2),
          y = c(2.75, 5.85),
          rot = 90,
          gp = gpar(cex = 1.25),
          default.units = "in")
dev.off()
}


## <<BY_RM_plot>>
## nrow(out_agg) = 800, so 300 should be BY or RM 
out_by_rm <- grepl(pattern = "[BR][YM]",
                   x = as.character(out_agg$strain))

out_br <- out_agg[out_by_rm, ]

out_br_cols <- vector()

for (i in 1:nrow(out_br)) {
out_br_cols[i] <- all_cols[out_br$strain[i] == names(all_cols)]
}

{
pdf(file = "~/Desktop/2022.01.20_by_rm_final.pdf", height = 14.5, width = 14)
print(
## need to make a TFT ratio amino acid index 
stripplot(TFT_ratio ~ strain | reporter,
          data = out_br,
          type = c("g", "p"),
          col = gray(0),
          fill = out_br_cols,
          pch = 21,
          cex = 1.1,
          scales = list(alternating = 3,
                        x = list(labels = c("BY", "RM"))),
          layout = c(4, 5),
          grid = T,
          main = list(label = "Arg/N-end Reporters"),
          sub = list(label = "Ac/N-end Reporters"),
          between = list(x = c(0, 0, 0),
                         y = c(0, 0, 3)),
          as.table = T,
          par.settings = list(strip.background = list(col = gray(0.9)),
                              clip = list(panel = FALSE),
                              par.main.text = list(font = 2,
                                                   cex = 1.25,
                                                   just = "center", 
                                                   x = grid::unit(7, "in"),
                                                   y = grid::unit(13, "in")),
                              par.sub.text = list(font = 2,
                                                  just = "center",
                                                  cex = 1.25,
                                                  x = grid::unit(7.05, "in"),
                                                  y = grid::unit(5.95, "in")),
                              axis.text = list(cex = 1),
                              par.ylab.text = list(cex = 1.25,
                                                   col = "white"),
                              par.xlab.text = list(cex = 1.25)),
          jitter.data = T,
          ## factor for jittering
          factor = 1.5,
          index.cond = list(aa_order),
          horizontal = F)
)
## main title, not sure why this doesn't show up otherwise
## probably some par setting re: going outside the grid....
grid.text(label = "Arg/N-end Reporters",
          x = 7.0,
          y = 14.2,
          default.units = "in",
          gp = gpar(cex = 1.25, font = 2))

## y axis labels for the strip
grid.text(label = "log2 TFT ratio",
          x = c(0.15, 0.15),
          y = c(7.45, 1.85),
          rot = 90,
          default.units = "in",
          gp = gpar(col = "black", cex = 1.25))
dev.off() 
}
