## -----
## load all required packages
source("~/emacs/R/functions/load_flow_packages.R")

needed_dirs <- c("/fcs", "/results", "/tables",
                 "/scripts", "/dataframes", "/sessions",
                 "/dataframes/gated", "/dataframes/ungated")

dir_maker <- function(x){
    ifelse(!dir.exists(paths = paste0(base_dir, x)),
           dir.create(path = paste0(base_dir, x)),
           paste0("dir ", paste0(base_dir, x), " exists."))
}

sapply(X = needed_dirs, FUN = dir_maker)
work_dir       <- paste0(base_dir, "/fcs")
results_dir    <- paste0(base_dir, "/results")
tables_dir     <- paste0(base_dir, "/tables")
sessions_dir   <- paste0(base_dir, "/sessions")
frame_dir      <- paste0(base_dir, "/dataframes")
gated_dir      <- paste0(frame_dir, "/gated/")
ungated_dir    <- paste0(frame_dir, "/ungated/")
out_log        <- paste0(results_dir, "/output_log")


## -----
## read in all the fcs files in a directory 
all_flow <- read.flowSet(files = NULL,
                         path = work_dir,
                         pattern = ".*.fcs",
                         alter.names = T,
                         min.limit = 2)

 
## get the number of cells per file; write output
cells_per_file <- fsApply(x = all_flow,
                          FUN = function(x) {
                              nrow(as.data.frame(exprs(x))) 
                          })

## write information string to the output log
read_string <- paste0("1. On ",
                      Sys.time(),
                      " files were read in via 'read.flowSet'.\n\nCells per file:\n")

write(x = read_string,
      file = out_log,
      append = F)

## for nicely aligned output in the resulting file
cat(capture.output(cells_per_file),
    file = out_log,
    append = T,
    sep = "\n")

all_filtered <- all_flow

## -----
## now convert flowframes to dataframes and merge to a single frame
## start w/ ungated cells
ungated_frames <- fsApply(x = all_filtered,
                          FUN = function(x) {
                              
                              ## extract file name; we'll
                              ## add this to the dataframe
                              file_name <- rep(x = x@description$GUID,
                                               times = nrow(as.data.frame(exprs(x))))
                              
                              ## extract the time the sample was run;
                              ## then convert to a numeric value
                              ## BTIM = "beginning time"
                              file_time <- rep(x = x@description$`$BTIM`,
                                               times = nrow(as.data.frame(exprs(x))))
                              time_conv <- HmsToSec(file_time)
                              
                              ## combine flow data and file name
                              ## into a single data frame 
                              cbind(as.data.frame(exprs(x)),
                                    file_name,
                                    file_time,
                                    time_conv)
})

ungated_final <- do.call("rbind", ungated_frames)


## -----
## log fluorescence values and compute TFT ratio
ungated_final$log_GFP   <- log10(x = ungated_final$eGFP.A)
ungated_final$log_RFP   <- log10(x = ungated_final$mCherry.A)
ungated_final$TFT_ratio <- -log(x = ungated_final$log_RFP / ungated_final$log_GFP,
                                base = 2)

## unique(is.na(ungated_final$TFT_ratio))


## -----
## now gate the cells to capture the haploid cell population.
## we identify haploids as a sharp peak in the lower end of
## the fsc density plot.  I take 10% above and below the
## max density value.
ungated_final$file_string <- as.character(ungated_final$file_name)

## ungated = all cells
ungated_final$ungated <- rep(T, nrow(ungated_final))

## gated = all cells w/in 10% +/- FSC max density
densities <- lapply(X = unique(ungated_final$file_string),
                    FUN = function(x) {
                        
                        ## file name for dataframe creation
                        f_name <- x

                        ## per file dataframe for calculating density
                        dat  <- ungated_final[ungated_final$file_name == f_name, ]
                        
                        ## get the density for FSC
                        dens <- density(dat$FSC.A)
                        
                        ## 10% +/- the max FSC density peak
                        dens_max <- dens$x[which.max(dens$y)]
                        dens_up <- (0.1 * dens_max) + dens_max
                        dens_down <- dens_max - (0.1 * dens_max)
                        data.frame(f_name, dens_up, dens_down, stringsAsFactors = F)
                    })

## bind each file's density estimate into a single dataframe
densities <- do.call("rbind", densities)

ungated_final$gated <- sapply(X = 1:nrow(ungated_final), 
                              FUN = function(x) {
                                  
                                  ## get file name for current row
                                  f_name <- ungated_final[x, "file_string"]

                                  ## get corresponding row of density estimates
                                  dens_row <- densities[densities$f_name == f_name, ]
                                  
                                  ## 10% +/- the max FSC density peak
                                  dens_up   <- dens_row$dens_up
                                  dens_down <- dens_row$dens_down
                                  
                                  ## create a gate and subset to keep cells in range
                                  ifelse(ungated_final[x, "FSC.A"] > dens_down &
                                         ungated_final[x, "FSC.A"] < dens_up,
                                         yes = T, no = F)

                     })


ungated_per_file <- tapply(X = ungated_final$FSC.A,
                           INDEX = ungated_final$file_name,
                           FUN = length)

gated_per_file <- tapply(X = ungated_final$FSC.A[ungated_final$gated == T],
                         INDEX = ungated_final$file_name[ungated_final$gated == T],
                         FUN = length)

fsc_gated_counts <- data.frame(kept = gated_per_file,
                               excluded = ungated_per_file - gated_per_file,
                               percent = gated_per_file / ungated_per_file)

gate_string <- paste0("\n3. On ", Sys.time(),
                      " Filtered cells to grab 10% +/- the central FSC peak\n\n",
                      "The following counts were obtained:\n")

write(x = gate_string,
      file = out_log,
      append = T)

cat(capture.output(fsc_gated_counts),
    file = out_log,
    append = T,
    sep = "\n")


## -----
## write out summary statistics for gated and ungated files
## extract only numeric parameters for summary statistics 
f_params <- colnames(ungated_final)
f_params <- f_params[sapply(X = f_params,
           FUN = function(x) {
               is.numeric(ungated_final[, x])
               })]

## -----
## write out summary statistics for all parameters 
## make sure all columns are printed together 
options(width = 200)
for(p in 1:length(f_params)) {

    gated_final    <- ungated_final[ungated_final$gated == T, ]

    ungated_mean   <- tapply(X = ungated_final[, f_params[p]],
                             INDEX = ungated_final$file_name,
                             FUN = mean)
    gated_mean     <- tapply(X = gated_final[, f_params[p]],
                             INDEX = gated_final$file_name,
                             FUN = mean)
    mean_frac      <- round(x = gated_mean / ungated_mean,
                            digits = 2)
    ungated_median <- tapply(X = ungated_final[, f_params[p]],
                             INDEX = ungated_final$file_name,
                             FUN = median)
    gated_median   <- tapply(X = gated_final[, f_params[p]],
                             INDEX = gated_final$file_name,
                             FUN = median)
    median_frac      <- round(x = gated_median / ungated_median,
                              digits = 2)
    ungated_sd     <- tapply(X = ungated_final[, f_params[p]],
                             INDEX = ungated_final$file_name,
                             FUN = sd)
    gated_sd       <- tapply(X = gated_final[, f_params[p]],
                             INDEX = gated_final$file_name,
                             FUN = sd)
        
    ungated_cv     <- ungated_mean / ungated_sd

    gated_cv       <- gated_mean / gated_sd

    sum_stats      <- data.frame(ungated_mean, gated_mean, mean_frac,
                                 ungated_median, gated_median, median_frac,
                                 ungated_sd, gated_sd,
                                 ungated_cv, gated_cv)

    sum_string     <- paste0("\nOn ", Sys.time(), " Obtained summary statistics for parameter: ",
                             f_params[p], "\n\n")

    write(x = sum_string,
      file = out_log,
      append = T)

    cat(capture.output(sum_stats),
        file = out_log,
        append = T,
        sep = "\n")
    
}


## -----
## extract strain, reporter, and replicate as factors

## now, extract the strain from the file name string
## this relies on naming the files w/ strain, reporter,
## and replicated separated by a '-' character, as in:
## "BY_full_gRNA_02-ODC_TFT-003.fcs" for the:
## strain - "BY full gRNA 02"
## reporter - ODC TFT
## replicate - 3

## because of 'strsplits' list output format, have to use
## 'sapply' to create the 'strain' variable.  Output is, e.g.,:
## [[1]]
## [1] "BY_full_gRNA_02" "ODC_TFT"         "003.fcs"
## so, take the 1st index of the appropriate list item
ungated_final$strain <- sapply(X = 1:length(ungated_final$file_string),
                               FUN = function(x){
                                   strsplit(x = ungated_final$file_string[x], split = "-")[[1]][1]
                               })
ungated_final$strain_factor <- as.factor(ungated_final$strain)

## now write the unique values we obtain for strain in ea. dataset
write(x = sprintf("%s", c("\nUnique values for 'strain' variable:\n",
                          unique(ungated_final$strain))),
      file = out_log,
      append = T)

## now for reporter
ungated_final$reporter <- sapply(X = 1:length(ungated_final$file_string),
                                 FUN = function(x){
                                     strsplit(x = ungated_final$file_string[x], split = "-")[[1]][2]
                                 })
ungated_final$reporter_factor <- as.factor(ungated_final$reporter)

## log the results
write(x = sprintf("%s", c("\nUnique values for 'reporter' variable for ungated set:\n",
                          unique(ungated_final$reporter))),
      file = out_log,
      append = T)


## now for replicates
ungated_final$replicate <- sapply(X = 1:length(ungated_final$file_string),
                                  FUN = function(x){
                                      as.numeric(
                                          gsub(pattern = ".fcs",
                                               replacement = "",
                                               x = strsplit(x = ungated_final$file_string[x],
                                                            split = "-")[[1]][3]))
                                  })
## for later use in plotting individual replicates
ungated_final$replicate_factor <- as.factor(ungated_final$replicate)


## log the results
write(x = sprintf("%s", c("\nUnique values for 'replicate' variable for ungated set:\n",
                          unique(ungated_final$replicate))),
      file = out_log,
      append = T)


## -----
## get plate information and create a variable that 
## expresses time relative to the first sample for ea. plate
ungated_final$plate <- ifelse(test = ungated_final$replicate <= 8,
                              yes  = "plate_01",
                              no   = "plate_02")

## ungated_final$plate <- as.numeric(sapply(X = sapply(X = ungated_final$file_string,
##                                                     FUN = function(x) {
##                                                         x = unlist(strsplit(x = x,
##                                                                             split = ".*-"))[2]}),
##                                          FUN = function(y) {
##                                              unlist(strsplit(x = y,
##                                                              split = ".fcs")[1])
##                                          }))

ungated_final$plate <- ifelse(test = ungated_final$replicate <= 8,
                              yes  = "plate_01",
                              no   = "plate_02")

plate_one_min <- min(ungated_final$time_conv[ungated_final$plate == "plate_01"])
plate_two_min <- min(ungated_final$time_conv[ungated_final$plate == "plate_02"])

## convert time to a relative value based on the first sample of ea. plate
ungated_final$time_rel <- sapply(X = 1:nrow(ungated_final),
                                 FUN = function(x) {
                                     ifelse(test = ungated_final[x, "plate"] == "plate_01",
                                            yes  = ungated_final[x, "time_conv"] - plate_one_min, 
                                            no   = ungated_final[x, "time_conv"] - plate_two_min)
                                 })

## log relative time values per file
time_string <- paste0("\nOn ", Sys.time(), "Extracted time values for each file/plate\n")

write(x = time_string,
      file = out_log,
      append = T) 

time_frame <- lapply(X = unique(ungated_final$file_string),
                     FUN = function(x) {

                         sub_set <- ungated_final$file_string == x
                         f_params <- c("file_string", "strain", "plate",
                                     "file_time", "time_conv", "time_rel")
                         dat <- ungated_final[sub_set, f_params]
                         dat_out <- data.frame(file = unique(dat$file_string),
                                               strain = unique(dat$strain),
                                               plate = unique(dat$plate),
                                               file_time = unique(dat$file_time),
                                               time_conv = unique(dat$time_conv),
                                               time_rel = unique(dat$time_rel))
                     })

time_frame <- do.call("rbind", time_frame)

cat(capture.output(time_frame),
    file = out_log,
    append = T,
    sep = "\n")


## -----
## save the final dataframe for future use:
save(ungated_final,
     file = paste0(frame_dir, "/ungated_final.R"))

write(x = paste0("\nOn ", Sys.time(), " saved final dataframe as: ",
                 paste0(frame_dir, "/", "ungated_final.R"),
                 "\n"),
                 file = out_log,
                 append = T)

save(work_dir, results_dir, tables_dir, sessions_dir,
     frame_dir, gated_dir, ungated_dir, out_log,
     file = paste0(frame_dir, "/dir_structure.R"))

write(x = paste0("\nOn ", Sys.time(), " saved dir structure as: ",
                 paste0(frame_dir, "/dir_structure.R"),
                 "\n"),
                 file = out_log,
                 append = T)

## -----
## re-load the data by uncommenting the following:
## load(file = paste0(frame_dir, "/", "ungated_final.R"))
## load(file = paste0(frame_dir, "/dir_structure.R"))
## source("~/emacs/R/functions/load_flow_packages.R")


## -----
## extract the median from each biological replicate and use this
## value to build a dataframe w/ n. replicate observations per strain
## per reporter.  This dataframe is what we'll use for stats and for
## creating stripcharts, boxplots, and heatmaps

## 'aggregate' creates a new dataframe from x by applying FUN to
## all unique combinations of the factors supplied to the 'by'
## argument - in this case, grab the mean of numeric data and
## keep everything else a factor
ungated_medians <- aggregate.data.frame(x = ungated_final,
                                        by = list(ungated_final$strain_factor,
                                                  ungated_final$replicate_factor,
                                                  ungated_final$reporter_factor),
                                        FUN = function(x) {
                                            ifelse(is.numeric(x), median(x), as.character(x))
                                        },
                                        ## simplify results to vector 
                                        simplify = T)


## -----
## 'aggregate' seems to strip the levels from factors, so add
## these back using the values present in the original dataframe
ungated_medians$strain_factor    <- as.factor(ungated_medians$strain)
ungated_medians$replicate_factor <- as.factor(ungated_medians$replicate)
ungated_medians$reporter_factor  <- as.factor(ungated_medians$reporter)
ungated_medians$gating           <- rep("ungated", nrow(ungated_medians))


## -----
## subtract the deletion strain value from each replicate
## the deletion strain TFT ratio will be ~0, but not exactly
## thus, we normalize to it before converting to z-scores 
## ungated_medians[which.min(ungated_medians$TFT_ratio), ]
TFT_min <- min(ungated_medians$TFT_ratio)
ungated_medians$TFT_sub <- ungated_medians$TFT_ratio - TFT_min
ungated_medians <- ungated_medians[ungated_medians$strain != "doa10_strain", ]
ungated_medians <- ungated_medians[ungated_medians$strain != "ubr1_strain", ]
ungated_medians <- ungated_medians[ungated_medians$strain != "nta1_strain", ]
ungated_medians <- ungated_medians[ungated_medians$strain != "rpn4_strain", ]
ungated_medians <- ungated_medians[ungated_medians$strain != "ubc6_strain", ]


## -----
## adjust for the effect of time on the TFT ratio,
## then convert it to a Z-score: 
ungated_loess <- loess(formula = TFT_sub ~ time_rel,
                       data = ungated_medians)
ungated_medians$TFT_loess <- ungated_loess$residuals + mean(ungated_medians$TFT_ratio)

## scale to wild-type median with sd = 1
by_median <- median(ungated_medians$TFT_loess[ungated_medians$strain == "BY_full"])
ungated_medians$TFT_scaled <- scale(x = ungated_medians$TFT_loess,
                                    center = by_median,
                                    scale = T)

## -----
## re-level prior to statistical analysis
ungated_medians$strain_ordered <- factor(ungated_medians$strain_factor,
                                         levels = levels(ungated_medians$strain_factor)[c(2, 1,
                                                                                          6, 5,
                                                                                          4, 7,
                                                                                          3)])

## drop the deletion strain level from the dataset
ungated_medians$strain_final <- droplevels(x = ungated_medians$strain_ordered)
## levels(ungated_medians$strain_final)

## -----
## have to extract medians from the gated data separately
gated_final <- ungated_final[ungated_final$gated == T, ]

gated_medians <- aggregate.data.frame(x = gated_final,
                                       by = list(gated_final$strain_factor,
                                                 gated_final$replicate_factor,
                                                 gated_final$reporter_factor),
                                       FUN = function(x) {
                                           ifelse(is.numeric(x), median(x), as.character(x))
                                       },
                                       ## simplify results to vector 
                                       simplify = T)

## 'aggregate' seems to strip the levels from factors, so add
## these back using the values present in the original dataframe
gated_medians$strain_factor    <- as.factor(gated_medians$strain)
gated_medians$replicate_factor <- as.factor(gated_medians$replicate)
gated_medians$reporter_factor  <- as.factor(gated_medians$reporter)
gated_medians$gating           <- rep("gated", nrow(gated_medians))

## -----
## subtract the deletion strain value from each replicate
## the deletion strain TFT ratio will be ~0, but not exactly
## thus, we normalize to it before converting to z-scores 
## gated_medians[which.min(gated_medians$TFT_ratio), ]
TFT_min <- min(gated_medians$TFT_ratio)
gated_medians$TFT_sub <- gated_medians$TFT_ratio - TFT_min
gated_medians <- gated_medians[gated_medians$strain != "doa10_strain", ]
gated_medians <- gated_medians[gated_medians$strain != "ubr1_strain", ]
gated_medians <- gated_medians[gated_medians$strain != "nta1_strain", ]
gated_medians <- gated_medians[gated_medians$strain != "rpn4_strain", ]
gated_medians <- gated_medians[gated_medians$strain != "ubc6_strain", ]


## -----
## adjust for the effect of time on the TFT ratio,
## then convert it to a Z-score: 
gated_loess <- loess(formula = TFT_sub ~ time_rel,
                     data = gated_medians)
gated_medians$TFT_loess <- gated_loess$residuals + mean(gated_medians$TFT_ratio)

## scale to wild-type median with sd = 1
by_median <- median(gated_medians$TFT_loess[gated_medians$strain == "BY_full"])
gated_medians$TFT_scaled <- scale(x = gated_medians$TFT_loess,
                                  center = by_median,
                                  scale = T)

## -----
## re-level prior to statistical analysis
gated_medians$strain_ordered <- factor(gated_medians$strain_factor,
                                       levels = levels(gated_medians$strain_factor)[c(2, 1,
                                                                                      6, 5,
                                                                                      4, 7,
                                                                                      3)])

gated_medians$strain_final <- droplevels(x = gated_medians$strain_ordered)
## levels(gated_medians$strain_final)


## -----
## re-load point 
## dump the gated and ungated cells into a single list
all_data <- list()
all_data[[1]] <- ungated_medians
all_data[[2]] <- gated_medians
names(all_data) <- c("ungated_medians", "gated_medians")
save(all_data,
     file = paste0(frame_dir, "/all_data.R"))


## -----
## ## re-load point 
## load(file = paste0(frame_dir, "/dir_structure.R"))
## load(file = paste0(frame_dir, "/all_data.R"))


## -----
## 'lapply' to get the statistics
f_params <- c("TFT_ratio", "TFT_loess", "TFT_scaled")

lapply(X = all_data,
       FUN = function(n) {
           
           sapply(X = f_params,
                  FUN = function(p) {

                      ## -----
                      ## write ungated aov output to file
                      write(x = paste0("-------------------------------------\n",
                                       "|", unique(n$gating),
                                       " cells: ANOVA for parameter:|\n",
                                       "-------------------------------------\n",
                                       p, " on", " ", Sys.time(), "\n"),
                            file = stats_log,
                            append = T)
                      
                      test_aov <- aov(n[, p] ~ ungated_medians[, "strain_final"])
                      
                      capture.output(summary(test_aov),
                                     file = stats_log,
                                     append = T,
                                     type = "output",
                                     split = F)
                      
                      write(x = "\n",
                            file = stats_log,
                            append = T)

                      ## -----
                      ## now the Tukey HSD posthoc test for the ungated set
                      test_posthoc <- TukeyHSD(test_aov)
                      test_posthoc_out <- as.data.frame(test_posthoc[[1]])
                      test_posthoc_out$reporter <- rep(unique(n$reporter),
                                                       nrow(test_posthoc_out))
                      
                      write(x = paste0("---------------------------------------\n",
                                       "|", unique(n$gating),
                                       " cells: Tukey HSD for parameter:|\n",
                                       "---------------------------------------\n",
                                       p, " on", " ", Sys.time(), "\n"),
                            file = stats_log,
                            append = T)
                      
                      capture.output(x = test_posthoc_out,
                                     file = stats_log,
                                     append = T,
                                     type = "output",
                                     split = F)
                      
                      write(x = "\n",
                            file = stats_log,
                            append = T)

                      ## -----
                      ## pairwise t-test posthoc (for comparison to Tukey) for the ungated set
                      write(x = paste0("----------------------------------------\n",
                                       "|", unique(n$gating),
                                       " cells: Pairwise t for parameter:|\n",
                                       "----------------------------------------\n",
                                       p, " on", " ", Sys.time(), "\n"),
                            file = stats_log,
                            append = T)
                      
                      capture.output(
                          pairwise.t.test(x = n[, p],
                                          g = n[, "strain_final"],
                                          p.adjust.method = "BH",
                                          alternative = "two.sided")$p.value,
                          file = stats_log,
                          append = T,
                          type = "output",
                          split = F)

                                            ## -----
                      ## pairwise t-test posthoc (for comparison to Tukey) for the ungated set
                      write(x = paste0("-------------------------------------------------\n",
                                       "|", unique(n$gating),
                                       " cells: Pairwise t w/ no pooled SD for parameter:|\n",
                                       "-------------------------------------------------\n",
                                       p, " on", " ", Sys.time(), "\n"),
                            file = stats_log,
                            append = T)
                      
                      capture.output(
                          pairwise.t.test(x = n[, p],
                                          g = n[, "strain_final"],
                                          p.adjust.method = "BH",
                                          pool.sd = F,
                                          alternative = "two.sided")$p.value,
                          file = stats_log,
                          append = T,
                          type = "output",
                          split = F)
                      
                      write(x = "\n",
                            file = stats_log,
                            append = T)

                  }
                , simplify = F)
})
