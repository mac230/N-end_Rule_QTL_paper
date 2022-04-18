## -----
## ## re-load point 
load(file = paste0(frame_dir, "/dir_structure.R"))
load(file = paste0(frame_dir, "/all_data.R"))

ungated_medians <- all_data[[1]]
gated_medians   <- all_data[[2]]

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
