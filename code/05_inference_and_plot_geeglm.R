source("code/03_plot_exploratory_analysis.R")
source("code/04_geeglm_full_run.R")

##### Results of ungrouped analysis #####
lapply(gee_mods_ungrouped, summary) # GEE summary
compute_all(gee_mods_ungrouped$coefs_vcov, c(0, 1, 0, 0)) # Comparision of Beta1 before
compute_all(gee_mods_ungrouped$coefs_vcov, c(0, 0, 0, 1)) # Comparision of Beta1 after
compute_all(gee_mods_ungrouped$coefs_vcov, c(0, -1, 0, 1)) # Comparision of Beta1 after - before
compute_all(gee_mods_ungrouped$coefs_vcov, c(-1, -log(60), 1, log(60))) # Comparision of log odds ratio

##### Results of grouped analysis #####
lapply(gee_mods_grouped, summary) # GEE summary
compute_all(gee_mods_grouped$coefs_vcov, c(0, 1, 0, 0)) # Comparision of Beta1 before
compute_all(gee_mods_grouped$coefs_vcov, c(0, 0, 0, 1)) # Comparision of Beta1 after
compute_all(gee_mods_grouped$coefs_vcov, c(0, -1, 0, 1)) # Comparision of Beta1 after - before
compute_all(gee_mods_grouped$coefs_vcov, c(-1, -log(60), 1, log(60))) # Comparision of log odds ratio

##### Sum of success by group condition with welch anova ###
anova_analysis_sumofsuccess$aov$gr[[1]]; build_multi_ci(anova_analysis_sumofsuccess$data$gr, 80)
anova_analysis_sumofsuccess$aov$gr[[2]]; build_multi_ci(anova_analysis_sumofsuccess$data$gr, 90)
anova_analysis_sumofsuccess$aov$gr[[3]]; build_multi_ci(anova_analysis_sumofsuccess$data$gr, 100)

t.test(success ~ cond, data = anova_analysis_sumofsuccess$data$gr[maxt == 80 & cond != "individual"])
t.test(success ~ cond, data = anova_analysis_sumofsuccess$data$gr[maxt == 90 & cond != "individual"])
t.test(success ~ cond, data = anova_analysis_sumofsuccess$data$gr[maxt == 100 & cond != "individual"])

##### Sum of success by group condition with model results ###
expected_sum_of_success_gee(80, gee_mods_ungrouped)
expected_sum_of_success_gee(90, gee_mods_ungrouped)
expected_sum_of_success_gee(100, gee_mods_ungrouped)

expected_sum_of_success_gee(80, gee_mods_grouped)
expected_sum_of_success_gee(90, gee_mods_grouped)
expected_sum_of_success_gee(100, gee_mods_grouped)

##### Plot Results #####
gee_mods <- gee_mods_ungrouped # gee_mods_grouped
predictions <- testData[,.(iscorrectGr = mean(iscorrectInd)), by = .(round, group_num, cond)]
predictions[cond == "individual", predictions := predict(gee_mods$ind, type = "response")]
predictions[cond == "nonfine", predictions := predict(gee_mods$nonfine, type = "response")]
predictions[cond == "fine", predictions := predict(gee_mods$fine, type = "response")]

mean_vs_predicted_plts <- lapply(c("props", "odds", "logodds"),
                                 function(scale) plot_diff_bingee(predictions, scale))

predictions2 <- testData[,.(iscorrectGr = mean(iscorrectInd)), by = .(round, group_num, cond)
                         ][
                           ,.(average = mean(iscorrectGr)), by = .(round, cond)]

predictions2[cond == "individual", c("pred", "sd") := get_preds_and_sd(gee_mods$ind)]
predictions2[cond == "nonfine", c("pred", "sd") := get_preds_and_sd(gee_mods$nonfine)]
predictions2[cond == "fine", c("pred", "sd") := get_preds_and_sd(gee_mods$fine)]
predictions2[, `:=` (lower = pred - qnorm(0.975)*sd, upper = pred + qnorm(0.975)*sd)]
predictions2[, `:=` (pred = plogis(pred), sd = NULL, lower = plogis(lower), upper = plogis(upper))]

mean_vs_predci_plt <- 
  ggplot(predictions2, aes(x = round)) +
  geom_ribbon(aes(ymin = lower, ymax = upper),
              alpha = 0.33, col = NA, fill = "#00CED1") +
  geom_line(aes(y = average), col = "#FF6347", size = 1) + 
  facet_grid(cond ~ .) + 
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + 
  labs(title = "Empiric Mean and 95% CI of Predicted Value",
       x = "Round", y = "Probability of Success")

inference_plots <- list(mean_vs_predicted_plts = mean_vs_predicted_plts,
                        mean_vs_predci_plt = mean_vs_predci_plt)
rm(gee_mods, predictions, predictions2, mean_vs_predicted_plts, mean_vs_predci_plt)