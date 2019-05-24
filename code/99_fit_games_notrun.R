source("code/analysis_tmp.R")

##### Real Results #####

gee_mods <- boot_geeglm(testData2, do.boot = FALSE, return.model = TRUE)
lapply(gee_mods, summary)

coefs_vcov <- get_coefs_vcov(lst = gee_mods)

log_odds <- compute_log_odds(coefs_vcov$coeffs,
                             coefs_vcov$vcov,
                             c(-59, -log(60), 59, log(60)))

diffs <- compute_all_contrasts(log_odds$log_odds, log_odds$var_log_odds)
compute_zvals(diffs$diffs, diffs$var_diffs)

##### Plot Results #####

predictions <- testData[,.(iscorrectGr = mean(iscorrectInd)), by = .(round, group_num, cond)]
predictions[cond == "individual", predictions := predict(gee_mods$ind, type = "response")]
predictions[cond == "nonfine", predictions := predict(gee_mods$nonfine, type = "response")]
predictions[cond == "fine", predictions := predict(gee_mods$fine, type = "response")]

p_scale_plt <- plot_diff_bingee(predictions, "props")
odds_scale_plt <- plot_diff_bingee(predictions, "odds")
log_odds_scale_plt <- plot_diff_bingee(predictions, "logodds")

p_scale_plt$vs
p_scale_plt$diff

odds_scale_plt$vs
odds_scale_plt$diff

log_odds_scale_plt$vs
log_odds_scale_plt$diff
