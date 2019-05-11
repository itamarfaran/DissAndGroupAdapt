source("code/analysis_tmp.R")

##### Functions #####
prepare_data <- function(dt){
  dt_new <- copy(dt)
  
  dt_new[, `:=` (intercept_ind_before = ifelse(cond == "individual" & round <= 60, 1, 0),
                 slopelog_ind_before = ifelse(cond == "individual" & round <= 60, log(round), 0),
                 intercept_ind_after = ifelse(cond == "individual" & round > 60, 1, 0),
                 slopelog_ind_after = ifelse(cond == "individual" & round > 60, log(round), 0)
  )]
  
  dt_new[, `:=` (intercept_nonfine_before = ifelse(cond == "nonfine" & round <= 60, 1, 0),
                 slopelog_nonfine_before = ifelse(cond == "nonfine" & round <= 60, log(round), 0),
                 intercept_nonfine_after = ifelse(cond == "nonfine" & round > 60, 1, 0),
                 slopelog_nonfine_after = ifelse(cond == "nonfine" & round > 60, log(round), 0)
  )]
  
  dt_new[, `:=` (intercept_fine_before = ifelse(cond == "fine" & round <= 60, 1, 0),
                 slopelog_fine_before = ifelse(cond == "fine" & round <= 60, log(round), 0),
                 intercept_fine_after = ifelse(cond == "fine" & round > 60, 1, 0),
                 slopelog_fine_after = ifelse(cond == "fine" & round > 60, log(round), 0)
  )]
  return(dt_new)
}
boot_geeglm <- function(dt, do.boot = TRUE, even.groups = TRUE, return.model = FALSE, seed){
  cdt <- copy(dt)
  
  cdt <- cdt[,.(success = sum(iscorrectInd), fail = .N - sum(iscorrectInd)),
             by = .(round, group_num, cond)]
  
  if(do.boot){
    if(!missing(seed)) set.seed(seed)
    groups2boot <- if(even.groups){
      sort(as.numeric(
        sapply(cdt[,unique(cond)], function(c) cdt[cond == c, sample(unique(group_num), replace = T)])
      ))
    } else cdt[,sort(sample(unique(group_num), replace = T))]
    cdt_booted <- 
      do.call(rbind,
              lapply(1:length(groups2boot),
                     function(i) cdt[group_num == groups2boot[i],
                                     .(round, group_num = i, cond, success, fail)]))
  } else cdt_booted <- cdt
  cdt2 <- prepare_data(cdt_booted)
  X <- as.matrix(cdt2[,6:ncol(cdt2)])
  y <- as.matrix(cdt2[,.(success, fail)])
  rounds <- cdt2[,round]
  group_num <- cdt2[,group_num]
  
  cond <- cdt2[,cond]
  individuals <- cond == "individual"
  nonfines <- cond == "nonfine"
  fines <- cond == "fine"
  
  gee_mod_ind <- geeglm(y[individuals,] ~ 0 + X[individuals, 1:4],
                        family = binomial(link = "logit"),
                        id = group_num[individuals], corstr = "ar1")
  gee_mod_nonfine <- geeglm(y[nonfines,] ~ 0 + X[nonfines, 5:8],
                            family = binomial(link = "logit"),
                            id = group_num[nonfines], corstr = "ar1")
  gee_mod_fine <- geeglm(y[fines,] ~ 0 + X[fines, 9:12],
                         family = binomial(link = "logit"),
                         id = group_num[fines], corstr = "ar1")
  
  names(gee_mod_ind$coefficients) <-
    paste0(c("intercept_before", "slopelog_before", "intercept_after", "slopelog_after"), "_ind")
  names(gee_mod_nonfine$coefficients) <-
    paste0(c("intercept_before", "slopelog_before", "intercept_after", "slopelog_after"), "_nonfine")
  names(gee_mod_fine$coefficients) <-
    paste0(c("intercept_before", "slopelog_before", "intercept_after", "slopelog_after"), "_fine")
  
  if(return.model) return(list(ind = gee_mod_ind, nonfine = gee_mod_nonfine, fine = gee_mod_fine))
  return(get_coefs_vcov(gee_mod_ind, gee_mod_nonfine, gee_mod_fine))
}
compute_log_odds <- function(coef, vcov, contrast){
  C1 <- rbind(c(contrast, rep(0, 4), rep(0, 4)),
              c(rep(0, 4), contrast, rep(0, 4)),
              c(rep(0, 4), rep(0, 4), contrast))
  rownames(C1) <- c("ind", "nonfine", "fine")
  return(list(log_odds = C1 %*% coef,
              var_log_odds = C1 %*% vcov %*% t(C1)))
}
compute_all_contrasts <- function(log_odds, var_log_odds){
  C2 <- rbind(c(-1, 0, 1),
              c(0, -1, 1),
              c(-1, 1, 0))
  rownames(C2) <- c("fine_ind", "fine_nonfine", "nonfine_ind")
  
  return(list(diffs = C2 %*% log_odds,
              var_diffs = C2 %*% var_log_odds %*% t(C2)))
}
compute_zvals <- function(est, vcov, sig.level = 0.05, adjust.ci = TRUE, p.adjust.method = "BH"){
  est
  sds <- sqrt(diag(vcov))
  q <- ifelse(adjust.ci,
              mvtnorm::qmvnorm(1 - sig.level/2, corr = cov2cor(vcov))$quantile,
              qnorm(1 - sig.level/2))
  z <- est/sds
  p <- 2*pnorm(abs(z), lower.tail = FALSE)
  p.adj <- p.adjust(p, p.adjust.method)
  
  res <- cbind(est, sds, est - q*sds, est + q*sds, z, p, p.adj)
  
  colnames(res) <- c("Estimate", "SD",
                     paste0("Lower", (1 - sig.level)*100, "%", ifelse(adjust.ci, "(", "(Not "), "Adjusted)"),
                     paste0("Upper", (1 - sig.level)*100, "%", ifelse(adjust.ci, "(", "(Not "), "Adjusted)"),
                     "Z-Value", "P-Value",
                     paste0("Adj.P-Value(", p.adjust.method, ")"))
  rownames(res) <- str_replace(rownames(res), "_", " - ")
  rownames(res) <- tools::toTitleCase(rownames(res))
  return(t(res))
}
get_coefs_vcov <- function(..., lst){
  models <- if(missing(lst)) list(...) else lst
  coeffs <- simplify(purrr::transpose(models)$coefficients)
  vcov <- as.matrix(bdiag(lapply(models, function(m) summary(m)$cov.scaled)))
  
  return(list(coeffs = coeffs, vcov = vcov))
}
plot_diff_bingee <- function(dt, scale = c("props", "odds", "logodds")){
  scale <- scale[1]
  if(!scale %in% c("props", "odds", "logodds")) stop("'scale' must be one of 'props', 'odds', 'logodds'")
  scale_fn <- switch(scale,
                     "props" = mean,
                     "odds" = function(x) mean(x)/(1 - mean(x)),
                     "logodds" = function(x) qlogis(mean(x)))
  scale_title <- switch(scale,
                        "props" = "(Proportion Scale)",
                        "odds" = "(Odds Scale)",
                        "logodds" = "(Log Odds Scale)")
  
  errors_dt <- dt[,.(diff = scale_fn(iscorrectGr) - scale_fn(predictions)),
                  by = .(round, cond)][,diff]
  
  error_metrics <- list(
    mean_error = mean(errors_dt),
    rmse = sqrt(mean(errors_dt^2)),
    median_error = median(errors_dt),
    mad_ = median(abs(errors_dt))
  )
  
  message("Summary of errors:")
  print(summary(errors_dt))
  cat(paste0("RMSE: ", round(error_metrics$rmse, 5), " || MAD: ", round(error_metrics$mad_, 5)))
  
  vs <- dt[,.(emp_mean = scale_fn(iscorrectGr),
              predicted_mean = scale_fn(predictions)),
           by = .(round, cond)] %>%
    gather(type, value, -round, -cond) %>%
    ggplot(aes(x = round, y = movavg(value, 4, "e"), col = cond, linetype = type)) +
    geom_line(size = 1) + geom_vline(xintercept = 0) +
    labs(title = paste("Empiric Mean vs. Predicted Value", scale_title),
         x = "Round", y = "Probability of Success", col = "Condition", type = "Type")
  
  vs <- vs + if(scale == "odds") scale_y_log10() else geom_hline(yintercept = 0)
  
  diff <- predictions[,.(diff = scale_fn(iscorrectGr) - scale_fn(predictions)),
                      by = .(round, cond)] %>%
    ggplot(aes(x = round, y = movavg(diff, 3, "e"), col = cond)) +
    geom_line(size = 1) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
    stat_summary(aes(group = 1), fun.y = mean, geom = "line", size = 1.2, linetype = 2) +
    labs(title = paste("Difference of Empiric Mean vs. Predicted Value", scale_title),
         x = "Round", y = "Differnce", col = "Condition", type = "Type")
  
  return(list(vs = vs,
              diff = diff,
              error_stats = error_metrics,
              errors = errors_dt))
}

##### Real Results #####

gee_mods <- boot_geeglm(testData2, do.boot = FALSE, return.model = TRUE)
lapply(gee_mods, summary)

coefs_vcov <- get_coefs_vcov(lst = gee_mods)

log_odds <- compute_log_odds(coefs_vcov$coeffs,
                             coefs_vcov$vcov,
                             c(-1, -log(60), 1, log(60)))

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

##### Bagging #####

ncores <- ifelse(.Platform$OS.type == "windows", 1, parallel::detectCores() - 2)
B <- 50*ncores

tt <- Sys.time()
res_boot <- pbmclapply(1:B, function(b) boot_geeglm(testData2), mc.cores = ncores)
tt <- Sys.time() - tt

coeffs_boot <- do.call(rbind, purrr::transpose(res_boot)$coeffs)
vcov_boot <- array(unlist(purrr::transpose(res_boot)$vcov),
                   dim = c(ncol(coeffs_boot), ncol(coeffs_boot), B))

drop_obs <- data.table(which(abs(coeffs_boot) > 10^2.5, arr.ind = T))
drop_obs[between(col, 1, 4), dropped := "ind"]
drop_obs[between(col, 5, 8), dropped := "nonfine"]
drop_obs[between(col, 9, 12), dropped := "fine"]
drop_obs_vec <- drop_obs[,sort(unique(row))]

if(length(drop_obs_vec) > 0){
  coeffs_boot_drp <- coeffs_boot[-drop_obs_vec,]
  vcov_boot_drp <- vcov_boot[,,-drop_obs_vec]
} else {
  coeffs_boot_drp <- coeffs_boot
  vcov_boot_drp <- vcov_boot
}
Btag <- B - length(drop_obs_vec)
message(paste("Bootstraps dropped: ", round(1 - Btag/B, 4)*100, "%" ))
unique(drop_obs[,.(row, dropped)])[,.N, by = dropped]

coeffs_boot_mean <- colMeans(coeffs_boot_drp)
vcov_boot_emp <- cov(coeffs_boot_drp)
vcov_boot_mean <- vcov_boot_drp[,,1]/Btag
for(b in 2:Btag) vcov_boot_mean <- vcov_boot_mean + vcov_boot_drp[,,b]/Btag

log_odds_boot <- compute_log_odds(coeffs_boot_mean,
                                  vcov_boot_emp,
                                  # vcov_boot_mean,
                                  c(-1, -log(60), 1, log(60)))
diffs_boot <- compute_all_contrasts(log_odds_boot$log_odds, log_odds_boot$var_log_odds)
compute_zvals(diffs$diffs, diffs$var_diffs)
compute_zvals(diffs_boot$diffs, diffs_boot$var_diffs)
tt

link <- paste0("Data/binomial_bagg_B", B, "_", Sys.time(), ".RData")
link <- gsub(":", "-", link)
link <- gsub(" ", "_", link)
save.image(link)
