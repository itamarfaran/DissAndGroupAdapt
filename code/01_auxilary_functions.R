ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

ipak(c("parallel", "pbmcapply", "readr", "tidyverse", "data.table", "dtplyr", "pracma",
       "geepack", "plotly", "Matrix", "mvtnorm"))
# ipak(c("corrplot", "msm", "markovchain", "lme4", "glmmTMB"))

createID <- function(x, factor = TRUE){
  x <- factor(x)
  x2 <- numeric(length(x))
  xlev <- levels(x)
  for(k in xlev){
    x2[x == k] <- which(xlev == k)
  }
  if(factor) return(factor(x2))
  return(x2)
}
makeAllNumeric <- function(df){
  for(i in 1:ncol(df)) df[,i] <- as.numeric(df[,i])
  return(df)
}
auto_cor <- Vectorize(function(x, lag) cor(x[1:(length(x) - lag)], x[(1 + lag):length(x)]), "lag")
nlargest <- function(x, k){
  n <- length(x)
  if(k > n) return(min(x))
  return(sort(x, decreasing = TRUE)[k])
}
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
build_multi_ci <- function(dt, t, sig.level = 0.05, correct = TRUE){
  means0 <- dt[maxt == t, .(mean = mean(success)), keyby = cond]
  means <- means0$mean
  names(means) <- means0$cond
  
  vars0 <- dt[maxt == t, .(var = var(success)/.N), keyby = cond]
  vars <- diag(vars0$var)
  rownames(vars) <- colnames(vars) <- vars0$cond
  
  C <- matrix(c(
    1, -1, 0,
    1, 0, -1,
    0, 1, -1),
    ncol = 3, byrow = T)
  
  
  diff <- as.vector(C %*% means)
  diff_vars <- C %*% vars %*% t(C)
  quantile_t <- qmvt(1 - sig.level/2,
                     df = 90 - 6,
                     corr = cov2cor(diff_vars))$quantile
  
  names(diff) <- c(
    paste0(names(means)[1], "_", names(means)[2]),
    paste0(names(means)[1], "_", names(means)[3]),
    paste0(names(means)[2], "_", names(means)[3])
  )
  
  
  cis <- (diff %o% rep(1, 3)) + (( quantile_t*sqrt(diag(diff_vars)) ) %o% c(-1, 0, 1))
  colnames(cis) <- c(paste0("lower", 100*(1 - sig.level), "%"),
                     "estimate",
                     paste0("upper", 100*(1 - sig.level), "%"))
  return(cis)
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
get_coefs_vcov <- function(..., lst){
  models <- if(missing(lst)) list(...) else lst
  coeffs <- simplify(purrr::transpose(models)$coefficients)
  vcov <- as.matrix(bdiag(lapply(models, function(m) summary(m)$cov.scaled)))
  
  return(list(coeffs = coeffs, vcov = vcov))
}
get_preds_and_sd <- function(mod){
  X <- rbind(
    cbind(rep(1, 60), log(1:60), rep(0, 60), rep(0, 60)),
    cbind(rep(0, 40), rep(0, 40), rep(1, 40), log(61:100))
  )
  
  return(list(
    X %*% mod$coefficients,
    sqrt(diag(X %*% summary(mod)$cov.scaled %*% t(X) ))
  ))
}
expected_sum_of_success_gee <- function(t, gee_mods){
  leng_sum <- t - 60 
  
  C1 <- matrix(c(1, -1, 0, 1, 0, -1, 0, 1, -1), ncol = 3, byrow = TRUE)
  
  C2 <- rbind(c(rep(1, leng_sum), rep(0, leng_sum), rep(0, leng_sum)),
              c(rep(0, leng_sum), rep(1, leng_sum), rep(0, leng_sum)),
              c(rep(0, leng_sum), rep(0, leng_sum), rep(1, leng_sum)))
  
  C3_t1 <- cbind(rep(1, leng_sum), log(61:t))
  C3_t2 <- matrix(0, nr = leng_sum, nc = 2)
  C3 <- rbind(cbind(C3_t1, C3_t2, C3_t2),
              cbind(C3_t2, C3_t1, C3_t2),
              cbind(C3_t2, C3_t2, C3_t1))
  
  C <- C1 %*% C2 %*% C3
  
  beta_ind <- c(3, 4, 7, 8, 11, 12)
  beta_coef <- gee_mods$coefs_vcov$coeffs[beta_ind]
  beta_var <- gee_mods$coefs_vcov$vcov[beta_ind, beta_ind]
  
  contrast_coef <- C %*% beta_coef
  contrast_var <- C %*% beta_var %*% t(C)
  
  names(contrast_coef) <- rownames(contrast_var) <- colnames(contrast_var) <-
    c("ind_nonfine", "ind_fine", "nonfine_fine")
  
  chi <- as.vector(t(contrast_coef[1:2]) %*% solve(contrast_var[1:2, 1:2]) %*% contrast_coef[1:2])
  
  sum_of_success_diff_mod <- list(
    coef = contrast_coef,
    vcov = contrast_var,
    chi = chi,
    p_chi = pchisq(chi, 2, lower.tail = F),
    pairwise = compute_zvals(contrast_coef, contrast_var)[-(4:6),]
  )
  
  return(sum_of_success_diff_mod)
}

compute_log_odds <- function(coef, vcov, contrast){
  l <- length(contrast)
  C1 <- rbind(c(contrast, rep(0, l), rep(0, l)),
              c(rep(0, l), contrast, rep(0, l)),
              c(rep(0, l), rep(0, l), contrast))
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
  sds <- sqrt(diag(vcov))
  q <- ifelse(adjust.ci,
              mvtnorm::qmvnorm(1 - sig.level/2, corr = cov2cor(vcov))$quantile,
              qnorm(1 - sig.level/2))
  z <- est/sds
  p <- 2*pnorm(abs(z), lower.tail = FALSE)
  p.adj <- p.adjust(p, p.adjust.method)
  
  res <- cbind(est, est - q*sds, est + q*sds,
               exp(est), exp(est - q*sds), exp(est + q*sds),
               sds, z, p, p.adj)
  
  colnames(res) <- c("Estimate",
                     paste0("Lower", (1 - sig.level)*100, "%", ifelse(adjust.ci, "(", "(Not "), "Adjusted)"),
                     paste0("Upper", (1 - sig.level)*100, "%", ifelse(adjust.ci, "(", "(Not "), "Adjusted)"),
                     "exp[Estimate]",
                     paste0("Lower", (1 - sig.level)*100, "%", ifelse(adjust.ci, "(", "(Not "), "Adjusted) (exp)"),
                     paste0("Upper", (1 - sig.level)*100, "%", ifelse(adjust.ci, "(", "(Not "), "Adjusted) (exp)"),
                     "SD", "Z-Value", "P-Value",
                     paste0("Adj.P-Value(", p.adjust.method, ")"))
  rownames(res) <- str_replace(rownames(res), "_", " - ")
  rownames(res) <- tools::toTitleCase(rownames(res))
  return(t(res))
}
compute_all <- function(coefs_vcov_lst, contrast, sig.level = 0.05, adjust.ci = TRUE, p.adjust.method = "BH"){
  contrast_per_group <- compute_log_odds(coefs_vcov_lst$coeffs,
                                         coefs_vcov_lst$vcov,
                                         contrast)
  diffs_contrast <- compute_all_contrasts(contrast_per_group$log_odds, contrast_per_group$var_log_odds)
  chisq_val <- t(diffs_contrast$diffs[1:2]) %*% solve(diffs_contrast$var_diffs[1:2, 1:2]) %*% diffs_contrast$diffs[1:2]
  
  pval_chisq <- pchisq(chisq_val, 2, lower.tail = FALSE)
  zvals <- compute_zvals(diffs_contrast$diffs, diffs_contrast$var_diffs,
                         sig.level = sig.level, adjust.ci = adjust.ci, p.adjust.method = p.adjust.method)
  
  return(list(chisq = c("Chisq Value" = chisq_val, "DF" = 2, "P-Value" = pval_chisq),
              z = zvals))
}

plot_diff_bingee <- function(dt, scale = c("props", "odds", "logodds")){
  scale <- scale[1]
  if(!scale %in% c("props", "odds", "logodds")) stop("'scale' must be one of 'props', 'odds', 'logodds'")
  scale_fn <- switch(scale,
                     "props" = mean,
                     "odds" = function(x) mean(x)/(1 - mean(x)),
                     "logodds" = function(x) qlogis(mean(x)))
  scale_title <- switch(scale,
                        "props" = "Proportion ",
                        "odds" = "Odds",
                        "logodds" = "Log Odds")
  
  errors_dt <- dt[,.(diff = scale_fn(iscorrectGr) - scale_fn(predictions)),
                  by = .(round, cond)][,diff]
  
  error_metrics <- list(
    mean_error = mean(errors_dt),
    rmse = sqrt(mean(errors_dt^2)),
    median_error = median(errors_dt),
    mad_ = median(abs(errors_dt))
  )
  
  cat("Summary of errors:\n")
  print(summary(errors_dt))
  cat(paste0("RMSE: ", round(error_metrics$rmse, 5), " || MAD: ", round(error_metrics$mad_, 5)))
  
  dt1 <- dt[,.(emp_mean = scale_fn(iscorrectGr),
               predicted_mean = scale_fn(predictions)),
            by = .(round, cond)]
  dt1[, `:=` (error = emp_mean - predicted_mean)]
  dt1[, `:=` (lower = predicted_mean + error*(error < 0),
              upper = predicted_mean + error*(error > 0)
  )]
  vs <-
    ggplot(dt1, aes(x = round, col = cond)) +
    geom_line(aes(y = predicted_mean), size = 1) +
    geom_vline(xintercept = 0) +
    geom_linerange(aes(ymin = lower, ymax = upper), alpha = 0.4) + 
    geom_point(aes(y = emp_mean), shape = 18, alpha = 0.6) + 
    labs(title = paste0("Empiric Mean vs. Predicted Value (", scale_title, " scale)"),
         x = "Round", y = paste0(scale_title, " of Success"), col = "Condition", type = "Type")
  
  vs <- vs + if(scale == "odds") scale_y_log10() else geom_hline(yintercept = 0)
  
  diff <- dt[,.(diff = scale_fn(iscorrectGr) - scale_fn(predictions)), by = .(round, cond)] %>%
    ggplot(aes(x = round, y = movavg(diff, 3, "e"), col = cond)) +
    geom_line(size = 1) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
    stat_summary(aes(group = 1), fun.y = mean, geom = "line", size = 1.2, linetype = 2) +
    labs(title = paste("Difference of Empiric Mean vs. Predicted Value (", scale_title, " scale)"),
         x = "Round", y = "Differnce", col = "Condition", type = "Type")
  
  return(list(vs = vs,
              diff = diff,
              error_stats = error_metrics,
              errors = errors_dt))
}
