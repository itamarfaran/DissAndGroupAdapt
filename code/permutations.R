source("code/analysis_tmp.R")

# Build coefficients for group slopes and intercepts
prepare_data <- function(dt){
  dt_new <- copy(dt[,.(round, subnum, group_num, cond, iscorrectInd)])
  
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

# Get data, subset randomly one user from each group, and compute all differences in log-odds-ratio
sample_from_groups <- function(dt, return_logodds = FALSE, round = c(60, 61)){
  if(length(round) == 1) round <- rep(round, 2)
  
  # sort users by groups
  group_num_indexes <- copy(
    unique(dt[,.(subnum, group_num)])[
      ,.(sub1 = round(nlargest(subnum, 1)),
         sub2 = round(nlargest(subnum, 2)),
         sub3 = round(nlargest(subnum, 3))),
      by = group_num]) 
  
  # Randomly select a user from each group
  which_obs_gee <- sample(1:3, group_num_indexes[,.N], T)
  which_obs_gee <- sapply(1:length(which_obs_gee),
                          function(i) as.matrix(group_num_indexes)[i, 1 + which_obs_gee[i]])
  
  # Build GEE model for subset of data without group interactions
  gee_mod <- geeglm(iscorrectInd ~ 0 +
                      intercept_ind_before + slopelog_ind_before +
                      intercept_ind_after + slopelog_ind_after + 
                      intercept_nonfine_before + slopelog_nonfine_before +
                      intercept_nonfine_after + slopelog_nonfine_after +
                      intercept_fine_before + slopelog_fine_before +
                      intercept_fine_after + slopelog_fine_after,
                    family = binomial(link = "logit"),
                    data = dt, subset = (subnum %in% which_obs_gee),
                    id = subnum, corstr = "ar1")
  
  if(!return_logodds) return(gee_mod$coefficients)
  
  # Build log odds ratio at rounds t by condition
  log_odds_ind <- 
    gee_mod$coefficients["intercept_ind_after"] +
    gee_mod$coefficients["slopelog_ind_after"]*log(round[2]) - 
    gee_mod$coefficients["intercept_ind_before"] -
    gee_mod$coefficients["slopelog_ind_before"]*log(round[1])
  
  log_odds_nonfine <- 
    gee_mod$coefficients["intercept_nonfine_after"] +
    gee_mod$coefficients["slopelog_nonfine_after"]*log(round[2]) - 
    gee_mod$coefficients["intercept_nonfine_before"] -
    gee_mod$coefficients["slopelog_nonfine_before"]*log(round[1])
  
  log_odds_fine <- 
    gee_mod$coefficients["intercept_fine_after"] +
    gee_mod$coefficients["slopelog_fine_after"]*log(round[2]) - 
    gee_mod$coefficients["intercept_fine_before"] -
    gee_mod$coefficients["slopelog_fine_before"]*log(round[1])
  
  # Return all differences in log odds ratio
  res <- 
    c(log_odds_fine - log_odds_ind,
      log_odds_fine - log_odds_nonfine,
      log_odds_nonfine - log_odds_ind)
  names(res) = c("fine_ind", "fine_nonefine", "nonfine_ind")
  return(res)
}

# Re-assign participants to groups and re-assign groups to conditions and average over B calls of sample_from_groups
reassign_groups_and_cond <- function(dt, return_logodds = FALSE, round = c(60, 61), B){
  subnums <- copy(dt[,unique(subnum)])
  groupnums <- copy(dt[,unique(group_num)])
  
  # Randomly assign groups to new conditions
  group_aloc <- data.table(cond = sample(rep(c("individual", "nonfine", "fine"), times = 30), length(groupnums)),
                           group_num = groupnums)[order(cond)]
  # Repeat non-ind conditioned groups 3 times
  group_aloc <- data.table(group_num = c(group_aloc[cond == "individual", group_num],
                                         rep(group_aloc[cond != "individual", group_num], each = 3)))[group_aloc, on = "group_num"]
  # Randomly assign participants to groups
  group_aloc[,subnum := sample(subnums, length(subnums))]
  
  # Join participants to their actual observations
  dt_new <- copy(group_aloc[dt[,.(round, subnum, iscorrectInd)], on = "subnum"])
  dt_new <- prepare_data(dt_new)
  
  # Run sample_from_groups B times and simple-average on the results
  res <- colMeans(do.call(
    rbind, parallel::mclapply(1:B, function(b) sample_from_groups(dt = dt_new,
                                                                  return_logodds = return_logodds,
                                                                  round = round),
                              mc.cores = 1)
    ))

  return(res)
}

rounds <- c(60, 61)
return_logodds <- FALSE
B_sample_groups <- 100
B_bootstrap <- 1000

tmp_dt <- prepare_data(testData2[,.(round, subnum, group_num, cond, iscorrectInd)])
real_res <- colMeans(do.call(
  rbind, parallel::mclapply(X = 1:B_sample_groups,
                            FUN = function(b) sample_from_groups(dt = tmp_dt,
                                                                 return_logodds = return_logodds,
                                                                 round = rounds),
                            mc.cores = 
                              ifelse(.Platform$OS.type == "windows", 1, parallel::detectCores() - 2)
                    )))

tt <- Sys.time()
boot_res <- do.call(
  rbind, pbmcapply::pbmclapply(X = 1:B_bootstrap,
                               FUN = function(k) reassign_groups_and_cond(dt = testData2,
                                                                          return_logodds = return_logodds,
                                                                          round = rounds,
                                                                          B = B_sample_groups),
                               mc.cores =
                                 ifelse(.Platform$OS.type == "windows", 1, parallel::detectCores() - 2)
                               ))
tt <- Sys.time() - tt

if(return_logodds){
  pvals <- sapply(1:3, function(i) mean(abs(boot_res[,i]) > abs(real_res[i])))
  names(pvals) <- names(real_res)
} else {
  diffs_of_coefs_real <- c(
    real_res[4] - real_res[2],
    real_res[8] - real_res[6],
    real_res[12] - real_res[10]
  )
  diffs_of_diffs_real <- c(
    diffs_of_coefs_real[3] - diffs_of_coefs_real[1],
    diffs_of_coefs_real[3] - diffs_of_coefs_real[2],
    diffs_of_coefs_real[2] - diffs_of_coefs_real[1]
  )
  names(diffs_of_diffs_real) <- c("fine_ind", "fine_nonefine", "nonfine_ind")
  
  diffs_of_coefs_boot <- cbind(
    boot_res[,4] - boot_res[,2],
    boot_res[,8] - boot_res[,6],
    boot_res[,12] - boot_res[,10]
  )
  diffs_of_diffs_boot <- cbind(
    diffs_of_coefs_boot[,3] - diffs_of_coefs_boot[,1],
    diffs_of_coefs_boot[,3] - diffs_of_coefs_boot[,2],
    diffs_of_coefs_boot[,2] - diffs_of_coefs_boot[,1]
  )
  colnames(diffs_of_diffs_boot) <- c("fine_ind", "fine_nonefine", "nonfine_ind")
  
  pvals <- sapply(1:3, function(i) mean(abs(boot_res[,i]) > abs(real_res[i])))
  names(pvals) <- names(diff_real_res)
}

max(sqrt(pvals*(1 - pvals)/B_bootstrap)) # max se
p.adjust(pvals, method = "none")
p.adjust(pvals, method = "BH")
tt

link <- paste0("Data/permute_run_", Sys.time(), ".RData")
link <- gsub(":", "-", link)
link <- gsub(" ", "_", link)
save.image(link)
###