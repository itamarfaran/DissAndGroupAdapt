source("code/02_get_data.R")

##### Build functions #####

# Print P_vals
print_pvals <- function(pvals, B, p.adjust.method = c("none", "BH")){
  max_se <- max(sqrt(pvals*(1 - pvals)/B))
  pval_matrix <- sapply(p.adjust.method, function(m) p.adjust(pvals, method = m))
  message("P-value Matrix:")
  print(pval_matrix)
  message(paste("Maximum SE of", max_se))
}

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
sample_from_groups <- function(dt){

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
  
  return(gee_mod$coefficients)
}

# Re-assign participants to groups and re-assign groups to conditions and average over B calls of sample_from_groups
reassign_groups_and_cond <- function(dt, B){
  subnums <- copy(dt[,unique(subnum)])
  groupnums <- copy(dt[,unique(group_num)])
  
  # Randomly assign groups to new conditions
  group_aloc <- data.table(cond = sample(rep(c("individual", "nonfine", "fine"), times = 30), length(groupnums)),
                           group_num = groupnums)[order(cond)]
  # Repeat non-ind conditioned groups 3 times
  group_aloc <- data.table(group_num = c(group_aloc[cond == "individual", group_num],
                                         rep(group_aloc[cond != "individual", group_num], each = 3)))[group_aloc, on = "group_num"]
  # Randomly assign participants to groups
  group_aloc[,subnum := sample(subnums, group_aloc[,.N], T)]
  
  # Join participants to their actual observations
  dt_new <- merge(dt[,.(round, subnum, iscorrectInd)], group_aloc, by = "subnum", allow.cartesian = TRUE)
  dt_new <- prepare_data(dt_new)
  
  # Run sample_from_groups B times and simple-average on the results
  res <- colMeans(do.call(
    rbind, parallel::mclapply(1:B, function(b) sample_from_groups(dt = dt_new),
                              mc.cores = 1)
    ))

  return(res)
}

##### Bootstrap #####

rounds <- c(60, 60)
B_sample_groups <- 100
B_bootstrap <- 1000

tmp_dt <- prepare_data(testData2[,.(round, subnum, group_num, cond, iscorrectInd)])
real_res <- colMeans(do.call(
  rbind, parallel::mclapply(X = 1:B_sample_groups,
                            FUN = function(b) sample_from_groups(dt = tmp_dt),
                            mc.cores = 
                              ifelse(.Platform$OS.type == "windows", 1, parallel::detectCores() - 2)
                    )))

tt <- Sys.time()
boot_res <- do.call(
  rbind, pbmcapply::pbmclapply(X = 1:B_bootstrap,
                               FUN = function(k) reassign_groups_and_cond(dt = testData2,
                                                                          B = B_sample_groups),
                               mc.cores =
                                 ifelse(.Platform$OS.type == "windows", 1, parallel::detectCores() - 2)
                               ))
tt <- Sys.time() - tt

##### Compute absolute diffs (beta1g - beta1g) #####
diffs_of_coefs_real <- c(
  real_res[4] - real_res[2], # ind
  real_res[8] - real_res[6], # nonfine
  real_res[12] - real_res[10] # fine
)
diffs_of_diffs_real <- c(
  diffs_of_coefs_real[3] - diffs_of_coefs_real[1], # fine - ind
  diffs_of_coefs_real[3] - diffs_of_coefs_real[2], # fine - nonfine
  diffs_of_coefs_real[2] - diffs_of_coefs_real[1] # nonfine - ind
)
names(diffs_of_diffs_real) <- c("fine_ind", "fine_nonefine", "nonfine_ind")

diffs_of_coefs_boot <- cbind(
  boot_res[,4] - boot_res[,2], # ind
  boot_res[,8] - boot_res[,6], # nonfine
  boot_res[,12] - boot_res[,10] # fine
)
diffs_of_diffs_boot <- cbind(
  diffs_of_coefs_boot[,3] - diffs_of_coefs_boot[,1], # fine - ind
  diffs_of_coefs_boot[,3] - diffs_of_coefs_boot[,2], # fine - nonfine
  diffs_of_coefs_boot[,2] - diffs_of_coefs_boot[,1] # nonfine - ind
)
colnames(diffs_of_diffs_boot) <- c("fine_ind", "fine_nonefine", "nonfine_ind")

pvals_diffs <- sapply(1:3, function(i) mean(abs(diffs_of_diffs_boot[,i]) > abs(diffs_of_diffs_real[i])))
names(pvals_diffs) <- names(diffs_of_diffs_real)

##### Compute log odds at time t #####
log_odds_ratio_real <- 
  sapply(c(1, 5, 9),
         function(i) real_res[i + 2] - real_res[i] + real_res[i + 3]*log(rounds[2]) - real_res[i + 1]*log(rounds[1]))
names(log_odds_ratio_real) <- c("ind", "nonfine", "fine")

all_diffs_real <- c(
  log_odds_ratio_real["fine"] - log_odds_ratio_real["ind"],
  log_odds_ratio_real["fine"] - log_odds_ratio_real["nonfine"],
  log_odds_ratio_real["nonfine"] - log_odds_ratio_real["ind"]
)
names(all_diffs_real) <- c("fine_ind", "fine_nonefine", "nonfine_ind")

log_odds_ratio_boot <- 
  sapply(c(1, 5, 9),
         function(i) boot_res[,i + 2] - boot_res[,i] + boot_res[,i + 3]*log(rounds[2]) - boot_res[,i + 1]*log(rounds[1]))
colnames(log_odds_ratio_boot) <- c("ind", "nonfine", "fine")

all_diffs_boot <- cbind(
  log_odds_ratio_boot[,"fine"] - log_odds_ratio_boot[,"ind"],
  log_odds_ratio_boot[,"fine"] - log_odds_ratio_boot[,"nonfine"],
  log_odds_ratio_boot[,"nonfine"] - log_odds_ratio_boot[,"ind"]
)
colnames(all_diffs_boot) <- c("fine_ind", "fine_nonefine", "nonfine_ind")

pvals_logodds <- sapply(1:3, function(i) mean(abs(all_diffs_boot[,i]) > abs(all_diffs_real[i])))
names(pvals_logodds) <- names(all_diffs_real)

##### Show and save results #####

print_pvals(pvals_diffs, B_bootstrap)
print_pvals(pvals_logodds, B_bootstrap)
tt

link <- paste0("Data/permute_run_", Sys.time(), ".RData")
link <- gsub(":", "-", link)
link <- gsub(" ", "_", link)
save.image(link)
###