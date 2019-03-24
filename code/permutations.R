source("code/00_baseFunctions.R")
source("code/analysis_tmp.R")

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
sample_from_groups <- function(dt, round){
  if(length(round) == 1) round <- rep(round, 2)
  
  group_num_indexes <- copy(
    unique(dt[,.(subnum, group_num)])[
      ,.(sub1 = round(nlargest(subnum, 1)),
         sub2 = round(nlargest(subnum, 2)),
         sub3 = round(nlargest(subnum, 3))),
      by = group_num])
  group_num_indexes[sub2 == 0 | sub3 == 0,`:=` (sub2 = sub1, sub3 = sub1)]
  
  which_obs_gee <- sample(1:3, group_num_indexes[,.N], T)
  which_obs_gee <- sapply(1:length(which_obs_gee),
                          function(i) as.matrix(group_num_indexes)[i, 1 + which_obs_gee[i]])
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
  
  res <- 
    c(log_odds_fine - log_odds_ind,
      log_odds_fine - log_odds_nonfine,
      log_odds_nonfine - log_odds_ind)
  names(res) = c("fine_ind", "fine_nonefine", "nonfine_ind")
  return(res)
}
reassign_groups_and_cond <- function(dt, round, B){
  subnums <- copy(dt[,unique(subnum)])
  groupnums <- copy(dt[,unique(group_num)])
  group_aloc <- data.table(cond = sample(rep(c("individual", "nonfine", "fine"), times = 30), length(groupnums)),
                           group_num = groupnums)[order(cond)]
  group_aloc <- data.table(group_num = c(group_aloc[cond == "individual", group_num],
                                         rep(group_aloc[cond != "individual", group_num], each = 3)))[group_aloc, on = "group_num"]
  group_aloc[,subnum := sample(subnums, length(subnums))]
  
  dt_new <- copy(group_aloc[testData2[,.(round, subnum, iscorrectInd)], on = "subnum"])
  dt_new <- prepare_data(dt_new)

  res <- try(colMeans(t(sapply(1:B, function(b) sample_from_groups(dt_new, round)))))
  if(class(res) != "numeric") res <- rep(NA, 3)
  
  return(res)
}

rounds <- 60
B_sample_groups <- 100
B_bootstrap <- 1000

testData3 <- prepare_data(testData2[,.(round, subnum, group_num, cond, iscorrectInd)])
real_res <- colMeans(t(sapply(1:B_sample_groups, function(b) sample_from_groups(testData3, 60))))

boot_res <- do.call(rbind, 
                    parallel::mclapply(X = 1:B_bootstrap,
                                       FUN = function(k) reassign_groups_and_cond(testData2, rounds, 1),
                                       mc.cores = 
                                         ifelse(.Platform$OS.type == "windows", 1, parallel::detectCores() - 2)
                                       ))

pvals <- sapply(1:3, function(i) mean(abs(boot_res[,i]) > abs(real_res[i])))
p.adjust(pvals, method = "BH")
###

