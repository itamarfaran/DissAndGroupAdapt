source("code/02_get_data.R")

gee_mods_ungrouped <- boot_geeglm(testData2, do.boot = FALSE, return.model = TRUE)
gee_mods_ungrouped$coefs_vcov <- get_coefs_vcov(lst = gee_mods_ungrouped)

testData2_grouped <- prepare_data(unique(testData2[,.(round, group_num, cond, iscorrectGr)]))
gee_mods_grouped <- list(
  ind = 
    geeglm(iscorrectGr ~ 0 + intercept_ind_before + slopelog_ind_before + 
             intercept_ind_after + slopelog_ind_after,
           data = testData2_grouped, family = binomial(link = "logit"),
           id = group_num, corstr = "ar1", subset = (cond == "individual")),
  nonfine = 
    geeglm(iscorrectGr ~ 0 + intercept_nonfine_before + slopelog_nonfine_before + 
             intercept_nonfine_after + slopelog_nonfine_after,
           data = testData2_grouped, family = binomial(link = "logit"),
           id = group_num, corstr = "ar1", subset = (cond == "nonfine")),
  fine = 
    geeglm(iscorrectGr ~ 0 + intercept_fine_before + slopelog_fine_before + 
             intercept_fine_after + slopelog_fine_after,
           data = testData2_grouped, family = binomial(link = "logit"),
           id = group_num, corstr = "ar1", subset = (cond == "fine"))
  )
gee_mods_grouped$coefs_vcov <- get_coefs_vcov(lst = gee_mods_grouped)

##### Sum of Successes T-test #####

t_stop_time <- c(80, 90, 100)

sum_success_dt_by_ind <-
  do.call(rbind,
          lapply(t_stop_time, function(t)
            testData2[between(round, 61, t)
                      ,.(N = .N, success = sum(iscorrectInd)),
                      by = .(round, group_num, cond)
                      ][,.(maxt = t, success = sum(success*3/N)),
                        by = .(group_num, cond)]
                 ))

sum_success_dt_by_group <-
  do.call(rbind,
          lapply(t_stop_time, function(t)
            unique(testData2[between(round, 61, t), .(round, group_num, cond, iscorrectGr)])[
              ,.(maxt = t, success = sum(iscorrectGr)), by = .(group_num, cond)]
          ))

aov_models_by_ind <- 
  lapply(t_stop_time,
         function(t) oneway.test(success ~ cond,
                         data = sum_success_dt_by_ind[maxt == t]
         ))

aov_models_by_group <- 
  lapply(t_stop_time,
         function(t) oneway.test(success ~ cond,
                                 data = sum_success_dt_by_group[maxt == t]
         ))

anova_analysis_sumofsuccess <- list(
  aov = list(ind = aov_models_by_ind, gr = aov_models_by_group),
  data = list(ind = sum_success_dt_by_ind, gr = sum_success_dt_by_group)
)
rm(aov_models_by_ind, aov_models_by_group,
   sum_success_dt_by_ind, sum_success_dt_by_group)
rm(t_stop_time)