source("code/02_get_data.R")

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
compute_zvals(diffs_boot$diffs, diffs_boot$var_diffs)
tt

link <- paste0("Data/binomial_bagg_B", B, "_", Sys.time(), ".RData")
link <- gsub(":", "-", link)
link <- gsub(" ", "_", link)
save.image(link)
