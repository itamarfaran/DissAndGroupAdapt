source("code/02_get_data.R")
# library(glmmTMB)
# library(vars)

nocorr_model <- glm(iscorrectInd ~ log(round)*cond*afterShock, data = testData2, family = binomial(link = "logit"))
summary(nocorr_model)
testData2[,`:=` (predict_link = predict(nocorr_model, type = "response"),
                 residuals = resid(nocorr_model, type = "pearson"))]

do_ma <- TRUE
B <- testData2[,length(unique(subnum))]
arima_coeffs <- matrix(0, nr = B, nc = ifelse(do_ma, 3, 2))
arima_vars <- array(0, dim = c(rep(ifelse(do_ma, 3, 2), 2), B))

colnames(arima_coeffs) <- colnames(arima_vars) <- row.names(arima_vars) <- 
  if(do_ma){ c("ar1", "ma1", "intercept") } else { c("ar1", "intercept") }

sapply(1:B, function(b) {
  tmp <- arima(testData2[subnum == b, residuals], order = c(1,0, 1*do_ma), include.mean = T)
  .GlobalEnv$arima_coeffs[b,] <- tmp$coef
  .GlobalEnv$arima_vars[,,b] <- tmp$var.coef
  return(".")
  })

hist(arima_coeffs[,1])
hist(p.adjust(
  2*pnorm(abs(arima_coeffs[,1]/sqrt(arima_vars[1,1,])), lower.tail = F),
  method = "BH"
))

hist(arima_coeffs[,2])
hist(p.adjust(
  2*pnorm(abs(arima_coeffs[,2]/sqrt(arima_vars[2,2,])), lower.tail = F),
  method = "BH"
))


non_fine_index <- testData2[cond == "nonfine", unique(group_num)]
nonfine_cor_arrays <- array(0, dim = c(3, 3, length(non_fine_index)))
yes_fine_index <- testData2[cond == "fine", unique(group_num)]
yesfine_cor_arrays <- array(0, dim = c(3, 3, length(yes_fine_index)))

sapply(1:length(non_fine_index), function(b){
  .GlobalEnv$nonfine_cor_arrays[,,b] <-
    split(testData2[group_num == non_fine_index[b], residuals],
          testData2[group_num == non_fine_index[b], subnum]) %>%
    unlist() %>% matrix(nc = 3) %>% cor()
  return(".")
})
sapply(1:length(yes_fine_index), function(b){
  .GlobalEnv$yesfine_cor_arrays[,,b] <-
    split(testData2[group_num == yes_fine_index[b], residuals],
          testData2[group_num == yes_fine_index[b], subnum]) %>%
    unlist() %>% matrix(nc = 3) %>% cor()
  return(".")
})

hist(c(yesfine_cor_arrays[1,2,],
       yesfine_cor_arrays[1,3,],
       yesfine_cor_arrays[2,3,]))
hist(c(nonfine_cor_arrays[1,2,],
       nonfine_cor_arrays[1,3,],
       nonfine_cor_arrays[2,3,]))

zvals_nonfine_cor <- nonfine_cor_arrays %>% (function(z) 0.5*log((1 + z)/(1 - z))*sqrt(97))
zvals_nonfine_cor <- c(zvals_nonfine_cor[1,2,],
                       zvals_nonfine_cor[1,3,],
                       zvals_nonfine_cor[2,3,])

zvals_yesfine_cor <- yesfine_cor_arrays %>% (function(z) 0.5*log((1 + z)/(1 - z))*sqrt(97))
zvals_yesfine_cor <- c(zvals_yesfine_cor[1,2,],
                       zvals_yesfine_cor[1,3,],
                       zvals_yesfine_cor[2,3,])

hist(p.adjust(
  2*pnorm(abs(zvals_nonfine_cor), lower.tail = F), 
  method = "BH"
))
hist(p.adjust(
  2*pnorm(abs(zvals_yesfine_cor), lower.tail = F),
  method = "BH"
))


data2 <- split(testData2[group_num == 2, residuals],
                 testData2[group_num == 2, subnum]) %>%
  unlist() %>% matrix(nc = 3)
colnames(data2) <- paste0("subject", 1:3)
