source("code/02_get_data.R")

testData2
individual_mod <- geeglm(iscorrectInd ~ log(round)*afterShock, family = binomial(link = "logit"),
                         data = testData2, subset = (cond == "individual"),
                         id = subnum, corstr = "ar1")
summary(individual_mod)


dt <- copy(testData2[cond == "fine"])
nonfine_mod <- geeglm(iscorrectInd ~ log(round)*afterShock, family = binomial(link = "logit"),
                      data = dt, id = subnum, corstr = "ar1")

dt[, `:=` (resid = residuals(nonfine_mod, type = "pearson"),
           index = 1:.N)]
ar_coefs <- sapply(dt[,unique(subnum)],
                   function(i)
                     arima(dt[subnum == i, resid], order = c(1,0,0), include.mean = TRUE)$coef["ar1"])
ar_coefs <- cbind(dt[,.N, by = subnum][,N], ar_coefs)

indexs <- dt[, unique(group_num)]
cor_arrays <- 
  t(sapply(1:length(indexs), function(b){
    corrmat <-
      split(dt[group_num == indexs[b], resid],
            dt[group_num == indexs[b], subnum]) %>%
      unlist() %>% matrix(nc = 3) %>% cor()
    return(corrmat[lower.tri(corrmat)])
}, simplify = "array"))

pvals_cor <- 2*pnorm(abs(cor_arrays*sqrt(100 - 3)), lower.tail = FALSE)
pvals_cor <- t(apply(pvals_cor, 1, p.adjust, method = "BH"))
cor_arrays[pvals_cor > 0.1] <- NA

corrmeans <- rowMeans(cor_arrays, na.rm = T)
corrmeans[is.na(corrmeans)] <- 0

build_ar1_matrix <- function(length, phi){
  mat <- matrix(0, nr = length, nc = length)
  return(phi^abs(row(mat) - col(mat)))
}

ar_matrices <- 
  lapply(1:nrow(ar_coefs), function(i) build_ar1_matrix(ar_coefs[i,1], ar_coefs[i,2]))
for(i in 1:length(ar_matrices)) ar_matrices[[i]][ar_matrices[[i]] < 0.001] < 0

corr_matrix <- matrix(0, nr = sum(ar_coefs[,1]), nc = sum(ar_coefs[,1]))

index_start <- 1
for(i in 1:nrow(ar_coefs)){
  index_end <- sum(ar_coefs[1:i,1])
  corr_matrix[index_start:index_end, index_start:index_end] <- ar_matrices[[i]]
  index_start <- index_end + 1
}

unique_groupnum <- dt[,unique(group_num)]
for(k in 1:length(unique_groupnum)) for(i in 1:100){
  indexs <- dt[(group_num == unique_groupnum[k]) & (round == i), index]
  corr_matrix[indexs[1], indexs[2]] <- corrmeans[k]
  corr_matrix[indexs[2], indexs[1]] <- corrmeans[k]
  corr_matrix[indexs[1], indexs[3]] <- corrmeans[k]
  corr_matrix[indexs[3], indexs[1]] <- corrmeans[k]
  corr_matrix[indexs[2], indexs[3]] <- corrmeans[k]
  corr_matrix[indexs[3], indexs[2]] <- corrmeans[k]
}

ipak(corrplot)
corrplot::corrplot(corr_matrix[50:150, 50:150])
dt

X3 <- testData2[,.( group_num, subnum, 
                    indBeforeInter = (cond == "individual")*(round <= 60),
                    fineBeforeInter = (cond == "fine")*(round <= 60),
                    nfineBeforeInter = (cond == "nonfine")*(round <= 60),
                    indAfterInter = (cond == "individual")*(round > 60),
                    fineAfterInter = (cond == "fine")*(round > 60),
                    nfineAfterInter = (cond == "nonfine")*(round > 60),
                    indBeforeSlope = (cond == "individual")*log(round)*(round <= 60),
                    fineBeforeSlope = (cond == "fine")*log(round)*(round <= 60),
                    nfineBeforeSlope = (cond == "nonfine")*log(round)*(round <= 60),
                    indAfterSlope = (cond == "individual")*log(round)*(round > 60),
                    fineAfterSlope = (cond == "fine")*log(round)*(round > 60),
                    nfineAfterSlope = (cond == "nonfine")*log(round)*(round > 60),
                    iscorrectInd)]

gee.subnums <- geeglm(iscorrectInd ~ 0 + indBeforeInter + fineBeforeInter + nfineBeforeInter +
                        indAfterInter + fineAfterInter + nfineAfterInter +
                        indBeforeSlope + fineBeforeSlope + nfineBeforeSlope +
                        indAfterSlope + fineAfterSlope + nfineAfterSlope,
                      family = binomial(link = "logit"),
                      data = X3, id = subnum, corstr = "ar1")

vcovv <- summary(gee.subnums)$cov.scaled
coeffs <- gee.subnums$coefficients
C1 <- rbind(c(-1, 0, 0, 1, 0, 0, -log(60), 0, 0, log(61), 0, 0),
            c(0, -1, 0, 0, 1, 0, 0, -log(60), 0, 0, log(61), 0),
            c(0, 0, -1, 0, 0, 1, 0, 0, -log(60), 0, 0, log(61)))
rownames(C1) <- c("ind", "fine", "nfine")
gamma <- C1 %*% coeffs
vargamma <- C1 %*% vcovv %*% t(C1)

C2 <- rbind(c(1, -1, 0),
            c(1, 0, -1),
            c(0, 1, -1))
rownames(C2) <- c("ind-fine", "ind-nfine", "fine-nfine")

diff <- C2 %*% gamma
vardiff <- C2 %*% vargamma %*% t(C2)
Zval <- diff/sqrt(diag(vardiff))
Pval <- 2*pnorm(abs(Zval), lower.tail = F)

results3 <- list(data = X1, model = mod.gee.ur, vcov = vcovv, coeffs = coeffs, C1 = C1, gamma = gamma, vargamma = vargamma,
                 diff = diff, vardiff = vardiff, Zval = Zval, Pval = Pval)

p.adjust(results3$Pval, "BH")
