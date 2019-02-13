source("code/00_baseFunctions.R")
source("code/group_anlysis.R")

testData4 <- unique(testData3[,.(round, group_num, cond, iscorrectGr)])
X1 <- testData4[,.( group_num,
                    indBeforeInter = (cond == "individual")*(round <= 60),
                    fineBeforeInter = (cond == "fine")*(round <= 60),
                    nfineBeforeInter = (cond == "nonfine")*(round <= 60),
                    indAfterInter = (cond == "individual")*(round > 60),
                    fineAfterInter = (cond == "fine")*(round > 60),
                    nfineAfterInter = (cond == "nonfine")*(round > 60),
                    indBeforeSlope = (cond == "individual")*round*(round <= 60),
                    fineBeforeSlope = (cond == "fine")*round*(round <= 60),
                    nfineBeforeSlope = (cond == "nonfine")*round*(round <= 60),
                    indAfterSlope = (cond == "individual")*round*(round > 60),
                    fineAfterSlope = (cond == "fine")*round*(round > 60),
                    nfineAfterSlope = (cond == "nonfine")*round*(round > 60),
                    iscorrectGr)]

mod1 <- geeglm(iscorrectGr ~ 0 + indBeforeInter + fineBeforeInter + nfineBeforeInter +
                 indAfterInter + fineAfterInter + nfineAfterInter +
                 indBeforeSlope + fineBeforeSlope + nfineBeforeSlope +
                 indAfterSlope + fineAfterSlope + nfineAfterSlope,
               family = binomial(link = "logit"),
               data = X1, id = group_num, corstr = "ar1")

####

X2 <- testData4[,.( group_num,
                    fineBeforeInter = (cond == "fine")*(round <= 60),
                    nfineBeforeInter = (cond != "fine")*(round <= 60),
                    fineAfterInter = (cond == "fine")*(round > 60),
                    nfineAfterInter = (cond != "fine")*(round > 60),
                    fineBeforeSlope = (cond == "fine")*round*(round <= 60),
                    nfineBeforeSlope = (cond != "fine")*round*(round <= 60),
                    fineAfterSlope = (cond == "fine")*round*(round > 60),
                    nfineAfterSlope = (cond != "fine")*round*(round > 60),
                    iscorrectGr)]

mod2 <- geeglm(iscorrectGr ~ 0 + fineBeforeInter + nfineBeforeInter + fineAfterInter + nfineAfterInter +
                 fineBeforeSlope + nfineBeforeSlope + fineAfterSlope + nfineAfterSlope,
               family = binomial(link = "logit"),
               data = X2, id = group_num, corstr = "ar1")

anova_gee <- anova(mod2, mod1)
anova_gee
# The data suggests that there is no big defference between individuals and non-fine groups.

vcovv <- summary(mod1)$cov.scaled
coeffs <- mod1$coefficients
C1 <- rbind(c(-1, 0, 0, 1, 0, 0, -60, 0, 0, 61, 0, 0),
            c(0, -1, 0, 0, 1, 0, 0, -60, 0, 0, 61, 0),
            c(0, 0, -1, 0, 0, 1, 0, 0, -60, 0, 0, 61))
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

results1 <- list(data = X1, model = mod1, vcov = vcovv, coeffs = coeffs, C1 = C1, gamma = gamma, vargamma = vargamma,
                 diff = diff, vardiff = vardiff, Zval = Zval, Pval = Pval)

###

vcovv <- summary(mod2)$cov.scaled
coeffs <- mod2$coefficients
C1 <- rbind(c(-1, 0, 1, 0, -60, 0, 61, 0),
            c(0, -1, 0, 1, 0, -60, 0, 61))
rownames(C1) <- c("fine", "nfine")
gamma <- C1 %*% coeffs
vargamma <- C1 %*% vcovv %*% t(C1)

C2 <- rbind(c(1, -1))
rownames(C2) <- c("diff")

diff <- C2 %*% gamma
vardiff <- C2 %*% vargamma %*% t(C2)
Zval <- diff/sqrt(diag(vardiff))
Pval <- 2*pnorm(abs(Zval), lower.tail = F)

results2 <- list(data = X2, model = mod2, vcov = vcovv, coeffs = coeffs, C1 = C1, gamma = gamma, vargamma = vargamma,
                 diff = diff, vardiff = vardiff, Zval = Zval, Pval = Pval)

results1$gamma
results1$diff
results1$Pval
results2$gamma
results2$diff
results2$Pval