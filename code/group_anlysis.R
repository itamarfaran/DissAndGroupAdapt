source("code/00_baseFunctions.R")
source("code/analysis_tmp.R")

testData3 <- unique(testData2[,.(round, group_num, cond, iscorrectGr)])
X1 <- testData3[,.( group_num,
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

mod.gee.ur <- geeglm(iscorrectGr ~ 0 + indBeforeInter + fineBeforeInter + nfineBeforeInter +
                     indAfterInter + fineAfterInter + nfineAfterInter +
                     indBeforeSlope + fineBeforeSlope + nfineBeforeSlope +
                     indAfterSlope + fineAfterSlope + nfineAfterSlope,
                   family = binomial(link = "logit"),
                   data = X1, id = group_num, corstr = "ar1")

mod.glm.ur <- glm(iscorrectGr ~ 0 + indBeforeInter + fineBeforeInter + nfineBeforeInter +
                  indAfterInter + fineAfterInter + nfineAfterInter +
                  indBeforeSlope + fineBeforeSlope + nfineBeforeSlope +
                  indAfterSlope + fineAfterSlope + nfineAfterSlope,
                family = binomial(link = "logit"), data = X1)

testData3[,`:=` (pred.gee.ur = predict(mod.gee.ur, type = "response"),
                 res.gee.ur = resid(mod.gee.ur, type = "response"),
                 pred.glm.ur = predict(mod.glm.ur, type = "response"),
                 res.glm.ur = resid(mod.glm.ur, type = "response")
                  )]

####

X2 <- testData3[,.( group_num,
                    fineBeforeInter = (cond == "fine")*(round <= 60),
                    nfineBeforeInter = (cond != "fine")*(round <= 60),
                    fineAfterInter = (cond == "fine")*(round > 60),
                    nfineAfterInter = (cond != "fine")*(round > 60),
                    fineBeforeSlope = (cond == "fine")*round*(round <= 60),
                    nfineBeforeSlope = (cond != "fine")*round*(round <= 60),
                    fineAfterSlope = (cond == "fine")*round*(round > 60),
                    nfineAfterSlope = (cond != "fine")*round*(round > 60),
                    iscorrectGr)]

mod.gee.r <- geeglm(iscorrectGr ~ 0 + fineBeforeInter + nfineBeforeInter + fineAfterInter + nfineAfterInter +
                       fineBeforeSlope + nfineBeforeSlope + fineAfterSlope + nfineAfterSlope,
                     family = binomial(link = "logit"),
                     data = X2, id = group_num, corstr = "ar1")

mod.glm.r <- glm(iscorrectGr ~ 0 + fineBeforeInter + nfineBeforeInter + fineAfterInter + nfineAfterInter +
                    fineBeforeSlope + nfineBeforeSlope + fineAfterSlope + nfineAfterSlope,
                  family = binomial(link = "logit"), data = X2)

testData3[,`:=` (pred.gee.r = predict(mod.gee.r, type = "response"),
                 res.gee.r = resid(mod.gee.r, type = "response"),
                 pred.glm.r = predict(mod.glm.r, type = "response"),
                 res.glm.r = resid(mod.glm.r, type = "response")
)]


anova_gee <- anova(mod.gee.r, mod.gee.ur)
anova_gee

anova_glm <- anova(mod.glm.r, mod.glm.ur)
anova_glm

# The data suggests that there is no big defference between individuals and non-fine groups.

vcovv <- summary(mod.gee.ur)$cov.scaled
coeffs <- mod.gee.ur$coefficients
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

results1 <- list(data = X1, model = mod.gee.ur, vcov = vcovv, coeffs = coeffs, C1 = C1, gamma = gamma, vargamma = vargamma,
                 diff = diff, vardiff = vardiff, Zval = Zval, Pval = Pval)

###

vcovv <- summary(mod.gee.r)$cov.scaled
coeffs <- mod.gee.r$coefficients
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

results2 <- list(data = X2, model = mod.gee.r, vcov = vcovv, coeffs = coeffs, C1 = C1, gamma = gamma, vargamma = vargamma,
                 diff = diff, vardiff = vardiff, Zval = Zval, Pval = Pval)

results1$gamma
results1$diff
results1$Pval
results2$gamma
results2$diff
results2$Pval

autocor <- copy(testData3[,.(corr = auto_cor(res.gee.ur, 1)), by = .(group_num, cond, round > 60)])
autocor[, fisher_z := 0.5*log((1 + corr) / (1 - corr))]
autocor
autocor[,.(mean = mean(fisher_z), sd = sd(fisher_z)), by = .(cond, round)]
anova(
  lmer(fisher_z ~ round + (1|group_num), data = autocor),
  lmer(fisher_z ~ cond + round + (1|group_num), data = autocor)
)

summary(aov(fisher_z ~ cond*round, data = autocor))
summary(aov(fisher_z ~ cond + round, data = autocor))
oneway.test(fisher_z ~ cond, data = autocor)
kruskal.test(fisher_z ~ cond, data = autocor)


testData3[,auto_cor(pred.gee.ur, 2), by = group_num][,V1] %>% hist()