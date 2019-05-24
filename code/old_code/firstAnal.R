source("code/01_auxilary_functions.R")

testData <- fread("data/decision_making_experiment.csv")
str(testData)

testData[,subnum := createID(subnum)]
testData[,group_num := createID(group_num)]
testData[,cond := factor(cond)]
testData[,firstchoice := factor(firstchoice)]
testData[,kind := factor(kind)]
testData[,finaldecision := factor(finaldecision)]

testData[firstchoice != finaldecision, .N]
testData[round <= 60, iscorrInd := (individualcoice == 2)]
testData[round > 60, iscorrInd := individualcoice != 2]
testData[round <= 60, iscorrGr := (finaldec == 2)]
testData[round > 60, iscorrGr := finaldec != 2]
testData[,individualcoice := as.integer(individualcoice - 1)]
testData[,finaldec := as.integer(finaldec - 1)]

testData[,indPrevCorr := c(FALSE, iscorrInd[1:(nrow(testData)-1)])]
testData[round == "1", indPrevCorr := iscorrInd]
testData[,grPrevCorr := c(FALSE, iscorrGr[1:(nrow(testData)-1)])]
testData[round == "1", grPrevCorr := iscorrInd]

str(testData)

testData[subnum == "17", 1*iscorrGr] %>% pacf

testData[,mean(iscorrGr), keyby = cond]
testData[,mean(iscorrInd), keyby = cond]

testData[,mean(iscorrGr), keyby = .(cond, round > 60)]
testData[,mean(iscorrInd), keyby = .(cond, round > 60)]

cor(testData[,.(grPrevCorr,iscorrGr)])
cor(testData[,.(indPrevCorr,iscorrInd)])


testData[,auto_cor(iscorrGr, 1), by = subnum][,2] %>% unlist %>% hist
testData[,auto_cor(iscorrGr, 2), by = subnum][,2] %>% unlist %>% hist
testData[,auto_cor(iscorrGr, 3), by = subnum][,2] %>% unlist %>% hist
testData[,auto_cor(iscorrGr, 7), by = subnum][,2] %>% unlist %>% hist
testData[,auto_cor(iscorrGr, 12), by = subnum][,2] %>% unlist %>% hist
summary(testData[,cor(grPrevCorr,iscorrGr), by = subnum])

testData[cond == 1,auto_cor(iscorrGr, 1), by = subnum][,2] %>% unlist %>% hist
testData[cond == 2,auto_cor(iscorrGr, 1), by = subnum][,2] %>% unlist %>% hist
testData[cond == 3,auto_cor(iscorrGr, 1), by = subnum][,2] %>% unlist %>% hist

plot_autocor <- function(lag){
  return(
    testData[,auto_cor(iscorrInd, lag), by = .(cond, subnum)] %>% ggplot(aes(V1)) +
      geom_histogram(aes(y=..density..), col = "white", fill = "lightblue", bins = 10) +
      facet_grid(.~cond, labeller = as_labeller(c("1" = "Individuals", "2" = "Non-Penalty Groups", "3" = "Penalty Groups"))) + 
      labs(x = paste0("Auto-Correlation, lag=", lag), y = "Density", title = "Individual Auto-Correlation Within Conditions")
  )
}

autocors <- 1:20

tmp <- testData[,auto_cor(iscorrInd, autocors), by = .(cond, subnum)]
tmp[,lag := rep(autocors, 208)]
tmp[,mean(V1), keyby = .(cond, lag)] %>% ggplot(aes(x = lag, y = V1, col = cond)) + geom_smooth(se = FALSE) + 
  labs(x = "Lag", y = "Average Auto-Correlation", col = "Condition", title = "Average Subject Auto-Correlation by Condition") + 
  scale_colour_discrete(labels = c("Individuals", "Non-Penalty Groups", "Penalty Groups")) +
  theme(legend.position="bottom") + scale_x_continuous(breaks = (0:10)*2) + geom_hline(yintercept = 0)
plot_autocor(5)

glm(iscorrInd ~ (round + (round > 60))*cond, testData, family = binomial(link = "logit")) %>% summary()
# mod <- glmer(iscorrInd ~ (round + (round > 60))*cond + (1|group_num), testData, family = binomial(link = "logit"))
# summary(mod)

X <- data.frame(Round = 1:100, Right = c(rep(0.9, 60), rep(0.5, 40)), Left = 0.7) %>% gather(key = Button, value = Probability, -Round)
ggplot(X, aes(x = Round, y = Probability, col = Button)) + geom_line(size = 1) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
  theme(legend.position = "bottom") + scale_y_continuous(breaks = 0:10/10, labels = 0:10/10) + ggtitle("Probability of Higher Gain")
