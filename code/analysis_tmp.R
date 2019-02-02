source("code/00_baseFunctions.R")

testData <- fread("data/decision_making_experiment.csv")
str(testData)

testData[, `:=` (subnum = createID(subnum, F), # Create ID numbers for subjects
                 group_num = createID(group_num, F), # Create ID numbers for subjects who played together
                 cond = factor(cond, levels = 1:3,
                               labels = c("individual", "nonfine", "fine")), # Factorize group type
                 firstchoice = factor(firstchoice), # Factor subject choices
                 kind = factor(kind), # Factor kind
                 individualcoice = as.integer(individualcoice - 1), # Factorize subject final decision as 0 / 1
                 finaldecision = factor(finaldecision), # Factor group choices
                 # logreaction = log(reaction), # Define log(time) for further anlysis if needed
                 finaldec = as.integer(finaldec - 1))] # Factorize group final decision as 0 / 1

testData[round <= 60, `:=` (iscorrectInd = 1*(individualcoice == 1), 
                            iscorrectGr = 1*(finaldec == 1))] # t<=60, the right choise is 1
testData[round > 60, `:=` (iscorrectInd = 1*(individualcoice == 0),
                           iscorrectGr = 1*(finaldec == 0))] # t>60, the right choise is 0
setorder(testData, cond, group_num, subnum, round)

# Who has missing observations?
missingSub <- testData[, .N, by = subnum][N < 100]
missingSub
testData[subnum %in% missingSub[,subnum],
          .(minRound = min(round),
            maxRound = max(round),
            sumMissingInMiddle = sum(diff(round) != 1)),
          by = subnum]
# We can see that missing values are only at the end of the game. So this won't interfere with our Markov model.

str(testData)

######

testData2 <- copy(testData)
testData2[, `:=` (afterShock = 1*(round > 60), # Remove unneeded columns
                  firstchoice = NULL,
                  kind = NULL,
                  individualcoice = NULL,
                  finaldec = NULL,
                  finaldecision = NULL,
                  payoff = NULL,
                  foregone = NULL,
                  reaction = NULL)]




genZcor(clusz = table(testData2$subnum), waves = testData2$round, corstrv = 3)
genZcor(clusz = table(testData2$group_num), waves = testData2$subnum, corstrv = 2)

mod1 <- geeglm(iscorrectInd ~ round*afterShock*cond, family = binomial(link = "logit"), data = testData2,
              id = testData2$subnum, waves = round, corstr = "ar1")
# mod2 <- geese(iscorrectInd ~ round*afterShock*cond, family = binomial(link = "logit"), data = testData2,
#               id = testData2$subnum, waves = round, corstr = "ar1")

summary(mod1)
summary(mod2)
