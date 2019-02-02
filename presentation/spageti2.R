source("presentation/firstAnal.R")

byTime <- testData[,.(round, group_num, cond, iscorrGr)] %>% unique() %>% spread(key = round, value = iscorrGr)
cummeans1 <- byTime[,3:62] %>% apply(1, cummean) %>% t
cummeans2 <- byTime[,63:102] %>% apply(1, cummean) %>% t
cummeans <- cbind(cummeans1, cummeans2)

averagecummeans <- cbind(unique(byTime[,cond]),
                         sapply(1:3, function(i) cummeans[byTime[,cond] == i,] %>%
                                  colMeans(na.rm = T)) %>% t) %>% as.data.frame()
colnames(averagecummeans) <- c("cond", gsub("V", "", colnames(byTime)[-(1:2)]))

averagecummeans <- gather(averagecummeans, key = round, value = cummean, -cond)
averagecummeans$round <- as.numeric(averagecummeans$round)
averagecummeans$group_num <- 0
averagecummeans$ismean <- TRUE
averagecummeans <- arrange(averagecummeans, cond, round)


byTime <- cbind(byTime[,1:2], cummeans)
colnames(byTime)[-(1:2)] <- gsub("V", "", colnames(byTime)[-(1:2)])

byTime <- gather(byTime, key = round, value = cummean, -cond, -group_num)
byTime$round <- as.numeric(byTime$round)
byTime <- arrange(byTime, cond, group_num, round)
byTime$ismean <- FALSE
averagecummeans <- select(averagecummeans, colnames(byTime))

ggplot() + geom_line(data = byTime, mapping = aes(x = round, y = cummean, col = factor(cond), group = group_num)) +
  geom_line(data = averagecummeans, mapping = aes(x = round, y = cummean, col = factor(cond)), size = 1) + 
  geom_vline(xintercept = 0) + geom_hline(yintercept = 0) + 
  geom_vline(xintercept = 60.5, size = 3) + geom_hline(yintercept = 0.5, col = "purple") +
  labs(title = "Cumulative Average by Group", x = "Round", y = "Cumulative Mean", col = "Group Type")
