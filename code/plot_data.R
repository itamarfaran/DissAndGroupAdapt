source("code/analysis_tmp.R")

compute_ci <- function(x, type, sig.level = 0.05, m_adjust = 1){
  mean(x) + switch(type, "l" = -1, "u" = 1)*qnorm(1 - sig.level/(2*m_adjust))*sd(x)/sqrt(length(x))
}

group_data <- testData2[,.(iscorrectGr = mean(iscorrectInd)), by = .(round, group_num, cond)]

cum_mean_plt <- 
  rbind(
    group_data[round <= 60, .(cond, round, cummean_correct = cummean(iscorrectGr)), by = .(group_num)],
    group_data[round > 60, .(cond, round, cummean_correct = cummean(iscorrectGr)), by = .(group_num)]
    ) %>%
  ggplot(aes(x = round, y = movavg(cummean_correct, 5, "e"), col = cond, by = factor(group_num))) +
  geom_line(alpha = 0.2) +
  stat_summary(aes(group = cond, col = cond), fun.y = mean, geom = "line", size = 1.5) +
  stat_summary(aes(group = cond, col = cond), fun.y = compute_ci, geom = "line", size = 1.2, linetype = 2,
               fun.args = list(type = "l", m_adjust = 3)) +
  stat_summary(aes(group = cond, col = cond), fun.y = compute_ci, geom = "line", size = 1.2, linetype = 2,
               fun.args = list(type = "u", m_adjust = 3)) +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
  labs(title = "Cummulative Mean of Group Responses",
       x = "Round", y = "Probability of Success", col = "Condition")

roll_mean_plt <-
  testData2[,.(iscorrectGr = mean(iscorrectInd)), by = .(round, group_num, cond)] %>%
  ggplot(aes(x = round, y = movavg(iscorrectGr, 10, "e"), col = cond, by = factor(group_num))) +
  geom_line(alpha = 0.15) +
  stat_summary(aes(group = cond, col = cond), fun.y = mean, geom = "line", size = 1.5) +
  stat_summary(aes(group = cond, col = cond), fun.y = compute_ci, geom = "line", size = 1.2, linetype = 2,
               fun.args = list(type = "l", m_adjust = 3)) +
  stat_summary(aes(group = cond, col = cond), fun.y = compute_ci, geom = "line", size = 1.2, linetype = 2,
               fun.args = list(type = "u", m_adjust = 3)) +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
  labs(title = "Rolling Mean of Group Responses",
       x = "Round", y = "Probability of Success", col = "Condition")
  
