source("code/02_get_data.R")

compute_ci <- function(x, type, sig.level = 0.05, m_adjust = 1){
  mean(x) + switch(type, "l" = -1, "u" = 1)*qnorm(1 - sig.level/(2*m_adjust))*sd(x)/sqrt(length(x))
}
plot_autocor <- function(lag){
  return(
    autocor_data[,auto_cor(iscorrectInd, lag), by = .(cond, subnum)] %>% ggplot(aes(V1)) +
      geom_histogram(aes(y=..density..), col = "white", fill = "lightblue", bins = 10) +
      facet_grid(.~cond) + 
      labs(x = paste0("Auto-Correlation, lag=", lag), y = "Density", title = "Individual Auto-Correlation Within Conditions")
  )
}

autocor_data <- copy(testData2)
autocor_data[,indPrevCorr := c(FALSE, iscorrectInd[1:(nrow(testData)-1)])]
autocor_data[round == "1", indPrevCorr := iscorrectInd]
autocor_data[,grPrevCorr := c(FALSE, iscorrectGr[1:(nrow(testData)-1)])]
autocor_data[round == "1", grPrevCorr := iscorrectInd]
autocors <- 1:20

autocor_dt <- testData2[,.(auto_cor = auto_cor(iscorrectInd, autocors)), by = .(cond, subnum)]
autocor_dt[,lag := rep(autocors, 208)]

auto_corr_by_lag <- 
  autocor_dt[,mean(auto_cor), keyby = .(cond, lag)] %>%
  ggplot(aes(x = lag, y = V1, col = cond)) + geom_smooth(se = FALSE) +
  labs(x = "Lag", y = "Average Auto-Correlation", col = "Condition",
       title = "Average Subject Auto-Correlation by Condition") + 
  scale_colour_discrete(labels = c("Individuals", "Non-Penalty Groups", "Penalty Groups")) +
  theme(legend.position="bottom") + scale_x_continuous(breaks = (0:10)*2) + geom_hline(yintercept = 0)

auto_cor_plots <- lapply(1:3, function(k) plot_autocor(k))

prob_of_higher_gain <- 
  data.frame(Round = 1:100, Right = c(rep(0.9, 60), rep(0.5, 40)), Left = 0.7) %>%
  gather(key = Button, value = Probability, -Round) %>% 
  ggplot(aes(x = Round, y = Probability, col = Button)) +
  geom_line(size = 1) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
  theme(legend.position = "bottom") + scale_y_continuous(breaks = 0:10/10, labels = 0:10/10) +
  ggtitle("Probability of Higher Gain")

expected_rewards <- 
  data.frame(Round = 1:100, Right = c(rep(5.8, 60), rep(1, 40)), Left = 3.4) %>%
  gather(key = Button, value = `Expected Reward`, -Round) %>% 
  ggplot(aes(x = Round, y = `Expected Reward`, col = Button)) +
  geom_line(size = 1) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
  theme(legend.position = "bottom") + scale_y_continuous(breaks = 0:10, labels = 0:10) +
  ggtitle("Probability of Higher Gain")


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
  
exploratory_plots <- list(auto_corr_by_lag = auto_corr_by_lag,
                          auto_cor_plots = auto_cor_plots,
                          prob_of_higher_gain = prob_of_higher_gain,
                          expected_rewards = expected_rewards,
                          cum_mean_plt = cum_mean_plt,
                          roll_mean_plt = roll_mean_plt)

rm(compute_ci, plot_autocor, autocor_data, autocors, autocor_dt,
   auto_corr_by_lag, auto_cor_plots, prob_of_higher_gain, expected_rewards,
   group_data, cum_mean_plt, roll_mean_plt)
