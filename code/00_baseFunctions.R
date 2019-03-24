ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}
ipak(c("parallel", "readr", "tidyverse", "data.table", "dtplyr", "pracma", "geepack"))
# ipak(c("corrplot", "msm", "markovchain", "lme4", "glmmTMB"))

createID <- function(x, factor = TRUE){
  x <- factor(x)
  x2 <- numeric(length(x))
  xlev <- levels(x)
  for(k in xlev){
    x2[x == k] <- which(xlev == k)
  }
  if(factor) return(factor(x2))
  return(x2)
}

makeAllNumeric <- function(df){
  for(i in 1:ncol(df)) df[,i] <- as.numeric(df[,i])
  return(df)
}

auto_cor <- Vectorize(function(x, lag) cor(x[1:(length(x) - lag)], x[(1 + lag):length(x)]), "lag")

nlargest <- function(x, k){
  n <- length(x)
  if(k > n) return(max(x))
  return(sort(x)[k])
}