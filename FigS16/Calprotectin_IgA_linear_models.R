data <- read.csv("Calprotectin_IgA.csv")
library(lme4)
head(data)

summary(lm(log2(Calprotectin_Day3) ~ Treatment, data = data))
summary(lm(log2(Calprotectin_Day7.8) ~ Treatment, data = data))
summary(lm(log2(Calprotectin_Day12) ~ Treatment, data = data))

summary(lm(log2(IgA_Day3) ~ Treatment, data = data))
summary(lm(log2(IgA_Day7.8) ~ Treatment, data = data))
summary(lm(log2(IgA_Day12) ~ Treatment, data = data))
