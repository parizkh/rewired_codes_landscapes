library(ggplot2)

data = read.delim("input/elife-16965-supp1-v4.csv")

WT_input = data$Count.input[1]
WT_selected = data$Count.selected[1]

data$LogFitness = log((data$Count.selected+0.5)/(WT_selected+0.5)) - 
  log((data$Count.input+0.5)/(WT_input[1]+0.5))
data$SE = sqrt(1/(data$Count.input+0.5) + 1/(WT_input+0.5) + 
                 1/(data$Count.selected+0.5) + 1/(WT_selected+0.5))
data$Var = data$SE^2

write.table(data, file="output/GB1_data_processed.csv", quote=FALSE, sep=",", row.names=FALSE)


# plots
ggplot(data, aes(x=LogFitness)) +
  geom_histogram(alpha=0.5, color="black") + theme_test() +
  labs(x="log binding affinity")

ggplot(data, aes(x=SE)) +
  geom_histogram(alpha=0.5, color="black") + theme_test() +
  labs(x="standard error")

ggplot(data, aes(x=LogFitness, y=SE)) +
  geom_point(alpha=0.2) + theme_test() +
  labs(x="log binding affinity", y="standard error")




