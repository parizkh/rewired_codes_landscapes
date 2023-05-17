library(ggplot2)

data = read.csv("input/220518allcodoncounts_wt.csv")
colnames(data)[1] = "Variant"

# the wildtype
WT = "DKE"

WT_input_rep1 = data$ParE3wt_GRE_rep1_t0[data$Variant==WT]
WT_input_rep2 = data$ParE3wt_GRE_rep2_t0[data$Variant==WT]
WT_selected_rep1 = data$ParE3wt_GRE_rep1_t600[data$Variant==WT]
WT_selected_rep2 = data$ParE3wt_GRE_rep2_t600[data$Variant==WT]

# estimate fitness replicate 1
data$LogFitness_rep1 = log((data$ParE3wt_GRE_rep1_t600+0.5)/(WT_selected_rep1+0.5)) - 
  log((data$ParE3wt_GRE_rep1_t0+0.5)/(WT_input_rep1[1]+0.5))
# estimate fitness replicate 2
data$LogFitness_rep2 = log((data$ParE3wt_GRE_rep2_t600+0.5)/(WT_selected_rep2+0.5)) - 
  log((data$ParE3wt_GRE_rep2_t0+0.5)/(WT_input_rep2[1]+0.5))
# variance replicate 1
data$Var_rep1 = 1/(data$ParE3wt_GRE_rep1_t0+0.5) + 1/(WT_input_rep1+0.5) + 
  1/(data$ParE3wt_GRE_rep1_t600+0.5) + 1/(WT_selected_rep1+0.5)
# variance replicate 2
data$Var_rep2 = 1/(data$ParE3wt_GRE_rep2_t0+0.5) + 1/(WT_input_rep2+0.5) + 
  1/(data$ParE3wt_GRE_rep2_t600+0.5) + 1/(WT_selected_rep2+0.5)

# combined fitness
data$LogFitness_mean = (data$LogFitness_rep1/data$Var_rep1 + data$LogFitness_rep2/data$Var_rep2)/(1/data$Var_rep1+1/data$Var_rep2)
# combined variance
data$Var = 1/(1/data$Var_rep1+1/data$Var_rep2)

# discard sequences containing stop codons
data = data[!grepl("*", data$Variant, fixed=TRUE),]

# save
data_toWrite = data[,c("Variant", "LogFitness_mean", "Var")]
write.table(data_toWrite, file="output/ParE3_data_processed.csv", quote=FALSE, sep=",", row.names=FALSE, col.names=FALSE)


# plots
ggplot(data, aes(x=LogFitness_mean)) +
  geom_histogram(alpha=0.5, color="black") + theme_test() +
  labs(x="Fitness")
#ggsave("hist_fitness.pdf")

ggplot(data, aes(x=Var)) +
  geom_histogram(alpha=0.5, color="black") + theme_test() +
  labs(x="Variance")
#ggsave("hist_var.pdf")

ggplot(data, aes(x=LogFitness_mean, y=Var)) +
  geom_point(alpha=0.2) + theme_test() +
  labs(x="Fitness", y="Variance")
#ggsave("corr_fitness_var.pdf")


