library(ggplot2)
library(tibble)
library(seqinr)
library(stringi)
library(dplyr)

data = readRDS("input/sequence_reads_counts.RDS")
rownames(data$before)[1:10]
df = t(data$before)
colnames(df) = paste0(colnames(df), "_pre")
df2 = t(data$`M9+TMP`)
colnames(df2) = paste0(colnames(df2), "_post")

df = merge(df, df2, by=0, all=TRUE)
df = rownames_to_column(df, "Variant_nuc")

# extract the relevant codons
df$Variant_nuc = substr(df$Row.names, 52, 60)
# translate
df$Variant_aa = sapply(df$Variant_nuc, function(x) stri_paste(translate(s2c(x)), collapse=""))

# sum reads for all nuc encodings of a given protein variant
df_aa = aggregate(.~Variant_aa, df[,3:15], sum)                      

# Compute fitness and variance for each replicate
WT="ADL"
for (rep in 1:6) {
  WT_pre = df_aa[df_aa==WT,rep+1]
  WT_post = df_aa[df_aa==WT,rep+7]
  # fitness
  df_aa[,paste0("Fitness_", rep)] = log((df_aa[,rep+7]+0.5)/(WT_post+0.5)) - 
    log((df_aa[,rep+1]+0.5)/(WT_pre+0.5))
  # var
  df_aa[,paste0("Var", rep)] = 1/(df_aa[,rep+7]+0.5) + 1/(WT_post+0.5) + 
    1/(df_aa[,rep+1]+0.5) + 1/(WT_pre+0.5)
}

# combine the replicates
fitness = 0
sum_inv_vars = 0
for (rep in 1:6) {
    fitness = fitness + df_aa[,paste0("Fitness_", rep)]/df_aa[,paste0("Var",rep)]
    sum_inv_vars = sum_inv_vars + 1/df_aa[,paste0("Var",rep)]
}
df_aa$Fitness = fitness/sum_inv_vars
df_aa$Var = 1/sum_inv_vars

# check correctness with the Papkou data
Papkou = readRDS("input/fitness_data_wt.rds")
Papkou$Variant_aa = paste0(Papkou$aa1,Papkou$aa2,Papkou$aa3)
tmp = merge(df_aa, Papkou, all=TRUE)
cor.test(tmp$Fitness, tmp$m)

# save the data
data_toWrite = df_aa[,c("Variant_aa", "Fitness", "Var")]
data_toWrite = data_toWrite[!grepl("*", data_toWrite$Variant_aa, fixed=TRUE),]
write.table(data_toWrite, file="output/DHFR_data_processed.csv", quote=FALSE, sep=",", row.names=FALSE, col.names=FALSE)


# what is the fitness cut-off, using our method of computing fitness?
for (rep in 1:6) {
  WT_pre = df_aa[df_aa==WT,rep+1]
  WT_post = df_aa[df_aa==WT,rep+7]
  # fitness
  df[,paste0("Fitness_", rep)] = log((df[,rep+8]+0.5)/(WT_post+0.5)) - 
    log((df[,rep+2]+0.5)/(WT_pre+0.5))
  # var
  df[,paste0("Var", rep)] = 1/(df[,rep+8]+0.5) + 1/(WT_post+0.5) + 
    1/(df[,rep+2]+0.5) + 1/(WT_pre+0.5)
}

# combine the replicates
fitness = 0
sum_inv_vars = 0
for (rep in 1:6) {
  fitness = fitness + df[,paste0("Fitness_", rep)]/df[,paste0("Var",rep)]
  sum_inv_vars = sum_inv_vars + 1/df[,paste0("Var",rep)]
}
df$Fitness = fitness/sum_inv_vars
df$Var = 1/sum_inv_vars

# find the 93% quantile
q=quantile(df$Fitness, 0.93, na.rm=TRUE)

# only run this after the smoothing
data_DHFR_smoothed = read.csv("../01_vcregression/output/map.txt")
sum(data_DHFR_smoothed<q)/8000
data_DHFR_smoothed[data_DHFR_smoothed<q] = -100
write.table(data_DHFR_smoothed, file="../01_vcregression/output/map_wNonfuctional.csv", sep=",", quote=FALSE, row.names=FALSE)
