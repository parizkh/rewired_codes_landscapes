library(ggplot2)

data = read.csv("input/NBS_parS_fitness_ALL.csv")

WT_input = data$X.reads.pre.selection.1[data$variants=="RTAG"]
WT_selected = data$X.reads.post.selection..parS.[data$variants=="RTAG"]

# NA = 0 reads in post selection
data$X.reads.post.selection..parS.[is.na(data$X.reads.post.selection..parS.)] = 0
# in pre-selection, copy NA values from NBS
for (i in 1:nrow(data)) {
  if (is.na(data$X.reads.pre.selection.1[i])) data$X.reads.pre.selection.1[i] = data$X.reads.pre.selection[i]
}

# add the missing variants (0 reads post-selection; unknown number of reads pre-selection, will use 10)
# function to generate all sequences over given alphabet
gen_all_seqs = function(L, alphabet) {
  if (L==1) {
    return(alphabet)
  } else {
    prefixes = gen_all_seqs(L-1, alphabet)
    res = c()
    for (p in prefixes) {
      for (c in alphabet) {
        res = c(res, paste0(p,c))
      }
    }
    return(res)
  }
}
seqs_all = gen_all_seqs(4, c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"))
tmp_data = data.frame(variants = seqs_all)
med_count = median(data$X.reads.pre.selection.1[data$X.reads.post.selection..parS.==0])
data = merge(data, tmp_data, all.y=TRUE)
for (i in 1:nrow(data)) {
  if (is.na(data$X.reads.pre.selection.1[i])) {
    data$X.reads.pre.selection.1[i] = med_count
    data$X.reads.post.selection..parS.[i] = 0
  }
}

# compute the fitness and variance
data$LogFitness = log((data$X.reads.post.selection..parS.+0.5)/(WT_selected+0.5)) - 
  log((data$X.reads.pre.selection.1+0.5)/(WT_input+0.5))
data$SE = sqrt(1/(data$X.reads.pre.selection.1+0.5) + 1/(WT_input+0.5) + 
                 1/(data$X.reads.post.selection..parS.+0.5) + 1/(WT_selected+0.5))
data$Var = data$SE^2

data_toPrint = data[,c("variants", "LogFitness", "Var")]
write.table(data_toPrint, file="output/parS_data_processed.csv", quote=FALSE, sep=",", row.names=FALSE)

ggplot(data, aes(x=LogFitness)) +
  geom_histogram(alpha=0.5, color="black") + theme_test() +
  labs(x="log binding affinity")

ggplot(data, aes(x=SE)) +
  geom_histogram(alpha=0.5, color="black") + theme_test() +
  labs(x="standard error")

ggplot(data, aes(x=LogFitness, y=SE)) +
  geom_point(alpha=0.2) + theme_test() +
  labs(x="log binding affinity", y="standard error")




