library(dplyr)
#the general difference matrix one for all data (sons and daughters)
differences <- (child-((mother+father)/2))

total_iterations <- nrow(differences)
pb <- txtProgressBar(min = 0, max = total_iterations, style = 3)

df_pvals <- data.frame(
  index = 1:nrow(differences),
  pvals_all = NA
)

number_eval_loci <- 0
cutOff <- 100 #to evaluate all with at least one data point
#for the whole differences matrix
for( i in 1:nrow(differences) ) #look at every row in the matrix (at every locus)
{
  if(sum(is.na(differences[i,])) <= length(differences)-cutOff) # don't consider loci where >500 families are NA
  {
    df_pvals[i,"pvals_all"]<-wilcox.test(as.numeric(differences[i,]))$p.value
    number_eval_loci <- number_eval_loci+1
  }
  setTxtProgressBar(pb, i)
}
close(pb)

df_pvals$row_means_all <- rowMeans(differences, na.rm = TRUE)
df_pvals$log_pval_all <- -log10(df_pvals$pvals_all)

p_cutoff <- (0.05/number_eval_loci) #Bonferroni correction
#How many Loci have p<0.05? Don't count NA
number_instable_loci_whole_data <- sum(as.numeric(df_pvals$pvals_all) < p_cutoff, na.rm = TRUE)

#list instable loci
instable_loci_whole_data <- df_pvals %>%
  select(index, pvals_all, row_means_all) %>%
  filter(pvals_all < p_cutoff)
instable_loci_whole_data <- instable_loci_whole_data %>%
  mutate(Locus = locus[index])

cat("Examining the whole dataset, there are ", number_instable_loci_whole_data, " instable loci out of ", number_eval_loci, " evaluated loci.")