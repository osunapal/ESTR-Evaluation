library(dplyr)

#the general difference matrix one for all data (sons and daughters)
differences <- (child-((mother+father)/2))

#create separate differences dataframes for male and female
differences_daughters <- data.frame(matrix(nrow = nrow(differences), ncol = 0))
differences_sons <- differences_daughters
for(i in 1:nrow(IDs))
{
  fam <- IDs[i,"FamilyID"]
  if(IDs[i,"Sex"] == 1)
  {
    differences_sons[[fam]] <- differences[[fam]]
  }
  if(IDs[i,"Sex"] == 2)
  {
    differences_daughters[[fam]] <- differences[[fam]]
  }
}

df_pvals_gender <- data.frame(
  index = 1:nrow(differences),
  pvals_sons = NA,
  pvals_daughters = NA
)

total_iterations <- nrow(differences)
pb <- txtProgressBar(min = 0, max = total_iterations, style = 3)

#for only sons
no_eval_sons <- 0
for( i in 1:nrow(differences_sons) ) #look at every row in the matrix (at every locus)
{
  if(sum(is.na(differences_sons[i,])) <= length(differences_sons)-cutOff) # don't consider loci where >180 families are NA (so at least 100 values)
  {
    df_pvals_gender[i,"pvals_sons"]<-wilcox.test(as.numeric(differences_sons[i,]))$p.value
    no_eval_sons <- no_eval_sons+1
  }
  setTxtProgressBar(pb, i)
}
#close(pb)

#for only daughters
no_eval_daughters <- 0
for( i in 1:nrow(differences_daughters) ) #look at every row in the matrix (at every locus)
{
  if(sum(is.na(differences_daughters[i,])) <= length(differences_daughters)-cutOff) # don't consider loci where >180 families are NA
  {
    df_pvals_gender[i,"pvals_daughters"]<-wilcox.test(as.numeric(differences_daughters[i,]))$p.value
    no_eval_daughters <- no_eval_daughters+1
  }
  setTxtProgressBar(pb, i)
}
close(pb)

df_pvals_gender$row_means_daughters <- rowMeans(differences_daughters, na.rm = TRUE)
df_pvals_gender$row_means_sons <- rowMeans(differences_sons, na.rm = TRUE)
df_pvals_gender$log_pval_daughters <- -log10(df_pvals_gender$pvals_daughters)
df_pvals_gender$log_pval_sons <- -log10(df_pvals_gender$pvals_sons)

sum(as.numeric(df_pvals_gender$pvals_sons) < (0.05/no_eval_sons), na.rm = TRUE)
sum(as.numeric(df_pvals_gender$pvals_daughters) < (0.05/no_eval_daughters), na.rm = TRUE)

#list instable loci for sons
sons_instable_loci <- df_pvals_gender %>%
  select(index, pvals_sons) %>%
  filter(pvals_sons < 0.05/no_eval_sons)
sons_instable_loci <- sons_instable_loci %>%
  mutate(Locus = locus[index])

#list instable loci for daughters
daughters_instable_loci <- df_pvals_gender %>%
  select(index, pvals_daughters) %>%
  filter(pvals_daughters < 0.05/no_eval_daughters)
daugthers_instable_loci <- daughters_instable_loci %>%
  mutate(Locus = locus[index])



#Transmissions from fathers to children:
differences_from_fathers <- (child-father)

#Transmissions from mothers to children:
differences_from_mothers <- (child-mother)

total_iterations <- nrow(differences)
pb <- txtProgressBar(min = 0, max = total_iterations, style = 3)

df_pvals_fa_mo <- data.frame(
  index = 1:nrow(differences),
  pvals_fa = NA,
  pvals_mo = NA,
  pvals_mo_da = NA,
  pvals_fa_so = NA,
  pvals_mo_so = NA,
  pvals_fa_da = NA
)

cutOff <- 100 #to evaluate all with at least 100 data points

number_eval_loci_fa <- 0
#for the fathers differences matrix
for( i in 1:nrow(differences_from_fathers) ) #look at every row in the matrix (at every locus)
{
  if(sum(is.na(differences_from_fathers[i,])) <= length(differences_from_fathers)-cutOff) # don't consider loci where >500 families are NA
  {
    df_pvals_fa_mo[i,"pvals_fa"]<-wilcox.test(as.numeric(differences_from_fathers[i,]))$p.value
    number_eval_loci_fa <- number_eval_loci_fa+1
  }
  setTxtProgressBar(pb, i)
}

#for the mothers differences matrix
number_eval_loci_mo <- 0
for( i in 1:nrow(differences_from_mothers) ) #look at every row in the matrix (at every locus)
{
  if(sum(is.na(differences_from_mothers[i,])) <= length(differences_from_mothers)-cutOff) # don't consider loci where >500 families are NA
  {
    df_pvals_fa_mo[i,"pvals_mo"]<-wilcox.test(as.numeric(differences_from_mothers[i,]))$p.value
    number_eval_loci_mo <- number_eval_loci_mo+1
  }
  setTxtProgressBar(pb, i)
}
#close(pb)


#create separate differences dataframes for male and female
differences_mother_daughter <- data.frame(matrix(nrow = nrow(differences), ncol = 0))
differences_father_son <- differences_mother_daughter
differences_father_daughter <- data.frame(matrix(nrow = nrow(differences), ncol = 0))
differences_mother_son <- data.frame(matrix(nrow = nrow(differences), ncol = 0))
for(i in 1:nrow(IDs))
{
  fam <- IDs[i,"FamilyID"]
  if(IDs[i,"Sex"] == 1)
  {
    differences_mother_son[[fam]] <- differences_from_mothers[[fam]]
    differences_father_son[[fam]] <- differences_from_fathers[[fam]]
  }
  if(IDs[i,"Sex"] == 2)
  {
    differences_father_daughter[[fam]] <- differences_from_fathers[[fam]]
    differences_mother_daughter[[fam]] <- differences_from_mothers[[fam]]
  }
}



#for father -> sons
no_eval_father_son <- 0
for( i in 1:nrow(differences_father_son) ) #look at every row in the matrix (at every locus)
{
  if(sum(is.na(differences_father_son[i,])) <= length(differences_father_son)-cutOff) # don't consider loci where >180 families are NA (so at least 100 values)
  {
    df_pvals_fa_mo[i,"pvals_fa_so"]<-wilcox.test(as.numeric(differences_father_son[i,]))$p.value
    no_eval_father_son <- no_eval_father_son+1
  }
  setTxtProgressBar(pb, i)
}
#close(pb)

#for mother -> daughters
no_eval_mother_daughter <- 0
for( i in 1:nrow(differences_mother_daughter) ) #look at every row in the matrix (at every locus)
{
  if(sum(is.na(differences_mother_daughter[i,])) <= length(differences_mother_daughter)-cutOff) # don't consider loci where >180 families are NA
  {
    df_pvals_fa_mo[i,"pvals_mo_da"]<-wilcox.test(as.numeric(differences_mother_daughter[i,]))$p.value
    no_eval_mother_daughter <- no_eval_mother_daughter+1
  }
  setTxtProgressBar(pb, i)
}

#for mothers -> sons
no_eval_mother_son <- 0
for( i in 1:nrow(differences_mother_son) ) #look at every row in the matrix (at every locus)
{
  if(sum(is.na(differences_mother_son[i,])) <= length(differences_mother_son)-cutOff) # don't consider loci where >180 families are NA
  {
    df_pvals_fa_mo[i,"pvals_mo_so"]<-wilcox.test(as.numeric(differences_mother_son[i,]))$p.value
    no_eval_mother_son <- no_eval_mother_son+1
  }
  setTxtProgressBar(pb, i)
}

#for fathers -> daughters
no_eval_father_daughter <- 0
for( i in 1:nrow(differences_father_daughter) ) #look at every row in the matrix (at every locus)
{
  if(sum(is.na(differences_father_daughter[i,])) <= length(differences_father_daughter)-cutOff) # don't consider loci where >180 families are NA
  {
    df_pvals_fa_mo[i,"pvals_fa_da"]<-wilcox.test(as.numeric(differences_father_daughter[i,]))$p.value
    no_eval_father_daughter <- no_eval_father_daughter+1
  }
  setTxtProgressBar(pb, i)
}


#How many Loci have p<0.05? Don't count NA
sum(as.numeric(df_pvals_fa_mo$pvals_fa) < (0.05/number_eval_loci_fa), na.rm = TRUE) #Bonferroni correction for fathers
sum(as.numeric(df_pvals_fa_mo$pvals_mo) < (0.05/number_eval_loci_mo), na.rm = TRUE) #Bonferroni correction for mothers
sum(as.numeric(df_pvals_fa_mo$pvals_fa_so) < (0.05/no_eval_father_son), na.rm = TRUE) #Bonferroni correction for mothers
sum(as.numeric(df_pvals_fa_mo$pvals_mo_da) < (0.05/no_eval_mother_daughter), na.rm = TRUE) #Bonferroni correction for mothers
sum(as.numeric(df_pvals_fa_mo$pvals_mo_so) < (0.05/no_eval_mother_son), na.rm = TRUE) #Bonferroni correction for mothers
sum(as.numeric(df_pvals_fa_mo$pvals_fa_da) < (0.05/no_eval_father_daughter), na.rm = TRUE) #Bonferroni correction for mothers


#list instable loci from fathers
from_fathers_instable_loci <- df_pvals_fa_mo %>%
  select(index, pvals_fa) %>%
  filter(pvals_fa < 0.05/number_eval_loci_fa)
from_fathers_instable_loci <- from_fathers_instable_loci %>%
  mutate(Locus = locus[index])

#list instable loci from mothers
from_mothers_instable_loci <- df_pvals_fa_mo %>%
  select(index, pvals_mo) %>%
  filter(pvals_mo < 0.05/number_eval_loci_mo)
from_mothers_instable_loci <- from_mothers_instable_loci %>%
  mutate(Locus = locus[index])

