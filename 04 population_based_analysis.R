#create seperate differences dataframes for superpopulations
differences <- (child-((mother+father)/2))

differences_sp_AMR <- data.frame(matrix(nrow = nrow(differences), ncol = 0))
differences_sp_EAS <- differences_sp_AMR
differences_sp_EUR <- differences_sp_AMR
differences_sp_AFR <- differences_sp_AMR
differences_sp_SAS <- differences_sp_AMR

for(i in 1:nrow(IDs))
{
  fam <- IDs[i,"FamilyID"]
  if(IDs[i,"Superpopulation"] == "AMR")
  {
    differences_sp_AMR[[fam]] <- differences[[fam]]
  }
  if(IDs[i,"Superpopulation"] == "EAS")
  {
    differences_sp_EAS[[fam]] <- differences[[fam]]
  }
  if(IDs[i,"Superpopulation"] == "EUR")
  {
    differences_sp_EUR[[fam]] <- differences[[fam]]
  }
  if(IDs[i,"Superpopulation"] == "AFR")
  {
    differences_sp_AFR[[fam]] <- differences[[fam]]
  }
  if(IDs[i,"Superpopulation"] == "SAS")
  {
    differences_sp_SAS[[fam]] <- differences[[fam]]
  }
}

total_iterations <- nrow(differences)
pb <- txtProgressBar(min = 0, max = total_iterations, style = 3)

df_pvals_superpop <- data.frame(
  index = 1:nrow(differences),
  pvals_AFR = NA,
  pvals_AMR = NA,
  pvals_EAS = NA,
  pvals_EUR = NA,
  pvals_SAS = NA
)

cutOff <- 50
#AFR
no_eval_AFR <- 0
for( i in 1:nrow(differences_sp_AFR) ) #look at every row in the matrix (at every locus)
{
  if(sum(is.na(differences_sp_AFR[i,])) <= (ncol(differences_sp_AFR)-cutOff)) # don't consider loci where <50 families are not NA
  {
    df_pvals_superpop[i,"pvals_AFR"]<-wilcox.test(as.numeric(differences_sp_AFR[i,]))$p.value
    no_eval_AFR <- no_eval_AFR+1
  }
  setTxtProgressBar(pb, i)
}

#AMR
no_eval_AMR <- 0
for( i in 1:nrow(differences_sp_AMR) ) #look at every row in the matrix (at every locus)
{
  if(sum(is.na(differences_sp_AMR[i,])) <= (ncol(differences_sp_AMR)-cutOff)) # don't consider loci where <50 families are not NA
  {
    df_pvals_superpop[i,"pvals_AMR"]<-wilcox.test(as.numeric(differences_sp_AMR[i,]))$p.value
    no_eval_AMR <- no_eval_AMR+1
  }
  setTxtProgressBar(pb, i)
}

#EAS
no_eval_EAS <- 0
for( i in 1:nrow(differences_sp_EAS) ) #look at every row in the matrix (at every locus)
{
  if(sum(is.na(differences_sp_EAS[i,])) <= (ncol(differences_sp_EAS)-cutOff)) # don't consider loci where <50 families are not NA
  {
    df_pvals_superpop[i,"pvals_EAS"]<-wilcox.test(as.numeric(differences_sp_EAS[i,]))$p.value
    no_eval_EAS <- no_eval_EAS+1
  }
  setTxtProgressBar(pb, i)
}

#EUR
no_eval_EUR <- 0
for( i in 1:nrow(differences_sp_EUR) ) #look at every row in the matrix (at every locus)
{
  if(sum(is.na(differences_sp_EUR[i,])) <= (ncol(differences_sp_EUR)-cutOff)) # don't consider loci where <50 families are not NA
  {
    df_pvals_superpop[i,"pvals_EUR"]<-wilcox.test(as.numeric(differences_sp_EUR[i,]))$p.value
    no_eval_EUR <- no_eval_EUR+1
  }
  setTxtProgressBar(pb, i)
}

#SAS
no_eval_SAS <- 0
for( i in 1:nrow(differences_sp_SAS) ) #look at every row in the matrix (at every locus)
{
  if(sum(is.na(differences_sp_SAS[i,])) <= (ncol(differences_sp_SAS)-cutOff)) # don't consider loci where <50 families are not NA
  {
    df_pvals_superpop[i,"pvals_SAS"]<-wilcox.test(as.numeric(differences_sp_SAS[i,]))$p.value
    no_eval_SAS <- no_eval_SAS+1
  }
  setTxtProgressBar(pb, i)
}

#For superpopulations: How many are below p-value?
sum(as.numeric(df_pvals_superpop$pvals_AFR) < (0.05/no_eval_AFR), na.rm = TRUE)
sum(as.numeric(df_pvals_superpop$pvals_AMR) < (0.05/no_eval_AMR), na.rm = TRUE)
sum(as.numeric(df_pvals_superpop$pvals_EAS) < (0.05/no_eval_EAS), na.rm = TRUE)
sum(as.numeric(df_pvals_superpop$pvals_EUR) < (0.05/no_eval_EUR), na.rm = TRUE)
sum(as.numeric(df_pvals_superpop$pvals_SAS) < (0.05/no_eval_SAS), na.rm = TRUE)