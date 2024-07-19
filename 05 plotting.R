#plotting the overall instabilities
library(ggplot2)
ggplot(df_pvals, aes(x = index, y = log_pval_all)) +
  geom_point() +
  geom_hline(yintercept = -log10(0.05/553), color = "red", linetype = "dashed") +
  annotate("text", x = 100, y = -log10(0.05/553), label = "Significance\nthreshold", color = "red", vjust = -1) +
  labs(
    #title = "Log p values of evaluated loci for whole dataset",
    x = "Index of loci",
    y = "-log10(p-value)"
  ) +
  theme_minimal()

#plotting depending on gender
ggplot(df_pvals_gender) +
  geom_point(aes(x = index, y = log_pval_daughters, color = "Daughters")) +
  geom_point(aes(x = index, y = log_pval_sons, color = "Sons")) +
  geom_hline(yintercept = -log10(0.000137), color = "red", linetype = "dashed") +
  annotate("text", x = 200, y = -log10(0.000137), label = "Significance\nthreshold", color = "red", vjust = -1) +
  labs(x = "Index of loci",
       y = "-log10(P-value)") +
  scale_color_manual(name = "Gender",
                     values = c("Daughters" = "pink", "Sons" = "skyblue")) +
  theme_minimal()
