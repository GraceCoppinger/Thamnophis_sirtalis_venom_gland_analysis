#Run Transcriptome.R first.
#Obtain just the pie graphs to generate male and femaale comparison
ToxinPie(TPM_df2, "Female_Average", class="class")
ToxinPie(TPM_df2, "Male_Average", class="class")
