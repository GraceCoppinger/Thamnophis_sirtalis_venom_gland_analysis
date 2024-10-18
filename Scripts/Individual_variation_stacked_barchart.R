
library(tidyr)
library(dplyr)
library(ggplot2)

# Read data
toxin_data <- read.csv("/Users/gracecoppinger/Desktop/Toxins_over500tpm.csv", stringsAsFactors = FALSE)

# Reshape data from wide to long format for each toxin family
toxin_data_long <- toxin_data %>%
  pivot_longer(cols = c(CTL, SVMP, FTx, CRISP, Other), names_to = "ToxinFamily", values_to = "Transcripts")

# Calculate Total Toxin Transcripts
toxin_data_long <- toxin_data_long %>%
  group_by(ID) %>%
  mutate(TotalTranscripts = sum(Transcripts)) %>%
  ungroup() %>%
  mutate(Percentage = Transcripts / TotalTranscripts * 100)

# Sort toxin_data_long by ID based on SVL from toxin_data
toxin_data_long <- toxin_data_long %>%
  arrange(SVL) %>%
  mutate(ID = factor(ID, levels = unique(ID[order(SVL)])))
# Specify the order of toxin families for stacking
toxin_order <- c("CTL", "SVMP", "FTx", "CRISP", "Other")

# Plot
ggplot(toxin_data_long, aes(x = factor(ID, levels = unique(ID)), y = Transcripts, fill = ToxinFamily)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#0ec434", "#2f2aa0", "#d30b94", "#772b9d", "#a4a4a4"),
                    breaks = toxin_order,
                    labels = c("CTL", "SVMP", "FTx", "CRISP", "Other")) +  # Adjusted labels
  labs(x = "Individual", y = "Percentage of Toxin Transcripts", title = "Stacked Bar Plot of Toxin Transcripts by Individual Sorted by SVL") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
