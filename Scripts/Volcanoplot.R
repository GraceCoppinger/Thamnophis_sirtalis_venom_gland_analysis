### DESEQ2 & Volcano Plots orginal script provided by Samuel Hirst with modifications by Grace Coppinger

# There will be a bit of repeat from the normalization

# Read the CSV file as a data frame without specifying row names
countData_df <- read.csv("gene_count_matrix_GC.csv", row.names = NULL)

# Extract unique row names from the "gene_id" column and make them unique
unique_row_names <- make.unique(countData_df$gene_id)

# Assign unique row names to the data frame
rownames(countData_df) <- unique_row_names
countData_df <- countData_df[-1]

# Convert the data frame to a matrix
countData <- as.matrix(countData_df)

colData <- read.csv("Tsirt_meta_data.csv", row.names = "Sample_ID")
colData$Sex <- factor(colData$Sex, levels = c("F", "M"))

all(rownames(colData) %in% colnames(countData))

countData <- countData[, rownames(colData)]

all(rownames(colData) == colnames(countData))

## This will be for labeling the volcano plots
# I have all my venom genes annotated so they begin with VENOM-. If you don't have that, you'll just have to search using a long character string using the toxin gene names
Venom <- c("VENOM-")

#Convert count data to integers if necessary
countData <- round(countData)

ddsPop <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ Sex)

ddsPop <- DESeq(ddsPop)

resPop <- results(ddsPop, contrast=c("Sex", "F", "M"))

all_genes <- data.frame(gene = rownames(resPop), pvalue = resPop$padj, log2FoldChange = resPop$log2FoldChange)

#This is will now give me a dataframe of the LFS and padj. 
write.csv(all_genes, "DESEQ2_pop.csv")

all_genes_plot <- all_genes
rownames(all_genes_plot) <- all_genes_plot$gene
all_genes_plot <- all_genes_plot[-1]

venom_genes <- all_genes[grepl(paste(Venom, collapse="|"), all_genes$gene), ]

significant_venom_genes <- venom_genes %>%
  filter(abs(log2FoldChange) >= 1, !is.na(pvalue), pvalue < 0.05)

Male_genes <- all_genes %>%
  filter(log2FoldChange <= -1, !is.na(pvalue), pvalue < 0.05)

Female_genes <- all_genes %>%
  filter(log2FoldChange >= 1, !is.na(pvalue), pvalue < 0.05)

keyvals <- ifelse(
  resPop$log2FoldChange < -1.0, '#2f2aa0',
  ifelse(resPop$log2FoldChange > 1.0, '#d30b94',
         'black'))
# Update 'keyvals' to set the color of selected genes to mAgenta
keyvals[rownames(resPop) %in% significant_venom_genes$gene] <- '#CCFF23'
keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == '#d30b94'] <- 'Female Biased'
names(keyvals)[keyvals == '#CCFF23'] <- 'Biased Venom Gene'
names(keyvals)[keyvals == '#2f2aa0'] <- 'Male Biased'

EnhancedVolcano(all_genes_plot,
                lab = rownames(all_genes_plot),
                x = 'log2FoldChange',
                y = 'pvalue',
                selectLab = c(significant_venom_genes$gene),
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 2.0,
                labSize = 3,
                labCol = 'black',
                labFace = 'bold',
                legendLabSize = 14,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                colConnectors = 'black', 
                max.overlaps = 100, 
                xlim = c(-14, 12),
                ylim = c(0, -log10(10e-12)),
                colCustom = keyvals)
