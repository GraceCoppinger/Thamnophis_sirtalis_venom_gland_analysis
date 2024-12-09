# Check out ToxcodAn as the majority of this script is based on this material
#Rhett Rautsaw, Andrew Mason, and Samual Hirst also provided originally scripts for which mine are modified from.

#Load libraries
library(tidyverse)
library(ggplot2)
library(pheatmap)

#Load in the csv
TPM_df <- read_csv("Path_to_toxin_nontoxin_expression_tpm.csv")
TPM_df <- TPM_df %>% mutate_if(is.character,as.factor)

#source the below script from the toxcodan guide https://github.com/pedronachtigall/ToxCodAn/tree/master/Guide
source("PlottingFunctions.R")
TPM_df <- df_clean(TPM_df, class="class", toxin_family="toxin_family", colors=toxin_colors)

# Impute 0 values
# Take numeric values (columns 4-12), transpose them, impute zeros, and re-transpose
TPM_df2 <- t(cmultRepl(t(TPM_df[,4:25]),output = "p-counts"))

# Bind back to original dataframe
TPM_df2 <- cbind(TPM_df[,1:3], TPM_df2)
rownames(TPM_df2) <- TPM_df2$gene_id

#Plot the average transcriptome
FancyFigure(TPM_df2,"Average",class="class",toxin_family="toxin_family", colors=toxin_colors)

#There are also options to plot individuals or male/female average
FancyFigure(TPM_df2,"JLS178",class="class",toxin_family="toxin_family", colors=toxin_colors)
FancyFigure(TPM_df2,"JLS181",class="class",toxin_family="toxin_family", colors=toxin_colors)
FancyFigure(TPM_df2,"JLS187",class="class",toxin_family="toxin_family", colors=toxin_colors)
FancyFigure(TPM_df2,"JLS202",class="class",toxin_family="toxin_family", colors=toxin_colors)
FancyFigure(TPM_df2,"JLS216",class="class",toxin_family="toxin_family", colors=toxin_colors)
FancyFigure(TPM_df2,"JLS223",class="class",toxin_family="toxin_family", colors=toxin_colors)
FancyFigure(TPM_df2,"JLS224",class="class",toxin_family="toxin_family", colors=toxin_colors)
FancyFigure(TPM_df2,"JLS228",class="class",toxin_family="toxin_family", colors=toxin_colors)
FancyFigure(TPM_df2,"JLS242",class="class",toxin_family="toxin_family", colors=toxin_colors)
FancyFigure(TPM_df2,"JLS247",class="class",toxin_family="toxin_family", colors=toxin_colors)
FancyFigure(TPM_df2,"JLS248",class="class",toxin_family="toxin_family", colors=toxin_colors)
FancyFigure(TPM_df2,"JLS249",class="class",toxin_family="toxin_family", colors=toxin_colors)
FancyFigure(TPM_df2,"JLS253",class="class",toxin_family="toxin_family", colors=toxin_colors)
FancyFigure(TPM_df2,"JLS254",class="class",toxin_family="toxin_family", colors=toxin_colors)
FancyFigure(TPM_df2,"JLS255",class="class",toxin_family="toxin_family", colors=toxin_colors)
FancyFigure(TPM_df2,"JLS256",class="class",toxin_family="toxin_family", colors=toxin_colors)
FancyFigure(TPM_df2,"JLS257",class="class",toxin_family="toxin_family", colors=toxin_colors)
FancyFigure(TPM_df2,"JLS260",class="class",toxin_family="toxin_family", colors=toxin_colors)
FancyFigure(TPM_df2,"JLS269",class="class",toxin_family="toxin_family", colors=toxin_colors)
FancyFigure(TPM_df2,"Female_Average",class="class",toxin_family="toxin_family", colors=toxin_colors)
FancyFigure(TPM_df2,"Male_Average",class="class",toxin_family="toxin_family", colors=toxin_colors)
