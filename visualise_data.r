#Below code - add sequence name to 1st and mapped_reads to 2nd column
{r}
parent_directory <- "/vol/data/smallRNAseq_12102023/X204SC23083269-Z01-F001/MasterThesis/results/idxstats"

# Get a list of all files in the parent directory and its subfolders
files <- list.files(parent_directory, pattern = "\\.log$", full.names = TRUE, recursive = TRUE)

# Loop through each file
for (file in files) {
  # Read the data from the file into a data frame
  data <- read.table(file, header = FALSE)
  
  # Extract the first and third columns
  table_data <- data[, c(1, 3)]
  
  # Add column names
  colnames(table_data) <- c("sequence_name", "mapped_reads")
  
  # Generate the output file path
  output_file <- sub("_hairpin_idxstats.log$", "_mapped_reads.txt", file)
  
  # Save the resulting table with headers
  write.table(table_data, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE)
  
  # Print a message indicating completion
  cat("Processed file:", file, "\n")
}

#Below code - combine all idxstats log file
{r}
library(readr)
library(dplyr)
library(tidyr)

# Define base directory and filenames
base_dir <- "/vol/data/smallRNAseq_12102023/X204SC23083269-Z01-F001/MasterThesis/results/idxstats/"
file_names <- c("Col0SR_1.sorted.bam.log", "Col0SR_2.sorted_bam.log", "Col0SR_3.sorted.bam.log", "u1c_1.sorted.bam.log", "u1c_2.sorted.bam.log", "u1c_3.sorted.bam.log", "u170k_1.sorted.bam.log", "u170k_2.sorted.bam.log", "u170k_3.sorted.bam.log")
full_paths <- paste0(base_dir, file_names)

# Read and process each file
data_frames <- lapply(full_paths, function(file) {
  df <- read_tsv(file, col_names = c("sequence", "mapped_reads"), show_col_types = FALSE)
  print(paste("Reading:", file, "with", nrow(df), "rows"))  # Debug: print rows read

  # Ensure mapped_reads is numeric
  df <- df %>%
    mutate(mapped_reads = as.numeric(as.character(mapped_reads)), 
           file_id = sub("\\.sorted.*", "", basename(file)))

  # Check for and report any conversion problems
  if (anyNA(df$mapped_reads)) {
    warning("NA introduced by coercion in file: ", file)
  }
  
  return(df)
})

# Combine and reshape data
combined_data <- bind_rows(data_frames) %>%
  pivot_wider(names_from = file_id, values_from = mapped_reads, names_prefix = "", values_fill = list(mapped_reads = 0))

# Save and report
output_path <- paste0(base_dir, "combined_log.tsv")
write_tsv(combined_data, output_path)
print(paste("Combined file created at", output_path, "with", nrow(combined_data), "total rows"))

#Differential analysis of gene for u1c
{r}

setwd("/vol/data/smallRNAseq_12102023/X204SC23083269-Z01-F001/MasterThesis/results/idxstats")

library(edgeR)

x <- read.delim("combined_log.tsv", header = TRUE)
colnames(x) <- c("sequence_name", "Col0SR_1", "Col0SR_2", "Col0SR_3", "u1c_1", "u1c_2", "u1c_3", "u170k_1", "u170k_2", "u170k_3")

x <- x[-1, ]
rownames(x) <- x[, 1]

x[, -1] <- sapply(x[, -1], as.numeric)


# Define groups for the samples
group <- factor(c(1, 1, 1, 2, 2, 2, 3, 3, 3))

#rownames(x) <- x[, 1]

x <- x[!rowSums(is.na(x)), ]

y <- DGEList(counts=x[, -1], group=group)

y <- DGEList(counts=x, group=group)

keep <- filterByExpr(y)
y <- y[keep, , keep.lib.sizes=FALSE]

y <- calcNormFactors(y)

design <- model.matrix(~group)
y <- estimateDisp(y, design)

fit <- glmQLFit(y, design)

qlf.2vs1 <- glmQLFTest(fit, coef = 2)

topTags(qlf.2vs1)

#qlf.3vs1 <- glmQLFTest(fit, coef=3)
#topTags(qlf.3vs1)

dfu1c <- topTags(qlf.2vs1, n = Inf)

write.csv(dfu1c, "dfu1c_significant.csv", row.names = T)

tmp <- read_csv("dfu1c_significant.csv")

de_u1c <- tmp[complete.cases(tmp), ]

write.csv(de_u1c, "volcano_data_u1c.csv", row.names = T)

de_u1c$diffexpressed <- "NO"
 de_u1c$diffexpressed [de_u1c$logFC > 1 & de_u1c$FDR < 0.05] <- "UP"
 de_u1c$diffexpressed [de_u1c$logFC < -1 & de_u1c$FDR < 0.05] <- "DOWN"
 
 p <- ggplot(data = de_u1c, aes(x = logFC, y = -log10(FDR), color = diffexpressed)) +
    geom_point() +
    theme_minimal() +
    labs(title = "Differential Expression in U1c") +
    theme(plot.title = element_text(hjust = 0.5)) +
    geom_vline(xintercept=c(-1, 1), color="red") +
    geom_hline(yintercept=-log10(0.05), color="red") +
    scale_color_manual(values=c("blue", "black", "red"))

print(p)

#differential gene expression analysis for u170k
{r}
setwd("/vol/data/smallRNAseq_12102023/X204SC23083269-Z01-F001/MasterThesis/results/idxstats")

library(edgeR)

x <- read.delim("combined_log.tsv", header = TRUE)
colnames(x) <- c("sequence_name", "Col0SR_1", "Col0SR_2", "Col0SR_3", "u1c_1", "u1c_2", "u1c_3", "u170k_1", "u170k_2", "u170k_3")

x <- x[-1, ]
rownames(x) <- x[, 1]

x[, -1] <- sapply(x[, -1], as.numeric)


# Define groups for the samples
group <- factor(c(1, 1, 1, 2, 2, 2, 3, 3, 3))

#rownames(x) <- x[, 1]

x <- x[!rowSums(is.na(x)), ]

y <- DGEList(counts=x[, -1], group=group)

y <- DGEList(counts=x, group=group)

keep <- filterByExpr(y)
y <- y[keep, , keep.lib.sizes=FALSE]

y <- calcNormFactors(y)

design <- model.matrix(~group)
y <- estimateDisp(y, design)

fit <- glmQLFit(y, design)


qlf.3vs1 <- glmQLFTest(fit, coef=3)
topTags(qlf.3vs1)

dfu170k <- topTags(qlf.3vs1, n = Inf)

write.csv(dfu170k, "dfu170k_significant.csv", row.names = T)

tmp1 <- read_csv("dfu170k_significant.csv")

de_u170k <- tmp1[complete.cases(tmp1), ]

de_u170k$diffexpressed <- "NO"
 de_u170k$diffexpressed [de_u170k$logFC > 1 & de_u170k$FDR < 0.05] <- "UP"
 de_u170k$diffexpressed [de_u170k$logFC < -1 & de_u170k$FDR < 0.05] <- "DOWN"
 
 write.csv(de_u170k, "volcano_data_u170k.csv", row.names = TRUE)

 
 p1 <- ggplot(data = de_u170k, aes(x = logFC, y = -log10(FDR), color = diffexpressed)) +
    geom_point() +
    theme_minimal() +
    labs(title = "Differential Expression in U170k") +
    theme(plot.title = element_text(hjust = 0.5)) +
    geom_vline(xintercept=c(-1, 1), color="red") +
    geom_hline(yintercept=-log10(0.05), color="red") +
    scale_color_manual(values=c("blue", "black", "red"))

 
print(p1)


#Code - PC plot

{r}

setwd("/vol/data/smallRNAseq_12102023/X204SC23083269-Z01-F001/MasterThesis/DATA/TEs_results/idxstats")

library(edgeR)
library(ggplot2)
library(ggrepel)

x <- read.delim("combined_log.tsv", header = TRUE)
group <- factor(c(1, 1, 1, 2, 2, 2, 3, 3, 3))
design <- model.matrix(~group)
x <- x[!rowSums(is.na(x)), ]

y <- DGEList(counts = x, group = group)
keep <- filterByExpr(y)
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
y <- estimateDisp(y, design)
tmm_log <- edgeR :: cpm(y, log = T, prior.count = 1)
tmm <- tmm_log

pca <- prcomp(t(tmm_log), center=TRUE, scale=TRUE)
pca[["rotation"]]
summary = summary(pca)
summary

pca = prcomp(t(tmm_log), scale=T)
scores = data.frame(pca$x)
rownames(scores)
loadings = pca$rotation
scores = data.frame(pca$x)
legend =gsub("_[1,2,3]","",row.names(scores))
scores$label = row.names(scores)
scores
scores$genotype
scores$genotype = substr(scores$label, start = 1, stop = nchar(scores$label)-2)
pca1 <- ggplot(data=scores,aes(x = PC1, y = PC2, color = legend)) +
geom_point(size=3)+
  theme_bw() +
  geom_text_repel(aes(label=label),size=2,color="black")+
  labs(x = paste0("PC1 (",round(summary$importance[2,1]*100),"%)"),
       y = paste0("PC2 (",round(summary$importance[2,2]*100),"%)"),
       title = "PC plot TE") 
