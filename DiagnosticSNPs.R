library(vcfR)
library(triangulaR)

setwd("D:/Analyses/Triturus/Diagnostic")

# Read your data
data <- read.vcfR("TivaTmac_allfilt.recode.vcf")
meta <- read.csv("metadata_TivaTmac.csv")

# Create population map
popmap <- data.frame(
  "id" = as.character(meta$Sample.ID[meta$use.in.data=="yes"]),
  "pop" = meta$Species[meta$use.in.data=="yes"]
)

# Get diagnostic SNPs (use your preferred threshold)
vcfR.diff <- alleleFreqDiff(
  vcfR = data, 
  pm = popmap, 
  p1 = "Tiva",
  p2 = "Tmac",
  difference = 1.0  # or 1.0 if you want fully diagnostic
)

n_diagnostic_snps <- nrow(vcfR.diff@fix)
cat("Number of diagnostic SNPs:", n_diagnostic_snps, "\n")

# Extract genotype matrix from the diagnostic SNPs
# vcfR.diff contains only the diagnostic SNPs
gt_matrix <- extract.gt(vcfR.diff, element = "GT")

# Check dimensions
cat("Genotype matrix dimensions:", nrow(gt_matrix), "SNPs x", ncol(gt_matrix), "individuals\n")

# Function to convert VCF genotypes to Genepop format
vcf_to_genepop <- function(gt) {
  # VCF format: "0/0", "0/1", "1/1", etc.
  # Genepop format: "001001", "001002", "002002", etc.
  
  if(is.na(gt) || gt == "./.") {
    return("000000")  # Missing data
  }
  
  # Split the genotype
  alleles <- unlist(strsplit(gt, "[/|]"))
  
  # Convert to numeric and add 1 (so 0 becomes 1, 1 becomes 2, etc.)
  allele1 <- as.numeric(alleles[1]) + 1
  allele2 <- as.numeric(alleles[2]) + 1
  
  # Format as 3-digit codes
  genepop_gt <- sprintf("%03d%03d", allele1, allele2)
  
  return(genepop_gt)
}

# Convert all genotypes to Genepop format
genepop_matrix <- matrix(NA, nrow = ncol(gt_matrix), ncol = nrow(gt_matrix))

for(i in 1:nrow(gt_matrix)) {  # For each SNP
  for(j in 1:ncol(gt_matrix)) {  # For each individual
    genepop_matrix[j, i] <- vcf_to_genepop(gt_matrix[i, j])
  }
}

# Add sample names
rownames(genepop_matrix) <- colnames(gt_matrix)

# Create locus names
locus_names <- paste0(vcfR.diff@fix[, "CHROM"], "_", vcfR.diff@fix[, "POS"])

# Get population assignments for each sample
sample_ids <- rownames(genepop_matrix)
sample_pops <- popmap$pop[match(sample_ids, meta$Sample.ID)]

# Create a dataframe with samples and their populations
sample_data <- data.frame(
  id = sample_ids,
  pop = sample_pops,
  stringsAsFactors = FALSE
)

# Add the genotype matrix
sample_data <- cbind(sample_data, genepop_matrix)

# Check population distribution
cat("\n=== POPULATION DISTRIBUTION ===\n")
print(table(sample_data$pop))

# Write Genepop file
output_file <- "your_output_genepop.txt"

# Open file connection
con <- file(output_file, "w")

# Write header
writeLines("Diagnostic SNPs - Genepop format", con)

# Write locus names
writeLines(locus_names, con)

# Get unique populations
unique_pops <- unique(sample_data$pop)
unique_pops <- unique_pops[!is.na(unique_pops)]  # Remove NAs if any

cat("\n=== WRITING POPULATIONS ===\n")

# Write each population separately
for(pop in unique_pops) {
  cat("Writing population:", pop, "\n")
  
  # Write "Pop" separator
  writeLines("Pop", con)
  
  # Get individuals from this population
  pop_samples <- sample_data[sample_data$pop == pop, ]
  
  # Write genotypes for each individual in this population
  for(i in 1:nrow(pop_samples)) {
    # Get individual ID
    ind_id <- pop_samples$id[i]
    
    # Get genotypes (skip the first two columns which are id and pop)
    genotypes <- paste(pop_samples[i, 3:ncol(pop_samples)], collapse = " ")
    
    # Write in Genepop format: "IndividualName , genotypes"
    writeLines(paste0(ind_id, " ,  ", genotypes), con)
  }
}

close(con)


cat("\n=== GENEPOP FILE CREATED ===\n")
cat("File:", output_file, "\n")
cat("Number of loci:", n_diagnostic_snps, "\n")
cat("Number of individuals:", nrow(genepop_matrix), "\n")
cat("\nYou can now use this file as input for recom-sim.py\n")


# Show first few lines of the file to verify
cat("\n=== FIRST FEW LINES OF OUTPUT ===\n")
cat(readLines(output_file, n = 10), sep = "\n")
