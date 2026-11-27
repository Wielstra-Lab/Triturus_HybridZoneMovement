#https://github.com/zgompert/bgc-hm/tree/main

library(devtools)
library(bgchm)

#devtools::install_github("zgompert/bgc-hm")

############################
# 1. Load packages
############################

library(bgchm)
library(readxl)

setwd("D:/Analyses/Triturus/Bgc")

# ============================================================
# bgchm genomic cline analysis: Tiva × Tana
# ============================================================


# ------------------------------------------------------------
# 1. Read genotype data
# ------------------------------------------------------------

dat <- read.csv(
  "TivaTcrigtHIMiss.csv",
  header = TRUE,
  check.names = FALSE
)

locus_names <- as.character(dat[[1]])

G <- as.matrix(dat[, -1])
mode(G) <- "numeric"
rownames(G) <- locus_names

# Save your existing genome-wide hybrid index
existing_HI <- as.numeric(G["Average", ])
names(existing_HI) <- colnames(G)

# Remove Average row
G <- G[rownames(G) != "Average", ]

# ------------------------------------------------------------
# 2. Read parental reference IDs
# ------------------------------------------------------------

refs <- read_excel("REFs.xlsx")

Tiva_ids <- na.omit(as.character(refs$Tiva))
Tcri_ids <- na.omit(as.character(refs$Tcri))

# Check that all reference IDs exist
stopifnot(all(Tiva_ids %in% colnames(G)))
stopifnot(all(Tcri_ids %in% colnames(G)))

hybrid_ids <- setdiff(colnames(G), c(Tiva_ids, Tcri_ids))

cat("Loci:", nrow(G), "\n")
cat("Individuals:", ncol(G), "\n")
cat("Tiva references:", length(Tiva_ids), "\n")
cat("Tcri references:", length(Tcri_ids), "\n")
cat("Hybrids:", length(hybrid_ids), "\n")

# ------------------------------------------------------------
# 3. Arrange data as individuals × loci
# ------------------------------------------------------------

G_Tiva <- t(G[, Tiva_ids, drop = FALSE])
G_Tcri <- t(G[, Tcri_ids, drop = FALSE])
G_Hyb  <- t(G[, hybrid_ids, drop = FALSE])

colnames(G_Tiva) <- rownames(G)
colnames(G_Tcri) <- rownames(G)
colnames(G_Hyb)  <- rownames(G)

# ------------------------------------------------------------
# 4. Estimate parental allele frequencies
# ------------------------------------------------------------

cat("\nEstimating parental allele frequencies...\n")

p_out <- est_p(
  G0 = G_Tiva * 2,
  G1 = G_Tcri * 2,
  model = "genotype",
  ploidy = "diploid",
  HMC = FALSE
)

p0 <- p_out$p0[, 1]
p1 <- p_out$p1[, 1]

# ------------------------------------------------------------
# 5. Convert hybrid genotypes
#    0 / 0.5 / 1  ->  0 / 1 / 2
# ------------------------------------------------------------

G_Hyb <- G_Hyb * 2

# Ploidy information:
# 2 = observed diploid genotype
# 0 = missing genotype
P_Hyb <- ifelse(is.na(G_Hyb), 0, 2)

# ------------------------------------------------------------
# 6. Estimate hybrid indices
# ------------------------------------------------------------

cat("\nEstimating hybrid indexes...\n")

h_out <- est_hi(
  Gx = G_Hyb,
  p0 = p0,
  p1 = p1,
  model = "genotype",
  ploidy = "mixed",
  pldat = P_Hyb
)

HI_bgchm <- h_out$hi[, 1]
names(HI_bgchm) <- hybrid_ids

# ------------------------------------------------------------
# 7. Compare with your existing Average HI
# ------------------------------------------------------------

HI_compare <- data.frame(
  Individual = hybrid_ids,
  HI_existing = existing_HI[hybrid_ids],
  HI_bgchm = HI_bgchm
)

HI_compare$Difference <-
  HI_compare$HI_bgchm - HI_compare$HI_existing

print(head(HI_compare))

cat(
  "\nCorrelation between existing and bgchm HI:",
  cor(
    HI_compare$HI_existing,
    HI_compare$HI_bgchm,
    use = "complete.obs"
  ),
  "\n"
)

# ------------------------------------------------------------
# 8. Estimate genomic clines
# ------------------------------------------------------------

cat("\nEstimating genomic clines...\n")

gc_out <- est_genocl(
  Gx = G_Hyb,
  p0 = p0,
  p1 = p1,
  H = HI_bgchm,
  model = "genotype",
  ploidy = "mixed",
  pldat = P_Hyb,
  hier = TRUE,
  sd0 = 2
)

# ------------------------------------------------------------
# 9. Save results
# ------------------------------------------------------------

write.csv(
  HI_compare,
  "bgchm_hybrid_indices_Tcri.csv",
  row.names = FALSE
)

save(
  h_out,
  gc_out,
  p_out,
  HI_compare,
  file = "bgchm_results_Tcri.RData"
)

cat("\nAnalysis complete.\n")

sz95 <- sum2zero(
  hmc = gc_out$gencline_hmc,
  transform = TRUE,
  ci = 0.95
)

alpha <- sz95$center

# Check the column names first
colnames(alpha)

# Flag alpha intervals that exclude 0.5
alpha_outliers <- alpha[, "2.5%"] > 0.5 |
  alpha[, "97.5%"] < 0.5

table(alpha_outliers)

# Direction of alpha deviations
alpha_tiva <- alpha[, "2.5%"] > 0.5
alpha_tcri <- alpha[, "97.5%"] < 0.5

c(
  total_loci = nrow(alpha),
  Tiva_direction = sum(alpha_tiva),
  Tcri_direction = sum(alpha_tcri),
  alpha_outliers = sum(alpha_outliers),
  not_flagged = sum(!alpha_outliers)
)

cvec <- ifelse(alpha_outliers, "red", "grey75")

cvec <- grDevices::hcl.colors(
  n = nrow(gc_out$center),
  palette = "Dynamic"
)

# Original plot
gencline_plot(
   center = gc_out$center, 
   v = gc_out$gradient, 
   pdf = TRUE, 
   outf  = "cline_TivaTcri_plot.pdf", 
   cvec = cvec, 
   xvec = NULL)

# Sum2zero adjusted plot
gencline_plot(
  center = sz95$center,
  v = sz95$gradient,
  pdf = TRUE,
  outf = "cline_TivaTcri_sum2zero_plot.pdf",
  cvec = cvec,
  xvec = NULL
)

# Remove the Average row, if present
locus_names <- locus_names[locus_names != "Average"]

# Check that the names match the alpha estimates
length(locus_names)
nrow(sz95$center)

alpha <- sz95$center

alpha_results <- data.frame(
  Locus = locus_names,
  Alpha_median = alpha[, "50%"],
  CI_lower = alpha[, "2.5%"],
  CI_upper = alpha[, "97.5%"],
  Alpha_outlier = alpha[, "2.5%"] > 0.5 |
    alpha[, "97.5%"] < 0.5
)

alpha_results$Direction <- ifelse(
  alpha_results$CI_lower > 0.5,
  "Above 0.5",
  ifelse(
    alpha_results$CI_upper < 0.5,
    "Below 0.5",
    "Not flagged"
  )
)

alpha_outlier_table <- subset(alpha_results, Alpha_outlier)

write.csv(
  alpha_outlier_table,
  "alpha_outlier_loci_Tcri.csv",
  row.names = FALSE
)
