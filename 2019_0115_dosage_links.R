#!/bin/Rscript
# Kirk Amundson
# 16 January 2019

# This is not a finished product, but I have to kick this to the server to have it run overnight.
# need eval runtime throughout

library(tidyverse)

# read in dosage genotpye data, obtained by mean shift clustering of relative coverage values in nonoverlapping bins
dosage_genos <- read_tsv("2019_0107_LOP_250k_dosage_genotypes_40percband.tsv", col_names = T) %>% 
  mutate(bin = paste(chrom, start, end, sep = "_")) # %>% 
  # filter(chrom %in% c("chr08", "chr09")) # test case using 2 chromosomes. Expect strong linkage within chromosomes but not between

# tally by dosage allele in each bin
# see here for inspiration: http://www.stat.wisc.edu/~yandell/R_for_data_sciences/curate/tidyverse.html
dosage_geno_counts <- dosage_genos %>% 
  count(bin, cluster, dosage.gt) %>% 
  filter(n > 3) %>% 
  filter(n < 164)

# set up all pairwise allele comparisons
# commented out distinct() call removes redundant comparisons
dosage_geno_comps <- dosage_geno_counts %>% 
  setNames(paste0("target.", names(.))) %>% 
  crossing(dosage_geno_counts) %>%
  filter(target.bin != bin & target.dosage.gt != dosage.gt) %>%
  mutate(bincomp = ifelse(target.bin > bin, paste(target.bin, target.dosage.gt, bin, dosage.gt, sep = "_"), paste(bin, dosage.gt, target.bin, target.dosage.gt, sep = "_"))) %>% 
  distinct(bincomp, .keep_all = T) # here, .keep_all=T specifies to keep all variables. Not doing this makes the computation take that much longer


dosage_geno_comps$inds.in.target <- rep(NA, nrow(dosage_geno_comps))
dosage_geno_comps$inds.out.target <- rep(NA, nrow(dosage_geno_comps))

# get counts of individuals matching allele under comparison and not matching allele under comparison in target bin.
# is a whale of a for loop
ptm <- proc.time()
for (i in 1:nrow(dosage_geno_comps)) {
  # for (i in 1:10) { # training
  source.bin <- dosage_geno_comps$bin[i]
  source.allele <- dosage_geno_comps$dosage.gt[i]
  target.bin <- dosage_geno_comps$target.bin[i]
  target.allele <- dosage_geno_comps$target.dosage.gt[i]
  inds <- filter(dosage_genos, bin == source.bin & dosage.gt == source.allele)$sample
  con <- filter(dosage_genos, sample %in% inds & bin == target.bin)
  dosage_geno_comps$inds.in.target[i] <- nrow(filter(con, dosage.gt == target.allele))
  dosage_geno_comps$inds.out.target[i] <- nrow(filter(con, dosage.gt != target.allele))
}
proc.time() - ptm

# get data frame of observed target allele counts in each bin
dosage_geno_counts_split <- split(dosage_geno_counts, f=dosage_geno_counts$bin)
bin_allele_counts <- unlist(lapply(dosage_geno_counts_split, function(x) nrow(x)))
bin_allele_counts <- data.frame("target.bin" = names(bin_allele_counts),
                                "target.allele.count" = bin_allele_counts)

# add target observed allele count and target expected allele count
fisher <- left_join(dosage_geno_comps, bin_allele_counts, by = "target.bin") %>% 
  mutate(target.expect = round((inds.in.target + inds.out.target) / target.allele.count))

# do fisher exact test
fisher$fisher.pval <- apply(fisher, 1, function(x) {
  tbl <- matrix(as.numeric(c(x[10], x[11], x[13], x[13])) , ncol = 2, byrow = T)
  fisher.test(tbl)$p.value
})

# finally, add numeric bin back
fisher$src.bin.num <- rep(NA, nrow(fisher))
fisher$target.bin.num <- rep(NA, nrow(fisher))
dg.bins <- unique(dosage_genos$bin)
for (i in 1:nrow(fisher)) {
  fisher$src.bin.num[i] <- which(dg.bins == fisher$bin[i])
  fisher$target.bin.num[i] <- which(dg.bins == fisher$target.bin[i])
}

# for graphing purposes, may want to check src.bin.num and target.bin.num, and reverse source-target so that all tiles show up on the same side of the diagonal
# tested this using a single chromosome, works. Only plots below the diagonal.
fisher <- fisher %>%
  mutate(src.bin.num.plot = ifelse(src.bin.num > target.bin.num, src.bin.num, target.bin.num)) %>%  
  mutate(target.bin.num.plot = ifelse(src.bin.num > target.bin.num, target.bin.num, src.bin.num)) # target plot always lesser of src.bin.num and target.bin.num

write_tsv(fisher, "2019_0118_fisher_pairwise_linkages.tsv", col_names = T, delim = "\t")

# draw plot
# p <- ggplot(fisher, aes(x = src.bin.num.plot, y = target.bin.num.plot, color = fisher.pval)) +
#   geom_tile()
# ggsave("2019_1116_fishermtx.pdf", width = 6, height = 6, units = "in", device = "pdf")