#!/bin/Rscript
# Kirk Amundson
# 16 January 2019

# This is not a finished product, but I have to kick this to the server to have it run overnight.
# need eval runtime throughout

library(tidyverse)
dosage_genos <- read_tsv("2019_0107_LOP_250k_dosage_genotypes_40percband.tsv", col_names = T) %>% 
  mutate(bin = paste(chrom, start, end, sep = "_"))
# filter(bin %in% c("chr01_0_250000", "chr01_250000_500000")) # training wheels, gradually ease off of them
head(dosage_genos)

# for all trisomics, filter out only the excess chromosome
dosage_genos_trifilt <-  dosage_genos %>% 
  filter(!(sample == "2x_LOP868_238" & chrom == "chr01")) %>%
  filter(!(sample == "2x_LOP868_259" & chrom == "chr02")) %>%
  filter(!(sample == "2x_LOP868_292" & chrom == "chr02")) %>%
  filter(!(sample == "2x_LOP868_452" & chrom == "chr02")) %>%
  filter(!(sample == "2x_LOP868_460" & chrom == "chr02")) %>%
  filter(!(sample == "2x_LOP868_1157" & chrom == "chr02")) %>%
  filter(!(sample == "2x_LOP868_1190" & chrom == "chr02")) %>%
  filter(!(sample == "2x_LOP868_373" & chrom == "chr03")) %>%
  filter(!(sample == "2x_LOP868_028" & chrom == "chr04")) %>%
  filter(!(sample == "2x_LOP868_433" & chrom == "chr04")) %>%
  filter(!(sample == "2x_LOP868_293" & chrom == "chr05")) %>%
  filter(!(sample == "2x_LOP868_363" & chrom == "chr05")) %>%
  filter(!(sample == "2x_LOP868_172" & chrom == "chr06")) %>%
  filter(!(sample == "2x_LOP868_065" & chrom == "chr07")) %>%
  filter(!(sample == "2x_LOP868_272" & chrom == "chr07")) %>%
  filter(!(sample == "2x_LOP868_108" & chrom == "chr08")) %>%
  filter(!(sample == "2x_LOP868_488" & chrom == "chr08")) %>%
  filter(!(sample == "2x_LOP868_397" & chrom == "chr10")) %>%
  filter(!(sample == "2x_LOP868_373" & chrom == "chr02")) %>%
  filter(!(sample == "2x_LOP868_128" & chrom == "chr12"))

# split by bin, did this to get contingency tables of dosage genotypes
dosage_genos_split <- split(dosage_genos_trifilt, f=dosage_genos_trifilt$bin)

# remove outlier individuals from data frame to be parsed while constructing fisher contingency tables
# this seems needlessly complicated...

# returns a list of named vectors. names of list are bin names. names of vectors are outlier dosage allele in that bin
outliers <- lapply(dosage_genos_split, function(x) table(x$dosage.gt)[which(table(x$dosage.gt) <= 3)])

# traverse outliers, parsing out bins and outlier dosage allele in each
z <- list()

# loop populates list of outlier bins and dosage alleles
for (i in 1:length(outliers)) {
  if (length(outliers[[i]]) >= 1) { # needs to skip bins without rows
    z[[i]] <- data.frame('dosage.gt' = names(outliers[[i]]),
                         'bin' = names(outliers)[[i]])
  }
}

# bind rows to get single data frame of bin-allele combinations of all outliers
z.done <- bind_rows(z) %>% 
  mutate(allele = paste(bin, dosage.gt, sep = "_"))

# using z.done, remove outlier bin-allele combinations from table of dosage alleles by individual.
dosage_genos_outlierfilt <- dosage_genos_trifilt %>% 
  mutate(allele = paste(bin, dosage.gt, sep = "_")) %>% 
  filter(!allele %in% z.done$allele)

# did we filter anything?
nrow(dosage_genos_trifilt)
nrow(dosage_genos_outlierfilt) # yes

# parsable data structure for nested loops is stored in dosage_genos_outlierfilt

# filter out other unique dosage variants. May need to adjust thresholding due to what I'm guessing are undetected tetraploids in the population
outliers_exclude <- lapply(dosage_genos_split, function(x) table(x$dosage.gt)[which(table(x$dosage.gt) > 3)]) # drops outliers, monomorphic still here
# outliers_exclude

# what's left? This provides the scope of the analysis: the bins that are polymorphic
# lifted from other scripts, may need to develop this further
test <- unlist(lapply(outliers_exclude, function(x) length(x)))
# head(test)
# length(unique(names(which(test >= 2))))
dosage_variable_bins <- unique(names(which(test >= 2)))
# head(dosage_variable_bins)
dosage_invariable_bins <- unique(names(which(test < 2)))
# head(dosage_invariable_bins)

# filter out monomorphic bins
monomorph_exclude <- outliers_exclude[which(names(outliers_exclude) %in% dosage_variable_bins)]

# Generalizable conversion of genotype summary data to list of summary data
geno.list <- list() # stores genotype count data for each nonoverlapping bin

# work in progress chunk for generating geno.list using lapply call or loop
# for (i in 1:length(geno_summaries)) { 
for (i in 1:length(monomorph_exclude)) { 
  geno.list[[i]] <- data.frame("bin" = names(monomorph_exclude)[i], monomorph_exclude[[i]]) %>% 
    mutate(Var1 = as.character(Var1)) %>% 
    mutate(bin_allele = paste(bin, Var1, sep = "_"))
}

# build fisher comparisons for all pairwise allele comparisons between two polymorphic loci
fisher.list <- list() # stores 2x2 contingency tables, lapply call fisher exact test at the end.

# nonredundant all by all pairwise comparison of data frames in list. each fisher contingency table saved as a named entry in the list.
for (x in seq(1, length(geno.list)-1)) { # loop through all except last bins. these are "source"
  # print(x)
  for (y in seq(x+1, length(geno.list))) { # all pairwise "target" bins, given the same source bin
    # print(y)
    for (i in 1:nrow(geno.list[[x]])) { # loop all alleles in source bin at hand
      for (j in 1:nrow(geno.list[[y]])) { # loop all alleles in target bin at hand
        source.bin <- geno.list[[x]][i,1]
        source.gt  <- geno.list[[x]][i,2]
        target.bin <- geno.list[[y]][j,1]
        target.gt  <- geno.list[[y]][j,2]
        # define individuals in source bin
        match.source.allele <- dosage_genos_outlierfilt %>%
          filter(bin == source.bin & dosage.gt == source.gt)
        target.allele.dist <- dosage_genos_outlierfilt %>%
          filter(bin == target.bin & sample %in% match.source.allele$sample) # why is this failing? It worked before!
        g <- table(target.allele.dist$dosage.gt)
        target.cluster.in.count <- g[which(names(g) == target.gt)]
        target.cluster.out.count <- sum(g) - target.cluster.in.count
        target.cluster.null <- round(mean(g))
        fisher.table <- matrix(c(target.cluster.in.count, target.cluster.out.count, target.cluster.null, target.cluster.null), nrow = 2)
        list.name <- paste(geno.list[[x]]$bin_allele[i], geno.list[[y]]$bin_allele[j], sep = " ")
        fisher.list[[list.name]] <- fisher.table
      }
    }
  }
}

# TODO: Parse Fisher exact output, return as flat table so I can finish making new correlation matrix
fisher.price <- lapply(fisher.list, function(x) fisher.test(x))
fisher.price[[1]]

# init empty data frame
fisher.df <- data.frame("src" = rep(NA, length(fisher.list)),
                        "trgt" = rep(NA, length(fisher.list)),
                        "pval" = rep(NA, length(fisher.list)))

# populate
for (i in 1:length(fisher.price)) {
  fisher.df$src[i]  <- str_split(names(fisher.price)[i], " ")[[1]][1]
  fisher.df$trgt[i] <- str_split(names(fisher.price)[i], " ")[[1]][2]
  fisher.df$pval[i] <- fisher.price[[i]]$p.value
}

# a few more cleanups
fisher.df.tidy <- fisher.df %>% 
  separate(src, into = paste0("src.", c("chrom", "start", "end", "gt")), remove = F, sep= "_") %>% 
  separate(trgt, into = paste0("trgt.", c("chrom", "start", "end", "gt")), remove = F, sep = "_") %>% 
  mutate(src.bin = str_replace_all(src, "(.+)_[0-9]$", "\\1")) %>% 
  mutate(trgt.bin = str_replace_all(trgt, "(.+)_[0-9]$", "\\1")) 

# arrange(src.chrom, src.start, src.end, trgt.chsrc.gt, trgt.gt)
head(fisher.df.tidy)

fisher.df.tidy$src.bin.num <- rep(NA, nrow(fisher.df.tidy))
fisher.df.tidy$trgt.bin.num <- rep(NA, nrow(fisher.df.tidy)) 

for (i in 1:nrow(fisher.df.tidy)) {
  fisher.df.tidy$src.bin.num[i] <- which(unique(dosage_genos$bin) == fisher.df.tidy$src.bin[i])
  fisher.df.tidy$trgt.bin.num[i] <- which(unique(dosage_genos$bin) == fisher.df.tidy$trgt.bin[i])
}

p <- ggplot(fisher.df.tidy, aes(x = src.bin.num, y = trgt.bin.num, color = pval)) +
  geom_tile()
ggsave("2019_1116_fishermtx.pdf", width = 6, height = 6, units = "in", device = "pdf")
