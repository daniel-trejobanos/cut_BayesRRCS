library(tidyverse)

height <- read_delim("../data/height_beta_means.txt",delim = "\t")

height %>% filter(means!=0)

sds <- readRDS("../data/ukb_imp_v3_UKB_EST_oct19_unrelated_full_stats.rds")
rs <-rownames(sds)
sds <- as_tibble(sds)
sds$rs <-rs
sds_tmp <- sds %>% tidyr::separate(col=rs,into=c("rs","allele"), sep="_")
mean_effects <- sds_tmp %>% dplyr::inner_join(height)

unscaled_effects <-mean_effects %>% mutate(beta=means/V1)

maf <- read_table("../data/ukb_imp_v3_UKB_EST_oct19_unrelated.frq")
maf <- maf %>% rename(rs = SNP)

unscaled_effects <- unscaled_effects %>% inner_join(maf)


#saveRDS(unscaled_effects,file = "effects_maf.rds")

data_ukb <- list(beta=unscaled_effects$beta, f=unscaled_effects$MAF, FST=0.01, J = nrow(unscaled_effects))
saveRDS(data_ukb, file="../data/data_mean_ukb_height.rds")

