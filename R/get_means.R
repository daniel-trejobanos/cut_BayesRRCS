library(tidyverse)
library(collapse)
library(logger)
logger::log_info("reading data")
samples <- read_tsv("groups78_mix4_cpus4_tasks8_nodes10_mpisync10_height_1_unrelated_148covadjusted_w35520NA.betMap", col_names=FALSE)
logger::log_info("computing mean betas")
mean_betas_height <- samples %>% collapse::fsubset(X1 >= 500) %>% get_vars(2:3) %>%  collapse::fgroup_by(X2) %>% fmean 
#save(list="mean_betas_height", file="mean_betas_height.RData")
var_height <- samples %>% collapse::fsubset(X1 >= 500) %>% get_vars(c(1,3)) %>% fgroup_by(X1) %>% fvar
nobs_height <- samples %>% collapse::fsubset(X1 >= 500) %>% get_vars(c(1,3)) %>%  fgroup_by(X1) %>% fNobs
save(list=c( "var_height" , "nobs_height"), file="summary_height.RData")
