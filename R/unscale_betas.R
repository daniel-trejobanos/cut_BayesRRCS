### This Rscript will read the betmap file for a given chain and then create the corresponding tibble of unscaled betas
require(tidyverse)
require(logger)
require(argparser)
#' unscale betas
#' based on a current chain file, loads posterior betas from iter_start to iter_end, and
#' multiplies them to the sds in the sds file. Result is a data frame with the sds and maf
#' for the corresponding un-scaled betas.
#' @param current_chain file and path for the chain to process
#' @param iter_start first iteration (inclusive) to process betas from
#' @param iter_end last iteration (inclusive) to process betas from
#' @param dest_file rds file to save the resulting tibble in.
#' @param maf_sds_file file containing a tible, (snp,sd,maf)
unscale_betas <- function(current_chain, iter_start, iter_end, maf_sds_file, dest_file){

    logger::log_info("reading {current_chain}")
    current_samples <- read_tsv(current_chain, col_names=FALSE)
    colnames(current_samples) <- c("iter", "rs","beta")

    logger::log_info("reading maf_sd file {maf_sds_file}")
    maf_sds <- readRDS(maf_sds_file)

    logger::log_info("filtering from iter {iter_start}")
    current_samples <- current_samples %>% filter(iter >= iter_start)
    logger::log_info("filtering to iter {iter_end}")
    current_samples <- current_samples %>% filter(iter <= iter_end)
    logger::log_info("merging samples and maf")
    current_samples <- current_samples %>% inner_join(maf_sds)

    logger::log_info("unscaling betas")
    current_samples <- current_samples %>% mutate(u_beta=beta/gen_sd)
    current_samples <- current_samples %>% mutate(beta=NULL)

    logger::log_info("saving unscaled betas {dest_file}")
    saveRDS(current_samples, dest_file)
}

### Create a parser
p <- arg_parser("Unscale betas")

### Add command line arguments
p <- add_argument(p, "--chain", help="betMap file with iter,snp,beta")
p <- add_argument(p, "--start", help="first iteration (inclusive)", type="numeric")
p <- add_argument(p, "--end", help="last iteration (inclusive)", type="numeric")
p <- add_argument(p, "--mafsd", help="file with the maf-sds tibble")
p <- add_argument(p, "--out", help="output file")

### Parse the command line arguments
argv <- parse_args(p)

### Call function
unscale_betas(current_chain=argv$chain, iter_start=argv$start, iter_end=argv$end, maf_sds_file=argv$mafsd, dest_file=argv$out  )
