### This script will read the sds and maf file (extension .frq) and create a
### a single tibble
require(tidyverse)
require(logger)
require(argparser)

get_maf_sds <- function(maf_file, sds_file, output_file){

    logger::log_info("reading stats file: {sds_file}")
    sds <- readRDS(sds_file)
    rs <-rownames(sds)
    sds <- as_tibble(sds)
    colnames(sds) <- "gen_sd"
    sds$rs <-rs
    sds <- sds %>% tidyr::separate(col=rs,into=c("rs","allele"), sep="_")

    logger::log_info("reading maf file: {maf_file}")
    maf <- read_table(maf_file)
    maf <- maf %>% rename(rs = SNP)
    maf <- maf %>% mutate(CHR=NULL, A1 = NULL, A2=NULL, NCHROBS=NULL)

    logger::log_info("merging sds and maf")
    sds_maf <- sds %>% inner_join(maf) %>% mutate(allele=NULL)
    saveRDS(sds_maf, output_file)
}

### Create a parser
p <- arg_parser("Create MAF sd file")

### Add command line arguments
p <- add_argument(p, "--maf", help=".frq file with the maf per snp")
p <- add_argument(p, "--sd", help=".rds file with the sd per snp")
p <- add_argument(p, "--out", help="output file")

### Parse the command line arguments
argv <- parse_args(p)


get_maf_sds(maf_file=argv$maf, sds_file=argv$sd, output_file=argv$out)
