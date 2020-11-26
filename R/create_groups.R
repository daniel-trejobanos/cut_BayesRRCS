require(tidyverse)
require(logger)
require(argparser)

create_groups_table  <- function(snp_id_file, groups_file,output){
    snps_id <- read_table(snp_id_file,col_names=F)
    groups <- read_table(groups_file,col_names=F)
    groups_table <- tibble(groups=groups$X1, rs = snps_id$X1)
    saveRDS(groups_table,output)
}

### Create a parser
p <- arg_parser("Create groups table file")

### Add command line arguments
p <- add_argument(p, "--snp_id", help=".txt file with  snps ordered as in the analysis")
p <- add_argument(p, "--groups", help=".group file with groups per snp")
p <- add_argument(p, "--out", help="output file")

### Parse the command line arguments
argv <- parse_args(p)


create_groups_table(snp_id_file=argv$snp_id, groups_file=argv$groups, output=argv$out)
