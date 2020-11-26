require(BEDMatrix)
require(logger)
get_sds(bed_file){
  b_matrix <- BEDMatrix::BEDMatrix(bed_file)
  sds <- matrix(0,nrow=ncol(bed_file), ncol=1)
  rownames(sds) <- colnames(b_matrix)
  logger::log_info("calulating sds")
  for(i in 1:ncol(b_matrix)){
    sds[i] <- sd(b_matrix[i],na.rm = T)
    if(i%%100 ==0)
      logger::log_info("{i}/{ncol(b_matrix)}")
  }
    
  out_file <- paste0(bed_file,"_stats.RData")
  save(list="b_matrix", file=out_file)
  logger::log_info("saved {out_file}")
}
args <- commandArgs()
if(length(args)>0){
  get_sds(args[1])
}