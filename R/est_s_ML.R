### here we will estimate s and sigma
require(cmdstanr)
require(logger)
require(tidyverse)
est_s_ML <- function(iters_file, model, output_file)
