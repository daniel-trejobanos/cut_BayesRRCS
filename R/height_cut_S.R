require(tidyverse)
require(cmdstanr)
require(logger)
require(parallel)
require(argparser)

## scaling according to maf
hweq <- function(x){sqrt(2*x*(1-x))}
logger::log_threshold(TRACE)
logger::log_info("setting cmdstan")
options(mc.cores = parallel::detectCores() -1)
set_cmdstan_path(path = "/work/ext-unil-ctgg/tdaniel/genomicarchitecture/cmdstan")

cut_S <- function(beta_file, hyper_file, chain, out_file, model="bayesS_sigmaG"){

    logger::log_info("openning: {beta_file}")
    ukb_data  <-  readRDS(beta_file)

    logger::log_info("calculating maf scaling")
    ukb_data <- ukb_data %>% mutate(hw = hweq(MAF))

    logger::log_info("reading hyper parameters file: {hyper_file}")
    hyper <- read_table(hyper_file, col_names=F)
    hyper$X1 <-  as.numeric(str_remove(hyper$X1, ","))

    hyper <- hyper %>% filter(X1>= 2500)
    ## we select the sigmaG for each group
    hyper_groups <- hyper[,3:80]
    colnames(hyper_groups) <- 1:ncol(hyper)
    n_groups <- ncol(hyper_groups)
    logger::log_trace("n groups: {n_groups}")
    hyper_iter  <-  hyper$X1/5

    hyper_groups$X1 <- hyper_iter

    print(hyper_iter)
    if(model=="bayesS_sigmaG"){
        sm <- cmdstan_model("/work/ext-unil-ctgg/tdaniel/genomicarchitecture/models/bayesS_sigmaG.stan")
    }else{
        if(model=="bayesS"){
            sm <- cmdstan_model("/work/ext-unil-ctgg/tdaniel/genomicarchitecture/models/cut_bayesS.stan")
        }else{
            stop("wrong model specified, otions are bayesS or bayesS_sigmaG")
        }

    }

    result_S  <- NULL
    for(i in hyper_iter ){
        logger::log_info("loading groups")
        snp_groups <- readRDS("../data/snp_groups.rds")
        logger::log_trace("reading betas for iter: {i}")
        cur_data <- ukb_data %>% filter(iter == i) %>% inner_join(snp_groups)
        logger::log_trace("betas read: {length(cur_data$iter)}")
        for(j in 1:77){
            cur_data_temp <- cur_data %>% filter(groups == j)
            logger::log_trace("reading betas for group {j}")
            sigma_g <- sqrt(as.numeric( hyper_groups %>% filter(X1==i) %>% select(j)))
            logger::log_trace("sigma_G: {sigma_g}")
            if(sigma_g == 0){stop("ERROR: sigma_G eq. 0")}

            cur_data_list  <- list(J= nrow(cur_data_temp), beta=cur_data_temp$u_beta, f=cur_data_temp$MAF, sigma_g=sigma_g)
            logger::log_trace("running stan")
            posterior<- sm$sample(data=cur_data_list,chains=4, iter_warmup=400, iter_sampling=400, refresh=10)

            S <- posterior$draws()[dim(posterior$draws())[1], , "S_coef" ]

            S_summary <- posterior$summary(c("S_coef"))
            S_summary$group <- j
            S_summary$var_g <- sigma_g^2
            S_summary$chain_sigmaG <- chain
            S_summary$Schain_1_S <- S[1]
            S_summary$Schain_2_S <- S[2]
            S_summary$Schain_3_S <- S[3]
            S_summary$Schain_4_S <- S[4]

            result_S  <- bind_rows(result_S,S_summary)
        }
        if(i%%2==0){
        	saveRDS(result_S,out_file)
        }
    }

    saveRDS(result_S,out_file)
}

### Create a parser
p <- arg_parser("Run cutBayesS")

### Add command line arguments
p <- add_argument(p, "--beta_file", help="file path of chain with unscaled betas (omitting the chain number)")
p <- add_argument(p, "--sigma_file", help= "file path of the hyperparameters (sigmaG and sigmaE)" )
p <- add_argument(p, "--chain", help="number of chain to be processed (HAS to be the same as the beta and hyper file)")
p <- add_argument(p, "--out", help="output file for the chain summary")
### Parse the command line arguments
argv <- parse_args(p)

cut_S(argv$beta_file, argv$sigma_file, argv$chain, argv$out)
