require(tidyverse)
require(cmdstanr)
require(logger)
require(parallel)
                                        #require(rstan)
options(mc.cores = detectCores()-1)
ukb_data <- readRDS("../test/unscaled_test.rds")
#ukb_data <- readRDS("")
ukb_data <- ukb_data %>% filter(iter==500)
hweq <- function(x){sqrt(2*x*(1-x))}
logger::log_info("max maf: {max( ukb_data$MAF)}")
ukb_data <- ukb_data %>% mutate(hw = hweq(MAF))
logger::log_info("sum of squares beta:{sum(ukb_data$u_beta^2)}")
logger::log_info("reading hyper parameters file")
hyper <- read_table("/work/ext-unil-ctgg/tdaniel/genomicarchitecture/data/groups78_mix4_cpus4_tasks8_nodes10_mpisync10_height_1_unrelated_148covadjusted_w35520NA.csvLong",col_names=F)
hyper$X1 <-  as.numeric(str_remove(hyper$X1, ","))
#hyper$X1
hyper<- hyper %>% filter(X1==500)
sigma_g = sqrt(sum(hyper[,3:80]))
sigma_g = sqrt(hyper$X82)
logger::log_info("sigma_g :{sigma_g}")
logger::log_info("J :{nrow(ukb_data)}")
logger::log_info("A : {sum(ukb_data$gen_sd^2)}")

logger::log_info("loading groups")
snp_groups <- readRDS("../data/snp_groups.rds")
ukb_data <- ukb_data %>% inner_join(snp_groups)
#ukb_data  <- ukb_data %>% filter(groups==12)
sigma_g = sqrt(hyper$X82)
iter=4000
options(mc.cores = parallel::detectCores() -1)
set_cmdstan_path(path = "/work/ext-unil-ctgg/tdaniel/genomicarchitecture/cmdstan")
sm <- cmdstan_model("/work/ext-unil-ctgg/tdaniel/genomicarchitecture/models/bayesS_sigmaG.stan")


ukb_list <- list(J= nrow(ukb_data), beta=ukb_data$u_beta, f=ukb_data$MAF, sigma_g=sigma_g)
height_simple<- sm$sa
height_simple<- sm$sample(data=ukb_list,chains=4, iter_warmup=200, iter_sampling=200, refresh=10)
print("S_coef global")
 tryCatch({print(height_simple$summary(c("S_coef")))},error=function(e){print("S_coef NA") })

S_tibble <- height_simple$summary(c("S_coef"))
S_tibble$group <- 79
S_tibble$sample <- 500
S_tibble$chain <- 1
S_tibble$var_g <- sigma_g^2
hyper <- hyper[,3:80]

#sm <- cmdstan_model("/work/ext-unil-ctgg/tdaniel/genomicarchitecture/models/simple.stan")
#sm <- stan_model("simple.stan")
#height_simple<- sm$sample(ukb_data,chains=2,iter_warmup =ceiling(iter/2), iter_sampling = ceiling(iter/2),thin=1, threads_per_chain = 10, refresh=10)
logger::log_info('optimizing')
print("S_coef local")
for(i in 0:77){
    ukb_data_temp <- ukb_data %>% filter(groups==i)
    sigma_g <- sqrt(as.numeric(hyper[,i+1]))
    print(sigma_g)
ukb_list <- list(J= nrow(ukb_data_temp), beta=ukb_data_temp$u_beta, f=ukb_data_temp$MAF, sigma_g=sigma_g)
 height_simple<- sm$sample(data=ukb_list,chains=4, iter_warmup=200, iter_sampling=200, refresh=10)
    tryCatch({print(height_simple$summary(c("S_coef")))},Error=function(e){print("S_coef NA") })
    summary_S <- height_simple$summary(c("S_coef"))
    summary_S$group <- i
    summary_S$sample <- 500
    summary_S$chain <- 1
    summary_S$var_g <- sigma_g^2
    S_tibble <- bind_rows(S_tibble, summary_S)

}

logger::log_info('summarising')

print(S_tibble)
logger::log_info('saving')
saveRDS(S_tibble, "test_500.rds")
