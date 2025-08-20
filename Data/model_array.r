#------------------------------------------------------------------------------#
# creating input file ID for submitting array jobs on HPC EVE
# require data stan_dat_1.RData as input file, 
# models coded in Stan and a parameter files 
# every job is a run on a model with specific structure 
# in total 18  models were tested 

#------------------------------------------------------------------------------#

rm(list=ls())

# load required libraries
library(data.table)
library(docopt)


library(rstan)
library(coda)
rstan_options(auto_write = TRUE)
options(mc.cores = 9)

load("stan_dat.RData")

# provide initial values for stan model
inits <- rep(list(list(sigma = 1.0)), 3)

doc <- "Usage: model_array.R <params>"
opts <- docopt(doc)
# data_path <- "/work/yuw/inter-div"


process_single <- function(params) {
  ## actual processing goes here
  ## you can access the parameters from the parameter file
  ## like you would normally with readCSV or data.table:
  
  # to recognize data.1 as data not character (no need here)
  # dat <- eval(parse(text = params$data)) 
  stan_model <- params$stan_model
  res <- params$result

  # compile model and sampling
  stan_fit <- sampling(  stan_model(file = stan_model),
                         data = stan_dat,
                         iter = 12000,
                         warmup = 6000,
                         thin = 1,
                         chains = 3, 
                         init = inits,
                         control = list(max_treedepth = 12,
                                        adapt_delta = 0.9))
  
  # save the output stan_fit object
  saveRDS(stan_fit, file = res)

}


## read parameter file
params <- fread(opts$params)

## try to get task id
task <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))

## check if task id exists
if (is.na(task)) {
  ## process everything
  result <- by(params, 1:nrow(params), process_single)
  attributes(result) <- NULL
  
  # saveRDS(result, file.path(opts$output_dir, "everything.rds"))
} else {
  ## process single item
  result <- process_single(params[task])
  
  # saveRDS(result, file.path(opts$output_dir, paste0("chunk-", task, ".rds")))
}
