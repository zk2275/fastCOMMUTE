.libPaths('library.path')
library(glmnet)

M <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID")) #read in array id
args = commandArgs(trailingOnly = T) #read in all items specified in the .sh file
outcome = as.character(args[[1]]) #set the first args item

print(outcome)

sim_setting <- expand.grid(sim = 1:100, method = 1:2)

sim = sim_setting[M, 'sim']
method_type = sim_setting[M, 'method']

print(sim)
print(method_type)

set.seed(sim)