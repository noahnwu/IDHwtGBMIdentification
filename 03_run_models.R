################### INSTALL AND ATTACH REQUIRED LIBRARIES ################### 
ini = readRDS("ini.RDS")
print("Running Regularized Models...")
print("Importing Required Libraries...")
#note the addition of the glmnet package into the required packages for this script. This package 
#comes with a lot of dependencies, including Rcppeigen, which requires Rtools for installation on a Windows machine
#contact Noah if there is any trouble with package installation and/or with the import into working environment
required.packages = c("dplyr", "tidyr", "ggplot2", "readr", 
                      "tibble", "purrr", "stringr", "forcats", "glmnet")

missing.packages = required.packages[!(required.packages %in% installed.packages())]
if(length(missing.packages) > 0) { 
  install.packages(missing.packages, dependencies = T)
}
lapply(as.list(required.packages), library, character.only = T, quietly = T)



############# IMPORT DATAFRAME ################
df = ini$dichot_df %>% 
  data.table::fread(.) 
#I had to split documents into smaller chunks for MetaMap due to memory constraints on my machine at NM
#The coefficients were then re-aggregated here by an identifier that I assigned. I'm assuming 
#you won't need this, so I commented it out
# %>%
#   mutate(file.id = str_remove(file.id, "_\\d+\\.txt\\.out")) %>%
#   group_by(file.id) %>%
#   summarise(across(.cols = !matches("file\\.id"), sum)) %>%
#   ungroup() %>%
#   mutate(across(.cols = !matches("file\\.id"), .fns = ~ifelse(.x >=1 , 1, 0)))

################# DEFINE FUNCTION FOR RUNNING MODEL ###############
run.mod = function(model.path, 
                   analytic.df = df, 
                   #probability threshold for model to determine positive class membership
                   prob.cut = 0.5, 
                   #this is hardcoded in previous script and so should be present as a column
                   id_col = "file.id", 
                   #positive class label
                   pos.class = "IDH WT GBM", 
                   neg.class = "NOT IDH WT GBM") { 
  mod = readRDS(model.path)
  #all variables used by the models (these need to be included in the dataframe for us to be able to use it)
  mod.vars = mod$glmnet.fit$beta %>% 
    as.matrix() %>%  
    rownames()
  #variables not currently included in our dataframe
  unident.vars = mod.vars[!(mod.vars %in% names(analytic.df))]
  #create a matrix of 0s to append to the current dataframe with appropriate column names
  unident.mat = matrix(data = 0, nrow = nrow(analytic.df), ncol = length(unident.vars)) %>% 
    as.data.frame() %>%  
    rename_with(.fn = ~unident.vars)
  #append the matrix of 0s to allow for prediction
  X = analytic.df %>% 
    column_to_rownames(id_col) %>% 
    bind_cols(unident.mat) %>% 
    select(mod.vars) %>% 
    as.matrix()
  prob = predict(mod, newx = X, type = "response") %>%  
    as.numeric()
  df.out = X %>%  
    as.data.frame() %>% 
    rownames_to_column() %>%  
    select(rowname) %>%  
    mutate(probability = prob, 
           prediction = ifelse(prob >= prob.cut, pos.class, neg.class)) %>%  
    rename(file.id = rowname)
  return(df.out)
}

################# RUN IDH WT IDENTIFICATION MODELS ON DATASET ##################
pred.gbm = run.mod(ini$idh_mod, 
                   analytic.df = df) 


############## RUN MGMT SUBTYPE IDENTIFICATION MODEL ON DATASET ############
#Step 1. Identify Notes in Which MGMT is Missing (ONLY RUN ON NOTES THAT WERE IDENTIFIED AS HAVING POSITIVE IDH WT)

idh.wt.gbm = df %>%  
  filter(file.id %in% (
    pred.gbm %>%  
      filter(prediction == "IDH WT GBM") %>%  
      pull(file.id)
  ))
if(nrow(idh.wt.gbm) > 0) { 
  pred.mgmt.missing = run.mod(ini$mgmt.miss_mod, 
                              idh.wt.gbm,
                              pos.class = "MISSING MGMT", 
                              neg.class = "NOT MISSING MGMT")
  mgmt.nomiss = idh.wt.gbm %>% 
    filter(file.id %in% (
      pred.mgmt.missing %>% 
        filter(prediction == "NOT MISSING MGMT") %>% 
        pull(file.id)
    ))
  #Step 2 Identify MGMT Subtype (If there are any that were predicted to be non-missing)
  if(nrow(mgmt.nomiss) > 0) { 
    pred.mgmt.subtype = run.mod(ini$mgmt.subtype_mod, 
                                analytic.df = mgmt.nomiss,
                                pos.class = "METHYLATED", 
                                neg.class = "UNMETHYLATED")
  }
  
} 

if(!exists("pred.mgmt.missing")) { 
  pred.mgmt.missing = data.frame(file.id = df$file.id, 
                                 probability = NA,
                                 prediction = NA)
}
if(!exists("pred.mgmt.subtype")) { 
  pred.mgmt.subtype = data.frame(file.id = df$file.id, 
                                 probability = NA, 
                                 prediction = NA)
}

################### JOIN ALL PREDICTIONS ######################
df.pred.all = df %>%  
  select(file.id) %>% 
  left_join(pred.gbm, 
            by = "file.id") %>%  
  left_join(pred.mgmt.missing, 
            by = "file.id", 
            suffix = c("_idhwt", "_MGMTmissing")) %>%  
  left_join(pred.mgmt.subtype %>%  
              rename_with(.cols = !matches("file.id"), ~paste0(.x, "_MGMTsubtype")))
#output file
write_csv(df.pred.all, ini$model.prediction_df)











