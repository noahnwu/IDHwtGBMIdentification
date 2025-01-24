ini = readRDS("ini.RDS")
print("Computing Model Performance...")
required.packages = c("dplyr", "tidyr", "ggplot2", "readr", 
                      "tibble", "purrr", "stringr", "forcats")
#glmnet not required by this script
missing.packages = required.packages[!(required.packages %in% installed.packages())]
if(length(missing.packages) > 0) { 
  install.packages(missing.packages, dependencies = T)
}
lapply(as.list(required.packages), library, character.only = T, quietly = T)

##################### READ IN DATA FRAMES #########################

#output from 03_run_models.R
pred = read_csv(ini$model.prediction_df)
#manually adjudicated labels from chart review
labels = read_csv(ini$gold.standard_df)
################### DECLARE FUNCTION ##########################3
calc_performance = function(
  #dataset produced by models in 03_run_models.R
  #class data.frame
  predictions, 
  #dataset with labels assigned through manual adjudication
  #class data.frame
  gold_standard, 
  #column that represents the column identifier (hardcoded from predictions dataset)
  id.col = "file.id", 
  #column from predictions dataset that contains model predictions on variable of interest
  pred.col, 
  #column from manually adjudicated gold standard representing the same variable of interest
  gs.col, 
  #positive class
  pos.class, 
  #negative class
  neg.class
  ) { 
  #ensures that patients are in the same order
  predictions[[id.col]] = str_remove(predictions[[id.col]], "\\.txt\\.out")
  labels[[id.col]] = str_remove(labels[[id.col]], "\\.txt\\.out")
  df.joined = predictions %>%  
    filter(!is.na(!!sym(pred.col))) %>% 
    #changed 1/16/2024 due to package version incompatibilities
    select(all_of(c(id.col, pred.col))) %>% 
    inner_join(gold_standard %>%  
                 filter(!is.na(!!sym(gs.col))) %>% 
                 #changed 1/16/2024 due to package version incompatibilities
                 select(all_of(c(id.col, gs.col))), 
               by = id.col) %>%  
    #changed 1/16/2024 due to package version incompatibilities
    select(all_of(c(gs.col, pred.col)))
  if(nrow(df.joined) == 0) { 
    stop("0 predictions made for this variable")
  }
  twobytwo = table(df.joined[[gs.col]], df.joined[[pred.col]]) %>%  
    as.data.frame() %>%  
    rename_with(.fn = ~c("gold", "model", "n")) %>%  
    mutate(across(gold:model, ~case_when(
      !(.x %in% c(pos.class, neg.class)) ~ NA, 
      .x == pos.class ~ 1, 
      .x == neg.class ~ 0)
    )
    )
  check.na = sum(apply(twobytwo, 2, function(x) sum(is.na(x))))
  if(check.na > 0) { 
    stop("Positive or Negative Class Misspecified (check for spelling of class labels)")
  }
  cells = list("a" = twobytwo$n[twobytwo$gold ==1 & twobytwo$model == 1],
                        "b" = twobytwo$n[twobytwo$gold ==0 & twobytwo$model == 1],
                        "c" = twobytwo$n[twobytwo$gold ==1 & twobytwo$model == 0],
                        "d" = twobytwo$n[twobytwo$gold ==0 & twobytwo$model == 0]) %>%  
    lapply(., function(x) { 
      ifelse(is_empty(x), 0, x)
      })
  performance = list(
    "Sensitivity" = cells$a/sum(cells$a, cells$c), 
    "Specificity" = cells$d/sum(cells$d, cells$b), 
    "PPV" = cells$a/sum(cells$a, cells$b), 
    "NPV" = cells$d/sum(cells$d, cells$c)
  )
  performance[["F-measure"]] = 2*((performance$Sensitivity * performance$PPV)/(performance$PPV + performance$Sensitivity))
  return(list(
    #keep twobytwo table in case there are any issues
    "twobytwo" = twobytwo, 
    "pos.class" = pos.class, 
    "neg.class" = neg.class,
    "performance" = performance
  ))
}

############################### CALCULATE PERFORMANCE IN FINDING IDHWT GBM ################################
idhwt.performance = calc_performance(predictions = pred, 
                 gold_standard = labels, 
                 pred.col = "prediction_idhwt", 
                 gs.col = "idhwt", 
                 pos.class = "IDH WT GBM", 
                 neg.class = "NOT IDH WT GBM")

############################# CALCULATE PERFORMANCE IN IDENTIFYING MISSING MGMT (IDHWT GBM ONLY) ##################################
mgmt.missing.performance = calc_performance(predictions = pred, 
                                            gold_standard = labels, 
                                            pred.col = "prediction_MGMTmissing", 
                                            gs.col = "MGMTmissing", 
                                            pos.class = "MISSING MGMT", 
                                            neg.class = "NOT MISSING MGMT")
########################### CALCULATE MODEL PERFORMANCE OF IDENTIFYING MGMT SUBTYPE #######################################
mgmt.subtype.performance = calc_performance(predictions = pred, 
                                            gold_standard = labels, 
                                            pred.col = "prediction_MGMTsubtype", 
                                            gs.col = "MGMTsubtype", 
                                            pos.class = "METHYLATED", 
                                            neg.class = "UNMETHYLATED")
####################### OUTPUT RESULTS ###########################
list(idhwt.performance, 
     mgmt.missing.performance, 
     mgmt.subtype.performance) %>%  
  lapply(., function(x) { 
    N = sum(x$twobytwo$n)
    x[["performance"]] %>%
      as.data.frame() %>% 
      mutate(N_documents = N, 
             `Positive Class` = x$pos.class, 
             `Negative Class` = x$neg.class)
    }) %>% 
  purrr::reduce(., bind_rows) %>%  
  mutate(Model = c("IDH WT GBM", 
                   "Missing MGMT", 
                   "MGMT Subtype")) %>%  
  relocate(Model) %>%  
  write_csv(., ini$model.performance_df)

######################## SCRAPPED CODE ##########################33
#### this was used to simulate a gold standard for performance calcualtion development
# set.seed(101397)
# labels = pred %>%  
#   select(file.id, matches("prediction")) %>%  
#   rename_with(.cols = matches("prediction"), ~str_remove(.x, "prediction_"))
# pred.cols = names(labels)[!str_detect(names(labels), "file\\.id")]
# pred.opt = list(unique(labels$idhwt), 
#                 c("NOT MISSING MGMT", "MISSING MGMT"), 
#                 c("METHYLATED", "UNMETHYLATED"))
# names(pred.opt) = c("idh", "mgmt.miss", "mgmt.subtype")
# prob = c(0.6, 0.4)
# rev.prob = c(0.4, 0.6)
# for(i in 1:nrow(labels)) { 
#   if(labels$idhwt[i] == "IDH WT GBM") { 
#     temp.idh = sample(pred.opt$idh, size = 1, prob = rev.prob)
#   } else { 
#     temp.idh = sample(pred.opt$idh, size = 1, prob = prob)
#   }
#   if(temp.idh == "NOT IDH WT GBM") { 
#     labels$idhwt[i] = temp.idh
#     labels$MGMTmissing[i] = NA
#     labels$MGMTsubtype[i] = NA
#   } else { 
#       if(labels$MGMTmissing[i] == "NOT MISSING MGMT" | is.na(labels$MGMTmissing[i])) { 
#         temp.mgmt.miss = sample(pred.opt$mgmt.miss, size = 1, prob = c(0.95, 0.05))
#       } else { 
#         temp.mgmt.miss = sample(pred.opt$mgmt.miss, size = 1, prob = c(0.05, 0.95))
#       }
#     if(temp.mgmt.miss == "MISSING MGMT") { 
#       labels$idhwt[i] = temp.idh
#       labels$MGMTmissing[i] = temp.mgmt.miss
#       labels$MGMTsubtype[i] = NA
#     } else { 
#       if(labels$MGMTsubtype[i] == "METHYLATED" | is.na(labels$MGMTsubtype[i])) { 
#         temp.mgmt = sample(pred.opt$mgmt.subtype, size = 1, prob = prob)
#       } else { 
#         temp.mgmt = sample(pred.opt$mgmt.subtype, size = 1, prob = rev.prob)  
#       }
#       labels$idhwt[i] = temp.idh
#       labels$MGMTmissing[i] = temp.mgmt.miss
#       labels$MGMTsubtype[i] = temp.mgmt
#       }
#     }
# }

