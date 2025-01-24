# This file parses the outputs from MetaMap and Creates Two 
# Versions of the Analytic DataFrame that Will be Utilized by 
# the Regularized Logistic Regression Models
ini = readRDS("ini.RDS")

########## DECLARING DIRECTORY ##########
#the directory in which the metamap outputs are saved
metamap.out.dir = ini$metmap.out
metamap.files = dir(metamap.out.dir, full.names = T)

#check that files are uniquely identified 
if(length(unique(metamap.files)) != length(metamap.files)) { 
  stop("Metamap Output Files Must have Unique Names in Order to Proceed")
}


################### INSTALL AND ATTACH REQUIRED LIBRARIES ################### 
print("Parsing MetaMap Output...")
print("Importing Required Libraries...")
required.packages = c("dplyr", "tidyr", "ggplot2", "readr", 
                      "tibble", "purrr", "stringr", "forcats")
missing.packages = required.packages[!(required.packages %in% installed.packages())]
if(length(missing.packages) > 0) { 
  install.packages(missing.packages)
}
lapply(as.list(required.packages), library, character.only = T, quietly = T)


################### PARSE NOTES AND CREATE ANALYTIC DATAFRAMES ###################

#we will create two versions of the analytic dataframe, 
#one that contains counts of the concepts within a document
#and another that is dichotomous (present/absent)
#this is because two versions of the models were created at NM
#that relied on two different versions of the dataset.

#this object contains the version with concept counts 
#will create dichotomous version at end of script
print("Parsing MetaMap Outputs")
df.count = vector("list", length = length(metamap.files))
for(i in 1:length(metamap.files)) { 
  print(i)
  df.count[[i]] = tryCatch(
    expr = {
      read_delim(metamap.files[i], delim = "|", 
                 col_names =F, 
                 show_col_types = F, 
                 progress = F
                 ) %>% 
        #keep only the columns that we will need for later analysis
        select(5, 7) %>% 
        rename_with(.fn = ~c("CUI", "trigger")) %>% 
        #there are some lines (usually 2-3 at the beginning of the output file) that do not represent
        #the identification of a CUI. These are removed here
        filter(str_detect(CUI, "^C")) %>%  
        #determine whether MetaMap identified term negation
        #1 if negated, 0 if not
        mutate(negation_flag = str_extract(trigger, "\\d\\]$"), 
               negation_flag = str_remove(negation_flag, "\\]"), 
               negation_flag = as.numeric(negation_flag), 
               CUI_final = ifelse(negation_flag == 1, paste0(CUI, "_neg"), CUI)) %>%  
        group_by(CUI_final) %>% 
        count() %>%  
        ungroup() %>% 
        #put note output in WIDE format
        pivot_wider(names_from = CUI_final, 
                    values_from = n) %>% 
        mutate(file.id = metamap.files[i]) %>%  
        relocate(file.id)
    }, 
    error=function(e) {
      print(e)
    }
  )
}
df.count = lapply(df.count, function(x) { 
  if("data.frame" %in% class(x)) { 
    return(x)
  } else{
    return(NULL)  
    }
  
  }) %>% 
  reduce(., bind_rows) %>%  
  mutate(across(!matches("file\\.id"), ~ifelse(is.na(.x), 0, .x)))

#if the count of a concept is greater than 0, then change to 1 (dichotomous version)
df.dichot = df.count %>% 
  mutate(across(!matches("file\\.id"), ~ifelse(.x > 0, 1, 0)))

write_csv(df.count, ini$counts_df)
write_csv(df.dichot, ini$dichot_df)

