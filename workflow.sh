### R scripts in the order they should be run in order for the pipleine to function properly. Can run this script to do everything in one go
#Run metamap on the raw notes
Rscript 01_run_metamap.R
Rscript 02_parse_metamap_outputs.R
Rscript 03_run_models.R
Rscript 04_calculate_model_performance.R
