# MODELS FOR IDENTIFYING IDHWT GBM AND MGMT PROMOTER METHYLATION FROM SURGICAL PATHOLOGY DOCUMENTS


Below is an explanation of the setup needed as well as files contained within the directory. Contact noah.forrest@northwestern.edu if any issues arise.

# SYSTEM INFORMATION

The current project was developed and tested on the following system: 
Description:	Ubuntu 22.04 LTS
Release:	22.04
Codename:	jammy

All R scripts were written and tested with: 
platform       x86_64-pc-linux-gnu         
arch           x86_64                      
os             linux-gnu                   
system         x86_64, linux-gnu           
status                                     
major          4                           
minor          3.1                         
year           2023                        
month          06                          
day            16                          
svn rev        84548                       
language       R                           
version.string R version 4.3.1 (2023-06-16)
nickname       Beagle Scouts        


# SETUP  

## Three setup steps are needed in order to run the scripts contained in this repository

### Step 1. Set the MetaMap bin path
	
	- All R scripts reference the 'ini.RDS' file for paths to files in the workflow
	- ini.RDS can be read by R with the readRDS() function, which will return a list
	- the metamap bin location is indexed by the name 'metamap.loc'
	- default value is ~/Desktop/metamap/bin
	- to change this value, enter the following code into R: 

	ini = readRDS("ini.RDS")
	ini[["metamap.loc"]] = "/your/path/here" #note forward slashes are preferred by R. You can also escape backslashes.
	saveRDS(ini, "ini.RDS")

### Step 2. Import your Gold Standard Dataset
	- This file is referened by default in the root project directory and should be named "adjudicated_gold_standard.csv". If you want to 
	  change the file that is referenced by the scripts, you can also change the ini.RDS, index 'gold.standard_df'.
	- File should be comma separated and contain the following columns 
		- file.id: the (unique) identifier for a file, extracted in later scripts from the name of the text file in input_files
		- idhwt: the labels for whether a ddocument was identified as positive for IDHWT GBM. Possible Values: "IDH WT GBM" or "NOT IDH WT GBM" (ALL CAPS)
		- MGMTmissing: the label for whether a document was missing information on the MGMT promoter methylation status. Possible Values: "MISSING MGMT" or "NOT MISSING MGMT"
		- MGMTsubtype: the label for the MGMT promoter methylation status. Possible Values: "METHYLATED" or "UNMETHYLATED"
### Step 3. Import your raw text files into input_files/ directory
	- MetaMap requires that all text be enclosed by single quotation marks. 
	- Each text file must be uniquely identified. R will throw an error if there are duplicate identifiers. 
	- Text file identifiers must match the file.id columnfound in the gold standard .csv mentioned above (no need to add the .txt.out at the end of the file.id

Following these setup steps, the R scripts can be run sequentially with the workflow.sh file (or equivalent .bat file if on a windows OS)

# SCRIPTS
All scripts are numbered by the order in which they were intended to be run
Any output file names or paths can be ammended by editing the "ini.R" file in a similar manner to that described in SETUP

### 01_run_metamap.R 

This script first identifies all raw text files in the input_files directory
MetaMap is then called from within the R script and runs CUI identification of each document in series
All MetaMap outputs are written  to the metmap_output/ directory (mispelled but this is correct) with the file naming as "text file name.txt.out"
See https://lhncbc.nlm.nih.gov/ii/tools/MetaMap/Docs/MMI_Output.pdf for details on the format for these output files, 
but this is hopefully irrelevant to you if the following script works as intended

### 02_parse_metamap_outputs.R 

This script parses the '|' delimited files written by MetaMap and extracts two pieces of information: 
	CUI: Concept Unique Idenntifer
	Negation: Whether the CUI identified was a positive or negative instance of the word
In the case where a negative instance of a CUI is detected, "_neg" is appended to the end of the column name (this was the convention I used when training the model, and
variable names need to be consistent when running the model on new data)

Two files are output: 
"analytic_df_counts.csv"
"analytic_df_dichotomous.csv"

These are two different versions of the same data parsed by metamap. The first contains the count of CUIs identified within a file and the latter is translated to 
dichotomous ("present/absent"). For our purposes, the dichotomous file is the only one we will use, as the best performing models all utilized dichotomous covariates.


### 03_run_models.R 

Runs the regularized models and outputs predictions to "model_performance.csv"
All models referenced are listed in the "ini.RDS" file
Note that the model predictions follow the same convention as the expected values for the Gold Standard Dataset (all caps)

### 04_calculate_model_performance.R 

Manual computation (to reduce dependencies) of model performance in identifying IDH WT GBM, MGMT Promoter Missingnes, and MGMT Promoter Subtype
Compares predictions made by the models to the gold standard dataset. Following statistics are output into the file "model_performance.csv": 

Sensitivity: TP/(TP+FN)
Specificity: TN(TN + FP)
PPV: TP/(TP + FP)
NPV: TN/(TN+FN)
F-measure: 2*[(Sensitivity * PPV)/(Sensitivity + PPV)]



# OTHER FILES & DIRECTORIES 

* idh_models: all models for identifying idhwt gbm
* mgmt_models: contains two subdirectories, one that contains models that identify misssing MGMT instances and the second that identifies MGMT promoter methylation subtype

**Note** All models in these directories are saved as RDS files and are of class cv.glmnet. Note that the glmnet package is required to do essentially anything with these once they 
are loaded into an R environment.

That's it (for now)! Please reach out if any questions come up.
