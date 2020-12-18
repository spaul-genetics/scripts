# scripts
Here I will keep and share scripts that I used to analyze genetics and genomics data. 

### `extract_ukb_fields.R` 

Given a list of UK Biobank field IDs in a text file, the script will extract all array indices. For example if f.20002.0 is in the `--fields` file then it will extract f.20002.0.0 to f.20002.0.33. The `--fields` should contain UKB Data Field ID in this format: `f.XXXX.X`. One ID in a row.  


### `qc_ukb_impute_release3.py`

A python script to perform QC steps on UK Biobank imputed data based on the `.mfi`, and PLINK output file from `--hardy`. For details on input parameters use `python qc_ukb_impute_release3.py --help`. 

