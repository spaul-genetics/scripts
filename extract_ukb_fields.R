# Given a list of UK Biobank field IDs in a text file, the script will extract all array indices
# For example if f.20002.0 is in the --fields file then it will extract f.20002.0.0 to f.20002.0.33
# Written by Subrata
# Date: Noveber 15, 2020
library(optparse)
library(data.table)
library(stringr)
parser <- OptionParser()
parser <- add_option(parser, c("-i", "--input"), type = "character",
                            default = '/work/IBG/UKBiobank/Keller-1665/phenotypes.all/ukb42049.tab',
                            help = 'The .tab file of UKB phenotypes')
parser <- add_option(parser, c("-o","--output"), type = "character",
                            default = 'output_phenotypes.tab', help = 'Output file name with location')
parser <- add_option(parser, c("-f", "--fields"), type = "character",
                            help = "Give a text file with the field ids one in a line in the format of f.XXXXX.X")
args = parse_args(parser)

all_fields = system (paste0("head -n 1 ",args$input), intern = T)
all_fields = str_split(all_fields, '\t', simplify = T)
print(paste0('There are ', length(all_fields), 'data fields in the input file'))
interest_fields = fread(args$fields, header = F)

extract_cols = grep(paste(interest_fields$V1, collapse = '|'), all_fields, value = T)
extract_cols = c('f.eid', extract_cols)
print(paste0('Number of fields in interest including f.eid = ', length(extract_cols)))

dat = fread(args$input, header = T, select = extract_cols)
fwrite(dat, args$output, row.names = F, col.names = T, sep = '\t', quote = F)
