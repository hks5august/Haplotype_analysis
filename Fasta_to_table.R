library (devtools)
library (tidyverse)
source_url("https://raw.githubusercontent.com/lrjoshi/FastaTabular/master/fasta_and_tabular.R")


#sequence to table
FastaToTabular("Popgendata_BegCape_combined.seq")

# Table to sequence 
TabularToFasta("all_333seq.csv")
