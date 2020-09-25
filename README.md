## diPeptides

A repository designed to search for and plot chromatograms of monopeptides (aka amino acids) 
and dipeptides in mass-spec data.

### Process

1. Run scrapeScript.R to download masses for amino acids and dipeptides

2. Edit ms1_data_maker.R to point to your mass-spec files (good to test with a subset first)

3. Run ms1_data_maker.R

4. Run makeEICs.R

5. View PDFs

### Inputs and outputs

#### scrapeScript.R

Requires nothing other than an internet connection and the R libraries `rvest` and `tidyverse`

Monopeptide source: http://education.expasy.org/student_projects/isotopident/htdocs/aa-list.html 
but those are hydrated masses so I subtract a water and the mass of an electron

Dipeptide source: https://www.ionsource.com/tutorial/DeNovo/denovo_tables.htm

Outputs two files: mono_peptides.csv and di_peptides.csv which can be reviewed and edited in Excel
if necessary

#### ms1_data_maker.R

Requires mono_peptides.csv and di_peptides.csv

Requires access to mass-spec files (e.g. mapped network drive)

Requires editing to suit your own computer - absolute file paths are used by `mzR` and should not
be assumed by the script

Requires libraries `tidyverse`, `data.table`, `pbapply`, **and `mzR`**

Outputs ms1_data.rds, a compressed version of the data extracted from each mass-spec file. 
Is essentially a data frame with *m/z*, retention time, and intensity columns plus the file source.
This file is ignored by Git because it can easily exceed Github's 25MB limit.

#### makeEICs.R

Requires mono_peptides.csv, di_peptides.csv, ms1_data.rds

Requires libraries `tidyverse` and `data.table`

Outputs two PDFs containing chromatograms for various mono- and di- peptides.
