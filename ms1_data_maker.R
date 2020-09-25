
# Setup things ----
library(tidyverse)
library(data.table)
library(pbapply)



# Functions ----
pmppm <- function(mass, ppm=4){
  if(mass<200){
    as.numeric(mass)+(c(-ppm, ppm)*200/1000000)
  } else {
    c(mass*(1-ppm/1000000), mass*(1+ppm/1000000))
  }
}
grabSingleFileData <- function(filename){
  msdata <- mzR:::openMSfile(filename)
  fullhd <- mzR::header(msdata)
  spectra_list <- lapply(seq_len(nrow(fullhd)), function(x){
    given_peaks <- mzR::peaks(msdata, x)
    rtime <- fullhd[x, "retentionTime"]
    return(cbind(rtime, given_peaks))
  })
  all_data <- `names<-`(as.data.frame(do.call(rbind, spectra_list)), 
                        c("rt", "mz", "int"))
  
  return(all_data)
}



# Find paths to files ----
ms_file_paths <- dir(r"(G:\My Drive\FalkorFactor\mzMLs\pos)",
                     pattern = ".mzML",
                     full.names = TRUE) %>%
  normalizePath() %>%
  `[`(1:3)



# Read in raw data after it's been created by scrapeScript.R ----
mono_peptides <- read.csv("mono_peptides.csv")
di_peptides <- read.csv("di_peptides.csv")
pep_masses <- c(mono_peptides$value, di_peptides$value)



# Grab data points within ppm of the peptide mass across all RTs
stan_MS1_data <- pblapply(ms_file_paths, function(file_path){
  # Open up the mzML as a data.frame with mz, rt, int columns
  raw_EIC <- grabSingleFileData(file_path)
  
  # Convert to data.table to (greatly) speed up subsetting
  raw_EIC_dt <- as.data.table(raw_EIC)
  
  # For each mono/di peptide, grab the parts of the mzML that have an m/z 
  # close to (within `ppm`) the peptide mass
  lapply(unique(pep_masses), function(mz_i){
    raw_EIC_dt[mz%between%pmppm(as.numeric(mz_i), ppm=5)]
  }) %>% do.call(what = rbind) %>%
    cbind(fileid=basename(file_path))
}) %>% do.call(what = rbind)



# Write out the data into an .rds object for high-fidelity encoding and speedy reading
saveRDS(stan_MS1_data, file = "stan_data.rds")
