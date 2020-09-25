library(tidyverse)
library(data.table)
pmppm <- function(mass, ppm=4){
  if(mass<200){
    c(mass-200*ppm/1000000, mass+200*ppm/1000000)
  } else {
    mass*c(1-ppm/1000000, 1+ppm/1000000)
  }
}

MS1_data <- readRDS("MS1_data_frame") %>% as.data.table()
metadata <- read.csv("falkor_metadata.csv")
mono_peptides <- read.csv("mono_peptides.csv", stringsAsFactors = FALSE)
di_peptides <- read.csv("di_peptides.csv", stringsAsFactors = FALSE)

pdf(file = "mono_peptide_EICs.pdf", height = 3)
for(i in seq_len(nrow(mono_peptides))){
  EIC <- MS1_data[mz%between%pmppm(mono_peptides[i,"value"]+1.007276, ppm = 5)]
  EIC_meta <- left_join(EIC, metadata, by="fileid")
  gp <- ggplot(EIC_meta) +
    geom_line(aes(x=rt, y=int, color=depth, group=fileid)) +
    scale_color_manual(values = c("25m"="blue", "DCM"="green", "Blank"="red")) +
    ggtitle(mono_peptides[i,"mono_name"]) +
    theme_bw()
  print(gp)
}
dev.off()

pdf(file = "di_peptide_EICs.pdf", height = 10, width = 20)
for(i in seq_along(unique(di_peptides$rowname))){
  peps_i <- di_peptides %>%
    filter(rowname==unique(di_peptides$rowname)[i])
  eics <- lapply(seq_len(nrow(peps_i)), function(j){
    EIC <- MS1_data[mz%between%pmppm(peps_i[j,"value"]+1.007276, ppm = 5)]
    cbind(EIC, sec_pep=peps_i[j,"di_name"])
  }) %>% do.call(what = rbind)
  EICs_meta <- left_join(eics, metadata, by="fileid")
  gp <- ggplot(EICs_meta) +
    geom_line(aes(x=rt, y=int, color=depth, group=fileid)) +
    scale_color_manual(values = c("25m"="blue", "DCM"="green", "Blank"="red")) +
    ggtitle(unique(di_peptides[,"rowname"])[i]) +
    theme_bw() + 
    facet_wrap(~sec_pep, scales = "free_y")
  print(gp)
}
dev.off()
