library(rvest)
library(tidyverse)

aa_raw <- read_html(x = "http://education.expasy.org/student_projects/isotopident/htdocs/aa-list.html")
html_nodes(aa_raw, "td")

table_titles <- html_nodes(aa_raw, "th") %>% html_text() %>% as.character()
aa_clean <- aa_raw %>%
  html_nodes("td") %>% 
  html_text() %>% 
  matrix(ncol=5, byrow=TRUE) %>%
  as.data.frame(stringsAsFactors=FALSE) %>%
  `names<-`(table_titles) %>%
  #select(c(2,4)) %>%
  mutate(Monoisotopic=as.numeric(Monoisotopic)+18.0105+0.00005)
aa_map <- pull(aa_clean, 2) %>% `names<-`(pull(aa_clean, 1))
mono_peptides <- select(aa_clean, mono_name=`3-letter code`, value=Monoisotopic)

di_peptides <- "https://www.ionsource.com/tutorial/DeNovo/denovo_tables.htm" %>%
  read_html() %>%
  html_nodes(xpath = '//table[@width="1154"]') %>%
  html_table(fill = TRUE) %>%
  as.data.frame() %>%
  slice(-2) %>%
  select(-2) %>%
  `colnames<-`(.[1,]) %>%
  select(-1) %>%
  slice(-1) %>%
  mutate(rowname=colnames(.)) %>%
  select(-"CMC") %>%
  filter(rowname!="CMC") %>%
  pivot_longer(cols = -last_col()) %>%
  mutate(rowname=recode(rowname, !!!aa_map)) %>%
  mutate(name=recode(name, !!!aa_map)) %>%
  mutate(value=as.numeric(value)+18.0105) %>%
  mutate(di_name=paste(rowname, name, sep = "-"))

write.csv(mono_peptides, file = "mono_peptides.csv", row.names = FALSE)
write.csv(di_peptides, file = "di_peptides.csv", row.names = FALSE)
