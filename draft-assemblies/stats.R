library(tidyverse)

a0 <- read_tsv(pipe("seqkit stats -T -a *.fa"))

a0 %>%
  mutate(
    assembler = case_when(
      str_sub(file, 6,6) == "c" ~ "Canu 1.8",
      str_sub(file, 6,6) == "f" ~ "Flye 2.3.7",
      str_sub(file, 6,6) == "w" ~ "Wtdbg 2.1",
      TRUE ~ "SPAdes 3.6.1")
  ) %>%
    arrange(file) %>%
    select(file, assembler, everything(), -format, -type, -avg_len, -starts_with("Q"), -sum_gap) %>%
    write_tsv("assembly-stats.tsv")
