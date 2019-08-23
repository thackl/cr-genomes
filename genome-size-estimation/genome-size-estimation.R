#+TITLE: Workpack genome-size-estimation of project cr-genomes
#+AUTHOR: Thomas Hackl
#+DATE: 2019-07-28
#+DESCRIPTION: Based on 19-kmer-size dist for Mavirus reactivation Nature manuscript ("CrA_01.raw-histo.tsv")

library(tidyverse)
library(thacklr)

d1 <- read_table2("CrE410P-ms-raw-kmer-histo.tsv", col_names=c("Coverage", "Frequency"));

peak_covs <- c(60,120,180)
peaks <- tibble(peaks=c("1n", "2n", "3n"), Coverage=peak_covs)

d1 %>% filter(Coverage >= 40 & Coverage <=2000) %>%
  mutate(bp = Coverage * Frequency) %>%
  pull(bp) %>% sum() / 120

gg <- ggplot(d1, aes(x=Coverage, y=Frequency)) +
  geom_line() + 
  xlim(limits=c(0, 300)) +
  ylim(limits=c(0, 1e+06)) +
  geom_text(aes(Coverage, y=300000, label=peaks), peaks) +
  annotate(geom="label", x=120, y=.8e6, label="40 Mbp") +
  theme_bw() + no_grid()
gg

ggsave("CrE410P-kmer-spectrum.pdf", width=6, height=3)
ggsave("CrE410P-kmer-spectrum.png", type="cairo", width=6, height=3)

