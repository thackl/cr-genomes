library(phytools)
library(phyloseq)
library(patchwork)
library(tidyverse)
library(ggtree)
#library(thacklr)
library(devtools)
devtools::load_all("~/Code/projects/thacklr", export_all=FALSE)
library(reshape2)

asms <- c("CrEa-c0c", "CrBV-c0c", "CrCf-c0c", "CrRC-c0c",
          "CrEa-c1c", "CrBV-c1c", "CrCf-c1c", "CrRC-c1c")

min_length <- 5000

dev.new()
asm <- asms[4+4]
path <- paste0("/MIT-Research/CroV-complex/cr-genomes/pacbio-based-assemblies/decontaminate/", asm,"/")
analyze_assembly(asm, min_length)



for (asm in asms){
  path <- paste0("/MIT-Research/CroV-complex/cr-genomes/pacbio-based-assemblies/decontaminate/", asm,"/")
  analyze_assembly(asm, min_length)
}
  

#------------------------------------------------------------------------------#

asm <- asms[4]
mt_file <- paste0(path, asm, "-mt-hits.paf");
mt0 <- read_paf(mt_file) %>%
  rename(seq_id=target_name)
mt0 %>% select(query_name)

asm <- asms[4]
path <- paste0("/MIT-Research/CroV-complex/cr-genomes/pacbio-based-assemblies/decontaminate/", asm,"/")
emale_file <- paste0(path, asm, "-emale-hits.paf");
emale0 <- read_paf(emale_file) %>%
  rename(seq_id=target_name)
emale0 %>% select(query_name)



#------------------------------------------------------------------------------#

clean_y_axis <- function(){
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())}

cor_dist <- function(m, ...) as.dist((1-(cor(t(m), ...)))/2)

analyze_assembly <- function(asm, min_length, diff_cov=FALSE){

  # bam-coverage now reports length
  #ln0 <- read_tsv(paste0(path, asm, ".fa.fai" ), qc(seq_id, length, i, dont, care)) %>%
  #  select(1:2) %>% filter(length >= min_length)

  # differential coverage -------------------------------------------------------#
  cov_files <- list.files(path, "*-cov.tsv", full.names=T)
  names(cov_files) <- cov_files %>% basename %>% str_match("(Cr...?-..)-cov.tsv") %>% `[`(,2)
  cov_files
  #
  cv_0 <- map_dfr(
    .id="sample", cov_files, read_tsv,
    col_names = c("seq_id", "length", "cov_mean", "cov_med"))
  cv_1 <- cv_0 %>% filter(length>=min_length)
  #cv_1 %>% arrange(-cov_med) %>% ggplot() + geom_point(aes(x=rank(cov_med), y=cov_med)) + ylim(0,200)
  cov_files
  # differential coverage ------------------------------------------------------#
  # needs multiple bam/cov.tsv files
  if(diff_cov){
    if(length(cov_files) < 2){
      warning("need at least 2 cov files for differential cov !!!!")}
    cd_1 <- cv_1 %>% select(-length, -cov_mean) %>% spread(sample, cov_med) %>%
      mutate_if(is.numeric, replace_na, 0) %>%
      rowwise() %>% do({
        samples <- names(.)[-1]
        tibble(
        seq_id=.$seq_id,
          samples=paste(samples[-length(samples)], samples[-1], sep=":"),
          cov_diff=diff(log(unlist(.[-1])+1)))
      })
    #
    cd_1
  }

  ##cv_1 %<>% arrange(-length) %>% mutate(length_cum=cumsum(length))
  ##cv_1 %<>% mutate(length_cat = case_when(
  ##               length_cum > sum(length*.90) ~ ">90%",
  ##               length_cum > sum(length*.75) ~ "<90%",
  ##               length_cum > sum(length*.50) ~ "<75%",
  ##               length_cum > sum(length*.25) ~ "<50%",
  ##               length_cum > sum(length*.1) ~ "<25%",
  ##               TRUE ~ "<10%"
  ## ))


  # # differential coverage
  # ggplot(cd_1) + geom_tile(aes(samples, seq_id, fill=cov_diff)) +
  #  scale_fill_distiller(palette="Spectral")

  # taxonomy from kaiju ----------------------------------------------------------#
  tax_file <- paste0(path, asm, "-tax-nr-euk.tsv" );
  ranks <- c('species', 'genus', 'family', 'order', 'class', 'phylum', 'kingdom', 'domain')

  tx_0 <- read_tsv(tax_file, col_names=c("UC", "seq_id", "kaiju", rev(ranks))) %>%
    extract(seq_id, c("seq_id", "window"),"^(.*)_(\\d+)", convert=TRUE)
  tx_0
  tx_1 <- tx_0 %>% count(seq_id, domain) %>% spread(domain,n)  %>%
    rename(Eukaryota=`2759`, Bacteria=`2`, Archaea=`2157`, Viruses=`10239`) %>%
    mutate(total = rowSums(select(., `<NA>`, Eukaryota, Bacteria, Archaea, Viruses), na.rm=TRUE)) %>%
    gather(domain, count, -seq_id, -total, -`<NA>`) %>%
    mutate_if(is.numeric, replace_na, 0) %>%
    mutate(count_rel = count/total) %>%
    arrange(`<NA>`/total) %>%
    mutate(seq_id = factor(.$seq_id, unique(.$seq_id)))

  tx_col <- gg_color_hue(5)
  names(tx_col) <- c("Viruses", "Mitochondrion", "Bacteria", "Eukaryota", "Archaea")

  # nucleotide composition -------------------------------------------------------#
  comp_file <- paste0(path, asm, "-comp.tsv" );
  cp_0 <- read_tsv(comp_file) %>%
    filter(window=="median" & seq_id %in% cv_1$seq_id)

  # GC
  gc_0 <- cp_0 %>% select(seq_id, A, C) %>% mutate(GC=C/(A+C))

  # tetra nucs
  c4_0 <- cp_0 %>% select(-window,-size, -A,-C,-AA,-AC,-AG,-AT,-CA,-CC,-CG,-GA,-GC,-TA)
  # c4_0
  kmers_ordered <- c4_0 %>% select(-seq_id) %>% colSums %>% sort(TRUE) %>% names
  c4_1 <- c4_0 %>% gather(kmer, freq, -seq_id) %>%
    group_by(seq_id)
  c4_1$kmer <- factor(c4_1$kmer, levels=kmers_ordered)
  # normalizations
  c4_1 %<>% mutate(
    freq_rel = freq/sum(freq),
    freq_log = log(freq+1) - log(sum(freq)) # with pseudo-count
  )
  #
  comp_file
  c4_1


  #ggplot(c4_1) + geom_tile(aes(kmer, seq_id, fill=freq_log)) +
  #  scale_fill_distiller(palette="Spectral")

  ## CLUSTERING ##################################################################
  # dist based ------------------------------------------------------------------#

  ## # aitchison (better for compositional data) - still not ideal
  ## # also doesn't work for freq==0, gives NA
  ## library(coda.base)
  ## set_attribute <- function(x, attribute, value){attr(x, attribute) <- value; x}
  ## dnm <- acast(dn1, label ~ dinuc, value.var="freq_rel")
  ## tr_freq_rel <- dnm %>% dist(method="aitchison")
  ##   # fix for coda.base::dist (fails to set names for method="aitchison")
  ## tr_freq_rel %<>% set_attribute("Labels", rownames(dnm))
  ## tr_freq_rel %<>% hclust("ward.D2") %>% as.phylo
  ## m <- matrix(1:10/10, byrow=T, nrow=2)
  ## dist(m, "aitchison")


  # euk dist
  ## dist_freq_rel <- acast(c4_1, seq_id ~ kmer, value.var="freq_log") %>%
  ##   dist(method="euc")
  ## dist_cov_diff <- acast(cd_1, seq_id ~ samples, value.var="cov_diff") %>%
  ##   dist(method="euc")

  # cor based - kendall > spearman > pearson -------------------------------#
  # usually not the best idea for compos. data, but does a good job here
  dist_freq_rel <- acast(c4_1, seq_id ~ kmer, value.var="freq_log") %>%
    cor_dist(method="spearman")
  #dist_cov_diff <- acast(cd_1, label ~ samples, value.var="cov_diff") %>%
  #  cor_dist()
  #dist_freq_rel %>% hclust("ward.D2") %>% cutree(plot
  tr_freq_rel <- dist_freq_rel %>% hclust("complete") %>% as.phylo
  ## tr_cov_diff <- dist_cov_diff %>% hclust("ward.D2") %>% as.phylo
  ## #
  ## dist_freq_cov <- foo <- 1 * dist_freq_rel + 0 * dist_cov_diff
  ## tr_freq_cov <- dist_freq_cov %>% hclust("ward.D2") %>% as.phylo

  gg_tr <- ggtree(tr_freq_rel) + scale_y_tree()


  # annotate groups -------------------------------------------------------#
  mt_file <- paste0(path, asm, "-mt-hits.paf");
  mt0 <- read_paf(mt_file) %>%
    rename(seq_id=query_name) %>%
    filter(seq_id %in% cv_1$seq_id) # ignore seqs we've already filtered out

  emale_file <- paste0(path, asm, "-emale-hits.paf");
  emale0 <- read_paf(emale_file) %>%
    rename(seq_id=query_name) %>%
    filter(seq_id %in% cv_1$seq_id) # ignore seqs we've already filtered out

  an0 <- cv_1 %>% select(seq_id) %>% unique %>%
    mutate(
      annot = case_when(
        seq_id %in% filter(tx_1, domain=="Bacteria" & count_rel >.5)$seq_id ~ "Bacterial contamination",
        seq_id %in% mt0$seq_id ~ "Mitochondrion",
        seq_id %in% emale0$seq_id ~ "EMALE",
        TRUE ~ "C.r."
      )) %>%
        rename(label=seq_id)

  mt_file
  
  an_col <- tribble(
    ~annot, ~color,
    "Bacterial contamination", tx_col["Bacteria"],
    "Mitochondrion", tx_col["Mitochondrion"],
    "EMALE", tx_col["Viruses"])

  an1 <- an0 %>% mutate(y=tree_y(gg_tr, .)) %>%
    filter(annot != 'C.r.') %>%
    left_join(an_col)

  anns <- annotate("tile", y=an1$y, x=0.5, width=Inf, height=1, fill=an1$color, alpha=.7)
  anns_flipped <- annotate("tile", x=an1$y, y=0.5, height=Inf, width=1, fill=an1$color, alpha=.7)


  bc_file <- paste0(path, asm, "-contaminated.tsv");

  cont <- an1 %>% filter(annot=="Bacterial contamination") %>%
    select(seq_id=label, reason=annot) %>%
    bind_rows(cv_0 %>% filter(length<min_length) %>%
    mutate(reason=paste0("too short (", length,")")) %>%
    select(seq_id, reason) %>% unique)

  write_tsv(cont, bc_file)


  # composite plot ---------------------------------------------------------#
  theme_set(theme_classic())
  # GC
  #gg_heat_gc <- ggtreeplot(gg_tr, ungroup(gc_0) %>% rename(label=seq_id), aes("GC")) +
  #  geom_tile(aes(fill=GC)) + scale_fill_distiller(palette="Spectral") +
  #  theme(legend.position="bottom") + no_y_axis()
  #
  gg_gc <- ggtreeplot(gg_tr, ungroup(gc_0) %>% rename(label=seq_id), aes(GC)) +
    anns + geom_point() + scale_fill_distiller(palette="Spectral") +
    theme(legend.position="bottom") + clean_y_axis()
  # tetra nuc
  gg_heat_freq <- ggtreeplot(gg_tr, ungroup(c4_1) %>% rename(label=seq_id), aes(kmer)) + geom_tile(aes(fill=freq_log)) + scale_fill_distiller(palette="Spectral") + no_legend() + no_y_axis() +
    theme(axis.text.x=element_blank())
  #
  #gg_heat_cov <- ggtreeplot(gg_tr, ungroup(cd_1) %>% rename(label=seq_id), aes(samples)) + geom_tile(aes(fill=cov_diff)) + scale_fill_distiller(palette="Spectral") + no_legend()
  #gg_heat_cov <- ggtreeplot(gg_tr, ungroup(cv_1) %>% rename(label=seq_id), aes(cov_sample)) + geom_tile(aes(fill=log(cov_med+1))) + scale_fill_distiller(palette="Spectral")
  #
  ## #gg_cov_len <- ggtreeplot(gg_tr, cv_1 %>% rename(label=seq_id), aes(x=1)) +
  ##   facet_wrap(~sample, nrow=1) +
  ##   geom_segment(aes(x=1, xend=cov_med, yend=y), color="grey50", size=1) +
  ##   geom_segment(aes(x=1, xend=length/2000, yend=y, color=length_cat), size=.5) +
  ##   geom_point(aes(x=cov_med), size=1, shape=19) +
  ##   geom_point(aes(x=length/2000,color=length_cat), size=1, shape=19) +
  ##   theme(legend.position="bottom") + coord_cartesian(xlim=c(0,500)) + xlim(1,NA) + theme_void() + anns
  #
  gg_len <- ggtreeplot(gg_tr, cv_1 %>% rename(label=seq_id) %>% filter(sample==sample[1]), aes(y=length), flip=TRUE) +
    anns_flipped +
    geom_col() + coord_flip() + clean_y_axis() +
    theme(legend.position="bottom") + scale_fill_grey()
  # gg_len
  anns_facet <- annotate("tile", x=rep(an1$y,2), y=0.5,, height=Inf, width=1, fill=rep(an1$color,2), alpha=.1)
  #anns_facet$data <- map_df(unique(cv_1$sample), ~mutate(anns_flipped$data, sample=.x))
  gg_cov <- ggtreeplot(gg_tr, cv_1 %>% rename(label=seq_id) %>% filter(grepl("-ms$",sample)), aes(y=cov_med), flip=TRUE) + anns_flipped +
    geom_col() + coord_flip(ylim=c(0,200)) + clean_y_axis() +
    theme(legend.position="bottom") + scale_fill_grey()
  # gg_cov
  gg_tax <- ggtreeplot(gg_tr, rename(tx_1, label=seq_id), aes(y=count_rel), flip=TRUE)  + anns_flipped +
    geom_col(aes(fill=domain)) + coord_flip() + clean_y_axis() +
    theme(legend.position="bottom") +
    scale_fill_manual(values=tx_col)
  #gg_tr + ggtitle(asm) + gg_heat_freq + gg_cov_len + gg_heat_cov + plot_layout(widths=c(2,3,3,.5), nrow=1)

  gg <- gg_tr + anns + gg_heat_freq + gg_len + gg_cov + gg_gc + gg_tax +
    plot_layout(widths=c(1,1,1,1,1,1), nrow=1) + plot_annotation(tag_levels="a", title=asm)

  ggsave(paste0(asm, "-tetranuc.pdf"), gg, width=20, height=20)
  gg
}

