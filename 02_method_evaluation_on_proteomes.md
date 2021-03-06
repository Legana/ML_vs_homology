
# BLAST vs.¬†classification models for finding AMPs in proteomes

Read in AMP database

``` r
uniprot_and_amp_dbs_amps <- readRDS("data/uniprot_amps_w_amp_dbsJuly21.rds") %>%
  rename(Entry_name = `Entry name`) %>% 
  mutate(Organism = str_remove(Organism, " \\(.*")) %>% 
  rename(Taxonomic_lineage = `Taxonomic lineage (ALL)`) %>% 
  rename(Order = `Taxonomic lineage (ORDER)`) %>% 
  rename(Phylum = `Taxonomic lineage (PHYLUM)`) %>%
  rename(Class = `Taxonomic lineage (CLASS)`) %>%
  mutate(Order = str_remove(Order, " \\(.*")) %>%
  mutate(Organism = str_replace_all(Organism, " ", "_")) %>%
  ungroup()
```

## Proteomes used

The selected organisms from were used to create different BLAST query
sets and classification models for each organism (see
05_amp_training_data_preparations.Rmd)

**Table 6.1:** Proteomes used with protein and gene count obtained from
UniProt and AMP count obtained from the AMP database

| Organism Name                         | Reference proteome ID                                        | Total proteins | AMPs | Gene count |
|---------------------------------------|--------------------------------------------------------------|----------------|------|------------|
| *Mus musculus* (mouse)                | [UP000000589](https://www.uniprot.org/proteomes/UP000000589) | 55,366         | 104  | 22,001     |
| *Homo sapiens* (human)                | [UP000005640](https://www.uniprot.org/proteomes/UP000005640) | 78,120         | 96   | 20,600     |
| *Bos taurus* (cattle)                 | [UP000009136](https://www.uniprot.org/proteomes/UP000009136) | 37,513         | 58   | 23,847     |
| *Oryctolagus cuniculus* (rabbit)      | [UP000001811](https://www.uniprot.org/proteomes/UP000001811) | 41,459         | 17   | 21,193     |
| *Ornithorhynchus anatinus* (platypus) | [UP000002279](https://www.uniprot.org/proteomes/UP000002279) | 32,824         | 11   | 17,390     |
| *Gallus gallus* (jungefowl)           | [UP000000539](https://www.uniprot.org/proteomes/UP000000539) | 27,535         | 25   | 18,113     |
| *Oncorhynchus mykiss* (trout)         | [UP000193380](https://www.uniprot.org/proteomes/UP000193380) | 46,447         | 12   | 46,405     |
| *Drosophila melanogaster* (fruitfly)  | [UP000000803](https://www.uniprot.org/proteomes/UP000000803) | 22,110         | 23   | 13,821     |
| *Penaeus vannamei* (shrimp)           | [UP000283509](https://www.uniprot.org/proteomes/UP000283509) | 25,399         | 18   | 25,399     |
| *Bombyx mori* (moth)                  | [UP000005204](https://www.uniprot.org/proteomes/UP000005204) | 14,776         | 15   | 14,773     |
| *Arabidopsis thaliana* (cress)        | [UP000006548](https://www.uniprot.org/proteomes/UP000006548) | 39,337         | 294  | 27,468     |
| *Lithobates catesbeianus* (frog)      | [UP000228934](https://www.uniprot.org/proteomes/UP000228934) | 28,218         | 13   | 28,218     |
| *Escherichia coli K-12* (bacteria)    | [UP000000625](https://www.uniprot.org/proteomes/UP000000625) | 4,438          | 29   | 4,392      |

## Read in proteomes from organisms

``` r
read_proteome_metadata <- function(path, organism) {
  read_tsv(path, col_types = cols(), show_col_types = FALSE) %>%
  rename("Entry_name" = `Entry name`) %>%
  mutate(Organism = organism) %>% 
  mutate(Label = case_when(str_detect(Keywords, "Antimicrobial") ~ "Pos", TRUE ~ "Neg"))
}

mouse_proteome_metadata <- read_proteome_metadata("data/proteomes/M_musculus-proteome-UP000000589.tab.gz", "Mus_musculus")
cow_proteome_metadata <- read_proteome_metadata("data/proteomes/B_taurus-proteome-UP000009136.tab.gz", "Bos_taurus")
human_proteome_metadata <- read_proteome_metadata("data/proteomes/H_sapiens-proteome-UP000005640.tab.gz", "Homo_sapiens")
rabbit_proteome_metadata <- read_proteome_metadata("data/proteomes/O_cuniculus-proteome-UP000001811.tab.gz", "Oryctolagus_cuniculus")
platypus_proteome_metadata <- read_proteome_metadata("data/proteomes/O_anatinus-proteome-UP000002279.tab.gz", "Ornithorhynchus_anatinus")
junglefowl_proteome_metadata <- read_proteome_metadata("data/proteomes/G_gallus-proteome-UP000000539.tab.gz", "Gallus_gallus")
trout_proteome_metadata <- read_proteome_metadata("data/proteomes/O_mykiss-proteome-UP000193380.tab.gz", "Oncorhynchus_mykiss")
fruitfly_proteome_metadata <- read_proteome_metadata("data/proteomes/D_melanogaster-proteome-UP000000803.tab.gz", "Drosophila_melanogaster")
shrimp_proteome_metadata <- read_proteome_metadata("data/proteomes/P_vannamei-proteome-UP000283509.tab.gz", "Penaeus_vannamei")
moth_proteome_metadata <- read_proteome_metadata("data/proteomes/B_mori-proteome-UP000005204.tab.gz", "Bombyx_mori")
cress_proteome_metadata <- read_proteome_metadata("data/proteomes/A_thaliana_proteome-UP000006548.tab.gz", "Arabidopsis_thaliana")
frog_proteome_metadata <- read_proteome_metadata("data/proteomes/L_catesbeianus-proteome-UP000228934.tab.gz", "Lithobates_catesbeianus")
bacteria_proteome_metadata <- read_proteome_metadata("data/proteomes/E_coli-proteome-UP000000625.tab.gz", "Escherichia_coli")
```

## Overlap between AMPs in AMP database and AMPs present in proteomes

``` r
amp_overlap_in_proteome_and_amp_db <- function(proteome_metadata, organism){
  proteome_metadata %>% 
    filter(Label == "Pos") %>% 
    mutate(in_amp_database = Entry %in% uniprot_and_amp_dbs_amps$Entry) %>%  
    summarise(AMPs_in_proteome = n(), AMPs_overlap_in_AMP_db = sum(in_amp_database)) %>% 
    mutate(Organism = organism)
}

mouse_amp_overlap <- amp_overlap_in_proteome_and_amp_db(mouse_proteome_metadata, "Mus musculus")
human_amp_overlap <- amp_overlap_in_proteome_and_amp_db(human_proteome_metadata, "Homo sapiens")
cow_amp_overlap <- amp_overlap_in_proteome_and_amp_db(cow_proteome_metadata, "Bos taurus")
rabbit_amp_overlap <- amp_overlap_in_proteome_and_amp_db(rabbit_proteome_metadata, "Oryctolagus cuniculus")
platypus_amp_overlap <- amp_overlap_in_proteome_and_amp_db(platypus_proteome_metadata, "Ornithorhynchus anatinus")
junglefowl_amp_overlap <- amp_overlap_in_proteome_and_amp_db(junglefowl_proteome_metadata, "Gallus gallus")
trout_amp_overlap <- amp_overlap_in_proteome_and_amp_db(trout_proteome_metadata, "Oncorhynchus mykiss")
fruitfly_amp_overlap <- amp_overlap_in_proteome_and_amp_db(fruitfly_proteome_metadata, "Drosophila melanogaster")
shrimp_amp_overlap <- amp_overlap_in_proteome_and_amp_db(shrimp_proteome_metadata, "Penaeus vannamei")
moth_amp_overlap <- amp_overlap_in_proteome_and_amp_db(moth_proteome_metadata, "Bombyx mori")
cress_amp_overlap <- amp_overlap_in_proteome_and_amp_db(cress_proteome_metadata, "Arabidopsis thaliana")
frog_amp_overlap <- amp_overlap_in_proteome_and_amp_db(frog_proteome_metadata, "Lithobates catesbeianus")
bacteria_amp_overlap <- amp_overlap_in_proteome_and_amp_db(bacteria_proteome_metadata, "Escherichia coli")

amp_overlap <- rbind(mouse_amp_overlap, human_amp_overlap, cow_amp_overlap, rabbit_amp_overlap, platypus_amp_overlap, junglefowl_amp_overlap, trout_amp_overlap, fruitfly_amp_overlap, shrimp_amp_overlap, moth_amp_overlap, cress_amp_overlap, frog_amp_overlap, bacteria_amp_overlap) %>% mutate(AMPs_in_AMP_db = c(104, 96, 58, 17, 11, 25, 12, 23 , 18, 15, 294, 13, 29))
```

**Table 6.2:** AMP count of 13 organisms in their respective proteomes
and in the AMP database, and the overlap between these AMPs

| Organism                 | AMPs_in_proteome | AMPs_overlap_in_AMP_db | AMPs_in_AMP_db |
|:-------------------------|-----------------:|-----------------------:|---------------:|
| Mus musculus             |              131 |                     99 |            104 |
| Homo sapiens             |              115 |                     95 |             96 |
| Bos taurus               |              116 |                     54 |             58 |
| Oryctolagus cuniculus    |               83 |                     17 |             17 |
| Ornithorhynchus anatinus |               27 |                     11 |             11 |
| Gallus gallus            |               29 |                     25 |             25 |
| Oncorhynchus mykiss      |               15 |                      0 |             12 |
| Drosophila melanogaster  |               30 |                     23 |             23 |
| Penaeus vannamei         |                3 |                      0 |             18 |
| Bombyx mori              |               25 |                     13 |             15 |
| Arabidopsis thaliana     |              294 |                    291 |            294 |
| Lithobates catesbeianus  |               11 |                      0 |             13 |
| Escherichia coli         |                4 |                      4 |             29 |

From the `amp_overlap` table, it can be seen that there are AMPs in the
proteomes that do not overlap with the AMPs in the AMP database. A
possible explanation for this is that the proteome AMPs were found via
the ‚ÄúAntimicrobial‚ÄĚ keyword in UniProt. These contain proteins not
experimentally verified to contain antimicrobial activity and are
therefore not included in the AMP database. Additionally, there are AMPs
in the AMP database not present in the proteomes.

To see if this is due to the presence of mature AMP sequences in the AMP
database, the AMPs in the AMP database corresponding to each organism
were matched to full length sequences, annotated as not antimicrobial,
in the organisms‚Äô respective proteomes. However, only a very few number
of AMPs were detected with this method. The absence of these AMP in the
proteome could be because these AMPs correspond to genes not yet
annotated

``` r
get_missing_amps <- function(proteome_metadata, organism) {
  seqs <- uniprot_and_amp_dbs_amps %>% filter(Organism == organism) %>% select(Sequence)
  proteome_metadata %>% filter(Label == "Neg") %>% filter(str_detect(Sequence, str_c(seqs$Sequence, collapse = "|")))
}

get_missing_amps(mouse_proteome_metadata, "Mus_musculus")
get_missing_amps(human_proteome_metadata, "Homo_sapiens")
get_missing_amps(cow_proteome_metadata, "Bos_taurus")
get_missing_amps(rabbit_proteome_metadata, "Oryctolagus_cuniculus")
get_missing_amps(platypus_proteome_metadata, "Ornithorhynchus_anatinus")
get_missing_amps(junglefowl_proteome_metadata, "Gallus_gallus")
get_missing_amps(trout_proteome_metadata, "Oncorhynchus_mykiss")
get_missing_amps(fruitfly_proteome_metadata, "Drosophila_melanogaster")
get_missing_amps(shrimp_proteome_metadata, "Penaeus_vannamei")
get_missing_amps(moth_proteome_metadata, "Bombyx_mori")
get_missing_amps(cress_proteome_metadata, "Arabidopsis_thaliana")
get_missing_amps(frog_proteome_metadata, "Lithobates_catesbeianus")
get_missing_amps(bacteria_proteome_metadata, "Escherichia_coli")
```

**Table 6.3:** AMPs found in proteomes via the matching of AMP sequences
in the AMP database to the proteomes and remaining AMPs (from the AMP
database) that were not found in the proteomes, despite belonging to the
same organism

| Organism Name              | Number of AMPs in proteome not annotated as AMP | AMPs not accounted for |
|----------------------------|-------------------------------------------------|------------------------|
| *Mus musculus*             | 1                                               | 4                      |
| *Homo sapiens*             | 1                                               | 0                      |
| *Bos taurus*               | 1                                               | 3                      |
| *Oryctolagus cuniculus*    | 0                                               | 0                      |
| *Ornithorhynchus anatinus* | 0                                               | 0                      |
| *Gallus gallus*            | 1                                               | 0                      |
| *Oncorhynchus mykiss*      | 2                                               | 10                     |
| *Drosophila melanogaster*  | 3                                               | 0                      |
| *Penaeus vannamei*         | 0                                               | 18                     |
| *Bombyx mori*              | 0                                               | 2                      |
| *Arabidopsis thaliana*     | 2                                               | 1                      |
| *Lithobates catesbeianus*  | 1                                               | 13                     |
| *Escherichia coli K-12*    | 0                                               | 25                     |

AMPs not accounted for was calculated by subtracting the
AMPs_overlap_in_AMP_db from AMPs_in_AMP_db (from the previous table 6.2)
and then substracting this number from the Number of AMPs in proteome
not annotated as AMP

*E. coli* issue

The current proteome used for *E. coli* is UP000000318 which refers to
*E. coli* with the organism ID 83333. However, this organism only
contains four AMPs in the AMP database. The organism with organism ID
562 includes the most AMPs (25) in the AMP database. Looking at the
[available proteomes for *E.
coli*](https://www.uniprot.org/proteomes/?query=taxonomy%3A%22Escherichia+coli+%5B562%5D%22&sort=score),
it can be seen that the 562 organism does not have a reference proteome.

``` r
uniprot_and_amp_dbs_amps %>% filter(Organism == "Escherichia_coli") %>% select(Entry, Proteomes, `Organism ID`)
```

    ## # A tibble: 29 √ó 3
    ##    Entry  Proteomes                                        `Organism ID`
    ##    <chr>  <chr>                                                    <dbl>
    ##  1 P77551 UP000000318: Chromosome, UP000000625: Chromosome         83333
    ##  2 P75719 UP000000318: Chromosome, UP000000625: Chromosome         83333
    ##  3 Q46971 <NA>                                                       562
    ##  4 P06716 <NA>                                                       562
    ##  5 P04479 <NA>                                                       562
    ##  6 Q47500 <NA>                                                       562
    ##  7 P17998 <NA>                                                       562
    ##  8 P18000 <NA>                                                       562
    ##  9 P17999 <NA>                                                       562
    ## 10 Q47125 <NA>                                                       562
    ## # ‚Ä¶ with 19 more rows

Change the AMPs labelled as non-AMPs to AMPs

``` r
mouse_proteome_metadata <- mouse_proteome_metadata %>% mutate(Label = ifelse(Entry == "Q5M9M1", "Pos", Label))
human_proteome_metadata <- human_proteome_metadata %>% mutate(Label = ifelse(Entry == "J3KNB4", "Pos", Label)) 
cow_proteome_metadata <- cow_proteome_metadata %>% mutate(Label = ifelse(Entry == "A0A452DIN2", "Pos", Label)) 
junglefowl_proteome_metadata <- junglefowl_proteome_metadata %>% mutate(Label = ifelse(Entry == "A0A3Q2TTA1", "Pos", Label))
trout_proteome_metadata <- trout_proteome_metadata %>% mutate(Label = case_when(Entry == "C1BHI8" ~ "Pos",
                                                                                Entry == "C1BFG2" ~ "Pos",
                                                                                TRUE ~ Label))
fruitfly_proteome_metadata <- fruitfly_proteome_metadata %>% mutate(Label = case_when(Entry == "X2J8E8" ~ "Pos",
                                                                                Entry == "A0A0B4K7Q5" ~ "Pos",
                                                                                Entry == "D3PFG8" ~ "Pos",
                                                                                TRUE ~ Label))
cress_proteome_metadata <- cress_proteome_metadata %>% mutate(Label = case_when(Entry == "A0A1P8B0X4" ~ "Pos",
                                                                                Entry == "A0A1P8ARX5" ~ "Pos",
                                                                                TRUE ~ Label))
frog_proteome_metadata <- frog_proteome_metadata %>% mutate(Label = ifelse(Entry == "A0A2G9Q9C7", "Pos", Label)) 
```

### Creating a more stringent test set

Table 6.2 shows that there are generally more AMPs in the proteomes for
each organism compared to the AMPs in the AMP database. To keep the
highest quality data (with more likely experimentally verified AMPs),
proteomes were selected that had more than 10 AMPs overlap with the AMP
database and proteins in these proteomes were only annotated as an AMP
if the proteins were also present in the AMP database. *O. mykiss*, *P.
vannamei*, *L. catesbeianus* and *E. coli* all had less than 10 AMPs
with the AMP database and were therefore excluded.

For the remaining nine organisms, a new ‚Äúpositive‚ÄĚ label was added to
the sequences, if these were present in the AMP database. The same
analysis were performed on these nine organisms as to the full 13
organisms that contain less stringent AMP labels.

``` r
mouse_proteome_metadata <- mouse_proteome_metadata %>% mutate(Label_strict = ifelse(Entry %in% uniprot_and_amp_dbs_amps$Entry, "Pos", "Neg"))
human_proteome_metadata <- human_proteome_metadata %>% mutate(Label_strict = ifelse(Entry %in% uniprot_and_amp_dbs_amps$Entry, "Pos", "Neg"))
cow_proteome_metadata <- cow_proteome_metadata %>% mutate(Label_strict = ifelse(Entry %in% uniprot_and_amp_dbs_amps$Entry, "Pos", "Neg"))
rabbit_proteome_metadata <- rabbit_proteome_metadata %>% mutate(Label_strict = ifelse(Entry %in% uniprot_and_amp_dbs_amps$Entry, "Pos", "Neg"))
platypus_proteome_metadata <- platypus_proteome_metadata %>% mutate(Label_strict = ifelse(Entry %in% uniprot_and_amp_dbs_amps$Entry, "Pos", "Neg"))
junglefowl_proteome_metadata <-junglefowl_proteome_metadata %>% mutate(Label_strict = ifelse(Entry %in% uniprot_and_amp_dbs_amps$Entry, "Pos", "Neg"))
fruitfly_proteome_metadata <- fruitfly_proteome_metadata %>% mutate(Label_strict = ifelse(Entry %in% uniprot_and_amp_dbs_amps$Entry, "Pos", "Neg"))
moth_proteome_metadata <- moth_proteome_metadata %>% mutate(Label_strict = ifelse(Entry %in% uniprot_and_amp_dbs_amps$Entry, "Pos", "Neg"))
cress_proteome_metadata <- cress_proteome_metadata %>% mutate(Label_strict = ifelse(Entry %in% uniprot_and_amp_dbs_amps$Entry, "Pos", "Neg"))

# add NA values column so datasets maintain same columns for combining them later
trout_proteome_metadata <- trout_proteome_metadata %>% mutate(Label_strict = NA)
shrimp_proteome_metadata <- shrimp_proteome_metadata %>% mutate(Label_strict = NA)
frog_proteome_metadata <- frog_proteome_metadata %>% mutate(Label_strict = NA)
bacteria_proteome_metadata <- bacteria_proteome_metadata %>% mutate(Label_strict = NA)
```

## BLAST results

The [BLAST+](https://pubmed.ncbi.nlm.nih.gov/20003500/) version used was
blast 2.11.0, build Nov 17 2020 for MacOS.

Each proteome was used to make a local BLAST database using
`makeblastdb`. This proteome database was then used to query the AMP
dataset with `blastp`, Protein-Protein BLAST 2.11.0+. The
[BLAST](BLAST.sh) script was used for this.

``` r
parse_blast_results <- function(blast_results_path, metadata) {
  
  blast_colnames <- c("qaccver","saccver","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore")
  
  read_tsv(blast_results_path, col_names = blast_colnames, show_col_types = FALSE) %>% 
  group_by(saccver) %>% 
  slice_max(n = 1, order_by = bitscore) %>%
  separate(saccver, into = c(NA, NA, "Entry_name"), sep = "\\|") %>%
  right_join(metadata, by = "Entry_name") %>% 
  mutate(bitscore = replace_na(bitscore, 0)) %>%
  distinct(Entry_name, .keep_all = TRUE)
}

mouse_blast <- parse_blast_results("data/blastp_results/Mus_musculus.blastp", mouse_proteome_metadata)
human_blast <- parse_blast_results("data/blastp_results/Homo_sapiens.blastp", human_proteome_metadata)
cow_blast <- parse_blast_results("data/blastp_results/Bos_taurus.blastp", cow_proteome_metadata)
rabbit_blast <- parse_blast_results("data/blastp_results/Oryctolagus_cuniculus.blastp", rabbit_proteome_metadata)
platypus_blast <- parse_blast_results("data/blastp_results/Ornithorhynchus_anatinus.blastp", platypus_proteome_metadata)
junglefowl_blast <- parse_blast_results("data/blastp_results/Gallus_gallus.blastp", junglefowl_proteome_metadata)
trout_blast <- parse_blast_results("data/blastp_results/Oncorhynchus_mykiss.blastp", trout_proteome_metadata)
fruitfly_blast <- parse_blast_results("data/blastp_results/Drosophila_melanogaster.blastp", fruitfly_proteome_metadata)
shrimp_blast <- parse_blast_results("data/blastp_results/Penaeus_vannamei.blastp", shrimp_proteome_metadata)
moth_blast <- parse_blast_results("data/blastp_results/Bombyx_mori.blastp", moth_proteome_metadata)
cress_blast <- parse_blast_results("data/blastp_results/Arabidopsis_thaliana.blastp", cress_proteome_metadata)
frog_blast <- parse_blast_results("data/blastp_results/Lithobates_catesbeianus.blastp", frog_proteome_metadata)
bacteria_blast <- parse_blast_results("data/blastp_results/Escherichia_coli.blastp", bacteria_proteome_metadata)

blast_results <- rbind(mouse_blast, human_blast, cow_blast, rabbit_blast, platypus_blast, junglefowl_blast, trout_blast, fruitfly_blast, shrimp_blast, moth_blast, cress_blast, frog_blast, bacteria_blast)
```

## Prediction with AMP classification model

Read in models

``` r
mouse_model <- readRDS("models/mouse_model.rds")
human_model <- readRDS("models/human_model.rds")
cow_model <- readRDS("models/cow_model.rds")
rabbit_model <- readRDS("models/rabbit_model.rds")
platypus_model <- readRDS("models/platypus_model.rds")
junglefowl_model <- readRDS("models/junglefowl_model.rds")
trout_model <- readRDS("models/trout_model.rds")
fruitfly_model <- readRDS("models/fruitfly_model.rds")
shrimp_model <- readRDS("models/shrimp_model.rds")
moth_model <- readRDS("models/moth_model.rds")
cress_model <- readRDS("models/cress_model.rds")
frog_model <- readRDS("models/frog_model.rds")
bacteria_model <- readRDS("models/bacteria_model.rds")
```

Predict AMPs in proteomes

``` r
mouse_pred <- read_faa("data/proteomes/M_musculus_UP000000589_10090.fasta.gz") %>% predict_amps(n_cores = 3, model = mouse_model)
human_pred <- read_faa("data/proteomes/H_sapiens_UP000005640_9606.fasta.gz") %>% predict_amps(n_cores = 4, model = human_model)
cow_pred <- read_faa("data/proteomes/B_taurus_UP000009136_9913.fasta.gz") %>% predict_amps(n_cores = 4, model = cow_model)
rabbit_pred <- read_faa("data/proteomes/O_cuniculus_UP000001811_9986.fasta.gz") %>% predict_amps(n_cores = 4, model = rabbit_model)
platypus_pred <- read_faa("data/proteomes/O_anatinus_UP000002279_9258.fasta.gz") %>% predict_amps(n_cores = 4, model = platypus_model)
junglefowl_pred <- read_faa("data/proteomes/G_gallus_UP000000539_9031.fasta.gz") %>% predict_amps(n_cores = 4, model = junglefowl_model)
trout_pred <- read_faa("data/proteomes/O_mykiss_UP000193380_8022.fasta.gz") %>% predict_amps(n_cores = 4, model = trout_model)
fruitfly_pred <- read_faa("data/proteomes/D_melanogaster_UP000000803_7227.fasta.gz") %>% predict_amps(n_cores = 4, model = fruitfly_model)
shrimp_pred <- read_faa("data/proteomes/P_vannamei_UP000283509_6689.fasta.gz") %>% predict_amps(n_cores = 4, model = shrimp_model)
moth_pred <- read_faa("data/proteomes/B_mori_UP000005204_7091.fasta.gz") %>% predict_amps(n_cores = 4, model = moth_model)
cress_pred <- read_faa("data/proteomes/A_thaliana_UP000006548_3702.fasta.gz") %>% predict_amps(n_cores = 4, model = cress_model)
frog_pred <- read_faa("data/proteomes/L_catesbeianus_UP000228934_8400.fasta.gz") %>% predict_amps(n_cores = 4, model = frog_model)
bacteria_pred <- read_faa("data/proteomes/E_coli_k12_UP000000625_83333.fasta.gz") %>% predict_amps(n_cores = 4, model = bacteria_model)
```

Combine the predictions with metadata

``` r
join_pred_with_metadata <- function(pred_data, metadata){
  pred_data %>%
  mutate(Entry_name = str_extract(seq_name, "(?<=\\|)[a-zA-Z0-9_]*(?=\\s)")) %>% 
  select(Entry_name, seq_aa, prob_AMP) %>% 
  right_join(metadata, by = "Entry_name") %>%
  mutate(prob_AMP = replace_na(prob_AMP, 0))
}
```

``` r
mouse_proteome_pred <- join_pred_with_metadata(mouse_pred, mouse_proteome_metadata)
human_proteome_pred <- join_pred_with_metadata(human_pred, human_proteome_metadata)
cow_proteome_pred <- join_pred_with_metadata(cow_pred, cow_proteome_metadata)
rabbit_proteome_pred <- join_pred_with_metadata(rabbit_pred, rabbit_proteome_metadata)
platypus_proteome_pred <- join_pred_with_metadata(platypus_pred, platypus_proteome_metadata)
junglefowl_proteome_pred <- join_pred_with_metadata(junglefowl_pred, junglefowl_proteome_metadata)
trout_proteome_pred <- join_pred_with_metadata(trout_pred, trout_proteome_metadata)
fruitfly_proteome_pred <- join_pred_with_metadata(fruitfly_pred, fruitfly_proteome_metadata)
shrimp_proteome_pred <- join_pred_with_metadata(shrimp_pred, shrimp_proteome_metadata)
moth_proteome_pred <- join_pred_with_metadata(moth_pred, moth_proteome_metadata)
cress_proteome_pred <- join_pred_with_metadata(cress_pred, cress_proteome_metadata)
frog_proteome_pred <- join_pred_with_metadata(frog_pred, frog_proteome_metadata)
bacteria_proteome_pred <- join_pred_with_metadata(bacteria_pred, bacteria_proteome_metadata)

proteome_predictions <- rbind(mouse_proteome_pred, human_proteome_pred, cow_proteome_pred, rabbit_proteome_pred, platypus_proteome_pred, junglefowl_proteome_pred, trout_proteome_pred, fruitfly_proteome_pred, shrimp_proteome_pred, moth_proteome_pred, cress_proteome_pred, frog_proteome_pred, bacteria_proteome_pred)  
```

## Correctly identified AMPs

``` r
blastAMPsidentified <- blast_results %>% 
  filter(bitscore >= 50 & Label == "Pos") %>%
  count(Organism, name = "AMPs_found_BLAST") 

classificationAMPsidentified <- proteome_predictions %>%
  filter(prob_AMP >= 0.5 & Label == "Pos") %>%
  count(Organism, name = "AMPs_found_Classification")

total_amp_count <- proteome_predictions %>% 
  filter(Label == "Pos") %>% 
  count(Organism, name = "Total_AMP_count")

amps_identified <- blastAMPsidentified %>%
  left_join(classificationAMPsidentified, by = "Organism") %>%
  left_join(total_amp_count, by = "Organism")
```

**Table 6.4:** Correctly identified AMPs in different proteomes with the
BLAST and classification methods.

| Organism                 | AMPs_found_BLAST | AMPs_found_Classification | Total_AMP_count |
|:-------------------------|-----------------:|--------------------------:|----------------:|
| Arabidopsis_thaliana     |               23 |                       172 |             296 |
| Bombyx_mori              |               19 |                        11 |              25 |
| Bos_taurus               |               85 |                        78 |             117 |
| Drosophila_melanogaster  |               23 |                        16 |              33 |
| Escherichia_coli         |                2 |                         2 |               4 |
| Gallus_gallus            |               14 |                        21 |              30 |
| Homo_sapiens             |               77 |                        62 |             116 |
| Lithobates_catesbeianus  |               11 |                         8 |              12 |
| Mus_musculus             |               96 |                        70 |             132 |
| Oncorhynchus_mykiss      |               13 |                        11 |              17 |
| Ornithorhynchus_anatinus |               11 |                        12 |              27 |
| Oryctolagus_cuniculus    |               54 |                        39 |              83 |
| Penaeus_vannamei         |                3 |                         1 |               3 |

None of the 11 AMPs in platypus were identified by BLAST in the stricter
AMP criteria. `count` does not count zeros in the current version so the
BLAST results for platypus was manually added as a row

``` r
organism_strict_selection <- c("Mus_musculus","Homo_sapiens","Bos_taurus","Oryctolagus_cuniculus","Ornithorhynchus_anatinus","Gallus_gallus", "Drosophila_melanogaster","Bombyx_mori","Arabidopsis_thaliana")

platypus_amp_count_found_by_blast <- blast_results %>% filter(Organism == "Ornithorhynchus_anatinus" & bitscore >= 50 & Label_strict == "Pos") %>% nrow() 

blastAMPsidentified2 <- blast_results %>% 
  filter(Organism %in% organism_strict_selection) %>%
  filter(bitscore >= 50 & Label_strict== "Pos") %>%
  count(Organism, name = "AMPs_found_BLAST") %>%
  add_row(Organism = "Ornithorhynchus_anatinus", AMPs_found_BLAST = platypus_amp_count_found_by_blast)

classificationAMPsidentified2 <- proteome_predictions %>%
  filter(prob_AMP >= 0.5 & Label_strict == "Pos") %>%
  count(Organism, name = "AMPs_found_Classification")

total_amp_count2 <- proteome_predictions %>% 
  filter(Organism %in% organism_strict_selection & Label_strict == "Pos") %>%
  count(Organism, name = "Total_AMPs_present")

amps_identified2 <- blastAMPsidentified2 %>%
  left_join(classificationAMPsidentified2, by = "Organism") %>%
  left_join(total_amp_count2, by = "Organism")
```

**Table 6.5:** Same as Table 6.4 but using the stricter AMP criteria,
i.e.¬†where a protein is considered to be an AMP if it overlaps with the
AMP database

| Organism                 | AMPs_found_BLAST | AMPs_found_Classification | Total_AMPs_present |
|:-------------------------|-----------------:|--------------------------:|-------------------:|
| Arabidopsis_thaliana     |               23 |                       172 |                291 |
| Bombyx_mori              |               10 |                         9 |                 13 |
| Bos_taurus               |               39 |                        37 |                 54 |
| Drosophila_melanogaster  |               17 |                        15 |                 23 |
| Gallus_gallus            |               11 |                        17 |                 25 |
| Homo_sapiens             |               73 |                        58 |                 95 |
| Mus_musculus             |               82 |                        61 |                 99 |
| Oryctolagus_cuniculus    |               12 |                        10 |                 17 |
| Ornithorhynchus_anatinus |                0 |                         7 |                 11 |

## Calculate the Precision Recall (PR) curve and the Area Under the Precision Recall Curve (AUPRC) for each method

### PR curve

**For all 13 organisms using the ‚ÄúAntimicrobial‚ÄĚ keyword as positive
AMP**

``` r
organisms <- unique(proteome_predictions$Organism)

source("scripts/calc_cm_metrics_from_bitscore.R")

get_blast_roc <- function(data, method){
  do.call(rbind,lapply(organisms,function(org){ 
    map_df(seq(0, 1000, 1), calc_cm_metrics_from_bitscore, data %>% filter(Organism==org)) %>%
    add_column(Organism = org)
  })) %>%   
  add_column(Method = method)
}

source("scripts/calc_cm_metrics_from_prob.R")

get_proteome_roc <- function(data, method){
  do.call(rbind,lapply(organisms,function(org){ 
    map_df(c(seq(0.01, 0.99, 0.01),seq(0.99, 0.990, 0.001)), calc_cm_metrics_from_prob, data %>% filter(Organism==org)) %>%
    add_column(Organism = org)
  })) %>%   
  add_column(Method = method)
}
```

``` r
blast_roc <- get_blast_roc(blast_results, "BLAST")

pred_roc <- get_proteome_roc(proteome_predictions, "Classification")


blast_and_pred_roc <- rbind(blast_roc, pred_roc)
```

``` r
saveRDS(blast_and_pred_roc, "cache/blast_and_pred_roc13.rds")
```

**For nine organisms using the AMPs that overlap with the AMP database
as positive AMP**

``` r
calc_pr_bitscore <- function(p_threshold, df) {
  
  TP <- df %>% filter((Label_strict =="Pos")) %>% filter(bitscore > p_threshold) %>% n_distinct()
  FP <- df %>% filter((Label_strict =="Neg")) %>% filter(bitscore > p_threshold) %>% n_distinct()
  TN <- df %>% filter((Label_strict =="Neg")) %>% filter(bitscore < p_threshold) %>% n_distinct()
  FN <- df %>% filter((Label_strict =="Pos")) %>% filter(bitscore < p_threshold) %>% n_distinct()
  
  Recall <- TP / (TP + FN) 
  Precision <- TP / (TP + FP) 

  cm <- c(Recall, Precision, p_threshold)
  names(cm) <-c("Recall", "Precision","p_threshold") 
  cm
}

calc_pr_prob <- function(p_threshold, df) {
  
  TP <- df %>% filter((Label_strict =="Pos")) %>% filter(prob_AMP > p_threshold) %>% n_distinct()
  FP <- df %>% filter((Label_strict =="Neg")) %>% filter(prob_AMP > p_threshold) %>% n_distinct()
  TN <- df %>% filter((Label_strict =="Neg")) %>% filter(prob_AMP < p_threshold) %>% n_distinct()
  FN <- df %>% filter((Label_strict =="Pos")) %>% filter(prob_AMP < p_threshold) %>% n_distinct()
  
  Recall <- TP / (TP + FN) 
  Precision <- TP / (TP + FP) 

  cm <- c(Recall, Precision, p_threshold)
  names(cm) <-c("Recall", "Precision","p_threshold") 
  cm
}

get_blast_roc_strict <- function(data, method){
  do.call(rbind,lapply(organism_strict_selection,function(org){ 
    map_df(seq(0, 1000, 1), calc_pr_bitscore, data %>% filter(Organism==org)) %>%
    add_column(Organism = org)
  })) %>%   
  add_column(Method = method)
}

get_proteome_roc_strict <- function(data, method){
  do.call(rbind,lapply(organism_strict_selection,function(org){ 
    map_df(c(seq(0.01, 0.99, 0.01),seq(0.99, 0.990, 0.001)), calc_pr_prob, data %>% filter(Organism==org)) %>%
    add_column(Organism = org)
  })) %>%   
  add_column(Method = method)
}
```

``` r
proteome_predictions_9 <- filter(proteome_predictions, Organism %in% organism_strict_selection)
blast_results_9 <- filter(blast_results, Organism %in% organism_strict_selection)

blast_roc_strict <- get_blast_roc_strict(blast_results_9, "BLAST")

pred_roc_strict <- get_proteome_roc_strict(proteome_predictions_9, "Classification")


blast_and_pred_roc_strict <- rbind(blast_roc_strict, pred_roc_strict)
```

``` r
saveRDS(blast_and_pred_roc_strict, "cache/blast_and_pred_roc_strict_9.rds")
```

### AUPRC

First made a function to calculate the AUPRC from a dataset, for the
classification and BLAST method, and then made a function to loop this
function to calculate the AUPRC for each organism in the dataset. This
step was done twice, once for the 13 organisms and once for the nine
with stricter AMP

``` r
calculate_auprc_probAMP <- function(df) {
  evalmod(scores = df[["prob_AMP"]], labels = df[["Label"]], mode = "rocprc") %>%
  precrec::auc() %>%
  select(curvetypes, aucs) %>%
  filter(curvetypes == "PRC") %>%
  pivot_wider(names_from = curvetypes, values_from = aucs) %>%
  rename(AUPRC = "PRC") %>%
  round(digits = 3)
}

organisms <- unique(proteome_predictions$Organism)

get_probAMP_auprc <- function(prediction_data, method_name){
  do.call(rbind,lapply(organisms,function(org){ 
    calculate_auprc_probAMP(prediction_data %>% filter(Organism==org)) %>%
    add_column(Organism = org)
  })) %>%   
  add_column(Method = method_name)
}

calculate_auprc_bitscore <- function(df) {
  evalmod(scores = df[["bitscore"]], labels = df[["Label"]], mode = "rocprc") %>%
  precrec::auc() %>%
  select(curvetypes, aucs) %>%
  filter(curvetypes == "PRC") %>%
  pivot_wider(names_from = curvetypes, values_from = aucs) %>%
  rename(AUPRC = "PRC") %>%
  round(digits = 3)
}

get_bitscore_auprc <- function(blast_results_data, method_name){
  do.call(rbind,lapply(organisms,function(org){ 
    calculate_auprc_bitscore(blast_results_data %>% filter(Organism==org)) %>%
    add_column(Organism = org)
  })) %>%   
  add_column(Method = method_name)
}


classification_auprc <- get_probAMP_auprc(proteome_predictions, "Classification")
blast_auprc <- get_bitscore_auprc(blast_results, "BLAST")


methods_auprc <- rbind(classification_auprc, blast_auprc)
```

``` r
proteome_predictions_9 <- filter(proteome_predictions, Organism %in% organism_strict_selection)
blast_results_9 <- filter(blast_results, Organism %in% organism_strict_selection)

organisms9 <- unique(proteome_predictions_9$Organism)

calculate_auprc_probAMP_9 <- function(df) {
  evalmod(scores = df[["prob_AMP"]], labels = df[["Label_strict"]], mode = "rocprc") %>%
  precrec::auc() %>%
  select(curvetypes, aucs) %>%
  filter(curvetypes == "PRC") %>%
  pivot_wider(names_from = curvetypes, values_from = aucs) %>%
  rename(AUPRC = "PRC") %>%
  round(digits = 3)
}

get_probAMP_auprc_9 <- function(prediction_data, method_name){
  do.call(rbind,lapply(organisms9,function(org){ 
    calculate_auprc_probAMP_9(prediction_data %>% filter(Organism==org)) %>%
    add_column(Organism = org)
  })) %>%   
  add_column(Method = method_name)
}

calculate_auprc_bitscore_9 <- function(df) {
  evalmod(scores = df[["bitscore"]], labels = df[["Label_strict"]], mode = "rocprc") %>%
  precrec::auc() %>%
  select(curvetypes, aucs) %>%
  filter(curvetypes == "PRC") %>%
  pivot_wider(names_from = curvetypes, values_from = aucs) %>%
  rename(AUPRC = "PRC") %>%
  round(digits = 3)
}

get_bitscore_auprc_9 <- function(blast_results_data, method_name){
  do.call(rbind,lapply(organisms9,function(org){ 
    calculate_auprc_bitscore_9(blast_results_data %>% filter(Organism==org)) %>%
    add_column(Organism = org)
  })) %>%   
  add_column(Method = method_name)
}

classification_auprc_9 <- get_probAMP_auprc_9(proteome_predictions_9, "Classification")
blast_auprc_9 <- get_bitscore_auprc_9(blast_results_9, "BLAST")


methods_auprc_9 <- rbind(classification_auprc_9, blast_auprc_9)
```

Change methods_auprc to wide format and add their respective AMP counts
used in performance evaluation to save for use in the next workflow. For
methods_auprc_13, the AMP count in the proteomes (found via the
‚ÄúAntimicrobial‚ÄĚ keyword) was used. For methods_auprc_9, the AMP count in
the proteomes that overlapped with the AMPs in the AMP database was used
(see Table 6.2)

``` r
methods_auprc_13_wide <- pivot_wider(methods_auprc, names_from = Method, values_from = AUPRC)

total_amps_in_proteomes_antimicrobial_kw_13 <- proteome_predictions %>% 
  filter(Label == "Pos") %>%
  count(Organism, name = "Total_AMPs_in_test") 

methods_auprc_13_wide_w_count <- methods_auprc_13_wide %>% left_join(total_amps_in_proteomes_antimicrobial_kw_13, by = "Organism")

methods_auprc_13_wide_w_count
```

    ## # A tibble: 13 √ó 4
    ##    Organism                 Classification BLAST Total_AMPs_in_test
    ##    <chr>                             <dbl> <dbl>              <int>
    ##  1 Mus_musculus                      0.365 0.305                132
    ##  2 Homo_sapiens                      0.286 0.391                116
    ##  3 Bos_taurus                        0.37  0.299                117
    ##  4 Oryctolagus_cuniculus             0.22  0.205                 83
    ##  5 Ornithorhynchus_anatinus          0.157 0.097                 27
    ##  6 Gallus_gallus                     0.435 0.121                 30
    ##  7 Oncorhynchus_mykiss               0.071 0.108                 17
    ##  8 Drosophila_melanogaster           0.053 0.193                 33
    ##  9 Penaeus_vannamei                  0.014 0.064                  3
    ## 10 Bombyx_mori                       0.065 0.14                  25
    ## 11 Arabidopsis_thaliana              0.318 0.041                296
    ## 12 Lithobates_catesbeianus           0.022 0.214                 12
    ## 13 Escherichia_coli                  0.227 0.5                    4

``` r
methods_auprc_9_wide <- pivot_wider(methods_auprc_9, names_from = Method, values_from = AUPRC)

total_amps_in_proteome_overlap_w_amp_db_9 <- proteome_predictions_9 %>%
  filter(Label_strict == "Pos") %>%
  count(Organism, name = "Total_AMPs_in_test") 

methods_auprc_9_wide_w_count <- methods_auprc_9_wide %>% left_join(total_amps_in_proteome_overlap_w_amp_db_9, by = "Organism")

methods_auprc_9_wide_w_count
```

    ## # A tibble: 9 √ó 4
    ##   Organism                 Classification BLAST Total_AMPs_in_test
    ##   <chr>                             <dbl> <dbl>              <int>
    ## 1 Mus_musculus                      0.398 0.332                 99
    ## 2 Homo_sapiens                      0.296 0.437                 95
    ## 3 Bos_taurus                        0.142 0.211                 54
    ## 4 Oryctolagus_cuniculus             0.043 0.087                 17
    ## 5 Ornithorhynchus_anatinus          0.149 0.017                 11
    ## 6 Gallus_gallus                     0.439 0.078                 25
    ## 7 Drosophila_melanogaster           0.058 0.186                 23
    ## 8 Bombyx_mori                       0.079 0.094                 13
    ## 9 Arabidopsis_thaliana              0.323 0.041                291

![](02_method_evaluation_on_proteomes_files/figure-gfm/unnamed-chunk-34-1.png)<!-- -->

**Figure 6.1:** The precision-recall curve (A) and the area under the
precision-recall curve (B) for each organism and AMP finding method
(BLAST and classification models)

![](02_method_evaluation_on_proteomes_files/figure-gfm/unnamed-chunk-37-1.png)<!-- -->

**Figure 6.2:** The precision-recall curve (A) and the area under the
precision-recall curve (B) for each organism and AMP finding method
(BLAST and classification models) for nine organisms with AMPs
identified from the AMP database
