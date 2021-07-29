
# BLAST vs. classification models for finding AMPs in proteomes

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

Table 6.1: Proteomes used with protein and gene count obtained from
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
  read_tsv(path, col_types = cols()) %>%
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

amp_overlap %>% select(Organism, AMPs_in_proteome, AMPs_overlap_in_AMP_db, AMPs_in_AMP_db)
```

    ## # A tibble: 13 x 4
    ##    Organism                AMPs_in_proteome AMPs_overlap_in_AMP_… AMPs_in_AMP_db
    ##    <chr>                              <int>                 <int>          <dbl>
    ##  1 Mus musculus                         131                    99            104
    ##  2 Homo sapiens                         115                    95             96
    ##  3 Bos taurus                           116                    54             58
    ##  4 Oryctolagus cuniculus                 83                    17             17
    ##  5 Ornithorhynchus anatin…               27                    11             11
    ##  6 Gallus gallus                         29                    25             25
    ##  7 Oncorhynchus mykiss                   15                     0             12
    ##  8 Drosophila melanogaster               30                    23             23
    ##  9 Penaeus vannamei                       3                     0             18
    ## 10 Bombyx mori                           25                    13             15
    ## 11 Arabidopsis thaliana                 294                   291            294
    ## 12 Lithobates catesbeianus               11                     0             13
    ## 13 Escherichia coli                       4                     4             29

From the `amp_overlap` table, it can be seen that there are AMPs in the
proteomes that do not overlap with the AMPs in the AMP database. A
possible explanation for this is that the proteome AMPs were found via
the “Antimicrobial” keyword in UniProt. These contain proteins not
experimentally verified to contain antimicrobial activity and are
therefore not included in the AMP database. Additionally, there are AMPs
in the AMP database not present in the proteomes.

To see if this is due to the presence of mature AMP sequences in the AMP
database, the AMPs in the AMP database corresponding to each organism
were matched to full length sequences, annotated as not antimicrobial,
in the organisms’ respective proteomes. However, only a very few number
of AMPs were detected with this method. The absence of these AMP in the
proteome

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
AMPs_overlap_in_AMP_db from AMPs_in_AMP_db and then substracting this
number from the Number of AMPs in proteome not annotated as AMP

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
  
  read_tsv(blast_results_path, col_names = blast_colnames) %>% 
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
