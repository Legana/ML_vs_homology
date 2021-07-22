
## Starting over with a new AMP database

A new AMP database was constructed which included:

-   the reviewed AMPs from UniProt (accessed 07 July 2021) - 3,242 AMPs
-   the unreviewed AMPs in the UniProt database, if present in the AMP
    databases APD, DRAMP or dbAMP (accessed 09 April 2021) - 170 AMPs
-   Total AMPs: 3,412

This approach is consistent with the creation of ampir’s models and
detailed in [01_collate_databases.md in the AMP_pub
repository](https://github.com/Legana/AMP_pub/blob/master/01_collate_databases.md).

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

Removing sequences containing anything other than the standard AA
results in 3,354 sequences

``` r
amps_standardaa <- uniprot_and_amp_dbs_amps %>% select("Entry_name", "Sequence") %>% as.data.frame() %>% remove_nonstandard_aa()
```

Use `cd-hit` to remove highly similar sequences

``` bash
cd-hit -i cache/amps_standardaa.fasta -o cache/amps_standardaa.fasta90.fasta -c 0.90 -g 1
```

``` r
amps_standardaa90 <- read_faa("cache/amps_standardaa.fasta90.fasta") %>% 
  left_join(uniprot_and_amp_dbs_amps, by = c("seq_name" = "Entry_name"))
```

`cd-hit` resulted in 2,339 representative sequences. 1635 of these
sequences are between 50 and 500 AA long, which is the length filter
that will be applied to train the classification models.

The phyla Chordata, Arthropoda and Streptophyta contain the most AMPs.
Therefore I want to pick the organism with the most AMPs in each
taxonomic group within these three phyla.

The phyla Evosea and Nematoda both have a single organism (*D.
discoideum* (amoeba) and *C. elegans* (roundworm), respectively) with
around 10 AMPs which could maybe be added to the test. The remaining
phyla only have organisms that have less than 5 AMPs each, so these are
not worth adding.

Notes: How many AMPs present for each *chosen* organism in new training
DB. Use range of taxonomic distances to stretch out distance x axis. 10
borderline, over 20 is good.

``` r
uniprot_and_amp_dbs_amps %>%
  mutate(in_90 = Entry_name %in% amps_standardaa90$seq_name) %>%
  filter(!grepl("Bacteria", Taxonomic_lineage)) %>% 
  filter(!grepl("Viruses", Taxonomic_lineage)) %>%
  group_by(Phylum) %>%
  summarise(AMP_count = n(), AMP_standardaa_90 = sum(in_90)) %>%
  arrange(.by_group = TRUE, desc(AMP_count))
```

    ## # A tibble: 14 x 3
    ##    Phylum                AMP_count AMP_standardaa_90
    ##    <chr>                     <int>             <int>
    ##  1 Chordata                   1747              1044
    ##  2 Arthropoda                  685               493
    ##  3 Streptophyta                556               470
    ##  4 Mollusca                     33                29
    ##  5 Evosea                       16                14
    ##  6 Nematoda (roundworms)        15                14
    ##  7 Annelida                     14                10
    ##  8 Ascomycota                   13                11
    ##  9 Cnidaria                     10                 9
    ## 10 Basidiomycota                 8                 4
    ## 11 Echinodermata                 4                 3
    ## 12 Euryarchaeota                 3                 3
    ## 13 Mucoromycota                  1                 0
    ## 14 Platyhelminthes               1                 1

## Chordata

### Organisms with most AMPs according to taxonomic **Class**

-   *Mus musculus* (mammal with 104 AMPs)
-   *Bombina maxima* (toad with 50 AMPs)
-   *Gallus gallus* (bird with 25 AMPs)
-   *Oncorhynchus mykiss* (fish with 12 AMPs)
-   *Styela clava* (tunicate with 11 AMPs)
-   *Crotalus durissus terrificus* (snake with 10 AMPs)

``` r
uniprot_and_amp_dbs_amps %>%
  mutate(in_90 = Entry_name %in% amps_standardaa90$seq_name) %>%
  filter(!grepl("Bacteria", Taxonomic_lineage)) %>% 
  filter(Phylum == "Chordata") %>%
  group_by(Order, Organism, Class) %>%
  summarise(AMP_count = n(), AMP_standardaa_90 = sum(in_90)) %>%
  ungroup() %>%
  group_by(Class) %>%
  slice_max(AMP_count, n = 1) %>%
  arrange(desc(AMP_count)) 
```

    ## # A tibble: 10 x 5
    ## # Groups:   Class [10]
    ##    Order         Organism            Class            AMP_count AMP_standardaa_…
    ##    <chr>         <chr>               <chr>                <int>            <int>
    ##  1 Rodentia      Mus_musculus        Mammalia               104               80
    ##  2 Anura         Bombina_maxima      Amphibia                50               11
    ##  3 Galliformes   Gallus_gallus       Aves                    25               20
    ##  4 Salmoniformes Oncorhynchus_mykiss Actinopteri             12                8
    ##  5 Stolidobranc… Styela_clava        Ascidiacea              11                5
    ##  6 Squamata      Crotalus_durissus_… Lepidosauria (l…        10                3
    ##  7 Testudines    Pelodiscus_sinensis <NA>                     2                2
    ##  8 Myliobatifor… Potamotrygon_cf._h… Chondrichthyes           1                1
    ##  9 Petromyzonti… Petromyzon_marinus  Hyperoartia              1                1
    ## 10 Amphioxiform… Branchiostoma_belc… Leptocardii              1                1

If taking the top organisms per **Order**, where each organism has at
least 10 AMPs, the following organisms would be candidates:

-   *Mus musculus* (mouse, mammal, 104 AMPs)
-   *Homo sapiens* (human, mammal, 96 AMPs)
-   *Bos taurus* (cattle, mammal, 58 AMPs)
-   *Oryctolagus cuniculus* (rabbit, mammal, 17 AMPs)
-   *Ornithorhynchus anatinus* (platypus, mammal, 11 AMPs )
-   *Bombina maxima* (toad, amphibian, 50 AMPs)
-   *Gallus gallus* (junglefowl, bird, 25 AMPs)
-   *Oncorhynchus mykiss* (fish with 12 AMPs)
-   *Styela clava* (tunicate with 11 AMPs)
-   *Crotalus durissus terrificus* (snake with 10 AMPs)

However, the toad *Bombina maxima*, does not have an available reference
proteome. The only [reference proteomes available for
anurans](https://www.uniprot.org/proteomes/?query=taxonomy%3A%22Anura+%289ANUR%29+%5B8342%5D%22&sort=score)
are: *Xenopus tropicalis* (*Silurana tropicalis*), *Xenopus leavis* and
*Lithobates catesbeianus* (*Rana catesbeiana*) which have 9, 10 and 13
AMPs, respectively. *Lithobates catesbeianus* has the most AMPs of these
three frogs but also has the most missing values (70%) in its reference
proteome BUSCO score.

The tunicate, *Styela clava* and the snake *Crotalus durissus
terrificus* do not have a proteome either and the other organisms within
their respective class have less than four known AMPs.

The majority of this list are mammals. The platypus however, is a
monotreme, and not a typical placental (eutherian) mammal, like the
mouse, human, cattle and rabbit are. Monotremes diverged earlier from
the mammalian common ancestor compared to the placental mammals and
subsequently have different characteristics [(Deakin, Graves and Rens
2012)](https://www.karger.com/Article/Fulltext/339433)

``` r
uniprot_and_amp_dbs_amps %>%
  mutate(in_90 = Entry_name %in% amps_standardaa90$seq_name) %>%
  filter(!grepl("Bacteria", Taxonomic_lineage)) %>% 
  filter(Phylum == "Chordata") %>%
  group_by(Order, Organism, Class) %>%
  summarise(AMP_count = n(), AMP_standardaa_90 = sum(in_90)) %>%
  ungroup() %>%
  group_by(Order) %>%
  slice_max(AMP_count, n = 1) %>%
  arrange(desc(AMP_count)) %>%
  head(10)
```

    ## # A tibble: 10 x 5
    ## # Groups:   Order [10]
    ##    Order       Organism             Class             AMP_count AMP_standardaa_…
    ##    <chr>       <chr>                <chr>                 <int>            <int>
    ##  1 Rodentia    Mus_musculus         Mammalia                104               80
    ##  2 Primates    Homo_sapiens         Mammalia                 96               64
    ##  3 Artiodacty… Bos_taurus           Mammalia                 58               41
    ##  4 Anura       Bombina_maxima       Amphibia                 50               11
    ##  5 Galliformes Gallus_gallus        Aves                     25               20
    ##  6 Lagomorpha  Oryctolagus_cunicul… Mammalia                 17               12
    ##  7 Salmonifor… Oncorhynchus_mykiss  Actinopteri              12                8
    ##  8 Monotremata Ornithorhynchus_ana… Mammalia                 11               11
    ##  9 Stolidobra… Styela_clava         Ascidiacea               11                5
    ## 10 Squamata    Crotalus_durissus_t… Lepidosauria (le…        10                3

``` r
frog_amps <- uniprot_and_amp_dbs_amps %>%
  filter(Order == "Anura") %>%
  count(Order, Organism, sort = TRUE, name = "AMP_count")
```

## Arthropoda

If choosing the organism with the most AMPs per **Order**, the following
organisms would be chosen:

-   *Lachesana tarabaevi* (ant spider, arachnid, 28 AMPs)
-   *Drosophila melanogaster* (fruitfly, insect, 23 AMPs)
-   *Penaeus vannamei* (shrimp, decapod, 18 AMPs)
-   *Neoponera goeldii* (ant, insect, 15 AMPs)
-   *Bombyx mori* (moth, insect, 15 AMPs)
-   *Tachypleus tridentatus* (horseshoe crab, 11 AMPs)
-   *Chaerilus_tricostatus* (scorpion, arachnid, 10 AMPs)

*Neoponera goeldii*, *Tachypleus tridentatus* and *Chaerilus
tricostatus* do not have a reference proteome and therefore could not be
included.

``` r
uniprot_and_amp_dbs_amps %>%
  mutate(in_90 = Entry_name %in% amps_standardaa90$seq_name) %>%
  filter(!grepl("Bacteria", Taxonomic_lineage)) %>% 
  filter(Phylum == "Arthropoda") %>%
  group_by(Order, Organism, Class) %>%
  summarise(AMP_count = n(), AMP_standardaa_90 = sum(in_90)) %>%
  ungroup() %>%
  group_by(Order) %>%
  slice_max(AMP_count, n = 1) %>%
  arrange(desc(AMP_count)) %>%
  head(10)
```

    ## # A tibble: 10 x 5
    ## # Groups:   Order [9]
    ##    Order     Organism           Class                 AMP_count AMP_standardaa_…
    ##    <chr>     <chr>              <chr>                     <int>            <int>
    ##  1 Araneae   Lachesana_tarabae… Arachnida                    28               11
    ##  2 Diptera   Drosophila_melano… Insecta                      23               17
    ##  3 Decapoda  Penaeus_vannamei   Malacostraca                 18                4
    ##  4 Hymenopt… Neoponera_goeldii  Insecta                      15               11
    ##  5 Lepidopt… Bombyx_mori        Insecta                      13                9
    ##  6 Xiphosura Tachypleus_triden… Merostomata (horsesh…        11                8
    ##  7 Scorpion… Chaerilus_tricost… Arachnida                    10                8
    ##  8 Hemiptera Palomena_prasina   Insecta                       5                4
    ##  9 Coleopte… Acrocinus_longima… Insecta                       3                2
    ## 10 Coleopte… Holotrichia_diomp… Insecta                       3                3

## Streptophyta

For plants only *Arabidopsis thaliana* would be chosen which has 284
AMPs (after `cd-hit`)

``` r
uniprot_and_amp_dbs_amps %>%
  mutate(in_90 = Entry_name %in% amps_standardaa90$seq_name) %>%
  filter(!grepl("Bacteria", Taxonomic_lineage)) %>% 
  filter(Phylum == "Streptophyta") %>%
  group_by(Order, Organism, Class) %>%
  summarise(AMP_count = n(), AMP_standardaa_90 = sum(in_90)) %>%
  ungroup() %>%
  group_by(Order) %>%
  slice_max(AMP_count, n = 1) %>%
  arrange(desc(AMP_count)) %>%
  head(10)
```

    ## # A tibble: 10 x 5
    ## # Groups:   Order [10]
    ##    Order          Organism               Class        AMP_count AMP_standardaa_…
    ##    <chr>          <chr>                  <chr>            <int>            <int>
    ##  1 Brassicales    Arabidopsis_thaliana   Magnoliopsi…       294              284
    ##  2 Poales         Triticum_kiharae       Magnoliopsi…         9                6
    ##  3 Caryophyllales Spinacia_oleracea      Magnoliopsi…         7                3
    ##  4 Fabales        Clitoria_ternatea      Magnoliopsi…         7                6
    ##  5 Malvales       Malva_parviflora       Magnoliopsi…         7                3
    ##  6 Asterales      Taraxacum_officinale   Magnoliopsi…         5                5
    ##  7 Ranunculales   Nigella_sativa         Magnoliopsi…         5                3
    ##  8 Solanales      Solanum_tuberosum      Magnoliopsi…         5                5
    ##  9 Proteales      Macadamia_integrifolia Magnoliopsi…         4                2
    ## 10 Arecales       Cocos_nucifera         Magnoliopsi…         3                1

## Bacteria

Bacteria contain AMPs called bacteriocins used for competitive
advantages. *Escherichia coli* contains the most known AMPs (29) and
have two reference proteomes corresponding to two different strains: K12
and 0157:H7 where the latter is considered to be the pathogenic strain,
which is inhibited by AMPs produced from the non-pathogenic *E. coli*
strains [(Askari and Ghanbarpour
2019)](https://bmcvetres.biomedcentral.com/articles/10.1186/s12917-018-1771-y)

``` r
bacteriocins <- uniprot_and_amp_dbs_amps %>%
 filter(grepl("Bacteria", Taxonomic_lineage)) %>%
 count(Organism, sort = TRUE, name = "AMP_count")
```

## Final organism selection

Table 5.1: Final organism selection

| Organism Name                         | Reference proteome ID                                        | Total proteins | AMPs | Gene count |
|---------------------------------------|--------------------------------------------------------------|----------------|------|------------|
| *Mus musculus* (mouse)                | [UP000000589](https://www.uniprot.org/proteomes/UP000000589) | 55,366         | 104  | 22,001     |
| *Homo sapiens* (human)                | [UP000005640](https://www.uniprot.org/proteomes/UP000005640) | 78,120         | 96   | 20,600     |
| *Bos taurus* (cattle)                 | [UP000009136](https://www.uniprot.org/proteomes/UP000009136) | 37,513         | 58   | 23,847     |
| *Oryctolagus cuniculus* (rabbit)      | [UP000001811](https://www.uniprot.org/proteomes/UP000001811) | 41,459         | 17   | 21,193     |
| *Ornithorhynchus anatinus* (platypus) | [UP000002279](https://www.uniprot.org/proteomes/UP000002279) | 32,824         | 11   | 17,390     |
| *Gallus gallus* (bird)                | [UP000000539](https://www.uniprot.org/proteomes/UP000000539) | 27,535         | 25   | 18,113     |
| *Oncorhynchus mykiss* (fish)          | [UP000193380](https://www.uniprot.org/proteomes/UP000193380) | 46,447         | 12   | 46,405     |
| *Drosophila melanogaster* (fruitfly)  | [UP000000803](https://www.uniprot.org/proteomes/UP000000803) | 22,110         | 23   | 13,821     |
| *Penaeus vannamei* (shrimp)           | [UP000283509](https://www.uniprot.org/proteomes/UP000283509) | 25,399         | 18   | 25,399     |
| *Bombyx mori* (moth)                  | [UP000005204](https://www.uniprot.org/proteomes/UP000005204) | 14,776         | 15   | 14,773     |
| *Arabidopsis thaliana* (plant)        | [UP000006548](https://www.uniprot.org/proteomes/UP000006548) | 39,337         | 294  | 27,468     |
| *Lithobates catesbeianus* (frog)      | [UP000228934](https://www.uniprot.org/proteomes/UP000228934) | 28,218         | 13   | 28,218     |
| *Escherichia coli K-12* (bacteria)    | [UP000228934](https://www.uniprot.org/proteomes/UP000228934) | 4,438          | 29   | 4,392      |
