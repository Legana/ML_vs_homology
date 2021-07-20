
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
Want to pick the organism with the most AMPs in each taxonomic group
within these three phyla.

The phyla Evosea and Nematoda both have a single organism (*D.
discoideum* (amoeba) and *C. elegans* (roundworm), respectively) with
around 10 AMPs which could maybe be added to the test. The remaining
phyla only have organisms that have less than 5 AMPs each, so these are
not worth adding.

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

Initially I thought to maybe take the top organism in every Chordata
Class (rather than from every Order, as that would have been a lot of
organisms), but after using `cd-hit` which removed many similar AMPs,
the AMP count in some classes/orders drastically changed. If I were to
take the organisms with the top AMP count, going by taxonomic **Class**,
I would only end up with:

-   *Mus musculus* (mammal with 80 AMPs)
-   *Bombina maxima* (toad with 11 AMPs)
-   *Gallus gallus* (bird with 20 AMPs)

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

How does representative sequence work for `cd-hit`? would it have a
different result each time (i.e remove all the similar sequences but
choose a different organism to keep with its representative sequence).
Is it safe to pick organisms based on AMP_standardaa_90?

If choosing the organism with the most (unique?) AMPs per **Order**, the
following organisms would be chosen:

-   *Mus musculus* (mouse, mammal, 80 AMPs)
-   *Homo sapiens* (human, mammal, 64 AMPs)
-   *Bos taurus* (cattle, mammal, 41 AMPs)
-   *Bombina maxima* (toad, amphibian, 11 AMPs)
-   *Gallus gallus* (junglefowl, bird, 20 AMPs)
-   *Oryctolagus cuniculus* (rabbit, mammal, 12 AMPs)
-   *Ornithorhynchus anatinus* (platypus, mammal, 11 AMPs )

This is still mostly mammals though, and the list only has 1 fish with
only 8 AMPs. I find this a bit odd, does this mean the fish AMPs were
mostly homologous to other (mammal?) AMPs and got cut out by `cd-hit`?

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

## Arthropoda

If choosing the organism with the most AMPs per **Order**, the following
organisms would be chosen:

-   *Drosophila melanogaster* (fruitfly, insect, 17 AMPs)
-   *Lachesana tarabaevi* (ant spider, arachnid, 11 AMPs)
-   *Neoponera_goeldii* (ant, insect, 11 AMPs)

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
