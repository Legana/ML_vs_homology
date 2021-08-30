
``` r
library(tidyverse)
library(ape)
library(treeio)
library(patchwork)
library(pals)
library(ggtext)
```

Read in AMP database to extract organisms from to submit to
[TimeTree](http://timetree.org/) for phylogenetic data.

``` r
amp_database <- readRDS("data/uniprot_amps_w_amp_dbsJuly21.rds") %>%
  rename(Entry_name = `Entry name`) %>% 
  mutate(Organism = str_remove(Organism, " \\(.*")) %>% 
  rename(Taxonomic_lineage = `Taxonomic lineage (ALL)`) %>% 
  rename(Order = `Taxonomic lineage (ORDER)`) %>% 
  rename(Phylum = `Taxonomic lineage (PHYLUM)`) %>%
  rename(Class = `Taxonomic lineage (CLASS)`) %>%
  mutate(Order = str_remove(Order, " \\(.*")) %>%
  filter(!grepl("Viruses", Taxonomic_lineage)) %>%
  filter(!grepl("\\.", Organism)) %>%
  mutate(Organism = word(Organism, 1, 2)) %>% 
  mutate(Organism = case_when(
    str_detect(Organism, "Phlyctimantis maculatus") ~ "Hylambates maculatus",
    str_detect(Organism, "Nyctimystes infrafrenatus") ~ "Litoria infrafrenata",
    str_detect(Organism, "Boana punctata") ~ "Hypsiboas punctatus",
                          TRUE ~ Organism)) %>%
  mutate(Organism = str_replace(Organism, "Ranoidea", "Litoria")) %>%
  mutate(Organism = str_replace_all(Organism, " ", "_")) %>%
  ungroup()


write_lines(unique(amp_database$Organism), "cache/amp_db_organism_list.txt")
```

There are a total of 788 organisms in this AMP database.

Read in the tree from [TimeTree](http://timetree.org/)

Normalise names: 221 organisms had unresolved names of which 81 were
replaced with different names and the remaining 140 organisms
potentially were not in the TimeTree database at the time.

``` r
timetree <- read.tree("data/amp_db_organism_list.nwk")

timetree_tibble <- as_tibble(timetree) %>%
   mutate(label = case_when(
    str_detect(label, "Halobacterium_salinarum") ~ "Haloarchaeon_S8a",
    str_detect(label, "_aerogenes") ~ "Klebsiella_pneumoniae",
    str_detect(label, "Achromobacter_denitrificans") ~ "Achromobacter_lyticus",
    str_detect(label, "Lactobacillus_casei") ~ "Lactobacillus_paracasei",
    str_detect(label, "Actinoplanes_missouriensis") ~ "Actinoplanes_garbadinensis",
    str_detect(label, "Agrocybe_praecox") ~ "Agrocybe_cylindracea",
    str_detect(label, "Coprinopsis_lagopus") ~ "Coprinopsis_cinerea",
    str_detect(label, "Ganoderma_sinense") ~ "Pleurotus_sajor-caju",
    str_detect(label, "Peziza_vesiculosa") ~ "Pseudoplectania_nigrella",
    str_detect(label, "Microascus_cirrosus") ~ "Pseudallescheria_apiosperma",
    str_detect(label, "Colletotrichum_acutatum") ~ "Colletotrichum_dematium",
    str_detect(label, "Penicillium_javanicum") ~ "Penicillium_rubens",
    str_detect(label, "Mucor_racemosus") ~ "Rhizomucor_pusillus",
    str_detect(label, "Schistosoma_japonicum") ~ "Schistosoma_mansoni",
    str_detect(label, "Molgula_tectiformis") ~ "Halocynthia_aurantium",
    str_detect(label, "Paralichthys_dentatus") ~ "Paralichthys_olivaceus",
    str_detect(label, "Siganus_vulpinus") ~ "Siganus_canaliculatus",
    str_detect(label, "Opsanus_tau") ~ "Thalassophryne_nattereri",
    str_detect(label, "Cynops_pyrrhogaster") ~ "Cynops_fudingensis ",
    str_detect(label, "Paracentrotus_lividus") ~ "Echinus_esculentus",
    str_detect(label, "Leptasterias_hexactis") ~ "Asterias_rubens",
    str_detect(label, "Scolopendra_cingulata") ~ "Scolopendra_subspinipes",
    str_detect(label, "Panulirus_interruptus") ~ "Panulirus_argus",
    str_detect(label, "Chionoecetes_opilio") ~ "Hyas_araneus",
    str_detect(label, "Anax_junius") ~ "Aeshna_cyanea",
    str_detect(label, "Anoplius_aethiops") ~ "Anoplius_samariensis",
    str_detect(label, "Osmia_lignaria") ~ "Osmia_rufa",
    str_detect(label, "Lasioglossum_calceatum") ~ "Lasioglossum_laticeps",
    str_detect(label, "Macropis_europaea") ~ "Macropis_fulvipes",
    str_detect(label, "Tetramorium_caespitum") ~ "Tetramorium_bicarinatum",
    str_detect(label, "Formica_moki") ~ "Formica_aquilonia",
    str_detect(label, "Odontomachus_clarus") ~ "Odontomachus_monticola",
    str_detect(label, "Polistes_tenebricosus") ~ "Polistes_hebraeus",
    str_detect(label, "Stizoides_foxi") ~ "Sphecius_speciosus ",
    str_detect(label, "Aleiodessp.BOLD") ~ "Pimpla_disparis",
    str_detect(label, "Circellium_bacchus") ~ "Copris_tripartitus",
    str_detect(label, "Osmoderma_eremita") ~ "Trypoxylus_dichotomus",
    str_detect(label, "Tetraopes_tetrophthalmus") ~ "Acalolepta_luxuriosa",
    str_detect(label, "Chrysolina_hyperici") ~ "Gastrophysa_atrocyanea",
    str_detect(label, "Culicoides_arakawae") ~ "Culicoides_sonorensis",
    str_detect(label, "Phlebotomus_papatasi") ~ "Phlebotomus_duboscqi",
    str_detect(label, "Sarcophaga_bullata") ~ "Sarcophaga_peregrina",
    str_detect(label, "Tinea_columbariella") ~ "Oiketicus_kirbyi",
    str_detect(label, "Manga_basilinea") ~ "Mamestra_brassicae",
    str_detect(label, "Heliothis_terracottoides") ~ "Heliothis_virescens",
    str_detect(label, "Saturnia_mendocino") ~ "Antheraea_mylitta",
    str_detect(label, "Accinctapubes_albifasciata") ~ "Galleria_mellonella",
    str_detect(label, "Aeschyntelus_notatus") ~ "Riptortus_clavatus",
    str_detect(label, "Brochymena_sp._WCW-2003") ~ "Podisus_maculiventris ",
    str_detect(label, "Platypleura_capensis") ~ "Cicada_flammata",
    str_detect(label, "Cryptotympana_atrata") ~ "Cryptotympana_dubia",
    str_detect(label, "Acyrthosiphon_pisum") ~ "Megoura_viciae",
    str_detect(label, "Ixodes_hexagonus") ~ "Ixodes_sinensis",
    str_detect(label, "Rhipicephalus_sanguineus") ~ "Rhipicephalus_microplus",
    str_detect(label, "Ornithodoros_moubata") ~ "Argas_monolakensis",
    str_detect(label, "Loxosceles_laeta") ~ "Loxosceles_gaucho",
    str_detect(label, "Lycosa_tarantula") ~ "Oxyopes_takobius",
    str_detect(label, "Vaejovis_bandido") ~ "Mesomexovis_subcristatus",
    str_detect(label, "Metaphire_sieboldi") ~ "Metaphire_guillelmi",
    str_detect(label, "Poecilobdella_manillensis") ~ "Hirudo_medicinalis",
    str_detect(label, "Nereis_denhamensis") ~ "Perinereis_aibuhitensis",
    str_detect(label, "Capitella_teleta") ~ "Arenicola_marina",
    str_detect(label, "Conus_purpurascens") ~ "Conus_mustelinus",
    str_detect(label, "Globba_radicalis") ~ "Curcuma_longa",
    str_detect(label, "Xanthosoma_helleborifolium") ~ "Xanthosoma_sagittifolium",
    str_detect(label, "Macadamia_tetraphylla") ~ "Macadamia_integrifolia",
    str_detect(label, "Brassica_rapa") ~ "Brassica_campestris",
    str_detect(label, "Cucurbita_pepo") ~ "Cucurbita_maxima",
    str_detect(label, "Inga_cf._edulis_Klitgaard_677") ~ "Inga_vera",
    str_detect(label, "Parkia_timoriana") ~ "Parkia_pendula",
    str_detect(label, "Cullen_australasicum") ~ "Cullen_corylifolium",
    str_detect(label, "Arachis_major") ~ "Arachis_hypogaea",
    str_detect(label, "Heuchera_merriamii") ~ "Heuchera_sanguinea",
    str_detect(label, "Portulaca_grandiflora") ~ "Basella_alba",
    str_detect(label, "Amaranthus_palmeri") ~ "Amaranthus_caudatus ",
    str_detect(label, "Bidens_cernua") ~ "Dahlia_merckii",
    str_detect(label, "Geophila_obvallata") ~ "Chassalia_chartacea ",
    str_detect(label, "Psychotria_kirkii") ~ "Psychotria_longipes",
    str_detect(label, "Lippia_javanica") ~ "Lippia_sidoides",
    str_detect(label, "Diospyros_cavalcantei") ~ "Diospyros_texana",
    str_detect(label, "_gnavus") ~ "Ruminococcus_gnavus",
    str_detect(label, "Mycobacterium_phlei") ~ "Mycolicibacterium_phlei",
    str_detect(label, "Kassina_maculata") ~ "Hylambates_maculatus",
    str_detect(label, "Rana_septentrionalis") ~ "Lithobates_septentrionalis",
    str_detect(label, "Rana_sevosa") ~ "Lithobates_sevosus",
    str_detect(label, "Rana_clamitans") ~ "Lithobates_clamitans",
    str_detect(label, "Rana_berlandieri") ~ "Lithobates_berlandieri",
    str_detect(label, "Rana_sphenocephala") ~ "Lithobates_sphenocephalus",
    str_detect(label, "Rana_palustris") ~ "Lithobates_palustris",
    str_detect(label, "Rana_pipiens") ~ "Lithobates_pipiens",
    str_detect(label, "Rana_sylvatica") ~ "Lithobates_sylvaticus",
    str_detect(label, "Rana_catesbeiana") ~ "Lithobates_catesbeianus",
    str_detect(label, "Rugosa_rugosa") ~ "Nidirana_pleuraden",
    str_detect(label, "Babina_okinavana") ~ "Glandirana_rugosa",
    str_detect(label, "Babina_pleuraden") ~ "Nidirana_pleuraden",
    str_detect(label, "Eupemphix_nattereri") ~ "Physalaemus_nattereri",
    str_detect(label, "Hypsiboas_albopunctatus") ~ "Boana_albopunctata",
    str_detect(label, "Phyllomedusa_sauvagii") ~ "Phyllomedusa_sauvagei",
    str_detect(label, "Phyllomedusa_hypochondrialis") ~ "Pithecopus_hypochondrialis",
    str_detect(label, "Phyllomedusa_oreades") ~ "Pithecopus_oreades",
    str_detect(label, "Pachymedusa_dacnicolor") ~ "Agalychnis_dacnicolor",
    str_detect(label, "Hylomantis_lemur") ~ "Agalychnis_lemur",
    str_detect(label, "Cercopithecus_preussi") ~ "Allochrocebus_preussi",
    str_detect(label, "Fenneropenaeus_merguiensis") ~ "Penaeus_merguiensis",
    str_detect(label, "Litopenaeus_vannamei") ~ "Penaeus_vannamei",
    str_detect(label, "Litopenaeus_setiferus") ~ "Penaeus_setiferus",
    str_detect(label, "Neobellieria_bullata") ~ "Sarcophaga_peregrina",
    str_detect(label, "Chaetomium_elatum") ~ "Chaetomium_globosum",
    str_detect(label, "Sylvirana_latouchii") ~ "Hylarana_latouchii",
    str_detect(label, "Heliophila_hurkana") ~ "Heliophila_coronopifolia",
                          TRUE ~ label)) %>%
  mutate(label = str_trim(label, side = "right")) # remove trailing whitespace
```

``` r
treelabel_notinAMPs <- timetree_tibble[!timetree_tibble$label %in% amp_database$Organism,]

AMPs_notintreelabel <- amp_database[!amp_database$Organism %in% timetree_tibble$label,]
```

After renaming the majority of species in the tree to match to the
species names used in the AMP database, 472 AMP entries corresponding to
140 organisms remained in SwissProt where the species could not be found
in the TimeTree database. These mostly includes arachnids (scorpions and
spiders), hymenopterans (wasps, bees, ants and sawflies) and anurans
(frogs).

``` r
AMPs_notintreelabel %>% count(Order, sort = TRUE)
```

    ## # A tibble: 37 × 2
    ##    Order               n
    ##    <chr>           <int>
    ##  1 Scorpiones         98
    ##  2 Hymenoptera        91
    ##  3 Anura              83
    ##  4 Araneae            67
    ##  5 Lepidoptera        15
    ##  6 Coleoptera         13
    ##  7 Stolidobranchia    13
    ##  8 Lactobacillales    11
    ##  9 Fabales             9
    ## 10 Poales              9
    ## # … with 27 more rows

Convert the tree tibble back to a tree

``` r
timetree_reworded <- as.treedata(timetree_tibble)

timetree_reworded_phylo <- timetree_reworded@phylo
```

## Pairwise distances between organisms

Calculate pairwise distances between pairs of tips from the tree using
the branch lengths with
[`cophenetic.phylo`](https://rdrr.io/cran/ape/man/cophenetic.phylo.html)

``` r
timetree_distance_matrix <- cophenetic.phylo(timetree_reworded_phylo)
```

Convert matrix to dataframe and left join to AMP database (by organism).
Then calculate the inverse pairwise distance (1/distance) for each
organism

``` r
timetree_dm_df <- as.data.frame(timetree_distance_matrix) %>% rownames_to_column("Organism")

amps_w_distance <- amp_database %>%
  left_join(timetree_dm_df, by = "Organism") %>% 
  mutate(across(.cols = Haloarchaeon_S8a:Cycas_revoluta, function(x) 1/x, .names = "{.col}_inverse" ))
```

![](07_taxonomic_distance_vs_performance_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

**Figure 4.1:** **A** Histogram of pairwise distance **B** Histogram of
the inverse pairwise distance between each faceted organism and 648 out
of the 788 other organisms present in the AMP dataset

``` r
min(amps_w_distance$Gloydius_halys, na.rm = TRUE)
```

    ## [1] 0

``` r
max(amps_w_distance$Gloydius_halys, na.rm = TRUE)
```

    ## [1] 8580

``` r
min(amps_w_distance$Arabidopsis_thaliana, na.rm = TRUE)
```

    ## [1] 0

``` r
max(amps_w_distance$Arabidopsis_thaliana, na.rm = TRUE)
```

    ## [1] 8580

``` r
amps_w_distance %>% 
  select(Arabidopsis_thaliana) %>% 
  count(Arabidopsis_thaliana, sort = TRUE)
```

    ## # A tibble: 33 × 2
    ##    Arabidopsis_thaliana     n
    ##                   <dbl> <int>
    ##  1                2992.   979
    ##  2                  NA    479
    ##  3                2992.   407
    ##  4                2992.   352
    ##  5                   0    294
    ##  6                2992.   221
    ##  7                8580.   188
    ##  8                2992.    58
    ##  9                2992.    54
    ## 10                 321.    42
    ## # … with 23 more rows

``` r
amps_w_distance %>% 
  select(Gloydius_halys) %>% 
  count(Gloydius_halys, sort = TRUE)
```

    ## # A tibble: 56 × 2
    ##    Gloydius_halys     n
    ##             <dbl> <int>
    ##  1           704.   549
    ##  2            NA    479
    ##  3          2992.   362
    ##  4          1593.   198
    ##  5           624.   197
    ##  6          8580    188
    ##  7           624.   149
    ##  8           624.   139
    ##  9          1593.   113
    ## 10           704.    74
    ## # … with 46 more rows

``` r
amps_w_distance %>% 
  select(Mus_musculus) %>% 
  count(Mus_musculus, sort = TRUE)
```

    ## # A tibble: 60 × 2
    ##    Mus_musculus     n
    ##           <dbl> <int>
    ##  1         704.   497
    ##  2          NA    479
    ##  3        2992.   362
    ##  4        8580.   188
    ##  5        1593.   156
    ##  6         180.   148
    ##  7         193.   118
    ##  8        1593.   113
    ##  9           0    104
    ## 10         704.    74
    ## # … with 50 more rows

``` r
amps_w_distance %>% 
  select(Bos_taurus) %>% 
  count(Bos_taurus, sort = TRUE)
```

    ## # A tibble: 58 × 2
    ##    Bos_taurus     n
    ##         <dbl> <int>
    ##  1       704.   497
    ##  2        NA    479
    ##  3      2992.   363
    ##  4      1593.   198
    ##  5      8580    188
    ##  6       193.   176
    ##  7       193.   148
    ##  8      1593.   113
    ##  9       193.    80
    ## 10       704.    74
    ## # … with 48 more rows

Examine values. There are 479 NA values for each column in both the
normal distance and inverse distance. The inverse distance columns
contain `Inf` values for each organism. The `Inf` values match up to the
number of AMPs for these organisms. In the original distance matrix
these values would have been 0, as an organism has 0 taxonomic distance
to itself. Because division by 0 is prohibited, the inverse distances
for these rows subsequently become `Inf` values.

``` r
amps_w_distance %>% 
  select(all_of(organism_selection)) %>%
  summarise(across(everything(), ~ sum(is.na(.)))) %>%
  pivot_longer(cols = everything(), names_to = "Organism", values_to = "NA values count")
```

    ## # A tibble: 13 × 2
    ##    Organism                 `NA values count`
    ##    <chr>                                <int>
    ##  1 Mus_musculus                           479
    ##  2 Homo_sapiens                           479
    ##  3 Bos_taurus                             479
    ##  4 Oryctolagus_cuniculus                  479
    ##  5 Ornithorhynchus_anatinus               479
    ##  6 Gallus_gallus                          479
    ##  7 Oncorhynchus_mykiss                    479
    ##  8 Drosophila_melanogaster                479
    ##  9 Penaeus_vannamei                       479
    ## 10 Bombyx_mori                            479
    ## 11 Arabidopsis_thaliana                   479
    ## 12 Lithobates_catesbeianus                479
    ## 13 Escherichia_coli                       479

``` r
amps_w_distance %>% 
  select(all_of(organism_selection_inverse)) %>%
  summarise(across(everything(), ~ sum(is.infinite(.)))) %>%
  pivot_longer(cols = everything(), names_to = "Organism", values_to = "Inf values count")
```

    ## # A tibble: 13 × 2
    ##    Organism                         `Inf values count`
    ##    <chr>                                         <int>
    ##  1 Mus_musculus_inverse                            104
    ##  2 Homo_sapiens_inverse                             96
    ##  3 Bos_taurus_inverse                               58
    ##  4 Oryctolagus_cuniculus_inverse                    17
    ##  5 Ornithorhynchus_anatinus_inverse                 11
    ##  6 Gallus_gallus_inverse                            25
    ##  7 Oncorhynchus_mykiss_inverse                      12
    ##  8 Drosophila_melanogaster_inverse                  23
    ##  9 Penaeus_vannamei_inverse                         18
    ## 10 Bombyx_mori_inverse                              13
    ## 11 Arabidopsis_thaliana_inverse                    294
    ## 12 Lithobates_catesbeianus_inverse                  13
    ## 13 Escherichia_coli_inverse                         30

## Summed inverse pairwise distance vs. AUPRC

Replace `Inf` values to `NA` and sum inverse distance values. Then
rename columns to original organism names and transform to longer format
so it can be more easily combined to the AUPRC results dataframe

``` r
summed_inverse_distance <- amps_w_distance %>%
  map_df(function(x) replace(x, is.infinite(x), NA)) %>%
  summarise(across(.cols = all_of(organism_selection_inverse), sum, na.rm = TRUE, .names = "{.col}_sum")) %>% 
  rename("Mus_musculus" = Mus_musculus_inverse_sum,
         "Homo_sapiens" = Homo_sapiens_inverse_sum,
         "Bos_taurus" = Bos_taurus_inverse_sum,
         "Oryctolagus_cuniculus" = Oryctolagus_cuniculus_inverse_sum,
         "Ornithorhynchus_anatinus" = Ornithorhynchus_anatinus_inverse_sum,
         "Gallus_gallus" = Gallus_gallus_inverse_sum,
         "Oncorhynchus_mykiss" = Oncorhynchus_mykiss_inverse_sum,
         "Drosophila_melanogaster" = Drosophila_melanogaster_inverse_sum,
         "Penaeus_vannamei" = Penaeus_vannamei_inverse_sum,
         "Bombyx_mori" = Bombyx_mori_inverse_sum,
         "Arabidopsis_thaliana" = Arabidopsis_thaliana_inverse_sum,
         "Lithobates_catesbeianus" = Lithobates_catesbeianus_inverse_sum,
         "Escherichia_coli" = Escherichia_coli_inverse_sum) %>%
  pivot_longer(cols = everything(), names_to = "Organism", values_to = "Inverse_distance_sum")
```

also get the normal distance sum to see how that works

``` r
summed_distance <- amps_w_distance %>%
  summarise(across(.cols = all_of(organism_selection), sum, na.rm = TRUE, .names = "{.col}_sum")) %>% 
  rename("Mus_musculus" = Mus_musculus_sum,
         "Homo_sapiens" = Homo_sapiens_sum,
         "Bos_taurus" = Bos_taurus_sum,
         "Oryctolagus_cuniculus" = Oryctolagus_cuniculus_sum,
         "Ornithorhynchus_anatinus" = Ornithorhynchus_anatinus_sum,
         "Gallus_gallus" = Gallus_gallus_sum,
         "Oncorhynchus_mykiss" = Oncorhynchus_mykiss_sum,
         "Drosophila_melanogaster" = Drosophila_melanogaster_sum,
         "Penaeus_vannamei" = Penaeus_vannamei_sum,
         "Bombyx_mori" = Bombyx_mori_sum,
         "Arabidopsis_thaliana" = Arabidopsis_thaliana_sum,
         "Lithobates_catesbeianus" = Lithobates_catesbeianus_sum,
         "Escherichia_coli" = Escherichia_coli_sum) %>%
  pivot_longer(cols = everything(), names_to = "Organism", values_to = "Distance_sum")
```

Read in previously calculated AUPRC values (see
03\_blast\_and\_prediction.Rmd) for the 13 organisms whose AMP count
refer to AMPs in their respective proteomes found via the UniProt
“Antimicrobial” keyword, and for the nine organisms, whose AMP count
refer to AMPs in their respective proteomes that overlap with the AMPs
in the AMP database.

``` r
methods_auprc_13_wide_w_count <- readRDS("cache/methods_auprc13_w_totalAMPs_wide.rds")
methods_auprc_9_wide_w_count <- readRDS("cache/methods_auprc9_w_totalAMPs_wide.rds")
```

Add the inverse distance metric to the performance and AMP count
datasets

``` r
auprc_and_distance_metric_wAMPcount_13 <- left_join(methods_auprc_13_wide_w_count, summed_inverse_distance, by = "Organism") %>% left_join(summed_distance)
auprc_and_distance_metric_wAMPcount_9 <- left_join(methods_auprc_9_wide_w_count, summed_inverse_distance, by = "Organism") %>% left_join(summed_distance)
```

**Table 7.1:** The summed inverse distance, summed distance and AMP
count for each organism and the AUPRC for the classification and BLAST
methods for 13 organisms

| Organism                  | Classification | BLAST | Total\_AMPs\_in\_test | Inverse\_distance\_sum | Distance\_sum |
|:--------------------------|---------------:|------:|----------------------:|-----------------------:|--------------:|
| Mus\_musculus             |          0.365 | 0.305 |                   132 |               6.218183 |       4773364 |
| Homo\_sapiens             |          0.286 | 0.391 |                   116 |               9.755967 |       4761142 |
| Bos\_taurus               |          0.370 | 0.299 |                   117 |               7.218126 |       4790223 |
| Oryctolagus\_cuniculus    |          0.220 | 0.205 |                    83 |               5.391011 |       4796781 |
| Ornithorhynchus\_anatinus |          0.157 | 0.097 |                    27 |               3.709931 |       4905419 |
| Gallus\_gallus            |          0.435 | 0.121 |                    30 |               3.239859 |       5038456 |
| Oncorhynchus\_mykiss      |          0.071 | 0.108 |                    17 |               2.439460 |       5376383 |
| Drosophila\_melanogaster  |          0.053 | 0.193 |                    33 |               3.136414 |       6237477 |
| Penaeus\_vannamei         |          0.014 | 0.064 |                     3 |               1.783726 |       6366313 |
| Bombyx\_mori              |          0.065 | 0.140 |                    25 |               1.953821 |       6250821 |
| Arabidopsis\_thaliana     |          0.318 | 0.041 |                   296 |               1.966165 |       8043766 |
| Lithobates\_catesbeianus  |          0.022 | 0.214 |                    12 |               6.894414 |       4792522 |
| Escherichia\_coli         |          0.227 | 0.500 |                     4 |             148.720591 |      23568582 |

**Table 7.2:** The summed inverse distance, summed distance and AMP
count for each organism and the AUPRC for the classification and BLAST
methods for 9 organisms

| Organism                  | Classification | BLAST | Total\_AMPs\_in\_test | Inverse\_distance\_sum | Distance\_sum |
|:--------------------------|---------------:|------:|----------------------:|-----------------------:|--------------:|
| Mus\_musculus             |          0.398 | 0.332 |                    99 |               6.218183 |       4773364 |
| Homo\_sapiens             |          0.296 | 0.437 |                    95 |               9.755967 |       4761142 |
| Bos\_taurus               |          0.142 | 0.211 |                    54 |               7.218126 |       4790223 |
| Oryctolagus\_cuniculus    |          0.043 | 0.087 |                    17 |               5.391011 |       4796781 |
| Ornithorhynchus\_anatinus |          0.149 | 0.017 |                    11 |               3.709931 |       4905419 |
| Gallus\_gallus            |          0.439 | 0.078 |                    25 |               3.239859 |       5038456 |
| Drosophila\_melanogaster  |          0.058 | 0.186 |                    23 |               3.136414 |       6237477 |
| Bombyx\_mori              |          0.079 | 0.094 |                    13 |               1.953821 |       6250821 |
| Arabidopsis\_thaliana     |          0.323 | 0.041 |                   291 |               1.966165 |       8043766 |

Change back to long format for plotting

``` r
auprc_and_distance_metric_long_13 <- auprc_and_distance_metric_wAMPcount_13 %>% 
  pivot_longer(cols = c(Classification, BLAST), names_to = "Method", values_to = "AUPRC")

auprc_and_distance_metric_long_9 <- auprc_and_distance_metric_wAMPcount_9 %>% 
  pivot_longer(cols = c(Classification, BLAST), names_to = "Method", values_to = "AUPRC")
```

``` r
ggplot(auprc_and_distance_metric_long_13, aes(x = Inverse_distance_sum, y = AUPRC)) +
  geom_line(aes(linetype = Method)) +
  geom_point(aes(colour = Organism, size = Total_AMPs_in_test)) +
  labs(x = "The sum of the inverse pairwise distance", colour = "", linetype = "") +
  theme_classic() +
  theme(legend.position = "bottom") +
  scale_colour_manual(values = watlington(13)) +
  guides(colour = guide_legend(label.theme = element_text(face = "italic", size = 9)))
```

![](07_taxonomic_distance_vs_performance_files/figure-gfm/unnamed-chunk-20-1.png)<!-- -->

**Figure 7.1:** **A** Scatter and line plot of the summed inverse
pairwise distance and the AUPRC for each AMP finding method for AMPs in
different organisms. The size of points depends on the number of AMPs in
the organism, represented by the AMP\_count.

``` r
ggplot(filter(auprc_and_distance_metric_long_13, Organism != "Escherichia_coli"), aes(x = Inverse_distance_sum, y = AUPRC)) +
  geom_line(aes(linetype = Method)) +
  geom_point(aes(colour = Organism, size = Total_AMPs_in_test)) +
  labs(x = "The sum of the inverse pairwise distance", colour = "", linetype = "") +
  theme_classic() +
  theme(legend.position = "bottom") +
  scale_colour_manual(values = watlington(12)) +
  guides(colour = guide_legend(label.theme = element_text(face = "italic", size = 9)))
```

![](07_taxonomic_distance_vs_performance_files/figure-gfm/unnamed-chunk-21-1.png)<!-- -->

**Figure 7.2:** Same as Figure 7.1, excluding *E. coli*

``` r
ggplot(auprc_and_distance_metric_long_9, aes(x = Inverse_distance_sum, y = AUPRC)) +
  geom_line(aes(linetype = Method)) +
  geom_point(aes(colour = Organism, size = Total_AMPs_in_test)) +
  labs(x = "The sum of the inverse pairwise distance", colour = "", linetype = "") +
  theme_classic() +
  theme(legend.position = "bottom") +
  scale_colour_manual(values = watlington(9)) +
  guides(colour = guide_legend(label.theme = element_text(face = "italic", size = 9)))
```

![](07_taxonomic_distance_vs_performance_files/figure-gfm/unnamed-chunk-22-1.png)<!-- -->

Images were obtained from [phylopic.org](http://phylopic.org/) all
images used were dedicated to the [public
domain](https://creativecommons.org/publicdomain/zero/1.0/) and are not
copyrighted.

``` r
link_to_img <- function(x, size = 25) {
  paste0("<img src='", x, "' width='", size, "'/>")
}
  
auprcplot_img_13 <- auprc_and_distance_metric_wAMPcount_13 %>%
  mutate(images = link_to_img(pics13)) %>% 
  filter(Organism != "Escherichia_coli") %>%
  filter(Organism != "Penaeus_vannamei") %>% 
  pivot_longer(cols = c(Classification, BLAST), names_to = "Method", values_to = "AUPRC") %>%
  ggplot(aes(x = Inverse_distance_sum, y = AUPRC, label = images)) +
  geom_richtext(fill = NA, label.color = NA) +
  geom_line(aes(linetype = Method)) +
  geom_point(aes(size = Total_AMPs_in_test), colour = "forestgreen", shape = 1) +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = "The sum of the inverse pairwise distance", linetype = "") +
  scale_x_continuous(breaks=c(2, 4, 6, 8, 10))


auprcplot_img_9 <- auprc_and_distance_metric_wAMPcount_9 %>%
  mutate(images = link_to_img(pics9)) %>% 
  pivot_longer(cols = c(Classification, BLAST), names_to = "Method", values_to = "AUPRC") %>%
  ggplot(aes(x = Inverse_distance_sum, y = AUPRC, label = images)) +
  geom_richtext(fill = NA, label.color = NA) +
  geom_line(aes(linetype = Method)) +
  geom_point(aes(size = Total_AMPs_in_test), colour = "forestgreen", shape = 1) +
  theme_classic() +
  theme(legend.position = "bottom") +
  labs(x = "The sum of the inverse pairwise distance", linetype = "", size = "AMP count")
```

![](07_taxonomic_distance_vs_performance_files/figure-gfm/unnamed-chunk-25-1.png)<!-- -->

**Figure 7.3** The performance of BLAST and classification models
measured in Area under the Precision-Recall curve (AUPRC) in finding
AMPs in the proteomes of **A)** 11 organisms where AMPs were labelled as
AMPs exclusively by using the “Antimicrobial” keyword from UniProt and
**B** nine organisms where AMPs were labelled as AMPs if annotated with
the UniProt “Antimicrobial” keyword **and** if these AMPs overlapped
with the AMP database generated from SwissProt and the APD, DRAMP or
dbAMP databases.

![](07_taxonomic_distance_vs_performance_files/figure-gfm/unnamed-chunk-28-1.png)<!-- -->

Figure \#: As figure above but using the sum of the normal pairwise
distance
