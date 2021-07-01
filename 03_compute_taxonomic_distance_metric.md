
``` r
library(tidyverse)
library(ape)
```

AMPs (3,371 sequences) were obtained from SwissProt on 01 July 2021
[UniProtKB 2021_03
results](https://www.uniprot.org/uniprot/?query=keyword%3A%22Antimicrobial%20%5BKW-0929%5D%22%20AND%20reviewed%3Ayes&columns=id%2Centry%20name%2Creviewed%2Cprotein%20names%2Cgenes%2Corganism%2Clength%2Ckeyword-id%2Ckeywords%2Cproteome%2Corganism-id%2Clineage(ORDER)%2Csequence%2Cexistence%2Clineage(ALL)&sort=sequence-modified).

UniProt contains a column called “Organism ID”. This ID refers to the
NCBI “taxid” which can be used in the [Taxonomy
Browser](https://www.ncbi.nlm.nih.gov/Taxonomy/CommonTree/wwwcmt.cgi) to
obtain a phylogenetic tree. The organism IDs, rather than the organism
names were used as some organism names listed in UniProt are under a
different name in NCBI, e.g. *Neoponera goeldii* is currently known as
*Pachycondyla goeldii* in the NCBI taxonomy database. However, the
taxonomic ID for this organism is 118888 in both databases.

Read in SwissProt AMPs and save the organism taxonomic IDs to submit to
the NCBI Taxonomy Browser to obtain a tree. Viruses were removed. This
resulted in 832 unique organisms.

``` r
swissprot_amps_july <- read_tsv("data/uniprot-keywordAntimicrobial+[KW-0929]+AND+reviewedyes01July2021.tab.gz") %>% 
  rename(Organism_ID = `Organism ID`) %>%
  rename(Taxonomic_lineage = `Taxonomic lineage (ALL)`) %>%
  filter(!grepl("Viruses", Taxonomic_lineage))


write_lines(unique(swissprot_amps_july$Organism_ID), "cache/organism_taxids.txt")
```

Read in the tree

``` r
tree_text <- readLines("data/phyliptree.phy") %>%
  paste0(collapse="")
tree <- read.tree(text = tree_text)

glimpse(tree)
```

    ## List of 6
    ##  $ edge       : int [1:941, 1:2] 692 693 694 694 695 696 696 696 697 697 ...
    ##  $ edge.length: num [1:941] 4 4 4 4 4 4 4 4 4 4 ...
    ##  $ Nnode      : int 251
    ##  $ node.label : chr [1:251] "Eukaryota" "Metazoa" "Chordata" "Mammalia" ...
    ##  $ tip.label  : chr [1:691] "Potamotrygoncf.henleiKC-2012" "Chinchillalanigera" "Caviaporcellus" "Rattusnorvegicus" ...
    ##  $ root.edge  : num 4
    ##  - attr(*, "class")= chr "phylo"
    ##  - attr(*, "order")= chr "cladewise"

Calculate branch lengths with
[`compute.brlen`](https://rdrr.io/cran/ape/man/compute.brlen.html)

``` r
tree_w_brlen <- compute.brlen(tree)

glimpse(tree_w_brlen)
```

    ## List of 6
    ##  $ edge       : int [1:941, 1:2] 692 693 694 694 695 696 696 696 697 697 ...
    ##  $ edge.length: num [1:941] 0.2014 0.3232 0.4754 0.3681 0.0986 ...
    ##  $ Nnode      : int 251
    ##  $ node.label : chr [1:251] "Eukaryota" "Metazoa" "Chordata" "Mammalia" ...
    ##  $ tip.label  : chr [1:691] "Potamotrygoncf.henleiKC-2012" "Chinchillalanigera" "Caviaporcellus" "Rattusnorvegicus" ...
    ##  $ root.edge  : num 4
    ##  - attr(*, "class")= chr "phylo"
    ##  - attr(*, "order")= chr "cladewise"

Calculate pairwise distances between pairs of tips from the tree using
the branch lengths with
[`cophenetic.phylo`](https://rdrr.io/cran/ape/man/cophenetic.phylo.html)

``` r
tree_distance_matrix <- cophenetic.phylo(tree_w_brlen)

glimpse(tree_distance_matrix)
```

    ##  num [1:691, 1:691] 0 0.951 0.951 0.951 0.951 ...
    ##  - attr(*, "dimnames")=List of 2
    ##   ..$ : chr [1:691] "Potamotrygoncf.henleiKC-2012" "Chinchillalanigera" "Caviaporcellus" "Rattusnorvegicus" ...
    ##   ..$ : chr [1:691] "Potamotrygoncf.henleiKC-2012" "Chinchillalanigera" "Caviaporcellus" "Rattusnorvegicus" ...
