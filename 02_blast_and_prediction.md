
# Using BLAST to find AMPs from AMP query database on proteomes

## Proteomes

The proteomes were downloaded from
[UniProt](https://www.uniprot.org/proteomes) on 10 June 2021.

**Table 1:** Proteome information from organisms used as subjects

| Organism Name  | Reference proteome ID                                        | Total proteins | Reviewed | Unreviewed |
|----------------|--------------------------------------------------------------|----------------|----------|------------|
| *Mus musculus* | [UP000000589](https://www.uniprot.org/proteomes/UP000000589) | 55,470         | 17,068   | 38,402     |
|                |                                                              |                |          |            |
|                |                                                              |                |          |            |
|                |                                                              |                |          |            |

## BLAST searches

The [BLAST+](https://pubmed.ncbi.nlm.nih.gov/20003500/) version used was
blast 2.11.0, build Nov 17 2020 for MacOS.

Each proteome was used to make a local BLAST database using
`makeblastdb`. This proteome database was then used to query the AMP
dataset with `blastp`, Protein-Protein BLAST 2.11.0+

``` bash
gunzip -dc data/proteomes/M_musculus-proteome-UP000000589.fasta.gz | makeblastdb -in - -title M_musculus-proteome-UP000000589 -dbtype prot -out cache/M_musculus-proteome-UP000000589.fasta
```

``` bash
blastp -db cache/M_musculus-proteome-UP000000589.fasta -query cache/Mus_musculus.fasta -outfmt 6 -max_target_seqs 5 -evalue=0.00001 > data/blastp_results/Mus_musculus.blastp
```

1.  qaccver - Query accession.version (AMPs list ID)
2.  saccver - Subject accession.version (Reference proteome ID)
3.  pident - Percentage of identical matches
4.  length - Alignment length (sequence overlap)
5.  mismatch - Number of mismatches
6.  gapopen - Number of gap openings
7.  qstart - Start of alignment in query
8.  qend - End of alignment in query
9.  sstart - Start of alignment in subject
10. send - End of alignment in subject
11. evalue - Expect value (the smaller the evalue, the better the
    homology match)
12. bitscore - Bit score (the higher the bitscore, the better the
    sequence similarity)

``` r
blast_colnames <- c("qaccver","saccver","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore")

mouse_amps_blast <- read_tsv("data/blastp_results/Mus_musculus.blastp", col_names = blast_colnames)

mouse_amps_blast %>% group_by(qaccver) %>% slice_max(n = 2, order_by = bitscore)
```

    ## # A tibble: 59 x 12
    ## # Groups:   qaccver [28]
    ##    qaccver  saccver     pident length mismatch gapopen qstart  qend sstart  send
    ##    <chr>    <chr>        <dbl>  <dbl>    <dbl>   <dbl>  <dbl> <dbl>  <dbl> <dbl>
    ##  1 CAMP_MA… sp|P51437|…   56.2    144       57       1     30   167     27   170
    ##  2 CAMP_MA… sp|O08692|…   34.4    122       78       1     11   132      3   122
    ##  3 CH3L1_B… sp|Q61362|…   69.9    382      113       2      1   381      9   389
    ##  4 CH3L1_B… sp|Q9D7Q1|…   49.6    387      184       5      1   380      1   383
    ##  5 CMGA_HU… sp|P26339|…   60.8    490      132      11      1   457      1   463
    ##  6 DB119_B… sp|Q8K3I8|…   48.5     68       35       0     15    82     15    82
    ##  7 DFP_LOC… tr|A0A0G2J…   28.2    163      106       6      1   159      8   163
    ##  8 DFP_LOC… sp|Q8K385|…   28.2    163      106       6      1   159      8   163
    ##  9 H2A_ONC… sp|Q8CGP6|…   94.5    128        7       0      1   128      1   128
    ## 10 H2A_ONC… sp|C0HKE8|…   94.5    128        7       0      1   128      1   128
    ## # … with 49 more rows, and 2 more variables: evalue <dbl>, bitscore <dbl>
