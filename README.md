
# ML vs. homology

Effectiveness of AMP prediction using machine learning versus sequence
similarity

-   [Find species represented with the most AMPs in
    UniProt](01_species_represented_with_most_amps_uniprot.md) : Rmd
    file
    [01\_species\_represented\_with\_most\_amps\_uniprot.Rmd](01_species_represented_with_most_amps_uniprot.Rmd)
-   [Create databases for model training and homology
    searches](02_create_databases.md) : Rmd file
    [02\_create\_databases.Rmd](02_create_databases.Rmd)
-   [Execute BLAST searches and AMP
    predictions](03_blast_and_prediction.md) : Rmd file
    [03\_blast\_and\_prediction.Rmd](03_blast_and_prediction.Rmd)
-   [Compute a taxonomic distance
    metric](04_compute_taxonomic_distance_metric.md) : Rmd file
    [04\_compute\_taxonomic\_distance\_metric.Rmd](04_compute_taxonomic_distance_metric.Rmd)

### The above workflow was repeated for a different AMP database with different organisms

-   [Find species represented with the most AMPs in a constructed AMP
    database and construct positive
    datasets](05_amp_training_data_preparations.md) : Rmd file
    [05\_amp\_training\_data\_preparations.Rmd](05_amp_training_data_preparations.Rmd)
-   [Test BLAST and classification
    methods](06_method_evaluation_on_proteomes.md) : Rmd file
    [06\_method\_evaluation\_on\_proteomes.Rmd](06_method_evaluation_on_proteomes.Rmd)
-   [Use phylogenetic data to assess performance over a taxonomic
    scale](07_taxonomic_distance_vs_performance.md) : Rmd file
    [07\_taxonomic\_distance\_vs\_performance.Rmd](07_taxonomic_distance_vs_performance.Rmd)

The files required to run the code in these Rmd files can be obtained by
clicking [here](https://cloudstor.aarnet.edu.au/plus/s/0QxalNJhHcbUIqk)
or by using the command:

``` bash
wget 'https://cloudstor.aarnet.edu.au/plus/s/0QxalNJhHcbUIqk/download' -O data.tgz
tar -zxvf data.tgz 
```
