# Mouse

gunzip -dk data/proteomes/M_musculus_UP000000589_10090.fasta.gz

makeblastdb -in data/proteomes/M_musculus_UP000000589_10090.fasta -dbtype 'prot'

blastp -db data/proteomes/M_musculus_UP000000589_10090.fasta -query cache/Mus_musculus.fasta -outfmt 6 -max_target_seqs 5 -evalue=10 > data/blastp_results/Mus_musculus.blastp


# Cow

gunzip -dk data/proteomes/B_taurus_UP000009136_9913.fasta.gz
 
makeblastdb -in data/proteomes/B_taurus_UP000009136_9913.fasta -dbtype 'prot' 

blastp -db data/proteomes/B_taurus_UP000009136_9913.fasta -query cache/Bos_taurus.fasta -outfmt 6 -max_target_seqs 5 -evalue=10 > data/blastp_results/Bos_taurus.blastp

# Human

gunzip -dk data/proteomes/H_sapiens_UP000005640_9606.fasta.gz 
 
makeblastdb -in data/proteomes/H_sapiens_UP000005640_9606.fasta -dbtype 'prot' 

blastp -db data/proteomes/H_sapiens_UP000005640_9606.fasta -query cache/Homo_sapiens.fasta -outfmt 6 -max_target_seqs 5 -evalue=10 > data/blastp_results/Homo_sapiens.blastp

# Rabbit

gunzip -dk data/proteomes/O_cuniculus_UP000001811_9986.fasta.gz 
 
makeblastdb -in data/proteomes/O_cuniculus_UP000001811_9986.fasta -dbtype 'prot' 

blastp -db data/proteomes/O_cuniculus_UP000001811_9986.fasta -query cache/Oryctolagus_cuniculus.fasta -outfmt 6 -max_target_seqs 5 -evalue=10 > data/blastp_results/Oryctolagus_cuniculus.blastp

# Platypus

gunzip -dk data/proteomes/O_anatinus_UP000002279_9258.fasta.gz 
 
makeblastdb -in data/proteomes/O_anatinus_UP000002279_9258.fasta -dbtype 'prot' 

blastp -db data/proteomes/O_anatinus_UP000002279_9258.fasta -query cache/Ornithorhynchus_anatinus.fasta -outfmt 6 -max_target_seqs 5 -evalue=10 > data/blastp_results/Ornithorhynchus_anatinus.blastp

# Junglefowl

gunzip -dk data/proteomes/G_gallus_UP000000539_9031.fasta.gz
 
makeblastdb -in data/proteomes/G_gallus_UP000000539_9031.fasta -dbtype 'prot' 

blastp -db data/proteomes/G_gallus_UP000000539_9031.fasta -query cache/Gallus_gallus.fasta -outfmt 6 -max_target_seqs 5 -evalue=10 > data/blastp_results/Gallus_gallus.blastp

# Trout

gunzip -dk data/proteomes/O_mykiss_UP000193380_8022.fasta.gz 
 
makeblastdb -in data/proteomes/O_mykiss_UP000193380_8022.fasta -dbtype 'prot' 

blastp -db data/proteomes/O_mykiss_UP000193380_8022.fasta -query cache/Oncorhynchus_mykiss.fasta -outfmt 6 -max_target_seqs 5 -evalue=10 > data/blastp_results/Oncorhynchus_mykiss.blastp

# Fruitfly 

gunzip -dk data/proteomes/D_melanogaster_UP000000803_7227.fasta.gz 
 
makeblastdb -in data/proteomes/D_melanogaster_UP000000803_7227.fasta -dbtype 'prot' 

blastp -db data/proteomes/D_melanogaster_UP000000803_7227.fasta -query cache/Drosophila_melanogaster.fasta -outfmt 6 -max_target_seqs 5 -evalue=10 > data/blastp_results/Drosophila_melanogaster.blastp

# Shrimp

gunzip -dk data/proteomes/P_vannamei_UP000283509_6689.fasta.gz 
 
makeblastdb -in data/proteomes/P_vannamei_UP000283509_6689.fasta -dbtype 'prot' 

blastp -db data/proteomes/P_vannamei_UP000283509_6689.fasta -query cache/Penaeus_vannamei.fasta -outfmt 6 -max_target_seqs 5 -evalue=10 > data/blastp_results/Penaeus_vannamei.blastp

# Silk moth

gunzip -dk data/proteomes/B_mori_UP000005204_7091.fasta.gz 
 
makeblastdb -in data/proteomes/B_mori_UP000005204_7091.fasta -dbtype 'prot' 

blastp -db data/proteomes/B_mori_UP000005204_7091.fasta -query cache/Bombyx_mori.fasta -outfmt 6 -max_target_seqs 5 -evalue=10 > data/blastp_results/Bombyx_mori.blastp

# Cress

gunzip -dk data/proteomes/A_thaliana_UP000006548_3702.fasta.gz 
 
makeblastdb -in data/proteomes/A_thaliana_UP000006548_3702.fasta -dbtype 'prot' 

blastp -db data/proteomes/A_thaliana_UP000006548_3702.fasta -query cache/Arabidopsis_thaliana.fasta -outfmt 6 -max_target_seqs 5 -evalue=10 > data/blastp_results/Arabidopsis_thaliana.blastp

# Bullfrog

gunzip -dk data/proteomes/L_catesbeianus_UP000228934_8400.fasta.gz 
 
makeblastdb -in data/proteomes/L_catesbeianus_UP000228934_8400.fasta -dbtype 'prot' 

blastp -db data/proteomes/L_catesbeianus_UP000228934_8400.fasta -query cache/Lithobates_catesbeianus.fasta -outfmt 6 -max_target_seqs 5 -evalue=10 > data/blastp_results/Lithobates_catesbeianus.blastp

# Bacteria

gunzip -dk data/proteomes/E_coli_k12_UP000000625_83333.fasta.gz 
 
makeblastdb -in data/proteomes/E_coli_k12_UP000000625_83333.fasta -dbtype 'prot' 

blastp -db data/proteomes/E_coli_k12_UP000000625_83333.fasta -query cache/Escherichia_coli.fasta -outfmt 6 -max_target_seqs 5 -evalue=10 > data/blastp_results/Escherichia_coli.blastp

find data/proteomes/ -type f -not -name '*.gz' -delete

