
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

# Rat

gunzip -dk data/proteomes/R_norvegicus_UP000002494_10116.fasta.gz 

makeblastdb -in data/proteomes/R_norvegicus_UP000002494_10116.fasta -dbtype 'prot' 

blastp -db data/proteomes/R_norvegicus_UP000002494_10116.fasta -query cache/Rattus_norvegicus.fasta -outfmt 6 -max_target_seqs 5 -evalue=10 > data/blastp_results/Rattus_norvegicus.blastp

# Chimpanzee

gunzip -dk data/proteomes/P_troglodytes_UP000002277_9598.fasta.gz 

makeblastdb -in data/proteomes/P_troglodytes_UP000002277_9598.fasta -dbtype 'prot' 

blastp -db data/proteomes/P_troglodytes_UP000002277_9598.fasta -query cache/Pan_troglodytes.fasta -outfmt 6 -max_target_seqs 5 -evalue=10 > data/blastp_results/Pan_troglodytes.blastp

# Pig

# Dog

# Bat (using general AMP set)

find data/proteomes/ -type f -not -name '*.gz' -delete
