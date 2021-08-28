
# Proteins

echo data/uniprot_amps_07Jul21_3371rev42126unrev.tab.gz >> data_list
echo data/uniprot-keywordAntimicrobial+[KW-0929]-filtered-reviewedyes24May21.tab >> data_list
echo data/uniprot-NOT+keyword_Antimicrobial+[KW-0929]+length[5+TO+500]24May21.tab >> data_list
echo data/uniprot_amps_w_amp_dbsJuly21.rds >> data_list
echo data/swissprot_amps_standardaa_90.fasta >> data_list
echo data/swissprot_nonamps_standardaa_90.fasta >> data_list
echo cache/amps_standardaa.fasta90.fasta >> data_list 

# Proteomes

echo data/proteomes/A_thaliana_UP000006548_3702.fasta.gz >> data_list
echo data/proteomes/A_thaliana_proteome-UP000006548.tab.gz >> data_list
echo data/proteomes/B_mori-proteome-UP000005204.tab.gz >> data_list
echo data/proteomes/B_mori_UP000005204_7091.fasta.gz >> data_list
echo data/proteomes/B_taurus-proteome-UP000009136.tab.gz >> data_list
echo data/proteomes/B_taurus_UP000009136_9913.fasta.gz >> data_list
echo data/proteomes/C_familiaris_UP000002254_9615.fasta.gz >> data_list
echo data/proteomes/C_familiaris_proteome-UP000002254.tab.gz >> data_list
echo data/proteomes/D_melanogaster-proteome-UP000000803.tab.gz >> data_list
echo data/proteomes/D_melanogaster_UP000000803_7227.fasta.gz >> data_list
echo data/proteomes/E_coli-proteome-UP000000625.tab.gz >> data_list
echo data/proteomes/E_coli_k12_UP000000625_83333.fasta.gz >> data_list
echo data/proteomes/G_gallus-proteome-UP000000539.tab.gz >> data_list
echo data/proteomes/G_gallus_UP000000539_9031.fasta.gz >> data_list
echo data/proteomes/H_sapiens-proteome-UP000005640.tab.gz >> data_list
echo data/proteomes/H_sapiens_UP000005640_9606.fasta.gz >> data_list
echo data/proteomes/L_catesbeianus-proteome-UP000228934.tab.gz >> data_list
echo data/proteomes/L_catesbeianus_UP000228934_8400.fasta.gz >> data_list
echo data/proteomes/M_musculus-proteome-UP000000589.tab.gz >> data_list
echo data/proteomes/M_musculus_UP000000589_10090.fasta.gz >> data_list
echo data/proteomes/O_anatinus-proteome-UP000002279.tab.gz >> data_list
echo data/proteomes/O_anatinus_UP000002279_9258.fasta.gz >> data_list
echo data/proteomes/O_cuniculus-proteome-UP000001811.tab.gz >> data_list
echo data/proteomes/O_cuniculus_UP000001811_9986.fasta.gz >> data_list
echo data/proteomes/O_mykiss-proteome-UP000193380.tab.gz >> data_list
echo data/proteomes/O_mykiss_UP000193380_8022.fasta.gz >> data_list
echo data/proteomes/P_troglodytes-proteome-UP000002277.tab.gz >> data_list
echo data/proteomes/P_troglodytes_UP000002277_9598.fasta.gz >> data_list
echo data/proteomes/P_vannamei-proteome-UP000283509.tab.gz >> data_list
echo data/proteomes/P_vannamei_UP000283509_6689.fasta.gz >> data_list
echo data/proteomes/R_ferrumequinum-proteome-UP000472240.tab.gz >> data_list
echo data/proteomes/R_ferrumequinum_UP000472240_59479.fasta.gz >> data_list
echo data/proteomes/R_norvegicus-proteome-UP000002494.tab.gz >> data_list
echo data/proteomes/R_norvegicus_UP000002494_10116.fasta.gz >> data_list
echo data/proteomes/S_scrofa-uniprot-proteome-UP000008227.tab.gz >> data_list
echo data/proteomes/S_scrofa_UP000008227_9823.fasta.gz >> data_list

# Models

echo models/Bos_taurus_model.rds >> data_list
echo models/Canis_familiaris_model.rds >> data_list
echo models/Homo_sapiens_model.rds >> data_list
echo models/Mus_musculus_model.rds >> data_list
echo models/Pan_troglodytes_model.rds >> data_list
echo models/Rattus_norvegicus_model.rds >> data_list
echo models/Sus_scrofa_model.rds >> data_list
echo models/bacteria_model.rds >> data_list
echo models/cow_model.rds >> data_list
echo models/cress_model.rds >> data_list
echo models/frog_model.rds >> data_list
echo models/fruitfly_model.rds >> data_list
echo models/general_model.rds >> data_list
echo models/human_model.rds >> data_list
echo models/junglefowl_model.rds >> data_list
echo models/moth_model.rds >> data_list
echo models/mouse_model.rds >> data_list
echo models/platypus_model.rds >> data_list
echo models/rabbit_model.rds >> data_list
echo models/shrimp_model.rds >> data_list
echo models/trout_model.rds >> data_list


# BLAST and prediction results

echo data/blastp_results/Arabidopsis_thaliana.blastp >> data_list
echo data/blastp_results/Bombyx_mori.blastp >> data_list
echo data/blastp_results/Bos_taurus.blastp >> data_list
echo data/blastp_results/Bos_taurus_proteome.blastp >> data_list
echo data/blastp_results/Canis_familiaris_proteome.blastp >> data_list
echo data/blastp_results/Canis_lupus_familiaris.blastp >> data_list
echo data/blastp_results/Drosophila_melanogaster.blastp >> data_list
echo data/blastp_results/Escherichia_coli.blastp >> data_list
echo data/blastp_results/Gallus_gallus.blastp >> data_list
echo data/blastp_results/Homo_sapiens.blastp >> data_list
echo data/blastp_results/Homo_sapiens_proteome.blastp >> data_list
echo data/blastp_results/Lithobates_catesbeianus.blastp >> data_list
echo data/blastp_results/Mus_musculus.blastp >> data_list
echo data/blastp_results/Mus_musculus_proteome.blastp >> data_list
echo data/blastp_results/Oncorhynchus_mykiss.blastp >> data_list
echo data/blastp_results/Ornithorhynchus_anatinus.blastp >> data_list
echo data/blastp_results/Oryctolagus_cuniculus.blastp >> data_list
echo data/blastp_results/Pan_troglodytes.blastp >> data_list
echo data/blastp_results/Pan_troglodytes_proteome.blastp >> data_list
echo data/blastp_results/Penaeus_vannamei.blastp >> data_list
echo data/blastp_results/R_ferrumequinum.blastp >> data_list
echo data/blastp_results/R_ferrumequinum_proteome.blastp >> data_list
echo data/blastp_results/Rattus_norvegicus.blastp >> data_list
echo data/blastp_results/Rattus_norvegicus_proteome.blastp >> data_list
echo data/blastp_results/Sus_scrofa.blastp >> data_list
echo data/blastp_results/Sus_scrofa_proteome.blastp >> data_list


echo cache/proteome_predictions13.rds >> data_list
echo cache/blast_and_pred_roc13.rds >> data_list
echo cache/blast_and_pred_roc_strict_9.rds >> data_list
echo cache/methods_auprc13_w_totalAMPs_wide.rds >> data_list
echo cache/methods_auprc9_w_totalAMPs_wide.rds >> data_list

# Trees

echo data/organism_list.nwk >> data_list
echo data/organism_list_w_bat.nwk >> data_list
echo data/amp_db_organism_list.nwk >> data_list

tar -zcvf data_amp_pred.tgz -T data_list

rm data_list