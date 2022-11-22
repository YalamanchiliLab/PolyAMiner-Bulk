# DEMO Contrast
## 3'UTR Only
## Softclipped + Apriori

python3 /mnt/belinda_local/venkata/data/PolyAMiner-Bulk/PolyA-miner.py -mode bam \
-fasta /mnt/belinda_local/venkata/data/Index_Files/Human/GenomeFasta_GTF/GRCh38.primary_assembly.genome.fa \
-gtf /mnt/belinda_local/venkata/data/Index_Files/Human/GenomeFasta_GTF/gencode.v33.primary_assembly.annotation.gtf \
-p 20 -a 0.65 -outPrefix 3UTROnly -expNovel 1 -t BB -s 2 \
-o /mnt/belinda_local/venkata/data/PolyAMiner-Bulk/Demo/Demo_Results/Demo_3UTROnly_Softclipped+APriori \
-c1 /mnt/belinda_local/venkata/data/PolyAMiner-Bulk/Demo/control1_-.subset.sorted.bam,\
/mnt/belinda_local/venkata/data/PolyAMiner-Bulk/Demo/control2_.subset.sorted.bam,\
/mnt/belinda_local/venkata/data/PolyAMiner-Bulk/Demo/control3_.subset.sorted.bam \
-c2 /mnt/belinda_local/venkata/data/PolyAMiner-Bulk/Demo/treatment1_.subset.sorted.bam,\
/mnt/belinda_local/venkata/data/PolyAMiner-Bulk/Demo/treatment2_.subset.sorted.bam,\
/mnt/belinda_local/venkata/data/PolyAMiner-Bulk/Demo/treatment3_.subset.sorted.bam \
-ignore UTR5,CDS,Intron,UN -apriori_annotations -modelOrganism human \
-visualizeTopNum 10 -visualizeCondition1Name Control -visualizeCondition2Name Treatment
