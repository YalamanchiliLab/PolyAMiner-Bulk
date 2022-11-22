# Contrast 1: Ctrl-siRNA vs Rbm17-siRNA
## 3'UTR Only
## Softclipped Assisted

# python3 /mnt/belinda_local/venkata/data/PolyAMiner-Bulk/PolyA-miner.py -mode bam \
# -fasta /mnt/belinda_local/venkata/data/Index_Files/Human/GenomeFasta_GTF/GRCh38.primary_assembly.genome.fa \
# -gtf /mnt/belinda_local/venkata/data/Index_Files/Human/GenomeFasta_GTF/gencode.v33.primary_assembly.annotation.gtf \
# -p 20 -pa_p 0.6 -pa_a 10 -pa_m 10 -ip_u 30 -ip_d 40 -a 0.65 -novel_d 5000 -outPrefix 3UTROnly -expNovel 1 -t BB -s 2 \
# -o /mnt/belinda_local/venkata/data/PolyAMiner-Bulk/Demo/Demo_Results/CtrlvsRBM17siRNA_3UTROnly_Softclipped_Assisted \
# -c1 /mnt/belinda_local/venkata/data/PolyAMiner-Bulk/Demo/HZ8169_.subset.sorted.bam,\
# /mnt/belinda_local/venkata/data/PolyAMiner-Bulk/Demo/HZ8170_.subset.sorted.bam,\
# /mnt/belinda_local/venkata/data/PolyAMiner-Bulk/Demo/HZ8171_.subset.sorted.bam \
# -c2 /mnt/belinda_local/venkata/data/PolyAMiner-Bulk/Demo/HZ8162_.subset.sorted.bam,\
# /mnt/belinda_local/venkata/data/PolyAMiner-Bulk/Demo/HZ8163_.subset.sorted.bam,\
# /mnt/belinda_local/venkata/data/PolyAMiner-Bulk/Demo/HZ8164_.subset.sorted.bam \
# -ignore UTR5,CDS,Intron,UN -apriori_annotations -modelOrganism human -softclippedNumReads 1 -softclippedNumSamples 1 -verboseLogging 
# -visualizeTopNum 100 -visualizeCondition1Name Control -visualizeCondition2Name RBM17-KD



# Contrast 1: Ctrl-siRNA vs Rbm17-siRNA
## 3'UTR Only
## Softclipped + Apriori

python3 /mnt/belinda_local/venkata/data/PolyAMiner-Bulk/PolyA-miner.py -mode bam \
-fasta /mnt/belinda_local/venkata/data/Index_Files/Human/GenomeFasta_GTF/GRCh38.primary_assembly.genome.fa \
-gtf /mnt/belinda_local/venkata/data/Index_Files/Human/GenomeFasta_GTF/gencode.v33.primary_assembly.annotation.gtf \
-p 20 -pa_p 0.6 -pa_a 10 -pa_m 10 -ip_u 30 -ip_d 40 -a 0.65 -novel_d 5000 -outPrefix 3UTROnly -expNovel 1 -t BB -s 2 \
-o /mnt/belinda_local/venkata/data/PolyAMiner-Bulk/Demo/Demo_Results/CtrlvsRBM17siRNA_3UTROnly_Softclipped+APriori \
-c1 /mnt/belinda_local/venkata/data/PolyAMiner-Bulk/Demo/HZ8169_.subset.sorted.bam,\
/mnt/belinda_local/venkata/data/PolyAMiner-Bulk/Demo/HZ8170_.subset.sorted.bam,\
/mnt/belinda_local/venkata/data/PolyAMiner-Bulk/Demo/HZ8171_.subset.sorted.bam \
-c2 /mnt/belinda_local/venkata/data/PolyAMiner-Bulk/Demo/HZ8162_.subset.sorted.bam,\
/mnt/belinda_local/venkata/data/PolyAMiner-Bulk/Demo/HZ8163_.subset.sorted.bam,\
/mnt/belinda_local/venkata/data/PolyAMiner-Bulk/Demo/HZ8164_.subset.sorted.bam \
-ignore UTR5,CDS,Intron,UN -apriori_annotations -modelOrganism human \
-visualizeTopNum 100 -visualizeCondition1Name Control -visualizeCondition2Name RBM17-KD
