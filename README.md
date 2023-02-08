# PolyAMiner-Bulk

This repository includes the PolyAMiner-Bulk computational tool from 'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'. Please cite our paper if you have used our computational tool, any of our machine learning models, or code snippets. Additionally, this repository is actively under development, so please kindly report any issues or feature requests.

In this package, we provide the following resources: 

(1) Source code of PolyAMiner-Bulk

(2) Test scripts delineating key usage scenarios

(3) Fine-tuned CPAS-BERT models for both human and mouse model organisms.

## Citation

If you have used PolyAMiner-Bulk in your research, please kindly cite the following publication:

XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


## Step 1: Setup the compute environment

### 1.1: Create and activate a new virtual environment
It is highly recommended that you use Anaconda to build a python virtual environment. Also, please make sure you have at least one NVIDIA GPU with Linux x86_64 Driver Version >= 410.48 (compatible with CUDA 10.0).

```
conda create -n cpasbert python=3.9
conda activate cpasbert
```

### 1.2: Install the necessary package dependencies

```
<!-- conda install pytorch torchvision cudatoolkit=10.0 -c pytorch -->
conda install pytorch==1.11.0 torchvision==0.12.0 torchaudio==0.11.0 cudatoolkit=10.0 -c pytorch

git clone https://github.com/venkatajonnakuti/PolyAMiner-Bulk
cd PolyAMiner-Bulk/lib/DNABERT
python3 -m pip install --editable .
Rscript installPkgs.R

conda install pandas statsmodels
pip3 install -U scikit-learn
pip install pysam tokenizers
pip install tensorboard
conda install -c bioconda pyfasta
conda install -c bioconda gtfparse
conda install -c bioconda pybedtools
conda install -c bioconda pybigwig
conda install -c bioconda subread
conda install -c bioconda samtools
conda install -c bioconda bedtools
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
sudo apt install build-essential
conda install -c anaconda seaborn 
sudo apt-get install gfortran
pip install deeptools
pip install pygenometracks==3.6

<!-- 
conda install pandas
conda install statsmodels
pip3 install -U scikit-learn
pip install pysam
conda install -c bioconda pyfasta
conda install -c bioconda gtfparse
pip install tokenizers
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
sudo apt install build-essential
pip install transformers==2.5.0
pip install tensorboard
conda install -c anaconda seaborn 
conda install -c bioconda pybedtools
conda install -c bioconda pybigwig
conda install -c bioconda subread
conda install -c bioconda samtools
conda install -c conda-forge r-kernsmooth
BiocManager::install("sva")
test again
 -->
```

## Step 2: Download trained ML-models + ML-dependencies + necessary reference files

### 2.1: Download "6-new-12w-0" from https://bcm.box.com/s/1lo85ig6eop2a5e12d52ytd6q98fg4xa and copy the extracted folder (NOTE: not the folder contents!) into the /lib folder

### 2.2: Download "CPASBERT_TrainedModels" from https://bcm.box.com/s/d1yt6iz95r18rfs2kf22thc58r3qipxg and copy the extracted folder (NOTE: not the folder contents!) into the /lib folder

### 2.3: Download organism-specific fasta and gtf reference files from https://bcm.box.com/s/vyze88bg9wr71kkxmz5yi4hvi8vstscd OR directly from Ensemble: http://uswest.ensembl.org/info/data/ftp/index.html/ and copy the extracted folder to a readable and writeable location on your computer. NOTE: Keep the locations of these fasta and gtf reference files on hand as they are required parameters for PolyAMiner-Bulk.

## Step 3: Understand PolyAMiner-Bulk Command-Line Parameters

#### Base Parameters
-mode = Run mode options: \'bam\' to start from mapped data, \'fastq\' to start from raw data' (string)

-index = Reference genome bowtie2 index. NOTE: Valid for -mode fastq ONLY! (string)

-d = Base directory of input fastq files. NOTE: Valid for -mode fastq ONLY! (string)

-o = Output directory; default = 'PolyAminer_OUT' (string)

-c1 = Comma-separated list of condition1 files. Full path for BAMs (index files are also expected) or just file names for fastq (string)

-c2 = Comma-separated list of condition2 files. Full path for BAMs (index files are also expected) or Just file names for fastq (string)

-s = Strand information. Use 0 for un-stranded, 1 for fwd-stranded, and 2 for rev-stranded (integer)

#### Required Ref. Files
-fasta = Reference fasta file

-gtf = Reference gtf file

-pa =PolyA annotations file standard 6 column bed format (string)

-apriori_annotations = Enable pre-loading of a priori PolyASite 2.0 and PolyADB 3.0 annotations (boolean toggle)

!Note: In general, between these -pa and -apriori_annotations options, use -apriori_annotations.

#### Tuning
-paired = Enable paired analyses where sample files are considered paired (i.e., pre-treatment vs post-treatment) for beta-binomial statistical test (boolean toggle)

## DEMO

Please refer to DEMO folder for PolyAMiner-Bulk demo command. 

```
python3 /mnt/belinda_local/venkata/data/PolyAMiner-Bulk/PolyA-miner.py -mode bam \
-fasta /mnt/belinda_local/venkata/data/Index_Files/Human/GenomeFasta_GTF/GRCh38.primary_assembly.genome.fa \
-gtf /mnt/belinda_local/venkata/data/Index_Files/Human/GenomeFasta_GTF/gencode.v33.primary_assembly.annotation.gtf \
-p 20 -a 0.65 -outPrefix 3UTROnly -expNovel 1 -s 2 \
-o /mnt/belinda_local/venkata/data/PolyAMiner-Bulk/Demo/Demo_Results/Demo_3UTROnly_Softclipped+APriori_ReRun \
-c1 /mnt/belinda_local/venkata/data/PolyAMiner-Bulk/Demo/control1_.subset.sorted.bam,\
/mnt/belinda_local/venkata/data/PolyAMiner-Bulk/Demo/control2_.subset.sorted.bam,\
/mnt/belinda_local/venkata/data/PolyAMiner-Bulk/Demo/control3_.subset.sorted.bam \
-c2 /mnt/belinda_local/venkata/data/PolyAMiner-Bulk/Demo/treatment1_.subset.sorted.bam,\
/mnt/belinda_local/venkata/data/PolyAMiner-Bulk/Demo/treatment2_.subset.sorted.bam,\
/mnt/belinda_local/venkata/data/PolyAMiner-Bulk/Demo/treatment3_.subset.sorted.bam \
-ignore UTR5,CDS,Intron,UN -apriori_annotations -modelOrganism human \
-visualizeTopNum 10 -visualizeCondition1Name Control -visualizeCondition2Name Treatment
```
Important notes: 
(1) User will have to adjust the file location parameters to their specific file system


## FAQ

Q1) I am unable to install the gplots package and its corresponding dependencies when running "Rscript installPkgs.R". 
A1) This usually occurs when OS-level dependencies are missing. If on Linux, try running "apt-get install libblas-dev liblapack-dev prior to running "Rscript installPkgs.R". 

Q2) I am unable to install the DESEQ2 package and its corresponding dependencies when running "Rscript installPkgs.R". 
A2) This usually occurs when OS-level dependencies are missing. Look through the error messages. Usually the console will request specific OS-level dependencies to be manually installed (see answer to Q1)

## Feedback?

Please complete the following form if you have any questions or feedback: https://forms.gle/8rF4TZcPoS15PEsB8. We will update this readme with answers to the most frequently asked questions. 