# PolyAMiner-Bulk

This repository includes the PolyAMiner-Bulk computational tool from 'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'. Please cite our paper if you have used our computational tool, any of our machine learning models, or code snippets. Additionally, this repository is actively under development, so please kindly report any issues or feature requests.

In this package, we provide the following resources: 

(1) Source code of PolyAMiner-Bulk

(2) Test scripts delineating key usage scenarios

(3) Fine-tuned CPAS-BERT models for both human and mouse model organisms.

## Citation

If you have used PolyAMiner-Bulk in your research, please kindly cite the following publication:

XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


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

### Cleaning

#### Optional Parameters
optional.add_argument('-ignore',help=' Comma-separated list of regions to igonore from APA analysis UTR5, Introns, CDS, UTR3', nargs='+',type=str,default='ZZZZZZ')
optional.add_argument('-apaBlock',help='Window size for annotated polyA sites',type=int, default=30)
optional.add_argument('-mdapa',help='Cluster distance for annotated polyA sites: Merge polyA sites with in this distance. ',type=int, default=0)
optional.add_argument('-md',help='Cluster distance for de-novo polyA sites: Merge polyA sites with in this distance',type=int, default=0)
optional.add_argument('-anchor',help='Overlap in "bp" for mapping de-novo polyA sites to annotated polyA sites ',type=int, default=1)
optional.add_argument('-softclippedNumReads',help='(Param 1 of 2 for Softclipped-Assisted Filter) Minimum # of required softclipped reads',type=int, default=0)
optional.add_argument('-softclippedNumSamples',help='(Param 2 of 2 for Softclipped-Assisted Filter) Minimum # of samples that need to meet softclippedNumReads requirement',type=int, default=0)
optional.add_argument('-slopDistanceParameter',help='Slop distance parameter',type=int, default = 25)
optional.add_argument('-clusterParameter',help='Cluster Parameter',type=int, default = 30)

#### Additional analysis 
optional.add_argument('-cluster_onGenes',help='Cluster samples based on gene counts - all polyA sites / 3utr polyA sites',choices=['all','3utr','none'],type=str,default='none')
optional.add_argument('-cluster_onPAsites',help='Cluster samples based on gene counts - all polyA sites / 3utr polyA sites',choices=['all','3utr','none'],type=str,default='none')
optional.add_argument('-DEG',help='Perform differential gene expression analysis - all polyA sites / 3utr polyA sites',choices=['all','3utr','none'],type=str,default='none')
optional.add_argument('-pa_usage',help='Perform differential polyA site usage analysis - all polyA sites / filtered polyA sites',choices=['all','filtered','none'],type=str,default='none')

#### Tuning 
optional.add_argument('-expNovel',help='Explore novel APA sites 0: only annotated sites 1: de-novo',choices=[1,0],type=int,default=0)
optional.add_argument('-novel_d',help='Distance from annotated TES to map novel pA sites',type=int, default=1000)
optional.add_argument('-p',help='No. of processors to use',type=int,default=4)
optional.add_argument('-ip_d',help='Downstream internal priming window',type=int, default=50)
optional.add_argument('-ip_u',help='Upstream internal priming window',type=int, default=50)
optional.add_argument('-a',help='Internal priming polyA fraction',type=float, default=0.65)
optional.add_argument('-pa_p',help='pOverA filter: P ',type=float, default=0.6)
optional.add_argument('-pa_a',help='pOverA filter: A ',type=int, default=5)
optional.add_argument('-pa_m',help='pOverA filter: M ',type=int, default=2)
optional.add_argument('-gene_min',help='Min counts per Gene',type=int, default=10)
optional.add_argument('-apa_min',help='Min. proportion per APA',type=float, default=0.05)
optional.add_argument('-batchCorrection',help='Comma-seperated list of batch membership for ComBat-seq correction',type=str, default='ZZZZZZ')
optional.add_argument('-modelOrganism',help='Model organism for CPAS-BERT Model',type=str, default="human")
optional.add_argument('-visualizeTopNum',help='Generate read density visualization plots for the top N genes, where N is an integer specified by the user ',type=int, default=0)
optional.add_argument('-visualizeCondition1Name',help='Name of Condition 1 when generating read density visualization plots',type=str, default="Control")
optional.add_argument('-visualizeCondition2Name',help='Name of Condition 2 when generating read density visualization plots',type=str, default="Treatment")
optional.add_argument('-visualizeCondition1NameHeatmap',help='Name of Condition 2 when generating read density visualization plots',type=str, default="CR")
optional.add_argument('-visualizeCondition2NameHeatmap',help='Name of Condition 2 when generating read density visualization plots',type=str, default="TR")
optional.add_argument('-outPrefix',help='Output file/s prefix', default="PolyAminer_Out",type=str)

## Step 4 (Optional): Run Test Script

## FAQ

Q1) I am unable to install the gplots package and its corresponding dependencies when running "Rscript installPkgs.R". 
A1) This usually occurs when OS-level dependencies are missing. If on Linux, try running "apt-get install libblas-dev liblapack-dev prior to running "Rscript installPkgs.R". 

Q2) I am unable to install the DESEQ2 package and its corresponding dependencies when running "Rscript installPkgs.R". 
A2) This usually occurs when OS-level dependencies are missing. Look through the error messages. Usually the console will request specific OS-level dependencies to be manually installed (see answer to Q1)

## Feedback?

Please complete the following form if you have any questions or feedback: https://forms.gle/8rF4TZcPoS15PEsB8. We will update this readme with answers to the most frequently asked questions. 