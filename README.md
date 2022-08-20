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
conda create -n cpasbert python=3.6
conda activate cpasbert
```

### 1.2: Install the necessary package dependencies

```
conda install pytorch torchvision cudatoolkit=10.0 -c pytorch

git clone https://github.com/venkatajonnakuti/PolyAMiner-Bulk
cd PolyAMiner-Bulk
python3 -m pip install -r requirements.txt
```

## Step 2: Download trained ML-models + ML-dependencies 

### 2.1: Download "6-new-12w-0" from https://bcm.box.com/s/1lo85ig6eop2a5e12d52ytd6q98fg4xa and extract it into the /lib folder

### 2.2: Download "CPASBERT_TrainedModels" from https://drive.google.com/file/d/1iMsBCbVSLPLjaHGHCiLMsxvqkzQgpk2Y/view?usp=sharing and extract it into the /lib folder

## Step 3: Download necessary reference files

### 3.1: Download organism-specific fasta and gtf reference files from Ensemble: http://uswest.ensembl.org/info/data/ftp/index.html/

## Step 4: Run PolyAMiner-Bulk

## Step 5 (Optional): Download RBM17 Test Data and Run Test Script

## Q&A
