# PolyAMiner-Bulk

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

git clone https://github.com/jerryji1993/DNABERT
cd DNABERT
python3 -m pip install --editable .
cd examples
python3 -m pip install -r requirements.txt
```

## Step 2: Download trained ML-models + ML-dependencies 

### 2.1: Download "6-new-12w-0" from https://drive.google.com/file/d/14bH1cxNX16bldIFIhkM33WFaNXp0XNPS/view?usp=sharing and extract it into the /lib folder

### 2.2: Download "CPASBERT_TrainedModels" from https://drive.google.com/file/d/1iMsBCbVSLPLjaHGHCiLMsxvqkzQgpk2Y/view?usp=sharing and extract it into the /lib folder

## Step 3: Download necessary reference files

### 3.1: Download organism-specific fasta and gtf reference files from Ensemble: http://uswest.ensembl.org/info/data/ftp/index.html/

## Step 4 (Optional): Download RBM17 Test Data and Run Test Script
