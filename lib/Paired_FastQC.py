#FastQC Script for NGS Pipeline

#sys.argv[1] = input directory where fastq.gz files are in the main directory (fastq_data/)
inputDir="/mnt/local_storage/venkata/data/ALS/RAW_Data/GSE74973"
#sys.argv[2] = output directory where FastQC reports for each sample and general MultiQC will be stored (FastQC/)
outputDir="/mnt/local_storage/venkata/data/ALS/GSE74973/FASTQC_Results"

import os,sys, glob
  

# print("Started FastQC analysis")
sys.argv[1]=sys.argv[1].rstrip("/")+"/"
sys.argv[2]=sys.argv[2].rstrip("/")+"/"

#grabbing paths of fastq.gz files
files=glob.glob(sys.argv[1]+"*.fastq.gz")
# print(files)

#empty list for processed samples
processed=[]

#creating directory for FastQC and MultiQC reports
os.system("mkdir "+sys.argv[2]+"FastQC")

file = files[0]
#sample name
sample=file.split("/")[-1].split("_")[0]
#sample="_".join(file.split("/")[-1].split("_")[:8])
#skip sample if already processed
if sample in processed:
    pass
#processing sample
else:
    processed.append(sample)
    print("\nProcessing:"+sample+"\n")
    #creating directory for FastQC reports for each sample
    os.system("mkdir "+sys.argv[2]+"FastQC/"+sample)
    #grabbing paths of R1/R2 files and sorting by ascending lanes
    r1=glob.glob(sys.argv[1]+sample+"_*R1*.fastq.gz")
    r2=glob.glob(sys.argv[1]+sample+"_*R2*.fastq.gz")
    print(r1)
    ##r1=glob.glob(sys.argv[1]+"*/"+sample+"_*R1*.fastq.gz")
    ##r2=glob.glob(sys.argv[1]+"*/"+sample+"_*R2*.fastq.gz")
    r1.sort();r2.sort();
    #multiple lanes for R1/R2
    if len(r1) >1:
        #merging lanes for R1
        cmd = "cat "+" ".join(r1)+" > "+sys.argv[2]+"FastQC/"+sample+"/"+sample+"_R1.fastq.gz"
        # print(cmd+"\n")
        os.system(cmd)
        #merging lanes for R2
        cmd = "cat "+" ".join(r2)+" > "+sys.argv[2]+"FastQC/"+sample+"/"+sample+"_R2.fastq.gz"
        # print(cmd+"\n")
        os.system(cmd)
       #creating FastQC reports
        cmd='fastqc -t 2 ' + sys.argv[2]+"FastQC/"+sample+"/"+sample+"_R1.fastq.gz "+sys.argv[2]+"FastQC/"+sample+"/"+sample+"_R2.fastq.gz -o "+sys.argv[2]+"FastQC/"+sample+"/"
        # print(cmd+"\n")
        os.system(cmd)
    #single lane for R1/R2    
    else:
        #grabbing paths for R1/R2 files
        r1=r1[0]
        r2=r2[0]
        #creating FastQC reports
        cmd='fastqc -t 6 ' +r1+" "+r2+" -o "+sys.argv[2]+"FastQC/"+sample
        # print(cmd+"\n")
        os.system(cmd)
