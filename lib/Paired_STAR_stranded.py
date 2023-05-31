#STAR Alignment Script for NGS Pipeline (Michelle)
  
#sys.argv[1] = input directory where trimmed fastq.gz files are in the Fastp directory (Fastp/Trim_#)
#sys.argv[2] = output directory where STAR alignment output will be stored in individual sample folders (STAR)
#sys.argv[3] = STAR genome directory for read length 
#sys.argv[4] = number of cores to use (usually 10)
#sys.argv[5] = overhang = readlength-1
##REMEMBER to change ohang (readLength-1) and ank to appropriate values 

import os, sys, glob

files=glob.glob(sys.argv[1]+"*.fastq.gz")
print(files)

processed=[]

#Generating STAR alignments
for file in files:

    sample=file.split("/")[-1].split("_")[0]
    #sample="_".join(file.split("/")[-1].split("_")[:2])

    
    if sample in processed:
        pass

    else:
        processed.append(sample)
        print("\nProcessing:"+sample+"\n")

        os.system("mkdir "+sys.argv[2]+sample)
    
        r1=glob.glob(sys.argv[1]+sample+"*R1*.fastq.gz")
        r2=glob.glob(sys.argv[1]+sample+"*R2*.fastq.gz")
        r1.sort();r2.sort();
        
        #multiple lanes for R1/R2
        if len(r1) >1:
            r1 = ",".join(r1)
            r2 = ",".join(r2)

            
        #single lane for R1/R2    
        else:
            r1=r1[0]
            r2=r2[0]

        #Alignment using STAR 
        genomeDir = sys.argv[3]
        ohang = sys.argv[5];
        # --outFilterScoreMinOverLread 0.2 --outFilterMatchNminOverLread 0.2
        cmd='STAR --runThreadN '+sys.argv[4]+' --genomeDir '+ genomeDir + ' --quantMode GeneCounts --outSJfilterOverhangMin 8 8 8 8 --outSJfilterCountUniqueMin 1 1 1 1 --outSJfilterCountTotalMin 1 1 1 1 --outSJfilterDistToOtherSJmin 10 0 5 10 --outSJfilterIntronMaxVsReadN 500000 500000 500000 --outSAMtype BAM SortedByCoordinate --sjdbOverhang '+ohang+' --outFileNamePrefix '+sys.argv[2]+sample+"/"+sample+'_ --readFilesCommand gunzip -c --readFilesIn '+r1+" "+r2
        print(cmd+"\n")
        os.system(cmd)

        # Command for Indexing  
#        cmd="samtools index "+sys.argv[2]+"*"+sample+"*/*.bam"
#        print(cmd+"\n")
#        os.system(cmd)

        #Command for renaming bam and bam index #
#        f=glob.glob(sys.argv[2]+"*"+sample+"*/*.bam")[0]
#        cmd="mv "+f+" "+f.replace("Aligned.sortedByCoord.out.bam",".sorted.bam")
#        print(cmd+"\n")
#        os.system(cmd)
#        f=glob.glob(sys.argv[2]+"*"+sample+"*/*.bam.bai")[0]
#        cmd="mv "+f+" "+f.replace("Aligned.sortedByCoord.out.bam.bai",".sorted.bam.bai")
#        print(cmd+"\n")
#        os.system(cmd)

        # Command for renaming counts file #
#        f=glob.glob(sys.argv[2]+"*"+sample+"*/*ReadsPerGene.out.tab")[0]
#        cmd="mv "+f+" "+f.replace("ReadsPerGene.out.tab",".counts.txt")
#        print(cmd+"\n")
#        os.system(cmd)

        print("Finished processing:"+sample)

