import os, sys, glob 
import pandas as pd 
import pysam
from CPASBERT import CPASBERT
from pybedtools import BedTool 
import gtfparse as gp   
import pybedtools as pb 

class ExtractPolyAsites4Bulk:
	def __init__(self, outDir, outPrefix, fasta, gtf, con1BAMFiles, con2BAMFiles, proportionA, modelOrganism, apriori_annotations, ignoreFeatures, clusterParameter, slopDistanceParameter, softclip_numSamples, softclip_numReads):
		#Housekeeping variables
		self.outDir = outDir.rstrip("/")+"/"
		self.outPrefix = outPrefix
		self.con1BAMFiles = "".join(con1BAMFiles).replace(" ","").split(",")
		self.con2BAMFiles = "".join(con2BAMFiles).replace(" ","").split(",")
		self.FASTA = fasta
		self.GTF = gtf
		
		#Variables that are transformations of inputs
		self.conditionBAMList = self.con1BAMFiles + self.con2BAMFiles

		#Tuning variables
		self.proportionA = proportionA #Proportion of A's in true polyA tails. Default: '0.90,0.85,0.80,0.75' ->0.90 (for a tail of 4nt long), 0.85 (for a tail of 4-8nt long), 0.80 (for a tail of 8-12nt long), 0.75 (for a tail > 12nt long) ''',nargs='+',type=str,default= "0.90,0.85,0.80,0.75"
		self.clusterParameter = clusterParameter
		self.slopDistanceParameter = slopDistanceParameter

		#Mode variables
		self.apriori_annotations = apriori_annotations
		self.softclip_numReads = int(softclip_numReads)
		self.softclip_numSamples = int(softclip_numSamples)
		self.ignoreFeatures = ignoreFeatures
		libPath = os.path.dirname(os.path.abspath(__file__)) + "/CPASBERT_TrainedModels"
		if modelOrganism == "mouse":
			self.modelPath = libPath + "/mm10_checkpoint-8000"
			self.POLYASITE = libPath + "/PolyASite_mouse_mm10.bed"
			self.POLYADB = libPath + "/PolyADB_mouse_mm10.bed"
			self.APRIORI = libPath + "/APrioriAnnotations_PolyADB_PolyASite_mouse_mm10.bed"

		elif modelOrganism == "human":
			self.modelPath = libPath + "/hg38_checkpoint-64000"
			self.POLYASITE = libPath + "/PolyASite_human_hg38.bed"
			self.POLYADB = libPath + "/PolyADB_human_hg38.bed"
			self.APRIORI = libPath + "/APrioriAnnotations_PolyADB_PolyASite_human_hg38.bed"

		#Global variables that are called between functions
		self.CONCATENATED_SOFTCLIPPED_BED_FILE = self.outDir + self.outPrefix + ".concatenated.softclipped.bed"
		self.CPAS_BED_FILE = self.outDir + self.outPrefix + ".CPAS.bed"
		self.GENE_BED = self.outDir + self.outPrefix + "Genes.bed"
		self.CDS_BED = self.outDir + self.outPrefix + "CDS.bed"
		self.UTR3_BED = self.outDir + self.outPrefix + "UTR3.bed"
		self.UTR5_BED = self.outDir + self.outPrefix + "UTR5.bed"
		self.INTRON_BED = self.outDir + self.outPrefix + "Intron.bed"
		self.FEATURES_BED = self.outDir + self.outPrefix + "Features.bed"
		self.CPAS_SAF_FILE = self.outDir + self.outPrefix + "_denovoAPAsites.saf"

	def _compute(self, string, key):
		string = string.strip().upper()
		return(string.count(key)/float(len(string)))

	# Extract soft-clipped (or polyA-capped) reads from each BAM file and save them as separate BAM and BED files	
	def _extractSoftClippedReads(self):
		prop = "".join(self.proportionA).replace(" ","").split(",")

		for BAM_FILE in self.conditionBAMList:
			CAPPED_BAM_FILE = self.outDir + self.outPrefix + os.path.basename(BAM_FILE).replace(".bam", ".capped.bam")
			CAPPED_BED_FILE = self.outDir + self.outPrefix + os.path.basename(BAM_FILE).replace(".bam", ".capped.bed")

			BED_object = open(CAPPED_BED_FILE, "w")
			bamfile = pysam.AlignmentFile(BAM_FILE, 'rb')
			polyAreads = pysam.AlignmentFile(CAPPED_BAM_FILE, "wb", template=bamfile)
			
			for read in bamfile:
				if (not read.is_unmapped) and ("S" in read.cigarstring):
					if (not read.is_reverse):
						if read.cigartuples[-1][0]== 4 and read.cigartuples[-1][1] >=4:
							count = self._compute(read.query_sequence[(-1*read.cigartuples[-1][1]):],"A")
							tail=read.cigartuples[-1][1]
							BED_Entry = read.reference_name + "\t" + str(read.reference_end - 1) + "\t" + str(read.reference_end + 1) + "\t\t\t+\n"		
					if (read.is_reverse):
						if read.cigartuples[0][0]== 4 and read.cigartuples[0][1] >=4:
							count = self._compute(read.query_sequence[:read.cigartuples[0][1]],"T")
							tail=read.cigartuples[0][1]
							BED_Entry = read.reference_name + "\t" + str(read.reference_start - 1) + "\t" + str(read.reference_start + 1) + "\t\t\t-\n"
					try:
						if tail > 0 and tail <=4 and count >= float(prop[0]):
							polyAreads.write(read)
							BED_object.write(BED_Entry)
						elif tail > 4 and tail <=8 and count >= float(prop[1]):
							polyAreads.write(read)
							BED_object.write(BED_Entry)
						elif tail > 8 and tail <=12 and count >= float(prop[2]):
							polyAreads.write(read)
							BED_object.write(BED_Entry)
						elif tail > 12 and count >= float(prop[3]):
							polyAreads.write(read)
							BED_object.write(BED_Entry)
						else:
							pass
					except:
						pass

			polyAreads.close()
			BED_object.close()
			os.system("samtools index " + CAPPED_BAM_FILE)

		CAPPED_BED_FILELIST = glob.glob(self.outDir + "*.capped.bed")
		CAPPED_BED_FILELIST_STRING = ' '.join(CAPPED_BED_FILELIST)
		cmd = "cat " + CAPPED_BED_FILELIST_STRING + " > " + self.CONCATENATED_SOFTCLIPPED_BED_FILE
		# print(cmd)
		os.system(cmd)
		BedTool(self.CONCATENATED_SOFTCLIPPED_BED_FILE).saveas(self.outDir + self.outPrefix + ".concatenated.softclipped_beforeCPASBERTFiltering.bed")

		#Remove duplicate BED entries
		TMP_BED_FILELOC = self.outDir + self.outPrefix + ".tmp.bed"
		cmd = "rm -f " + TMP_BED_FILELOC
		# print(cmd)
		os.system(cmd)
		cmd = "awk \'!x[$0]++\' " + self.CONCATENATED_SOFTCLIPPED_BED_FILE + " > " + TMP_BED_FILELOC
		# print(cmd)
		os.system(cmd)
		cmd = "mv " + TMP_BED_FILELOC + " " + self.CONCATENATED_SOFTCLIPPED_BED_FILE
		# print(cmd)
		os.system(cmd)
		BedTool(self.CONCATENATED_SOFTCLIPPED_BED_FILE).saveas(self.outDir + self.outPrefix + ".concatenated.softclipped_beforeCPASBERTFilteringAndAfterDeDup.bed")

	# Append filtered APriori Annotations to BED File containing filtered soft-clipped reads
	def _loadFilteredAPrioriAnnotations(self):
		APRIORI_DF = pd.read_csv(self.APRIORI, sep = "\t", header = None, names = ["Chr", "Start", "End", "Feature1", "Feature2", "Strand"])
		APRIORI_DF[['Start', 'End']] = APRIORI_DF[['Start', 'End']].astype("int")
		APRIORI_DF["Midpoint"] = APRIORI_DF[['Start', 'End']].mean(axis=1)
		APRIORI_DF['Midpoint'] = APRIORI_DF['Midpoint'].apply(lambda x: round(x, 0))
		APRIORI_DF = APRIORI_DF.astype({"Midpoint": int})
		APRIORI_DF["Start"] = APRIORI_DF["Midpoint"] - 1
		APRIORI_DF["End"] = APRIORI_DF["Midpoint"] + 1
		APRIORI_DF = APRIORI_DF.drop('Midpoint', axis=1)
		modifiedAPrioriFileLoc = self.outDir + self.outPrefix + ".modified.aprioriAnnotations.bed"
		APRIORI_DF.to_csv(modifiedAPrioriFileLoc, sep='\t', index=False, header=None)

		cmd = "cat " + modifiedAPrioriFileLoc + " " + self.CONCATENATED_SOFTCLIPPED_BED_FILE + " > " + self.CPAS_BED_FILE
		# print(cmd)
		os.system(cmd)

	def _performSoftClippedAssistedFiltering(self):
		if not self.apriori_annotations:
			self.softclip_numSamples = 1
			self.softclip_numReads = 1 

		cappedCountsFileList = []
		KEEP_CLUSTERID_LIST = []
		TMP_BED_FILELOC = self.outDir + self.outPrefix + ".tmp.bed"
		
		#Remove duplicate BED entries
		cmd = "awk \'!x[$0]++\' " + self.CPAS_BED_FILE + " > " + TMP_BED_FILELOC
		# print(cmd)
		os.system(cmd)
		cmd = "mv " + TMP_BED_FILELOC + " " + self.CPAS_BED_FILE
		# print(cmd)
		os.system(cmd)
		
		#Import BED file as a BedTool
		inputBEDObject = BedTool(self.CPAS_BED_FILE)

		#Perform clustering
		CLUSTERED_BED_FILELOC = self.outDir + self.outPrefix + "CPASwithClusteredColumn.BED"
		clusteredBEDObject = inputBEDObject.sort().cluster(s = True, d = self.clusterParameter)
		clusteredBEDObject.saveas(CLUSTERED_BED_FILELOC)
		# df = pd.read_csv(CLUSTERED_BED_FILELOC, sep = "\t", header = None)
		# df[2] = df[1] - 5
		# df[2] = df[2] + 5
		# df.to_csv(CLUSTERED_BED_FILELOC, sep = "\t", index = False, header = None)
		clusteredBEDObject = BedTool(CLUSTERED_BED_FILELOC)

		##KEEP ONLY CLUSTERS SUPPORTED BY A SOFT-CLIPPED READ
		for BAM_FILE in self.conditionBAMList:
			CAPPED_BED_FILE = self.outDir + self.outPrefix + os.path.basename(BAM_FILE).replace(".bam", ".capped.bed")

			cmd = "rm -f " + TMP_BED_FILELOC
			# print(cmd)
			os.system(cmd)

			df = pd.read_csv(CAPPED_BED_FILE, sep = "\t", header = None)
			df = df[df[1] != -1]
			df = df[df[2] != -1]
			df.to_csv(CAPPED_BED_FILE, sep = "\t", header = None, index = None)

			cappedBEDObject = BedTool(CAPPED_BED_FILE)
			
			tmpBEDObject = clusteredBEDObject.intersect(cappedBEDObject, c = True, s = True, a = clusteredBEDObject, b = cappedBEDObject, output = TMP_BED_FILELOC)
			BedTool(TMP_BED_FILELOC).saveas(self.outDir + self.outPrefix + "Intersect." + str(self.conditionBAMList.index(BAM_FILE)) + ".TMP.bed")
			# cmd = "bedtools intersect -wao -s -a " + CLUSTERED_BED_FILELOC + " -b " + CAPPED_BED_FILE + " > " + TMP_BED_FILELOC
			# print(cmd)
			# os.system(cmd)

			df = pd.read_csv(TMP_BED_FILELOC, sep = "\t", header = None)
			df["Cluster"] = df[6]
			cols = [0, 1, 2, 5]
			df[os.path.basename(BAM_FILE).replace(".bam", ".softclipped_ct")] = df[7]
			df['BED_ID'] = df[cols].apply(lambda row: '_'.join(row.values.astype(str)), axis=1)
			# df_drop = df[df[7] != -1]

			# df2 = df_drop.groupby("BED_ID").count()[[0]]
			# df2.columns = [os.path.basename(BAM_FILE).replace(".bam", ".softclipped_ct")]
			# df = pd.merge(df, df2, on='BED_ID', how = "outer")
			# df[os.path.basename(BAM_FILE).replace(".bam", ".softclipped_ct")] = df[os.path.basename(BAM_FILE).replace(".bam", ".softclipped_ct")].fillna(0)
			# df[os.path.basename(BAM_FILE).replace(".bam", ".softclipped_ct")] = pd.to_numeric(df[os.path.basename(BAM_FILE).replace(".bam", ".softclipped_ct")], downcast='integer')
			# cols = ["Cluster", "BED_ID", os.path.basename(BAM_FILE).replace(".bam", ".softclipped_ct")]
			# df = df[cols]
			df = df[["BED_ID", "Cluster", os.path.basename(BAM_FILE).replace(".bam", ".softclipped_ct")]]
			
			CAPPED_COUNTS_FILE = self.outDir + self.outPrefix + os.path.basename(BAM_FILE).replace(".bam", ".softclipped.counts")
			cappedCountsFileList.append(CAPPED_COUNTS_FILE)
			df.to_csv(CAPPED_COUNTS_FILE, sep = "\t", index = None)

			cmd = "awk \'!x[$0]++\' " + CAPPED_COUNTS_FILE + " > " + TMP_BED_FILELOC
			# print(cmd)
			os.system(cmd)
			cmd = "mv " + TMP_BED_FILELOC + " " + CAPPED_COUNTS_FILE
			# print(cmd)
			os.system(cmd)

		cappedCountsFileMergedDataframe = pd.read_csv(cappedCountsFileList[0], sep = "\t")

		for i in range(1,len(cappedCountsFileList)):
			data = pd.read_csv(cappedCountsFileList[i], sep = "\t")
			df_inter = pd.DataFrame(data)
			df_inter = df_inter.drop(['Cluster'], axis=1)
			cappedCountsFileMergedDataframe = pd.merge(cappedCountsFileMergedDataframe, df_inter, on='BED_ID')

		cappedCountsFileMergedDataframe.to_csv(TMP_BED_FILELOC, sep = "\t", index = None)
		#####Change bottom value to 2 and not 1!!
		final_dataframe = cappedCountsFileMergedDataframe.iloc[:, 1:]
		final_dataframe["Cluster"] = cappedCountsFileMergedDataframe["Cluster"]

		SUMMARIZED_INTER_SOFTCLIPPED_COUNTS = self.outDir + self.outPrefix + ".summarized.softclipped.counts"
		final_dataframe.to_csv(SUMMARIZED_INTER_SOFTCLIPPED_COUNTS, sep = "\t", index = False)
		final_dataframe = final_dataframe.set_index('Cluster')

		for index, row in final_dataframe.iterrows():
			samples = [num for num in row[1:].values if num >= self.softclip_numReads]
			if (len(samples) >= self.softclip_numSamples):
				KEEP_CLUSTERID_LIST.append(index)

		CPAS_df = pd.read_csv(CLUSTERED_BED_FILELOC, sep = "\t", header = None, names = ["Chr", "Start", "End", "Feature1", "Feature2", "Strand", "Cluster"])
		CPAS_df = CPAS_df[CPAS_df.Cluster.isin(KEEP_CLUSTERID_LIST)]
		CPAS_df.to_csv(self.CPAS_BED_FILE, sep = "\t", index = None, header = None)
		BedTool(self.CPAS_BED_FILE).saveas(self.outDir + self.outPrefix + ".concatenated.softclipped_afterSoftClippedAssistedFilteringModule.bed")

		# cmd = "rm " + TMP_BED_FILELOC
		# print(cmd)
		# os.system(cmd)



		# cmd = "mv " + TMP_BED_FILELOC + " " + self.CPAS_BED_FILE
		# print(cmd)
		# os.system(cmd)


	# Given a BED file, keep only the most proximal read per cluster. Then, extend (or slop) the read by a user-determined distance on both ends. 
	def _transformCPASBED(self):
		#Keep only the most proximal read per cluster
		CPAS_df = pd.read_csv(self.CPAS_BED_FILE, sep = "\t", header = None, names = ["Chr", "Start", "End", "Feature", "Feature2", "Strand", "Cluster"])
		#Below is a modification that I made on 10.14.22
		CPAS_df = CPAS_df[~CPAS_df.Chr.str.contains(" ")].reset_index(drop=True)
		Filtered_DF_posStrand = CPAS_df.query('Strand == "+"').reset_index(drop=True)
		Filtered_DF_negStrand = CPAS_df.query('Strand == "-"').reset_index(drop=True)
		Filtered_DF_posStrand = Filtered_DF_posStrand.loc[Filtered_DF_posStrand.groupby('Cluster').Start.idxmin()]
		Filtered_DF_negStrand = Filtered_DF_negStrand.loc[Filtered_DF_negStrand.groupby('Cluster').End.idxmax()]
		pieces = (Filtered_DF_posStrand, Filtered_DF_negStrand)
		Filtered_DF = pd.concat(pieces, ignore_index = True)
		Filtered_DF = Filtered_DF.drop(["Cluster"], axis=1)
		Filtered_DF.to_csv(self.CPAS_BED_FILE, sep = "\t", index = None, header = None)
		BedTool(self.CPAS_BED_FILE).saveas(self.outDir + self.outPrefix + ".concatenated.softclipped_inTransformCPASBEDModule_AfterChoosingOneReadPerCluster.bed")

		#Slop the BED file by a user-determined distance on both ends
		cmd = "samtools faidx " + self.FASTA
		# print(cmd)
		os.system(cmd)

		FASTA_INDEX = self.FASTA.replace(".fa", ".fa.fai")
		CHROM_SIZES = self.outDir + self.outPrefix + "chrom.sizes"
		cmd = "cut -f 1,2 " + FASTA_INDEX + " > " + CHROM_SIZES
		# print(cmd)
		os.system(cmd)

		BedTool(self.CPAS_BED_FILE).slop(b = self.slopDistanceParameter, g = CHROM_SIZES).sort().saveas(self.CPAS_BED_FILE)
		BedTool(self.CPAS_BED_FILE).saveas(self.outDir + self.outPrefix + ".concatenated.softclipped_inTransformCPASBEDModule_AfterSlopping.bed")

	def _makeGeneBed(self):
		df = gp.read_gtf(self.GTF)
		df = df[df['seqname'].astype(str).str.contains('chr')]
		
		# Gene table in bed format #
		df_gene = df[df['feature'] == 'gene']
		df_gene = df_gene[['seqname', 'start', 'end', 'gene_id', 'gene_name', 'strand']]
		df_gene = df_gene.sort_values(["seqname", "start", "end"], ascending = (True, True,True))
		
		df_gene.to_csv(self.GENE_BED, sep='\t', index=False, header=False)

		return None

	def _addGeneNames(self):
		self._makeGeneBed()

		TEMPFILE = self.outDir + self.outPrefix + "tempFile.bed"
		cmd = 'bedtools closest -nonamecheck -a ' + self.CPAS_BED_FILE + ' -b ' + self.GENE_BED + ' -s -id -D a -t first -k 1 > ' + TEMPFILE
		os.system(cmd)
		tempdf = pd.read_csv(TEMPFILE, sep='\t', index_col=None, header=None)
		tempdf = tempdf.iloc[:, [0, 1, 2, 5, 10, 12]]
		tempdf.columns = ['Chr', 'Start', 'End', 'Strand', 'Gene', 'Distance']
		tempdf = tempdf[['Chr', 'Start', 'End', 'Gene', 'Distance', 'Strand']]
		# os.system('rm ' + TEMPFILE)
		tempdf.to_csv(self.CPAS_BED_FILE, sep='\t', index=False, header=None)

	def _mapAPA2Features(self):
		APAfile = self.CPAS_BED_FILE

		PA_df = pd.read_csv(APAfile,sep="\t",header=None,index_col=None)
		PA_df.columns = ['Chr', 'Start', 'End', 'GeneID', 'feature', 'Strand']
		PA = pb.BedTool.from_dataframe(PA_df)

		FT = pb.BedTool(self.FEATURES_BED)
		PA_FT_bed = PA.intersect(FT, nonamecheck=True, s=True, wao=True)
		PA_FT_df = PA_FT_bed.to_dataframe(disable_auto_names=True, header=None)
		PA_FT_df.columns = ['Chr', 'Start', 'End', 'GeneID', 'APAID', 'Strand','Chr_2', 'Start_2', 'End_2', 'GeneID_2', 'APAID_2', 'Strand_2','overlap']	
		
		PA_FT_df_nomap = PA_FT_df[PA_FT_df['Chr_2']=="."].copy()
		PA_FT_df_nomap = PA_FT_df_nomap[['Chr', 'Start', 'End', 'GeneID', 'APAID', 'Strand']]
		PA_FT_df_nomap['APAID'] = "UN"
			
		PA_nomap_FT = pb.BedTool.from_dataframe(PA_FT_df_nomap).intersect(FT, nonamecheck=True, s=False, wao=True)
		PA_nomap_FT_df = PA_nomap_FT.to_dataframe(disable_auto_names=True, header=None)
		PA_nomap_FT_df.columns = ['Chr', 'Start', 'End', 'GeneID', 'APAID', 'Strand','Chr_2', 'Start_2', 'End_2', 'GeneID_2', 'APAID_2', 'Strand_2','overlap']
		UTR3_EX_df = PA_nomap_FT_df[PA_nomap_FT_df['Start_2']==-1].copy()
		UTR3_EX_df = UTR3_EX_df[['Chr', 'Start', 'End', 'GeneID', 'APAID', 'Strand']]
		UTR3_EX_df['APAID'] = "UTR3-EX"

		PA_FT_df = PA_FT_df[(((PA_FT_df['Strand'] =="+") & (PA_FT_df['APAID_2'] !="UTR3")) & (PA_FT_df['End'] >= PA_FT_df['Start_2']) & (PA_FT_df['End'] <= PA_FT_df['End_2'])) | (((PA_FT_df['Strand'] =="+") & (PA_FT_df['APAID_2'] =="UTR3")) & (PA_FT_df['End'] >= PA_FT_df['Start_2'])) | (((PA_FT_df['Strand'] =="-") & (PA_FT_df['APAID_2'] !="UTR3")) & (PA_FT_df['Start'] <= PA_FT_df['End_2']) & (PA_FT_df['Start'] >= PA_FT_df['Start_2'])) | (((PA_FT_df['Strand'] =="-") & (PA_FT_df['APAID_2'] =="UTR3")) & (PA_FT_df['Start'] <= PA_FT_df['End_2']))]
		PA_FT_df['APAID'] = PA_FT_df['APAID_2']
		PA_FT_df = PA_FT_df[['Chr', 'Start', 'End', 'GeneID', 'APAID', 'Strand']]
		PA_FT_bed = pb.BedTool.from_dataframe(PA_FT_df)
		PA_FT_bed = PA_FT_bed.sort()
		PA_FT_bed = PA_FT_bed.merge(c=[4,5,6],s=True,o="distinct")
		PA_FT_df = PA_FT_bed.to_dataframe(disable_auto_names=True, header=None)
		PA_FT_df.columns = ['Chr', 'Start', 'End', 'GeneID', 'APAID', 'Strand']
		
		PAmap_df = pd.concat([PA_FT_df,UTR3_EX_df], axis=0,ignore_index=True)
		PA_df = PA_df.merge(PAmap_df,on=['Chr', 'Start', 'End', 'GeneID', 'Strand'],how='left').fillna("UN")
		PA_df = PA_df.drop(columns=['feature'])
		PA_df = PA_df[['Chr', 'Start', 'End', 'GeneID', 'APAID', 'Strand']]
		PA_df.to_csv(APAfile, sep='\t', index=False, header=False)

	def _checkUTR (self, m):
		return_list = []

		for i in range (0, len(m)):
			status="UN"
			
			if m[i][5] =="+":
				if 	m[i][1] >= m[i][7]:
					status="3UTR"
				elif m[i][2] <= m[i][10]:
					status="5UTR"
				else:
					pass
			
			if m[i][5] =="-":
				if m[i][2] <= m[i][8]:
					status="3UTR"
				elif m[i][1] >= m[i][10]:
					status="5UTR"
				else:
					pass
			
			return_list.append(status)

		return(return_list)

	def _getTables(self, df):
		# CDS/Coding Exon table #
		df_cds=df[df['feature']=='CDS']
		df_cds=df_cds[['seqname','start', 'end','gene_id','exon_id','strand']]
		df_cds=df_cds[df_cds['seqname'].str.contains("chr")]
		df_cds['exon_id']="CDS"
		df_cds=df_cds.sort_values(["seqname", "start","end"], ascending = (True, True,True))

		# For non zfish GTF -ENCODE #
		features_list=list(set(df['feature']))
		if "three_prime_utr" in features_list and "five_prime_utr" in features_list:
			df_utr5=df[df['feature']=='CDS']
			df_utr5=df_utr5[['seqname','start', 'end','gene_id','exon_id','strand']]
			df_utr5=df_utr5[df_utr5['seqname'].str.contains("chr")]
			df_utr5['transcript_id']="UTR5"
			df_utr5=df_utr5.sort_values(["seqname", "start","end"], ascending = (True, True,True))

			df_utr3=df[df['feature']=='CDS']
			df_utr3=df_utr3[['seqname','start', 'end','gene_id','exon_id','strand']]
			df_utr3=df_utr3[df_utr3['seqname'].str.contains("chr")]
			df_utr3['transcript_id']="UTR3"
			df_utr3=df_utr3.sort_values(["seqname", "start","end"], ascending = (True, True,True))
		
		else:
			# UTR table in bed format #
			df_utr=df[df['feature']=='UTR']
			df_utr=df_utr[['seqname','start', 'end','transcript_id','gene_id','strand']]
			df_utr=df_utr[df_utr['seqname'].str.contains("chr")]
			df_utr=df_utr.sort_values(["seqname", "start","end"], ascending = (True, True,True))

			# Stop table#
			df_stp=df[df['feature']=='stop_codon']
			df_stp=df_stp[['seqname','start', 'end','transcript_id','gene_id','strand']]
			df_stp=df_stp.sort_values(["seqname", "start","end"], ascending = (True, True,True))
			df_stp['transcripts_stop']=df_stp['transcript_id']+"_"+df_stp['start'].astype(str)+"_"+df_stp['end'].astype(str)
			df_stp_temp=df_stp[['seqname','start', 'end','transcript_id']]
			df_stp_temp.columns=['stp_seqname','stp_start', 'stp_end','transcript_id']
			df_utr=df_utr.merge(df_stp_temp,on=['transcript_id'],how='inner')
			
			# Start table#
			df_str=df[df['feature']=='start_codon']
			df_str=df_str[['seqname','start', 'end','transcript_id','gene_id','strand']]
			df_str=df_str.sort_values(["seqname", "start","end"], ascending = (True, True,True))
			df_str['transcripts_start']=df_str['transcript_id']+"_"+df_str['start'].astype(str)+"_"+df_str['end'].astype(str)
			df_str_temp=df_str[['seqname','start', 'end','transcript_id']]
			df_str_temp.columns=['str_seqname','str_start', 'str_end','transcript_id']
			df_utr=df_utr.merge(df_str_temp,on=['transcript_id'],how='inner')


			# Mark 5'UTR and 3'UTR #
			df_utr['status']=self._checkUTR(df_utr.values)
			df_utr5=df_utr[df_utr['status']=="5UTR"]
			df_utr5=df_utr5[['seqname','start','end','gene_id','transcript_id','strand']]
			df_utr5['transcript_id']='UTR5'
			df_utr3=df_utr[df_utr['status']=="3UTR"]
			df_utr3=df_utr3[['seqname','start','end','gene_id','transcript_id','strand']]
			df_utr3['transcript_id']='UTR3'

		return(df_cds,df_utr5,df_utr3)

	def _markFeatures(self):
		genes_bed = pb.BedTool(self.GENE_BED)

		df = gp.read_gtf(self.GTF)
		df = df[df['seqname'].astype(str).str.contains('chr')]
		cds_df, utr5_df, utr3_df = self._getTables(df)

		utr5_bed = pb.BedTool.from_dataframe(utr5_df)
		utr5_bed = utr5_bed.sort()
		utr5_merge = utr5_bed.merge(s=True, c=[4,5,6], o='distinct', d=1)
		utr5_df = utr5_merge.to_dataframe(disable_auto_names=True,header=None)
		utr5_df.columns=['chr','start','end','GeneID','Feature','strand']

		utr3_bed = pb.BedTool.from_dataframe(utr3_df)
		utr3_bed = utr3_bed.sort()
		utr3_merge = utr3_bed.merge(s=True, c=[4,5,6], o='distinct', d=1)
		utr3_df = utr3_merge.to_dataframe(disable_auto_names=True,header=None)
		utr3_df.columns = ['chr','start','end','GeneID','Feature','strand']

		cds_bed = pb.BedTool.from_dataframe(cds_df)
		cds_bed = cds_bed.sort()
		cds_merge = cds_bed.merge(s=True, c=[4,5,6], o='distinct', d=1)
		cds_merge = cds_merge.subtract(utr5_merge, s=True)
		cds_merge = cds_merge.subtract(utr3_merge, s=True)
		cds_df = cds_merge.to_dataframe(disable_auto_names=True,header=None)
		cds_df.columns=['chr','start','end','GeneID','Feature','strand']

		intr_bed = genes_bed.subtract(cds_merge, s=True)
		intr_bed = intr_bed.subtract(utr5_merge, s=True)
		intr_bed = intr_bed.subtract(utr3_merge, s=True)
		intr_df = intr_bed.to_dataframe(disable_auto_names=True,header=None)
		intr_df.columns=['chr','start','end','GeneID','Feature','strand']
		intr_df['Feature']='Intron'
		
		cds_df.to_csv(self.CDS_BED, sep="\t", header=False, index=False)
		utr3_df.to_csv(self.UTR3_BED, sep="\t", header=False, index=False)
		utr5_df.to_csv(self.UTR5_BED, sep="\t", header=False, index=False)
		intr_df.to_csv(self.INTRON_BED, sep="\t", header=False, index=False)

		feature_df = pd.concat([intr_df, cds_df, utr5_df, utr3_df], axis=0, ignore_index=True)
		feature_df.to_csv(self.FEATURES_BED, sep="\t", header=False, index=False)

		self._mapAPA2Features()
		
		return(1)

	def _annotateCPASBED(self):
		self._addGeneNames()
		self._markFeatures()

	def _filterIgnoredFeatures(self):
		ignore_features="".join(self.ignoreFeatures).replace(" ","").split(",")
		
		if len(ignore_features) >=1:
			PA_FT_df=pd.read_csv(self.CPAS_BED_FILE,sep="\t",header=None,index_col=None)
			PA_FT_df.columns=['chr','start','end','GeneID','Feature','strand']
			for ig in ignore_features:
				PA_FT_df=PA_FT_df[~PA_FT_df['Feature'].str.contains(ig)]
			PA_FT_df.to_csv(self.CPAS_BED_FILE,sep="\t",header=False,index=False)

	def _generateSAF(self):
		df = pd.read_csv(self.CPAS_BED_FILE, sep='\t', index_col=None, names=['Chr', 'Start', 'End', 'Gene', 'APA', 'Strand'])
		df['GeneID'] = df['Gene'] + '@' + df['Chr'] + '_' + df['Start'].apply(str) + '_' + df['End'].apply(str) + '_' + df['Strand']+ '@' + df['APA']
		df = df[['GeneID', 'Chr', 'Start', 'End', 'Strand']]
		df.to_csv(self.CPAS_SAF_FILE, sep='\t', index=False, header=None)

	def extractPolyA(self):
		# Extract soft-clipped (or polyA-capped) reads from each BAM file and save them as separate BAM and BED files; concatenated capped BED files and save as Global Variable "self.CONCATENATED_SOFTCLIPPED_BED_FILE"
		self._extractSoftClippedReads()

		# Filter extracted soft-clipped reads for only reads that contain a C/PAS site
		CPASBERT1 = CPASBERT(outDir = self.outDir,
		outPrefix = self.outPrefix, 
		modelPath = self.modelPath,
		fasta = self.FASTA,
		gtf = self.GTF,
		predictions = "",
		denovoapasiteBED = self.CONCATENATED_SOFTCLIPPED_BED_FILE,
		geneBED = "",
		polyADB = self.POLYADB,
		polyASite = self.POLYASITE
		)
		CPASBERT1.filter_CPAS_Sites()
		BedTool(self.CONCATENATED_SOFTCLIPPED_BED_FILE).saveas(self.outDir + self.outPrefix + ".concatenated.softclipped_afterCPASBERTFiltering.bed")

		# If user requests, append filtered APriori Annotations to BED File containing filtered soft-clipped reads. 
		if (self.apriori_annotations):
			self._loadFilteredAPrioriAnnotations()
			BedTool(self.CPAS_BED_FILE).saveas(self.outDir + self.outPrefix + ".concatenated.softclipped_AfterAprioriAnnotations.bed")
		# Otherwise, continue using filtered soft-clipped BED File (with modified file name)
		else:
			cmd = "cat " + self.CONCATENATED_SOFTCLIPPED_BED_FILE + " > " + self.CPAS_BED_FILE
			# print(cmd)
			os.system(cmd)
			BedTool(self.CPAS_BED_FILE).saveas(self.outDir + self.outPrefix + ".concatenated.softclipped_AfterNoAprioriAnnotations.bed")

		# If user requests, cluster on a user-determined distance and perform soft-clipped assisted filtering, where we only keep clusters that are supported by a soft-clipped read.
		self._performSoftClippedAssistedFiltering()
		# # Otherwise, cluster on a user-determined distance.
		# else:
		# 	BedTool(self.CPAS_BED_FILE).sort().cluster(s = True, d = self.clusterParameter).saveas(self.CPAS_BED_FILE)
		# 	BedTool(self.CPAS_BED_FILE).saveas(self.outDir + self.outPrefix + ".concatenated.softclipped_Clustered_BeforeTransformations.bed")

		# Given a BED file, keep only the most proximal read per cluster. Then, extend (or slop) the read by a user-determined distance on both ends. 
		self._transformCPASBED()

		# exit()
		# Given a BED file of the desired C/PAS coordinates, annotate C/PAS coordinate BED File with gene name and features. 
		self._annotateCPASBED()
		BedTool(self.CPAS_BED_FILE).saveas(self.outDir + self.outPrefix + ".concatenated.softclipped_AfterAnnotations.bed")

		# Remove (or ignore) any user-specified annotations from the C/PAS coordinate BED File
		self._filterIgnoredFeatures()
		BedTool(self.CPAS_BED_FILE).saveas(self.outDir + self.outPrefix + ".concatenated.softclipped_AfterFilteringIgnoredFeatures.bed")

		# Generate final SAF file for downstream APA calculations.
		self._generateSAF()

		return True

	def _mergePolyADBandPolyASite(self):
		#Make PolyADB and PolyASite 2.0 into similarly structured dataframe with a PAS_ID column
		df1 = pd.read_csv(self.POLYADB, sep = "\t", header = None, names = ["Chromosome", "Start", "End", "Feature1", "Feature2", "Strand"])
		df1['Start'] = df1['Start'].astype(str)
		df1['End'] = df1['End'].astype(str)
		df1["PAS_ID"] = df1["Chromosome"]+":"+df1["Start"]+"-"+df1["End"]+":"+df1["Strand"]
		cols = ['Chromosome', 'Start', 'End', 'Feature1', 'PAS_ID', 'Strand']
		df1 = df1[cols]

		df2 = pd.read_csv(self.POLYASITE, sep = "\t", header = None)
		df2 = df2.iloc[: , :6]
		df2.columns = ['Chromosome', 'Start', 'End', 'Feature1', 'Feature2', 'Strand']
		df2['Chromosome'] = "chr"+df2['Chromosome'].astype(str)
		df2[['Start', 'End']] = df2[['Start', 'End']].astype("int")
		df2["Midpoint"] = df2[['Start', 'End']].mean(axis=1).astype(int)
		df2['Start'] = df2['Midpoint'] - 1
		df2['End'] = df2['Midpoint']
		df2['Start'] = df2['Start'].astype(str)
		df2['End'] = df2['End'].astype(str)
		df2["PAS_ID"] = df2["Chromosome"]+":"+df2["Start"]+"-"+df2["End"]+":"+df2["Strand"]
		cols = ['Chromosome', 'Start', 'End', 'Midpoint', 'PAS_ID', 'Strand']
		df2 = df2[cols]

		pieces = (df1, df2)
		df_final = pd.concat(pieces, ignore_index = True)
		df_final = df_final.drop(columns=['Midpoint'])
		APRIORI_BED_FILE = self.outDir + self.outPrefix + "APrioriAnnotations.capped.bed"
		df_final.to_csv(APRIORI_BED_FILE, sep='\t', index=False, header=None)

		# cmd = "samtools faidx " + self.FASTA
		# print(cmd)
		# os.system(cmd)

		# FASTA_INDEX = self.FASTA.replace(".fa", ".fa.fai")
		# CHROM_SIZES = self.outDir + self.outPrefix + "chrom.sizes"
		# cmd = "cut -f 1,2 " + FASTA_INDEX + " > " + CHROM_SIZES
		# print(cmd)
		# os.system(cmd)

		# df3 = pd.read_csv(CHROM_SIZES, sep = "\t", header = None, names = ["Chromosome", "Size"])
		# ChromosomeNames = df3['Chromosome'].tolist()
		# df_final = df_final[df_final.Chromosome.isin(ChromosomeNames)]
		# df_final.to_csv(CAPPED_BED_FILE, sep='\t', index=False, header=None)

		# SLOPPED_CAPPED_BED_FILE = self.outDir + os.path.basename(CAPPED_BED_FILE).replace(".capped.bed", ".slopped.capped.bed")
		# cmd = "bedtools slop -i " + CAPPED_BED_FILE + " -b " + self.slopDistanceParameter + " -g " + CHROM_SIZES + " > " + SLOPPED_CAPPED_BED_FILE
		# print(cmd)
		# os.system(cmd)

		# SORTED_SLOPPED_CAPPED_BED_FILE = self.outDir + os.path.basename(CAPPED_BED_FILE).replace(".capped.bed", ".sorted.slopped.capped.bed")
		# cmd = "bedtools sort -i " + SLOPPED_CAPPED_BED_FILE + " > " + SORTED_SLOPPED_CAPPED_BED_FILE
		# print(cmd)
		# os.system(cmd)

	def _loadFilteredPolyADBandPolyASite(self):
		#Generate final Filtered BED File
		CONCATENATED_APRIORI_SOFTCLIPPED_BED_FILE = self.outDir + self.outPrefix + ".concatenated.apriori.softclipped.bed"
		APRIORI_DF = pd.read_csv(self.APRIORI, sep = "\t", header = None, names = ["Chr", "Start", "End", "Feature1", "Feature2", "Strand"])
		APRIORI_DF[['Start', 'End']] = APRIORI_DF[['Start', 'End']].astype("int")
		APRIORI_DF["Midpoint"] = APRIORI_DF[['Start', 'End']].mean(axis=1)
		APRIORI_DF['Midpoint'] = APRIORI_DF['Midpoint'].apply(lambda x: round(x, 0))
		APRIORI_DF = APRIORI_DF.astype({"Midpoint": int})
		APRIORI_DF["Start"] = APRIORI_DF["Midpoint"] - 1
		APRIORI_DF["End"] = APRIORI_DF["Midpoint"]
		APRIORI_DF = APRIORI_DF.drop('Midpoint', axis=1)
		modifiedAPrioriFileLoc = self.outDir + self.outPrefix + ".modified.aprioriAnnotations.bed"
		APRIORI_DF.to_csv(modifiedAPrioriFileLoc, sep='\t', index=False, header=None)

		cmd = "cat " + modifiedAPrioriFileLoc + " " + self.CONCATENATED_CAPPED_BED_FILE + " > " + self.FINAL_SOFTCLIPPED_BED_FILE
		# print(cmd)
		os.system(cmd)

		# self._clusterAndSortBED(inputBED = CONCATENATED_APRIORI_SOFTCLIPPED_BED_FILE)

	def filterPolyADBandPolyASite(self):
		self._mergePolyADBandPolyASite()
		APRIORI_BED_FILE = self.outDir + self.outPrefix + "APrioriAnnotations.capped.bed"

		CPASBERT1 = CPASBERT(outDir = self.outDir,
		outPrefix = self.outPrefix, 
		modelPath = self.modelPath,
		fasta = self.FASTA,
		gtf = self.GTF,
		predictions = "",
		denovoapasiteBED = APRIORI_BED_FILE,
		geneBED = "",
		polyADB = self.POLYADB,
		polyASite = self.POLYASITE
		)

		CPASBERT1.filter_CPAS_Sites()

		# SORTED_SLOPPED_CAPPED_APRIORIANNOTATIONS_BED_FILE = self.outDir + self.outPrefix + "APrioriAnnotations.sorted.slopped.capped.bed"
		self._clusterAndSortBED(inputBED = APRIORI_BED_FILE)

	def _clusterAndSortBED(self, inputBED):
		#Look at _generateFinalSoftClippedBEDFile and filterPolyADBandPolyASite functions; take the common code and refactor into this function
		
		CLUSTERED_BED_FILE = self.outDir + self.outPrefix + ".clustered.bed"
		cmd = "bedtools cluster -s -d " + self.clusterParameter + " -i " + inputBED + " > " + CLUSTERED_BED_FILE
		# print(cmd)
		os.system(cmd)

		df = pd.read_csv(CLUSTERED_BED_FILE, sep = "\t", header = None, names = ["Chr", "Start", "End", "Feature", "Feature2", "Strand", "Cluster"])
		Filtered_DF_posStrand = df.query('Strand == "+"').reset_index(drop=True)
		Filtered_DF_negStrand = df.query('Strand == "-"').reset_index(drop=True)
		Filtered_DF_posStrand = Filtered_DF_posStrand.loc[Filtered_DF_posStrand.groupby('Cluster').Start.idxmin()]
		Filtered_DF_negStrand = Filtered_DF_negStrand.loc[Filtered_DF_negStrand.groupby('Cluster').End.idxmax()]

		pieces = (Filtered_DF_posStrand, Filtered_DF_negStrand)
		Filtered_DF = pd.concat(pieces, ignore_index = True)

		CLEANED_BED_FILE = self.outDir + self.outPrefix + ".cleaned.bed"
		Filtered_DF.to_csv(CLEANED_BED_FILE, sep='\t', index=False, header=None)
		Filtered_DF = Filtered_DF.drop(["Cluster"], axis=1)
		Filtered_DF.to_csv(self.FINAL_SOFTCLIPPED_BED_FILE, sep='\t', index=False, header=None)

		# os.system("rm " + self.FINAL_SOFTCLIPPED_BED_FILE)
		# cmd = "bedtools sort -i " + CLEANED_BED_FILE + " > " + self.FINAL_SOFTCLIPPED_BED_FILE
		# print(cmd)
		# os.system(cmd)

		self._addGeneNames()
		self._markFeatures()

	def testExtractPolyA(self):
		self.CPAS_BED_FILE = "/mnt/belinda_local/venkata/data/Project_Human_RBM17_HEK/PolyAMiner_Results/Run5/00TESTING/TestingPlayground/3UTROnly.concatenated.softclipped_afterCPASBERTFiltering.bed"
		# If user requests, cluster on a user-determined distance and perform soft-clipped assisted filtering, where we only keep clusters that are supported by a soft-clipped read.
		# if (self.softclip_numReads > 0 and self.softclip_numSamples > 0):
		self._performSoftClippedAssistedFiltering()
		# else:
		# 	BedTool(self.CPAS_BED_FILE).sort().cluster(s = True, d = self.clusterParameter).saveas(self.CPAS_BED_FILE)

		# # Given a BED file, keep only the most proximal read per cluster. Then, extend (or slop) the read by a user-determined distance on both ends. 
		# self._transformCPASBED()

		# # Given a BED file of the desired C/PAS coordinates, annotate C/PAS coordinate BED File with gene name and features. 
		# self._annotateCPASBED()

def main():
	print("Testing ExtractPolyAsites4Bulk.py module (Filtering new BAM files portion)... ")
	ExtractPolyAsites4Bulk1 = ExtractPolyAsites4Bulk(
		outDir = "/mnt/belinda_local/venkata/data/Project_Human_RBM17_HEK/PolyAMiner_Results/Run5/00TESTING/TestingPlayground", 
		outPrefix = "3UTROnly", 
		fasta = "/mnt/belinda_local/venkata/data/Index_Files/Human/GenomeFasta_GTF/GRCh38.primary_assembly.genome.fa", 
		gtf = "/mnt/belinda_local/venkata/data/Index_Files/Human/GenomeFasta_GTF/gencode.v33.primary_assembly.annotation.gtf", 
		con1BAMFiles = "/mnt/belinda_local/venkata/data/Project_Human_RBM17_HEK/BAM/HZ8169/HZ8169_.sorted.bam,/mnt/belinda_local/venkata/data/Project_Human_RBM17_HEK/BAM/HZ8170/HZ8170_.sorted.bam,/mnt/belinda_local/venkata/data/Project_Human_RBM17_HEK/BAM/HZ8171/HZ8171_.sorted.bam", 
		con2BAMFiles = "/mnt/belinda_local/venkata/data/Project_Human_RBM17_HEK/BAM/HZ8162/HZ8162_.sorted.bam,/mnt/belinda_local/venkata/data/Project_Human_RBM17_HEK/BAM/HZ8163/HZ8163_.sorted.bam,/mnt/belinda_local/venkata/data/Project_Human_RBM17_HEK/BAM/HZ8164/HZ8164_.sorted.bam", 
		proportionA = "0.90,0.85,0.80,0.75", 
		modelOrganism = "human", 
		apriori_annotations = False, 
		ignoreFeatures = "",
		softclip_numSamples = 1,
		softclip_numReads = 1,
		clusterParameter = 30,
		slopDistanceParameter = 25
		)

	ExtractPolyAsites4Bulk1.testExtractPolyA()

if __name__ == "__main__":
	main()

