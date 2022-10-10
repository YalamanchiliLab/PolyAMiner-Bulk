import pandas as pd
import os, sys


class SoftclippedAssistedFiltering:
	def __init__(self, outDir, outPrefix, con1BAMFiles, con2BAMFiles, softclip_numSamples, softclip_numReads, clusterParameter, fasta):
		self.outDir = outDir.rstrip("/")+"/"
		self.outPrefix = outPrefix

		self.FINAL_SAF_FILE = self.outDir + self.outPrefix + "_denovoAPAsites.saf"
		self.SAF_BED_FILE = self.outDir + self.outPrefix + "_denovoAPAsites.SAF.bed"
		self.SUMMARIZED_SOFTCLIPPED_COUNTS = self.outDir + self.outPrefix + "summarized.softclipped.counts"

		self.con1BAMFiles = "".join(con1BAMFiles).replace(" ","").split(",")
		self.con2BAMFiles = "".join(con2BAMFiles).replace(" ","").split(",")
		self.conditionBAMList = self.con1BAMFiles + self.con2BAMFiles

		self.softclip_numReads = softclip_numReads
		self.softclip_numSamples = softclip_numSamples
		self.clusterParameter = clusterParameter

		self.FASTA = fasta

	def _convertSAFtoBED(self, inputSAF, outputBED):
		FullSAFDF = pd.read_csv(inputSAF, sep = "\t", header = None)
		FullSAFDF = FullSAFDF[[0]]
		FullSAFDF[['Gene', "BEDEntry",'Feature']] = FullSAFDF[0].astype("string").str.split('@',expand=True)
		FullSAFDF[['Chr', "Start",'End', "Strand"]] = FullSAFDF["BEDEntry"].astype("string").str.split('_',expand=True)
		FullSAFDF["PAS_ID"] = FullSAFDF['BEDEntry'] + "@" + FullSAFDF['Feature']
		#take midpoint of featureID
		FullSAFDF[['Start', 'End']] = FullSAFDF[['Start', 'End']].astype("int")
		FullSAFDF["Midpoint"] = FullSAFDF[['Start', 'End']].mean(axis=1)
		FullSAFDF['Midpoint'] = FullSAFDF['Midpoint'].apply(lambda x: round(x, 0))
		FullSAFDF = FullSAFDF.astype({"Midpoint": int})
		FullSAFDF["Start"] = FullSAFDF["Midpoint"] - 1
		FullSAFDF["End"] = FullSAFDF["Midpoint"]
		FullSAFDF = FullSAFDF.drop('Midpoint', axis=1)
		# FullSAFDF = FullSAFDF.drop('index', axis=1)
		FullSAFDF.reset_index(inplace=True)
		FullSAFDF = FullSAFDF[["Chr", "Start", "End", "Gene", "PAS_ID", "Strand"]]

		intermediateBED =  self.outDir + self.outPrefix + "intermediate.pre_slop.saf.bed"
		FullSAFDF.to_csv(intermediateBED, sep = "\t", index=False, header=None)		

		#slop by cluster parameter
		cmd = "samtools faidx " + self.FASTA
		print(cmd)
		os.system(cmd)

		FASTA_INDEX = self.FASTA.replace(".fa", ".fa.fai")
		CHROM_SIZES = self.outDir + self.outPrefix + "chrom.sizes"
		cmd = "cut -f 1,2 " + FASTA_INDEX + " > " + CHROM_SIZES
		print(cmd)
		os.system(cmd)

		cmd = "bedtools slop -i " + intermediateBED + " -b " + str(self.clusterParameter) + " -g " + CHROM_SIZES + " > " + outputBED
		print(cmd)
		os.system(cmd)

	def _createFilteredSoftclippedCountsBED_Deprecated(self):
		conditionBEDList = []
		cappedCountsFileList = []

		INTERMEDIATE_STEP3_SOFTCLIPPED_BED_FILE = self.outDir + self.outPrefix + "Intermediate.Step3.SoftClipped.bed"

		for BAM_FILE in self.conditionBAMList:
			cmd = "rm " + INTERMEDIATE_STEP3_SOFTCLIPPED_BED_FILE
			print(cmd)
			os.system(cmd)

			CAPPED_BED_FILE = self.outDir + self.outPrefix + os.path.basename(BAM_FILE).replace(".bam", ".capped.bed")

			cmd = "bedtools intersect -wao -s -a " + self.SAF_BED_FILE + " -b " + CAPPED_BED_FILE + " > " + INTERMEDIATE_STEP3_SOFTCLIPPED_BED_FILE
			print(cmd)
			os.system(cmd)

			DF = pd.read_csv(INTERMEDIATE_STEP3_SOFTCLIPPED_BED_FILE, sep = "\t", header = None)
			# cols = [0, 1, 2, 5]
			# DF['inter_1'] = DF[cols].apply(lambda row: '_'.join(row.values.astype(str)), axis=1)
			# cols = ["inter_1", 4]
			# DF['PAS_ID'] = DF[cols].apply(lambda row: '@'.join(row.values.astype(str)), axis=1)
			DF['PAS_ID'] = DF[4]
			DF = DF.groupby("PAS_ID").count()[[0]] - 1
			DF.columns = [os.path.basename(BAM_FILE).replace(".bam", ".softclipped_ct")]

			CAPPED_COUNTS_FILE = self.outDir + self.outPrefix + os.path.basename(BAM_FILE).replace(".bam", ".softclipped.counts")
			cappedCountsFileList.append(CAPPED_COUNTS_FILE)
			DF.to_csv(CAPPED_COUNTS_FILE, sep = "\t")

		main_dataframe = pd.DataFrame(pd.read_csv(cappedCountsFileList[0], sep = "\t"))

		for i in range(1,len(cappedCountsFileList)):
			data = pd.read_csv(cappedCountsFileList[i], sep = "\t")
			df_inter = pd.DataFrame(data)
			main_dataframe = pd.concat([main_dataframe, df_inter],axis=1)

		final_dataframe = main_dataframe.iloc[:, 1::2]
		final_dataframe["PAS_ID"] = main_dataframe.iloc[:, 0]

		SUMMARIZED_INTER_SOFTCLIPPED_COUNTS = self.outDir + self.outPrefix + "intermediate.summarized.softclipped.counts"
		final_dataframe.to_csv(SUMMARIZED_INTER_SOFTCLIPPED_COUNTS, sep = "\t", index = False)
		final_dataframe = final_dataframe.set_index('PAS_ID')

		KEEP_PASID_LIST = []

		for index, row in final_dataframe.iterrows():
			samples = [num for num in row.values if num >= self.softclip_numReads]
			if (len(samples) >= self.softclip_numSamples):
				KEEP_PASID_LIST.append(index)

		final_dataframe.reset_index(inplace=True)
		final_dataframe = final_dataframe[final_dataframe.PAS_ID.isin(KEEP_PASID_LIST)]
		final_dataframe.to_csv(self.SUMMARIZED_SOFTCLIPPED_COUNTS, sep = "\t", index = False)

	def _createFilteredSoftclippedCountsBED(self):
		conditionBEDList = []
		cappedCountsFileList = []

		INTERMEDIATE_STEP3_SOFTCLIPPED_BED_FILE = self.outDir + self.outPrefix + "Intermediate.Step3.SoftClipped.bed"
		CLUSTERED_BED_FILE = self.outDir + self.outPrefix + ".clustered.bed"

		for BAM_FILE in self.conditionBAMList:
			cmd = "rm " + INTERMEDIATE_STEP3_SOFTCLIPPED_BED_FILE
			print(cmd)
			os.system(cmd)

			CAPPED_BED_FILE = self.outDir + self.outPrefix + os.path.basename(BAM_FILE).replace(".bam", ".capped.bed")

			cmd = "bedtools intersect -wao -s -a " + CLUSTERED_BED_FILE + " -b " + CAPPED_BED_FILE + " > " + INTERMEDIATE_STEP3_SOFTCLIPPED_BED_FILE
			print(cmd)
			os.system(cmd)
			print("Terminated the program at line 132 in SoftclippedAssistedFilter.py")
			quit()

			DF = pd.read_csv(INTERMEDIATE_STEP3_SOFTCLIPPED_BED_FILE, sep = "\t", header = None)
			cols = [0, 1, 2, 5]
			DF['inter_1'] = DF[cols].apply(lambda row: '_'.join(row.values.astype(str)), axis=1)
			cols = ["inter_1", 4]
			DF['PAS_ID'] = DF[cols].apply(lambda row: '@'.join(row.values.astype(str)), axis=1)
			# DF['PAS_ID'] = DF[4]
			DF = DF.groupby("PAS_ID").count()[[0]] - 1
			DF.columns = [os.path.basename(BAM_FILE).replace(".bam", ".softclipped_ct")]

			CAPPED_COUNTS_FILE = self.outDir + self.outPrefix + os.path.basename(BAM_FILE).replace(".bam", ".softclipped.counts")
			cappedCountsFileList.append(CAPPED_COUNTS_FILE)
			DF.to_csv(CAPPED_COUNTS_FILE, sep = "\t")

		main_dataframe = pd.DataFrame(pd.read_csv(cappedCountsFileList[0], sep = "\t"))

		for i in range(1,len(cappedCountsFileList)):
			data = pd.read_csv(cappedCountsFileList[i], sep = "\t")
			df_inter = pd.DataFrame(data)
			main_dataframe = pd.concat([main_dataframe, df_inter],axis=1)

		final_dataframe = main_dataframe.iloc[:, 1::2]
		final_dataframe["PAS_ID"] = main_dataframe.iloc[:, 0]

		SUMMARIZED_INTER_SOFTCLIPPED_COUNTS = self.outDir + self.outPrefix + "intermediate.summarized.softclipped.counts"
		final_dataframe.to_csv(SUMMARIZED_INTER_SOFTCLIPPED_COUNTS, sep = "\t", index = False)
		final_dataframe = final_dataframe.set_index('PAS_ID')

		KEEP_PASID_LIST = []

		for index, row in final_dataframe.iterrows():
			samples = [num for num in row.values if num >= self.softclip_numReads]
			if (len(samples) >= self.softclip_numSamples):
				KEEP_PASID_LIST.append(index)

		final_dataframe.reset_index(inplace=True)
		final_dataframe = final_dataframe[final_dataframe.PAS_ID.isin(KEEP_PASID_LIST)]
		final_dataframe.to_csv(self.SUMMARIZED_SOFTCLIPPED_COUNTS, sep = "\t", index = False)

	def _updateMasterSAF(self):
		filteredSoftclippedDF = pd.read_csv(self.SUMMARIZED_SOFTCLIPPED_COUNTS, sep = "\t")
		KEEP_PASID_LIST = filteredSoftclippedDF["PAS_ID"]

		FullSAFDF = pd.read_csv(self.FINAL_SAF_FILE, sep = "\t", header = None, names = ["Full_SAF_Feature", "Chr", "Start", "End", "Strand"])
		BEFORE_FILTERING_SAF = self.outDir + self.outPrefix + "_denovoAPAsites_beforeSoftClippedAssitedFiltering.saf"
		FullSAFDF.to_csv(BEFORE_FILTERING_SAF, sep = "\t", index = False)
		FullSAFDF[['Gene', "BEDEntry",'Feature']] = FullSAFDF["Full_SAF_Feature"].astype("string").str.split('@',expand=True)
		FullSAFDF["PAS_ID"] = FullSAFDF['BEDEntry'] + "@" + FullSAFDF['Feature']
		FullSAFDF = FullSAFDF[FullSAFDF.PAS_ID.isin(KEEP_PASID_LIST)]
		FullSAFDF = FullSAFDF.iloc[: , :-4]

		#Comment out later
		# self.FINAL_SAF_FILE = self.outDir + self.outPrefix + "_denovoAPAsites_testingRun1.saf"

		FullSAFDF.to_csv(self.FINAL_SAF_FILE, sep = "\t", index = False, header = False)

	def performSoftclippedAssistedFiltering(self):
		self._convertSAFtoBED(inputSAF = self.FINAL_SAF_FILE, outputBED = self.SAF_BED_FILE)
		self._createFilteredSoftclippedCountsBED()
		self._updateMasterSAF()

	def test(self):
		self.SAF_BED_FILE = "/mnt/belinda_local/venkata/data/Project_Human_RBM17_HEK/PolyAMiner_Results/Run5/CtrlvsRBM17siRNA_3UTROnly_SoftClipped+Annotations_Run6.5_TESTRUN_fixedSoftClippedFiltering/intermediate_test1.post_slop.bed"
		
		self._createFilteredSoftclippedCountsBED()
		self._updateMasterSAF()

def main():
	SoftclippedAssistedFiltering1 = SoftclippedAssistedFiltering(outDir = "/mnt/belinda_local/venkata/data/Project_Human_RBM17_HEK/PolyAMiner_Results/Run5/CtrlvsRBM17siRNA_3UTROnly_SoftClipped+Annotations_Run6.5_TESTRUN_fixedSoftClippedFiltering",
		outPrefix = "3UTROnly",
		con1BAMFiles = "/mnt/belinda_local/venkata/data/Project_Human_RBM17_HEK/BAM/HZ8169/HZ8169_.sorted.bam,/mnt/belinda_local/venkata/data/Project_Human_RBM17_HEK/BAM/HZ8170/HZ8170_.sorted.bam,/mnt/belinda_local/venkata/data/Project_Human_RBM17_HEK/BAM/HZ8171/HZ8171_.sorted.bam",
		con2BAMFiles = "/mnt/belinda_local/venkata/data/Project_Human_RBM17_HEK/BAM/HZ8162/HZ8162_.sorted.bam,/mnt/belinda_local/venkata/data/Project_Human_RBM17_HEK/BAM/HZ8163/HZ8163_.sorted.bam,/mnt/belinda_local/venkata/data/Project_Human_RBM17_HEK/BAM/HZ8164/HZ8164_.sorted.bam",
		softclip_numReads = 1,
		softclip_numSamples = 1,
		clusterParameter = 30,
		fasta = ""
		)

	SoftclippedAssistedFiltering1.test()

if __name__ == "__main__":
    main()
