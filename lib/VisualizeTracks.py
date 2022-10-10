import os, sys, glob
import argparse
import pandas as pd

class VisualizeTracks:
	def __init__ (self, outDir, outPrefix, gtf, polyAResults, condition1SamplesBAM, condition2SamplesBAM, condition1Name, condition2Name, numTop):
		self.outDir = outDir.rstrip("/")+"/"
		self.outPrefix = outPrefix
		self.GTF = gtf
		self.polyAResults = polyAResults
		self.condition1SamplesBAM = "".join(condition1SamplesBAM).replace(" ","").split(",")
		self.condition2SamplesBAM = "".join(condition2SamplesBAM).replace(" ","").split(",")

		self.condition1Name = condition1Name
		self.condition2Name = condition2Name

		self.numTop = numTop

		self.condition1SamplesBW_FORWARD = []
		self.condition1SamplesBW_REVERSE = []

		self.condition2SamplesBW_FORWARD = []
		self.condition2SamplesBW_REVERSE = []

		self.CONFIG_FILEPATH_FORWARD = self.outDir + self.outPrefix + "forward.config.ini"
		self.CONFIG_FILEPATH_REVERSE = self.outDir + self.outPrefix + "reverse.config.ini"

	def _checkDir (self, Dir):
		DirNoSlash = Dir.rstrip("/")
		if (os.path.isdir(DirNoSlash)):
			os.system("rm -R " + DirNoSlash)
		os.system("mkdir " + DirNoSlash)

	def _convertBam2BW(self):
		counter = 1
		total = (len(self.condition1SamplesBAM) * 2) + (len(self.condition2SamplesBAM)*2)

		for file in self.condition1SamplesBAM:
			BASENAME_FORWARD = os.path.basename(file).replace(".bam","_forward_.bw")
			BASENAME_REVERSE = os.path.basename(file).replace(".bam","_reverse_.bw")
			OUTPUT_DIR = self.outDir + "Stranded_BW"
			if (os.path.isdir(OUTPUT_DIR) == False):
				os.system("mkdir " + OUTPUT_DIR)
			OUTPUT_FORWARD = OUTPUT_DIR + "/" + BASENAME_FORWARD
			OUTPUT_REVERSE = OUTPUT_DIR + "/" + BASENAME_REVERSE

			if os.path.exists(OUTPUT_FORWARD) == False:
				cmd = ("bamCoverage -b "+file+" -bs 5 -p 20 --normalizeUsing CPM --skipNonCoveredRegions --smoothLength 15 --centerReads --filterRNAstrand forward -o "+OUTPUT_FORWARD)
				print("Converting BAMs to BWs: " +  str(counter) + " of " + str(total))
				counter += 1
				print(cmd)
				os.system(cmd)
			if os.path.exists(OUTPUT_REVERSE) == False: 
				cmd = ("bamCoverage -b "+file+" -bs 5 -p 20 --normalizeUsing CPM --skipNonCoveredRegions --smoothLength 15 --centerReads --filterRNAstrand reverse -o "+OUTPUT_REVERSE)
				print("Converting BAMs to BWs: " +  str(counter) + " of " + str(total))
				counter += 1
				print(cmd)
				os.system(cmd)

			self.condition1SamplesBW_FORWARD.append(OUTPUT_FORWARD)
			self.condition1SamplesBW_REVERSE.append(OUTPUT_REVERSE)

		for file in self.condition2SamplesBAM:
			BASENAME_FORWARD = os.path.basename(file).replace(".bam","_forward_.bw")
			BASENAME_REVERSE = os.path.basename(file).replace(".bam","_reverse_.bw")
			OUTPUT_FORWARD = OUTPUT_DIR + "/" + BASENAME_FORWARD
			OUTPUT_REVERSE = OUTPUT_DIR + "/" + BASENAME_REVERSE

			if os.path.exists(OUTPUT_FORWARD) == False:
				cmd = ("bamCoverage -b "+file+" -bs 5 -p 20 --normalizeUsing CPM --skipNonCoveredRegions --smoothLength 15 --centerReads --filterRNAstrand forward -o "+OUTPUT_FORWARD)
				print("Converting BAMs to BWs: " +  str(counter) + " of " + str(total))
				counter += 1
				print(cmd)
				os.system(cmd)
			if os.path.exists(OUTPUT_REVERSE) == False: 
				cmd = ("bamCoverage -b "+file+" -bs 5 -p 20 --normalizeUsing CPM --skipNonCoveredRegions --smoothLength 15 --centerReads --filterRNAstrand reverse -o "+OUTPUT_REVERSE)
				print("Converting BAMs to BWs: " +  str(counter) + " of " + str(total))
				counter += 1
				print(cmd)
				os.system(cmd)

			self.condition2SamplesBW_FORWARD.append(OUTPUT_FORWARD)
			self.condition2SamplesBW_REVERSE.append(OUTPUT_REVERSE)

	def _makeConfigFile(self, strand):		
		if (strand == "forward"):
			CONFIG_FILEPATH = self.CONFIG_FILEPATH_FORWARD
			condition1SamplesBW = self.condition1SamplesBW_FORWARD
			condition2SamplesBW = self.condition2SamplesBW_FORWARD
		
		elif (strand == "reverse"):
			CONFIG_FILEPATH = self.CONFIG_FILEPATH_REVERSE
			condition1SamplesBW = self.condition1SamplesBW_REVERSE
			condition2SamplesBW = self.condition2SamplesBW_REVERSE

		fw = open(CONFIG_FILEPATH, 'w')

		if (len(condition1SamplesBW) > 4) or (len(condition1SamplesBW) > 4) :
			for i in range (0, len(condition1SamplesBW)):
				fw.write("[bigwig control file]\n")
				fw.write("file = " + condition1SamplesBW[i] + "\n")
				fw.write("height = 4\n")
				fw.write("color = green\n")
				fw.write("nans_to_zeros = true\n")
				fw.write("summary_method = mean\n")
				fw.write("show_data_range = true\n")
				fw.write("alpha = 0.5\n")
				fw.write("title = Control BigWigs (x"+ str(len(condition1SamplesBW)) +")\n")
				fw.write("min_value = 0\n")
				fw.write("overlay_previous = share-y\n")
				# fw.write("max_value = 0.5\n")
				
			fw.write("[spacer]\n")
			fw.write("height = 4\n")

			for i in range (0, len(condition2SamplesBW)):
				fw.write("[bigwig treatment file]\n")
				fw.write("file = " + condition2SamplesBW[i] + "\n")
				fw.write("height = 4\n")
				fw.write("color = red\n")
				fw.write("nans_to_zeros = true\n")
				fw.write("summary_method = mean\n")
				fw.write("show_data_range = true\n")
				fw.write("alpha = 0.5\n")
				fw.write("title = Treatment BigWigs (x"+ str(len(condition2SamplesBW)) +")\n")
				fw.write("min_value = 0\n")
				fw.write("overlay_previous = share-y\n")
				# fw.write("max_value = 0.5\n")
		else:
			for i in range (0, len(condition1SamplesBW)):
				fw.write("[bigwig control file]\n")
				fw.write("file = " + condition1SamplesBW[i] + "\n")
				fw.write("height = 4\n")
				fw.write("color = green\n")
				fw.write("nans_to_zeros = true\n")
				fw.write("summary_method = mean\n")
				fw.write("show_data_range = true\n")
				fw.write("alpha = 0.5\n")
				fw.write("title = " + self.condition1Name + " BigWig\n")
				fw.write("min_value = 0\n")
				# fw.write("overlay_previous = share-y\n")
				# fw.write("max_value = 0.5\n")
				
			# fw.write("[spacer]\n")
			# fw.write("height = 4\n")

			for i in range (0, len(condition2SamplesBW)):
				fw.write("[bigwig treatment file]\n")
				fw.write("file = " + condition2SamplesBW[i] + "\n")
				fw.write("height = 4\n")
				fw.write("color = red\n")
				fw.write("nans_to_zeros = true\n")
				fw.write("summary_method = mean\n")
				fw.write("show_data_range = true\n")
				fw.write("alpha = 0.5\n")
				fw.write("title = " + self.condition2Name + " BigWig\n")
				fw.write("min_value = 0\n")
				# fw.write("overlay_previous = share-y\n")
				# fw.write("max_value = 0.5\n")


		fw.write("[spacer]\n")
		# fw.write("height = 2\n")

		fw.write("[test gtf collapsed]\n")
		fw.write("file = " + self.GTF + "\n")
		fw.write("height = 3\n")
		fw.write("merge_transcripts = true\n")
		# fw.write("merge_overlapping_exons = true\n")
		fw.write("color_utr = purple\n")
		fw.write("labels = true\n")
		fw.write("height_utr = 0.4\n")
		fw.write("style = UCSC\n")
		fw.write("arrow_interval = 10\n")
		fw.write("arrowhead_included = true\n")
		fw.write("gene_rows = 1\n")
		fw.write("prefered_name = gene_name\n")
		fw.write("fontsize = 12\n")
		fw.write("display = stacked\n")
		fw.write("labels_in_margin = true\n")
		fw.write("file_type = gtf\n")

		fw.write("[C/PAS BED]\n")
		fw.write("file = " + self.formattedCPAS_BED_FileLoc + "\n")
		fw.write("height = 0\n")
		fw.write("fontsize = 12\n")
		fw.write("file_type = bed\n")
		fw.write("merge_transcripts = true\n")
		fw.write("overlay_previous = share-y\n")
		fw.write("color = #e3dc62\n")
		fw.write("fontstyle = italic\n")
		fw.write("orientation = inverted\n")
		fw.write("labels = false\n")
		fw.write("file_type = bed\n")

		fw.write("[Label for C/PAS BED]\n")
		fw.write("file = " + self.formattedLabeledCPAS_BED_FileLoc + "\n")
		fw.write("height = 0\n")
		fw.write("fontsize = 12\n")
		fw.write("file_type = bed\n")
		fw.write("merge_transcripts = true\n")
		# fw.write("overlay_previous = share-y\n")
		fw.write("color = #FFFFFF\n")
		fw.write("border_color = #FFFFFF\n")
		fw.write("fontstyle = italic\n")
		fw.write("orientation = inverted\n")
		fw.write("file_type = bed\n")

		fw.write("[vlines]\n")
		fw.write("file = " + self.formattedCPAS_BED_FileLoc + "\n")
		fw.write("type = vlines")

		# fw.write("[x-axis]\n")
		# fw.write("fontsize = 10\n")

		fw.close()

	def _parseResults(self, numTop, NegOrPosPolyAIndex):
		resultsDF = pd.read_csv(self.polyAResults, sep = "\t")
		resultsDF = resultsDF[resultsDF['AdjG_Pval'] < 0.05]

		if NegOrPosPolyAIndex == "Positive": 
			resultsDF = resultsDF[resultsDF['PolyAIndex'] > 0]
			resultsDF = resultsDF.sort_values(by='PolyAIndex', ascending=False)
		else:
			resultsDF = resultsDF[resultsDF['PolyAIndex'] < 0]
			resultsDF = resultsDF.sort_values(by='PolyAIndex', ascending=True)
		
		resultsDF = resultsDF.reset_index()
		resultsDF = resultsDF.iloc[:numTop]

		return resultsDF

	def _generatePlots(self, resultsDF):
		for index, row in resultsDF.iterrows():
			Gene = row["Symbol"]
			self.formattedCPAS_BED_FileLoc = self.outDir + self.outPrefix + str(Gene) + "_CPASdb.bed"
			self.formattedLabeledCPAS_BED_FileLoc = self.outDir + self.outPrefix + str(Gene) + "_LabeledCPASdb.bed"

			PolyASiteList = row["PolyASites"].split(",c")
			regionList = [i.split('@')[0].split('_') for i in PolyASiteList] 
			featureList = [i.split('@')[1] for i in PolyASiteList] 
			chrList = list(list(zip(*regionList))[0])
			chrList = ['c' + x if x[0] != 'c' else x for x in chrList]
			startList = list(list(zip(*regionList))[1])
			endList = list(list(zip(*regionList))[2])
			dummyList = [1000 for i in range(len(endList))]
			strandList = list(list(zip(*regionList))[3])
			CPAS_BED_DF = pd.DataFrame(list(zip(chrList, startList, endList, featureList, dummyList, strandList)))
			CPAS_BED_DF.to_csv(self.formattedCPAS_BED_FileLoc, sep = "\t", header = None, index = False)

			CPAS_BED_DF[1] = pd.to_numeric(CPAS_BED_DF[1]) - 160
			CPAS_BED_DF[2] = pd.to_numeric(CPAS_BED_DF[2]) - 160
			CPAS_BED_DF.to_csv(self.formattedLabeledCPAS_BED_FileLoc, sep = "\t", header = None, index = False)
			#Maybe slop CPAS BED by X coordinates on both sides???

			self._makeConfigFile(strand = "forward")
			self._makeConfigFile(strand = "reverse")

			OUTPUT_FILEPATH = self.outDir + self.outPrefix + str(index+1) + "_" + str(Gene) +".DAG_Track_WholeGeneView.png"
			chromosome = CPAS_BED_DF[0][0]
			start = int(CPAS_BED_DF[1].min()) - 2000
			end = int(CPAS_BED_DF[2].max()) + 2000
			strand = CPAS_BED_DF[5][0]
		
			if (start > end):
				temp1 = start
				temp2 = end
				end = temp1
				start = temp2
			region = chromosome + ":" + str(start) + "-" + str(end)
			
			if (strand == "+"):
				cmd = "pyGenomeTracks --tracks " + self.CONFIG_FILEPATH_FORWARD + " --region " + region + " --dpi 150 --fontSize 14 --trackLabelFraction 0 --width 50 --outFileName " + OUTPUT_FILEPATH
			elif (strand == "-"):
				cmd = "pyGenomeTracks --tracks " + self.CONFIG_FILEPATH_REVERSE + " --region " + region + " --dpi 150 --fontSize 14 --trackLabelFraction 0 --width 50 --outFileName " + OUTPUT_FILEPATH

			os.system(cmd)
			
			try:
				OUTPUT_FILEPATH = self.outDir + self.outPrefix + str(index+1) + "_" + str(Gene) +".DAG_Track_3UTRView.png" 
				CPAS_BED_DF = CPAS_BED_DF[CPAS_BED_DF[3].str.contains("UTR3") | CPAS_BED_DF[3].str.contains("UN")]
				start = int(CPAS_BED_DF[1].min()) - 2000
				end = int(CPAS_BED_DF[2].max()) + 2000
				if (start > end):
					temp1 = start
					temp2 = end
					end = temp1
					start = temp2
				region = chromosome + ":" + str(start) + "-" + str(end)

				if (strand == "+"):
					cmd = "pyGenomeTracks --tracks " + self.CONFIG_FILEPATH_FORWARD + " --region " + region + " --dpi 150 --fontSize 14 --trackLabelFraction 0 --width 50 --outFileName " + OUTPUT_FILEPATH
				elif (strand == "-"):
					cmd = "pyGenomeTracks --tracks " + self.CONFIG_FILEPATH_REVERSE + " --region " + region + " --dpi 150 --fontSize 14 --trackLabelFraction 0 --width 50 --outFileName " + OUTPUT_FILEPATH
				os.system(cmd)
			except:
				print(str(Gene) + "has no 3'UTR or UN C/PASs.....")

	def visualizeTopDAGs(self):
		self._convertBam2BW()

		commonBase = self.outDir

		self.outDir = commonBase + "Graphics_NegPolyAIndex/"
		self._checkDir(self.outDir)
		resultsDF = self._parseResults(numTop = self.numTop, NegOrPosPolyAIndex = "Negative")
		self._generatePlots(resultsDF)

		self.outDir = commonBase + "Graphics_PosPolyAIndex/"
		self._checkDir(self.outDir)
		resultsDF = self._parseResults(numTop = self.numTop, NegOrPosPolyAIndex = "Positive")
		self._generatePlots(resultsDF)

def main ():
	parser = argparse.ArgumentParser(description='''PolyAMiner-Bulk Visualization Module: Visualize PolyAMiner-Bulk Results - Venkata Jonnakuti et al., \n''',formatter_class=argparse.ArgumentDefaultsHelpFormatter)

	required = parser.add_argument_group('Required arguments')

	required.add_argument("-o", help = 'Output directory', type = str, default = 'PolyAminer_OUT')
	required.add_argument("-outPrefix", help = 'Output file/s prefix', default = "PolyAminer_Out", type = str)
	required.add_argument('-gtf',help = 'Reference gtf file',required = 'True', type = str)
	required.add_argument("-polyAResults", help = "PolyAMiner-Bulk Results File", required = "True", type = str)
	required.add_argument('-c1',help='Comma-separated list of condition1 BAM files in full path format. Index files are also expected', nargs='+',required='True',type=str)
	required.add_argument('-c2',help='Comma-separated list of condition2 BAM files in full path format. Index files are also expected', nargs='+',required='True',type=str)
	required.add_argument("-c1Name", help = "Condition 1 Sample Name", type = str, default = "Control")
	required.add_argument("-c2Name", help = "Condition 2 Sample Name", type = str, default = "Treatment")
	required.add_argument("-numTop", help = "Number of significant DAGs to visualize", type = int, default = 100)

	args = parser.parse_args()
	args_dict = vars(args)
	args.o=args.o.rstrip("/")

	VisualizeTracks1 = VisualizeTracks(outDir = args.o,
		outPrefix = args.outPrefix,
		gtf = args.gtf,
		polyAResults = args.polyAResults,
		condition1SamplesBAM = args.c1,
		condition2SamplesBAM = args.c2,
		condition1Name = args.c1Name,
		condition2Name = args.c2Name,
		numTop = args.numTop
		)

	# VisualizeTracks1 = VisualizeTracks(outDir = "/mnt/belinda_local/venkata/data/Project_Human_RBM17_HEK/TEST_VISUALIZATIONS",
	# 	outPrefix = "TEST_",
	# 	gtf = "/mnt/belinda_local/venkata/data/Index_Files/Human/GenomeFasta_GTF/gencode.v33.primary_assembly.annotation.gtf",
	# 	polyAResults = "/mnt/belinda_local/venkata/data/Project_Human_RBM17_HEK/PolyAMiner_Results/CtrlvsRBM17siRNA_3UTROnly_SoftClipped+Annotations_Run4/3UTROnly_PolyA-miner.Results.txt",
	# 	# condition1SamplesBW = "/mnt/belinda_local/venkata/data/Project_Human_RBM17_HEK/BW/HZ8169_.sorted.bw,/mnt/belinda_local/venkata/data/Project_Human_RBM17_HEK/BW/HZ8170_.sorted.bw,/mnt/belinda_local/venkata/data/Project_Human_RBM17_HEK/BW/HZ8171_.sorted.bw",
	# 	# condition2SamplesBW = "/mnt/belinda_local/venkata/data/Project_Human_RBM17_HEK/BW/HZ8162_.sorted.bw,/mnt/belinda_local/venkata/data/Project_Human_RBM17_HEK/BW/HZ8163_.sorted.bw,/mnt/belinda_local/venkata/data/Project_Human_RBM17_HEK/BW/HZ8164_.sorted.bw",
	# 	condition1SamplesBAM = "/mnt/belinda_local/venkata/data/Project_Human_RBM17_HEK/BAM/HZ8169/HZ8169_.sorted.bam,/mnt/belinda_local/venkata/data/Project_Human_RBM17_HEK/BAM/HZ8170/HZ8170_.sorted.bam,/mnt/belinda_local/venkata/data/Project_Human_RBM17_HEK/BAM/HZ8171/HZ8171_.sorted.bam",
	# 	condition2SamplesBAM = "/mnt/belinda_local/venkata/data/Project_Human_RBM17_HEK/BAM/HZ8162/HZ8162_.sorted.bam,/mnt/belinda_local/venkata/data/Project_Human_RBM17_HEK/BAM/HZ8163/HZ8163_.sorted.bam,/mnt/belinda_local/venkata/data/Project_Human_RBM17_HEK/BAM/HZ8164/HZ8164_.sorted.bam",
	# 	numTop = 100
	# 	)

	# VisualizeTracks1 = VisualizeTracks(outDir = "/mnt/belinda_local/venkata/data/Project_Meningioma_AkashPatel_NSG/PolyAMiner_Results/SubtypeCvsA_3UTROnly_Run6_Visualizations",
	# 	outPrefix = "",
	# 	gtf = "/mnt/belinda_local/venkata/data/Index_Files/Human/GenomeFasta_GTF/gencode.v33.primary_assembly.annotation.gtf",
	# 	polyAResults = "/mnt/belinda_local/venkata/data/Project_Meningioma_AkashPatel_NSG/PolyAMiner_Results/SubtypeCvsA_3UTROnly_Run6/3UTROnly_PolyA-miner.Results.txt",
	# 	condition1SamplesBAM = "/mnt/belinda_local/venkata/data/Project_Meningioma_AkashPatel_NSG/hari_APA_Akash/02_BAM/TL-21-VZKP229D/TL-21-VZKP229D.sorted.bam,/mnt/belinda_local/venkata/data/Project_Meningioma_AkashPatel_NSG/hari_APA_Akash/02_BAM/TL-21-QGUU886F/TL-21-QGUU886F.sorted.bam,/mnt/belinda_local/venkata/data/Project_Meningioma_AkashPatel_NSG/hari_APA_Akash/02_BAM/TL-20-DF4101/TL-20-DF4101.sorted.bam,/mnt/belinda_local/venkata/data/Project_Meningioma_AkashPatel_NSG/hari_APA_Akash/02_BAM/TL-20-36A961/TL-20-36A961.sorted.bam",
	# 	condition2SamplesBAM = "/mnt/belinda_local/venkata/data/Project_Meningioma_AkashPatel_NSG/hari_APA_Akash/02_BAM/TL-20-56286D/TL-20-56286D.sorted.bam,/mnt/belinda_local/venkata/data/Project_Meningioma_AkashPatel_NSG/hari_APA_Akash/02_BAM/TL-20-2F7E80/TL-20-2F7E80.sorted.bam,/mnt/belinda_local/venkata/data/Project_Meningioma_AkashPatel_NSG/hari_APA_Akash/02_BAM/TL-19-EBC5FF/TL-19-EBC5FF.sorted.bam,/mnt/belinda_local/venkata/data/Project_Meningioma_AkashPatel_NSG/hari_APA_Akash/02_BAM/TL-19-C46B1C/TL-19-C46B1C.sorted.bam",
	# 	numTop = 1000
	# 	)

	VisualizeTracks1.visualizeTopDAGs()

if __name__ == "__main__":
	main()