import os, sys
import pandas as pd

class VisualizeTracks:
	def __init__ (self, outDir, outPrefix, gtf, polyAResults, condition1SamplesBW, condition2SamplesBW, numTop):
		self.outDir = outDir.rstrip("/")+"/"
		self.outPrefix = outPrefix
		self.GTF = gtf
		self.polyAResults = polyAResults
		self.condition1SamplesBW = condition1SamplesBW.split(",")
		self.condition2SamplesBW = condition2SamplesBW.split(",")
		self.numTop = numTop

		self.CONFIG_FILEPATH = self.outDir + self.outPrefix + ".config.ini"


	def _convertBam2BW(self):
		return None

	def _makeConfigFile(self):
		fw = open(self.CONFIG_FILEPATH, 'w')
		for i in range (0, len(self.condition1SamplesBW)):
			fw.write("[bigwig control file]\n")
			fw.write("file = " + self.condition1SamplesBW[i] + "\n")
			fw.write("height = 4\n")
			fw.write("color = green\n")
			fw.write("nans_to_zeros = true\n")
			fw.write("summary_method = mean\n")
			fw.write("show_data_range = true\n")
			fw.write("alpha = 0.5\n")
			fw.write("title = Control BigWigs (x"+ str(len(self.condition1SamplesBW)) +")\n")
			fw.write("min_value = 0\n")
			fw.write("overlay_previous = share-y\n")
			fw.write("max_value = 0.5\n")
			
		fw.write("[spacer]\n")
		fw.write("height = 4\n")

		for i in range (0, len(self.condition2SamplesBW)):
			fw.write("[bigwig treatment file]\n")
			fw.write("file = " + self.condition2SamplesBW[i] + "\n")
			fw.write("height = 4\n")
			fw.write("color = red\n")
			fw.write("nans_to_zeros = true\n")
			fw.write("summary_method = mean\n")
			fw.write("show_data_range = true\n")
			fw.write("alpha = 0.5\n")
			fw.write("title = Treatment BigWigs (x"+ str(len(self.condition2SamplesBW)) +")\n")
			fw.write("min_value = 0\n")
			fw.write("overlay_previous = share-y\n")
			fw.write("max_value = 0.5\n")

		fw.write("[spacer]\n")
		# fw.write("height = 2\n")

		fw.write("[test gtf collapsed]\n")
		fw.write("file = " + self.GTF + "\n")
		fw.write("height = 4\n")
		fw.write("merge_transcripts = true\n")
		fw.write("merge_overlapping_exons = true\n")
		fw.write("color_utr = purple\n")
		fw.write("labels = true\n")
		fw.write("height_utr = 0.4\n")
		fw.write("style = flybase\n")
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
		fw.write("height = 7\n")
		fw.write("fontsize = 10\n")
		fw.write("file_type = bed\n")
		fw.write("merge_transcripts = true\n")
		fw.write("overlay_previous = yes\n")
		fw.write("color = #e3dc62\n")
		fw.write("all_labels_inside = true\n")
		fw.write("file_type = bed\n")

		fw.write("[x-axis]\n")
		fw.write("fontsize = 10\n")

		fw.write("[vlines]\n")
		fw.write("file = " + self.formattedCPAS_BED_FileLoc + "\n")
		fw.write("type = vlines")

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
			self.formattedCPAS_BED_FileLoc = self.outDir + self.outPrefix + Gene + "_CPASdb.bed"

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
			#Maybe slop CPAS BED by X coordinates on both sides???

			self._makeConfigFile()

			OUTPUT_FILEPATH = self.outDir + self.outPrefix + str(index+1) + "_" + Gene +".DAG_Track_WholeGeneView.png"
			chromosome = CPAS_BED_DF[0][0]
			start = int(CPAS_BED_DF[1].min()) - 10000
			end = int(CPAS_BED_DF[2].max()) + 10000
			if (start > end):
				temp1 = start
				temp2 = end
				end = temp1
				start = temp2
			region = chromosome + ":" + str(start) + "-" + str(end)
			cmd = "pyGenomeTracks --tracks " + self.CONFIG_FILEPATH + " --region " + region + " --dpi 150 --fontSize 14 --trackLabelFraction 0 --width 50 --outFileName " + OUTPUT_FILEPATH
			os.system(cmd)
			
			try:
				OUTPUT_FILEPATH = self.outDir + self.outPrefix + str(index+1) + "_" + Gene +".DAG_Track_3UTRView.png" 
				CPAS_BED_DF = CPAS_BED_DF[CPAS_BED_DF[3].str.contains("UTR3") | CPAS_BED_DF[3].str.contains("UN")]
				start = int(CPAS_BED_DF[1].min()) - 10000
				end = int(CPAS_BED_DF[2].max()) + 10000
				if (start > end):
					temp1 = start
					temp2 = end
					end = temp1
					start = temp2
				region = chromosome + ":" + str(start) + "-" + str(end)
				cmd = "pyGenomeTracks --tracks " + self.CONFIG_FILEPATH + " --region " + region + " --dpi 150 --fontSize 14 --trackLabelFraction 0 --width 50 --outFileName " + OUTPUT_FILEPATH
				os.system(cmd)
			except:
				print(Gene + "has no 3'UTR or UN C/PASs.....")

	def visualizeTopDAGs(self):
		commonBase = self.outDir
		self.outDir = commonBase + "Graphics_PosPolyAIndex/"
		outDirNoSlash = self.outDir.rstrip("/")
		if (os.path.isdir(outDirNoSlash)):
			os.system("rm -R " + outDirNoSlash)
		os.system("mkdir " + outDirNoSlash)

		resultsDF = self._parseResults(numTop = self.numTop, NegOrPosPolyAIndex = "Positive")
		self._generatePlots(resultsDF)

		self.outDir = commonBase + "Graphics_NegPolyAIndex/"
		outDirNoSlash = self.outDir.rstrip("/")
		if (os.path.isdir(outDirNoSlash)):
			os.system("rm -R " + outDirNoSlash)
		os.system("mkdir " + outDirNoSlash)

		resultsDF = self._parseResults(numTop = self.numTop, NegOrPosPolyAIndex = "Negative")
		self._generatePlots(resultsDF)

def main ():
	VisualizeTracks1 = VisualizeTracks(outDir = "/mnt/belinda_local/venkata/data/Project_Human_RBM17_HEK/TEST_VISUALIZATIONS",
		outPrefix = "TEST_",
		gtf = "/mnt/belinda_local/venkata/data/Index_Files/Human/GenomeFasta_GTF/gencode.v33.primary_assembly.annotation.gtf",
		polyAResults = "/mnt/belinda_local/venkata/data/Project_Human_RBM17_HEK/PolyAMiner_Results/CtrlvsRBM17siRNA_3UTROnly_SoftClipped+Annotations_Run4/3UTROnly_PolyA-miner.Results.txt",
		condition1SamplesBW = "/mnt/belinda_local/venkata/data/Project_Human_RBM17_HEK/BW/HZ8169_.sorted.bw,/mnt/belinda_local/venkata/data/Project_Human_RBM17_HEK/BW/HZ8170_.sorted.bw,/mnt/belinda_local/venkata/data/Project_Human_RBM17_HEK/BW/HZ8171_.sorted.bw",
		condition2SamplesBW = "/mnt/belinda_local/venkata/data/Project_Human_RBM17_HEK/BW/HZ8162_.sorted.bw,/mnt/belinda_local/venkata/data/Project_Human_RBM17_HEK/BW/HZ8163_.sorted.bw,/mnt/belinda_local/venkata/data/Project_Human_RBM17_HEK/BW/HZ8164_.sorted.bw",
		numTop = 100
		)

	VisualizeTracks1.visualizeTopDAGs()

if __name__ == "__main__":
	main()