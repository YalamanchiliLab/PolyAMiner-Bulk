import os, sys
import pandas as pd

class DEGAnalyzer:
	def __init__(self, outDir, outPrefix, bed, gtf, condition1Samples, condition2Samples):
		self.outDir = outDir
		self.outPrefix = outPrefix
		self.BED = bed
		self.GTF = gtf

		self.APACOUNTMATRIX4DEGS = outDir + outPrefix + '_APACountMatrix4DEGs.txt'
		self.APACOUNTMATRIX4DEGSDF = pd.read_csv(self.APACOUNTMATRIX4DEGS,comment='#',sep="\t",index_col=None)

		self.DESIGNFILE = self.outDir+self.outPrefix+"_Design.txt"
		self.APACOUNTMATRIXFORPCA = self.APACOUNTMATRIX4DEGS.replace(".txt","_APA_PCA.txt")
		self.APACOUNTMATRIX_GENELEVEL = self.APACOUNTMATRIX4DEGS.replace(".txt","_Genes.txt")

		self.OVERALL_DESEQ2_RESULTS_FILE = outDir + outPrefix + "_DEG-results.txt"
		self.SIMPLE_OVERALL_DESEQ2_RESULTS_FILE = outDir + outPrefix + "_simple_DEG-results.txt"
		self.APAFACTOR_DESEQ2_RESULTS_FILE = outDir + outPrefix + "_APAFACTOR_DEG-results.txt"

		self.OVERALL_HEATMAP_FILE = outDir + outPrefix + "_overallDEGHeatmap.pdf"
		self.OVERALL_VOLCANOPLOT_FILE = outDir + outPrefix + "_overallDEGVolcanoPlot.tiff"

		self.APAFACTOR_HEATMAP_FILE = outDir + outPrefix + "_APAFactor_DEGHeatmap.pdf"
		self.APAFACTOR_VOLCANOPLOT_FILE = outDir + outPrefix + "_APAFactor_DEGVolcanoPlot.tiff"

		self.condition1Samples = condition1Samples
		self.condition2Samples = condition2Samples
		self.cond1NumSamples = self.condition1Samples[0].count(",") + 1
		self.cond2NumSamples = self.condition2Samples[0].count(",") + 1
		self.conditionSamplesList = condition1Samples + condition2Samples
		self.conditionSamples = " ".join(str(e) for e in self.conditionSamplesList)
		self.conditionSamples = self.conditionSamples.replace(",", " ")

		return None

	def _createDesignFile(self):
		# Write design file -APA -DEG #
		fw=open(self.DESIGNFILE,"w")
		fw.write("\tGenotype\n")
		for s in self.APACOUNTMATRIX4DEGSDF.columns[2:2+self.cond1NumSamples]:
			fw.write(s+"\tCR\n")
		for s in self.APACOUNTMATRIX4DEGSDF.columns[self.cond1NumSamples+2:self.cond1NumSamples+2+self.cond2NumSamples]:
			fw.write(s+"\tTR\n")
		fw.close()

	def _plotPCAandTSNE(self, inputMatrixFile, outputFileName):
		ClusteringRScriptLoc = os.path.dirname(__file__) + "/PolyA-miner_Clustering.r"
		os.system("Rscript " + ClusteringRScriptLoc + " " + inputMatrixFile + " " + self.DESIGNFILE + " " + outputFileName)

	def _getGeneLevelCounts(self):
		gdf=self.APACOUNTMATRIX4DEGSDF.drop(columns=['feature_id']).copy()
		gdf[gdf.columns[1:]]=gdf[gdf.columns[1:]].apply(pd.to_numeric)
		gdf = gdf.groupby(['gene_id'], as_index=False, sort=False).aggregate('sum')
		gdf=gdf[(self.APACOUNTMATRIX4DEGSDF.gene_id != ".")]
		gdf.to_csv(self.APACOUNTMATRIX_GENELEVEL, sep="\t", index=False)

	def _calculateDEGs(self, inputMatrixFile):
		DEGTestRScriptLoc = os.path.dirname(__file__) + "/PolyA-miner_DEGtest.r"
		os.system("Rscript " + DEGTestRScriptLoc + " " + inputMatrixFile + " " + self.DESIGNFILE + " " + self.outDir + self.outPrefix)
	
	def _generateHeatmap(self, DESEQ2ResultsFile, heatmapFileName):
		endPos = str(9 + self.cond1NumSamples + self.cond2NumSamples)
		cmd = "cat " + DESEQ2ResultsFile + " | cut -f 1,9-" + endPos + " > " + self.SIMPLE_OVERALL_DESEQ2_RESULTS_FILE
		os.system(cmd)

		DrawHeatmapRScriptLoc = os.path.dirname(__file__) + "/draw-heatmap.r"
		cmd = "cat " + self.SIMPLE_OVERALL_DESEQ2_RESULTS_FILE + " | Rscript " + DrawHeatmapRScriptLoc + " > " + heatmapFileName
		os.system(cmd)

	def _generateDEGVolcanoPlot(self, DESEQ2ResultsFile, volcanoPlotFileName, plotTitle):
		data=pd.read_csv(DESEQ2ResultsFile,sep="\t",header=0,index_col=None)
		
		genes=pd.read_csv(self.BED,sep="\t",header=None,index_col=None)
		genes.columns=["Chr","Start","End","Gene","Symbol","Strand"]
		if "." in data['Gene'].iloc[0]:
			data=pd.merge(data,genes,on=["Gene"],how="left")
			data=data.drop(columns=["Chr","Start","End","Strand"])
		else:
			try:
				genes[['Gene','Version']] = genes['Gene'].str.split('.',expand=True)
				data=pd.merge(data,genes,on=["Gene"],how="left")
				data=data.drop(columns=["Chr","Start","End","Strand","Version"])
			except:
				data=pd.merge(data,genes,on=["Gene"],how="left")
				data=data.drop(columns=["Chr","Start","End","Strand"])

		DESEQ2ResultsFileSymbol = DESEQ2ResultsFile[:len(DESEQ2ResultsFile)-4]+"_Symbol.txt"
		data.to_csv(DESEQ2ResultsFileSymbol, index=False, sep ="\t")

		temp1 = data[data['padj']<=0.05]
		temp2 = data[(data['padj']<=0.05) & ((data['log2FoldChange']>=0.263) | (data['log2FoldChange']<=-0.263))]
		sorteddf=temp2.sort_values(by=['log2FoldChange'])
		
		newlist = [x for x in list(sorteddf.head(5)['Symbol']) if pd.isnull(x) == False and x != 'nan']
		x="_".join(newlist)
		#print(x)
		#print(sorteddf.head(5))
		#print(sorteddf.tail(5))
		
		newlist = [x for x in list(sorteddf.tail(5)['Symbol']) if pd.isnull(x) == False and x != 'nan']
		y="_".join(newlist)
		#print(y)
		temp3 = data[(data['padj']<=0.05) & (data['log2FoldChange']>=0.263)]
		temp4 = data[(data['padj']<=0.05) & (data['log2FoldChange']<=-0.263)]
		
		ndeg=str(temp1.shape[0])
		fcdeg=str(temp2.shape[0])
		updeg=str(temp3.shape[0])
		dndeg=str(temp4.shape[0])

		DEGVolcanoPlotRScriptLoc = os.path.dirname(__file__) + "/DEGvolcanoPlot.r"
		os.system("Rscript " + DEGVolcanoPlotRScriptLoc + " "+DESEQ2ResultsFileSymbol+" "+volcanoPlotFileName+" "+ndeg+" "+fcdeg+" "+updeg+" "+dndeg+" "+x+" "+y+ " "+plotTitle)
		# os.system("Rscript DEGvolcanoPlot.r "+degfile+" "+degfile.split("_DEG")[0]+" "+ndeg+" "+fcdeg+" "+updeg+" "+dndeg+" "+x+" "+y)
		

	def _filterAPAFactors(self):
		OverallDEGResultsDF = pd.read_csv(self.OVERALL_DESEQ2_RESULTS_FILE, sep = "\t")
		OverallDEGResultsDF['Gene'] = OverallDEGResultsDF.Gene.map(lambda x: x[0: x.find('.')] if '.' in x else x)

		Ensembl_Human_APAFactorList = ["ENSG00000071894","ENSG00000165934","ENSG00000119203","ENSG00000160917","ENSG00000136709","ENSG00000145216","ENSG00000101138","ENSG00000101811","ENSG00000177613","ENSG00000176102","ENSG00000167005","ENSG00000149532","ENSG00000111605","ENSG00000125755","ENSG00000165494","ENSG00000172409","ENSG00000090060","ENSG00000115421","ENSG00000172531","ENSG00000213639","ENSG00000070756","ENSG00000090621","ENSG00000122257","ENSG00000100836"]
		OverallDEGResultsDF = OverallDEGResultsDF[OverallDEGResultsDF['Gene'].isin(Ensembl_Human_APAFactorList)]

		if OverallDEGResultsDF.empty:
			print("Entered filtering overall DEG Results by ensemble mouse APA Factor List")
			OverallDEGResultsDF = pd.read_csv(self.OVERALL_DESEQ2_RESULTS_FILE, sep = "\t")
			OverallDEGResultsDF['Gene'] = OverallDEGResultsDF.Gene.map(lambda x: x[0: x.find('.')] if '.' in x else x)
			Ensembl_Mouse_APAFactorList = ["ENSMUSG00000034022", "ENSMUSG00000041781", "ENSMUSG00000054309", "ENSMUSG00000029625","ENSMUSG00000024400", "ENSMUSG00000029227", "ENSMUSG00000027498", "ENSMUSG00000031256", "ENSMUSG00000053536", "ENSMUSG00000027176", "ENSMUSG00000031754", "ENSMUSG00000034820", "ENSMUSG00000055531", "ENSMUSG00000023118", "ENSMUSG00000041328", "ENSMUSG00000027079", "ENSMUSG00000021111", "ENSMUSG00000020273", "ENSMUSG00000040385", "ENSMUSG00000014956", "ENSMUSG00000022283", "ENSMUSG00000011257", "ENSMUSG00000030779", "ENSMUSG00000022194"]
			OverallDEGResultsDF = OverallDEGResultsDF[OverallDEGResultsDF['Gene'].isin(Ensembl_Mouse_APAFactorList)]
		if OverallDEGResultsDF.empty:
			print("Entered filtering overall DEG Results by symbol mouse APA Factor List")
			OverallDEGResultsDF = pd.read_csv(self.OVERALL_DESEQ2_RESULTS_FILE, sep = "\t")
			OverallDEGResultsDF['Gene'] = OverallDEGResultsDF.Gene.map(lambda x: x[0: x.find('.')] if '.' in x else x)
			Symbol_Mouse_APAFactorList = ["Cpsf1", "Cpsf2", "Cpsf3", "Cpsf4", "Wdr33", "Fip1l1", "Cstf1", "Cstf2", "Cstf2t", "Cstf3", "Nudt21", "Cpsf7", "Cpsf6", "Sympk", "Pcf11", "Clp1" "Papola", "Papolg", "Ppp1ca", "Ppp1cb", "Pabpc1", "Pabpc4", "Rbbp6", "Pabpn1"]
			OverallDEGResultsDF = OverallDEGResultsDF[OverallDEGResultsDF['Gene'].isin(Symbol_Mouse_APAFactorList)]
		if OverallDEGResultsDF.empty:
			print("Entered filtering overall DEG Results by symbol human APA Factor List")
			OverallDEGResultsDF = pd.read_csv(self.OVERALL_DESEQ2_RESULTS_FILE, sep = "\t")
			OverallDEGResultsDF['Gene'] = OverallDEGResultsDF.Gene.map(lambda x: x[0: x.find('.')] if '.' in x else x)
			Symbol_Human_APAFactorList = ["CPSF1", "CPSF2", "CPSF3", "CPSF4", "WDR33", "FIP1L1", "CSTF1", "CSTF2", "CSTF2T", "CSTF3", "NUDT21", "CPSF7", "CPSF6", "SYMPK", "PCF11", "CLP1", "PAPOLA", "PAPOLG", "PPP1CA", "PPP1CB", "PABPC1", "PABPC4", "RBBP6", "PABPN1"]
			OverallDEGResultsDF = OverallDEGResultsDF[OverallDEGResultsDF['Gene'].isin(Symbol_Human_APAFactorList)]

		OverallDEGResultsDF.to_csv(self.APAFACTOR_DESEQ2_RESULTS_FILE, index=False, sep ="\t")

		if OverallDEGResultsDF.empty:
			print("Dataset not sufficient depth")
			exit()
			return False

		return True

	def clusterAndAnalyzeDEGs(self):
		self._createDesignFile()

		# Prepare APA counts for PCA #
		adf=self.APACOUNTMATRIX4DEGSDF.drop(columns=['gene_id']).copy()
		adf.to_csv(self.APACOUNTMATRIXFORPCA,sep="\t", index=False)

		self._getGeneLevelCounts()

		# Plot PCA & t-SNE #
		self._plotPCAandTSNE(inputMatrixFile = self.APACOUNTMATRIXFORPCA, outputFileName = self.outDir+self.outPrefix+"_PA")
		self._plotPCAandTSNE(inputMatrixFile = self.APACOUNTMATRIX_GENELEVEL, outputFileName = self.outDir+self.outPrefix+"_Gene")

		#Overall DEGs
		self._calculateDEGs(inputMatrixFile = self.APACOUNTMATRIX_GENELEVEL)
		self._generateHeatmap(DESEQ2ResultsFile = self.OVERALL_DESEQ2_RESULTS_FILE, heatmapFileName = self.OVERALL_HEATMAP_FILE)
		self._generateDEGVolcanoPlot(DESEQ2ResultsFile = self.OVERALL_DESEQ2_RESULTS_FILE, volcanoPlotFileName = self.OVERALL_VOLCANOPLOT_FILE, plotTitle = "DEG_Volcano_Plot")

		#Core APA Factor DEGs
		if (self._filterAPAFactors()):
			self._generateHeatmap(DESEQ2ResultsFile = self.APAFACTOR_DESEQ2_RESULTS_FILE, heatmapFileName = self.APAFACTOR_HEATMAP_FILE)
			self._generateDEGVolcanoPlot(DESEQ2ResultsFile = self.APAFACTOR_DESEQ2_RESULTS_FILE, volcanoPlotFileName = self.APAFACTOR_VOLCANOPLOT_FILE, plotTitle = "APA_Factor_Volcano_Plot")
		