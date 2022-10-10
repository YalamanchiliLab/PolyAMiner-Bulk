import os, sys
import pandas as pd
import numpy as np
from pyfasta import Fasta
import gtfparse as gp
import torch
from transformers import BertTokenizer, BertModel, DNATokenizer
import seaborn as sns
import matplotlib.pyplot as plt

class CPASBERT:
	def __init__(self, outDir, outPrefix, modelPath, fasta, gtf, predictions, denovoapasiteBED, geneBED, polyADB, polyASite):
		self.outDir = outDir.rstrip("/")+"/"
		self.outPrefix = outPrefix

		self.MODEL_PATH = modelPath
		self.DATA_PATH = self.outDir
		self.PREDICTION_PATH = self.outDir
		self.FASTA = fasta
		self.GTF = gtf
		self.PREDICTIONS = predictions 

		self.GENE_BED = geneBED
		self.POLYADB = polyADB
		self.POLYASITE = polyASite
		self.DENOVOAPASITES_BED = denovoapasiteBED
		self.DENOVOAPASITES_SEQUENCE_ATLAS = self.outDir + self.outPrefix + '_denovoAPAsitesSeqAtlas.tsv'
		self.DEVFILE = self.outDir + "dev.tsv"
		self.PREDRESULTS =self.outDir + "pred_results.npy"

		# self.DENOVOAPASITES_BED = "/mnt/localstorage/venkata/data/Project_AD_YLAB/PolyAMiner_Output/CtrlvsADCerebellum_AllRegionsPAusage/AllRegions_PAusage_denovoAPAsites.bed"

	def _makeGeneBed(self):
		df = gp.read_gtf(self.GTF)
		df = df[df['seqname'].astype(str).str.contains('chr')]
		
		# Gene table in bed format #
		df_gene = df[df['feature'] == 'gene']
		df_gene = df_gene[['seqname', 'start', 'end', 'gene_id', 'gene_name', 'strand']]
		df_gene = df_gene.sort_values(["seqname", "start", "end"], ascending = (True, True,True))
		
		df_gene.to_csv(self.GENE_BED, sep='\t', index=False, header=False)

	def _convertBEDtoSequenceFile(self):
		DENOVOAPASITES_BED_DF = pd.read_csv(self.DENOVOAPASITES_BED, sep = "\t", header = None, names=['Chr', 'Start', 'End', 'gene_id', 'Feature', 'Strand'])
		DENOVOAPASITES_BED_DF.to_csv(self.outDir + self.outPrefix + '_denovoAPAsitesORIGINAL.bed', sep="\t", header=False, index=False)

		DENOVOAPASITES_BED_DF["Midpoint"] = DENOVOAPASITES_BED_DF[['Start', 'End']].mean(axis=1)
		DENOVOAPASITES_BED_DF['Midpoint'] = DENOVOAPASITES_BED_DF['Midpoint'].apply(lambda x: round(x, 0))
		DENOVOAPASITES_BED_DF = DENOVOAPASITES_BED_DF.astype({"Midpoint": int})
		DENOVOAPASITES_BED_DF["Start"] = DENOVOAPASITES_BED_DF["Midpoint"] - 1
		DENOVOAPASITES_BED_DF["End"] = DENOVOAPASITES_BED_DF["Midpoint"]
		DENOVOAPASITES_BED_DF = DENOVOAPASITES_BED_DF.drop('Midpoint', axis=1)
		DENOVOAPASITES_BED_DF.to_csv(self.outDir + self.outPrefix + '_denovoAPAsitesMidpoint.bed', sep="\t", header=False, index=False)

		cmd = "samtools faidx " + self.FASTA
		print(cmd)
		os.system(cmd)

		FASTA_INDEX = self.FASTA.replace(".fa", ".fa.fai")
		CHROM_SIZES = self.outDir + self.outPrefix + "chrom.sizes"
		cmd = "cut -f 1,2 " + FASTA_INDEX + " > " + CHROM_SIZES
		print(cmd)
		os.system(cmd)

		DENOVOAPASITES_MIDPOINT_BED = self.outDir + self.outPrefix + '_denovoAPAsitesMidpoint.bed'
		df = pd.read_csv(DENOVOAPASITES_MIDPOINT_BED, sep = "\t", header = None, names = ['Chr', 'Start', 'End', 'gene_id', 'Feature', 'Strand'])
		# df = df[df["Chr"].str.contains("_")==False]
		# df = df[df["Chr"].str.contains("chrGL")==False]
		df_ChromSizes = pd.read_csv(CHROM_SIZES, sep = "\t", header = None)
		df = df[df['Chr'].isin(df_ChromSizes[0])]
		df.to_csv(DENOVOAPASITES_MIDPOINT_BED, sep="\t", header=False, index=False)

		cmd = "bedtools slop -i " + DENOVOAPASITES_MIDPOINT_BED + " -g " + CHROM_SIZES + " -b 50 > " + self.DENOVOAPASITES_SEQUENCE_ATLAS
		print(cmd)
		os.system(cmd)

		DENOVOAPASITES_SEQUENCE_ATLAS_DF = pd.read_csv(self.DENOVOAPASITES_SEQUENCE_ATLAS, sep = "\t", header = None, names=['Chr', 'Start', 'End', 'gene_id', 'Feature', 'Strand'])
		DENOVOAPASITES_SEQUENCE_ATLAS_DF['Chr'] = DENOVOAPASITES_SEQUENCE_ATLAS_DF['Chr'].str[3:]
		DENOVOAPASITES_SEQUENCE_ATLAS_DF["Chr"] = "chr" + DENOVOAPASITES_SEQUENCE_ATLAS_DF["Chr"] + " " + DENOVOAPASITES_SEQUENCE_ATLAS_DF["Chr"]
		sequenceList = []
		f = Fasta(self.FASTA)
		for index, row in DENOVOAPASITES_SEQUENCE_ATLAS_DF.iterrows():
			try:
				sequence = f.sequence({'chr': row["Chr"], 'start': row["Start"], 'stop': row["End"], 'strand': row["Strand"]}, one_based = False)
			except:
				sequence = "NOT_FOUND"
				# print(sequence)
			sequenceList.append(sequence)

		DENOVOAPASITES_SEQUENCE_ATLAS_DF["Sequence"] = sequenceList
		DENOVOAPASITES_SEQUENCE_ATLAS_DF["SequenceLen"] = DENOVOAPASITES_SEQUENCE_ATLAS_DF['Sequence'].apply(lambda x: len(x))
		DENOVOAPASITES_SEQUENCE_ATLAS_DF = DENOVOAPASITES_SEQUENCE_ATLAS_DF[DENOVOAPASITES_SEQUENCE_ATLAS_DF['Sequence'] != "ERROR"]
		DENOVOAPASITES_SEQUENCE_ATLAS_DF.to_csv(self.DENOVOAPASITES_SEQUENCE_ATLAS, sep = "\t", index = False)

	def _getKmerSentence(self, original_string, kmer = 1, stride = 1):
		if kmer == -1:
			return original_string

		sentence = ""
		original_string = original_string.replace("\n", "")

		i = 0
		while i < len(original_string)-kmer:
			sentence += original_string[i:i+kmer] + " "
			i += stride

		return sentence[:-1].strip("\"")

	def _convertSequenceFiletoDevFile(self):
		ToBFilteredDB = pd.read_csv(self.DENOVOAPASITES_SEQUENCE_ATLAS, sep = "\t")
		# ToBFilteredDB = ToBFilteredDB.iloc[:20,:] #DeleteThisLine

		kmer = []
		for i in ToBFilteredDB["Sequence"]:
			kmer.append(self._getKmerSentence(i, kmer = 6) + "\t1")
		    
		with open(self.DEVFILE, "w") as fileObj:
			fileObj.write("sequence\tlabel\n")
			for element in kmer:
				fileObj.write(element + "\n")

		cmd = "sed '$ s/.$/0/' " + self.DEVFILE + " > " + self.DEVFILE.replace(".tsv", ".intermediate.tsv")
		print(cmd)
		os.system(cmd)

		cmd = "rm " + self.DEVFILE
		print(cmd)
		os.system(cmd)

		cmd = "mv " + self.DEVFILE.replace(".tsv", ".intermediate.tsv") + " " + self.DEVFILE
		print(cmd)
		os.system(cmd)

	def _predict_CPAS_Probabilty(self):
		libPath = os.path.dirname(os.path.abspath(__file__))
		RUN_FINETUNE_PROGRAM = libPath + "/DNABERT/examples/run_finetune.py "

		cmd = "python " + RUN_FINETUNE_PROGRAM + " --model_type dna --tokenizer_name=dna6 --model_name_or_path " + self.MODEL_PATH + \
		" --task_name dnaprom --do_predict --data_dir " + self.DATA_PATH + \
		" --max_seq_length 510 --per_gpu_pred_batch_size=10 --output_dir " + self.MODEL_PATH + \
		" --predict_dir " + self.PREDICTION_PATH + \
		" --n_process 40 --overwrite_cache"

		print(cmd) #Make sure code is written properly
		os.system(cmd)

	def _load_CPAS_Probabilty(self):
		self.PREDRESULTS = self.PREDICTIONS 

	def _createFiltered_CPAS_BED(self):
		DENOVOAPASITES_SEQUENCE_ATLAS_DF = pd.read_csv(self.DENOVOAPASITES_SEQUENCE_ATLAS, sep = "\t")
		DENOVOAPASITES_SEQUENCE_ATLAS_DF["Chr"] = DENOVOAPASITES_SEQUENCE_ATLAS_DF["Chr"].str[:-2]
		DENOVOAPASITES_SEQUENCE_ATLAS_DF['Chr'] = DENOVOAPASITES_SEQUENCE_ATLAS_DF['Chr'].str.strip()

		predictions = np.load(self.PREDRESULTS)
		DENOVOAPASITES_SEQUENCE_ATLAS_DF['Predictions'] = predictions.tolist()
		# DENOVOAPASITES_SEQUENCE_ATLAS_DF['Predictions'] = DENOVOAPASITES_SEQUENCE_ATLAS_DF['Predictions'].astype(float)
		DENOVOAPASITES_SEQUENCE_ATLAS_DF.to_csv(self.DENOVOAPASITES_SEQUENCE_ATLAS, sep = "\t", index = False)
		originalSize = len(DENOVOAPASITES_SEQUENCE_ATLAS_DF.index)

		DENOVOAPASITES_SEQUENCE_ATLAS_DF = DENOVOAPASITES_SEQUENCE_ATLAS_DF[DENOVOAPASITES_SEQUENCE_ATLAS_DF['Predictions'] > 0.5]

		DENOVOAPASITES_SEQUENCE_ATLAS_DF["Midpoint"] = DENOVOAPASITES_SEQUENCE_ATLAS_DF[['Start', 'End']].mean(axis=1)
		DENOVOAPASITES_SEQUENCE_ATLAS_DF['Midpoint'] = DENOVOAPASITES_SEQUENCE_ATLAS_DF['Midpoint'].apply(lambda x: round(x, 0))
		DENOVOAPASITES_SEQUENCE_ATLAS_DF = DENOVOAPASITES_SEQUENCE_ATLAS_DF.astype({"Midpoint": int})
		DENOVOAPASITES_SEQUENCE_ATLAS_DF["Start"] = DENOVOAPASITES_SEQUENCE_ATLAS_DF["Midpoint"] - 1
		DENOVOAPASITES_SEQUENCE_ATLAS_DF["End"] = DENOVOAPASITES_SEQUENCE_ATLAS_DF["Midpoint"]
		
		DENOVOAPASITES_SEQUENCE_ATLAS_DF = DENOVOAPASITES_SEQUENCE_ATLAS_DF.drop('Midpoint', axis=1)
		DENOVOAPASITES_SEQUENCE_ATLAS_DF = DENOVOAPASITES_SEQUENCE_ATLAS_DF.drop('Predictions', axis=1)
		DENOVOAPASITES_SEQUENCE_ATLAS_DF = DENOVOAPASITES_SEQUENCE_ATLAS_DF.drop('Sequence', axis=1)
		DENOVOAPASITES_SEQUENCE_ATLAS_DF = DENOVOAPASITES_SEQUENCE_ATLAS_DF.drop('SequenceLen', axis=1)
		finalSize = len(DENOVOAPASITES_SEQUENCE_ATLAS_DF.index)

		keptPercentage = finalSize/originalSize * 100 
		keptPercentage = str(round(keptPercentage, 2))
		print("After filtering with CPAS-BERT, " + keptPercentage + "%/ of the original " + self.DENOVOAPASITES_BED +" file remains...")

		DENOVOAPASITES_SEQUENCE_ATLAS_DF.to_csv(self.DENOVOAPASITES_BED, sep = "\t", index = False, header = None)

	def filter_CPAS_Sites(self):
		self._convertBEDtoSequenceFile()
		self._convertSequenceFiletoDevFile()
		self._predict_CPAS_Probabilty()
		self._createFiltered_CPAS_BED()

	def filter_CPAS_Sites_fast(self):
		self._convertBEDtoSequenceFile()
		self._convertSequenceFiletoDevFile()
		self._load_CPAS_Probabilty()
		self._createFiltered_CPAS_BED()

	def _upgradePolyADB(self):
		df = pd.read_csv(self.POLYADB, sep = "\t")
		df = df[['Chromosome', 'Position', "Strand"]]
		df = df.rename(columns={"Position": "End"})
		df["Start"] = df["End"] - 1
		df["Feature1"] = "-"
		df["Feature2"] = "-"
		cols = ["Chromosome", "Start", "End", "Feature1", "Feature2", "Strand"]
		df = df[cols]

		CLEANED_POLYADB= self.outDir + self.outPrefix + os.path.basename(self.POLYADB).replace(".txt", ".bed")
		df.to_csv(CLEANED_POLYADB, sep = "\t", index = False, header = None) #Pass to liftover program in UCSC genome browswer

	def _generatePositiveLabels(self):
		#Make PolyADB and PolyASite 2.0 into similarly structured dataframe with a PAS_ID column
		df1 = pd.read_csv(self.POLYADB, sep = "\t", header = None, names = ["Chromosome", "Start", "End", "Feature1", "Feature2", "Strand"])
		df1['Start'] = df1['Start'].astype(str)
		df1['End'] = df1['End'].astype(str)
		df1["PAS_ID"] = df1["Chromosome"]+":"+df1["Start"]+"-"+df1["End"]+":"+df1["Strand"]
		cols = ['Chromosome', 'Start', 'End', 'Feature1', 'PAS_ID', 'Strand']
		df1 = df1[cols]
		FORMATTED_POLYADB = self.outDir + self.outPrefix + os.path.basename(self.POLYADB).replace(".bed", ".formatted.bed")
		df1.to_csv(FORMATTED_POLYADB, sep = "\t", index = False, header = None)

		df2 = pd.read_csv(self.POLYASITE, sep = "\t", header = None)
		df2 = df2.iloc[: , :6]
		df2.columns = ['Chromosome', 'Start', 'End', 'Feature1', 'Feature2', 'Strand']
		df2['Chromosome'] = "chr"+df2['Chromosome'].astype(str)
		df2["Midpoint"] = df2[['Start', 'End']].mean(axis=1).astype(int)
		df2['Start'] = df2['Midpoint'] - 1
		df2['End'] = df2['Midpoint']
		df2['Start'] = df2['Start'].astype(str)
		df2['End'] = df2['End'].astype(str)
		df2["PAS_ID"] = df2["Chromosome"]+":"+df2["Start"]+"-"+df2["End"]+":"+df2["Strand"]
		cols = ['Chromosome', 'Start', 'End', 'Midpoint', 'PAS_ID', 'Strand']
		df2 = df2[cols]
		FORMATTED_POLYASITE = self.outDir + self.outPrefix + os.path.basename(self.POLYASITE).replace(".bed", ".formatted.bed")
		df2.to_csv(FORMATTED_POLYASITE, sep = "\t", index = False, header = None)

		WINDOWED_CPAS_BED = self.outDir + self.outPrefix + "windowed.CPAS.bed"
		cmd = "bedtools window -a " + FORMATTED_POLYADB + " -b " + FORMATTED_POLYASITE + " -w 24 -sm > " + WINDOWED_CPAS_BED
		print(cmd)
		os.system(cmd)

		df = pd.read_csv(WINDOWED_CPAS_BED, sep = "\t", header = None, names = ["AChromosome", "AStart", "AEnd", "AFeature1", "AFeature2", "AStrand", "BChromosome", "BStart", "BEnd", "BFeature1", "BFeature2", "BStrand"])
		PAS_ID1 = df["AFeature2"].tolist()
		PAS_ID2 = df["BFeature2"].tolist()
		PAS_ID = PAS_ID1 + PAS_ID2
		PAS_ID = [*set(PAS_ID)]

		df = pd.DataFrame(PAS_ID, columns = ["PAS_ID"])
		df[["Chromosome", "Positions", "Strand"]] = df['PAS_ID'].str.split(":", expand=True)
		df[["Start", "End"]] = df['Positions'].str.split("-", expand=True)
		df = df.drop(['Positions'], axis=1)
		df["Feature1"] = "-"
		cols = ["Chromosome", "Start", "End", "Feature1","PAS_ID","Strand"]
		df = df[cols]
		EXPANEDED_WINDOWED_CPAS_BED = self.outDir + self.outPrefix + "expanded.windowed.CPAS.bed"
		df.to_csv(EXPANEDED_WINDOWED_CPAS_BED, sep = "\t", index = False, header = None)

		FASTA_INDEX = self.FASTA.replace(".fa", ".fa.fai")
		CHROM_SIZES = self.outDir + self.outPrefix + "chrom.sizes"
		cmd = "cut -f 1,2 " + FASTA_INDEX + " > " + CHROM_SIZES
		print(cmd)
		os.system(cmd)

		SLOPPED_EXPANEDED_WINDOWED_CPAS_BED = self.outDir + self.outPrefix + "slopped.expanded.windowed.CPAS.bed"
		cmd = "bedtools slop -b 50 -i " + EXPANEDED_WINDOWED_CPAS_BED + " -g " + CHROM_SIZES +  " > " + SLOPPED_EXPANEDED_WINDOWED_CPAS_BED
		print(cmd)
		os.system(cmd)

		SEQ_SLOPPED_EXPANEDED_WINDOWED_CPAS_BED = self.outDir + self.outPrefix + "sequence.slopped.expanded.windowed.CPAS.bed"
		cmd = "bedtools getfasta -s -bedOut -fi " + self.FASTA + " -bed " + SLOPPED_EXPANEDED_WINDOWED_CPAS_BED + " > " + SEQ_SLOPPED_EXPANEDED_WINDOWED_CPAS_BED
		print(cmd)
		os.system(cmd) 

		df = pd.read_csv(SEQ_SLOPPED_EXPANEDED_WINDOWED_CPAS_BED, sep = "\t", header = None, names = ["Chr", "Start", "End", "Feature1", "Feature2", "Strand", "Sequence"])
		df["Label"] = 1
		df = df[df["Sequence"].str.contains("N")==False]
		df.to_csv(SEQ_SLOPPED_EXPANEDED_WINDOWED_CPAS_BED, sep = '\t', index = False, header = None)
		return df


	def _generateNegativeLabels(self):
		self._makeGeneBed()

		FASTA_INDEX = self.FASTA.replace(".fa", ".fa.fai")
		CHROM_SIZES = self.outDir + self.outPrefix + "chrom.sizes"
		cmd = "cut -f 1,2 " + FASTA_INDEX + " > " + CHROM_SIZES
		print(cmd)
		os.system(cmd)

		SLOPPED_GENE_BED = self.outDir + self.outPrefix + os.path.basename(self.GENE_BED).replace(".bed", ".slopped.bed")
		cmd = "bedtools slop -b 3000 -i " + self.GENE_BED + " -g " + CHROM_SIZES + " > " + SLOPPED_GENE_BED
		print(cmd)
		os.system(cmd)

		RANDOM_BED = self.outDir + self.outPrefix + '_random101bpRegions.bed'
		cmd = "bedtools random -l 101 -n 525000 -seed 2022 -g " + CHROM_SIZES + " > " + RANDOM_BED
		print(cmd)
		os.system(cmd)

		FILTERED_RANDOM_BED = self.outDir + self.outPrefix + '_filtered_random101bpRegions.bed'
		cmd = "bedtools subtract -s -A -a " + RANDOM_BED + " -b " + SLOPPED_GENE_BED + " > " + FILTERED_RANDOM_BED
		print(cmd)
		os.system(cmd)

		SEQ_FILTERED_RANDOM_BED = self.outDir + self.outPrefix + '_sequence_filtered_random101bpRegions.bed'
		cmd = "bedtools getfasta -s -bedOut -fi " + self.FASTA + " -bed " + FILTERED_RANDOM_BED + " > " + SEQ_FILTERED_RANDOM_BED
		print(cmd)
		os.system(cmd) 

		df = pd.read_csv(SEQ_FILTERED_RANDOM_BED, sep = "\t", header = None, names = ["Chr", "Start", "End", "Feature1", "Feature2", "Strand", "Sequence"])
		df["Label"] = 0
		df = df[df["Sequence"].str.contains("N")==False]
		df.to_csv(SEQ_FILTERED_RANDOM_BED, sep = '\t', index = False, header = None)
		return df

	def _getKmerSentence(self, original_string, kmer = 1, stride = 1):
		if kmer == -1:
			return original_string

		sentence = ""
		original_string = original_string.replace("\n", "")

		i = 0
		while i < len(original_string)-kmer:
			sentence += original_string[i:i+kmer] + " "
			i += stride

		return sentence[:-1].strip("\"")

	def generateTrainingData(self):
		posDF = self._generatePositiveLabels()
		negDF = self._generateNegativeLabels()

		trainingDF = pd.concat([posDF,negDF],ignore_index=True)
		trainingDF = trainingDF.sample(frac=1)

		TRAINING_DF = self.outDir + self.outPrefix + '_trainingDF.bed'
		trainingDF.to_csv(TRAINING_DF, sep = '\t', index = False)
		trainingDF['Label'] = trainingDF['Label'].astype(str)
		kmer = []
		for index, row in trainingDF.iterrows():
			kmer.append(self._getKmerSentence(row["Sequence"], kmer = 6) + "\t"+ row["Label"])

		num_kmer_test = int(len(kmer)*0.1)
		kmer_test = kmer[:num_kmer_test]

		train = kmer[num_kmer_test:]
		test = kmer_test

		trainFile = self.outDir + "train.tsv"
		with open(trainFile, "w") as fileObj:
			fileObj.write("sequence\tlabel\n")
			for element in train:
				fileObj.write(element + "\n")

		devFile = self.outDir + "dev.tsv"
		with open(devFile, "w") as fileObj:
			fileObj.write("sequence\tlabel\n")
			for element in test:
				fileObj.write(element + "\n")

	def trainModel(self):
		libPath = os.path.dirname(os.path.abspath(__file__))
		PRETRAINED_DNABERT6 = libPath + "/6-new-12w-0"
		RUN_FINETUNE_PROGRAM = libPath + "/DNABERT/examples/run_finetune.py "
		self.DATA_PATH = "/mnt/atlas_local/venkata/data/PolyA-miner_v1.5/CPAS_BERT_TESTING/Human_hg38 "

		cmd = "python " + RUN_FINETUNE_PROGRAM + \
		" --model_type dna --tokenizer_name=dna6 --model_name_or_path " + PRETRAINED_DNABERT6 + " --task_name dnaprom " + \
		" --do_train --do_eval --data_dir " + self.DATA_PATH + \
		" --max_seq_length 510 --per_gpu_eval_batch_size=10 --per_gpu_train_batch_size=10 --num_train_epochs 10 --output_dir " + self.DATA_PATH + \
		" --evaluate_during_training --logging_steps 100 --save_steps 4000 --warmup_percent 0.1 --hidden_dropout_prob 0.1 --overwrite_output " + \
		" --weight_decay 0.01 --n_process 8 --overwrite_cache --should_continue"

		print(cmd) #Make sure code is written properly
		os.system(cmd)

	def _formatAttention(self, attention):
		squeezed = []
		for layer_attention in attention:
			# 1 x num_heads x seq_len x seq_len
			if len(layer_attention.shape) != 4:
				raise ValueError("The attention tensor does not have the correct number of dimensions. Make sure you set "
					"output_attentions=True when initializing your model.")
			squeezed.append(layer_attention.squeeze(0))
		# num_layers x num_heads x seq_len x seq_len
		return torch.stack(squeezed)

	def _getAttentionDNA(self, model, tokenizer, sentence_a, start, end):
		inputs = tokenizer.encode_plus(sentence_a, sentence_b=None, return_tensors='pt', add_special_tokens=True)
		input_ids = inputs['input_ids']
		attention = model(input_ids)[-1]
		input_id_list = input_ids[0].tolist() # Batch index 0
		tokens = tokenizer.convert_ids_to_tokens(input_id_list) 
		attn = self._formatAttention(attention)
		attn_score = []
		for i in range(1, len(tokens)-1):
			attn_score.append(float(attn[start:end+1,:,0,i].sum()))
		return attn_score

	def _getRealScore(self, attention_scores, kmer, metric):
		counts = np.zeros([len(attention_scores)+kmer-1])
		real_scores = np.zeros([len(attention_scores)+kmer-1])

		if metric == "mean":
			for i, score in enumerate(attention_scores):
				for j in range(kmer):
					counts[i+j] += 1.0
					real_scores[i+j] += score

			real_scores = real_scores/counts
		else:
			pass

		return real_scores

	def generateAttentionHeatmap(self, SequenceBEDFileLoc, imagePrefix, cmapPalette):
		tokenizer_name = 'dna6'
		model = BertModel.from_pretrained(self.MODEL_PATH, output_attentions=True)
		tokenizer = DNATokenizer.from_pretrained(tokenizer_name, do_lower_case=False)

		#might have to play around with the size of master_scores; either self.intervalLength, self.intervalLength-1, or self.intervalLength+1
		master_scores = np.empty((0, 100), int)

		SequenceBed = pd.read_csv(SequenceBEDFileLoc, sep = "\t", header = None, names = ["Chr", "Start", "End", "rnd", "Length", "Strand", "Sequence", "Label"])
		SequenceBed = SequenceBed.loc[SequenceBed['Strand'] == "+"]
		SequenceList = SequenceBed.Sequence.values.tolist()
		LabelList = SequenceBed.Label.values.tolist()
		SequenceList = SequenceList[0:300000]
		print(LabelList[0])

		x = 0
		for raw_sentence in SequenceList:
			raw_sentence = str(raw_sentence)
			sentence_a = self._getKmerSentence(raw_sentence, kmer = 6)
			tokens = sentence_a.split()

			#might have to change args.start_layer and args.end_layer parameter
			attention = self._getAttentionDNA(model, tokenizer, sentence_a, start=9, end=12)
			print(attention)
			attention_scores = np.array(attention).reshape(np.array(attention).shape[0],1)

			real_scores = self._getRealScore(attention_scores, kmer = 6, metric = "mean")
			scores = real_scores.reshape(1, real_scores.shape[0])

			master_scores = np.concatenate((master_scores, scores))
			x += 1
			print(str(x) + " sequences analyzed!")

		master_scores_npy_saveLoc = self.outDir + "master_attention_scores.npy"
		np.save(master_scores_npy_saveLoc, master_scores)

		# plot        
		sns.set()
		# ax = sns.heatmap(master_scores, cmap='YlGnBu', vmin=0, vmax=2)
		# ax = sns.clustermap(master_scores, cmap='YlGnBu', vmin=0, vmax=2, row_cluster=True, col_cluster = False, yticklabels=False)
		ax = sns.clustermap(master_scores, cmap=cmapPalette, vmin=0, vmax=2, row_cluster=True, col_cluster = False, yticklabels=False)
		ax.ax_row_dendrogram.set_visible(False)
		ax.ax_row_dendrogram.set_xlim([0,0])
		x0, y0, w, h = ax.cbar_pos
		ax.ax_cbar.set_position([x0, 0.4, w, h])

		# x_values = []
		# x= -50
		# while x < 50:
		# 	x_values.append(x)
		# 	x += 5
		# ax.set_xticklabels(x_values)
		plt.show()
		plt.savefig(self.outDir + imagePrefix + "heatmap_visualization.png")
		plt.figure().clear()
		plt.close()
		plt.cla()
		plt.clf()

def main():
	libPath = os.path.dirname(os.path.abspath(__file__)) + "/CPASBERT_TrainedModels"
	modelPath = libPath + "/hg38_checkpoint-64000"
	print(modelPath)
	CPASBERT1 = CPASBERT(outDir = "/mnt/belinda_local/venkata/data/PolyAMiner-Bulk/TestFiles_Human_APriori",
		outPrefix = "Human_hg38_", 
		modelPath = modelPath,
		fasta = "",
		gtf = "",
		predictions = "",
		denovoapasiteBED = "",
		geneBED = "",
		polyADB = "",
		polyASite = ""
		)


	NegLabeledSequenceBedLoc = "/mnt/belinda_local/venkata/data/PolyAMiner-Bulk/TestFiles_Human_APriori/Human_hg38_NegLabeledSequence.bed"
	CPASBERT1.generateAttentionHeatmap(NegLabeledSequenceBedLoc, "NegSequence", cmapPalette = "YlOrRd")

	PosLabeledSequenceBedLoc = "/mnt/belinda_local/venkata/data/PolyAMiner-Bulk/TestFiles_Human_APriori/Human_hg38_PosLabeledSequence.bed"
	CPASBERT1.generateAttentionHeatmap(PosLabeledSequenceBedLoc, "PosSequence", cmapPalette = "YlGnBu")
	# CPASBERT1.trainModel()

	cmd = "python3 /mnt/belinda_local/venkata/data/venkata/DNABERT/motif/find_motifs.py  --min_len 6 --window_size 6  --align_all_ties  --save_file_dir /mnt/belinda_local/venkata/data/PolyAMiner-Bulk/TestFiles_Human_APriori --predict_dir /mnt/belinda_local/venkata/data/PolyAMiner-Bulk/TestFiles_Human_APriori  --data_dir /mnt/belinda_local/venkata/data/PolyAMiner-Bulk/TestFiles_Human_APriori --verbose --min_n_motif 3"
	print(cmd)
	os.system(cmd)

if __name__ == "__main__":
	main()

