# v1.2
# Extract PolyA sites from soft-clipping and CPAS-BERT Filtered Databases #

import uuid, math
import pandas as pd,subprocess
from pybedtools import BedTool
import pybedtools as pb, gtfparse as gp
import os, sys, glob, time, argparse, pysam
import concurrent.futures as cf
import numpy as np
from CPASBERT import CPASBERT

class ExtractPolyAsites4Bulk:
	def __init__(self, outDir, outPrefix, fasta, gtf, con1BAMFiles, con2BAMFiles, proportionA, modelOrganism, apriori_annotations, ignoreFeatures):
		self.outDir = outDir.rstrip("/")+"/"
		self.outPrefix = outPrefix

		self.apriori_annotations = apriori_annotations
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

		self.FASTA = fasta
		self.GTF = gtf

		self.con1BAMFiles = "".join(con1BAMFiles).replace(" ","").split(",")
		self.con2BAMFiles = "".join(con2BAMFiles).replace(" ","").split(",")
		self.conditionBAMList = self.con1BAMFiles + self.con2BAMFiles

		self.proportionA = proportionA
		self.slopDistanceParameter = "25"
		#Proportion of A's in true polyA tails. Default: '0.90,0.85,0.80,0.75' ->0.90 (for a tail of 4nt long), 0.85 (for a tail of 4-8nt long), 0.80 (for a tail of 8-12nt long), 0.75 (for a tail > 12nt long) ''',nargs='+',type=str,default= "0.90,0.85,0.80,0.75"
		
		self.GENE_BED = self.outDir + self.outPrefix + "Genes.bed"
		self.CDS_BED = self.outDir + self.outPrefix + "CDS.bed"
		self.UTR3_BED = self.outDir + self.outPrefix + "UTR3.bed"
		self.UTR5_BED = self.outDir + self.outPrefix + "UTR5.bed"
		self.INTRON_BED = self.outDir + self.outPrefix + "Intron.bed"
		self.FEATURES_BED = self.outDir + self.outPrefix + "Features.bed"

		self.FINAL_SOFTCLIPPED_BED_FILE = self.outDir + self.outPrefix + "Final.SoftClipped.bed"
		self.FINAL_SAF_FILE = self.outDir + self.outPrefix + "_denovoAPAsites.saf"

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
		cmd = 'bedtools closest -nonamecheck -a ' + self.FINAL_SOFTCLIPPED_BED_FILE + ' -b ' + self.GENE_BED + ' -s -id -D a -t first -k 1 > ' + TEMPFILE
		os.system(cmd)
		tempdf = pd.read_csv(TEMPFILE, sep='\t', index_col=None, header=None)
		tempdf = tempdf.iloc[:, [0, 1, 2, 5, 10, 13]]
		tempdf.columns = ['Chr', 'Start', 'End', 'Strand', 'Gene', 'Distance']
		tempdf = tempdf[['Chr', 'Start', 'End', 'Gene', 'Distance', 'Strand']]
		# os.system('rm ' + TEMPFILE)
		tempdf.to_csv(self.FINAL_SOFTCLIPPED_BED_FILE, sep='\t', index=False, header=None)

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

	def _mapAPA2Features(self):
		APAfile = self.FINAL_SOFTCLIPPED_BED_FILE

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

	def _compute(self, string, key):
		string = string.strip().upper()
		return(string.count(key)/float(len(string)))

	def _truePolyA(self):
		prop = "".join(self.proportionA).replace(" ","").split(",")

		for BAM_FILE in self.conditionBAMList:
			bamfile = pysam.AlignmentFile(BAM_FILE, 'rb')
			CAPPED_BAM_FILE = self.outDir + self.outPrefix + os.path.basename(BAM_FILE).replace(".bam", ".capped.bam")
			CAPPED_BED_FILE = self.outDir + self.outPrefix + os.path.basename(BAM_FILE).replace(".bam", ".capped.bed")

			BED_object = open(CAPPED_BED_FILE, "w")

			polyAreads = pysam.AlignmentFile(CAPPED_BAM_FILE, "wb", template=bamfile)
			for read in bamfile:
				if (not read.is_unmapped) and ("S" in read.cigarstring):
					if (not read.is_reverse):
						#print(read.cigartuples)
						#print(read.cigarstring)
						#print(read.cigartuples[-1])
						#print(read.query_sequence)
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

			# index #
			os.system("samtools index "+CAPPED_BAM_FILE)

			cmd = "samtools faidx " + self.FASTA
			print(cmd)
			os.system(cmd)

			FASTA_INDEX = self.FASTA.replace(".fa", ".fa.fai")
			CHROM_SIZES = self.outDir + self.outPrefix + "chrom.sizes"
			cmd = "cut -f 1,2 " + FASTA_INDEX + " > " + CHROM_SIZES
			print(cmd)
			os.system(cmd)

			SLOPPED_CAPPED_BED_FILE = self.outDir + self.outPrefix + os.path.basename(BAM_FILE).replace(".bam", ".slopped.capped.bed")
			cmd = "bedtools slop -i " + CAPPED_BED_FILE + " -b " + self.slopDistanceParameter + " -g " + CHROM_SIZES + " > " + SLOPPED_CAPPED_BED_FILE
			print(cmd)
			os.system(cmd)

			SORTED_SLOPPED_CAPPED_BED_FILE = self.outDir + self.outPrefix + os.path.basename(BAM_FILE).replace(".bam", ".sorted.slopped.capped.bed")
			cmd = "bedtools sort -i " + SLOPPED_CAPPED_BED_FILE + " > " + SORTED_SLOPPED_CAPPED_BED_FILE
			print(cmd)
			os.system(cmd)

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
		CAPPED_BED_FILE = self.outDir + self.outPrefix + "APrioriAnnotations.capped.bed"
		df_final.to_csv(CAPPED_BED_FILE, sep='\t', index=False, header=None)

		cmd = "samtools faidx " + self.FASTA
		print(cmd)
		os.system(cmd)

		FASTA_INDEX = self.FASTA.replace(".fa", ".fa.fai")
		CHROM_SIZES = self.outDir + self.outPrefix + "chrom.sizes"
		cmd = "cut -f 1,2 " + FASTA_INDEX + " > " + CHROM_SIZES
		print(cmd)
		os.system(cmd)

		df3 = pd.read_csv(CHROM_SIZES, sep = "\t", header = None, names = ["Chromosome", "Size"])
		ChromosomeNames = df3['Chromosome'].tolist()
		df_final = df_final[df_final.Chromosome.isin(ChromosomeNames)]
		df_final.to_csv(CAPPED_BED_FILE, sep='\t', index=False, header=None)

		SLOPPED_CAPPED_BED_FILE = self.outDir + os.path.basename(CAPPED_BED_FILE).replace(".capped.bed", ".slopped.capped.bed")
		cmd = "bedtools slop -i " + CAPPED_BED_FILE + " -b " + self.slopDistanceParameter + " -g " + CHROM_SIZES + " > " + SLOPPED_CAPPED_BED_FILE
		print(cmd)
		os.system(cmd)

		SORTED_SLOPPED_CAPPED_BED_FILE = self.outDir + os.path.basename(CAPPED_BED_FILE).replace(".capped.bed", ".sorted.slopped.capped.bed")
		cmd = "bedtools sort -i " + SLOPPED_CAPPED_BED_FILE + " > " + SORTED_SLOPPED_CAPPED_BED_FILE
		print(cmd)
		os.system(cmd)

	def _generateFinalSoftClippedBEDFile (self):
		SORTED_SLOPPED_CAPPED_BED_FILELIST = glob.glob(self.outDir + "*.sorted.slopped.capped.bed")
		SORTED_SLOPPED_CAPPED_BED_FILELIST_STRING = ' '.join(SORTED_SLOPPED_CAPPED_BED_FILELIST)

		CONCATENATED_SORTED_SLOPPED_CAPPED_BED_FILE = self.outDir + self.outPrefix + ".concatenated.sorted.slopped.capped.bed"
		cmd = "cat " + SORTED_SLOPPED_CAPPED_BED_FILELIST_STRING + " > " + CONCATENATED_SORTED_SLOPPED_CAPPED_BED_FILE
		print(cmd)
		os.system(cmd)

		self._clusterAndSortBED(inputBED = CONCATENATED_SORTED_SLOPPED_CAPPED_BED_FILE)
		
	def _clusterAndSortBED(self, inputBED):
		#Look at _generateFinalSoftClippedBEDFile and filterPolyADBandPolyASite functions; take the common code and refactor into this function
		
		CLUSTERED_BED_FILE = self.outDir + self.outPrefix + ".clustered.bed"
		cmd = "bedtools cluster -s -d 24 -i " + inputBED + " > " + CLUSTERED_BED_FILE
		print(cmd)
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

		cmd = "bedtools sort -i " + CLEANED_BED_FILE + " > " + self.FINAL_SOFTCLIPPED_BED_FILE
		print(cmd)
		os.system(cmd)

		self._addGeneNames()
		self._markFeatures()

	def _loadFilteredPolyADBandPolyASite(self):
		#Generate final Filtered BED File
		CONCATENATED_APRIORI_SOFTCLIPPED_BED_FILE = self.outDir + self.outPrefix + ".concatenated.apriori.softclipped.bed"
		cmd = "cat " + self.APRIORI + " " + self.FINAL_SOFTCLIPPED_BED_FILE + " > " + CONCATENATED_APRIORI_SOFTCLIPPED_BED_FILE
		print(cmd)
		os.system(cmd)

		self._clusterAndSortBED(inputBED = CONCATENATED_APRIORI_SOFTCLIPPED_BED_FILE)

	def _generateSAF(self):
		df = pd.read_csv(self.FINAL_SOFTCLIPPED_BED_FILE, sep='\t', index_col=None, names=['Chr', 'Start', 'End', 'Gene', 'APA', 'Strand'])

		# Uncomment bottom line once you add Gene and APA features!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		df['GeneID'] = df['Gene'] + '@' + df['Chr'] + '_' + df['Start'].apply(str) + '_' + df['End'].apply(str) + '_' + df['Strand']+ '@' + df['APA']

		# Change bottom command from Gene to GeneID!!!!!!!!!!!!!!!!!!!!!!!!!!!
		df = df[['GeneID', 'Chr', 'Start', 'End', 'Strand']]

		df.to_csv(self.FINAL_SAF_FILE, sep='\t', index=False, header=None)

	def _filterIgnoredFeatures(self):
		ignore_features="".join(self.ignoreFeatures).replace(" ","").split(",")
		
		if len(ignore_features) >=1:
			PA_FT_df=pd.read_csv(self.FINAL_SOFTCLIPPED_BED_FILE,sep="\t",header=None,index_col=None)
			PA_FT_df.columns=['chr','start','end','GeneID','Feature','strand']
			for ig in ignore_features:
				PA_FT_df=PA_FT_df[~PA_FT_df['Feature'].str.contains(ig)]
			PA_FT_df.to_csv(self.FINAL_SOFTCLIPPED_BED_FILE,sep="\t",header=False,index=False)

	def extractPolyA(self):
		self._truePolyA()
		self._generateFinalSoftClippedBEDFile()

		CPASBERT1 = CPASBERT(outDir = self.outDir,
		outPrefix = self.outPrefix, 
		modelPath = self.modelPath,
		fasta = self.FASTA,
		gtf = self.GTF,
		predictions = "",
		denovoapasiteBED = self.FINAL_SOFTCLIPPED_BED_FILE,
		geneBED = "",
		polyADB = self.POLYADB,
		polyASite = self.POLYASITE
		)
		
		CPASBERT1.filter_CPAS_Sites()

		if (self.apriori_annotations):
			self._loadFilteredPolyADBandPolyASite()

		self._filterIgnoredFeatures()
		self._generateSAF()


	def filterPolyADBandPolyASite(self):
		self._mergePolyADBandPolyASite()

		SORTED_SLOPPED_CAPPED_APRIORIANNOTATIONS_BED_FILE = self.outDir + self.outPrefix + "APrioriAnnotations.sorted.slopped.capped.bed"
		self._clusterAndSortBED(inputBED = SORTED_SLOPPED_CAPPED_APRIORIANNOTATIONS_BED_FILE)

		CPASBERT1 = CPASBERT(outDir = self.outDir,
		outPrefix = self.outPrefix, 
		modelPath = self.modelPath,
		fasta = self.FASTA,
		gtf = self.GTF,
		predictions = "",
		denovoapasiteBED = self.FINAL_SOFTCLIPPED_BED_FILE,
		geneBED = "",
		polyADB = self.POLYADB,
		polyASite = self.POLYASITE
		)

		CPASBERT1.filter_CPAS_Sites()

def main():
	ExtractPolyASites4Bulk1 = ExtractPolyAsites4Bulk(outDir = "/mnt/belinda_local/venkata/data/PolyAMiner-Bulk/TestFiles_Human_APriori",
		outPrefix = "human_", 
		fasta = "/mnt/belinda_local/venkata/data/Index_Files/Human/GenomeFasta_GTF/GRCh38.primary_assembly.genome.fa",
		gtf = "/mnt/belinda_local/venkata/data/Index_Files/Human/GenomeFasta_GTF/gencode.v33.primary_assembly.annotation.gtf",
		con1BAMFiles = "/mnt/belinda_local/venkata/data/Project_Meningioma_AkashPatel_NSG/hari_APA_Akash/02_BAM/TL-21-VZKP229D/TL-21-VZKP229D.sorted.bam,/mnt/belinda_local/venkata/data/Project_Meningioma_AkashPatel_NSG/hari_APA_Akash/02_BAM/TL-21-QGUU886F/TL-21-QGUU886F.sorted.bam,/mnt/belinda_local/venkata/data/Project_Meningioma_AkashPatel_NSG/hari_APA_Akash/02_BAM/TL-20-DF4101/TL-20-DF4101.sorted.bam,/mnt/belinda_local/venkata/data/Project_Meningioma_AkashPatel_NSG/hari_APA_Akash/02_BAM/TL-20-879184/TL-20-879184.sorted.bam,/mnt/belinda_local/venkata/data/Project_Meningioma_AkashPatel_NSG/hari_APA_Akash/02_BAM/TL-20-664423/TL-20-664423.sorted.bam,/mnt/belinda_local/venkata/data/Project_Meningioma_AkashPatel_NSG/hari_APA_Akash/02_BAM/TL-20-36A961/TL-20-36A961.sorted.bam,/mnt/belinda_local/venkata/data/Project_Meningioma_AkashPatel_NSG/hari_APA_Akash/02_BAM/TL-20-9A6CA1/TL-20-9A6CA1.sorted.bam,/mnt/belinda_local/venkata/data/Project_Meningioma_AkashPatel_NSG/hari_APA_Akash/02_BAM/TL-20-7E5975/TL-20-7E5975.sorted.bam", 
		con2BAMFiles = "/mnt/belinda_local/venkata/data/Project_Meningioma_AkashPatel_NSG/hari_APA_Akash/02_BAM/TL-21-ZVRBQA9H/TL-21-ZVRBQA9H.sorted.bam,/mnt/belinda_local/venkata/data/Project_Meningioma_AkashPatel_NSG/hari_APA_Akash/02_BAM/TL-21-E6727A/TL-21-E6727A.sorted.bam,/mnt/belinda_local/venkata/data/Project_Meningioma_AkashPatel_NSG/hari_APA_Akash/02_BAM/TL-21-62427B/TL-21-62427B.sorted.bam,/mnt/belinda_local/venkata/data/Project_Meningioma_AkashPatel_NSG/hari_APA_Akash/02_BAM/TL-21-84U67PPG/TL-21-84U67PPG.sorted.bam,/mnt/belinda_local/venkata/data/Project_Meningioma_AkashPatel_NSG/hari_APA_Akash/02_BAM/TL-20-56286D/TL-20-56286D.sorted.bam,/mnt/belinda_local/venkata/data/Project_Meningioma_AkashPatel_NSG/hari_APA_Akash/02_BAM/TL-20-2F7E80/TL-20-2F7E80.sorted.bam,/mnt/belinda_local/venkata/data/Project_Meningioma_AkashPatel_NSG/hari_APA_Akash/02_BAM/TL-19-EBC5FF/TL-19-EBC5FF.sorted.bam,/mnt/belinda_local/venkata/data/Project_Meningioma_AkashPatel_NSG/hari_APA_Akash/02_BAM/TL-19-C46B1C/TL-19-C46B1C.sorted.bam",
		proportionA = "0.90,0.85,0.80,0.75",
		modelOrganism = "human",
		apriori_annotations = True,
		ignoreFeatures = "" 
		)

	ExtractPolyASites4Bulk1.filterPolyADBandPolyASite()

if __name__ == "__main__":
    main()

