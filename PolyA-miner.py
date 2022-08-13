# PolyA-miner.py -help # 
# Venkata Jonnakuti et al., last update 01/21/2022 #
# V1.5
# Tested on Python 3.9.5; look out for rstrip if using newer versions (alternative is v3.9.8 without rstrip)

import os, sys
import time,argparse,subprocess, pandas as pd
sys.path.append("/".join(os.path.abspath(sys.argv[0]).split("/")[0:-1])+"/lib")
import CheckDependency, DataProcessing
import MakeAPAMatrix, GenePolyAIndex, STest, PAusage
from DEGAnalyzer import DEGAnalyzer
from ExtractPolyAsites4Bulk import ExtractPolyAsites4Bulk
from PolyASafety import PolyASafety


def check_files (checkfiles,logfile):
	for checkfile in checkfiles:
		if os.path.exists(checkfile):
			pass
		else:
			with open(logfile, "a") as fileObj:
				fileObj.write("\nError cannot locate "+checkfile+"\n")
			print("\nError cannot locate "+checkfile+"\n")
			exit()
	return(1)

def logEvent(logfile, event):
	localdate = time.strftime('%a %m/%d/%Y')
	localtime = time.strftime('%H:%M:%S')
	with open(logfile, "a") as fileObj:
		fileObj.write('# '+ event + ' on: ' + localdate + ' at: ' + localtime + ' \n')
		print('# '+ event + ' on: ' + localdate + ' at: ' + localtime)

def main():
	parser = argparse.ArgumentParser(description='''PolyA-miner v1.5: Inferring alternative poly-adenylation changes
	from bulk-RNAseq data  - Venkata Jonnakuti et al., \n''',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	optional = parser._action_groups.pop()
	required = parser.add_argument_group('Required arguments')
	parser._action_groups.append(optional)

	required.add_argument('-mode',help='Run mode options: \'bam\' to start from mapped data, \'fastq\' to start from raw data',choices=['bam','fastq'],default='bam',required='True',type=str)
	optional.add_argument('-d',help='Base directory of input fastq files. Valid for -mode fastq ',type=str)
	optional.add_argument('-o',help='Out put directory',type=str,default='PolyAminer_OUT')
	required.add_argument('-c1',help='Comma-separated list of condition1 files. Full path for BAMs (index files are also expected) or Just file names for fastq', nargs='+',required='True',type=str)
	required.add_argument('-c2',help='Comma-separated list of condition2 files. Full path for BAMs (index files are also expected) or Just file names for fastq ', nargs='+',required='True',type=str)
	parser.add_argument('-s',help='Strand information 0: un-stranded 1: fwd-strand 2:rev-strand. ',choices=[0,1,2],type=int,default=0)
	
	# Ref. files
	optional.add_argument('-index',help='Reference genome bowtie2 index. Valid for -mode fastq',type=str)
	required.add_argument('-fasta',help='Reference fasta sequence',required='True',type=str)
	required.add_argument('-gtf',help='Reference gtf file',required='True',type=str)
	required.add_argument('-pa',help='PolyA annotations file standard 6 column bed format',type=str)
	required.add_argument('-apriori_annotations',help='Use pre-loaded a priori PolyASite 2.0 and PolyADB 3.0 annotations', action=argparse.BooleanOptionalAction)
	
	# Optional #
	optional.add_argument('-umi',help='Length of UMIs, 0 if not used', type=int,default=0)
	optional.add_argument('-ignore',help=' Comma-separated list of regions to igonore from APA analysis UTR5, Introns, CDS, UTR3', nargs='+',type=str,default='ZZZZZZ')
	optional.add_argument('-apaBlock',help='Window size for annotated polyA sites',type=int, default=30)
	optional.add_argument('-mdapa',help='Cluster distance for annotated polyA sites: Merge polyA sites with in this distance. ',type=int, default=0)
	optional.add_argument('-md',help='Cluster distance for de-novo polyA sites: Merge polyA sites with in this distance',type=int, default=0)
	optional.add_argument('-anchor',help='Overlap in "bp" for mapping de-novo polyA sites to annotated polyA sites ',type=int, default=1)
	
	# Additional analysis 
	optional.add_argument('-cluster_onGenes',help='Cluster samples based on gene counts - all polyA sites / 3utr polyA sites',choices=['all','3utr','none'],type=str,default='none')
	optional.add_argument('-cluster_onPAsites',help='Cluster samples based on gene counts - all polyA sites / 3utr polyA sites',choices=['all','3utr','none'],type=str,default='none')
	optional.add_argument('-DEG',help='Perform differential gene expression analysis - all polyA sites / 3utr polyA sites',choices=['all','3utr','none'],type=str,default='none')
	optional.add_argument('-pa_usage',help='Perform differential polyA site usage analysis - all polyA sites / filtered polyA sites',choices=['all','filtered','none'],type=str,default='none')

	# Tuning #
	optional.add_argument('-expNovel',help='Explore novel APA sites 0: only annotated sites 1: de-novo',choices=[1,0],type=int,default=0)
	optional.add_argument('-novel_d',help='Distance from annotated TES to map novel pA sites',type=int, default=1000)
	optional.add_argument('-p',help='No. of processors to use',type=int,default=4)
	optional.add_argument('-ip_d',help='Downstream internal priming window',type=int, default=50)
	optional.add_argument('-ip_u',help='Upstream internal priming window',type=int, default=50)
	optional.add_argument('-a',help='Internal priming polyA fraction',type=float, default=0.65)
	optional.add_argument('-pa_p',help='pOverA filter: P ',type=float, default=0.6)
	optional.add_argument('-pa_a',help='pOverA filter: A ',type=int, default=5)
	optional.add_argument('-pa_m',help='pOverA filter: M ',type=int, default=2)
	optional.add_argument('-gene_min',help='Min counts per Gene',type=int, default=10)
	optional.add_argument('-apa_min',help='Min. proportion per APA',type=float, default=0.05)
	optional.add_argument('-batchCorrection',help='Comma-seperated list of batch membership for ComBat-seq correction',type=str, default='ZZZZZZ')
	optional.add_argument('-modelOrganism',help='Model organism for CPAS-BERT Model',type=str, default="human")

	#
	optional.add_argument('-t',help='Statistical Test- BB: for beta-binomial or iNMF: for iterative NMF. For small sample size BB is recommended ',choices=['BB','iNMF'],type=str,default="BB")
	optional.add_argument('-i',help='No. of NMF iterations. Valid only for -t iNMF',type=int,default=100)
	optional.add_argument('-outPrefix',help='Output file/s prefix', default="PolyAminer_Out",type=str)

	if len(sys.argv)==1:
		parser.print_help(sys.stderr)
		sys.exit(1)
	args=parser.parse_args()
	args_dict = vars(args)
	args.o=args.o.rstrip("/")+"/"

	############################################################################################################
	# Module 0: Safety Checks -- Check output directory, start logging, check dependencies, check input files  #
	############################################################################################################

	if args.mode =='bam' and args.umi >1:
		print("\nUMI option is only valid in fastq mode ...\n")
		exit()
	try:
		if args.mode =='bam' and len(args.index) >= 1:
			print("\nindex option is only valid in fastq mode. Skipping ...\n")
	except:
		pass

	# checkfiles = args.c1 + args.c2 + [args.fasta] + [args.gtf] + [args.pa]
	# checkfiles = [args.fasta] + [args.gtf] + [args.pa]
	checkfiles = args.c1 + args.c2 + [args.fasta] + [args.gtf]
	logfile = args.o + args.outPrefix + '.PolyA-miner.log.txt'

	PolyASafety1 = PolyASafety(argsDict = args_dict, outDir = args.o, inputFiles = checkfiles, logfile = logfile)
	if (PolyASafety1.runInitialSafetyChecks() == 0):
		exit()

	logBook = open(logfile,'w')

	# if no ref. polyA use dummy #
	localdate = time.strftime('%a_%m_%d_%Y')
	localtime = time.strftime('%H_%M_%S')
	if args.pa is None:
		if args.expNovel == 1:
			dummy_refPA=args.o+"DummyRefPA_"+localdate+"_"+localtime+".bed"
			fw=open(dummy_refPA,"w")
			fw.write("chr1	634823	634823	D000001	DummyGene	+\n")
			fw.close()
			args.pa=dummy_refPA
		else:
			logBook.write("Err: Reference PolyA file is needed if expNovel is set to 0 ...\n")
			print("Err: Reference PolyA file is needed if expNovel is set to 0 ...\n")
			exit()

	###################################################
	# Module 1A: Arguments check and Data preprocessing 
	###################################################
	# Number of samples #
	controls="".join(args.c1).replace(" ","").split(",")
	nc=len(controls)
	treated="".join(args.c2).replace(" ","").split(",")
	nt=len(treated)
	args.bed = args.o + args.outPrefix + "Genes.bed"
	
	if args.mode=='bam':
		logEvent(logfile = logfile, event = 'Run mode: from preprocessed bam files')
		# check files #
		checkfiles=controls+treated+[args.fasta]+[args.pa]
		if check_files(checkfiles,logfile):
			logEvent(logfile = logfile, event = 'Arguments checked')
		
	if args.mode=='fastq':
		logEvent(logfile = logfile, event = 'Run mode: from raw fastq data')
		# check basedir #:
		if os.path.isdir(args.d):
			pass
		else:
			logBook.write("\nError cannot locate "+args.d+"\n")
			exit()
		
		# Check Ref #:
		rg="/".join(args.index.split('/')[:-1])
		if os.path.isdir(rg):
			pass
		else:
			logBook.write("\nError cannot locate "+rg+"\n")
			exit()
		
		# Check files #
		cck=[]
		for c in controls:
			cck.append(args.d+"/"+c)
		tck=[]
		for t in treated:
			tck.append(args.d+"/"+t)
		checkfiles=cck+tck+[args.fasta]+[args.pa]
		if check_files(checkfiles,logfile):
			logEvent(logfile = logfile, event = 'Arguments checked')
		try:
			samples = controls + treated
			nc = len(controls)
			nt = len(treated)
			baseDir =args.d; outDir=args.o; np=args.p; ref_genome=args.index; fkey=args.outPrefix;
			for s in samples:
				DataProcessing.process_rawfastq(args.d,args.o, s, args.umi,args.p)
				logEvent(logfile = logfile, event = 'Finished fastp/umi')
				DataProcessing.mapping_bowtie2(args.o, s, args.p, args.index,args.umi)
				logEvent(logfile = logfile, event = 'Finished mapping')
			logEvent(logfile = logfile, event = 'Completed data processing')
			pass

		except:
			logEvent(logfile = logfile, event = 'Error in data processing module')
			logfile.close()
			exit()

	args.o=args.o.rstrip("/")+"/"

	################################
	# Module2: Extract PolyA sites #
	################################

	ExtractPolyASites4Bulk1 = ExtractPolyAsites4Bulk(outDir = args.o,
		outPrefix = args.outPrefix, 
		fasta = args.fasta,
		gtf = args.gtf,
		con1BAMFiles = args.c1, 
		con2BAMFiles = args.c2,
		proportionA = "0.90,0.85,0.80,0.75",
		modelOrganism = args.modelOrganism,
		apriori_annotations = args.apriori_annotations
		)

	if args.expNovel == 1:
		if ExtractPolyASites4Bulk1.extractPolyA():
			logEvent(logfile = logfile, event = 'Completed extracting soft-clipped (Poly(A)-Capped) polyadenylation sites')
			pass
		else:
			logEvent(logfile = logfile, event = "Error in extracting soft-clipped (Poly(A)-Capped) polyadenylation sites")
			exit()
	

	###################################
	# Module 3: Make APA count matrix #
	###################################

	if MakeAPAMatrix.MakeMatrix(args.o, args.p, args.outPrefix, args.pa_p, args.pa_a, args.pa_m, controls, treated, args.apa_min, args.gene_min, args.mode,controls+treated,logBook,args.pa_usage) == 1:
		logEvent(logfile = logfile, event = 'Completed abstracting APA proportions')
		pass
	else:
		logEvent(logfile = logfile, event = "Error in abstracting APA proportions")
		exit()

	#######################################################################
	# Module 3.1: RUN BATCH CORRECTION USING COMBAT-SEQ IF USER SPECIFIED #
	#######################################################################

	if args.batchCorrection != "ZZZZZZ":
		try:
			CombatSeqCorrectionRScriptLoc = os.path.dirname(__file__) + "/lib/ComBatSeqBatchCorrection.r"
			APAMatrixLoc = args.o.rstrip("/")+"/"+args.outPrefix + '_APA.CountMatrix.GFil.PA.PR.txt'
			DEGMatrixLoc = args.o.rstrip("/")+"/"+args.outPrefix + '_APACountMatrix4DEGs.txt'
			
			cmd = "Rscript " + CombatSeqCorrectionRScriptLoc + " " + APAMatrixLoc + " " + args.batchCorrection
			print(cmd)
			os.system(cmd)

			cmd = "Rscript " + CombatSeqCorrectionRScriptLoc + " " + DEGMatrixLoc + " " + args.batchCorrection
			print(cmd)
			os.system(cmd)

			logEvent(logfile = logfile, event = 'Completed performing batch corrections')
			pass
		except:
			logEvent(logfile = logfile, event = "Error in performing batch correction!")
			exit()

	################################
	# Module 3.5: RUN DEG ANALYSIS #
	################################

	DEGAnalyzer1 = DEGAnalyzer(outDir = args.o, outPrefix = args.outPrefix, bed = args.bed, gtf = args.gtf, condition1Samples = args.c1, condition2Samples = args.c2)
	DEGAnalyzer1.clusterAndAnalyzeDEGs()
	logEvent(logfile = logfile, event = 'Completed DEG (Overall and Core APA Factor) Analysis')
	
	###################################
	# Module 4: Gene level PolyA Index 
	###################################
	if GenePolyAIndex.VectorPro(args.o.rstrip("/")+"/", args.outPrefix, nc, nt, args.p, args.o.rstrip("/")+"/"+args.outPrefix+ '_APA.CountMatrix.GFil.PA.PR.txt', args.o.rstrip("/")+"/LibSize.txt"):
		logEvent(logfile = logfile, event = 'Completed computing PolyA Index')
		pass
	else:
		logEvent(logfile = logfile, event = "Error in computing PolyA Index")
		exit()

	
	###################################
	# Module 5: Stat BetaBinomial iNMF 
	###################################
	if args.t == "BB":
		logEvent(logfile = logfile, event = 'Using BB model')
		if STest.runBBtest(args.o.rstrip("/")+"/"+args.outPrefix+ '_APA.CountMatrix.GFil.PA.PR.txt', nc, nt, args.o, args.outPrefix, args.p,logBook):
			logEvent(logfile = logfile, event = 'Completed beta-binomial testing')
		else:
			logEvent(logfile = logfile, event = "Error in beta-binomial testing")
			exit()

	if args.t =='iNMF':
		logEvent(logfile = logfile, event = 'Using iNMF model')
		if STest.iNMFtest(args.o,args.outPrefix, nc, nt, args.i,args.p, args.o.rstrip("/")+"/"+args.outPrefix+ '_APA.CountMatrix.GFil.PA.PR.txt', logBook):
			logEvent(logfile = logfile, event = 'Completed iNMF testing')
		else:
			logEvent(logfile = logfile, event = "Error in iNMF testing")
			exit()

	
	# Merge 4 and 5 add gene symbols #
	pdata=pd.read_csv(args.o.rstrip("/")+"/"+args.outPrefix+'_PolyA-miner.Results.txt',sep="\t",header=0,index_col=None)
	stat_data=pd.read_csv(args.o.rstrip("/")+"/"+args.outPrefix+"_Gene_Stats.txt",sep="\t",header=0,index_col=None)
	pdata=pd.merge(pdata,stat_data,left_on=["Gene"],right_on=["Gene"],how="outer")
	
	genes=pd.read_csv(args.bed,sep="\t",header=None,index_col=None)
	genes.columns=["Chr","Start","End","Gene","Symbol","Strand"]
	if "." in pdata['Gene'].iloc[0]:
		pdata=pd.merge(pdata,genes,on=["Gene"],how="left")
		pdata=pdata.drop(columns=["Chr","Start","End","Strand"])
	else:
		try:
			genes[['Gene','Version']] = genes['Gene'].str.split('.',expand=True)
			pdata=pd.merge(pdata,genes,on=["Gene"],how="left")
			pdata=pdata.drop(columns=["Chr","Start","End","Strand","Version"])
		except:
			pdata=pd.merge(pdata,genes,on=["Gene"],how="left")
			pdata=pdata.drop(columns=["Chr","Start","End","Strand"])
	pdata.to_csv(args.o.rstrip("/")+"/"+args.outPrefix+'_PolyA-miner.Results.txt',sep="\t",header=True,index=False)
	
	# Summary #
	logEvent(logfile = logfile, event = 'PolyA-miner results summary')
	nsig_c=pdata[pdata['AdjG_Pval']<=0.05].shape[0]
	nsig_s=pdata[(pdata['AdjG_Pval']<=0.05) & (pdata['PolyAIndex'] <0)].shape[0]
	nsig_l=pdata[(pdata['AdjG_Pval']<=0.05) & (pdata['PolyAIndex'] >0)].shape[0]
	logBook.write('# Significant PolyA changes:\t'+str(nsig_c)+"\n"+'# Significant shortening changes:\t'+str(nsig_s)+"\n"+'# Significant lengthening changes:\t'+str(nsig_l)+"\n")
	logEvent(logfile = logfile, event = 'Finished PolyA-miner APA analysis')

	#Create APA Volcano Plot
	apafile = args.o.rstrip("/")+"/"+args.outPrefix+'_PolyA-miner.Results.txt'
	data=pd.read_csv(apafile,sep="\t",header=0,index_col=None)

	temp1 = data[data['AdjG_Pval']<=0.05]
	temp2 = data[(data['AdjG_Pval']<=0.05) & ((data['PolyAIndex']>=0.5) | (data['PolyAIndex']<=-0.5))]
	sorteddf=temp2.sort_values(by=['PolyAIndex'])
	newlist = [x for x in list(sorteddf.head(8)['Gene']) if pd.isnull(x) == False and x != 'nan']
	x="_".join(newlist)
	print(x)
	#print(sorteddf.head(5))
	#print(sorteddf.tail(5))
	newlist = [x for x in list(sorteddf.tail(8)['Gene']) if pd.isnull(x) == False and x != 'nan']
	y="_".join(newlist)
	print(y)
	temp3 = data[(data['AdjG_Pval']<=0.05) & (data['PolyAIndex']>=0.5)]
	temp4 = data[(data['AdjG_Pval']<=0.05) & (data['PolyAIndex']<=-0.5)]
	#temp3 = data[(data['adj_pvalue']<=0.05) & (data['PolyAIndex']>=0)]
	#temp4 = data[(data['adj_pvalue']<=0.05) & (data['PolyAIndex']<=0)]

	ndeg=str(temp1.shape[0])
	fcdeg=str(temp2.shape[0])
	updeg=str(temp3.shape[0])
	dndeg=str(temp4.shape[0])

	DrawDAGVolcanoPlotRScriptLoc = os.path.dirname(__file__) + "/lib/APAvolcanoPlot.r"
	DAGvolcanoPlotFileName = args.o.rstrip("/")+"/"+args.outPrefix+'_DAGVolcanoPlot.tiff'

	cmd = "Rscript " + DrawDAGVolcanoPlotRScriptLoc + " " + apafile+" "+DAGvolcanoPlotFileName+" "+ndeg+" "+fcdeg+" "+updeg+" "+dndeg+" "+x+" "+y+" "+"DAG_Volcano_Plot"
	print(cmd)
	os.system(cmd)		
	

	#######################################
	# Module 6: PolyA Site Usage + Tidy Up
	#######################################
	if args.pa_usage !="none":
		if(PAusage.PAusage(args.o.rstrip("/")+"/"+args.outPrefix+"Genes.bed",args.o.rstrip("/")+"/"+args.outPrefix+"_PAusage.CountMatrix.txt",args.o,args.outPrefix,args.p,logBook)):
			logEvent(logfile = logfile, event = 'Finished PolyA site usage analysis')
		
		else:
			logEvent(logfile = logfile, event = 'Failed PolyA site usage analysis')
	
	DAGresultFile = args.o.rstrip("/")+"/"+args.outPrefix+'_PolyA-miner.Results.txt'
	DEGresultFile = args.o.rstrip("/")+"/"+args.outPrefix+"_DEG-results.txt"
	DEGHeatmap = args.o.rstrip("/")+"/"+args.outPrefix+"_overallDEGHeatmap.pdf"
	DEGVolcanoPlot = args.o.rstrip("/")+"/"+args.outPrefix+"_overallDEGVolcanoPlot.tiff"
	APAFactorHeatmap = args.o.rstrip("/")+"/"+args.outPrefix+"_APAFactor_DEGHeatmap.pdf"
	APAFactorVolcanoPlot = args.o.rstrip("/")+"/"+args.outPrefix+"_APAFactor_DEGVolcanoPlot.tiff"
	DAGVolcanoPlot = args.o.rstrip("/")+"/"+args.outPrefix+'_DAGVolcanoPlot.tiff'
	PA_PCAPlot = args.o.rstrip("/")+"/"+args.outPrefix+"_PA.PCA.tif"
	PA_tSNEPlot = args.o.rstrip("/")+"/"+args.outPrefix+"_PA.t-SNE.tif"
	Gene_PCAPlot = args.o.rstrip("/")+"/"+args.outPrefix+"_Gene.PCA.tif"
	Gene_tSNEPlot = args.o.rstrip("/")+"/"+args.outPrefix+"_Gene.t-SNE.tif"
	CDS_PAusage = args.o.rstrip("/")+"/"+args.outPrefix+"_CDS_PAusage_Results.txt"
	Intron_PAusage = args.o.rstrip("/")+"/"+args.outPrefix+"_Intron_PAusage_Results.txt"
	UTR3_PAusage = args.o.rstrip("/")+"/"+args.outPrefix+"_UTR3_PAusage_Results.txt"
	UTR5_PAusage = args.o.rstrip("/")+"/"+args.outPrefix+"_UTR5_PAusage_Results.txt"

	keepFiles = [logfile, DAGresultFile, DEGresultFile, DEGHeatmap, DEGVolcanoPlot, APAFactorHeatmap, APAFactorVolcanoPlot, DAGVolcanoPlot, PA_PCAPlot, PA_tSNEPlot, Gene_PCAPlot, Gene_tSNEPlot, CDS_PAusage, Intron_PAusage, UTR3_PAusage, UTR5_PAusage]

	# PolyASafety1.tidyUp(keep = keepFiles)

	logBook.close()

if __name__ == "__main__":
	main()
