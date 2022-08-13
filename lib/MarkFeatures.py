#v1.2
import os, sys, glob
import pandas as pd, pybedtools as pb, gtfparse as gp



def makeGeneBed(args_gtf,args_fkey,args_o):
	args_o=args_o.strip("/")+"/"+args_fkey
	df=gp.read_gtf(args_gtf)
	df=df[df['seqname'].astype(str).str.contains('chr')]
	# Gene table in bed format #
	df_gene=df[df['feature']=='gene']
	df_gene=df_gene[['seqname','start', 'end','gene_id','gene_name','strand']]
	df_gene=df_gene.sort_values(["seqname", "start","end"], ascending = (True, True,True))
	df_gene.to_csv(args_o+'.Genes.bed', sep='\t', index=False,header=False)
	return(args_o+'.Genes.bed')



def getTables(df):
	# Gene table in bed format #
	#df_gene=df[df['feature']=='gene']
	#df_gene=df_gene[['seqname','start', 'end','gene_id','gene_name','strand']]
	#df_gene=df_gene.sort_values(["seqname", "start","end"], ascending = (True, True,True))

	# Exon table #
	#df_exn=df[df['feature']=='exon']
	#df_exn=df_exn[['seqname','start', 'end','gene_id','exon_id','strand']]
	#df_exn=df_exn[df_exn['seqname'].str.contains("chr")]
	#df_exn=df_exn.sort_values(["seqname", "start","end"], ascending = (True, True,True))

	# CDS/Coding Exon table #
	df_cds=df[df['feature']=='CDS']
	df_cds=df_cds[['seqname','start', 'end','gene_id','exon_id','strand']]
	df_cds=df_cds[df_cds['seqname'].str.contains("chr")]
	df_cds['exon_id']="CDS"
	df_cds=df_cds.sort_values(["seqname", "start","end"], ascending = (True, True,True))
	
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
	df_utr['status']=checkUTR(df_utr.values)
	df_utr5=df_utr[df_utr['status']=="5UTR"]
	df_utr5=df_utr5[['seqname','start','end','gene_id','transcript_id','strand']]
	df_utr5['transcript_id']='UTR5'
	df_utr3=df_utr[df_utr['status']=="3UTR"]
	df_utr3=df_utr3[['seqname','start','end','gene_id','transcript_id','strand']]
	df_utr3['transcript_id']='UTR3'

	return(df_cds,df_utr5,df_utr3)

def checkUTR(m):
	return_list=[]
	for i in range(0,len(m)):
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


def mapAPA2Features(APAfile,Features):
	PA_df=pd.read_csv(APAfile,sep="\t",header=0,index_col=None)
	PA_df[['Chr', 'Start', 'End',"Strand"]]=PA_df['feature_id'].str.split("_", expand=True)
	PA_df=PA_df[['Chr', 'Start', 'End', 'gene_id', 'feature_id', 'Strand']]
	PA=pb.BedTool.from_dataframe(PA_df)

	FT = pb.BedTool(Features)
	PA_FT_bed = PA.intersect(FT, nonamecheck=True, s=True, wao=True)
	PA_FT_df=PA_FT_bed.to_dataframe(disable_auto_names=True, header=None)
	PA_FT_df.columns=['Chr', 'Start', 'End', 'GeneID', 'APAID', 'Strand','Chr_2', 'Start_2', 'End_2', 'GeneID_2', 'APAID_2', 'Strand_2','overlap']	
	
	PA_FT_df_nomap=PA_FT_df[PA_FT_df['Chr_2']=="."].copy()
	PA_FT_df_nomap=PA_FT_df_nomap[['Chr', 'Start', 'End', 'GeneID', 'APAID', 'Strand']]
	PA_FT_df_nomap['APAID']="NA"
		
	PA_nomap_FT=pb.BedTool.from_dataframe(PA_FT_df_nomap).intersect(FT, nonamecheck=True, s=False, wao=True)
	PA_nomap_FT_df=PA_nomap_FT.to_dataframe(disable_auto_names=True, header=None)
	PA_nomap_FT_df.columns=['Chr', 'Start', 'End', 'GeneID', 'APAID', 'Strand','Chr_2', 'Start_2', 'End_2', 'GeneID_2', 'APAID_2', 'Strand_2','overlap']
	UTR3_EX_df=PA_nomap_FT_df[PA_nomap_FT_df['Start_2']==-1].copy()
	UTR3_EX_df=UTR3_EX_df[['Chr', 'Start', 'End', 'GeneID', 'APAID', 'Strand']]
	UTR3_EX_df['APAID']="UTR3_EX"

	PA_FT_df=PA_FT_df[(((PA_FT_df['Strand'] =="+") & (PA_FT_df['APAID_2'] !="UTR3")) & (PA_FT_df['End'] >= PA_FT_df['Start_2']) & (PA_FT_df['End'] <= PA_FT_df['End_2'])) | (((PA_FT_df['Strand'] =="+") & (PA_FT_df['APAID_2'] =="UTR3")) & (PA_FT_df['End'] >= PA_FT_df['Start_2'])) | (((PA_FT_df['Strand'] =="-") & (PA_FT_df['APAID_2'] !="UTR3")) & (PA_FT_df['Start'] <= PA_FT_df['End_2']) & (PA_FT_df['Start'] >= PA_FT_df['Start_2'])) | (((PA_FT_df['Strand'] =="-") & (PA_FT_df['APAID_2'] =="UTR3")) & (PA_FT_df['Start'] <= PA_FT_df['End_2']))]
	PA_FT_df['APAID']=PA_FT_df['APAID_2']
	PA_FT_df=PA_FT_df[['Chr', 'Start', 'End', 'GeneID', 'APAID', 'Strand']]
	PA_FT_bed=pb.BedTool.from_dataframe(PA_FT_df)
	PA_FT_bed = PA_FT_bed.sort()
	PA_FT_bed=PA_FT_bed.merge(c=[4,5,6],s=True,o="distinct")
	PA_FT_df=PA_FT_bed.to_dataframe(disable_auto_names=True, header=None)
	PA_FT_df.columns=['Chr', 'Start', 'End', 'GeneID', 'APAID', 'Strand']
	
	PAmap_df=pd.concat([PA_FT_df,UTR3_EX_df], axis=0,ignore_index=True)
	PAmap_df['feature_id']=PAmap_df['Chr']+"_"+PAmap_df['Start'].map(str)+"_"+PAmap_df['End'].map(str)+"_"+PAmap_df['Strand']
	PAmap_df=PAmap_df[['feature_id','APAID']]
	
	PA_df=pd.read_csv(APAfile,sep="\t",header=0,index_col=None)
	PA_df=PA_df.merge(PAmap_df,on=['feature_id'],how='left').fillna("NA")
	PA_df['feature_id']=PA_df['feature_id']+"_"+PA_df['APAID']
	PA_df=PA_df.drop(columns=['APAID'])
	PA_df.to_csv(APAfile, sep='\t', index=False, header=True)
	

def MarkFeatures(args_gtf,args_fkey,args_o,APAfile):
	args_o=args_o.strip("/")+"/"+args_fkey
	genes_bed=pb.BedTool(args_o+'.Genes.bed')

	df=gp.read_gtf(args_gtf)
	df=df[df['seqname'].astype(str).str.contains('chr')]
	cds_df,utr5_df,utr3_df=getTables(df)

	#gene_df.to_csv(args_o+'.Genes.bed',sep="\t",header=False,index=False)
	#genes_bed=pb.BedTool(args_o+'.Genes.bed')

	utr5_bed=pb.BedTool.from_dataframe(utr5_df)
	utr5_bed = utr5_bed.sort()
	utr5_merge=utr5_bed.merge(s=True, c=[4,5,6], o='distinct', d=1)
	utr5_df=utr5_merge.to_dataframe(disable_auto_names=True,header=None)
	utr5_df.columns=['chr','start','end','GeneID','Feature','strand']

	utr3_bed=pb.BedTool.from_dataframe(utr3_df)
	utr3_bed = utr3_bed.sort()
	utr3_merge=utr3_bed.merge(s=True, c=[4,5,6], o='distinct', d=1)
	utr3_df=utr3_merge.to_dataframe(disable_auto_names=True,header=None)
	utr3_df.columns=['chr','start','end','GeneID','Feature','strand']

	cds_bed=pb.BedTool.from_dataframe(cds_df)
	cds_bed = cds_bed.sort()
	cds_merge=cds_bed.merge(s=True, c=[4,5,6], o='distinct', d=1)
	cds_merge=cds_merge.subtract(utr5_merge, s=True)
	cds_merge=cds_merge.subtract(utr3_merge, s=True)
	cds_df=cds_merge.to_dataframe(disable_auto_names=True,header=None)
	cds_df.columns=['chr','start','end','GeneID','Feature','strand']

	'''
	exn_bed=pb.BedTool.from_dataframe(exon_df)
	exn_bed = exn_bed.sort()
	exn_merge=exn_bed.merge(s=True, c=[4,5,6], o='distinct', d=1)
	exn_merge=exn_merge.subtract(utr5_merge, s=True)
	exn_merge=exn_merge.subtract(utr3_merge, s=True)
	exn_df=exn_merge.to_dataframe(disable_auto_names=True,header=None)
	exn_df.columns=['chr','start','end','GeneID','Feature','strand']
	exn_df.to_csv("EXN.bed",sep="\t",header=False,index=False)
	'''

	intr_bed = genes_bed.subtract(cds_merge, s=True)
	intr_bed = intr_bed.subtract(utr5_merge, s=True)
	intr_bed = intr_bed.subtract(utr3_merge, s=True)
	intr_df=intr_bed.to_dataframe(disable_auto_names=True,header=None)
	intr_df.columns=['chr','start','end','GeneID','Feature','strand']
	intr_df['Feature']='Intron'
	
	feature_df=pd.concat([intr_df,cds_df,utr5_df,utr3_df], axis=0,ignore_index=True)
	feature_df.to_csv(args_o+".Features.bed",sep="\t",header=False,index=False)


	mapAPA2Features(APAfile,args_o+".Features.bed")
	return(1)
