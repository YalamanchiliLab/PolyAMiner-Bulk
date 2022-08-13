#  PolyA gene feature analysis  # 
#  Reads filtered Matrix from PolyA-miner #

import os, sys, glob, time, math, warnings
import numpy as np, pandas as pd, statistics
import statsmodels.stats.multitest as multi
import concurrent.futures as cf
from numpy import array, dot, arccos, clip
from scipy import stats

def runBBtest(apa_df, feature, nc, nt, out, key, npc,logfile):
	from rpy2.robjects.packages import importr
	import rpy2.robjects as ro
	from rpy2.robjects import pandas2ri
	pandas2ri.activate()
	ro.r['options'](warn=-1)
	countdata=importr('countdata')
	cols=apa_df.columns
	# Adjust counts #
	apa_mat=apa_df.values
	for row in apa_mat:
		ck=row[1:1+nc]; tk=row[1+nc:1+nc+nt]
		cn=row[1+nc+nt:1+nc+nt+nc]; tn=row[1+nc+nt+nc:1+nc+nt+nc+nt]
		if (ck==cn).all():
			cn=cn+1
			row[1+nc+nt:1+nc+nt+nc]=cn
		if (tk==tn).all():
			tn=tn+1
			row[1+nc+nt+nc:1+nc+nt+nc+nt]=tn
	apa_df2 = pd.DataFrame(apa_mat, columns = cols)
	
	ns=nc+nt
	K=cols[1:(ns+1)];
	N=cols[(ns+1):((ns*2)+1)]
	Kdf=apa_df2[K] # pad 1 ?
	Ndf=apa_df2[N] # pad 1 ?

	Kdf_matrix=Kdf.values
	Kdf = pd.DataFrame(Kdf_matrix, columns=K)
	Ndf_matrix=Ndf.values
	Ndf = pd.DataFrame(Ndf_matrix, columns=N)
	grp=[]

	for i in range(0,nc):
		grp.append("C1")
	for i in range(0,nt):
		grp.append("C2")
	
	Kdf=Kdf.astype('int64')
	Ndf=Ndf.astype('int64')
	Ndf=Ndf.replace(0, 1) # tx positive 
	Kdf_r=ro.conversion.py2rpy(Kdf)
	Ndf_r=ro.conversion.py2rpy(Ndf)
	n_threads=ro.conversion.py2rpy(npc)
	alternative = ro.conversion.py2rpy("two.sided")
	verbose = ro.vectors.BoolVector([False])
	results=countdata.bb_test(Kdf_r,Ndf_r,ro.StrVector(grp),alternative,n_threads,verbose)

	localdate = time.strftime('%a %m/%d/%Y')
	localtime = time.strftime('%H:%M:%S')
	logfile.write('# Finished BB test : ' + localdate + ' at: ' + localtime + ' \n')
	apa_df['P-value']=results[results.names.index('p.value')][:].tolist()
	apa_df["P-value"] = pd.to_numeric(apa_df["P-value"])
	apa_df['Adj-Pval'] = (multi.multipletests(apa_df["P-value"], method='fdr_bh', is_sorted=False, returnsorted=False))[1].tolist()		
	apa_df.to_csv(out.rstrip("/")+"/"+key+"_"+feature+"_PAusage_Results.txt",sep="\t",header=True,index=False)
	return(1)


def addTotalSum(sdf,sample):
	# add sample level counts #
	sdf.loc['sum'] = sdf.sum()
	as_list = sdf.loc['sum'].tolist()
	as_list[0]="Sample"
	sdf.loc['sum']=as_list
	# add column #
	sdf[sample+'_Total']=sdf.sum(axis=1)
	return(sdf)

def dropGenes(m,nc,nt):
	returnS=[]
	for i in range(0, len(m)):
		if np.sum(m[i][1:nc+nt]) == 0:
			returnS.append(0)
		elif 0 in list(m[i][1:nc]) and 0 in list(m[i][nc:nc+nt]):
			returnS.append(0)
		elif np.mean(m[i][1:nc]) < 10 and np.mean(m[i][nc:nc+nt]) < 10:
			returnS.append(0)
		elif min(list(m[i][nc+nt:nc+nt+nc])) < 10 and min(list(m[i][nc+nt+nc:nc+nt+nc+nt])) < 10:
			returnS.append(0)
		else:
			returnS.append(1)
	return(returnS)


def PAusage(gene_bed,APAcounts,outDir,prefix,np,logfile):
	out=outDir.strip("/")+"/"+prefix
	# gene id to symbol #
	gene_map=pd.read_csv(gene_bed,sep="\t",header=None,index_col=None)
	gene_map.columns=['chr','start','end','gene_id',"Symbol","strand"]
	gene_map=gene_map[['gene_id',"Symbol"]]
	
	# split filteres APA counts to regions #
	df=pd.read_csv(APAcounts,sep="\t",header=0,index_col=None)
	Mcols=df.columns[2:]
	features=['UTR3','UTR5','Intron','CDS']
	fcounts=[]
	for f in features:
		dftemp=df[df['feature_id'].str.contains(f)].copy()
		dftemp=dftemp.drop(columns=['feature_id'])
		dftemp[dftemp.columns[1:]]=dftemp[dftemp.columns[1:]].apply(pd.to_numeric)
		dftemp = dftemp.groupby(['gene_id'], as_index=False, sort=False).aggregate('sum')
		cols=dftemp.columns[1:]
		cols=[x + "_"+f for x in cols]
		dftemp.columns=['gene_id']+cols
		fcounts.append(dftemp)

	df=fcounts[0].copy()
	for dftemp in fcounts[1:]:
		df=df.merge(dftemp,on=['gene_id'],how='outer')
	df=df.fillna(0)


	# get sample proportions #
	all_cols=df.columns
	sdfs=[]
	for sample in Mcols:

		scols=['gene_id']
		for i in range (0, len(all_cols)):
			if sample in all_cols[i]:
				scols.append(all_cols[i])
		sdf=df[scols].copy()
		sdf=addTotalSum(sdf,sample)
		sdfs.append(sdf)
	

	# process count matrix #
	df = sdfs[0].copy()
	for s in sdfs[1:]:
		df=df.merge(s,on=['gene_id'],how='outer')
	df=df.fillna(0)


	# grouping and test BB#
	controls=Mcols[:int(len(Mcols)/2)]
	treated=Mcols[int(len(Mcols)/2):]
	all_cols=df.columns
	for feature in features:
		control_cols=[];control_tot=[];treated_cols=[];treated_tot=[]
		for i in range (0, len(all_cols)):
			for j in range(0,len(controls)):
				if controls[j] in all_cols[i] and feature in all_cols[i]:
					control_cols.append(all_cols[i])
				if controls[j] in all_cols[i] and "Total" in all_cols[i]:
					control_tot.append(all_cols[i])
			for j in range(0,len(treated)):
				if treated[j] in all_cols[i] and feature in all_cols[i]:
					treated_cols.append(all_cols[i])
				if treated[j] in all_cols[i] and "Total" in all_cols[i]:
					treated_tot.append(all_cols[i])
		# sort #
		control_cols.sort();control_tot.sort();treated_cols.sort();treated_tot.sort()
		sdf=df[['gene_id']+control_cols+treated_cols+control_tot+treated_tot].copy()
		#sdf.to_csv("check1.txt",sep="\t",header=True,index=False)
		sdf['Status']=dropGenes(sdf.values,len(controls),len(treated))
		#sdf.to_csv("check2.txt",sep="\t",header=True,index=False)
		sdf=sdf[sdf['Status']==1]
		sdf=sdf[['gene_id']+control_cols+treated_cols+control_tot+treated_tot]
		sdf["Control_mean"]=sdf[control_cols].mean(axis=1)
		sdf["Treated_mean"]=sdf[treated_cols].mean(axis=1)
		sdf["Control_Tmean"]=sdf[control_tot].mean(axis=1)
		sdf["Treated_Tmean"]=sdf[treated_tot].mean(axis=1)
		sdf["Control"]=sdf["Control_mean"].div(sdf["Control_Tmean"],axis=0)
		sdf["Treated"]=sdf["Treated_mean"].div(sdf["Treated_Tmean"],axis=0)
		sdf[feature+"_FC"]=sdf["Treated"].div(sdf["Control"],axis=0)
		sdf=sdf.merge(gene_map,on=["gene_id"],how="left")
		sdf=sdf[['gene_id']+control_cols+treated_cols+control_tot+treated_tot+['Symbol','Control','Treated',feature+"_FC"]]
		#sdf.to_csv(out+"_"+feature+"_BBinput.txt",sep="\t",index=None) # can go
		s=0
		s=runBBtest(sdf, feature,len(controls), len(treated), outDir,prefix,np,logfile)
		#os.system("rm "+out+"_"+feature+"_BBinput.txt")
		if s!=1:
			return (0)
	return (1)


