# v1.2
# Compute Gene level PolA #

import numpy as np, pandas as pd
import math,concurrent.futures as cf
from numpy import array, dot
from numpy.linalg import norm

def count2prop(sdf):
	X = sdf.values
	X_sums = np.sum(X, axis=0)
	for i in range(X.shape[0]):
		for j in range(X.shape[1]):
			X[(i, j)] = X[(i, j)] / X_sums[j]
	return (X)

def count2proplist(lst):
	tot = sum(lst)
	if tot > 0:
		for i in range(0, len(lst)):
			lst[i] = lst[i] / tot

		return (lst)
	temp = []
	for i in range(0, len(lst)):
		temp.append(0)

	return (temp)

def ComputeLengtheningScoreDeprecated(W_g, strand):
	vd = np.asarray((np.zeros(shape=(1, W_g.shape[0])))[0])
	vs = np.asarray((np.zeros(shape=(1, W_g.shape[0])))[0])
	#for al in range(0,len(vs)):
	#	vs[al]=1-((1/float(len(vs)))*al)
	if strand == '+':
		vs[0] = 1
		vs[-1] = 0
		vd[0] = 0
		vd[-1] = 1
		#vd = vs[::-1]
	if strand == '-':
		vs[0] = 0
		vs[-1] = 1
		vd[0] = 1
		vd[-1] = 0
		#vd=vs
		#vs=vd[::-1]	
	try:
		a = count2proplist(list(W_g[:, 0].flat))
		b = count2proplist(list(W_g[:, 1].flat))
	except:
		exit()

	cos_theta_cp = dot(vs, a) / (norm(vs) * norm(a))
	cos_theta_cd = dot(vd, a) / (norm(vd) * norm(a))
	cos_theta_tp = dot(vs, b) / (norm(vs) * norm(b))
	cos_theta_td = dot(vd, b) / (norm(vd) * norm(b))
	tp = norm(b) * cos_theta_tp
	cp = norm(a) * cos_theta_cp
	td = norm(b) * cos_theta_td
	cd = norm(a) * cos_theta_cd
	try:
		paindex=math.log((cp/(cp+cd))/(tp/(tp+td)),2)
	except:
		paindex="-inf"	
	return (cp, tp, cd, td, paindex)

def ComputeLengtheningScore(W_g, strand, APAsites, optionVectorToZero, optionPAIndexDivideByDenominator):
	APAsitesdf = pd.DataFrame(APAsites)
	vd = np.asarray((np.zeros(shape=(1, W_g.shape[0])))[0])
	vs = np.asarray((np.zeros(shape=(1, W_g.shape[0])))[0])

	# if strand == '+':
	# 	vs[0] = 1
	# 	for i in range(1,vs.size):
	# 		vs[i] = vs[i-1] - 1/(vs.size)
	# 	vd[-1] = 1
	# 	for i in range(0,vd.size):
	# 		if i == 0:
	# 			vd[i] = 1/(vd.size)
	# 		else:
	# 			vd[i] = vd[i-1] + 1/(vd.size)
	# if strand == '-':
	# 	vs[-1] = 1
	# 	for i in range(0,vs.size):
	# 		if i == 0:
	# 			vs[i] = 1/(vs.size)
	# 		else:
	# 			vs[i] = vs[i-1] + 1/(vs.size)
	# 	vd[0] = 1
	# 	for i in range(1,vd.size):
	# 		vd[i] = vd[i-1] - 1/(vd.size)

	if strand == "+":
		APAsitesdf[['Chr','Start', "End", "Strand"]] = APAsitesdf[0].str.split('_',expand=True)
		CPAS_MagnitudeOfLength = abs(int(APAsitesdf["Start"].iat[0]) - int(APAsitesdf["End"].iat[-1]))
		APAsitesdf["CPAS_Percentile_vd"] = abs(APAsitesdf["End"].astype(float) - float(APAsitesdf["Start"].iat[0])) / CPAS_MagnitudeOfLength
		APAsitesdf["CPAS_Percentile_vs"] = 1 - abs(APAsitesdf["Start"].astype(float) - float(APAsitesdf["Start"].iat[0])) / CPAS_MagnitudeOfLength
		vd = APAsitesdf.CPAS_Percentile_vd.tolist()
		vs = APAsitesdf.CPAS_Percentile_vs.tolist()
	if strand == "-":
		APAsitesdf[['Chr','Start', "End", "Strand"]] = APAsitesdf[0].str.split('_',expand=True)
		CPAS_MagnitudeOfLength = abs(int(APAsitesdf["Start"].iat[0]) - int(APAsitesdf["End"].iat[-1]))
		APAsitesdf["CPAS_Percentile_vs"] = abs(APAsitesdf["End"].astype(float) - float(APAsitesdf["Start"].iat[0])) / CPAS_MagnitudeOfLength
		APAsitesdf["CPAS_Percentile_vd"] = 1 - abs(APAsitesdf["Start"].astype(float) - float(APAsitesdf["Start"].iat[0])) / CPAS_MagnitudeOfLength
		vd = APAsitesdf.CPAS_Percentile_vd.tolist()
		vs = APAsitesdf.CPAS_Percentile_vs.tolist()

	try:
		a = count2proplist(list(W_g[:, 0].flat))
		b = count2proplist(list(W_g[:, 1].flat))
	except:
		exit()

	cos_theta_cp = dot(vs, a) / (norm(vs) * norm(a))
	cos_theta_cd = dot(vd, a) / (norm(vd) * norm(a))
	cos_theta_tp = dot(vs, b) / (norm(vs) * norm(b))
	cos_theta_td = dot(vd, b) / (norm(vd) * norm(b))
	tp = norm(b) * cos_theta_tp
	cp = norm(a) * cos_theta_cp
	td = norm(b) * cos_theta_td
	cd = norm(a) * cos_theta_cd
	try:
		paindex = math.log((cp/tp),2)
		# paindex=math.log((cp/(cp+cd))/(tp/(tp+td)),2)
	except:
		paindex = float('-inf')
	a = ",".join(list(map(str,list(a))))
	b = ",".join(list(map(str,list(b))))	
	return (a, b, cp, tp, cd, td, paindex)

def ComputeVectorPro(x):
	sdf = x[1][x[1]['gene_id'] == x[0]]
	x[1] = 'NULL'
	if sdf.shape[0] > 1:
		APAsites = list(sdf['feature_id'])
		strand = APAsites[0].split('_')[-1]
		Samples = list(sdf.columns[2:])
		sdf = sdf[sdf.columns[2:]]
		osdf = sdf.copy()
		sdf = count2prop(sdf)
		sdf[np.isnan(sdf)] = 0.0
		APAmean = np.hstack((sdf[:, 0:x[2]].sum(axis=1, keepdims=1), sdf[:, x[2]:x[2] + x[3]].sum(axis=1, keepdims=1)))
		#Try 4 iterations: (1) optionVectorToZero=F, optionPAIndexDivideByDenominator=F, (2) optionVectorToZero=F, optionPAIndexDivideByDenominator=T, (3) optionVectorToZero=T, optionPAIndexDivideByDenominator=F, (4) optionVectorToZero=T, optionPAIndexDivideByDenominator=T
		vector_a, vector_b, CP, TP, CD, TD, PolyAIndex_length = ComputeLengtheningScore(APAmean, strand, APAsites, optionVectorToZero = False, optionPAIndexDivideByDenominator = False)
		return ([x[0], str(vector_a), str(vector_b), str(CP), str(TP), str(CD), str(TD), str(PolyAIndex_length)])
			 
def VectorPro(outDir, fkey, nc, nt, npc, matrix, LS):
	df = pd.read_csv(matrix, sep='\t', index_col=None)
	df[['feature_id','FM']]=df.feature_id.str.split('@', expand=True)
	df= df.drop(columns=['FM'])
	
	genes = list(set(list(df['gene_id'])))
	dflb = pd.read_csv(LS, sep='\t', index_col=None)
	cols = df.columns[2:]

	for i in range(0, len(cols)):
		df[cols[i]] = df[cols[i]] / float(dflb[cols[i]][0]) * 1000000.0

	passlist = []
	#genes = ['VMA21', 'LAMC1', 'PAK1', 'PAK2']
	for gene in genes:
		passlist.append([gene, df, nc, nt])

	with cf.ProcessPoolExecutor(max_workers=npc) as (executor):
		result = list(executor.map(ComputeVectorPro, passlist))
	result=list(result) 
	fw = open(outDir + '/' + fkey + '_PolyA-miner.Results.txt', 'w')
	fw.write('Gene\tVector_a\tVector_b\tCP\tTP\tCD\tTD\tPolyAIndex\n')
	for p in result:
		try:
			fw.write(('\t').join(p) + '\n')
		except:
			pass
	fw.close()
	return (1)