# PolyASafety Class #

import time, subprocess, os, sys, glob
import pkg_resources,subprocess

class PolyASafety:
	def __init__(self, argsDict, outDir, inputFiles, logfile):
		self.argsDict = argsDict
		self.outDir = outDir
		self.inputFiles = inputFiles
		self.logfile = logfile

	def _logMessage(self, message):
		localdate = time.strftime('%a %m/%d/%Y')
		localtime = time.strftime('%H:%M:%S')
		with open(self.logfile, "a") as fileObj:
			fileObj.write("# " + message + " on: " + localdate + ' at: ' + localtime + ' \n')
			print("# " + message + " on: " + localdate + ' at: ' + localtime)

	def _checkOutDir(self):
		outDirNoSlash = self.outDir.rstrip("/")
		if (os.path.isdir(outDirNoSlash)):
			os.system("rm -R " + outDirNoSlash)
			# subprocess.run(['rm','-R',outDirNoSlash], stderr=subprocess.DEVNULL, shell=False)
		# subprocess.run(['mkdir',outDirNoSlash], stderr=subprocess.DEVNULL, shell=False)
		os.system("mkdir " + outDirNoSlash)

	def _checkInputFiles(self):
		for file in self.inputFiles:
			if os.path.exists(file):
				pass
			else:
				self._logMessage("Error: Cannot Find Input Files")
				return(0)
		self._logMessage("Arguments checked")
		return(1)

	def _startLog(self):
		self._logMessage("Starting PolyA-miner")
		with open(self.logfile, "a") as fileObj:
			fileObj.write('# *********** Arguments ************** \n')
			for k in self.argsDict:
				if k in ["c1","c2"]:
					fileObj.write("\t"+k+' : ' + "\t".join(self.argsDict[k])+'\n')
				else:
					fileObj.write("\t"+k+' : ' + str(self.argsDict[k])+'\n')
			fileObj.write('# *********** ********* ************** \n')

	def _checkDependency(self):
		dependency=['pandas','numpy','Cython','pybedtools','scipy','sklearn','statsmodels','uuid','rpy2']
		installed_packages = " ".join([str(d) for d in pkg_resources.working_set])
		for package in dependency:
			if package in installed_packages:
				pass
			else:
				self._logMessage("Package " + package + " not found. Installing ....")
				try:		
					os.system("pip3 install " + package)
				except:
					self._logMessage("Failed installing " + package + " ....")
					return(0)

		if subprocess.call(['which','featureCounts'], stderr=subprocess.PIPE, stdout=subprocess.PIPE, shell=False) ==1:
			self._logMessage("Err: Missing featureCounts.")
			return(0)
		
		if subprocess.call(['which','bedtools'], stderr=subprocess.PIPE, stdout=subprocess.PIPE, shell=False) == 1:
			self._logMessage("Err: Missing bedtools. Install bedtools.")
			return(0)
			
		import rpy2
		from rpy2.robjects.packages import importr
		rpy2.robjects.r['options'](warn=-1)

		try:
			countdata=importr('countdata')
		except:
			self._logMessage("Err: Missing R package 'countdata'.")
			return(0)

		return(1)

	def runInitialSafetyChecks(self):
		self._checkOutDir()
		self._startLog()
		if (self._checkDependency() == 0):
			self._logMessage("Err: Fix dependencies and run ...")
			return(0)
		# if (self._checkInputFiles() == 0):
		# 	return(0)

	def tidyUp(self, keep):
		cleanupFiles = glob.glob(self.outDir + "/*")
		cleanupFiles = [i for i in cleanupFiles if i not in keep]
		for cleanupFile in cleanupFiles:
			if "Stranded" not in os.path.basename(cleanupFile):
				if "Graphics" not in os.path.basename(cleanupFile):
					os.remove(cleanupFile) 
		self._logMessage("Finished PolyAMiner")
		# Tidy up  #
		# files=[args.o.rstrip("/")+"/LibSize.txt",args.o.rstrip("/")+"/"+args.outPrefix+"_Gene_Stats.txt",args.o.rstrip("/")+"/"+args.outPrefix+'.APSitesDB.bed',args.o.rstrip("/")+"/"+args.outPrefix+"_APA.CountMatrix.GFil.PA.PR.txt",args.o.rstrip("/")+"/"+args.outPrefix+"_denovoAPAsites.bed"]
		# log=subprocess.run(['rm','-R']+files,stderr=subprocess.DEVNULL,shell=False)
		# try:
		# 	subprocess.run(['rm',dummy_refPA]+files,stderr=subprocess.DEVNULL,shell=False)
		# except: 
		# 	pass



		
