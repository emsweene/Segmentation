#!/usr/bin/env python

import subprocess
import os
import shutil
import sys
import glob

class msg(object):
	__quiet = False
	__debug = False
	OK = 0
	WARNING = 1
	ERROR = 2
	SKIPPED = 3
	defaultColor = "\033[0m"
	okColor = "\033[32m"
	skippedColor = '\033[34m'
	warningColor = "\033[35m"
	errorColor = "\033[1;31m"
	bold = '\033[1;34m'
	
	@classmethod
	def printMsg(cls,status,message):
		if status == cls.OK:
			if cls.__quiet:
				return
			else:
				print "".join((cls.okColor,"[OK] ",cls.defaultColor,message))
		elif status == cls.WARNING:
			if cls.__quiet:
				return
			else:
				print "".join((cls.warningColor,"[WARNING] ",cls.defaultColor,message))
		elif status == cls.ERROR:
			print "".join((cls.errorColor,"[ERROR] ",cls.defaultColor,message))
		elif status == cls.SKIPPED:
			if cls.__quiet:
				return
			else:
				print "".join((cls.skippedColor,"[SKIPPED] ",cls.defaultColor,message))
		else:
			return
	
	@classmethod
	def printDebug(cls,message):
		if cls.__debug:
			print "".join((cls.bold,"[DEBUG] ",cls.defaultColor,message))
		else:
			return
	
	@classmethod
	def debugOn(cls):
		cls.__debug = True

class outputLog(object):
	stdout = []
	stderr = []
	
	def writeSTDOUT(self):
		pass
	
	def writeSTDERR(self):
	    pass
	
	def writeLog(self):
	    pass

def fslmaths(*args):
	cmd = ["fslmaths"]
	for arg in args:
	    cmd.append(arg)
	subprocess.call(cmd,stdout=logFile,stderr=errorFile)
	
def fslmask(inputFile,maskFile,outputFile):
	fslmaths(inputFile,"-mas",maskFile,outputFile)

def fslnorm(inputFile,maskFile,outputFile):
	mean = subprocess.Popen(["fslstats",inputFile,"-M"],stdout=subprocess.PIPE,stderr=errorFile)
	mean = mean.stdout.read()
	stddev = subprocess.Popen(["fslstats",inputFile,"-S"],stdout=subprocess.PIPE,stderr=errorFile)
	stddev = stddev.stdout.read()
	fslmaths(inputFile,"-sub",str(mean).rstrip(),"-div",str(stddev).rstrip(),"-mas",maskFile,outputFile)

def fslsmooth(inputFile,maskFile,outputFile,size):
    fslmaths(inputFile,"-mas",maskFile,"-s",str(size),outputFile)

def fslsmoothdiv(inputFile,divFile,maskFile,outputFile):
    fslmaths(inputFile,"-div",divFile,"-mas",maskFile,outputFile)

def gzipOutput():
	cmd = ["gzip","-f"]
	files = glob.glob(outputDir+"*.nii")
	for fileName in files:
	    cmd.append(fileName)
	subprocess.call(cmd,stdout=logFile,stderr=errorFile)

	
if(__name__ == "__main__"):	
	msg.debugOn()
	
	pipelineHome = "/Users/deweybe/Desktop/PyOASIS/"
	MRN = "04881278"
	scanID = "20120822-05845"
	inputDir = pipelineHome + "Input/" + MRN + "/" + scanID +"/"
	outputDir = pipelineHome + "Output/" + MRN + "/" + scanID +"/"
	logDir = outputDir+"logs/"
	OASIS_Step1 = pipelineHome + "PyOASIS_Step1.R"
	OASIS_Step2 = pipelineHome + "PyOASIS_Step2.R"
	
	try:
		os.makedirs(outputDir)
	except OSError:	
		overwrite = raw_input("Output folder already exists. Would you like to overwrite? (y,n) ")
		if str(overwrite).upper() == "Y":
			shutil.rmtree(outputDir)
			os.makedirs(outputDir)
			os.makedirs(logDir)
			logFile = open(logDir+"stdOut.txt","w")
			errorFile = open(logDir+"stdErr.txt","w")
		else:
			check = raw_input("Would you like to resume a previous session? (y,n) ")
			if str(check).upper() == "Y":
				logFile = open(logDir+"stdOut.txt","a")
				errorFile = open(logDir+"stdErr.txt","a")
			else:
				msg.printMsg(msg.ERROR,"Unable to create output directory. Directory already exists")
				sys.exit()
	else:
		os.makedirs(logDir)
		logFile = open(logDir+"stdOut.txt","w")
		errorFile = open(logDir+"stdErr.txt","w")
	
	try:
		#Filename Variables
		CSFMask = outputDir+"OASIS_CSFMask.nii.gz"
		
		scanTypes = ["FLAIR","PD","T2","VolumetricT1"]
		kernelSizes = [10,20]
		
		# Erodes the brain mask
		if os.path.exists(outputDir+"OASIS_DuraMask.nii.gz"):
			msg.printMsg(msg.SKIPPED,"Eroded DuraMask already exists.")
		else:
			fslmaths(inputDir+"DuraMask.nii.gz","-kernel","box","5x5x5","-ero",outputDir+"OASIS_DuraMask.nii.gz")
		
		# Mask the FLAIR
		if os.path.exists(outputDir+"OASIS_FLAIRStrip.nii.gz"):
			msg.printMsg(msg.SKIPPED,"DuraMasked FLAIR already exists.")
		else:
			fslmask(inputDir+"FLAIRStrip.nii.gz",outputDir+"OASIS_DuraMask.nii.gz",outputDir+"OASIS_FLAIRStrip.nii.gz")
			
		## Exclude the bottom 10 percent of the FLAIR (to remove CSF)
		if os.path.exists(CSFMask):
			msg.printMsg(msg.SKIPPED,"CSF Mask already exists.")
		else:
			thres = subprocess.Popen(["fslstats",outputDir+"OASIS_FLAIRStrip.nii.gz","-P","15"],stdout=subprocess.PIPE,stderr=errorFile)
			thres = thres.stdout.read()
			fslmaths(outputDir+"OASIS_FLAIRStrip.nii.gz","-thr",str(thres).rstrip(),"-bin",CSFMask,"-odt","int")
		msg.printMsg(msg.OK,"CSF Mask Creation Complete")
		
		# Full brain normalization of the images
		for scan in scanTypes:
			if os.path.exists(outputDir+"OASIS_"+scan+"Strip_masked.nii.gz"):
				msg.printMsg(msg.SKIPPED,"CSF Masked "+scan+" already exists.")
			else:
				fslmask(inputDir+scan+"Strip.nii.gz",CSFMask,outputDir+"OASIS_"+scan+"Strip_masked.nii.gz")
		
		for scan in scanTypes:
			if os.path.exists(outputDir+"OASIS_"+scan+"Norm.nii.gz"):
				msg.printMsg(msg.SKIPPED,"Normalized "+scan+" already exists.")
			else:
				fslnorm(outputDir+"OASIS_"+scan+"Strip_masked.nii.gz",CSFMask,outputDir+"OASIS_"+scan+"Norm.nii.gz")
		msg.printMsg(msg.OK,"Full Brain Normalization Complete")
	
		## Create first smoothed volumes (kernel window size of 10 and 20)
		for kernel in kernelSizes:
			if os.path.exists(outputDir+"OASIS_CSFMask_blur1_"+str(kernel)+".nii.gz"):
				msg.printMsg(msg.SKIPPED,"Blurred CSF Mask with kernel size "+str(kernel)+" already exists.")
			else:
				fslmaths(CSFMask,"-s",str(kernel),outputDir+"OASIS_CSFMask_blur1_"+str(kernel)+".nii.gz")
	
			for scan in scanTypes:
				if os.path.exists(outputDir+"OASIS_"+scan+"Norm_blur1_"+str(kernel)+".nii.gz"):
					msg.printMsg(msg.SKIPPED,"First Blurred "+scan+" with kernel size "+str(kernel)+" already exists.")
				else:
					fslsmooth(outputDir+"OASIS_"+scan+"Norm.nii.gz",CSFMask,outputDir+"OASIS_"+scan+"Norm_blur1_"+str(kernel)+".nii.gz",kernel)
				if os.path.exists(outputDir+"OASIS_"+scan+"Norm_blur1_"+str(kernel)+"_div.nii.gz"):
					msg.printMsg(msg.SKIPPED,"First Divided and Blurred "+scan+" with kernel size "+str(kernel)+" already exists.")
				else:
					fslsmoothdiv(outputDir+"OASIS_"+scan+"Norm_blur1_"+str(kernel)+".nii.gz",outputDir+"OASIS_CSFMask_blur1_"+str(kernel)+".nii.gz",CSFMask,outputDir+"OASIS_"+scan+"Norm_blur1_"+str(kernel)+"_div.nii.gz")
		msg.printMsg(msg.OK,"First Image Smoothing Complete")
		
		if(os.path.exists(outputDir+"OASIS_prob_map_1.nii.gz") and os.path.exists(outputDir+"OASIS_segment_run1.nii.gz")):
			msg.printMsg(msg.SKIPPED,"Probability Map and Segmentation from Step 1 already exist.")
		else:
			Rcmd1 = ["Rscript",OASIS_Step1,CSFMask]
			fileTypes = ["Norm.nii.gz","Norm_blur1_10_div.nii.gz","Norm_blur1_20_div.nii.gz"]
			for fileType in fileTypes:
				for scan in scanTypes:
					Rcmd1.append(outputDir+"OASIS_"+scan+fileType)
			Rcmd1.append(outputDir)
			subprocess.call(Rcmd1,stdout=logFile,stderr=errorFile)
			gzipOutput()
		msg.printMsg(msg.OK,"First Segmentation Step Complete")
		
		## Dilate Step 1 Lesion Segmentation
		if os.path.exists(outputDir+"OASIS_segment_dilate.nii.gz"):
			msg.printMsg(msg.SKIPPED,"Dilated Lesion Segmentation already exists.")
		else:	
			fslmaths(outputDir+"OASIS_segment_run1.nii.gz","-kernel","boxv","5","-dilM",outputDir+"OASIS_segment_dilate.nii.gz")
		
		## Apply CSF Mask to Dilated Lesion Segmentation
		LesionMask = outputDir+"OASIS_lesion_segment_dilate.nii.gz"
		if os.path.exists(LesionMask):
			msg.printMsg(msg.SKIPPED,"Masked Dilated Lesion Segmentation already exists.")
		else:	
			fslmaths(CSFMask,"-sub",outputDir+"OASIS_segment_dilate.nii.gz","-mas",CSFMask,LesionMask)
		msg.printMsg(msg.OK,"First Leasion Mask Creation Complete")
		
		## Create second smoothed volumes (kernel window size of 10 and 20)
		for kernel in kernelSizes:
			if os.path.exists(outputDir+"OASIS_SegmentDilate_blur2_"+str(kernel)+".nii.gz"):
				msg.printMsg(msg.SKIPPED,"Second Blurred CSF Mask with kernel size "+str(kernel)+" already exists.")
			else:
				fslmaths(LesionMask,"-s",str(kernel),outputDir+"OASIS_SegmentDilate_blur2_"+str(kernel)+".nii.gz")
	
			for scan in scanTypes:
				if os.path.exists(outputDir+"OASIS_"+scan+"Norm_blur2_"+str(kernel)+".nii.gz"):
					msg.printMsg(msg.SKIPPED,"Blurred "+scan+" with kernel size "+str(kernel)+" already exists.")
				else:
					fslsmooth(outputDir+"OASIS_"+scan+"Norm.nii.gz",LesionMask,outputDir+"OASIS_"+scan+"Norm_blur2_"+str(kernel)+".nii.gz",kernel)
				if os.path.exists(outputDir+"OASIS_"+scan+"Norm_blur2_"+str(kernel)+"_div.nii.gz"):
					msg.printMsg(msg.SKIPPED,"Second Divided and Blurred "+scan+" with kernel size "+str(kernel)+" already exists.")
				else:
					fslsmoothdiv(outputDir+"OASIS_"+scan+"Norm_blur2_"+str(kernel)+".nii.gz",outputDir+"OASIS_SegmentDilate_blur2_"+str(kernel)+".nii.gz",CSFMask,outputDir+"OASIS_"+scan+"Norm_blur2_"+str(kernel)+"_div.nii.gz")
		msg.printMsg(msg.OK,"Second Image Smoothing Complete")
		
		if(os.path.exists(outputDir+"OASIS_segment_run2_1.nii.gz") and os.path.exists(outputDir+"OASIS_segment_run2_2.nii.gz") and os.path.exists(outputDir+"OASIS_segment_run2_3.nii.gz") and os.path.exists(outputDir+"OASIS_segment_run2_4.nii.gz") and os.path.exists(outputDir+"OASIS_segment_run2_5.nii.gz")):
			msg.printMsg(msg.SKIPPED,"Segmentations from Step 2 already exist.")
		else:
			Rcmd2 = ["Rscript",OASIS_Step2,CSFMask]
			fileTypes = ["Norm.nii.gz","Norm_blur2_10_div.nii.gz","Norm_blur2_20_div.nii.gz"]
			for fileType in fileTypes:
				for scan in scanTypes:
					Rcmd2.append(outputDir+"OASIS_"+scan+fileType)
			Rcmd2.append(outputDir)
			subprocess.call(Rcmd2,stdout=logFile,stderr=errorFile)
			gzipOutput()
		msg.printMsg(msg.OK,"Second Segmentation Step Complete")
		
		shutil.rmtree(outputDir+"temp/")
		msg.printMsg(msg.OK,"OASIS Segmentation Complete")
		logFile.close()
		errorFile.close()
	except KeyboardInterrupt:
		logFile.close()
		errorFile.close()
		msg.printMsg(msg.ERROR,"OASIS Program Terminated")