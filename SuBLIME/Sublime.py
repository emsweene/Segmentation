import numpy as np
import nibabel as nib
import argparse
import multiprocessing
import config
import utilities
import sys
from scipy import ndimage

LOG = config.getLogger(__name__,"Sublime")

class Sublime(utilities.CaseObject):
	
	def __init__(self,inputDir=None,outputDir=None):
		utilities.CaseObject.__init__(self,inputDir,outputDir)
		if self.__class__.__name__ is "Sublime":
			self.cases = multiprocessing.Queue()
			self.mgr = multiprocessing.Manager()
			self.completedCases = self.mgr.list()
	
	def run(self,inputCases):
		self.cases = self.getInputCases(inputCases)
		self.cases.put(None)
		while True:
			currentCase= self.cases.get()
			if currentCase is None:
				break
				self.DB.close()
			try:
				if len(currentCase) == 2:
					self.runByCasePair(currentCase[0],currentCase[1],self.inputRoot,self.outputRoot)
				else:
					raise TypeError("Array must contain two cases")
			except KeyboardInterrupt:
				sys.exit(1)
				raise
			except:
				LOG.exception(currentCase)
			else:
				self.completedCases.append(currentCase)
		config.defaultOutput("Sublime",self.completedCases)
		return self.completedCases
	
	@classmethod
	def runByCasePair(self,baseline,followup,inputDir,outputDir):
		start = baseline.find("/")
		scanID1 = baseline[start+1:]
		start = followup.find("/")
		MRN = followup[:start]
		scanID2 = followup[start+1:]
		outputDir = outputDir+MRN+"/"+scanID1+"_"+scanID2+"/"
		
		baseline = inputDir+baseline+"/"
		followup = inputDir+followup+"/"
		
		FLAIR_B = baseline+"FLAIRNorm.nii.gz"
		PD_B = baseline+"PDNorm.nii.gz"
		T1_B = baseline+"VolumetricT1Norm.nii.gz"
		T2_B = baseline+"T2Norm.nii.gz"
		FLAIR_F = followup+"FLAIRNorm.nii.gz"
		PD_F = followup+"PDNorm.nii.gz"
		T1_F = followup+"VolumetricT1Norm.nii.gz"
		T2_F = followup+"T2Norm.nii.gz"
		
		self.createProbablilityMap(FLAIR_B,FLAIR_F,PD_B,PD_F,T2_B,T2_F,T1_B,T1_F,outputDir)
		
	@classmethod
	def createProbablilityMap(self,FLAIR_B,FLAIR_F,PD_B,PD_F,T2_B,T2_F,T1_B,T1_F,outputDir,FOV_B=None,FOV_F=None):
		utilities.mkdir_p(outputDir)
		
		intercept = -9.1008
		T1 = 0.5098
		FLAIR = 0.7388
		PD = 0.5606
		T2 = -0.2531
		T1sub = -0.8282
		FLAIRsub = 0.0540
		PDsub = 0.3959
		T2sub = 0.6503
		
		FLAIR_B = nib.load(FLAIR_B)
		FLAIR_F = nib.load(FLAIR_F)
		PD_B = nib.load(PD_B)
		PD_F = nib.load(PD_F)
		T1_B = nib.load(T1_B)
		T1_F = nib.load(T1_F)
		T2_B = nib.load(T2_B)
		T2_F = nib.load(T2_F)
				
		FLAIR_Bdata = FLAIR_B.get_data()
		FLAIR_Fdata = FLAIR_F.get_data()
		PD_Bdata = PD_B.get_data()
		PD_Fdata = PD_F.get_data()
		T1_Bdata = T1_B.get_data()
		T1_Fdata = T1_F.get_data()
		T2_Bdata = T2_B.get_data()
		T2_Fdata = T2_F.get_data()
				
		FLAIR_Sdata = FLAIR_Fdata-FLAIR_Bdata
		PD_Sdata = PD_Fdata-PD_Bdata
		T1_Sdata = T1_Fdata-T1_Bdata
				
		if FOV_B is None:
			T2_Bmask = np.array(T2_Bdata>T2_Bdata[0,0,0],dtype=int)
			FLAIR_Bmask = np.array(FLAIR_Bdata>FLAIR_Bdata[0,0,0],dtype=int)
			PD_Bmask = np.array(PD_Bdata>PD_Bdata[0,0,0],dtype=int)
			T1_Bmask = np.array(T1_Bdata>T1_Bdata[0,0,0],dtype=int)
			FOV_Bdata = T2_Bmask*FLAIR_Bmask*PD_Bmask*T1_Bmask
		else:
			FOV_B = nib.load(FOV_B)
			FOV_Bdata = FOV_B.get_data()
					
		if FOV_F is None:
			T2_Fmask = np.array(T2_Fdata>T2_Fdata[0,0,0],dtype=int)
			FLAIR_Fmask = np.array(FLAIR_Fdata>FLAIR_Fdata[0,0,0],dtype=int)
			PD_Fmask = np.array(PD_Fdata>PD_Fdata[0,0,0],dtype=int)
			T1_Fmask = np.array(T1_Fdata>T1_Fdata[0,0,0],dtype=int)
			FOV_Fdata = T2_Fmask*FLAIR_Fmask*PD_Fmask*T1_Fmask
		else:
			FOV_F = nib.load(FOV_F)
			FOV_Fdata = FOV_F.get_data()
				
		FOV_data = FOV_Fdata*FOV_Bdata
		T2_Sdata = (T2_Fdata-T2_Bdata)*FOV_data
				
		T2_SBlur = ndimage.filters.convolve(T2_Sdata,self.gauss3D(5,3),mode='constant',cval=T2_Sdata[0,0,0])
				
		T2_SBlur_stdev = T2_SBlur.std()
		voxel_select_mask_data = np.array(T2_SBlur>T2_SBlur_stdev,dtype=int)

		voxel_select_mask = nib.Nifti1Image(voxel_select_mask_data,None,header=T2_F.get_header())
		voxel_select_mask.to_filename(outputDir+"SublimeVoxelSelectMask.nii.gz")

		prob_map_raw_data = intercept + T1*T1_Fdata + FLAIR*FLAIR_Fdata + PD*PD_Fdata + T2*T2_Fdata + \
							T1sub*T1_Sdata + FLAIRsub*FLAIR_Sdata + PDsub*PD_Sdata + T2sub*T2_Sdata
		prob_map_raw_data = np.exp(prob_map_raw_data)
		prob_map_raw_data = prob_map_raw_data/(1.0+prob_map_raw_data)

		prob_map_raw = nib.Nifti1Image(prob_map_raw_data,None,header=T2_F.get_header())
		prob_map_raw.to_filename(outputDir+"SublimeProbMapRaw.nii.gz")
		prob_map_data = prob_map_raw_data*voxel_select_mask_data
		prob_map_data = ndimage.filters.convolve(prob_map_data,self.gauss3D(3,3),mode='constant',cval=T2_Sdata[0,0,0])

		prob_map = nib.Nifti1Image(prob_map_data,None,header=T2_F.get_header())
		prob_map.to_filename(outputDir+"SublimeProbMap.nii.gz")


	@classmethod
	def gauss3D(self,kernelsize,sigma):
		center = (kernelsize+1)/2
		kernel = np.zeros((kernelsize,kernelsize,kernelsize))
		covarray = np.array([[sigma,0,0],[0,sigma,0],[0,0,sigma]])
		covdet = abs(np.linalg.det(covarray))
		covinv = np.linalg.inv(covarray)
		for i in xrange(1,kernelsize+1):
			for j in xrange(1,kernelsize+1):
				for k in xrange(1,kernelsize+1):
					point = np.array([i-center,j-center,k-center])
					kernel[i-1,j-1,k-1] = ((2.0*np.pi)**(-3.0/2.0))*np.exp(-0.5*(np.dot(np.dot(np.transpose(point),covinv),point)))/np.sqrt(covdet)
		kernel /= np.sum(kernel)
		return kernel
	
		
if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument("--FLAIRB",type=str,required=True)
	parser.add_argument("--FLAIRF",type=str,required=True)
	parser.add_argument("--T1B",type=str,required=True)
	parser.add_argument("--T1F",type=str,required=True)
	parser.add_argument("--PDB",type=str,required=True)
	parser.add_argument("--PDF",type=str,required=True)
	parser.add_argument("--T2B",type=str,required=True)
	parser.add_argument("--T2F",type=str,required=True)
	parser.add_argument("--FOVB",type=str,default=None)
	parser.add_argument("--FOVF",type=str,default=None)
	parser.add_argument("-o","--outputDir",type=str,required=True)
	args = parser.parse_args()
	S = Sublime()
	S.createProbablilityMap(args.FLAIRB,args.FLAIRF,args.PDB,args.PDF,args.T2B,args.T2F,args.T1B,args.T1F,args.outputDir,args.FOVB,args.FOVF)
