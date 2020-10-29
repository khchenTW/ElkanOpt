import sys
from sklearn.datasets import load_svmlight_file
from sklearn.datasets import load_svmlight_files

def libsvmPreprocess(pathInput, pathOutput, name):

	# download the raw dataset from libsvm website, e.g., HIGGS.bz2

	# at first load the raw dataset
	x1 = load_svmlight_file(pathInput+name)

	B = x1[0].toarray()
	#for toarray
	with open(pathOutput+name+'.txt', 'w') as myfile:
		for i,line in enumerate(B):
		#myfile.write("%s\n" % line)
		
			line=str(line)+"\n"	
			myfile.writelines(line.replace('[','').replace(']',''))

	#save memory
	f = open(pathOutput+name+'.txt','r')
	temp = f.read()
	f.close()

	f = open(pathOutput+name+'_clean.txt', 'w')
	# KHCHEN
	# this is not automatic yet
	f.write("11000000 28\n")

	f.write(temp)
	f.close()

def insert(originalfile,string):
	with open(originalfile,'r') as f:
		with open('/home/qiao/Downloads/datasets/gassensor_clean.txt','w') as f2: 
			f2.write(string)
			f2.write(f.read())

def datasetPreprocess(pathInput, pathOutput, name):
	
	# if the raw dataset is downloaded from , e.g., gassensor
	x1,y1,x2,y2,x3,y3,x4,y4,x5,y5,x6,y6,x7,y7,x8,y8,x9,y9,x10,y10=load_svmlight_files(("/home/qiao/Desktop/fast-kmeans-master/datasets/gassensor/batch1.dat",
											   "/home/qiao/Desktop/fast-kmeans-master/datasets/gassensor/batch2.dat",
											  "/home/qiao/Desktop/fast-kmeans-master/datasets/gassensor/batch3.dat",
											  "/home/qiao/Desktop/fast-kmeans-master/datasets/gassensor/batch4.dat",
											  "/home/qiao/Desktop/fast-kmeans-master/datasets/gassensor/batch5.dat",
											  "/home/qiao/Desktop/fast-kmeans-master/datasets/gassensor/batch6.dat",
											  "/home/qiao/Desktop/fast-kmeans-master/datasets/gassensor/batch7.dat",
											  "/home/qiao/Desktop/fast-kmeans-master/datasets/gassensor/batch8.dat",
											  "/home/qiao/Desktop/fast-kmeans-master/datasets/gassensor/batch9.dat",
											  "/home/qiao/Desktop/fast-kmeans-master/datasets/gassensor/batch10.dat"))
	B=x1.toarray().tolist()+x2.toarray().tolist()+x3.toarray().tolist()+x4.toarray().tolist()+x5.toarray().tolist()+x6.toarray().tolist()+x7.toarray().tolist()+x8.toarray().tolist()+x9.toarray().tolist()+x10.toarray().tolist()


	#for tolist
	import re
	#B = x1[0].toarray().tolist()

	with open('/home/qiao/Desktop/fast-kmeans-master/datasets/gassensor.txt', 'w') as myfile:
		for i,line in enumerate(B):#myfile.write("%s\n" % line)
			line=str(line)+"\n"
			myfile.writelines(line.replace('[','').replace(']','').replace(',',' '))

		#os.rename('newfile.txt',originalfile)
	file=insert('/home/qiao/Desktop/fast-kmeans-master/datasets/gassensor.txt',"14000 128\n")

def main():
	args = sys.argv
	if len(args) < 4:
		print ("Usage: python dataPreprocessing.py <flag> <path/of/inputDataset> <path/of/outputDataset> <name>, flag 0 is for LIBSVM")
		return
	mode = int(args[1])
	pathI = str(args[2])
	pathO = str(args[3])
	name = str(args[4])
	
	if mode == 0:
		libsvmPreprocess(pathI, pathO, name)
	else:
		datasetPreprocess(pathI, pathO, name)

if __name__=="__main__":
	main()
