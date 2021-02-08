import sys
from sklearn.datasets import load_svmlight_file
from sklearn.datasets import load_svmlight_files

def libsvmPreprocess(pathInput, pathOutput, name):

	# download the raw dataset from libsvm website, e.g., skin_nonskin, HIGGS.bz2

	# at first load the raw dataset
	x1 = load_svmlight_file(pathInput+name)

	B = x1[0].toarray()
	#for toarray
	line_count = 0
	with open(pathOutput+name+'.txt', 'w') as myfile:
		myfile.write(str(len(B))+" "+str(len(B[0]))+"\n") # number data points / number of features
		for i,line in enumerate(B):

			line=str(line)+"\n"
			myfile.writelines(line.replace('[','').replace(']',''))
	myfile.close()


def insert(originalfile,string):
	with open(originalfile,'r') as f:
		with open('./datasets/gassensor_clean.txt','w') as f2:
			f2.write(string)
			f2.write(f.read())

def datasetPreprocess(pathInput, pathOutput, name):

	# if the raw dataset is downloaded from UCI, e.g., gassensor
	x1,y1,x2,y2,x3,y3,x4,y4,x5,y5,x6,y6,x7,y7,x8,y8,x9,y9,x10,y10=load_svmlight_files(("./datasets/gassensor/batch1.dat",
											   "./datasets/gassensor/batch2.dat",
											  "./datasets/gassensor/batch3.dat",
											  "./datasets/gassensor/batch4.dat",
											  "./datasets/gassensor/batch5.dat",
											  "./datasets/gassensor/batch6.dat",
											  "./datasets/gassensor/batch7.dat",
											  "./datasets/gassensor/batch8.dat",
											  "./datasets/gassensor/batch9.dat",
											  "./datasets/gassensor/batch10.dat"), multilabel=True)
	B=x1.toarray().tolist()+x2.toarray().tolist()+x3.toarray().tolist()+x4.toarray().tolist()+x5.toarray().tolist()+x6.toarray().tolist()+x7.toarray().tolist()+x8.toarray().tolist()+x9.toarray().tolist()+x10.toarray().tolist()
	# x1,y1,x2,y2=load_svmlight_files(("datasets/gassensor/batch1.dat","datasets/gassensor/batch2.dat"))
	# B=x1.toarray().tolist()+x2.toarray().tolist()
	print(test)

	#for tolist
	import re
	#B = x1[0].toarray().tolist()

	with open('./datasets/gassensor.txt', 'w') as myfile:
		for i,line in enumerate(B):#myfile.write("%s\n" % line)
			line=str(line)+"\n"
			myfile.writelines(line.replace('[','').replace(']','').replace(',',' '))

		#os.rename('newfile.txt',originalfile)
	file=insert('./datasets/gassensor.txt',"14000 128\n")

def main():
	args = sys.argv
	if len(args) < 5:
		print ("Usage: python dataPreprocessing.py <flag> <path/of/inputDataset> <path/of/outputDataset> <name>, flag 0 for LIBSVM; 1 for UCI repository ")
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
