# Optimized Elkan K-means 

## HOW TO RUN THE SOFTWARE:
The driver is designed to take commands from standard input, usually a file
that's been redirected as input:

./driver-standalone <algorithm> <dataset> <number of centers> <output>

You can read the source to find all the possible commands, but here is a summary:
- algorithm -- the name of algorithm, eg. elkan, FB-elkan, MO-elkan to cluster
- dataset --the given path to a file as the dataset 
- number of centers --the given k as initial number of centers
- output --output the result of clustering eg. <centers> outputs results of all centers. Otherwise, output the results of assigned center of each point 
For example:
```
make all
./bin/driver-standalone MO_elkan ../datasets/skin_nonskin.txt 50 centers
```
Please note that yinyang_k-means is incompleted version.

Some notes for using jupyter notebook:
- Data preprocess for livsvmdataset
- data preprocessing.ipythonnotebook
- jupyter notebook > data prepeocessing 

