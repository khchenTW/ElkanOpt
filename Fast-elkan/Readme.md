# Optimized Elkan K-means 

## HOW TO RUN THE SOFTWARE:
The driver is designed to take commands from standard input, usually a file
that's been redirected as input:

./driver-standalone <algorithm> <dataset> <number of centers> <output>

You can read the source to find all the possible commands, but here is a summary:
- algorithm --use the name of algorithm, eg. elkan, FB-elkan, MO-elkan to cluster
- dataset --use the given path name to a file as the dataset 
- number of centers --use the given number as initial number of centers
- output --output the result of clustering eg. <centers> outputs results of all centers. Otherwise, output the results of assigned center of each point 
For example:
```
make all
cd bin
./driver-standalone MO_elkan ../datasets/skin_nonskin_clean.txt 50 centers
```

Some notes for the past:
- Data preprocess for livsvmdataset
- data preprocessing.ipythonnotebook
- jupyter notebook > data prepeocessing 
- Noted that yinyang_k-means is incompleted version
   


