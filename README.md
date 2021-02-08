# Placeholder for SISAP 2020 paper
Yu et al. Using a Set of Triangle Inequalities to Accelerate K-means Clustering, SISAP 2020

##File structures:

    .
    ├── Fast-elkan                       # Placeholder for outputs
    │   ├── Makefile          
    │   ├── Readme.md    # Detailed readme for the framework
    |   ├── src                # source code (*.cpp)
    │   └── include        # header (*.h)
    ├── datasets                    # Placeholder for dataset files    
    ├── datapreprocessing.py     # Preprocessing python script
    └── README.md
    
The preprocessing python script works for the following tested datasets at least:
- USPS: [https://www.csie.ntu.edu.tw/~cjlin/libsvmtools/datasets/multiclass.html#usps]
- Skin_nonskin: [https://www.csie.ntu.edu.tw/~cjlin/libsvmtools/datasets/binary.html#skin_nonskin]
- gas sensor: [https://archive.ics.uci.edu/ml/datasets/Gas+Sensor+Array+Drift+Dataset+at+Different+Concentrations]
    -  (semi-auto) please manually remove semicolon in the raw data
    
Usage of preprocessing script:
```
Usage: python dataPreprocessing.py <flag> <path/of/inputDataset> <path/of/outputDataset> <name>, flag 0 for LIBSVM; 1 for UCI repository 
for example: python dataPreprocessing.py 0 datasets/ datasets/ skin_nonskin
```
