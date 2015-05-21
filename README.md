# MotifExtraction(https://github.com/THUIR/MotifExtraction)

MotifExtraction is the implementation of motif extraction algorithm: Yiqun Liu, Ye Chen, Jinhui Tang, Jiashen Sun, Min Zhang, Shaoping Ma, Xuan Zhu. Different Users, Different Opinions: Predicting Search Satisfaction with Mouse Movement Information. SIGIR (2015)

A *motif* is a frequently-appeared sequence of mouse positions, which can be used for satisfaction prediction.

This project is aimed to implement motif extraction algorithm and intended to be easy-to-read and easy-to-modify. If it's not, please let me know how to improve it :)

# Files
## README.md
This file.
 
## bin/
Directory with the Java class file.

## dataset/
Directory with the sample dataset.

- MouseLog: all search session mouse logs used in paper
- Satisfaction_Scores: satisfaction scores from both users and external assessors used in paper

## src/
Directory with the Java code.

# Format of the Input Data 
Examples can be found under `dataset/` (tab-separated).

- for Data Pre-Processing:
  - raw_data: organized as in MouseLog, each file contains the mouse log of a search session
  - sat_data: organized as SampleSatData, each line contains the file name and its satisfaction score, 2 means satisfactory and 1 means dissatisfactory.

- for Motif Extraction
  - train_data, test_data: organized as in CandidateData, each file contains the motif candidates extracted from the mouse log of a search session. In each file, different lines represent different motifs and its satisfaction scores, separated by '\t'. 2 means satisfactory and 1 means dissatisfactory.
  - e.g.: x1,y1;x2,y2;x3,y3;...;x10,y10;  2.0

# Usage
- for Data Pre-Processing:
  - java SatPredict DataProcess dataDirectory SatFile outputDirectory windowlen
  - e.g.: java SatPredict DataProcess ../dataset/MouseLog/ ../dataset/SampleSatData ../dataset/CandidateData/ 3000

- for Motif Extraction:
    - java SatPredict MotifExtraction method ['Frequency', 'Distance', 'Distribution'] R trainDataDirectory testDataDirectory motifNum resultFile
    - e.g.: java SatPredict MotifExtraction Frequency 3000 ../dataset/CandidateData/ ../dataset/CandidateData/ 5 ../dataset/res

# Output
- for Data Pre-Processing:
  - result files will be generated in outputDirectory, which can be directly used as train_data or test_data in Motif Extraction

- for Motif Extraction:
  - result file will be automatically generated, each line represent a search session instance, containing filename, n motif features and satisfaction score, separated by '\t'
  - e.g.: "01_1  1.0 1.0 1.0 1.0 1.0 2.0" means filename is 01_1, all five motif features are 1.0 and the satisfaction score is 2.0





