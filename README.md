# Masters-Thesis-Project
This thesis proposes a general-purpose data science framework using a supervised learning technique.
# Exploring Vector Attributes of Maxtree, Component tree using Supervised Machine Learning Technique

Description: In this master thesis we focus on exploring vector attributes to detect lung tumors using supervised machine learning technique namely Random Forest. You can find complete work described in the final report "mCS_2023_PoojaG" uploaded here.

In this repository you can find following files in respective directories :
| Items (Scripts/data/results) | Filename/directory |
| ------ | ------ |
| Supervised Trained Model      | Model/randomforest.m      |
|   Testing the Model using 22 whole test set volumes  |  Model/Testrnd1.m      | 
|   Testing the Model using single test set volume  |  Model/Testrnd2.m      | 
|  Filtering of classified training set and test volumes  |    Performance Evaluation/nifti1_Segment_filter.c  | 
|  Reading and Comparing Classified Volumes with Ground Truth  |   Performance Evaluation/nifti1_read_compare.c  | 
| Results Obtained  |    Results  | 
| Master's Defense Slides  |   S4410963_MasterThesis_Defense.pdf  | 
| Final Report  |   mCS_2023_PoojaG.pdf  | 

# Usage
1.  Compile script "nifti1_Segment_filter.c" to filter all the classified volumes(Training or Test volumes) at once.
```bash
gcc -o nifti1_segmentts nifti1_segmentimp.c nifti1_read_write.c -lpthread
```
2. Run the script as follows, here "input directory" directory which consists of respective nifti volumes and "output directory" consists of classified volumes by the trained classifier:
```bash
 ./nifti1_segmentts <input directory> <output directory> <nthreads>
 ```
1.  Compile script "nifti1_read_compare.c" to read all the ground truth volumes and the classified volumes stored in dir1_path and dir2_path (to be specified in the script)

```bash
gcc -o nifti_to_avss nifti1_read_compare.c avs_io.c -lpthread
```
2. Run
```bash
./nifti_to_avss
```
To run matlab files use matlab version of: R2021a (MATLAB 9.10) or any newer version of matlab will be supported.



