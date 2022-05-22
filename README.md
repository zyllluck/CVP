#CVP
This is source code for inferring causality by the CVP algorithm in python. The CVP is an algorithm for causal inference, the full name of which is "Causality inference based on cross-validation prediction for non-time series data".

The 'CVP_inferring_causality.py' can be used to build the causality network by the CVP algorithm, and details of the parameters are shown below.

##Dependencies

python 3.7 or later
numpy 1.2 or later
pingouin 0.3.10 or later
scipy 1.6.1 or later
pandas 1.0.3 or later
scikit-learn 1.0.2 or later

##Usage
### Run the main program
python CVP_inferring_causality.py -m=m -infile=infile_name -threshold=thres  -outfile=outfile_name

-m: The parameter 'm' indicates that the group for cross-validation. The default value is 3.

-threshold: The parameter 'threshold' is a correlation threshold produced by program 'pcc_threshold.py', the fault threshold value is 0.

-infile: The parameter 'infile' indicates the file name for input. The input file is a data matrix that each row represents a sample and each column represents a variable. The file 'data_example.txt' is a example file for input file.

-outfile: The parameter 'outfile' indicates the name of output file to store the causality with format is '.txt'. The default value is 'CVP_result.txt'. The outputfile has three columns, the first column indicates the dependent variable, the second column indicates the independent variable and the third column indicates the causal intensity.

###Calculate the correlation threshold 

python pcc_threshold.py -infile=data_pcc_example.txt -gn=number -outfile=outfile_name

-infile: The parameter 'infile' indicates the file name for input. The input file is a data matrix that each row represents a sample and each column represents a variable. The file 'data_pcc_example.txt' is a example file for input file.

-gn: the number of variables

-outfile: The parameter 'outfile' indicates the name of output file to store threshold parameter for program 'CVP_inferring_causality.py' with format is '.txt'. The default value is 'threshold_result.txt'. The outfile contains two lines, the first line is 100, indicating the average degree, and the second line corresponds to the correlation threshold.




##Example 1: Using the CVP program to calculate the causality between variables

This file 'data_example.txt' is an example to infer causality for 9 variables with 8 samples. Each column in this file represents a variable, and each row in this file represents a sample.
We can infer the causality of the 9 variables by executing below:

python CVP_inferring_causality.py -infile=data_example.txt -threshold=0 -m=3 -outfile=outfile_name.txt


This 'CVP_inferring_causality.py' is used directly as long as the number of variables is less than 1000, otherwise the corresponding correlation threshold needs to be calculated by the 'pcc_threshold.py'.

Using the 'CVP_inferring_causality.py', a causal network of nine variables is constructed and the output is 'CVP_result.txt'. The output file has three columns, the first column indicates the dependent variable, the second column indicates the independent variable and the third column indicates the causal intensity.
Each row is a causality between two variables in the 9 test variables in 'data_example.txt'.

##Example 2: Calculate the correlation threshold when the number of variables is greater than 1000.

python pcc_threshold.py -infile=data_pcc_example.txt -gn=1137 -outfile=outfile.txt

This file 'data_pcc_example.txt' is matrix data, with columns representing variables and rows representing samples. The parameter 'gn' is the number of variables.

Using the 'pcc_threshold.py',the outputfile contains two lines, the first line is 100, indicating the average degree, and the second line corresponds to the correlation threshold. The default output is the average degree of 100 and the corresponding correlation threshold. Correlation thresholds can be calculated for different average degrees, depending on the amount of data.
