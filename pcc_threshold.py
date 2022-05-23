import numpy as np
from scipy.stats import pearsonr
import pandas as pd
import os, sys

def text_save(filename, data):
	file = open(filename,'w')
	for i in range(len(data)):
		s = str(data[i]).replace('[','').replace(']','')
		s = s.replace("'",'').replace(',','') +'\n'
		file.write(s)
	file.close()

def text(file_name):
	with open (file_name,'r+') as f:
		line = f.readline()
		while line:
			yield line.split()
			line = f.readline()

def montage(file_name):
	a=[]
	for i in text(file_name):
		a.append(i)
	return a

def  lexicon1(file_name):
	data=montage(file_name)
	data_dic={}
	for i in range(len(data[0])):
		data_dic[data[0][i]]=[data[j][i] for j in range(1,len(data))]
	gene_name=list(data_dic.keys())
	return data_dic,gene_name


param={}
for i in range(1,len(sys.argv)):
	t=sys.argv[i].split("=")
	param[t[0].lower()]=t[1]
	
help_msg="""
usage: python pcc.py -infile = data_example.txt
Options and arguments:
-infile: The parameter 'infile' indicates the file name for input. The input file is a data matrix that each row represents a sample and each column represents a variable. The file 'data_example.txt' is an example file for input file.
"""

if "-help" in param.keys() or "-h" in param.keys():
	print(help_msg)
	
if "-infile" not in param.keys():
	print("Parameter missing!")
	print(help_msg)
	exit()
else:
	path_name=[param["-infile"]]

#---changed----
path0='./CVP-main/'
path0='.'+os.sep
#path_name=['data_example.txt']
path=path0+path_name[0]
data_dic,gene_name=lexicon1(path)
data_df=pd.read_csv(path,sep='\t')
ps=data_df.corr('pearson')
ps = ps.mask(np.triu(np.ones(ps.shape, dtype=np.bool_)))
ps=np.array(ps)
ps=abs(ps.flatten())
ps=-np.sort(-ps)
ps=np.delete(ps,np.where(np.isnan(ps))[0],axis=0)

if "-outfile" not in param.keys():	
	fold="threshold_result.txt"
else:
	fold=param["-outfile"]

if len(gene_name)>1000:   
	#gene_tho=[100,ps[int(50*len(gene_name))]]
	gene_tho=[ps[int(50*len(gene_name))]]
	text_save(path0+fold,gene_tho)
	print("The correlation threshold for CVP algorithm is:",gene_tho[0])

else:
    print("The correlation threshold for CVP algorithm is: 0")
