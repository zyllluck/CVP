# -*- coding: utf-8 -*-
"""
Created on Mon Mar  8 16:10:44 2021

@author: zyll1
"""
import numpy as np
import itertools, time
from scipy import stats
from sklearn.model_selection import KFold
from sklearn.linear_model import LinearRegression
import copy
import multiprocessing
from multiprocessing import Pool
from scipy.stats import pearsonr
import pandas as pd
import random
import pingouin
import codecs
import scipy
from scipy import sparse
import scipy.sparse
import os, sys

def text_save(filename, data):
	file = open(filename,'w')
	for i in range(len(data)):
		s = str(data[i]).replace('[','').replace(']','')
		s = s.replace("'",'').replace(',','') +'\n'
		file.write(s)
	file.close()


def ind(gene,gene_name):
	xx=[gene_name.index(i) for i in gene]
	return xx


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
usage: python CVP_inferring_causality.py -m=m -infile=infile_name -threshold=thres -outfile=outfile_name
Options and arguments:
-m: The parameter 'm' indicates the grouping for m-fold cross-validation. The default value is 3.
-threshold: The parameter 'threshold' is a correlation threshold produced by program 'pcc_threshold.py', the fault threshold value is 0.
-infile: The parameter 'infile' indicates the file name for input. The input file is a data matrix that each row represents a sample and each column represents a variable. The file 'data_example.txt' is an example file for input file.
-outfile: The parameter 'outfile' indicates the name of output file to store the causality with format '.txt'. The default value is 'CVP_result.txt'. The output file has three columns, the first column indicates the cause variable, the second column indicates the effect variable and the third column indicates the causal strength.
"""


if "-help" in param.keys() or "-h" in param.keys():
    print(help_msg)

if "-infile" not in param.keys():
    print("Parameter missing!")
    print(help_msg)
    exit()
else:
    path_name=[param["-infile"]]


if "-threshold" not in param.keys():
    print("Parameter missing!")
    print(help_msg)
    exit()

if "-m" not in param.keys():
    print("Parameter missing!")
    print(help_msg)
    exit()

if "-outfile" not in param.keys():    
    fold="CVP_result.txt"
else:
    fold=param["-outfile"]

if "-threshold" in param.keys():
    threshold=float(param["-threshold"])
    if threshold != 0:
        print("Please set the correct threshold of Correlation threshold")
        print(help_msg)


if "-m" in param.keys():
    m=int(float(param["-m"]))
    if m != 3:
        print("Please set the correct group for cross-validation")
        print(help_msg)





path0='./CVP-main/'
path0='.'+os.sep
path=path0+path_name[0]
data_dic,gene_name=lexicon1(path)
data_example_pd=pd.read_csv(path,sep='\s+')
data_0=np.array(list(data_dic.values()),dtype=float).T
pcc0_matrix=np.corrcoef(data_0.T)
pcc_0=pcc0_matrix
ptcc_tho=0.1


def pcc0_screen(gene1,pcc_0,threshold):
	gene1_index=ind([gene1],gene_name)
	pcc0_gene1=pcc_0[gene1_index[0]]	
	pcc0_pair=[{gene1:gene_name[i],'pcc0':pcc0_gene1[i]} for i in range(len(gene_name))]
	pcc0_pair=sorted(pcc0_pair,key=lambda x:abs(x['pcc0']),reverse=False)
	pcc0_pairs=[item for item in pcc0_pair if not item[gene1]==gene1]
	pcc0_pairs=[item for item in pcc0_pairs if abs(item['pcc0'])>=threshold]
	#pcc_gene1=[list(pcc0_pairs[i].values()) for i in range(len(pcc0_pairs))]
	pcc_gene1=[list(pcc0_pairs[i].values()) for i in range(len(pcc0_pairs))]
	co_gene=[lis[0] for lis in pcc_gene1]
	return co_gene
	
def ttr(gene1,gene_1,gn,m):
	if len(gene_1)>=1:
		erro=[]
		k=0
		kf = KFold(n_splits=m,shuffle=True,random_state=None)
		y_index=ind([gene1],gn)
		x_index=ind(gene_1,gn)
		for train_index, test_index in kf.split(data_0):
			k=k+1
			train = data_0[list(train_index)]
			test= data_0[list(test_index)]
			reg= LinearRegression().fit(train[:,x_index],train[:,y_index])#通过训练集模拟回归模型
			erro_1= np.sum((reg.predict(test[:,x_index])-test[:,y_index])**2)#计算测试集的残差平方和
			erro.append(erro_1)
		erro_avg=np.mean(erro)
	else:
		erro_avg=0		
	return erro_avg

def regulon(gene1):
	gene2=pcc0_screen(gene1,pcc_0,threshold) 
	end=[]
	end_1=[gene1]
	end_2=[0]
	if len(gene2)>=2:
		erro_1=ttr(gene1,gene2,gene_name,m)
		##########
		gene_b=copy.deepcopy(gene2)
		gene2_len=list(range(len(gene2)))
		for j in gene2_len:
			gene_b.remove(gene2[j])
			erro_2=ttr(gene1,gene_b,gene_name,m)
			erro_j=erro_1-erro_2
			if erro_j<0:
				gene_b.append(gene2[j])
				end_1.append(gene2[j])
				end_2.append(erro_j)
			else:
				erro_1=erro_2
		end.append([end_1,end_2])
	else:
		end.append([end_1,end_2])
	return end

#A single process can run directly
if __name__=="__main__":
	result=[]
	for i in gene_name:
		regu=regulon(i)
		result.append(regu)
	text_save(path0+path_name[0][0:-4]+'_pcc0_r.txt',result)
'''
#Multi-process can be run directly
if __name__=="__main__":
	pool = Pool()
	result =pool.map(regulon,gene_name)
	text_save(path0+path_name[0][0:-4]+'_pcc0_r.txt',result)
'''

f = codecs.open(path0+path_name[0][0:-4]+'_pcc0_r.txt', mode='r', encoding='utf-8')
line = f.readline() 
targets = []
regulators=[]
while line:
  a = line.split()
  b = a[0:1] 
  b1=a[0:int(len(a)/2)]
  targets.append(b[0])
  regulators.append(b1)
  line = f.readline()
f.close()


def gene_list(gene_all,gene):
	gene_pcc=gene_all[1:]
	gene_pcc.remove(gene)
	return gene_pcc

def partial_cc(gene_all,ptcc_tho):
	partial=[data_example_pd[gene_all].partial_corr(x=gene_all[0], y=gene_all[i], covar=gene_list(gene_all,gene_all[i]),
			 method='pearson').round(3) for i in range(1,len(gene_all))]
	partial_gene=[[gene_all[0],gene_all[i],list(partial[i-1]['p-val'])[0]] for i in range(1,len(gene_all))]
	partial_gene1= [l[0:3] for l in partial_gene if l[2] <=ptcc_tho]
	partial_gene1.sort(key= lambda x:x[2],reverse=True)
	co_pcc_gene=[gene_all[0]]+[x[1] for x in partial_gene1]
	return co_pcc_gene

def regulon_partial(gene_all):
	gene2=partial_cc(gene_all,ptcc_tho) 
	gene2=gene2[1:]
	end=[]
	end_1=[gene_all[0]]
	end_2=[0]
	if len(gene2)>=2:
		erro_1=ttr(gene_all[0],gene2,gene_name,m)
		gene_b=copy.deepcopy(gene2)
		gene2_len=list(range(len(gene2)))
		for j in gene2_len:
			gene_b.remove(gene2[j])
			erro_2=ttr(gene_all[0],gene_b,gene_name,m)
			erro_j=erro_1-erro_2
			if erro_j<0:
				gene_b.append(gene2[j])
				end_1.append(gene2[j])
				end_2.append(erro_j)
			else:
				erro_1=erro_2
		end.append([end_1,end_2])
	else:
		end.append([end_1,end_2])
	return end


#Single process can be run directly
if __name__=="__main__":
	result1=[]
	for i in regulators:
		#print(i)
		regu=regulon_partial(i)
		result1.append(regu)		   
	text_save(path0+path_name[0][0:-4]+'_pcc1_r.txt',result1)
'''
#Multi-process can be run directly
if __name__=="__main__":
	pool = Pool()
	result1 =pool.map(regulon_partial,regulators)
	text_save(path0+path_name[0][0:-4]+'_pcc1_r.txt',result1)
'''


filename=path0+path_name[0][0:-4]+'_pcc1_r.txt'
def  pairs1(filename):
	data=montage(filename)
	gstr=[data[i][0:int(len(data[i])/2)] for i in range(len(data))]
	rdata=[data[i][int(len(data[i])/2)+1:len(data[i])] for i in range(len(data))]
	rdata= [list(map(float, rdata[i])) for i in range(len(rdata))]
	z=np.concatenate((rdata), axis=0)
	z_tho=0
	g=[]
	g0=[]
	for i in range(len(gstr)):
		gstr0=[gstr[i][0]]
		for j in range(len(rdata[i])):
			if abs(rdata[i][j])>z_tho:
				gstr0.append(gstr[i][j+1])
				gname=[gstr[i][j+1],gstr[i][0],-rdata[i][j]]
				g0.append(gname)
		if len(gstr0)>1:
			g.append(gstr0)
	return g0
predict_pairs1=pairs1(path0+path_name[0][0:-4]+'_pcc1_r.txt')
predict_pairs1.sort(key= lambda x:abs(x[2]),reverse=True)
text_save(path0+fold,predict_pairs1)

os.remove(path0+path_name[0][0:-4]+'_pcc0_r.txt')
os.remove(path0+path_name[0][0:-4]+'_pcc1_r.txt')




