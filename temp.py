import requests
import csv
import os
import h5py
import pickle
import numpy as np
import pandas as pd
from joblib import Parallel, delayed 

def open_gvm(fname):
	#Open the file.
	gvm = pd.read_csv(fname, keep_default_na = False, na_values=('',), sep='\t', low_memory=False, encoding='Latin-1', index_col=0)
	##Workaroud from bug which raised error if keep_default_na=False and index_col=0 are both used.
	#gvm = pd.read_csv(fname, keep_default_na = False, na_values=('',), sep='\t', low_memory=False, encoding='Latin-1')
	#lib_name = fname.partition('_gvm.csv')[0]
	#gvm.set_index(gvm[lib_name], inplace=True)
	#gvm.drop([lib_name], axis=1, inplace=True)

	#Convert blank cells to False, and ensure bool type. 
	gvm = gvm.fillna(False).astype(bool)
	return gvm

os.chdir('original_drug-gene_libs')

#for filename in os.listdir(os.getcwd()):
#	os.rename(filename, filename.replace('_transformed', '_gvm'))

dmu = open_gvm('DrugMatrix_Union_gvm.csv')
dmu = dmu.replace(to_replace=False, value='')
dmu.index.name = 'DrugMatrix_Union'
dmu.to_csv('DrugMatrix_Union_gvm.csv',sep='\t')

os.chdir('..')