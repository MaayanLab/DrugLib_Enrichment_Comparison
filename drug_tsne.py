import os
import numpy as np
import pandas as pd
import scipy.sparse
import matplotlib.pyplot as plt
from sklearn.manifold import TSNE
from sklearn.decomposition import TruncatedSVD
from scipy.spatial.distance import euclidean
import time



def open_gvm(fname, already_have_LINCS=True):
	print('opening', fname)
	#Read in chunks if the file is too big.
	if 'interactions' in fname or 'LINCS' in fname:
		file_chunks = pd.read_csv(fname, keep_default_na = False, sep='\t', 
			low_memory=False, encoding='Latin-1', index_col=0, chunksize=1000)
		gvm = pd.concat(file_chunks)
		gvm = gvm.replace(to_replace='',value=False).astype(bool)
	else:
		gvm = pd.read_csv(fname, keep_default_na = False, na_values=('',), sep='\t', 
			low_memory=False, encoding='Latin-1', index_col=0)
		##Workaroud from bug which raised error if keep_default_na=False and index_col=0 are both used.
		#gvm = pd.read_csv(fname, keep_default_na = False, na_values=('',), sep='\t', low_memory=False, encoding='Latin-1')
		#lib_name = fname.partition('_gvm.csv')[0]
		#gvm.set_index(gvm[lib_name], inplace=True)
		#gvm.drop([lib_name], axis=1, inplace=True)

		#Convert blank cells to False, and ensure bool type. 
		gvm = gvm.fillna(False).astype(bool)
	return gvm

GET_DTARGET_LIBNAME = {
	'1_DrugBank_Column_5OrMoreTargets_MissingChEMBL_IDsRemoved_1-14-18':'DrugBank',
	'2_TargetCentral_Column_5OrMoreTargets_MissingChEMBL_IDsRemoved_1-14-18':'TargetCentral', 
	'3_RepurposeHub_Column_5OrMoreTargets_MissingChEMBL_IDsRemoved_1-14-18':'RepurposeHub', 
	'4_DGIdb_Column_5OrMoreTargets_MissingChEMBL_IDsRemoved_1-14-18':'DGIdb', 
	'5b_DrugCentral_Column_5OrMoreTargets_MissingChEMBL_IDsRemoved_Human_1-14-18':'DrugCentral',
	'3_EdgeLists_Union_10-05-17':'DBTC Union',
	'4_EdgeLists_Intersection_10-05-17':'DBTC Intersect'
	}

dtarget_libs = (
	'1_DrugBank_Column_5OrMoreTargets_MissingChEMBL_IDsRemoved_1-14-18',
	'2_TargetCentral_Column_5OrMoreTargets_MissingChEMBL_IDsRemoved_1-14-18', 
	'3_RepurposeHub_Column_5OrMoreTargets_MissingChEMBL_IDsRemoved_1-14-18', 
	'4_DGIdb_Column_5OrMoreTargets_MissingChEMBL_IDsRemoved_1-14-18', 
	'5b_DrugCentral_Column_5OrMoreTargets_MissingChEMBL_IDsRemoved_Human_1-14-18',
	)
ppi_libs = (
	'hu.MAP',
	'BioGRID',
	'ARCHS4'
	)

dtarget_libnames = ['expanded_drug-gene_libs\\' + dtarget + '_expanded_with_' + ppi + '_gvm2.csv' 
	for dtarget in dtarget_libs for ppi in ppi_libs]

x = dtarget_libnames[0]

x = open_gvm(x)

for i in x.columns.values:
	print(i, end='')
	time.sleep(.003)
