import csv
import os
import pickle
import numpy as np
import pandas as pd
from joblib import Parallel, delayed
import enrichment_functions as m
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier, \
	RandomTreesEmbedding, AdaBoostClassifier, ExtraTreesClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.svm import LinearSVC
from xgboost import XGBClassifier
from random import uniform as rand
import scipy.stats as stats
#from get_classifiers import get_classifiers #REMOVED. the script was put in the old or unused scripts folder.
import h5py

def open_gvm(fname):
	print('opening', fname)
	#Read in chunks if the file is too big.
	if 'interactions' in fname:
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

def clean(annot):
	'''Extracts drug synonym from annotation.'''
	return str(annot).partition('|||')[0]

def get_overlapping_ilib_annots(ilib_name, ilib_annots, slib_name, slib_annots):
	'''Return the annotations in the input library which have matches in the search library.'''
	cleaned_ilib_annots = {clean(annot) for annot in ilib_annots}
	cleaned_slib_annots = {clean(annot) for annot in slib_annots}
	cleaned_overlaps = cleaned_ilib_annots & cleaned_slib_annots
	overlapping_ilib_annots = [annot for annot in ilib_annots if clean(annot) in cleaned_overlaps]
	return overlapping_ilib_annots

def get_enrichment_algorithms(ilib_gvm, slib_gvm, ilib_name, slib_name):
	'''
	Returns a dataframe with methods and their parameters.
	ilib_gvm : pandas.DataFrame
		the input library's gene vector matrix, from which column vectors are being used as input gene sets
	slib_gvm: pandas.DataFrame
		the search library's gvm, whose annotations are being ranked by the enrichment algorithm based on their gene sets
	ilib_name : str
		the name of ilib_gvm, for example "ChEA_2016"
	slib_name: str
		the name of slib_gvm, for example "CREEDS"
	'''
	#=================================================
	#Define any other variables, as necessary, here.
	#=================================================
	train_group = slib_gvm
	features = slib_gvm.columns.values
	#REMOVED-DOES NOT WORK. see get_classifiers.py in the old or unused scripts folder.
	##Create a classifier - only necessary for ML_fisher_features (REMOVED).
	##Otherwise, comment out. 
	#classifier = get_classifiers(ilib_name, slib_name)

	#======================================================================================================================
	#This is where you specify which enrichment methods to run.
	#Specify each method as a column in a dataframe, where the first row is the enrichment function,
	#	the second row is a tuple of parameters, and the column name is the chosen name of the enrichment method.
	#For example, for a method using the RandomForestClassifier with max_depth=3, you might say:
	#	df['RF_md3'] = (m.ML_wrapper, (RandomForestClassifier, train_group, features, 10304, max_depth=3))
	#See enrichment_functions.py for the available methods and necessary params.
	#IMPORTANT: DO NOT specify the input_geneset parameter:
	#	this is created and called later, in perform_enrichment().
	#======================================================================================================================
	enrichment_algorithms = pd.DataFrame(index=['func', 'params'])
	enrichment_algorithms['Fisher'] = [m.Fisher, [slib_gvm]] 
	#enrichment_algorithms['RandomForest'] = [m.ML_wrapper, [RandomForestClassifier, train_group, features, 101317]]
	return enrichment_algorithms
	#======================================================================================================================

def perform_enrichment(pair):
	'''This function is called for each lib pair, and iterates over each method and each tf. 
	pair : tuple of str
		(input library fname, search library fname)
	'''

	#Get the algorithms with which to perform enrichment. 
	algorithm_name = 'Fisher'

	ilib_name, slib_name = [lib.partition('\\')[2].partition('_gvm')[0] for lib in pair]
	prefix = 'input_' + ilib_name + '_into_' + slib_name
	output_fnames = ('results\\' + prefix + '_' + algorithm_name + '.csv',)

	#TEMPORARY
	if 'interactions' in ilib_name and 'interactions' in slib_name:
		print('skipping')
		return

	if os.path.isfile(output_fnames[0]): 
		print('score file already created for', algorithm_name, prefix)
		return

	ilib, slib = (open_gvm(lib) for lib in pair)
	func = m.Fisher
	params = [slib,]
	print('Beginning enrichment analysis inputting', ilib_name, 'into', slib_name)

	#Get the input library annotations whose corresponding tf/drug also corresponds to 
	#	at least one search library annotation.
	overlapping_ilib_annots = get_overlapping_ilib_annots(ilib_name, ilib.columns.values, 
		slib_name, slib.columns.values)
	print(str(len(overlapping_ilib_annots)), 'overlaps')

	#Results will be stored after each tf iteration.
	score_dfs = [pd.DataFrame() for n in range(len(output_fnames))]
	#Iterate over each overlapping ilib annotation.
	for annot in overlapping_ilib_annots:
		print(algorithm_name, annot) #for diagnostics.
		input_geneset = ilib.index[ilib[annot]]
		#Get scores for all the slib annotations.
		result = func(input_geneset, *params)
		for x in range(len(score_dfs)): 
			df = score_dfs[x]
			#Store this result as a column in the score df.
			if len(score_dfs) == 1: df[annot] = result
			else: df[annot] = result[x]
	for x in range(len(score_dfs)): 
		#Save the score_dfs as csv files.
		df = score_dfs[x]
		df.index = slib.columns
		df.to_csv(output_fnames[x], sep='\t')
	return

if __name__ == '__main__':
	
	expanded_drug_libs = ['repurposing_drugs_20170327','interactions', 
		'1_DrugBank_Edgelist_10-05-17', '2_TargetCentral_Edgelist_10-05-17',
		'3_Edgelists_Union_10-05-17', '4_EdgeLists_Intersection_10-05-17']
	expansion_ppi_libs = ['hu.MAP','BioGRID','ARCHS4']
	expanded_gvms = ['expanded_drug-gene_libs\\' + drug_lib + '_expanded_with_' + ppi_lib + 
		'_gvm.csv' for drug_lib in expanded_drug_libs for ppi_lib in expansion_ppi_libs]

	non_expanded_drug_libs = ['CREEDS_Drugs','DrugMatrix_Union']
	non_expanded_gvms = ['original_drug-gene_libs\\' + drug_lib + 
		'_gvm2.csv' for drug_lib in non_expanded_drug_libs]

	#======================================================
	#Choose the libraries with which to perform enrichment.
	#======================================================
	libs = expanded_gvms + non_expanded_gvms
	#======================================================

	if not os.path.isdir('results'): os.makedirs('results')

	#Iterate over each gmt pair.
	lib_df_pairs = [(a,b) for a in libs for b in libs if a != b]
	Parallel(n_jobs=1, verbose=0)(delayed(perform_enrichment)(pair) for pair in lib_df_pairs)