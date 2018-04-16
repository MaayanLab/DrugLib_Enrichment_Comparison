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

def list_of_drug_libs_fnames():
	'''Generates and returns a tuple of drug library file names.'''
	expanded_drug_libs = (
		'stitch_500cutoff_noENSP',
		'stitch_600cutoff_noENSP',
		'stitch_700cutoff_noENSP',
		'stitch_800cutoff_noENSP',
		'drugtargetcommons_5ormoretargets')
	ppi_libs = (
		'huMAP',
		'BioGRID',
		'ARCHS4')

	non_expanded_drug_libs = (
		'CREEDS',
		'LINCS_L1000_Chem_Pert')

	drug_libs = tuple(
		'expanded_drug-gene_libs\\' + edl + '_expanded_with_' + ppi + '_gvm.csv' 
		for edl in expanded_drug_libs for ppi in ppi_libs) + tuple(
		'original_drug-gene_libs\\' + drug_lib + '_pcids_gvm.csv' 
		for drug_lib in non_expanded_drug_libs)

	drug_libs2 = tuple(
		'original_drug-gene_libs\\' + edl + '_gvm.csv' 
		for edl in expanded_drug_libs) + tuple(
		'original_drug-gene_libs\\' + drug_lib + '_pcids_gvm.csv' 
		for drug_lib in non_expanded_drug_libs)

	return drug_libs2

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

def clean(annot):
	'''Extracts drug synonym from annotation.'''
	#This handles duplicate column names.
	return str(annot).partition('.')[0]

def get_common_ilib_annots(ilib_name, ilib_annots, slib_name, slib_annots):
	'''Return the annotations in the input library which have matches in the search library.'''
	cleaned_ilib_annots = {clean(annot) for annot in ilib_annots}
	cleaned_slib_annots = {clean(annot) for annot in slib_annots}
	cleaned_overlaps = cleaned_ilib_annots & cleaned_slib_annots
	common_ilib_annots = [annot for annot in ilib_annots if clean(annot) in cleaned_overlaps]
	return common_ilib_annots

def get_enrichment_algorithms(ilib_gvm, ilib_name, slib_gvm, slib_name):
	'''
	Returns a dataframe with methods and their parameters.
	ilib_gvm : pandas.DataFrame
		the input library's gene vector matrix, from which column vectors are being used as input gene sets
	slib_gvm: pandas.DataFrame
		the search library's gvm, whose annotations are being ranked by the enrichment algorithm based on their gene sets
	ilib_name : str
		the name of ilib_gvm, for example "ChEA_2016"
	slib_name : str
		the name of slib_gvm, for example "CREEDS_TFs"
	algorithm_name_only : bool
		If True, returns only the names of the algorithms. 
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
	#	this is created and called later, in enrichment().
	#======================================================================================================================
	enrichment_algorithms = pd.DataFrame(index=['func', 'params'])
	enrichment_algorithms['Fisher'] = [m.Fisher, [slib_gvm]] 
	#enrichment_algorithms['RandomForest'] = [m.ML_wrapper, [RandomForestClassifier, train_group, features, 101317]]
	#======================================================================================================================
	return enrichment_algorithms

def enrichment(pair):
	'''This function is called for each lib pair, and iterates over each method and each tf. 
	pair : tuple of str
		(input library fname, search library fname)
	'''
	ilib_fname, slib_fname = pair
	ilib_name, slib_name = (fname.rpartition('\\')[2].partition('_gvm')[0] for fname in pair)
	prefix = 'input_' + ilib_name + '_into_' + slib_name

	script_dir = os.path.dirname(os.path.abspath(__file__))
	results_dir = os.path.join(script_dir, 'results\\')

	#==============================================================================================================
	#OPTIONAL: To speedily allow this script to continue working from where it left off when it was previously run,
	#	specify ALL the enrichment algorithm names again here.
	#==============================================================================================================
	algorithm_names = ('Fisher',) #e.g. `('Fisher', 'RandomForest')` or `None`. 
	#==============================================================================================================

	#TEMP
	if 'interactions' in ilib_name and 'interactions' in slib_name:
		print('skipping', prefix)
		return

	#Exit if enrichment between this library pair has already been done.
	#See above chunk: `algorithm_names` must be specified.
	if algorithm_names is not None:
		output_fnames = tuple(results_dir + prefix + '_' + name + '.csv' for name in algorithm_names)
		if all((os.path.isfile(fname) for fname in output_fnames)):
			print('Already done enrichment for', prefix)
			return

	#Otherwise, begin by loading the gvms.
	print('Beginning enrichment analysis inputting', ilib_name, 'into', slib_name)
	ilib_gvm, slib_gvm = (open_gvm(fname) for fname in pair)

	#Get the algorithms with which to perform enrichment. 
	enrichment_algorithms = get_enrichment_algorithms(ilib_gvm, ilib_name, slib_gvm, slib_name)

	#Get the input library annotations whose corresponding tf/drug also corresponds to 
	#	at least one search library annotation.
	common_ilib_annots = get_common_ilib_annots(
		ilib_name, ilib_gvm.columns.values, 
		slib_name, slib_gvm.columns.values)
	print(str(len(common_ilib_annots)), 'overlaps')

	ilib_gvm = ilib_gvm.loc[:,common_ilib_annots]

	#Iterate over each algorithm (i.e. each column).
	for algorithm_name in enrichment_algorithms:
		algorithm = enrichment_algorithms[algorithm_name]

		#Some methods actually return multiple results. These will need multiple score files.
		if algorithm_name == 'ZAndCombined': output_fnames = (
			results_dir + prefix + '_Z.csv', 
			results_dir + prefix + '_Combined.csv')
		#========================================================================================
		#If using a algorithm which returns multiple results, add the appropriate elif statement here.
		#========================================================================================
		#elif algorithm_name == 'Pairwise_Gini': output_fnames = (prefix + '_Pair_Gini_ltf100_1.csv', 
		#	prefix + '_Pair_Gini_ltf100_5.csv', prefix + '_Pair_Gini_ltf100_10.csv', prefix + '_Pair_Gini_ltf100_25.csv')
		#========================================================================================
		else: output_fnames = (results_dir + prefix + '_' + algorithm_name + '.csv',)

		#Exit if this algorithm has already been run for this library pair.
		if os.path.isfile(output_fnames[0]): 
			print('score file already created for', algorithm_name)
			return

		#Otherwise, run enrichment.
		#Results will be stored after iteration over the common ilib annotations.
		score_dfs = [pd.DataFrame() for n in range(len(output_fnames))]
		for annot in common_ilib_annots:
			print(algorithm_name, annot) #for checking progress.
			input_geneset = ilib_gvm.index[ilib_gvm[annot]]
			#Get scores for all the slib annotations.
			result = algorithm['func'](input_geneset, *algorithm['params'])
			for x in range(len(score_dfs)): 
				df = score_dfs[x]
				#Store this result as a column in the score df.
				if len(score_dfs) == 1: df[annot] = result
				else: df[annot] = result[x]
		for x in range(len(score_dfs)): 
			#Save the score_dfs as csv files.
			df = score_dfs[x]
			df.index = slib_gvm.columns
			df.to_csv(output_fnames[x], sep='\t')
	return

if __name__ == '__main__':

	libs = list_of_drug_libs_fnames()
	target_libs = libs[:-2]
	perturb_libs = libs[-2:]

	#fwd_pairs = [(a,b) for a in target_libs for b in perturb_libs]
	#bck_pairs = [(a,b) for a in perturb_libs for b in target_libs]

	fwd_pairs = [(a,b) for a in target_libs for b in perturb_libs]
	bck_pairs = [(b,a) for a in target_libs for b in perturb_libs]

	if not os.path.isdir('results'): os.makedirs('results')

	#Iterate over each gmt pair.
	lib_df_pairs = fwd_pairs + bck_pairs
	#Parallel(n_jobs=1, verbose=0)(delayed(enrichment)(pair) for pair in lib_df_pairs)
	for pair in reversed(lib_df_pairs): enrichment(pair)