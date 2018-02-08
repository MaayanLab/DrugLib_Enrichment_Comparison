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

def file_exists(fname):
	'''Checks if a file exists in the directory, printing a statement if so.'''
	if os.path.isfile(fname):
		print(fname, 'has already been created.')
		return True
	else: return False

def convert_gmt(gmt_fname, output_type='gvm'):
	'''
	Converts a gmt file to either a dict or a gvm dataframe, and then returns it.
	If converting to a gvm, the resulting dataframe is also saved as a csv file in the current working directory.
	gmt_fname : str
		The gmt filename, i.e. "ChEA_2016.txt". The file must be in the current working directory. 
	output_type : str
		Is either 'gvm' or 'dict'.
	'''
	#The first segment of this list, ending at 'nan', comes from pandas.read_csv(na_values) documentation.
		#From the original list, 'NA' was removed because it is, in fact, a gene. 
	#The second segment of this list, beginning with '---', was added according to my own observations.
	MY_NA_VALS = {'', '#N/A', '#N/A N/A', '#NA', '-1.#IND', '-1.#QNAN',
		'-NaN', '-nan', '1.#IND', '1.#QNAN', 'N/A', 'NULL', 'NaN', 'nan', 
		'---', '[NULL]'}

	def gmt_to_gvm(reader, output_fname):
		'''
		Transforms a gmt .txt file to a gvm dataframe, saves it, and returns it.
		reader : reader object
		output_fname : str
			The output filename, i.e. "ChEA_2016_gvm.csv".
		'''
		lib_name = output_fname.partition('_gvm')[0]

		gvm = pd.DataFrame(dtype=bool)
		#Each row corresponds to an annotation (e.g. a tf or drug). So, for each annotation:
		for row in reader:
			annotation = row[0]
			print(annotation)
			if annotation not in MY_NA_VALS:
				#Obtain the gene set, removing delimiting strings and not accepting null values.
				genes = [str(x).replace(',1.0', '') for x in row[2:] if str(x).replace(',1.0', '') not in MY_NA_VALS]
				#Remove repeats. To save time, comment out this below line if gene sets are known to be unique.
				genes = set(genes)
				#Create a column gene vector for this annotation...
				vec = pd.DataFrame(True, index = genes, columns = [annotation], dtype=bool)
				#...and concatenate it with the growing dataframe. 
				gvm = pd.concat([gvm,vec], axis=1)
		gvm.index.name = lib_name
		#Save file to the current working directory, with values True and '' (instead of False, to save space).
		gvm.to_csv(output_fname, sep='\t')
		#Return the dataframe in sparse form, with values True and False.
		gvm = gvm.fillna(False).to_sparse()
		return gvm

	def gmt_to_dict(reader):
		'''
		Converts a gmt .txt file to a dict, and returns it.
		reader : reader object
		'''
		d = {}
		#Each row corresponds to an annotation (e.g. a tf or drug). So, for each annotation:
		for row in reader:
			annotation = row[0]
			if annotation not in MY_NA_VALS:
				#Obtain the gene set, removing delimiting strings and not accepting null values.
				#Store the row data in the dict with the annotation as the key and the gene set as the value.
				d[annotation] = {str(x).replace(',1.0', '') for x in row[2:] if x not in MY_NA_VALS}
		return d

	#If the gvm is requested, check to see if it has already been created. If so, simply load it and return it.
	output_fname = gmt_fname.partition('.')[0] + '_gvm.csv'
	if output_type == 'gvm' and os.path.isfile(output_fname): 
		print('will use old gvm file for', gmt_fname.partition('.')[0])
		result = open_gvm(output_fname)
		#Ensure the index name is the name of the gmt.
		result.index.name = gmt_fname.partition('.')[0]
		return result

	#Else, open the gmt file.
	print('getting', gmt_fname, 'as', output_type)
	if not os.path.isfile(gmt_fname): raise ValueError(gmt_fname, 'is not in the current working directory.')	
	with open(gmt_fname, 'r') as f:
		reader = csv.reader(f, delimiter = '\t')
		#Use either gmt_to_gvm or gmt_to_dict to create and return the desired data structure.
		if output_type == 'gvm': return gmt_to_gvm(reader, output_fname)
		elif output_type == 'dict': return gmt_to_dict(reader)
		else: raise ValueError(output_type, 'must be either gvm or dict.')

def combine_gmts(gmt_fnames, output_fname, merge_type='union'):
	'''
	Creates a gvm with merged gene sets from two gmt files with corresponding annotations in the same order.
	E.g. gmt file A has down-regulated genes, and gmt file B has up-regulated genes for the same annotations. 
	gmt_fnames : list-like
		Contains the names of the two gmt files, which must be in the current working directory. 
	output_fname : str
		Name of the output file.
	merge_type : str
		'union' or 'intersection'
	'''
	def dict_to_gvm(d, output_fname, ROWS_PER_CHUNK=2500):
		chunk_fname_prefix = output_fname.partition('_gvm.csv')[0] + '_gvm_chunk_'

		#Convert each key to a gene vector column in the matrix.
		df = pd.DataFrame(dtype=bool)
		it = 0
		subit = 0
		for k in d:
			print(it, subit, k)
			if not os.path.isfile(chunk_fname_prefix + str(it) + '.csv'):
				vec = pd.DataFrame(True, index = d[k], columns = [k], dtype=bool)
				df = pd.concat([df,vec], axis=1)
				df.index.name = output_fname.partition('_gvm.csv')[0]
			subit += 1
			#Switch to a new data.frame 'chunk' every `ROWS_PER_CHUNK` rows.
			if subit == ROWS_PER_CHUNK:
				if not df.empty: df.to_csv(chunk_fname_prefix + str(it) + '.csv', sep='\t')
				df = pd.DataFrame(dtype=bool)
				subit = 0
				it += 1
		if (subit != 0) and (not df.empty): df.to_csv(chunk_fname_prefix + str(it) + '.csv', sep='\t')

		#Merge all the "chunks" together.
		gvm = pd.DataFrame(dtype=bool)
		for fname in os.listdir(os.getcwd()):
			if chunk_fname_prefix in fname:
				print(fname)
				chunk = open_gvm(fname)
				gvm = pd.concat([gvm,chunk], axis=1)
		gvm = gvm.replace(to_replace=False, value='')
		gvm.to_csv(output_fname, sep='\t')

		print('Done. Delete the chunk files manually after verifying accuracy of merged gvm.')

		return df

	#Exit if already completed.
	if file_exists(output_fname): return
	print('creating', output_fname)

	#Otherwise, first convert each gmt to a dict.
	print('getting gmts as dicts')
	dicts = [convert_gmt(x, 'dict') for x in gmt_fnames]

	#Then, merge the dicts according to the specified method, union or intersection.
	print('merging dicts')
	combined_dict = {}
	if merge_type == 'union':
		for k in dicts[0]: combined_dict[k] = dicts[0][k] | dicts[1][k]
	elif merge_type == 'intersection': 
		for k in dicts[0]: combined_dict[k] = dicts[0][k] & dicts[1][k]

	#Finally, convert the merged dict to a gvm matrix.
	#(Also save it in the current working directory.)
	print('converting merged dict to df')
	return dict_to_gvm(combined_dict, output_fname)

os.chdir('original_drug-gene_libs')

combine_gmts(['Drug_Perturbations_from_GEO_down.txt', 'Drug_Perturbations_from_GEO_up.txt'], 
	'CREEDS_Drugs_gvm.csv', merge_type='union')
combine_gmts(('LINCS_L1000_Chem_Pert_up.txt','LINCS_L1000_Chem_Pert_down.txt'), 
	'LINCS_L1000_Chem_Pert_gvm.csv', merge_type='union')

os.chdir('..')