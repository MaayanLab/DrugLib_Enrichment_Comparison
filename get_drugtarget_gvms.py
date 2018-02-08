import requests
import csv
import os
import h5py
import pickle
import numpy as np
import pandas as pd
from joblib import Parallel, delayed 

def get_interactionlist(fname):
	if '1-14-18' in fname:
		interactions = pd.read_csv(fname, sep=',', encoding='latin1')
		interactions = interactions[['TargetGeneSymbol_Entrez','DrugName_Final']]

	elif '_10-05-17.csv' in fname: 
		interactions = pd.read_csv(fname, sep=',', encoding='latin1')
		interactions = interactions[['TargetGeneSymbol_Entrez','DrugName']]

	elif fname == 'interactions.tsv':
		interactions = pd.read_csv(fname, sep='\t', encoding='latin1')
		interactions.loc[interactions['gene_name'].isnull().values, 'gene_name'] = interactions.loc[
			interactions['gene_name'].isnull().values, 'gene_claim_name']
		interactions = interactions[['gene_name','drug_claim_primary_name']]

	elif fname == 'repurposing_drugs_20170327.txt':
		f = pd.read_csv(fname,sep='\t',skiprows=9,encoding='latin1')
		f = f.loc[~f['target'].isnull().values,]
		interactions = []
		for row in f.index:
			geneset = f.at[row,'target'].split('|')
			interactions = interactions + [np.column_stack([geneset,[f.at[row,'pert_iname']] * len(geneset)])]
		interactions = pd.DataFrame(np.vstack(interactions))

	else: raise ValueError(fname, ' was not recognized.')

	interactions.columns = ('gene','annotation')
	interactions = interactions.loc[~interactions['gene'].isnull().values,]
	interactions = interactions.loc[~interactions['annotation'].isnull().values,]
	interactions = interactions.drop_duplicates()
	return interactions

def interactionlist_to_gvm(interactionlist_fname):
	#Check if this has already been done. 
	output_fname = interactionlist_fname.partition('.')[0]+ '_gvm.csv'
	if os.path.isfile(output_fname): 
		#print(output_fname, 'already created.')
		return

	#Otherwise, proceed.
	#print(interactionlist_fname)
	interactionlist = get_interactionlist(interactionlist_fname)

	#The first segment of this list, ending at 'nan', comes from pandas.read_csv(na_values) documentation.
		#From the original list, 'NA' was removed because it is, in fact, a gene. 
	#The second segment of this list, beginning with '---', was added according to my own observations.
	MY_NA_VALS = {'', '#N/A', '#N/A N/A', '#NA', '-1.#IND', '-1.#QNAN',
		'-NaN', '-nan', '1.#IND', '1.#QNAN', 'N/A', 'NULL', 'NaN', 'nan', 
		'---', '[NULL]'}

	#Initialize the gvm with annotations as indices.
	gvm = pd.DataFrame(index=set(interactionlist['annotation']), dtype=bool)
	#For each gene:
	##(I think this is faster than iterating over annotations because there are less unique genes than annots.)
	genes = set(interactionlist['gene'])
	for gene in genes:
		#print(gene)
		#Get this gene's annotation vector: the collection of annotations with this gene in their set.
		annotations = interactionlist.loc[interactionlist['gene'] == gene, 'annotation'].values
		#annotations = set(annotations)
		#Add this annotation vector to the gvm.
		vec = pd.DataFrame(True, index=annotations, columns=[gene], dtype=bool)
		gvm = pd.concat([gvm,vec], axis=1)
	#Transpose such that genes are indices and annotations are columns.
	gvm = gvm.transpose()

	#Save the results.
	gvm = gvm.replace(to_replace=False, value='')
	gvm.to_csv(output_fname, sep='\t')
	return

os.chdir('original_drug-gene_libs')

interactionlists = (
	'1_DrugBank_Column_5OrMoreTargets_MissingChEMBL_IDsRemoved_1-14-18.csv',
	'2_TargetCentral_Column_5OrMoreTargets_MissingChEMBL_IDsRemoved_1-14-18.csv',
	'3_RepurposeHub_Column_5OrMoreTargets_MissingChEMBL_IDsRemoved_1-14-18.csv',
	'4_DGIdb_Column_5OrMoreTargets_MissingChEMBL_IDsRemoved_1-14-18.csv',
	'5b_DrugCentral_Column_5OrMoreTargets_MissingChEMBL_IDsRemoved_Human_1-14-18.csv')

for interactionlist in interactionlists: interactionlist_to_gvm(interactionlist)

os.chdir('..')