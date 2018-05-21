import pandas as pd
import numpy as np
import time
import sys
import os
import csv
import timeit
from convert_files_scripts import *

lib_to_convert = sys.argv[1]

if lib_to_convert == 'fivedtlibs':
	interactionlist_fnames = (
		'intermediate_files/1_DrugBank_Column_5OrMoreTargets_MissingChEMBL_IDsRemoved_1-14-18_PCIDs.csv',
		'intermediate_files/2_TargetCentral_Column_5OrMoreTargets_MissingChEMBL_IDsRemoved_1-14-18_PCIDs.csv',
		'intermediate_files/3_RepurposeHub_Column_5OrMoreTargets_MissingChEMBL_IDsRemoved_1-14-18_PCIDs.csv',
		'intermediate_files/4_DGIdb_Column_5OrMoreTargets_MissingChEMBL_IDsRemoved_1-14-18_PCIDs.csv',
		'intermediate_files/5b_DrugCentral_Column_5OrMoreTargets_MissingChEMBL_IDsRemoved_Human_1-14-18_PCIDs.csv'
		)
	libnames = ('DrugBank','TargetCentral','RepurposeHub','DGIdb','DrugCentral')
	gslists = [get_genesetlist(get_interactionlist(i_fname), 'ilist') for i_fname in interactionlist_fnames]
	for i in range(len(gslists)): 
		gslist = gslists[i]
		libname = libnames[i] + '_5OrMoreTargets'
		gmt_fname = 'gmts/' + libname + '_gmt.csv'
		gvm_fname = 'gvms/' + libname + '_gvm.csv'
		print('getting files for ' + libname)
		get_gmt_and_gvm(gslist, gmt_fname, gvm_fname)

elif lib_to_convert == 'STITCH':
	CONFIDENCE_SCORES = [800, 700, 600, 500]

	CONFIDENCE_SCORES = sorted(CONFIDENCE_SCORES)
	STITCH = get_interactionlist('9606.protein_chemical.links.v5.0.tsv', CONFIDENCE_SCORES)

	for conf_score in CONFIDENCE_SCORES:
		print('creating STITCH files for interactions at least ' + str(conf_score) + ' confidence.')
		prefix = 'STITCH_' + str(conf_score)
		STITCH = STITCH.loc[STITCH['combined_score'] >= conf_score]

		gslist = get_genesetlist(STITCH[['gene','annotation']], 'ilist')
		gslist_noENSP = get_genesetlist(STITCH.loc[ ~STITCH['gene'].str.match('ENSP') ][['gene','annotation']], 'ilist')

		#THE GVM FILES ARE TOO BIG IF DRUGS WITH LESS THAN FIVE TARGETS ARE KEPT.
		#SO, JUST GET THE GMT.
		gmt_fname = 'gmts/' + prefix + 'cutoff_gmt.csv'
		convert_genesetlist(gslist, to='gmt', output_fname=gmt_fname)
		#get_gmt_and_gvm(gslist, gmt_fname)
		#No ENSP.
		gmt_fname = 'gmts/' + prefix + 'cutoff_noENSP_gmt.csv'
		convert_genesetlist(gslist_noENSP, to='gmt', output_fname=gmt_fname)
		#get_gmt_and_gvm(gslist_noENSP, gmt_fname)

		#5 or more targets.
		gmt_fname = 'gmts/' + prefix + 'cutoff_5OrMoreTargets_gmt.csv'
		get_gmt_and_gvm(remove_small_genesets(gslist), gmt_fname)
		#5 or more targets and no ENSP.
		gmt_fname = 'gmts/' + prefix + 'cutoff_5OrMoreTargets_noENSP_gmt.csv'
		get_gmt_and_gvm(remove_small_genesets(gslist_noENSP), gmt_fname)

elif lib_to_convert == 'DTCommons':
	gslist = get_genesetlist(
		get_interactionlist('intermediate_files/DTCommons_interactionlist.txt'), item_type='ilist')
	convert_genesetlist(gslist, to='gmt', output_fname='gmts/DTCommons_gmt.csv')
	get_gmt_and_gvm(remove_small_genesets(gslist), 'gmts/DTCommons_5OrMoreTargets_gmt.csv')

elif lib_to_convert == 'pertlibs':
	CREEDS_up = get_genesetlist(
		'original_drug-gene_libs/Drug_Perturbations_from_GEO_up.txt', item_type='gmt_fname')
	CREEDS_down = get_genesetlist(
		'original_drug-gene_libs/Drug_Perturbations_from_GEO_down.txt', item_type='gmt_fname')
	CREEDS = combine_genesetlists(CREEDS_up, CREEDS_down)
	get_gmt_and_gvm(CREEDS, 'intermediate_files/CREEDS_noPCIDs_gmt.csv')
	annotations = pd.Series(CREEDS.index.values)
	annotations.name = 'annotations'
	annotations.to_csv('intermediate_files/CREEDS_original_annotations.csv', header=True, index=False, sep='\n')

	LINCS_up = get_genesetlist(
		'original_drug-gene_libs/LINCS_L1000_Chem_Pert_up.txt', item_type='gmt_fname')
	LINCS_down = get_genesetlist(
		'original_drug-gene_libs/LINCS_L1000_Chem_Pert_down.txt', item_type='gmt_fname')
	LINCS = combine_genesetlists(LINCS_up, LINCS_down)
	get_gmt_and_gvm(LINCS,'intermediate_files/LINCS_noPCIDs_gmt.csv')
	annotations = pd.Series(LINCS.index.values)
	annotations.name = 'annotations'
	annotations.to_csv('intermediate_files/LINCS_original_annotations.csv', header=True, index=False, sep='\n')

elif lib_to_convert == 'LINCS2':
	LINCS = get_genesetlist('intermediate_files/LINCS_noPCIDs_gmt.csv', item_type='gmt_fname')
	LINCS_identifiers = pd.read_csv('synonyms/LINCS_identifiers.csv', sep='\t')
	if not all(LINCS.index.values == LINCS_identifiers['sample']): raise ValueError('Mismatched annotations.')
	LINCS.index = LINCS_identifiers['gvm']
	get_gmt_and_gvm(LINCS, 'gmts/LINCS_gmt.csv')

elif lib_to_convert == 'expanded_libs':
	gvm_fnames = [fname for fname in os.listdir('gvms') if 'gvm.' in fname]
	gvm_fnames = [fname for fname in gvm_fnames if (('CREEDS' not in fname) and ('LINCS' not in fname))]
	expansion_fnames = ['ARCHS4_human_top_100.csv', 'BioGRID_top_100.csv', 'huMAP_top_100.csv']

	for gvm_fname in gvm_fnames:
		gvm = open_gvm('gvms//' + gvm_fname)
		for expansion_fname in expansion_fnames:
			expansion = pd.read_csv('ppi-coexp_libs//' + expansion_fname, sep='\t')
			expansion = expansion.set_index('gene').transpose()
			output_fname = 'gvms//expanded//' + gvm_fname.partition('_gvm.')[0] + '_expanded_with_' + expansion_fname.partition('_top')[0] + '_gvm.csv'
			if 'csv' not in gvm_fname: output_fname = output_fname.replace('.csv','.pkl')
			print(output_fname)
			expand_gvm(gvm, expansion, output_fname)
			convert_genesetlist(get_genesetlist(output_fname, 'gvm_fname'), to='gmt', output_fname=output_fname.replace('gvm','gmt'))

else: raise ValueError('Invalid library to convert argument.')