import os
import numpy as np
import pandas as pd
from convert_files_scripts import *

libs = ['DrugBank_5OrMoreTargets', 'TargetCentral_5OrMoreTargets', 'DGIdb_5OrMoreTargets', 
	'RepurposeHub_5OrMoreTargets', 'DrugCentral_5OrMoreTargets', 'DTCommons_5OrMoreTargets',
	'STITCH_500cutoff_5OrMoreTargets_noENSP', 'STITCH_600cutoff_5OrMoreTargets_noENSP', 
	'STITCH_700cutoff_5OrMoreTargets_noENSP', 'STITCH_800cutoff_5OrMoreTargets_noENSP', 'CREEDS', 'LINCS']

names = ['DrugBank','Target Central Resource Database', 'Drug Gene Interaction Database', 
	'Drug Repurposing Hub', 'Drug Central', 'Drug Target Commons',
	'Search Tool for Interacting Chemicals (STITCH), 500 cutoff', 'STITCH, 600 cutoff',
	'STITCH, 700 cutoff', ' STITCH, 800 cutoff', 
	'CRowd Extracted Expression of Differential Signatures',
	'Library of Integrated Network-Based Cellular Signatures']

abbrevs = [i.replace('_5OrMoreTargets','').replace('_noENSP','') for i in libs]

link = ['drugbank.ca', 'juniper.health.unm.edu/tcrd', 'dgidb.org', 'clue.io/repurposing',
	'drugcentral.org', 'drugtargetcommons.fimm.fi', 'STITCH.embl.de', 'STITCH.embl.de', 
	'STITCH.embl.de', 'STITCH.embl.de', 'amp.pharm.mssm.edu/CREEDS', 'lincsproject.org']

lib_type = ['Drug Target'] * 10 + ['Perturbation Signature'] * 2

description = ['Encyclopedia-like drug database compiled from scientific publications.',
  'Part of the "Illuminating the Druggable Genome" project, which compiles data from multiple sources and emphasizes GPCR, kinase, ion channel and nuclear receptor targets.',
  'Drug-gene interactions compiled by curation and text-mining from multiple sources.',
  'Drug screening library containing known clinical drugs and their targets as reported by scientific publications.',
  'Database of drugs approved or regulated by government agencies.',
  'Drug bioactivity database extracted and organized from multiple sources by crowdsource.',
  'Chemical-protein interaction network obtained from computational predictions, inter-organism knowledge, and other databases.',
  'See above.', 'See above.', 'See above.',
  'Human, mouse and rat drug perturbation signatures extracted from the Gene Expression Omnibus by crowdsource.',
  'Gene expression profiles from flow cytometry experiments with different perturbagens, time points, doses, and cell lines.']

n_samples = []
n_unique_drugs = []
n_genes = []
avg_n_genes_per_annot = []

preprocessing = ['Edgelist. Removed samples missing the ChEMBL_ID, or with less than five targets.'] * 5 + ['Converted UniProt protein IDs to HGNC symbols. Removed samples with blank, negative, and inconclusive activities, mutant cell lines, missing an HGNC symbol, or with less than five targets.',
	'Converted ENSP IDs to HGNC symbols. Used an interaction score cutoff of 500. Removed samples with less than five targets.',
	'Interaction score cutoff of 600', 'Interaction score cutoff of 700', 'Interaction score cutoff of 800',
	'Took the union of up and down DE gene sets from Enrichr gmts.',
	'Took the union of up and down DE gene sets from Enrichr gmts.']

all_drugs = pd.Series(index=abbrevs)

for i in range(len(libs)):
	lib = libs[i]
	print(lib)
	fname = 'gvms/' + lib + '_gvm.csv'
	if not os.path.isfile(fname): fname = fname.replace('.csv','.pkl')
	if not os.path.isfile(fname): raise ValueError('Unknown file', fname)
	gvm = open_gvm(fname)
	samples = [x.partition('.')[0] for x in gvm.columns.values]
	n_samples.append(str(len(samples)))
	drugs = set(samples)
	all_drugs.iloc[i] = drugs
	n_unique_drugs.append(str(len(drugs)))
	n_genes.append(str(gvm.shape[0]))
	avg_n_genes_per_annot.append(str(round(sum(sum(gvm.values))/int(n_samples[-1]), 2)))

table = pd.DataFrame([names, link, lib_type, description, n_samples, n_unique_drugs, n_genes, avg_n_genes_per_annot, preprocessing]).transpose()
table.columns = ['name', 'link', 'lib_type', 'description', 'n_samples', 'n_unique_drugs', 'n_genes', 'avg_n_genes_per_annot', 'preprocessing']
table.to_csv('intermediate_files/lib_summary_table.csv', sep='\t', index=False)

matches = pd.DataFrame(index = abbrevs, columns = abbrevs)
for a in abbrevs:
	for b in abbrevs:
		matches.at[a,b] = len(all_drugs[a].intersection(all_drugs[b]))
matches.to_csv('intermediate_files/lib_summary_matches.csv', sep='\t', index=False)

ppi_coexp_libs = ['Original','ARCHS4_human','BioGRID','huMAP']
for lib in libs:
	if lib == 'CREEDS' or lib == 'LINCS': continue
	geneset_size_fname = 'intermediate_files/lib_summary_' + lib + '_geneset_sizes.csv'
	for ppi_coexp_lib in ppi_coexp_libs:
		if ppi_coexp_lib == 'Original': 
			fname = 'gvms/' + lib + '_gvm.csv'
		else: fname = 'gvms/expanded/' + lib + '_expanded_with_' + ppi_coexp_lib + '_gvm.csv'
		if not os.path.isfile(fname): fname = fname.replace('.csv','.pkl')
		if not os.path.isfile(fname): raise ValueError('Unknown file', fname)
		gvm = open_gvm(fname)
		if ppi_coexp_lib == 'Original': geneset_sizes = pd.DataFrame(index = ppi_coexp_libs, columns = gvm.columns)
		if type(gvm) == pd.core.sparse.frame.SparseDataFrame: geneset_sizes.loc[ppi_coexp_lib,:] = gvm.to_dense().sum()
		else: geneset_sizes.loc[ppi_coexp_lib,:] = gvm.sum()
	geneset_sizes.to_csv('intermediate_files/lib_summary_geneset_sizes_' + lib + '.csv', sep='\t', index=False)

