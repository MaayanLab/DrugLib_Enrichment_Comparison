

		gvm = pd.DataFrame(index=all_genes, dtype=bool)
		gvm_chunk = pd.DataFrame(dtype=bool)

		it = 0
		for annotation in annotations:
			if verbose: print(annotation)
			#Get this annotation's gene vector: the collection of genes with this annotation in their set.
			genes = genesetlist[annotation]
			#Add this annotation's gene vector to the gvm if there are at least five genes.
			vec = pd.DataFrame(True, index=genes, columns=[annotation], dtype=bool)
			gvm_chunk = pd.concat([gvm_chunk,vec], axis=1)

			it += 1
			if it >= 150: 
				gvm = pd.concat([gvm, gvm_chunk], axis=1)
				gvm_chunk = pd.DataFrame(dtype=bool)
				print(str(len(gvm.columns)) + ' out of ' + str(len(annotations)))
				it = 0

		if it != 0:
			gvm = pd.concat([gvm, gvm_chunk], axis=1)

		#Save the results.
		gvm.index.name = lib_name

def convert_gmt(gmt_fname, output_type='gvm', verbose = False, output_fname = None):
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

	def gmt_to_gvm(genesetlist, output_fname):
		'''
		Transforms a gmt .txt file to a gvm dataframe, saves it, and returns it.
		reader : reader object
		output_fname : str
			The output filename, i.e. "ChEA_2016_gvm.csv".
		'''
		lib_name = output_fname.partition('_gvm')[0].rsplit('/', maxsplit=1)[1]

		return genesetlist_to_gvm(genesetlist, output_fname=output_fname)

def genesetlist_to_gvm(genesetlist, output_fname = None, verbose = False, GVM_CHUNKSIZE = 200):

		gvm = pd.DataFrame(dtype=bool)
		gvm_chunk = pd.DataFrame(dtype=bool)
		it = 0
		#Each row corresponds to an annotation (e.g. a tf or drug). So, for each annotation:
		for i in interactions.index:
			annotation = row[0]
			if verbose: print(annotation)
			if annotation in MY_NA_VALS: next
			#Obtain the gene set, removing delimiting strings and not accepting null values.
			genes = interactions[i]
			#Remove repeats. To save time, comment out this below line if gene sets are known to be unique.
			genes = set(genes)
			#Create a column gene vector for this annotation...
			vec = pd.DataFrame(True, index = genes, columns = [annotation], dtype=bool)
			#...and concatenate it with the growing dataframe. 
			gvm_chunk = pd.concat([gvm_chunk,vec], axis=1)

			it += 1
			if it > 200: 
				print(annotation)
				gvm = pd.concat([gvm, gvm_chunk], axis=1)
				gvm_chunk = pd.DataFrame(dtype=bool)
				it = 0

		if it != 0:
			gvm = pd.concat([gvm, gvm_chunk], axis=1)

		gvm.index.name = lib_name
		#Save file to the current working directory, with values True and '' (instead of False, to save space).
		gvm.to_csv(output_fname, sep='\t')
		#Return the dataframe in sparse form, with values True and False.
		gvm = gvm.fillna(False).to_sparse()
		return gvm

	def gmt_to_dict(reader):
		#TO FIX: REMOVE TOO FEW GENES THINGY U KNO
		'''
		Converts a gmt .txt file to a dict, and returns it.
		reader : reader object
		'''
		d = {row[0]:[str(g).replace(',1.0','') for g in row[2:] if g not in MY_NA_VALS] for r in reader}
		d = {k:v for (k,v) in d.items() if len(v) > 4}
		return d

	#If the gvm is requested, check to see if it has already been created. If so, simply load it and return it.
	if output_type == 'gvm':
		if output_fname is None: output_fname = gmt_fname.replace('gmt','gvm').replace('.txt','.csv')
		if file_exists(output_fname): 
			gvm = open_gvm(output_fname)
			gvm.index.name = output_fname.rpartition('/')[2].partition('_gvm.csv')[0]
			return gvm

	#Else, open the gmt file.
	if not os.path.isfile(gmt_fname): raise ValueError(gmt_fname, 'could not be found.')	
	with open(gmt_fname, 'r') as f:
		reader = csv.reader(f, delimiter = '\t')
		out = gmt_to_dict(reader)
		if output_type == 'dict': return out

		out = pd.Series(out)
		if output_type == 'gvm': return gmt_to_gvm(out, output_fname)
		elif output_type == 'dict': return gmt_to_dict(reader)
		else: raise ValueError(output_type, 'must be either gvm or dict.')

		mask = [len(x) > 4 for x in genesetlist]
		genesetlist = genesetlist[mask]

def gvm_to_gmt(gvm_fname, gmt_fname = None, duplicates_might_exist=True):
	if gmt_fname is None: gmt_fname = gvm_fname.replace('gvm', 'gmt')
	if file_exists(gmt_fname): return
	gvm = open_gvm(gvm_fname)
	gmt_rows = [gvm.index.values[list(compress(range(len(gvm[t])), gvm[t]))] for t in gvm.columns.values]

	#If duplicates might exist, we will use the string before "." as the annotation name. 
	#This will cause problems if "." might actually be part of the annotation names!
	#Because of this, gvms with "." in column names and some duplicate column names are not supported.
	if duplicates_might_exist: colnames = [c.partition('.')[0] for c in gvm.columns.values]
	else: colnames = gvm.columns.values

	gmt_rows = [[t, None] + list(r) for (t, r) in zip(colnames, gmt_rows)]
	with open(gmt_fname, "w") as out:
		wr = csv.writer(out, delimiter='\t', lineterminator='\n')
		for r in gmt_rows:
			wr.writerow(r)