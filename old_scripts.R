

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

def subset_gmt_and_gvm(filter, axis, gmt_fname, new_gmt_fname, gvm_fname = None, new_gvm_fname = None):
	'''
	
	'''
	gvm = open_gvm(gvm_fname)
	if axis == 'genes':
		#gvm.
		keep = filter(gvm.index.values)
		genes_to_keep = set(gvm.index.values[keep])
		gvm = gvm[keep]
		gvm.to_csv(new_gvm_fname, sep='\t')
		#gmt.
		with open(gmt_fname, 'r') as f:
			reader = csv.reader(f, delimiter = '\t')
			gmt = [row[0:2] + [str(g) for g in row[2:] if g in genes_to_keep] for row in reader]

		with open(new_gmt_fname, 'w', newline='') as f:
			writer = csv.writer(f, delimiter='\t')
			for geneset in gmt: writer.writerow(geneset)

	elif axis == 'annotations':
		print('b')
	else: raise ValueError('axis must be genes or annotations')

dtc = fread('original_drug-gene_libs/DtcDrugTargetInteractions.csv')
dtc$activity_comment = tolower(dtc$activity_comment)
dtc = dtc[grepl('inconclusive|not determined|not active'), dtc$activity_comment]
dtc = subset(dtc, activity_comment 'inconclusive')
dtc = subset(dtc, activity_comment!='not determined')
dtc = dtc[,c('compound_id','target_id','target_pref_name','gene_names','pubmed_id')]

#only need compound id and gene names

lowercase the comment first


not active|compound not obtained|could not be accurately determined due to low solubility|inactive|inconclusive|ineffective|is not a...|no action|no activity|no activity detected|no binding|no block|no change|no data|no delta cs|no detectable activity|no displacement|no effect|no effect at 1 mm|no enzyme activity|no evidence for binding|no inactivation|no incorporation|no increase|no inhibition...|no inhibitory action|no interactions|no loss in activity|no measureable reactivation observed|no reaction|no reactivation|no reproducible apparent ki data|no significant...|no substrate detected|no time-dependent inhibition...|no turn over|non valid test|none|none observed|not a substrate|not active...|not detectable|not detected|not determined|not estimated accurately|not evaluated|not examined|not measureable |not...|unstable|unspecified|unable to be measured|the compound was too insoluble to obtain reliable datas

up = UniProt.ws(taxId=9606)
hgnc_to_uniprotkb = select(up, keys=x, columns='HGNC', keytype='UNIPROTKB')

ds <- parse.gctx("original_drug-gene_libs/CD_signatures_full_42809x22268.gctx", rid=1:10, cid=1:10)

Below is a comparison of how the gene set sizes changed after expansion.

```{r size_comparison, include=TRUE}
gvm_fnames_cleaned = c('DrugBank','TargetCentral','DGIdb','DrugRepHub','DrugCentral')

for(i in 1:length(gvm_fnames)){
  old_gvm_fname = gvm_fnames[[i]]
  old_gvm = fread(paste0('original_drug-gene_libs/', old_gvm_fname, '_gvm.csv'))
  
  ngenes = data.frame(Original=log(colSums(old_gvm==TRUE, na.rm=TRUE), base=10)[-1])
  
  for(j in 1:length(expansion_libs)){
    expansion_lib_name = expansion_libs_names[[j]]
    new_gvm_fname = paste0('expanded_drug-gene_libs/', gvm_libnames[i], '_expanded_with_', expansion_lib_name, '_gvm.csv')
    new_gvm = fread(new_gvm_fname, header=FALSE, showProgress=FALSE)
    ngenes[[expansion_lib_name]] = log(colSums(new_gvm==TRUE, na.rm=TRUE), base=10)[-1]
  }
  ngenes[ngenes==-Inf]=0

  print(ggplot(data=melt(ngenes, measure.vars=c('Original', expansion_libs_names)), aes(x=value, fill=variable)) + geom_histogram(bins=150) + facet_grid(variable~.) + labs(x='log 10 of gene set size', y='count', title=gvm_fnames_cleaned[[i]]))
  
}
```

Here is a summary of how each of the libraries changed.

```{r summarize_expanded_gvms}
libs = c('DrugBank_5OrMoreTargets', 'TargetCentral_5OrMoreTargets', 'DGIdb_5OrMoreTargets', 
	'RepurposeHub_5OrMoreTargets', 'DrugCentral_5OrMoreTargets', 'DTCommons_5OrMoreTargets',
	'STITCH_500cutoff_5OrMoreTargets_noENSP', 'STITCH_600cutoff_5OrMoreTargets_noENSP', 
	'STITCH_700cutoff_5OrMoreTargets_noENSP', 'STITCH_800cutoff_5OrMoreTargets_noENSP', 'CREEDS',
	'LINCS')

for(lib in libs){
  geneset_sizes = read.csv('intermediate_files/lib_summary_geneset_sizes_' + lib + '.csv')
  ggplot(gather(geneset_sizes), aes(value)) + 
    geom_histogram(bins = 50) + 
    facet_wrap(~key, scales = 'free_x')
}
```
expand_gvm = function(gvm, expansion, output_fname){
  #This function creates and saves, if not already done, a gvm file with expanded gene sets
  ##based on an expansion file with gene-gene coexpression data. 
  
  if(file.exists(output_fname)){
    print(paste0(output_fname, ' already created.'))
    return(NULL)
  }
  
  #Get all the genes in either the gvm or the expansion file (or both).
  gene_union = union(as.character(unlist(gvm[,1])), as.character(unlist(expansion)))
  
  #Initialize a matrix to store the expanded gvm.
  expanded_gvm = matrix(NA, ncol=ncol(gvm)-1, nrow=length(gene_union))
  colnames(expanded_gvm) = colnames(gvm)[-1]
  rownames(expanded_gvm) = gene_union
  
  #For each gene set from the old gvm,
  for(c in colnames(expanded_gvm)[-1]){
    #print(c)
    
    #Get all the genes in the original gene set...
    gene_indices = which(gvm[[c]]==TRUE)
    genes = as.character(unlist(gvm[gene_indices,1]))
    
    #and for each of these genes, lookup the corresponding top-n-coexpressed-genes. Combine these sets into a single set.
    e_genes_with_na = as.character(unlist(expansion[genes]))
    
    #Remove NA values.
    e_genes = e_genes_with_na[!is.na(e_genes_with_na)]
    
    #Combine this set with the original genes (some might not have been a lookup value in the expansion file)
    genes_to_add = union(genes, e_genes)

    #And now save this set to the expanded gvm.
    expanded_gvm[genes_to_add, c] = 'True'
  }
  #Save the expanded gvm. 
  write.table(expanded_gvm, output_fname, na='', sep='\t', quote=FALSE) #, row.names=FALSE
}

compare_expanded_gmt = function(old_gmt_fname, expanded_gmt_fname){
  old_gmt = read.csv(old_gmt_fname, sep='\t', colClasses='character', row.names=1)
  expanded_gmt = read.csv(expanded_gmt_fname, sep='\t', colClasses='character', row.names=1)
  
  #Create histograms comparing the distributions of gene set sizes before and after expansion. 
  ngenes = data.frame(old_ngenes=colSums(old_gmt=='True'))
  ngenes$expanded_ngenes = colSums(expanded_gmt=='True')
  print(ggplot(data=melt(ngenes, measure.vars=c('old_ngenes','expanded_ngenes')), aes(x=value, fill=variable)) + geom_histogram(bins=150) + facet_grid(variable~.) + labs(x='gene set size', y='count', title=strsplit(expanded_gmt_fname, '_transformed')[[1]]))
}

# setwd('C:\\Users\\damon\\Desktop\\drug benchmarking\\ppi_libs')
# for(file in c('human_correlation.rda', 'mouse_correlation.rda', 'genename_pairsWprob.txt',
#   'BIOGRID-ORGANISM-Homo_sapiens-3.4.156.tab2.txt')){
#   print(file)
#   get_top_n_coexpressed(file, paste0(strsplit(file, '.', fixed=TRUE)[[1]], '_top_100.csv'))
# }
# setwd('..')

# library(ggplot2)
# library(reshape2)

# for(gmt in c('1_DrugBank_EdgeList_10-05-17', '2_TargetCentral_EdgeList_10-05-17', '3_EdgeLists_Union_10-05-17', '4_EdgeLists_Intersection_10-05-17')){
#   for(expansion in c('genename_pairsWprob_top_100', 'human_correlation_top_100')){
#     if(expansion=='genename_pairsWprob_top_100'){exp_libname='hu.MAP'}
#     else{exp_libname='ARCHS4'}
#     output_fname = paste0('libs\\', gmt, '_expanded_with_', exp_libname, '_transformed.csv')
#     print(output_fname)
#     expand_gmt(gmt_fname=paste0('libs\\', gmt, '_transformed.csv'), expansion_fname=paste0('ppi_libs\\', expansion, '.csv'), output_fname=output_fname)
#     compare_expanded_gmt(old_gmt_fname=paste0('libs\\', gmt, '_transformed.csv'), expanded_gmt_fname=output_fname)
#   }
# }