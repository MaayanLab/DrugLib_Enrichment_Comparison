#================================================================
#This script expands the gene sets in a gmt file.
#================================================================ 
get_top_n_coexpressed = function(expansion_library_fname, output_fname, n=100){
  #This function creates and saves, if not already done, and returns a lookup table 
  ##whose rows are the genes in the input expansion library, and whose columns
  ##are the top n genes, in order, which are coexpressed with the corresponding row gene. 
  ##The expansion libraries which are supported are the gene-gene correlation datasets from ARCHS4 and hu.MAP.
  print(output_fname)
  
  if(file.exists(output_fname)){
  	print('already created.')
  	return(NULL)
  }
  
  #Different expansion library files are formatted differently. So, process accordingly. 
  
  #For the ARCHS4 files:
  if(expansion_library_fname %in% c('human_correlation.rda', 'mouse_correlation.rda')){
    #This file is large, so break the it down into ten chunks.
    ##Each row can be computed independently, so take the top tenth rows, and then the next tenth, and so on. 
    #Begin by initializing a list to store the results from each chunk.
    rows_list = vector('list', length=10)
    for(i in 0:9){
      #Load the expansion library data.table. 
      print(i)
      temp_space = new.env()
      f = load(expansion_library_fname, temp_space)
      expansion_library = get(f, temp_space)
      rm(temp_space)
      
      #Take the ith chunk. 
      start = floor(nrow(expansion_library)/10*i) + 1
      if(i==9){
        end=nrow(expansion_library)} else {
          end = floor(nrow(expansion_library)/10*(i + 1))}
      expansion_library = expansion_library[start:end,]
      
      #Define a function that returns the top n genes
      ##with the highest coexpression values, for a given row.
      get_top_n = function(row){
        return(names(row)[order(row, decreasing=TRUE)[1:n]])
      }
      
      #Get the top n genes for each row.
      rows_list[[i+1]] = t(apply(expansion_library, 1, get_top_n))
    }
    
    #Combine the results of each chunk.
    top_n = do.call(rbind, rows_list)
    print(dim(top_n))
    
    #Remove the (now duplicated) list of chunk results to save memory.
    rm(rows_list)
    
    #For the hu.MAP files:
  } else if(expansion_library_fname %in% c('genename_pairsWprob.txt')) {
    expansion_library = read.table(expansion_library_fname)
    colnames(expansion_library) = c('gene1', 'gene2', 'p')
    
    #Remove rows with p-value == 0.
    expansion_library = subset(expansion_library, p > 0)
    
    #Get a vector with all the genes in this expansion library...
    genes = union(expansion_library$gene1, expansion_library$gene2)
    
    #...and then create a list for each of these genes which we will use to store
    ##the genes which are coexpressed.
    coexpressed = as.list(genes)
    names(coexpressed) = genes
    
    #For each coexpression pair, add gene A to the coexpression list of gene B
    ##and add gene B to the coexpression list of gene A.
    ##Since the expansion library is already sorted by ascending p value,
    ##The first n genes in these lists will be the top n coexpressed.
    for(i in 1:nrow(expansion_library)){
      g1 = as.character(expansion_library[[i,'gene1']])
      g2 = as.character(expansion_library[[i,'gene2']])
      
      coexpressed[[g1]] = append(coexpressed[[g1]], g2)
      coexpressed[[g2]] = append(coexpressed[[g2]], g1)
    }
    
    #Only keep the top n coexpressed.
    for(i in names(coexpressed)){
      if(length(coexpressed[[i]]) > n){
        coexpressed[[i]]=coexpressed[[i]][1:n]
      }
      length(coexpressed[[i]]) = 100
    }
    
    #Combine the lists into a data.frame
    top_n = do.call(rbind, coexpressed)
  }
  
  #Save the coexpression data.frame to file and return.
  write.table(top_n, output_fname, sep='\t')
  return(top_n)
  
}

expand_gmt = function(gmt_fname, expansion_fname, output_fname){
  #This function creates and saves, if not already done, a gmt file with expanded gene sets
  ##Based on an expansion file with gene-gene coexpression data. 
  
  if(file.exists(output_fname)){
    print(paste0(output_fname, ' already created.'))
    return(NULL)
  }
  
  gmt=read.csv(gmt_fname, sep='\t', colClasses='character', row.names=1)
  expansion = read.csv(expansion_fname, sep='\t', colClasses='character', row.names=1)
  
  #Get all the genes in either the gmt or the expansion file (or both).
  gene_union = union(rownames(gmt), rownames(expansion))
  
  #Initialize a matrix to store the expanded gmt.
  expanded_gmt = matrix(NA, ncol=ncol(gmt), nrow=length(gene_union))
  colnames(expanded_gmt) = colnames(gmt)
  rownames(expanded_gmt) = gene_union
  
  #For each gene set from the old gmt,
  for(c in colnames(expanded_gmt)){
    #print(c)
    
    #Get all the genes in the original gene set,
    genes = as.character(rownames(gmt[gmt[c]=='True',]))
    
    #And for each of these genes, lookup the corresponding top-n-coexpressed-genes. Combine these sets into a single set.
    e_genes_with_na = as.character(unlist(expansion[genes,]))
    
    #Remove NA values.
    e_genes = e_genes_with_na[!is.na(e_genes_with_na)]
    
    #Combine this set with the original genes (some might not have been a lookup value in the expansion file)
    genes_to_add = union(genes, e_genes)
    
    #And now save this set to the expanded gmt.
    expanded_gmt[genes_to_add, c] = 'True'
  }
  #Save the expanded gmt. 
  write.table(expanded_gmt, output_fname, na='', sep='\t', quote=FALSE)
}

compare_expanded_gmt = function(old_gmt_fname, expanded_gmt_fname){
  old_gmt = read.csv(old_gmt_fname, sep='\t', colClasses='character', row.names=1)
  expanded_gmt = read.csv(expanded_gmt_fname, sep='\t', colClasses='character', row.names=1)
  
  #Create histograms comparing the distributions of gene set sizes before and after expansion. 
  ngenes = data.frame(old_ngenes=colSums(old_gmt=='True'))
  ngenes$expanded_ngenes = colSums(expanded_gmt=='True')
  print(ggplot(data=melt(ngenes, measure.vars=c('old_ngenes','expanded_ngenes')), aes(x=value, fill=variable)) + geom_histogram(bins=150) + facet_grid(variable~.) + labs(x='gene set size', y='count', title=strsplit(expanded_gmt_fname, '_transformed')[[1]]))
  }

setwd('C:\\Users\\damon\\Desktop\\drug benchmarking\\gene_coexp_data')
for(file in c('human_correlation.rda', 'mouse_correlation.rda', 'genename_pairsWprob.txt')){
  print(file)
  get_top_n_coexpressed(file, paste0(strsplit(file, '.', fixed=TRUE)[[1]], '_top_100.csv'))
}
setwd('..')

library(ggplot2)
library(reshape2)

for(gmt in c('1_DrugBank_EdgeList_10-05-17', '2_TargetCentral_EdgeList_10-05-17', '3_EdgeLists_Union_10-05-17', '4_EdgeLists_Intersection_10-05-17')){
  for(expansion in c('genename_pairsWprob_top_100', 'human_correlation_top_100')){
    if(expansion=='genename_pairsWprob_top_100'){exp_libname='hu.MAP'}
    else{exp_libname='ARCHS4'}
    output_fname = paste0('libs\\', gmt, '_expanded_with_', exp_libname, '_transformed.csv')
    print(output_fname)
    expand_gmt(gmt_fname=paste0('libs\\', gmt, '_transformed.csv'), expansion_fname=paste0('gene_coexp_data\\', expansion, '.csv'), output_fname=output_fname)
    compare_expanded_gmt(old_gmt_fname=paste0('libs\\', gmt, '_transformed.csv'), expanded_gmt_fname=output_fname)
  }
}