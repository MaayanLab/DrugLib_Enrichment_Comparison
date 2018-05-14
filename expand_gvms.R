#================================================================
#This script expands the gene sets in a gmt file.
#================================================================ 

get_top_n_coexpressed = function(expansion_library_fname, output_fname, n=100){
  #This function creates and saves, if not already done, and returns a lookup table 
  ##whose rows are the genes in the input expansion library, and whose columns
  ##are the top n genes, in order, which are coexpressed with the corresponding row gene. 
  ##The expansion libraries which are supported are the gene-gene correlation datasets from ARCHS4, hu.MAP and BioGRIDs.
  
  if(file.exists(output_fname)){
    print(paste0(output_fname, ' is already created.'))
    return(NULL)
  }
  
  #Different expansion library files are formatted differently. So, process accordingly. 
  
  ARCHS4_top_n = function(fname, n){
    #This file is large, so break the it down into ten chunks.
    ##Each row can be computed independently, so take the top tenth rows, and then the next tenth, and so on. 
    #Begin by initializing a list to store the results from each chunk.
    rows_list = vector('list', length=10)
    for(i in 0:9){
      #Load the expansion library data.table. 
      print(i)
      temp_space = new.env()
      f = load(fname, temp_space)
      ARCHS4 = get(f, temp_space)
      rm(temp_space)
      
      #Take the ith chunk. 
      start = floor(nrow(ARCHS4)/10*i) + 1
      if(i==9){
        end=nrow(ARCHS4)} else {
          end = floor(nrow(ARCHS4)/10*(i + 1))}
      ARCHS4 = ARCHS4[start:end,]
      
      #Define a function that returns the top n genes
      ##with the highest coexpression values, for a given row.
      get_top_n = function(row){
        return(names(row)[order(row, decreasing=TRUE)[1:(n+1)]])
      }
      
      #Get the top n genes for each row.
      rows_list[[i+1]] = t(apply(ARCHS4, 1, get_top_n))
    }
    
    #Combine the results of each chunk.
    top_n = do.call(rbind, rows_list)
    print(dim(top_n))
    
    #Remove the (now duplicated) list of chunk results to save memory.
    rm(rows_list)

    return(top_n)
  }

  hu.MAP_top_n = function(fname, n){
    hu.MAP = read.table(fname)

    colnames(hu.MAP) = c('gene1', 'gene2', 'p')
    hu.MAP['gene1'] = toupper(hu.MAP[['gene1']])
    hu.MAP['gene2'] = toupper(hu.MAP[['gene2']])
    
    #Remove rows with probability score  == 0.
    hu.MAP = subset(hu.MAP, p > 0)
    
    #Get a vector with all the genes in this expansion library...
    genes = union(hu.MAP$gene1, hu.MAP$gene2)
    
    #...and then create a list for each of these genes which we will use to store
    ##the genes which are coexpressed.
    coexpressed = as.list(genes)
    names(coexpressed) = genes
    
    #For each coexpression pair, add gene A to the coexpression list of gene B
    ##and add gene B to the coexpression list of gene A.
    ##Since the expansion library is already sorted by ascending p value,
    ##The first n genes in these lists will be the top n coexpressed.
    for(i in 1:nrow(hu.MAP)){
      g1 = as.character(hu.MAP[[i,'gene1']])
      g2 = as.character(hu.MAP[[i,'gene2']])
      
      coexpressed[[g1]] = append(coexpressed[[g1]], g2)
      coexpressed[[g2]] = append(coexpressed[[g2]], g1)
    }
    
    #Only keep the top n coexpressed.
    for(i in names(coexpressed)){
      if(length(coexpressed[[i]]) > n){
        coexpressed[[i]]=coexpressed[[i]][1:(n+1)]
      }
      length(coexpressed[[i]]) = n+1
    }
    
    #Combine the lists into a data.frame
    top_n = do.call(rbind, coexpressed)

    return(top_n)
  }

  BioGRID_top_n = function(fname, n){
    biogrid = fread(fname)
    biogrid = biogrid[Throughput == 'Low Throughput', c('Official Symbol Interactor A','Official Symbol Interactor B')]
    colnames(biogrid) = c('gene1','gene2')
    biogrid[,'gene1'] = toupper(unlist(biogrid[,'gene1']))
    biogrid[,'gene2'] = toupper(unlist(biogrid[,'gene2']))

    #Get a vector with all the genes in this expansion library...
    genes = union(biogrid[['gene1']], biogrid[['gene2']])
    
    #...and then create a list for each of these genes which we will use to store
    ##the genes which are coexpressed.
    coexpressed = as.list(genes)
    names(coexpressed) = genes
    
    #For each coexpression pair, add gene A to the coexpression list of gene B
    ##and add gene B to the coexpression list of gene A.
    ##Since the expansion library is already sorted by ascending p value,
    ##The first n genes in these lists will be the top n coexpressed.
    for(i in 1:nrow(biogrid)){
      g1 = as.character(biogrid[[i,'gene1']])
      g2 = as.character(biogrid[[i,'gene2']])
      
      coexpressed[[g1]] = append(coexpressed[[g1]], g2)
      coexpressed[[g2]] = append(coexpressed[[g2]], g1)
    }
    
    #BioGRID has repeats, so remove them.
    coexpressed = sapply(coexpressed, unique)

    #Only keep the top n coexpressed.
    for(i in names(coexpressed)){
      if(length(coexpressed[[i]]) > n){
        coexpressed[[i]]=coexpressed[[i]][1:(n+1)]
      }
      length(coexpressed[[i]]) = n+1
    }
    
    #Combine the lists into a data.frame
    top_n = do.call(rbind, coexpressed)

    return(top_n)
  }

  if(expansion_library_fname %in% c('ppi-coexp_libs\\human_correlation.rda', 'ppi-coexp_libs\\mouse_correlation.rda')){
    top_n = ARCHS4_top_n(expansion_library_fname, n)
  } else if(expansion_library_fname == 'ppi-coexp_libs\\genename_pairsWprob.txt') {
    top_n = hu.MAP_top_n(expansion_library_fname, n)
  } else if(expansion_library_fname == 'ppi-coexp_libs\\BIOGRID-ORGANISM-Homo_sapiens-3.4.160.tab2.txt'){
    top_n = BioGRID_top_n(expansion_library_fname, n)
  } else { stop(c('Unknown expansion library file name ', expansion_library_fname)) }

  colnames(top_n) = c('gene', paste0('int-coexp_', as.character(1:n)))

  #Save the coexpression data.frame to file and return.
  write.table(top_n, output_fname, sep = '\t', quote = FALSE, row.names = FALSE)
  return(top_n)
  
}

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