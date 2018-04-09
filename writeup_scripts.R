clean_drug_names = function(drug_names){
	#Only keep alphanumeric characters. Lowercase all letters.
	return(gsub('[^[:alnum:]]', '', tolower(drug_names)))
}


extract_embedded_synonyms = function(drug_names, left_paren_split = ' \\(', no_syns_substrs = NULL){
	#Identifies drug names containing a synonym enclosed in parentheses. Replaces the drug name with the first synonym, and creates a table with the first synonym in one column and the second synonym in another. The first synonym is defined as the substring to the left of `left_paren_split`, and the second synonym is defined as the substring between `left_paren_split` and the right-most ")". For example, the drug name "chembl1255588 ((+)-butaclamol)" will become "chembl1255588" and ('chembl1255588', '(+)-butaclamol') will be added to the synonym table. 
	#Returns the input `drug_names` vector, but with second synonyms removed, and a data.table of the synonyms. 
	#left_paren_split : the substring used to identify the left parenthesis enclosing the second synonym. Usually, these are separated from the actual drug name by a space, so searching for ' \\(' instead of '\\(' prevents names like "4-(benzyloxy)phenol" from being detected as having a synonym.
	#no_syns_substrs : a vector of substrings. If any are present in a drug name, the drug name will not be identified as having a synonym. For example, 'methyl \\(' will prevent "methyl (2z)-..." from being identified. One problem is that it will also prevent "methyl (2z)-... (true-synonym)" from being identified, even though it does contain a synonym. Luckily, this case does not occur in the given libraries. 
	
	#Identify names with potential synonyms.
	if(is.null(no_syns_substrs)){
		might_have_synonym = grepl(left_paren_split, drug_names)
	} else {
		might_have_synonym = grepl(left_paren_split, drug_names) & !grepl(paste(no_syns_substrs, collapse='|'), drug_names)
	}
	
	#Exit if no potential synonyms are detected.
	if(!any(might_have_synonym)){
		return(list(drug_names, NA))
	}
	
	#Confirm potential synonyms by checking each, one at a time.
	drug_names.might_have_synonym = drug_names[might_have_synonym]
	synonym_list = matrix(data=NA,nrow=length(drug_names.might_have_synonym), ncol=2)
	for(i in 1:length(drug_names.might_have_synonym)){
		name = drug_names.might_have_synonym[i]
		name_split = strsplit(name, split=left_paren_split)[[1]] #Split on first occurence of `left_paren_split`
		syn1 = name_split[1] #The first synonym is the substring to the left of `left_paren`split`.
		syn2 = paste0(name_split[-1]) #The second synonym is the substring to the right...
		syn2 = gsub("(.*)\\)(.*)", "\\1\\2", syn2) #...without the right-most right parenthesis.
		
		if(length(syn2) > 1){
			pick = menu(syn2, graphics = FALSE, 'Select one:')
			syn2 = syn2[pick]
		}
		
		#Currently, the only other condition for a potential synonym to be a true synonym is
		#that both the true drug name and its synonym must have three or more alphanumeric characters.
		cleaned_syns = clean_drug_names(c(syn1, syn2))
		if(all(nchar(cleaned_syns) > 2)){
			#If true synonym, put the first synonym in the first column, and the second in the second.
			synonym_list[i,] = c(syn1, syn2)
		} else {
			#If false synonym, put the original name in the first column, and NA in the second.
			synonym_list[i,] = c(name, NA)
		}
	}
	
	#Replace the potential synonyms.
	drug_names[might_have_synonym] = synonym_list[,1]
	return(list(drug_names, unique(na.omit(synonym_list))))
}


do_chunkwise = function(FUN, x, out_fname = 'temp.txt', CHUNKSIZE = 1000, ATTEMPTS = 100, DELETE_OUTFILE = TRUE, ...){
	#if(is.null(out_fname)){stop('out_fname must be given.')}
	
	chunks = split(x, ceiling(seq_along(x)/CHUNKSIZE))

	next_chunk = 1
	if(file.exists(out_fname)){
		print(paste0('using saved results from ', out_fname))
		temp_out = read.table(out_fname, sep='\t', quote = '', header = FALSE, comment.char = "")
		nrows = nrow(temp_out)
		if(nrows == length(x)){return(temp_out[,2])}
		stopifnot(nrows %% CHUNKSIZE == 0)
		next_chunk = as.integer(nrows/CHUNKSIZE) + 1 
	}
	
	for(i in next_chunk:length(chunks)){
		print(paste0(i, ' out of ', length(chunks), ' chunks'))
		query = chunks[[i]]
	
		out=NULL
		attempt = 0
		while(is.null(out) && attempt <= ATTEMPTS) {
			attempt = attempt + 1
			try(
				out <- FUN(query, ...)
			)
		}

		if(is.null(out)){stop('Too many failed attempts.')}
		
		out.table = data.frame(firstcol_input = query, secondcol_output = out)
		write.table(x=out.table, file = out_fname, sep = '\t', append = TRUE, quote = FALSE,
								col.names = FALSE, row.names = FALSE)
	}
	
	out.table.complete = read.table(out_fname, sep='\t', quote = '', header = FALSE, comment.char = "")

	if(nrow(out.table.complete) != length(x)){stop('An error occured.')}
	if(DELETE_OUTFILE){file.remove(out_fname)}

	return(out.table.complete[,2])
}

#ADAPTED FROM WEBCHEM. THE CHANGE IS THAT JSON QUERIES WILL CONTAIN MULTIPLE COMPOUNDS IF POSSIBLE. THIS IS FASTER.
pc_synonyms <- function(query, from = 'name', interactive = 0, verbose = TRUE, CHUNKSIZE=10, arg = NULL, ...) {
	# from can be cid | name | smiles | inchi | sdf | inchikey | formula
	# query <- c('Aspirin')
	# from = 'name'
	
	foo <- function(query, from, verbose, ...){
		prolog <- 'https://pubchem.ncbi.nlm.nih.gov/rest/pug'
		input <- paste0('/compound/', from)
		output <- '/synonyms/JSON'
		if (!is.null(arg))
			arg <- paste0('?', arg)
		qurl <- paste0(prolog, input, output, arg)
		if (verbose)
			message(qurl)
		Sys.sleep(0.2)
		
		query_commasep = paste0(query,collapse=',')
		
		cont <- try(content(POST(qurl,
														 body = paste0(from, '=', query_commasep)
		)), silent = TRUE
		)
		
		if (inherits(cont, "try-error")) {
			if (length(query) > 1) {
				warning('Problem with web service encountered... Trying chemicals individually.')
				cont <- lapply(query, foo, from=from, verbose=verbose)
			} else {
				warning('Problem with web service encountered... Returning NA.')
				return(NA)
			}
		} else if (names(cont) == 'Fault') {
			if (length(query) > 1) {
				warning(cont$Fault$Details, '. Trying chemicals individually.')
				cont <- lapply(query, foo, from=from, verbose=verbose)
			} else {
				warning(cont$Fault$Details, '. Returning NA.')
				return(NA)
			}
		} else if (length(query) > 1) {
			cont <- cont$InformationList$Information
			if (length(cont) > length(query)) {
				#One or more chemicals map to multiple CIDs.
				#The current results cannot distinguish which, so redo each chemical individually.
				cont <- lapply(query, foo, from=from, verbose=verbose)
			}
		}
		
		out <- lapply(lapply(cont, unlist), unname)
		names(out) <- query

		if (interactive > 0) {
			for (i in 1:length(out)) {
				results <- out[[i]]
				if (length(results) > 1) {
					pick <- menu(results[seq_len(interactive)], graphics = FALSE, 'Select one:')
					out[[i]] <- results[pick]          
				}
			}     
		}
		return(out)
	}
	
	if(from %in% c('name','smiles','inchi','inchikey','sdf')){
		#Only single strings per request are suppoprted for these identifiers.
		out <- lapply(query, foo, from=from, verbose=verbose)
		out <- do.call(c, out)
	} else if (length(query) > CHUNKSIZE) {
		query_chunks <- split(query, ceiling(seq_along(query)/CHUNKSIZE))
		out <- lapply(query_chunks, foo, from=from, verbose=verbose)
		out <- do.call(c, out)
	} else {
		out <- foo(query, from = from, verbose = verbose)
	}
	
	out <- setNames(out, query)
	return(out)
}

read_table = function(fname){
	read.table(fname, sep='\t', quote='', header=TRUE, comment.char='')
}

read_old_table = function(fname){
	read.table(fname, sep='\t', quote='', header=FALSE, comment.char='')	
}

write_table = function(table, fname, ...){
	write.table(table, file = fname, sep = '\t', quote = FALSE, row.names = FALSE, ...)
}

write_gvm = function(gvm, fname){
	write.table(gvm, file = fname, sep = '\t', na='', quote = FALSE, row.names = FALSE)
}