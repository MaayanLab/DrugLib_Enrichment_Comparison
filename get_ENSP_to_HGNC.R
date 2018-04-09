list.of.packages = c('knitr', 'data.table', 'biomaRt', 'webchem', 'httr', 'jsonlite')
new.packages = list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, require, character.only = TRUE)
rm(list.of.packages, new.packages)

setwd('/Users/damon/Desktop/drug benchmarking')
#Read STITCH file to get a list of ENSP IDs to convert.
stitch = fread('original_drug-gene_libs/9606.protein_chemical.links.v5.0.tsv')

#Get ENSP-to-HGNC lookup table. 
gene_symbol_lookup_fname = 'synonyms/ENSP_to_HGNC.csv'
ensembl = useMart('ENSEMBL_MART_ENSEMBL',dataset='hsapiens_gene_ensembl')
gene_symbol_lookup = getBM(filters= "ensembl_peptide_id",
                           attributes= c("ensembl_peptide_id","hgnc_symbol"),
                           #Ask for HGNC symbols for each ENSP ID in STITCH.
                           values=unique(stitch$protein),
                           mart=ensembl)
#Replace blank symbol values with the ENSP.
gene_symbol_lookup[gene_symbol_lookup$hgnc_symbol=='','hgnc_symbol'] =
  gene_symbol_lookup[gene_symbol_lookup$hgnc_symbol=='','ensembl_peptide_id']
#Save lookup table to csv.
write.table(gene_symbol_lookup, file=gene_symbol_lookup_fname, row.names=FALSE, quote=FALSE)