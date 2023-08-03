# convert the refseq id to HSGC gene symbol
RefSeq2HSGC <- function(refseq_id){
  require(biomaRt)
  mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  
  mapping <- getBM(
    attributes = c("ensembl_gene_id", "hgnc_symbol", "refseq_mrna", "refseq_ncrna"),
    mart=mart,
    filters = "refseq_mrna",
    values = refseq_id,
    uniqueRows = T
  )
  remaining_id <- refseq_id[!(refseq_id %in% mapping$refseq_mrna)]
  mapping <- rbind(mapping,
                   getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", "refseq_mrna", "refseq_ncrna"),
                         mart=mart,
                         filters = "refseq_ncrna",
                         values = remaining_id,
                         uniqueRows = T)
  )
  return(mapping$hgnc_symbol)
}

Ens2HSGC <- function(ensembl_id){
  require(biomaRt)
  mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  
  mapping <- getBM(
    attributes = c("ensembl_gene_id_version", "hgnc_symbol"),
    mart=mart,
    filters = "ensembl_gene_id_version",
    values = ensembl_id,
    uniqueRows = T
  )
  return(mapping$hgnc_symbol)
}

HGNC2Des <- function(hgnc_name){
  require(biomaRt)
  mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  
  mapping <- getBM(
    attributes = c("hgnc_symbol", "description"),
    mart=mart,
    filters = "hgnc_symbol",
    values = hgnc_name,
    uniqueRows = T
  )
  return(mapping$description)
}

Ens2Des <- function(ensembl_id){
  require(biomaRt)
  mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  
  mapping <- getBM(
    attributes = c("ensembl_gene_id_version", "description"),
    mart=mart,
    filters = "ensembl_gene_id_version",
    values = ensembl_id,
    uniqueRows = T
  )
  return(mapping$description)
}

# GTEx variant IDs to rsIDs

GTEx2rs <- function(variant_ids, reference){

  #check the input data is provided
  if(missing(variant_ids)){
    stop("Please specify a vector of GTEx variant IDs.")
  }
  if(missing(reference)){
    stop("Please provide a reference table. Reference can be downloaded here: https://storage.googleapis.com/gtex_analysis_v8/reference/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.txt.gz")
  }
  stopifnot(('variant_id' %in% names(reference)) & ('rs_id' %in% names(reference)))
  #mapping variant ids to rs ids
  variant_ids.df = data.frame('variant_id' = variant_ids)
  map = merge(variant_ids.df, reference, by = "variant_id", all.x = TRUE)
  return(map$rs_id)
}

rs2GTEx <- function(rsids, reference){
  #check the input data is provided
  if(missing(rsids)){
    stop("Please specify a vector of rs IDs.")
  }
  if(missing(reference) ){
    stop("Please provide a reference table. Reference can be downloaded here: https://storage.googleapis.com/gtex_analysis_v8/reference/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.txt.gz")
  }
  stopifnot(('variant_id' %in% names(reference)) & ('rs_id' %in% names(reference)))

  #mapping rs ids to variant ids
  rsids.df = data.frame('rs_id' = rsids)
  map = merge(rsids.df, reference, by = "rs_id", all.x = TRUE)
  return(map$variant_id)
}

# from Human gene symbol to mouse
Human2Mouse <- function(x){
  require("biomaRt")
  human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/") 
  mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
  
  genesV2 <- getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", 
                    values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  
  mousex <- unique(genesV2[, 2]) # unique the 2nd column values
  
  human_genes_number <- length(x)
  mouse_genes_number <- length(mousex)
  
  if(human_genes_number != mouse_genes_number){
    genes_not_trans <- setdiff(x, genesV2$HGNC.symbol)
    print("These genes could not be translated:")
    print(genes_not_trans)
    print(paste("A total number of ",length(genes_not_trans),"genes could not be translated!"),sep=" ")
  }else{
    print("All genes were translated successfully!")
  }
  
  return(mousex)
}
