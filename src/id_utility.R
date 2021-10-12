

entrez_to_symbol <- function(vec,hgnc){
  hgnc %<>% dplyr::select(.,symbol,entrez_id)
  hgnc$entrez_id %<>% as.character()
  hgnc %<>% subset(.,!is.na(entrez_id))
  return(hgnc$symbol[match(vec,hgnc$entrez_id)])
}

entrez_to_cds <- function(vec,hgnc){
  hgnc %<>% dplyr::select(.,symbol,entrez_id)
  hgnc$entrez_id %<>% as.character()
  hgnc %<>% subset(.,!is.na(entrez_id))
  hgnc$cds <- paste0(hgnc$symbol," (",hgnc$entrez_id,")")
  return(hgnc$cds[match(vec,hgnc$entrez_id)])
}

extract_entrez <- function(vec){
  parens <- unlist(regmatches(vec, gregexpr("\\(.*?\\)", vec)))
  ids <- gsub("[\\(\\)]", "", parens)
  
  if (all(is.na(ids))){
    stop("Gene names do not appear to be in CDS format, all entries are NA")
  }
  return(ids)
}