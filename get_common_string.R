get_common_string = function(a, b){
  sb = stri_sub(b, 1, 1:nchar(b))
  sstr = na.omit(stri_extract_all_coll(a, sb, simplify=TRUE))
  sstr = sstr[which.max(nchar(sstr))]
  return(sstr)
}
