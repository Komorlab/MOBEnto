read_sam = function(file_rx){
  f = readLines(con = file_rx)
  f = f[nchar(f)>100]
  f = f[-1]
  f = strsplit(f ,split = "\t") %>%
    map(function(x){
      t(x)[c(1,10)]
    }) 
  f = tibble(barcode = unlist(lapply(f, function(x){x[1]})),
             aligned_sequence = unlist(lapply(f, function(x){x[2]}))) %>%
    filter(nchar(aligned_sequence)>0)
  return(f)
}