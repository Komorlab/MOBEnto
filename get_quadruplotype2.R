library(dplyr)
library(purrr)
library(knitr)
library(writexl)

get_quadruplotype_from_snippet = function(meta_snippet){
  
  cat("-----------------------------------------------------\n")
  cat(paste0("Quantifying Quadruplotypes for: \n"))
  filename = gsub(pattern = "_R1_", "_RX_", meta_snippet$r1)
  cat(paste0(filename,"...\n"))
  
  # Read san file
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
  
  # Combine read 1 and read 2
  file_r1 = read_sam(file_rx = paste0(meta_snippet$path, meta_snippet$r1))
  file_r2 = read_sam(file_rx = paste0(meta_snippet$path, meta_snippet$r2))
  sam_file = inner_join(file_r1, file_r2, by = "barcode", suffix = c("_r1", "_r2")) %>%
    mutate(aligned_sequence = paste0(aligned_sequence_r1, aligned_sequence_r2))
  
  # Calculate total reads
  total_read = dim(sam_file)[1]
  
  # Format protospacer sequence
  proto1 = toupper(meta_snippet$proto1)
  proto2 = toupper(meta_snippet$proto2)
  
  # Filter sam files based on common protospacer start
  proto_start1 = strsplit(sam_file$aligned_sequence, split = proto1) %>%
    lapply(function(x){x[1]}) %>%
    unlist() %>%
    lapply(nchar) %>%
    unlist()
  proto_start1_demax = proto_start1[proto_start1!=max(proto_start1)]
  proto_common_start1 = names(table(proto_start1_demax)[which.max(table(proto_start1_demax))]) %>%
    as.numeric()
  sam_file_filtered1 = sam_file[!(proto_start1 %in% c((proto_common_start1-5):(proto_common_start1-1),
                                                      (proto_common_start1+1):(proto_common_start1+5))),]
  
  proto_start2 = strsplit(sam_file$aligned_sequence, split = proto2) %>%
    lapply(function(x){x[1]}) %>%
    unlist() %>%
    lapply(nchar) %>%
    unlist()
  proto_start2_demax = proto_start2[proto_start2!=max(proto_start2)]
  proto_common_start2 = names(table(proto_start2_demax)[which.max(table(proto_start2_demax))]) %>%
    as.numeric()
  sam_file_filtered2 = sam_file[!(proto_start2 %in% c((proto_common_start2-5):(proto_common_start2-1),
                                                      (proto_common_start2+1):(proto_common_start2+5))),]
  
  # Retrieve genotypes
  gt1 = tibble(
    proto1_be_gt = sapply(sam_file_filtered1$aligned_sequence, strsplit, split = "") %>%
      map(function(x){
        x[c(proto_common_start1+meta_snippet$proto1_be_index)]
      }) %>%
      unlist(),
    proto1_ct_gt = sapply(sam_file_filtered1$aligned_sequence, strsplit, split = "") %>%
      map(function(x){
        x[c(proto_common_start1+meta_snippet$proto1_ct_index)]
      }) %>%
      unlist()
  ) %>%
    mutate(barcode = sam_file_filtered1$barcode)
  
  gt2 = tibble(
    proto2_be_gt = sapply(sam_file_filtered2$aligned_sequence, strsplit, split = "") %>%
      map(function(x){
        x[c(proto_common_start2+meta_snippet$proto2_be_index)]
      }) %>%
      unlist(),
    proto2_ct_gt = sapply(sam_file_filtered2$aligned_sequence, strsplit, split = "") %>%
      map(function(x){
        x[c(proto_common_start2+meta_snippet$proto2_ct_index)]
      }) %>%
      unlist()
  ) %>%
    mutate(barcode = sam_file_filtered2$barcode)
  
  # Generate quadruplotype summary
  quadruplotype_summary = inner_join(gt1, gt2, by = "barcode") %>%
    select(barcode, everything()) %>%
    filter(proto1_be_gt==toupper(meta_snippet$proto1_be_from)|proto1_be_gt==toupper(meta_snippet$proto1_be_to)) %>%
    filter(proto2_be_gt==toupper(meta_snippet$proto2_be_from)|proto2_be_gt==toupper(meta_snippet$proto2_be_to)) %>%
    filter(proto1_ct_gt==toupper(meta_snippet$proto1_ct_from)|proto1_ct_gt==toupper(meta_snippet$proto1_ct_to)) %>%
    filter(proto2_ct_gt==toupper(meta_snippet$proto2_ct_from)|proto2_ct_gt==toupper(meta_snippet$proto2_ct_to)) %>%
    mutate(proto1_be = toupper(meta_snippet$proto1_be_from)) %>%
    mutate(proto1_ct = toupper(meta_snippet$proto1_ct_from)) %>%
    mutate(proto2_be = toupper(meta_snippet$proto2_be_from)) %>%
    mutate(proto2_ct = toupper(meta_snippet$proto2_ct_from)) %>%
    mutate(proto1_be_status = ifelse(proto1_be_gt==proto1_be,"nbe","be")) %>%
    mutate(proto1_ct_status = ifelse(proto1_ct_gt==proto1_ct,"nct","ct")) %>%
    mutate(proto2_be_status = ifelse(proto2_be_gt==proto2_be,"nbe","be")) %>%
    mutate(proto2_ct_status = ifelse(proto2_ct_gt==proto2_ct,"nct","ct")) %>%
    mutate(quadruplotype = paste(proto1_be_gt, proto1_ct_gt, proto2_be_gt, proto2_ct_gt, sep = "")) %>%
    mutate(quadruplotype_status = paste0("p1",proto1_be_status,proto1_ct_status,"_","p2",proto2_be_status,proto2_ct_status))
  
  quadruplotype_status = table(quadruplotype_summary$quadruplotype_status)
  quadruplotype_status = tibble(quadruplotype = names(quadruplotype_status),
                                count = quadruplotype_status)
  names(quadruplotype_status)[2] = filename
  print(kable(quadruplotype_status))
  
  return(quadruplotype_status)
  
}


get_quadruplotype_from_meta = function(meta_dir = "emx1_meta.csv"){
  
  meta = read.csv(file = meta_dir)
  
  quadruplotypes = tibble(
    quadruplotype = c("p1nbenct_p2nbenct")
  )
  
  for(i in 1:dim(meta)[1]){
    try({
      quadruplotypes = quadruplotypes %>%
        full_join(get_quadruplotype_from_snippet(meta[i,]), by = "quadruplotype")
    })
  }
  
  table_name = paste0("quadruplotypes_from_", meta_dir)
  table_name = gsub(pattern = ".csv", replacement = ".xlsx", table_name)
  write_xlsx(quadruplotypes, path = table_name)
  
  return(quadruplotypes)
  
}



