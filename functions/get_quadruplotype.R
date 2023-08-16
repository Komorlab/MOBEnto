get_quadruplotype = function(read1_sam_filename, read2_sam_filename, meta_snippet){
  
  # Load aligned files
  read1_sequence = read_sam(read1_sam_filename)
  read2_sequence = read_sam(read2_sam_filename)
  
  # Merge pair-end reads
  sam_file = inner_join(read1_sequence, read2_sequence, by = "barcode", suffix = c("_r1", "_r2")) %>%
    mutate(aligned_sequence = paste0(aligned_sequence_r1, aligned_sequence_r2))
  
  # Filter sam files based on common protospacer index
  proto_start1 = strsplit(sam_file$aligned_sequence, split = meta_snippet$proto1) %>%
    lapply(function(x){x[1]}) %>%
    unlist() %>%
    lapply(nchar) %>%
    unlist()
  if(length(unique(proto_start1))!=1){
    proto_start1_demax = proto_start1[proto_start1!=max(proto_start1)]
  }else{
    proto_start1_demax = max(proto_start1)
  }
  proto_common_start1 = names(table(proto_start1_demax)[which.max(table(proto_start1_demax))]) %>%
    as.numeric()
  sam_file_filtered1 = sam_file[!(proto_start1 %in% c((proto_common_start1-5):(proto_common_start1-1),
                                                      (proto_common_start1+1):(proto_common_start1+5))),]
  
  proto_start2 = strsplit(sam_file$aligned_sequence, split = meta_snippet$proto2) %>%
    lapply(function(x){x[1]}) %>%
    unlist() %>%
    lapply(nchar) %>%
    unlist()
  if(length(unique(proto_start2))!=1){
    proto_start2_demax = proto_start2[proto_start2!=max(proto_start2)]
  }else{
    proto_start2_demax = max(proto_start2)
  }
  proto_common_start2 = names(table(proto_start2_demax)[which.max(table(proto_start2_demax))]) %>%
    as.numeric()
  sam_file_filtered2 = sam_file[!(proto_start2 %in% c((proto_common_start2-5):(proto_common_start2-1),
                                                      (proto_common_start2+1):(proto_common_start2+5))),]
  
  # Retrieve genotypes at each specified sites
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
    filter(proto1_be_gt==meta_snippet$proto1_be_from|proto1_be_gt==meta_snippet$proto1_be_to) %>%
    filter(proto2_be_gt==meta_snippet$proto2_be_from|proto2_be_gt==meta_snippet$proto2_be_to) %>%
    filter(proto1_ct_gt==meta_snippet$proto1_ct_from|proto1_ct_gt==meta_snippet$proto1_ct_to) %>%
    filter(proto2_ct_gt==meta_snippet$proto2_ct_from|proto2_ct_gt==meta_snippet$proto2_ct_to) %>%
    mutate(proto1_be = meta_snippet$proto1_be_from) %>%
    mutate(proto1_ct = meta_snippet$proto1_ct_from) %>%
    mutate(proto2_be = meta_snippet$proto2_be_from) %>%
    mutate(proto2_ct = meta_snippet$proto2_ct_from) %>%
    mutate(proto1_be_status = ifelse(proto1_be_gt==proto1_be,"nbe","be")) %>%
    mutate(proto1_ct_status = ifelse(proto1_ct_gt==proto1_ct,"nct","ct")) %>%
    mutate(proto2_be_status = ifelse(proto2_be_gt==proto2_be,"nbe","be")) %>%
    mutate(proto2_ct_status = ifelse(proto2_ct_gt==proto2_ct,"nct","ct")) %>%
    mutate(quadruplotype = paste(proto1_be_gt, proto1_ct_gt, proto2_be_gt, proto2_ct_gt, sep = "")) %>%
    mutate(quadruplotype_status = paste0("p1",proto1_be_status,proto1_ct_status,"_","p2",proto2_be_status,proto2_ct_status))
  
  quadruplotype_status = table(quadruplotype_summary$quadruplotype_status)
  quadruplotype_status = tibble(quadruplotype = names(quadruplotype_status),
                                count = as.numeric(quadruplotype_status))
  
  return(quadruplotype_status)
  
}