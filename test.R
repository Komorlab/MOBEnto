library(dplyr)
library(purrr)
library(knitr)
library(Rsubread)

read_filelist = tibble(
  read1 = c("Q-Rn-A-OBE1-QTC281-Tfx-210909_S37_L001_R1_001.fastq.gz"),
  read2 = c("Q-Rn-A-OBE1-QTC281-Tfx-210909_S37_L001_R2_001.fastq.gz")
)

spacies = "Human"
chromosome = 1
proto1 = "ttcctgacttctgtatgttg"
proto1_be_index	= 16
proto1_be_from = "t"
proto1_be_to = "c"
proto1_ct_index	= 17
proto1_ct_from = "g"
proto1_ct_to = "a"
proto2 = "caggtaatgactaagatgac"
proto2_be_index	= 15
proto2_be_from = "g"
proto2_be_to = "a"
proto2_ct_index	= 17
proto2_ct_from = "t"
proto2_ct_to = "c"

proto1 = toupper(proto1)
proto2 = toupper(proto2)
proto1_be_from = toupper(proto1_be_from)
proto1_be_to = toupper(proto1_be_to)
proto1_ct_from = toupper(proto1_ct_from)
proto1_ct_to = toupper(proto1_ct_to)
proto2_be_from = toupper(proto2_be_from)
proto2_be_to = toupper(proto2_be_to)
proto2_ct_from = toupper(proto2_ct_from)
proto2_ct_to = toupper(proto2_ct_to)

meta_snippet = tibble(
  spacies,
  chromosome,
  proto1,
  proto1_be_index,	
  proto1_be_from, 
  proto1_be_to, 
  proto1_ct_index,	
  proto1_ct_from,
  proto1_ct_to, 
  proto2,
  proto2_be_index,	
  proto2_be_from,
  proto2_be_to,
  proto2_ct_index,
  proto2_ct_from,
  proto2_ct_to
)


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

get_quadruplotype = function(read1, read2, meta_snippet){
  
  # Generate alignment parameters
  reference_file = paste("ref/", meta_snippet$spacies, "/chr", meta_snippet$chromosome,"_ref", sep = "")
  read1_suffix = gsub(".fastq.gz", "", read1)
  read2_suffix = gsub(".fastq.gz", "", read2)
  read1_sam_filename = paste0(read1_suffix, ".sam")
  read2_sam_filename = paste0(read2_suffix, ".sam")
  
  # Perform alignment
  align(index = reference_file, readfile1 = read1, output_format = "SAM", type = "dna", output_file = read1_sam_filename)
  align(index = reference_file, readfile1 = read2, output_format = "SAM", type = "dna", output_file = read2_sam_filename)
  
  # Load aligned files
  read1_sequence = read_sam(read1_sam_filename)
  read2_sequence = read_sam(read2_sam_filename)
  
  # Merge pair-end reads
  sam_file = inner_join(read1_sequence, read2_sequence, by = "barcode", suffix = c("_r1", "_r2")) %>%
    mutate(aligned_sequence = paste0(aligned_sequence_r1, aligned_sequence_r2))
  
  # Count total reads
  total_read = dim(sam_file)[1]
  
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
                                count = quadruplotype_status)
  
  return(quadruplotype_status)
}

merge_quadruplotype_result = function(read_filelist, meta_snippet){
  
  quadruplotypes = tibble(
    quadruplotype = c("p1nbenct_p2nbenct")
  )
  
  for(i in 1:dim(read_filelist)[1]){
    try({
      read1 = read_filelist$read1[i]
      read2 = read_filelist$read2[i]
      quadruplotype_snippet = get_quadruplotype(read1 = read1, read2 = read2,
                                                meta_snippet = meta_snippet)
      quadruplotypes = quadruplotypes %>%
        full_join(quadruplotype_snippet, by = "quadruplotype")
    })
  }
  
  return(quadruplotypes)
}

merge_quadruplotype_result(read_filelist, meta_snippet)
