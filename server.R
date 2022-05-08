server = function(input, output, session) {
  
  values = reactiveValues()
  
  output$ref_used = renderText({
    req(input$species)
    req(input$chromosome)
    reference_file = paste("ref/", input$species, "/chr", input$chromosome,"_ref", sep = "")
    values[["reference_file"]] = reference_file
    paste0("Reference file used: \n", reference_file)
  })
  
  output$uploaded_fastq_files = renderDataTable({
    
    req(input$read1)
    req(input$read2)
    
    fastq_names1 = input$read1[["name"]]
    fastq_names2 = input$read2[["name"]]
    
    fastq_sizes1 = input$read1[["size"]]
    fastq_sizes2 = input$read2[["size"]]
    
    fastq_paths1 = input$read1[["datapath"]]
    fastq_paths2 = input$read2[["datapath"]]
    
    uploaded_fastq_info = tibble(
      `R1` = fastq_names1,
      `R1 File Size (Kb)` = fastq_sizes1,
      `R1 Path` = fastq_paths1,
      `R2` = fastq_names2,
      `R2 File Size (Kb)` = fastq_sizes2,
      `R2 Path` = fastq_paths2,
    ) %>%
      mutate(`R1 File Size (Kb)` = round(`R1 File Size (Kb)`/1024, 0)) %>%
      mutate(`R2 File Size (Kb)` = round(`R2 File Size (Kb)`/1024, 0))
    
    values[["uploaded_fastq_info"]] = uploaded_fastq_info
    
    brief_uploaded_fastq_info = uploaded_fastq_info %>%
      select(`R1`, `R1 File Size (Kb)`, `R2`, `R2 File Size (Kb)`)
    brief_uploaded_fastq_info
    
  })
  
  observeEvent(input$start_alignment,{
    
    file.remove(file.path("temp", list.files(path = "temp/")))
    
    withProgress(message = "Performing Short Sequence Alignment...", value = 0, {
      
      uploaded_fastq_info = values[["uploaded_fastq_info"]]
      reference_file = values[["reference_file"]]
      
      read_summary = vector()
      sam_paths = vector()
      for(i in 1:dim(uploaded_fastq_info)[1]){
        
        read1_path = uploaded_fastq_info[["R1 Path"]][i]
        read2_path = uploaded_fastq_info[["R2 Path"]][i]
        
        read1_suffix = gsub(".fastq.gz", "", uploaded_fastq_info[["R1"]][i])
        read2_suffix = gsub(".fastq.gz", "", uploaded_fastq_info[["R2"]][i])
        read1_sam_filename = paste0(read1_suffix, ".sam")
        read2_sam_filename = paste0(read2_suffix, ".sam")
        
        if(all(
          file.exists(paste0(reference_file, ".00.b.array")),
          file.exists(paste0(reference_file, ".00.b.tab")),
          file.exists(paste0(reference_file, ".files")),
          file.exists(paste0(reference_file, ".reads")),
          file.exists(read1_path),
          file.exists(read2_path)
        )){
          align(index = reference_file, readfile1 = read1_path, output_format = "SAM", type = "dna", output_file = file.path("temp", read1_sam_filename))
          align(index = reference_file, readfile1 = read2_path, output_format = "SAM", type = "dna", output_file = file.path("temp", read2_sam_filename))
        }else{
          cat("Reference File(s) and/or Sequencing File(s) Missing...\n")
        }
        
        read1_summary = read.delim(paste0(file.path("temp",read1_suffix), ".sam.summary"),
                                   sep = "\t", header = FALSE)[,2]
        read1_summary = tibble(
          `File` = read1_suffix,
          `Total Reads` = read1_summary[1],
          `Mapped Reads` = read1_summary[2],
          `Uniquely Mapped Reads` = read1_summary[3],
          `Multiple Mapped Reads` = read1_summary[4],
          `Unmapped Reads` = read1_summary[5],
          `Indels` = read1_summary[6]
        )
        
        read2_summary = read.delim(paste0(file.path("temp",read2_suffix), ".sam.summary"),
                                   sep = "\t", header = FALSE)[,2]
        read2_summary = tibble(
          `File` = read2_suffix,
          `Total Reads` = read2_summary[1],
          `Mapped Reads` = read2_summary[2],
          `Uniquely Mapped Reads` = read2_summary[3],
          `Multiple Mapped Reads` = read2_summary[4],
          `Unmapped Reads` = read2_summary[5],
          `Indels` = read2_summary[6]
        )
        
        read_summary = rbind.data.frame(read_summary, read1_summary, read2_summary)
        sam_paths = rbind.data.frame(
          sam_paths,
          tibble(
            r1 = read1_suffix,
            r2 = read2_suffix,
            r1_sam = file.path("temp", read1_sam_filename),
            r2_sam = file.path("temp", read2_sam_filename)
          )
        )
        
        incProgress(1/dim(uploaded_fastq_info)[1], detail = paste("Aligning", uploaded_fastq_info[["R1"]][i],
                                                                  "and", uploaded_fastq_info[["R2"]][i]))
      }
      
    })
    
    values[["read_alignment_summary"]] = read_summary
    values[["sam_paths"]] = sam_paths
    cat("Alignment Complete...\n")
    print(kable(read_summary))
    
  })
  
  output$read_summary_table = renderDataTable({
    
    summary_table = tibble(
      `File` = NULL,
      `Total Reads` = NULL,
      `Mapped Reads` = NULL,
      `Uniquely Mapped Reads` = NULL,
      `Multiple Mapped Reads` = NULL,
      `Unmapped Reads` = NULL,
      `Indels` = NULL
    )
    
    try({
      summary_table = values[["read_alignment_summary"]]
    }, silent = TRUE)
    
    summary_table
  })
  
  output$aligned_sam_table = renderDataTable({
    sam_table = tibble(
      r1 = NULL,
      r2 = NULL,
      r1_sam = NULL,
      r2_sam = NULL
    )
    try({
      sam_table = values[["sam_paths"]]
    }, silent = TRUE)
    
    sam_table
  })
  
  observeEvent(input$update_meta,{
    output$meta_table = renderText({
      
      species = isolate(input$species)
      chromosome = isolate(input$chromosome)
      proto1 = toupper(isolate(input$proto1))
      proto1_be_index = isolate(input$proto1_be_index)
      proto1_be_from = toupper(isolate(input$proto1_be_from))
      proto1_be_to = toupper(isolate(input$proto1_be_to))
      proto1_ct_index = isolate(input$proto1_ct_index)
      proto1_ct_from = toupper(isolate(input$proto1_ct_from))
      proto1_ct_to = toupper(isolate(input$proto1_ct_to))
      
      proto2 = toupper(isolate(input$proto2))
      proto2_be_index = isolate(input$proto2_be_index)
      proto2_be_from = toupper(isolate(input$proto2_be_from))
      proto2_be_to = toupper(isolate(input$proto2_be_to))
      proto2_ct_index = isolate(input$proto2_ct_index)
      proto2_ct_from = toupper(isolate(input$proto2_ct_from))
      proto2_ct_to = toupper(isolate(input$proto2_ct_to))
      
      meta_snippet = tibble(
        species,
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
      
      brief_meta_snippet = tibble(
        `Species` = meta_snippet$species,
        `Chr` = meta_snippet$chromosome,
        `Protospacer 1` = meta_snippet$proto1,
        `Index 1` = paste0("BE: ", meta_snippet$proto1_be_index,"; CT: ", meta_snippet$proto1_ct_index),
        `R1 BE; CT` = paste0(meta_snippet$proto1_be_from,"->",meta_snippet$proto1_be_to,"; ",meta_snippet$proto1_ct_from,"->",meta_snippet$proto1_ct_to),
        `Protospacer 2` = meta_snippet$proto2,
        `Index 2` = paste0("BE: ", meta_snippet$proto2_be_index,"; CT: ", meta_snippet$proto2_ct_index),
        `R2 BE; CT` = paste0(meta_snippet$proto2_be_from,"->",meta_snippet$proto2_be_to,"; ",meta_snippet$proto2_ct_from,"->",meta_snippet$proto2_ct_to)
      )
      
      values[["meta_snippet"]] = meta_snippet
      
      paste0("Meta Parameters: \n",
             paste(unlist(brief_meta_snippet),collapse = "\n"))
      
    })
  })
  
  observeEvent(input$generate_haplotype_table, {
    output$haplotype_table = renderDataTable({
      
      sam_paths = values[["sam_paths"]]
      read1_headers = sam_paths$r1
      read2_headers = sam_paths$r2
      read1_sam_filenames = sam_paths$r1_sam
      read2_sam_filenames = sam_paths$r2_sam
      meta_snippet = values[["meta_snippet"]]
      
      print(kable(sam_paths))
      
      quadruplotypes = tibble(
        quadruplotype = c("p1nbenct_p2nbenct")
      )
      
      withProgress(message = "Performing Short Sequence Alignment...", value = 0, {
        
        for(i in 1:dim(sam_paths)[1]){
          table_header = get_common_string(read1_headers[i], read2_headers[i])
          quadruplotype_snippet = get_quadruplotype(read1_sam_filenames[i],
                                                    read2_sam_filenames[i],
                                                    meta_snippet)
          names(quadruplotype_snippet)[2] = table_header
          quadruplotypes = full_join(quadruplotypes, quadruplotype_snippet,
                                     by = "quadruplotype")
        }
        
        incProgress(1/dim(sam_paths)[1], detail = paste0("Processing ", table_header))
        
      })
      
      quadruplotypes
      
    })
  })
  
  session$onSessionEnded(function() {
    stopApp()
  })
}