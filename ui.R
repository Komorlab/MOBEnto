ui = fluidPage(
  title = "MOBEnto v0.1",
  titlePanel("MOBEnto"),
  
  theme = bs_theme(version = 4, bootswatch = "cosmo"),
  
  tabsetPanel(
    
    tabPanel("Generate Reference", fluid = TRUE,
             br(),
             sidebarLayout(
               sidebarPanel = sidebarPanel(
                 p("Create Customized Reference Sequence Functionality is not Available Yet...\n")
               ),
               mainPanel = mainPanel(
                 reactableOutput("file_tree")
               )
             )
    ),
    
    tabPanel("Perform Alignment", fluid = TRUE,
             br(),
             sidebarLayout(
               sidebarPanel = sidebarPanel(
                 
                 pickerInput(inputId = 'species',
                             label = 'Species',
                             choices = c("Human", "Mouse"),
                             selected = "Human"),
                 
                 pickerInput(inputId = 'chromosome',
                             label = 'Chromosome',
                             choices = c(seq(1,22), "X", "Y", "M"),
                             selected = "1"),
                 
                 fileInput(inputId = "read1", label = "Upload R1 Files (fastq,gz)", 
                           multiple = TRUE, accept = ".fastq.gz"),
                 
                 fileInput(inputId = "read2", label = "Upload R2 Files (fastq,gz)", 
                           multiple = TRUE, accept = ".fastq.gz"),
                 
                 actionButton(inputId = "start_alignment", label = "Run Alignment")
                 
               ),
               mainPanel = mainPanel(
                 verbatimTextOutput("ref_used"),
                 reactableOutput("uploaded_fastq_files"),
                 reactableOutput("read_summary_table")
               )
             )
    ),
    
    tabPanel("Create Metafile", fluid = TRUE,
             br(),
             sidebarLayout(
               sidebarPanel = sidebarPanel(
                 
                 textInput(inputId = "proto1",
                           label = "Protospacer 1",
                           value = "TTCCTGACTTCTGTATGTTG",
                           placeholder = "TTCCTGACTTCTGTATGTTG"),
                 
                 numericInput("proto1_be_index", "Base Edit Site Index",
                              value = 16, min = 1, max = 100),
                 
                 fluidRow(
                   column(width = 6,
                          pickerInput(inputId = "proto1_be_from",
                                      label = "From",
                                      choices = c("A", "T", "C", "G"),
                                      selected = "T")
                   ),
                   
                   column(width = 6,
                          
                          pickerInput(inputId = "proto1_be_to",
                                      label = "To",
                                      choices = c("A", "T", "C", "G"),
                                      selected = "C")
                   )
                 ),
                 
                 numericInput("proto1_ct_index", "Crosstalk Site Index",
                              value = 17, min = 1, max = 100),
                 
                 fluidRow(
                   column(width = 6,
                          pickerInput(inputId = "proto1_ct_from",
                                      label = "From",
                                      choices = c("A", "T", "C", "G"),
                                      selected = "G")
                   ),
                   
                   column(width = 6,
                          
                          pickerInput(inputId = "proto1_ct_to",
                                      label = "To",
                                      choices = c("A", "T", "C", "G"),
                                      selected = "A")
                   )
                 ),
                 
                 br(),
                 
                 textInput(inputId = "proto2",
                           label = "Protospacer 2",
                           value = "CAGGTAATGACTAAGATGAC",
                           placeholder = "CAGGTAATGACTAAGATGAC"),
                 
                 numericInput("proto2_be_index", "Base Edit Site Index",
                              value = 15, min = 1, max = 100),
                 
                 fluidRow(
                   column(width = 6,
                          pickerInput(inputId = "proto2_be_from",
                                      label = "From",
                                      choices = c("A", "T", "C", "G"),
                                      selected = "G")
                   ),
                   
                   column(width = 6,
                          
                          pickerInput(inputId = "proto2_be_to",
                                      label = "To",
                                      choices = c("A", "T", "C", "G"),
                                      selected = "A")
                   )
                 ),
                 
                 numericInput("proto2_ct_index", "Crosstalk Site Index",
                              value = 17, min = 1, max = 100),
                 
                 fluidRow(
                   column(width = 6,
                          pickerInput(inputId = "proto2_ct_from",
                                      label = "From",
                                      choices = c("A", "T", "C", "G"),
                                      selected = "T")
                   ),
                   
                   column(width = 6,
                          
                          pickerInput(inputId = "proto2_ct_to",
                                      label = "To",
                                      choices = c("A", "T", "C", "G"),
                                      selected = "C")
                   )
                 ),
                 
                 actionButton(inputId = "update_meta", "Update Meta")
                 
                 
               ),
               
               
               mainPanel = mainPanel(
                 
                 verbatimTextOutput("meta_table"),
                 reactableOutput("aligned_sam_table")
                 
               )
             )
             
    ),
    
    tabPanel("Quantify Haplotype", fluid = TRUE,
             br(),
             sidebarLayout(
               sidebarPanel = sidebarPanel(
                 
                 actionButton(inputId = "generate_haplotype_table", label = "Quantify Haplotype")
                 
               ),
               mainPanel = mainPanel(
                 
                 reactableOutput("haplotype_table")
                 
               )
             )
    ),
    
  )
  
)
