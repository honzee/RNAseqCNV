####Multi-tab user interface####
shinyAppUi <- navbarPage("RNAseq CNA analysis", id = "tabs",

                           #Input tab

                           tabPanel("Input",

                                    fluidPage(

                                      titlePanel("Input"),

                                      sidebarLayout(

                                        sidebarPanel(
                                          fileInput("config", "config file"),
                                          htmlOutput("mess_config"),
                                          br(),
                                          fileInput("metadata", "metadata files"),
                                          htmlOutput("mess_metadata"),
                                          br(),
                                          radioButtons("snv_format", "Select input data format for SNV information", choiceNames = c("vcf", "custom table"), choiceValues = c("vcf", "custom"), inline = TRUE),
                                          br(),
                                          radioButtons("genome_version", "Select genome version according to the input files", choiceNames = c("hg38", "hg19"), choiceValues = c("hg38", "hg19"), inline = TRUE),
                                          br(),
                                          checkboxInput("adjust_in", "Apply diploid adjustement", value = TRUE),
                                          br(),
                                          checkboxInput("batch", "Analyze samples as a batch", value = FALSE),
                                          br(),
                                          checkboxInput("generate_weights", "Generate gene weights from the analyzed samples", value = FALSE),
                                          br(),
                                          checkboxInput("arm_lvl", "Generate arm level figures (increases run-time)", value = TRUE),
                                          br(),
                                          checkboxInput("estimate_lab", "Plot with estimate labels", value = TRUE),
                                          br(),
                                          actionButton("preview", "Analyze first sample"),
                                          br(),
                                          br(),
                                          actionButton("analyze", "Analyze all samples"),
                                          br(),
                                          br(),
                                          shinyDirButton("dir_button", "Mock analysis", "Please select output directory for mock analysis")
                                        ),

                                        mainPanel(
                                          imageOutput("main_fig_prev", width = "100%", height = "auto"),
                                          conditionalPanel(condition = "output.chr_fig_prev != null",
                                            fluidRow(
                                              column(2, actionButton("prev_butt_chr_prev", "Previous", width = "100%")
                                              ),
                                              column(2, offset = 8, actionButton("next_butt_chr_prev", "Next", width = "100%"))
                                            )
                                          ),
                                          imageOutput("chr_fig_prev",  width = "100%", height = "auto")

        )
      )
    )
  )
)

