####Multi-tab user interface####
shinyAppUi <- navbarPage("RNAseq CNA analysis", id = "tabs",

                         #Input tab

                         tabPanel("Input",

                                  fluidPage(

                                    titlePanel("Input"),

                                    sidebarLayout(

                                      sidebarPanel(
                                        fileInput("metadata", "metadata files"),
                                        htmlOutput("mess_metadata"),
                                        br(),
                                        fileInput("config", "config file"),
                                        htmlOutput("mess_config"),
                                        checkboxInput("adjust_in", "Apply diploid adjustement", value = TRUE),
                                        br(),
                                        checkboxInput("arm_lvl", "Generate arm level figures (increases run-time)", value = TRUE),
                                        br(),
                                        checkboxInput("estimate", "Plot with estimate labels", value = TRUE),
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
                         ),
                         tabPanel("Manual CNV analysis",

                                  fluidPage(

                                    tags$style(type = "text/css",
                                               "label {font-size: 16px;}"
                                    ),

                                    titlePanel("Manual Analysis"),
                                    br(),
                                    br(),
                                    fluidRow(
                                      column(3,
                                             wellPanel(
                                               h3("Select Figures"),
                                               br(),
                                               uiOutput("figure_select"),
                                               br(),
                                               uiOutput("chr_sel"),
                                               br(),
                                               h3("Estimation correction"),
                                               br(),
                                               uiOutput("type_select"),
                                               br(),
                                               br(),
                                               uiOutput("gender_select"),
                                               br(),
                                               br(),
                                               uiOutput("chromn_text"),
                                               br(),
                                               br(),
                                               uiOutput("alt_text"),
                                               htmlOutput("examp"),
                                               htmlOutput("war_message"),
                                               br(),
                                               br(),
                                               uiOutput("comments_text"),
                                               br(),
                                               br(),
                                               fluidRow(
                                                 column(6,
                                                        uiOutput("default")
                                                 ),
                                                 column(3, offset = 3,
                                                        uiOutput("save_butt")
                                                 )
                                               ),
                                               fluidRow(
                                                 column(3, offset = 9,
                                                        htmlOutput("status"))
                                               )
                                             )
                                      ),
                                      column(9,
                                             fluidRow(
                                               column(2,
                                                      uiOutput("prev_butt")
                                               ),
                                               column(8,
                                                      htmlOutput("sample_num")
                                               ),
                                               column(2,
                                                      uiOutput("next_butt"))
                                             ),
                                             fluidRow(
                                               column(12,
                                                      imageOutput("main_fig", width = "100%", height = "auto")
                                               )
                                             ),
                                             conditionalPanel(condition = "output.chr_choices != null",
                                                              fluidRow(
                                                                column(2,
                                                                       uiOutput("prev_butt_chr")
                                                                ),
                                                                column(2, offset = 8,
                                                                       uiOutput("next_butt_chr"))
                                                              ),
                                                              fluidRow(
                                                                column(12,
                                                                       imageOutput("chr_fig", width = "100%", height = "auto")
                                                                )
                                                              )
                                             )
                                      )
                                    )
                                  )
                         ),
                         tabPanel("Export",

                                  fluidPage(

                                    titlePanel("Export analyzed table"),

                                    sidebarLayout(

                                      sidebarPanel(

                                        uiOutput("columns"),
                                        br(),
                                        br(),
                                        uiOutput("format"),
                                        br(),
                                        br(),
                                        shinyDirButton("export", label = "Export to selected directory", title = "Select directory")

                                      ),

                                      mainPanel(

                                        h2("Output preview"),
                                        br(),
                                        br(),
                                        tableOutput("prev_tab")

                                      )

                                    )
                                  ))
)
