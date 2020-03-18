####Server####
shinyAppServer <- function(input, output, session) {

  #change the size of possible upload
  options(shiny.maxRequestSize=30*1024^2)
  volumes <- c(Home = fs::path_home(), "R Installation" = R.home(), getVolumes()())

  #do the directories from the input exist?
  observe({
  if (!is.null(input$config)) {
    source(input$config$datapath)

      #check count_dir
      if (!dir.exists(count_dir)) {
        count <- "not_found"
      } else {
        count <- count_dir
      }

      #check snv_dir
      if (!dir.exists(snv_dir)) {
        snv <- "not_found"
      } else {
        snv <- snv_dir
      }

      #check out_dir
      if (!dir.exists(out_dir)) {
        out <- "not_found"
      } else {
        out <- out_dir
      }

      react_val$config <- c(count_dir = count, snv_dir = snv, out_dir = out)

    } else {
      react_val$config <- "no_input"
    }
  })


  #render warning messages if directories don't exist
  output$mess_config <- renderText({
    if(any(react_val$config %in% "not_found")) {

      message <- NULL

      for (i in 1:length(react_val$config)) {
        if (react_val$config[i] == "not_found") {
          message <- paste0(message, '<br><font size="1px"><p style="text-align:left;"><font color="red">', names(react_val$config)[i], " was not found")
        }
      }

      return(message)

    } else {
      NULL
    }

  })


  #check metadata file for three columns and read it
  observe({

    if (is.null(input$metadata)) return("no_input")

    metadata_tab = fread(input$metadata$datapath, header = FALSE)

    if (ncol(metadata_tab) == 3) {
      react_val$metadata <- metadata_tab
    } else {
      react_val$metadata <- "incorrect_format"
    }
  })

  #render warning messages if metadata table does not have three columns
  output$mess_metadata <- renderText({

    if (all( react_val$metadata == "no_input")) return(NULL)

    if(all( react_val$metadata == "incorrect_format")) {

      return('<font size="1px"><p style="text-align:left;"><font color="red">metadata table does not have three columns')

    }
  })


  #create a sample table
  sample_table <- reactive({

    if (all(react_val$metadata == "no_input")) return(NULL)
    if (all(react_val$config == "no_input")) return(NULL)


    if (any(react_val$config == "not_found")) return(NULL)
    if (all(react_val$metadata == "incorrect_format")) return(NULL)

    sample_table <-  react_val$metadata

    HTSeq_f = pull(sample_table, 2)
    snv_f = pull(sample_table, 3)

    #create paths to the files
    sample_table$count_path <- file.path(react_val$config["count_dir"], HTSeq_f)
    sample_table$snv_path <-  file.path(react_val$config["snv_dir"], snv_f)
    return(sample_table)
  })


  #check which files are available
  avail <- eventReactive(c(input$preview,
                           input$analyze), {
                             if (is.null(input$metadata) | is.null(input$config)) return(NULL)
                             count_ex <- sample_table()$count_path[!file.exists(sample_table()$count_path)]
                             snv_ex <- sample_table()$snv_path[!file.exists(sample_table()$snv_path)]
                             if (length(count_ex) > 0 | length(snv_ex) > 0) {
                               return(paste0("Missing files: ", paste(paste0(count_ex, collapse = ", "), paste(snv_ex, collapse = ", "), sep = ", ")))
                             } else {
                               "all_present"
                             }
                           })

  # ####Generate a figure from the first file in sample table######
  figures <- eventReactive(input$preview, {

    gen_fig_wrapper(react_val$config,  react_val$metadata, snv_format = input$snv_format, avail(), sample_table(), to_analyse = 1, adjust = input$adjust_in, arm_lvl = input$arm_lvl, estimate_lab = input$estimate_lab,
                    refDataExp, keepSNP, par_reg, centr_ref, weight_table, model_gender, model_dipl, model_alt, model_noSNV, chrs,
                    diploid_standard, scaleCols, dpRatioChrEdge)

    chr_figs <- file.path(react_val$config["out_dir"], sample_table()[1, 1], paste0("chromosome_", c(1:22, "X"), ".png"))
    main_fig <- file.path(react_val$config["out_dir"], sample_table()[1, 1], paste0(sample_table()[1, 1], "_CNV_main_fig.png"))

    figures <- list(chr_figs = chr_figs, main_fig = main_fig)

    def_table <- read.table(file = paste0(react_val$config["out_dir"], "/", "manual_an_table.tsv"), stringsAsFactors = FALSE, sep = "\t")
    react_val$man_table <- def_table
    react_val$def_table <- def_table
    react_val$check <- TRUE

    return(figures)
  })

  # Render main fig image for the first sample analysis
  output$main_fig_prev <- renderImage({

    width  <- session$clientData$output_myImage_width
    height <- session$clientData$output_myImage_height

    list(src = figures()$main_fig,
         contentType = "image/png",
         width = "100%",
         height = "auto"
    )
  }, deleteFile = FALSE)

  # Render chromosome figure based on the selected chromosome
  output$chr_fig_prev <- renderImage({
    list(src = figures()$chr_fig[chromosome()],
         contentType = "image/png",
         width = "100%",
         height = "auto"
    )
  }, deleteFile = FALSE)

  # Reactive value for keeping track of selected chromosome
  chromosome <- eventReactive(c(input$next_butt_chr_prev, input$prev_butt_chr_prev), {

    sel_chrom = input$next_butt_chr_prev - input$prev_butt_chr_prev

    if (sel_chrom >= 0) {
      return(sel_chrom %% 23 + 1)
    } else {
      return(23 - (abs(sel_chrom) %% 23) +1)
    }

  })

  # Create reactiveValues variable in order to load a estimation table
  react_val <- reactiveValues()
  # Default check reactive value as FALSE
  react_val$check <- FALSE

  ####Analyze all samples and save figures####
  observeEvent(input$analyze, {

    if (!is.null(react_val$config) & !is.null(react_val$metadata)) {

      if (!is.null(sample_table())) {

        gen_fig_wrapper(react_val$config,  react_val$metadata, snv_format = input$snv_format, avail(), sample_table(), to_analyse = nrow( react_val$metadata), adjust = input$adjust_in, arm_lvl = input$arm_lvl, estimate_lab = input$estimate_lab,
                        refDataExp, keepSNP, par_reg, centr_ref, weight_table, model_gender, model_dipl, model_alt, model_noSNV, chrs,
                        diploid_standard, scaleCols, dpRatioChrEdge)

        def_table <- read.table(file = paste0(react_val$config["out_dir"], "/", "manual_an_table.tsv"), stringsAsFactors = FALSE, sep = "\t")
        react_val$man_table <- def_table
        react_val$def_table <- def_table
        react_val$check <- TRUE

      }

    }

  })


  #render select button for scrolling through jpg files
  output$figure_select <- renderUI({

    react_val$config["out_dir"]

    if (react_val$check == FALSE) return(NULL)


    selectInput("sel_sample", "Select sample to visualize",
                choices = react_val$man_table$sample, selected = react_val$man_table$sample[1])
  })

  #read default estimation table in case config file is changed
  observeEvent(react_val$config, {

    table_exists <- file.exists(paste0(react_val$config["out_dir"], "/", "estimation_table.tsv"))

    if (table_exists == TRUE) {
      est_def <- read.table(file = paste0(react_val$config["out_dir"], "/", "estimation_table.tsv"), stringsAsFactors = FALSE, sep = "\t")
      est_def <- cbind(est_def, status = "not checked", comments = "none")

      react_val$def_table <- est_def
    } else {
      NULL
    }

  })

  #read estimation table for manual curration in case config file is changed
  observeEvent(react_val$config, {

    table_exists <- file.exists(paste0(react_val$config["out_dir"], "/", "manual_an_table.tsv"))

    if (table_exists == TRUE) {
      est_man <- read.table(file = paste0(react_val$config["out_dir"], "/", "manual_an_table.tsv"), stringsAsFactors = FALSE, sep = "\t")

      react_val$man_table <- est_man

    } else {
      NULL
    }

  })


  # list all figures in output directory
  fig_sam <- reactive({

    input$config
    react_val$check

    figs = sub(".*/", "", list.files(path = react_val$config["out_dir"], pattern = "_CNV_main_fig.png", recursive = TRUE, full.names = FALSE))

    if(length(figs) == 0) return(NULL)

    fig_sam <- sub("_CNV_main_fig.png", "", figs)
    return(fig_sam)

  })




  # check if the figures and estimation table match
  observe({

    #check whether files are present
    if (is.null(fig_sam()) | is.null(react_val$def_table) | is.null(react_val$man_table) | is.null(sample_table())) return(FALSE)

    figures <- fig_sam()

    #the check is completed only if the samples in the given directory match with samples in estimation tables and and samples from the input

      if (all(figures %in% react_val$def_table$sample) & all(figures %in% react_val$man_table$sample) & all(figures %in% pull(sample_table(), 1))) {
        react_val$check <- TRUE
      } else {
        react_val$check <- FALSE
      }

  })

  #Hide tabs in default
  hideTab(inputId = "tabs", target = "Manual CNV analysis")
  hideTab(inputId = "tabs", target = "Export")

  observe({

    if (react_val$check == TRUE) {
      showTab(inputId = "tabs", target = "Manual CNV analysis")
      showTab(inputId = "tabs", target = "Export")
    } else {
      hideTab(inputId = "tabs", target = "Manual CNV analysis")
      hideTab(inputId = "tabs", target = "Export")

    }
  })


  #render image based on the select button####
  output$main_fig <- renderImage({
    list(src = paste0(react_val$config["out_dir"], "/", input$sel_sample, "/", input$sel_sample,  "_CNV_main_fig.png"),
         contentType = "image/png",
         width = "100%",
         height = "auto"
    )
  }, deleteFile = FALSE)


  #render choices for selectInput from generated chromosome figures
  chr_choices <- eventReactive(input$sel_sample, {

    chromosomes <- list.files(path = paste0(react_val$config["out_dir"], "/", input$sel_sample, "/"), pattern = "^chromosome_.*.png$")

    if (length(chromosomes) == 0) {
      return(NULL)
    } else if (length(chromosomes) > 0) {
      choices <-  factor(gsub("_|.png", " ", chromosomes), levels = paste0("chromosome ", c(1:22, "X"), " "))
      choices <- choices[order(choices)]
      return(choices)
    }
  })

  #pass it to output variable for conditional ui
  output$chr_choices <- reactive({
    chr_choices()
  })

  outputOptions(output, "chr_choices", suspendWhenHidden = FALSE)


  #render select input widget for chromosomes
  output$chr_sel <- renderUI({

    input$sel_sample

    if (is.null(chr_choices())) {
      return(NULL)
    } else {
      selectInput("chr_sel", "Select chromosome to view in detail", choices = chr_choices())
    }

  })


  #create a filename for selected chromosome
  image <- reactive({
    if(is.null(input$chr_sel)) return(NULL)
    image <- sub(" ", "_", sub(" $", ".png", input$chr_sel))
    return(image)
  })

  #render chromosome image
  output$chr_fig <- renderImage({

    list(src = paste0(react_val$config["out_dir"], "/", input$sel_sample, "/", image()),
         contentType = "image/png",
         width = "100%",
         height = "auto")

  }, deleteFile = FALSE)


  #Update the ui when sample is changed
  observe({
    if(is.null(input$sel_sample)) return(NULL)

    #prevent crashing when changing output directories and updating selectInput and other inputs
    cur_sel <- input$sel_sample

    if (!cur_sel %in% react_val$man_table$sample) {
      sel <- react_val$man_table$sample[1]
    } else {
      sel <- cur_sel
    }

    #update gender
    gender <- react_val$man_table$gender[react_val$man_table$sample == sel]
    output$gender_select <- renderUI({
      sexes <- c("male", "female")
      selectInput(inputId = "gender_select", choices = sexes, selected = gender, label = "Select gender of the sample")
    })

    #update chromosome number
    chromn <- react_val$man_table$chrom_n[react_val$man_table$sample == sel]
    output$chromn_text <- renderUI({
      textInput(inputId = "chromn_text", value = chromn, label = "Type in the number of chromosomes")
    })

    #update alterations
    alt <- react_val$man_table$alterations[react_val$man_table$sample == sel]
    output$alt_text <- renderUI({
      textInput(inputId = "alt_text", value = alt, label = "Type in the chromosomal alterations")
    })

    #update comments
    comments <- react_val$man_table$comments[react_val$man_table$sample == sel]
    output$comments_text <- renderUI({
      textInput(inputId = "comments_text", value = comments, label = "Comments")
    })

  })

  #render next button and previous button and connect them with selected sample
  output$next_butt <- renderUI({

    actionButton("next_butt", "Next", width = "100%")

  })

  observeEvent(input$next_butt, {

    samples <- unlist(fig_sam())

    selected <- which(input$sel_sample == samples)

    if (selected == length(samples)) {
      next_sam <- samples[1]
    } else {
      next_sam <- samples[selected + 1]
    }

    updateSelectInput(session, inputId = "sel_sample", selected = next_sam)


  })

  output$prev_butt <- renderUI({
    if(is.null(input$sel_sample)) return(NULL)

    actionButton("prev_butt", "Previous", width = "100%")
  })

  observeEvent(input$prev_butt, {

    samples <- unlist(react_val$man_table$sample)

    selected <- which(input$sel_sample == samples)

    if (selected == 1) {
      prev_sam <- samples[length(samples)]
    } else {
      prev_sam <- samples[selected -1]
    }

    updateSelectInput(session, inputId = "sel_sample", selected = prev_sam)


  })

  #render save button
  output$save_butt <- renderUI({
    if(is.null(react_val$man_table)) return(NULL)

    actionButton("save_changes", "Save", width = "100%")
  })


  #save changes in the sample when button is clicked
  observeEvent(input$save_changes, {

    #read in table for manual analysis
    man_an <- read.table(paste0(react_val$config["out_dir"], "/", "manual_an_table.tsv"), stringsAsFactors = FALSE, sep = "\t")

    #change values from text input
    sam_ind <- which(man_an$sample == input$sel_sample)
    man_an[sam_ind, "gender"] <- input$gender_select
    man_an[sam_ind, "chrom_n"] <- input$chromn_text
    man_an[sam_ind, "alterations"] <- input$alt_text
    man_an[sam_ind, "comments"] <- input$comments_text
    man_an[sam_ind, "status"] <- "checked"

    write.table(man_an, file = paste0(react_val$config["out_dir"], "/", "manual_an_table.tsv"), sep = "\t", quote = FALSE)
    react_val$man_table <- man_an

  })

  #keep track of which files have and havn't been checked and of unsaved changes
  observeEvent(c(input$gender_select,
                 input$chromn_text,
                 input$alt_text,
                 input$comments_text,
                 input$save_changes
  ), {
    #do not run the code if input values are null
    if (is.null(input$gender_select) |
        is.null(input$chromn_text) |
        is.null(input$alt_text) |
        is.null(input$comments_text)
    ) {
      return(NULL)
    } else {
      #check whetter any text input changed
      if(input$gender_select != as.character(react_val$man_table$gender[react_val$man_table$sample == input$sel_sample]) |
         input$chromn_text != as.character(react_val$man_table$chrom_n[react_val$man_table$sample == input$sel_sample]) |
         input$alt_text != as.character(react_val$man_table$alterations[react_val$man_table$sample == input$sel_sample]) |
         input$comments_text != as.character(react_val$man_table$comments[react_val$man_table$sample == input$sel_sample])
      ) {

        status <- "unsaved changes"

      } else {
        status <- as.character(react_val$man_table$status[react_val$man_table$sample == input$sel_sample])
      }


      output$status <- renderText({
        if (status == "checked") {
          paste0('<p style="text-align:right;"><font color="green">', "<b>", status)
        } else if (status == "not checked") {
          paste0('<p style="text-align:right;"><font color="red">', "<b>", status)
        } else if (status == "unsaved changes") {
          paste0('<p style="text-align:right;"><font color="blue">', "<b>", status)
        }

      })
    }
  })

  #next and previous buttons for chromosome selection
  output$next_butt_chr <- renderUI({
    if (is.null(input$chr_sel)) return(NULL)

    actionButton("next_chr", "Next", width = "100%")

  })

  observeEvent(input$next_chr, {

    choices <- as.character(chr_choices())

    selected <- which(input$chr_sel == choices)

    if (selected == length(choices)) {
      next_chr <- choices[1]
    } else {
      next_chr <- choices[(selected + 1)]
    }

    updateSelectInput(session, inputId = "chr_sel", selected = next_chr)
  })


  output$prev_butt_chr <- renderUI({
    if (is.null(input$chr_sel)) return(NULL)

    actionButton("prev_chr", "Previous", width = "100%")

  })

  observeEvent(input$prev_chr, {

    choices <- as.character(chr_choices())

    selected <- which(input$chr_sel == choices)

    if (selected == 1) {
      prev_chr <- choices[length(choices)]
    } else {
      prev_chr <- choices[(selected - 1)]
    }

    updateSelectInput(session, inputId = "chr_sel", selected = prev_chr)
  })

  #render default ui button and allow to get back to default estimation
  output$default <- renderUI({
    if (is.null(react_val$def_table)) return(NULL)

    actionButton("default", "Default estimate")
  })

  observeEvent(input$default,{

    #update gender
    def_gender <- react_val$def_table$gender[react_val$def_table$sample == input$sel_sample]
    updateTextInput(session, inputId = "gender_select", value = def_gender)
    #update chromosome number
    def_chromn <- react_val$def_table$chrom_n[react_val$def_table$sample == input$sel_sample]
    updateTextInput(session, inputId = "chromn_text", value = def_chromn)
    #update alterations
    def_alt <- react_val$def_table$alterations[react_val$def_table$sample == input$sel_sample]
    updateTextInput(session, inputId = "alt_text", value = def_alt)
    #update comments
    def_comments <- react_val$def_table$comments[react_val$def_table$sample == input$sel_sample]
    updateTextInput(session, inputId = "comments_text", value = def_comments)
  })


  #render Text to represent which sample in row is being analyzed
  output$sample_num <- renderText({

    samples = unlist(react_val$man_table$sample)
    num = which(samples == input$sel_sample)

    paste0("<font size='5px'><p align='center'>Sample <b>", num, "</b> out of <b>", length(samples), "</b>")
  })

  #render checkbox user interface for export tab
  output$columns <- renderUI({

    variables = colnames(react_val$man_table)[-1]
    checkboxGroupInput("columns", "Choose columns to export", choices = variables, selected = variables)

  })

  #render radio buttons for export format
  output$format <- renderUI({

    radioButtons("format", "Choose format to export the file in", choices = c(".csv", ".tsv"), selected = ".csv")
  })


  #render preview for table to export
  output$prev_tab <- renderTable({

    to_show <- react_val$man_table %>% select("sample", input$columns) %>% head(20)
    return(to_show)

  })

  #Save the manually currated table after user selected a directory
  shinyDirChoose(input, "export", roots = volumes, session = session)

  observeEvent(input$export, {

     if (!is.integer(input$export) & length(input$export) > 1) {

      table <- react_val$man_table %>% select("sample", input$columns)
      exp_dir <- parseDirPath(volumes, input$export)

      if (input$format == ".csv") {
        write.csv(x = table, file = paste0(exp_dir, "/RNAseqCNA_an_table.csv"), row.names = FALSE, quote = FALSE)
      }

      if (input$format == ".tsv") {
        write.table(x = table, file = paste0(exp_dir, "/RNAseqCNA_an_table.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
      }
    }

  })


  #create a vector of input alterations
  alt_vec <- reactive({
    str_split(input$alt_text, ",")[[1]] %>% str_trim
  })

  #parse the alt_text text input and check whether the format is acceptable
  acceptable <- eventReactive(input$alt_text, {
    if (any(str_detect(string = alt_vec(), pattern = "^([1-9]|^1[0-9]|^2[0-2]|^X){1}[p-q]{0,1}(\\+|-)$|^none$") == FALSE)) {
      return(FALSE)
    } else {
      return(TRUE)
    }

  })

  #update chromosome number base on the alteration input
  observeEvent(input$alt_text, {


    plus = sum(str_count(alt_vec(), "^([1-9]|^1[0-9]|^2[0-2]|^X){1}\\+$"))
    minus = sum(str_count(alt_vec(), "^([1-9]|^1[0-9]|^2[0-2]|^X){1}-$"))

    updateTextInput(session, "chromn_text", value = (46 + plus - minus))

  })

  #render warning message if the alt_text is not in acceptable format
  output$war_message <- renderText({

    if (acceptable() == FALSE) {
      '<font color="red">This alteration format is not supported. Chromosome number will not be adjusted'
    } else {
      NULL
    }

  })

  #render prefered format hint for alteration input
  output$examp <- renderText({

    if (!is.null(input$alt_text)) {
      '<font size="2px"><p style="text-align:left;">Prefered format: 1+,3-,6p+,Xq-'
    }
  })

  #render button for choosing mock output directory
  shinyDirChoose(input, "dir_button", roots = volumes, session = session)

 #perform mock analysis when mock output directory is selected
  observeEvent(input$dir_button, {

    if (!is.integer(input$dir_button)) {

      source(system.file(package = "RNAseqCNVapp", "inst/extdata/config_mock.txt"))
      react_val$config <- c(count_dir = count_dir, snv_dir = snv_dir, out_dir = parseDirPath(volumes, input$dir_button))
      react_val$metadata <- fread(system.file(package = "RNAseqCNVapp", "inst/extdata/metadata_mock"), header = FALSE)
      sample_table <-  react_val$metadata %>% mutate(count_path = file.path(react_val$config["count_dir"], pull(., 2)), snv_path = file.path(react_val$config["snv_dir"], pull(., 3)))

          gen_fig_wrapper(react_val$config,  react_val$metadata, snv_format = input$snv_format, avail = "all_present", sample_table = sample_table, to_analyse = nrow( react_val$metadata), adjust = input$adjust_in, arm_lvl = input$arm_lvl, estimate_lab = input$estimate_lab,
                          refDataExp, keepSNP, par_reg, centr_ref, weight_table, model_gender, model_dipl, model_alt, model_noSNV, chrs,
                          diploid_standard, scaleCols, dpRatioChrEdge)

          def_table <- read.table(file = paste0(react_val$config["out_dir"], "/", "manual_an_table.tsv"), stringsAsFactors = FALSE, sep = "\t")
          react_val$man_table <- def_table
          react_val$def_table <- def_table
          react_val$check <- TRUE
    }
  })

}
