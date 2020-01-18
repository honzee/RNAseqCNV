####Server####
shinyAppServer <- function(input, output, session) {

  #change the size of possible upload
  options(shiny.maxRequestSize=30*1024^2)

  #do the directories from the input exist?
  config <- reactive({
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

      return(c(count_dir = count, snv_dir = snv, out_dir = out))

    } else {
      return("no_input")
    }
  })


  #render warning messages if directories don't exist
  output$mess_config <- renderText({
    if(any(config() %in% "not_found")) {

      message <- NULL

      for (i in 1:length(config())) {
        if (config()[i] == "not_found") {
          message <- paste0(message, '<br><font size="1px"><p style="text-align:left;"><font color="red">', names(config())[i], " was not found")
        }
      }

      return(message)

    } else {
      NULL
    }

  })


  #check metadata file for three columns and read it
  metadata <- reactive({

    if (is.null(input$metadata)) return("no_input")

    metadata_tab = fread(input$metadata$datapath, header = FALSE)

    if (ncol(metadata_tab) == 3) {
      return(metadata_tab)
    } else {
      return("incorrect_format")
    }
  })

  #render warning messages if metadata table does not have three columns
  output$mess_metadata <- renderText({

    if (all(metadata() == "no_input")) return(NULL)

    if(all(metadata() == "incorrect_format")) {

      return('<font size="1px"><p style="text-align:left;"><font color="red">Metadata table does not have three columns')

    }
  })


  #create a sample table
  sample_table <- reactive({

    if (all(metadata() == "no_input")) return(NULL)
    if (all(config() == "no_input")) return(NULL)


    if (any(config() == "not_found")) return(NULL)
    if (all(metadata() == "incorrect_format")) return(NULL)

    sample_table <- metadata()

    HTSeq_f = pull(sample_table, 2)
    snv_f = pull(sample_table, 3)
    #create paths to the files
    sample_table$count_path <- file.path(config()["count_dir"], HTSeq_f)
    sample_table$snv_path <-  file.path(config()["snv_dir"], snv_f)
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
  observeEvent(input$preview, {

    prev_fig <- gen_fig_wrapper(config(), metadata(), avail(), sample_table(), preview = TRUE, prev_chr = 1, adjust = input$adjust_in, arm_lvl = input$arm_lvl, estimate = input$estimate,
                                refDataExp, keepSNP, par_reg, centr_ref, weight_table, model_gender, model_dipl, model_alt, chrs,
                                base_matr, base_col, scaleCols, dpRatioChrEdge)

    if (!is.null(prev_fig)) {
      output$plot_arm <- renderPlot({
        prev_fig$gg_arm
      })

      output$arr <- renderPlot({
        prev_fig$fig
      }, height = (function() {
        session$clientData$output_arr_width/2
      }) )
    }

  })

  ####Analyze all samples and save figures####
  observeEvent(input$analyze, {

    gen_fig_wrapper(config(), metadata(), avail(), sample_table(), preview = FALSE, prev_chr = 1, adjust = input$adjust_in, arm_lvl = input$arm_lvl, estimate = input$estimate,
                    refDataExp, keepSNP, par_reg, centr_ref, weight_table, model_gender, model_dipl, model_alt, chrs,
                    base_matr, base_col, scaleCols = scaleCols_DES_norm, dpRatioChrEdge)


    #Reload the directory for the the analysis tab to refresh
    updateTextInput(session, inputId = "out_dir", value = "")
    updateTextInput(session, inputId = "out_dir", value = config()["out_dir"])

  })


  #render select button for scrolling through jpg files
  output$figure_select <- renderUI({

    config()["out_dir"]

    if (check() == FALSE) return(NULL)


    selectInput("sel_sample", "Select sample to visualize",
                choices = est_table_man()$sample, selected = est_table_man()$sample[1])
  })

  #read default estimation table if possible
  est_table_def <- reactive({

    table_exists <- file.exists(paste0(config()["out_dir"], "/", "estimation_table.tsv"))

    if (table_exists == TRUE) {
      est_def <- read.table(file = paste0(config()["out_dir"], "/", "estimation_table.tsv"), stringsAsFactors = FALSE)
      est_def <- cbind(est_def, status = "not checked", comments = "none")
      return(est_def)
    } else {
      NULL
    }

  })

  #read estimation table for manual corrections
  est_table_man <- reactive({

    input$save_changes

    table_exists <- file.exists(paste0(config()["out_dir"], "/", "manual_an_table.tsv"))

    if (table_exists == TRUE) {
      est_man <- read.table(file = paste0(config()["out_dir"], "/", "manual_an_table.tsv"), stringsAsFactors = FALSE)
      return(est_man)
    } else {
      NULL
    }

  })


  #list all figures in output directory
  fig_sam <- reactive({

    figs = sub(".*/", "", list.files(path = config()["out_dir"], pattern = "_CNV_main_fig.png", recursive = TRUE, full.names = FALSE))

    if(length(figs) == 0) return(NULL)

    fig_sam <- sub("_CNV_main_fig.png", "", figs)
    return(fig_sam)

  })


  #check if the figures and estimation table match
  check <- reactive({

    #check whether files are present
    if (is.null(fig_sam()) | is.null(est_table_def()) | is.null(est_table_man())) return(FALSE)

    figures <- fig_sam()

    #the check is completed only if the samples in the table and figure samples are identical

    if (length(figures) == nrow(est_table_def()) & length(figures) == nrow(est_table_man())) {

      if (all(figures %in% est_table_def()$sample) & all(figures %in% est_table_man()$sample)) {
        return(TRUE)
      } else {
        return(FALSE)
      }
    } else {
      return(FALSE)
    }

  })


  #show/hide analysis and export tab based upon the table for manual curration and check value
  observe({

    if(is.null(est_table_man()) | check() == FALSE) {
      hideTab(inputId = "tabs", target = "Analysis")
      hideTab(inputId = "tabs", target = "Export")
    } else {
      showTab(inputId = "tabs", target = "Analysis")
      showTab(inputId = "tabs", target = "Export")

    }

  })

  #render image based on the select button####
  output$main_fig <- renderImage({
    list(src = paste0(config()["out_dir"], "/", input$sel_sample, "/", input$sel_sample,  "_CNV_main_fig.png"),
         contentType = "image/png",
         width = "100%",
         height = "auto"
    )
  }, deleteFile = FALSE)


  #render choices for selectInput from generated chromosome figures
  chr_choices <- eventReactive(input$sel_sample, {

    chromosomes <- list.files(path = paste0(config()["out_dir"], "/", input$sel_sample, "/"), pattern = "^chromosome_.*.png$")

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

    list(src = paste0(config()["out_dir"], "/", input$sel_sample, "/", image()),
         contentType = "image/png",
         width = "100%",
         height = "auto")

  }, deleteFile = FALSE)


  #Update the ui when sample is changed
  observe({
    if(is.null(input$sel_sample)) return(NULL)

    #prevent crashing when changing output directories and updating selectInput and other inputs
    cur_sel <- input$sel_sample

    if (!cur_sel %in% est_table_man()$sample) {
      sel <- est_table_man()$sample[1]
    } else {
      sel <- cur_sel
    }

    #update selectinput
    updateSelectInput(session, "sel_sample", choices = fig_sam(), selected = sel)

    #update type
    type <- est_table_man()$type[est_table_man()$sample == sel]
    output$type_text <- renderUI({
      textInput(inputId = "type_text", value = type, label = "Type in karyotypic subtype")
    })

    #update gender
    gender <- est_table_man()$gender[est_table_man()$sample == sel]
    output$gender_text <- renderUI({
      textInput(inputId = "gender_text", value = gender, label = "Type in the gender of the sample")
    })

    #update chromosome number
    chromn <- est_table_man()$chrom_n[est_table_man()$sample == sel]
    output$chromn_text <- renderUI({
      textInput(inputId = "chromn_text", value = chromn, label = "Type in the number of chromosomes")
    })

    #update alterations
    alt <- est_table_man()$alterations[est_table_man()$sample == sel]
    output$alt_text <- renderUI({
      textInput(inputId = "alt_text", value = alt, label = "Type in the chromosomal alterations")
    })

    #update comments
    comments <- est_table_man()$comments[est_table_man()$sample == sel]
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

    samples <- unlist(est_table_man()$sample)

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
    if(is.null(est_table_man())) return(NULL)

    actionButton("save_changes", "Save", width = "100%")
  })


  #save changes in the sample when button is clicked
  observeEvent(input$save_changes, {

    #read in table for manual analysis
    man_an <- read.table(paste0(config()["out_dir"], "/", "manual_an_table.tsv"), stringsAsFactors = FALSE)

    #change values from text input
    sam_ind <- which(man_an$sample == input$sel_sample)
    man_an[sam_ind, "type"] <- input$type_text
    man_an[sam_ind, "gender"] <- input$gender_text
    man_an[sam_ind, "chrom_n"] <- input$chromn_text
    man_an[sam_ind, "alterations"] <- input$alt_text
    man_an[sam_ind, "comments"] <- input$comments_text
    man_an[sam_ind, "status"] <- "checked"

    write.table(man_an, file = paste0(config()["out_dir"], "/", "manual_an_table.tsv"), sep = "\t")

  })

  #keep track of which files have and havn't been checked and of unsaved changes
  observeEvent(c(input$type_text,
                 input$gender_text,
                 input$chromn_text,
                 input$alt_text,
                 input$comments_text,
                 input$save_changes
  ), {
    #do not run the code if input values are null
    if (is.null(input$type_text) |
        is.null(input$gender_text) |
        is.null(input$chromn_text) |
        is.null(input$alt_text) |
        is.null(input$comments_text)
    ) {
      return(NULL)
    } else {
      #check whether any text input changed
      if(input$type_text != as.character(est_table_man()$type[est_table_man()$sample == input$sel_sample]) |
         input$gender_text != as.character(est_table_man()$gender[est_table_man()$sample == input$sel_sample]) |
         input$chromn_text != as.character(est_table_man()$chrom_n[est_table_man()$sample == input$sel_sample]) |
         input$alt_text != as.character(est_table_man()$alterations[est_table_man()$sample == input$sel_sample]) |
         input$comments_text != as.character(est_table_man()$comments[est_table_man()$sample == input$sel_sample])
      ) {

        status <- "unsaved changes"

      } else {
        status <- as.character(est_table_man()$status[est_table_man()$sample == input$sel_sample])
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
    if (is.null(est_table_def())) return(NULL)

    actionButton("default", "Default estimate")
  })

  observeEvent(input$default,{

    #update type
    def_type <- est_table_def()$type[est_table_def()$sample == input$sel_sample]
    updateTextInput(session, inputId = "type_text", value = def_type)
    #update gender
    def_gender <- est_table_def()$gender[est_table_def()$sample == input$sel_sample]
    updateTextInput(session, inputId = "gender_text", value = def_gender)
    #update chromosome number
    def_chromn <- est_table_def()$chrom_n[est_table_def()$sample == input$sel_sample]
    updateTextInput(session, inputId = "chromn_text", value = def_chromn)
    #update alterations
    def_alt <- est_table_def()$alterations[est_table_def()$sample == input$sel_sample]
    updateTextInput(session, inputId = "alt_text", value = def_alt)
    #update comments
    def_comments <- est_table_def()$comments[est_table_def()$sample == input$sel_sample]
    updateTextInput(session, inputId = "comments_text", value = def_comments)
  })


  #render Text to represent which sample in row is being analyzed
  output$sample_num <- renderText({

    samples = unlist(est_table_man()$sample)
    num = which(samples == input$sel_sample)

    paste0("<font size='5px'><p align='center'>Sample <b>", num, "</b> out of <b>", length(samples), "</b>")
  })

  #render checkbox user interface for export tab
  output$columns <- renderUI({

    variables = colnames(est_table_man())[-1]
    checkboxGroupInput("columns", "Choose columns to export", choices = variables, selected = variables)

  })

  #render textInput for export directory
  output$exp_out_dir <- renderUI({

    textInput("exp_out_dir", "Directory to export the file to (defaults to the output directory from the input tab)")
  })

  #check export directory
  check_exp_dir <- reactive({
    dir.exists(input$exp_out_dir) | input$exp_out_dir == ""
  })

  #render warning message if export directory does not exist
  output$mess_exp_out <- renderText({
    if (check_exp_dir() == FALSE) {
      return('<font size="1px"><p style="text-align:left;"><font color="red">Could not find this directory')
    }
  })


  #render radio buttons for export format
  output$format <- renderUI({

    radioButtons("format", "Choose format to export the file in", choices = c(".csv", ".tsv"), selected = ".csv")
  })

  #render the export button
  output$export <- renderUI({

    actionButton("export", "Export the table")
  })

  #render preview for table to export
  output$prev_tab <- renderTable({

    to_show <- est_table_man() %>% select("sample", input$columns) %>% head(20)
    return(to_show)

  })

  #export table when actionButton is clicked
  observeEvent(input$export, {

    if(check_exp_dir() == FALSE) return(showNotification("Could not find this directory", duration = 5, type = "message"))

    out_dir <- if (input$exp_out_dir == "") {
      config()["out_dir"]
    } else {
      input$exp_out_dir
    }

    table <- est_table_man() %>% select("sample", input$columns)

    #save a format based on the input from radio_buttons

    if (input$format == ".csv") {
      write.csv(x = table, file = paste0(out_dir, "/RNAseqCNA_an_table.csv"), row.names = FALSE)
    }

    if (input$format == ".tsv") {
      write.csv(x = table, file = paste0(out_dir, "/RNAseqCNA_an_table.tsv"), sep = "\t", row.names = FALSE)
    }

    output$export_mess <- renderText({
      if (file.exists(paste0(out_dir, "/RNAseqCNA_an_table", input$format))) {
        paste0('<font color="green">', "<b>", "saved")
      } else {
        paste0('<font color="red">', "<b>", "failed")
      }
    })

  })

  #hide the export message if any value changes
  observeEvent(c(input$columns,
                 input$exp_out_dir,
                 input$format), {

                   output$export_mess <- renderText({
                     return("")
                   })

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

}
