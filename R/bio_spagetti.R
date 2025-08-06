#' bio_spagetti
#'
#' @description This function generates a spagetti plot of biomarker readouts.
#' The output is stored in 1 RTF file under REPORT/OUTPUT/subFolder under eWISE environment,
#' or current working path/subFolder under local environment.
#'
#' @note In the input dataset indt, in addition to the required parameters, there must also be the variable USUBJID.
#'
#' @param indt Required, a R dataset name which is used to create the spagetti, such as adbm.
#' @param pageV A character string representing the page variable name, such as 'PARAM'.
#' Each value of page variable will have a spagetti on one page.
#' The default value is NULL, which will create only 1 spagetti.
#' @param xV Required, a character string representing the X variable name of spagetti, such as 'VISIT'.
#' @param xOrdV A character string representing the sorting variable name for X axis, such as 'VISITNUM'.
#' The default value is NULL, it will be displayed in the order of xV.
#' @param yV Required, a character string representing the Y variable name of spagetti, such as 'AVAL'.
#' @param columnV a character string representing the column variable name for facet_grid, such as 'TRT01A'.
#' The default value is NULL.
#' @param columnOrdV a character string representing sorting variable name for columnV, such as 'TRT01AN'.
#' The default value is NULL, it will be displayed in the order of columnV.
#' @param rowV a character string representing the column variable name for facet_grid, such as 'AGEGR1'.
#' The default value is NULL.
#' @param rowOrdV a character string representing sorting variable name for rowV, such as 'AGEGR1N'.
#' The default value is NULL, it will be displayed in the order of rowV
#' @param gmColor A character vector representing the color for geometric mean by columnV,
#' so the number of gmColor should be equal to the number of different value of columnV.
#' The default value is NULL, it will color randomly.
#' @param tit A character string representing the title of output, the default value is 'Biomarker Spagetti'.
#' @param footn A character string representing the nth footnote displayed as footnote, n = 1, ..., 8,
#' the default value is ''.
#' @param pgName A character string representing the program name displayed in the last footnote,
#' the default value is 'bio_spagetti'.
#' @param subFolder A character string representing sub folder name where RTF file locates, such as 'Primary', 'CMI', ...
#' the default value is NULL.
#'
#' @examples
#' \dontrun{
#'
#' adbm <- readRDS(system.file("extdata", "example.data.rds", "VacBioTLF")) %>%
#' filter(!is.na(AVAL)) %>%
#' mutate(VIS = str_split_i(VISIT, '/', 2),
#' sexn = ifelse(SEX == 'F', 2, 1))
#'
#' bio_spagetti(
#' indt = adbm,
#' pageV = 'PARAM',
#' xV = 'VIS',
#' yV = 'AVAL',
#' columnV = 'TRT01A',
#' columnOrdV = 'TRT01AN',
#' rowV = 'SEX',
#' rowOrdV = 'sexn',
#' gmColor = c('red', 'yellow')
#' )
#'
#' }
#'
#' @export
bio_spagetti <- function(indt,
                         pageV = NULL,
                         xV,
                         xOrdV = NULL,
                         yV,
                         columnV = NULL,
                         columnOrdV = NULL,
                         rowV = NULL,
                         rowOrdV = NULL,
                         gmColor = NULL,
                         tit = 'Biomarker Spagetti with GM and 95% CI',
                         foot1 = '',
                         foot2 = '',
                         foot3 = '',
                         foot4 = '',
                         foot5 = '',
                         foot6 = '',
                         foot7 = '',
                         foot8 = '',
                         pgName = 'bio_spagetti',
                         subFolder = NULL) {

  # preprocessing
  if (!is.null(columnV)){
    columnOrdV <- NULL
    gmColor <- NULL
  }

  if (!is.null(rowV)){
    rowOrdV <- NULL
  }


  # prepare color for GM group
  if (!is.null(gmColor) & !is.null(columnV)){

    tempColorDT <- indt %>%
      ungroup() %>%
      distinct(!!!rlang::syms(c(columnV, columnOrdV))) %>%
      arrange(!!!rlang::syms(c(columnOrdV, columnV)))


    # in case the number of color is not equal to the number of treatment
    if (length(gmColor) < nrow(tempColorDT) ){
      gmColor <- rep(gmColor, ceiling(nrow(tempColorDT)/length(gmColor)))[1:nrow(tempColorDT)]
    }else if (length(gmColor) > nrow(tempColorDT)) {
      gmColor <- gmColor[1:nrow(tempColorDT)]
    }

    # name gmColor
    names(gmColor) <- tempColorDT[[1]]
  }



  if (is.null(pageV)){

    temp2 <- indt %>%
      ungroup() %>%
      mutate(self_page = '')
  }else{

    temp2 <- indt %>%
      ungroup() %>%
      mutate(self_page = get(pageV))
  }



  final <- temp2 %>%
    ungroup() %>%
    filter(!is.na(get(yV))) %>%
    mutate(self_x = get(xV),
           self_y = get(yV),
           self_log = log10(get(yV))
    ) %>%
    select(USUBJID, starts_with('self_'), !!!rlang::syms(c(columnOrdV, columnV, rowOrdV, rowV)))



  # define the display order for x, column and row variables
  display_order <- function (inVar, inOrdVar){

    if (!is.null(inVar)){

      tempOrdDT <- final %>%
        distinct(!!!rlang::syms(c(inVar, inOrdVar))) %>%
        arrange(!!!rlang::syms(c(inOrdVar, inVar)))

      tempOrdDTLevel <- tempOrdDT[[inVar]]
    }
  }

  xOrdDTLevel <- display_order('self_x', xOrdV)
  columnOrdDTLevel <- display_order(columnV, columnOrdV)
  rowOrdDTLevel <- display_order(rowV, rowOrdV)

  final$self_x <- factor(final$self_x, levels = xOrdDTLevel)

  if(!is.null(columnV)){

    final[[columnV]] <- factor(final[[columnV]], levels = columnOrdDTLevel)
  }

  if(!is.null(rowV)){
    final[[rowV]] <- factor(final[[rowV]], levels = rowOrdDTLevel)
  }



  # output to RTF
  footVars <- unlist(mget(paste0("foot", 1:8)))

  outT <- 'rtf'
  dtsource <- 'adbm'
  heightP <- 5.4 - 0.2*sum(footVars != '')
  widthP <- 9
  tits <- tit


  # create folder to store the RTF file and dataset
  if (!exists("REPO")){
    REPO <- here::here()

    W_STUDY <- ''
  }

  if (!exists("REPD")){
    REPD <- here::here()
  }


  if (is.null(subFolder)){

    fname <- paste0(pgName, '.', outT)

  }else{

    folder_path <- file.path(REPO, subFolder)
    if (!file.exists(folder_path)) {
      dir.create(folder_path, recursive = TRUE)
    }

    fname <- paste0(subFolder, '/', pgName, '.', outT)

  }




  # footnotes
  foot9 <- paste0('Study: ', W_STUDY,' Program: ', tolower(pgName), '.R Datasets=', dtsource,
                  ' Output:/', str_replace_all(REPO,'~/wise/',''), '/', fname,
                  ' (', format(Sys.time(),'%d%b%Y %H:%M'),')')


  # Create temporary path
  tmp <- file.path(REPO, fname)


  # Add plot to report
  rpt <- create_report(tmp, output_type = toupper(outT), font = 'Times', font_size = 9) %>%
    set_margins(top = 1, bottom = 0.8) %>%
    options_fixed(font_size = 9)




  # circulation vars
  PAGE.Vector <- unique(final$self_page)

  allTmpPic <- NULL
  for (p in 1:length(PAGE.Vector)){

    tempDT <- final %>%
      filter(self_page == PAGE.Vector[p])


    # calculate GM and 95CI%
    # calculate 95% CI based on t distribution
    alpha <- 0.05
    tempDT2 <- tempDT %>%
      group_by(self_page, self_x, !!!rlang::syms(c(columnV, rowV))) %>%
      summarise(m = n(), mean = mean(self_log), sd = sd(self_log), df = m-1)

    quanT <- qt(1-alpha/2, tempDT2$df)
    tempDT2$quanT <- quanT

    tempDT3 <- tempDT2 %>%
      mutate(ci.lower = 10**(mean-quanT*sd/sqrt(m)),
             ci.upper = 10**(mean+quanT*sd/sqrt(m)),
             self_y = 10**(mean),
             USUBJID = 'MEAN') %>%
      select(self_page, self_x, !!!rlang::syms(c(columnV, rowV)),
             USUBJID, ci.lower, ci.upper, self_y)



    tempP <- ggplot(NULL, aes(self_x, self_y) ) +
      geom_line(data = tempDT, aes(group = USUBJID), color = '#636363', alpha = 0.5) +
      geom_point(data = tempDT, size = 1, alpha = 0.5)


    # Geometric Mean and CI
    if (!is.null(columnV)){

      tempP <- tempP +
        geom_line(data = tempDT3, aes(self_x, self_y, group = 1, color = get(columnV))) +
        geom_errorbar(data = tempDT3, aes(ymin  = ci.lower, ymax = ci.upper, color = get(columnV)), width = 0.2)


      if (!is.null(gmColor)){

        tempP <- tempP +
          scale_color_manual(values = gmColor)
      }

    }else {

      tempP <- tempP +
        geom_line(data = tempDT3, aes(self_x, self_y, group = 1)) +
        geom_errorbar(data = tempDT3, aes(ymin  = ci.lower, ymax = ci.upper), width = 0.2)
    }




    tempP <- tempP +
      scale_y_continuous(trans='log10',
                         breaks=scales::trans_breaks('log10', function(x) 10^x),
                         labels=scales::trans_format('log10', scales::math_format(10^.x))) +
      labs(title = paste0(PAGE.Vector[p]),
           x = '',
           y = paste('GM and 95% CI'),
           color = '') +
      theme_light() +
      theme(text = element_text(size = 9, family='serif'),
            plot.title = element_text(size = 10),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            strip.background = element_rect(fill = 'white',colour = 'white'),
            strip.text.x = element_text(color = 'black', size = 6),
            strip.text.y = element_text(color = 'black', size = 8),
            legend.position = 'none',
            legend.title = element_text(size = 8),
            legend.text = element_text(size = 7))


    if (!is.null(columnV) & !is.null(rowV)){

      tempP <- tempP +
        facet_grid( get(rowV) ~ get(columnV) )

    }else if (!is.null(rowV)){

      tempP <- tempP +
        facet_grid(get(rowV) ~ . )

    }else if (!is.null(columnV)){

      tempP <- tempP +
        facet_grid(. ~ get(columnV) )
    }



    if (is.null(subFolder)){
      tmpPic <- file.path(REPO, paste0(pgName, p, '.png'))
    }else{
      tmpPic <- file.path(REPO, paste0(subFolder, '/', pgName, p, '.png'))
    }

    ggsave(tmpPic,
           plot = tempP,
           # type='cairo',
           width = widthP,
           height = heightP,
           units = 'in')


    # Define plot object
    tempP <- create_plot(tmpPic,
                         height = heightP,
                         width = widthP,
                         borders = 'none')


    rpt <- rpt %>%
      add_content(tempP, page_break = T)


    rm(tempP)
    allTmpPic <- c(allTmpPic, tmpPic)
    assign(x = 'allTmpPic', value = allTmpPic, envir = .GlobalEnv)


  }


  footVars2 <- footVars[footVars != '']
  rpt <- rpt %>%
    title_header(' ', tit, right = c('Page [pg] of [tpg]',' '), borders = 'bottom', blank_row = 'below') %>%
    footnotes(paste0(c(footVars2, foot9), collapse = '\n'),
              borders = 'top',
              blank_row = 'none')

  # Write out report
  res <- write_report(rpt)


  # remove all PNG
  print(length(allTmpPic))
  print(allTmpPic)
  file.remove(allTmpPic)

}




