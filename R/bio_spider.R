#' bio_spider
#'
#' @description This function generates a spider plot of mean of biomarker readouts.
#' The output is stored in 1 RTF file under REPORT/OUTPUT/Primary under eWISE environment,
#' or current working path/Primary under local environment.
#'
#' @param indt Required, a R dataset name which is used to create the spider, such as adbm.
#' @param pageV A character string representing the page variable name, such as 'VISIT'.
#' Each value of page variable will have a spider on one page.
#' The default value is NULL, which will create only 1 spider.
#' @param dimV Required, a character string representing the dimension variable name of spider, such as 'PARAM'.
#' @param dimOrdV A character string representing the sorting variable name for dimension variable, such as 'PARAMN'.
#' The default value is NULL, it will be displayed in the order of dimV.
#' @param valueV Required, a character string representing the Value variable name of spider, such as 'AVAL'.
#' @param groupV a character string representing the group variable name of spider, such as 'TRT01A'.
#' The default value is NULL.
#' @param groupOrdV a character string representing the sorting order for groupColor.
#' The default value is NULL.
#' @param groupColor a character vector representing the color for group,
#' so the number of groupColor should be equal to the number of different value of groupV or groupOrdV.
#' The default value is NULL, it will color randomly.
#' @param refValue a floating number representing the reference circle of the spider, the default value is 0.
#' @param pValueDT a R dataset name contains the p-values of pairwise comparisons between groupV.
#' In addition to other necessary variables, variable estimates must be included in the pValueDT.
#' For example, dimV and groupV, or pageV (if any).
#' @param tit A character string representing the title of output, the default value is 'Biomarker Spider'.
#' @param footn A character string representing the nth footnote displayed as footnote, n = 1, ..., 8,
#' the default value is ''.
#' @param pgName A character string representing the program name displayed in the last footnote,
#' the default value is 'bio_spider'.
#' @param subFolder A character string representing sub folder name where RTF file locates, such as 'Primary', 'CMI', ...
#' the default value is NULL.
#'
#' @examples
#' \dontrun{
#'
#' adbm <- readRDS(here::here("data", "example.data.rds")) %>%
#' filter(!is.na(R2BASE) & VISITNUM > 10) %>%
#' mutate(VIS = str_split_i(VISIT, '/', 2),
#' sexn = ifelse(SEX == 'F', 2, 1))
#'
#' pvalue.data <- readRDS(here::here("data", "pvalue.data.rds"))
#'
#' bio_spider(
#' indt = adbm,
#' pageV = 'VIS',
#' dimV = 'PARAM',
#' dimOrdV = NULL,
#' valueV = 'R2BASE',
#' groupV = 'TRT01A',
#' groupOrdV = 'TRT01AN',
#' refValue = 1,
#' pValueDT = pvalue.data
#' )
#'
#' }  # Dimension Variable, Value Variable, Grouping Variable
#'
#' @export
bio_spider <- function(indt,
                       pageV = NULL,
                       dimV,
                       dimOrdV = NULL,
                       valueV,
                       groupV = NULL,
                       groupOrdV = NULL,
                       groupColor = c( "#FFEC8B", "#C4B045", "#A89222", "#8B7500", 'gray35', 'gray'),
                       refValue = NULL,
                       pValueDT = NULL,
                       tit = 'Spider plot of mean of biomarker readouts',
                       foot1 = '',
                       foot2 = '',
                       foot3 = '',
                       foot4 = '',
                       foot5 = '',
                       foot6 = '',
                       foot7 = '',
                       foot8 = '',
                       pgName = 'bio_spider',
                       subFolder = NULL) {


  if (is.null(groupV)){
    groupOrdV <- NULL
    pValueDT <- NULL
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




  # number of page variable value
  pageCntDT <- temp2 %>%
    distinct(self_page)

  pageCnt <- nrow(pageCntDT)

  # 5 common markers
  commonMarker <- temp2 %>%
    ungroup() %>%
    distinct(self_page, !!sym(dimV)) %>%
    group_by(!!sym(dimV)) %>%
    summarise(cnt = n()) %>%
    filter(cnt == pageCnt) %>%
    ungroup()


  # self_theta is the angle for dimension variable on the plot
  commonMarkerCnt <- length(unique(commonMarker[[dimV]]))


  final <- inner_join(temp2, commonMarker) %>%
    arrange(!!!syms(c(dimOrdV, dimV))) %>%
    group_by(!!!syms(c(dimOrdV, dimV))) %>%
    mutate(self_dimVid = cur_group_id(),
           self_theta = 2*pi/commonMarkerCnt * (self_dimVid-1) + pi/2) %>%
    ungroup()  %>%
    arrange(self_page)


  thetaDT <- final %>%
    distinct(self_theta, self_dimVid, !!sym(dimV))



  if (!is.null(pValueDT)){
    tempPvalue <- pValueDT  %>%
      left_join(thetaDT)

    if (is.null(pageV)){
        tempPvalue2 <- tempPvalue %>%
          ungroup() %>%
          mutate(self_page = '')
    }else{
        tempPvalue2 <- tempPvalue %>%
          ungroup() %>%
          mutate(self_page = get(pageV))
    }
  }






  # circulation vars
  PAGE.Vector <- unique(final$self_page)


  # Define the maximum and minimum values of dimension variable
  # Keep same maximum and minimum values for each spider plot
  maxminDT <- final %>%
    group_by(self_page, !!!syms(c(groupV, dimV))) %>%
    summarise(self_mean = mean(!!sym(valueV))) %>%
    pivot_wider(id_cols = c(self_page, all_of(groupV)), names_from = dimV, values_from = self_mean) %>%
    ungroup() %>%
    select(-self_page, -all_of(groupV)) %>%
    select_if(~ !any(is.na(.)))


  maxDim <- max(maxminDT, na.rm = T)
  minDim <- min(maxminDT, na.rm = T)



  if (!is.null(refValue)){
    # Create the max and min row
    maxDT <- maxminDT %>%
      slice_head() %>%
      mutate(across(everything(), ~ max(maxDim, refValue)))

    minDT <- maxDT %>%
      mutate(across(everything(), ~ min(minDim, refValue)))

    refCircle <- maxDT %>%
      mutate(across(where(is.numeric), ~ refValue))


    ySeq <- seq(min(minDim, refValue), max(maxDim, refValue), length.out = 5 )
    ySeqC <- formatC(ySeq, format = 'f',digits = 2)
  }else if (is.null(refValue)){
    # Create the max and min row
    maxDT <- maxminDT %>%
      slice_head() %>%
      mutate(across(everything(), ~ maxDim))

    minDT <- maxDT %>%
      mutate(across(everything(), ~ minDim))


    ySeq <- seq(minDim, maxDim, length.out = 5 )
    ySeqC <- formatC(ySeq, format = 'f',digits = 2)
  }


  if (!is.null(groupV)){
    # prepare color
    colorOrder <- final %>%
      distinct(!!!syms(c(groupOrdV, groupV))) %>%
      arrange(!!!syms(c(groupOrdV, groupV)))


    # in case the number of color is not equal to the number of treatment
    if (length(groupColor) < nrow(colorOrder) ){
      groupColor <- rep(groupColor, ceiling(nrow(colorOrder)/length(groupColor)))[1:nrow(colorOrder)]
    }else if (length(groupColor) > nrow(colorOrder)) {
      groupColor <- groupColor[1:nrow(colorOrder)]
    }

    names(groupColor) <- colorOrder[[groupV]]
  }



  # output to RTF
  footVars <- unlist(mget(paste0("foot", 1:8)))

  outT <- 'rtf'
  dtsource <- 'adbm'
  heightP <- 5 - 0.2*sum(footVars != '')
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




  allTmpPic <- NULL
  for (p in 1:length(PAGE.Vector)){

    tempDT <- final %>%
      filter(self_page == PAGE.Vector[p]) %>%
      ungroup() %>%
      arrange(!!!syms(c('self_page', groupV)))



    tempDT2 <- tempDT %>%
      group_by(self_page, !!!syms(c(groupV, dimV))) %>%
      summarise(self_mean = mean(!!sym(valueV))) %>%
      pivot_wider(id_cols = c(self_page, all_of(groupV)), names_from = dimV, values_from = self_mean) %>%
      ungroup() %>%
      arrange(!!!syms(c('self_page', groupV))) %>%
      select_if(~ !any(is.na(.))) %>%
      select(-self_page, -all_of(groupV))


    tempDT10 <- bind_rows(maxDT, minDT, tempDT2)


    if (is.null(groupV)){
      groupColor3 <- 'red'
      groupColor3Legend <- ''
      groupty <- 1
    }else if (!is.null(groupV)){
      # subset color for each page
      groupVector <- unique(tempDT[[groupV]])
      groupColor2 <- groupColor[groupVector]

      if (!is.null(refValue)){
        tempDT10 <- bind_rows(tempDT10, refCircle)

        groupColor3 <- c(groupColor2, 'red')
        groupColor3Legend <- c(groupVector, paste0('ref = ', refValue))
        groupty <- c(1: length(groupColor2), 1)
      }else if(is.null(refValue)){

        groupColor3 <- groupColor2
        groupColor3Legend <- groupVector
        groupty <- c(1: length(groupColor2))
      }
    }


    if (is.null(subFolder)){
      tmpPic <- file.path(REPO, paste0(pgName, p, '.png'))
    }else{
      tmpPic <- file.path(REPO, paste0(subFolder, '/', pgName, p, '.png'))
    }


    png(tmpPic,
        width = widthP,
        height = heightP,
        units = 'in',
        res = 800,
        pointsize = 3,
        type='cairo'
    )


    radarchart(tempDT10,
               axistype = 1,
               pcol = groupColor3, pfcol = scales::alpha(groupColor3, 0.1), plwd = 1.5, plty = groupty,
               cglcol = 'grey', cglty = 3, cglwd = 0.8,
               vlcex = 2,
               calcex = 3, caxislabels = ySeqC, axislabcol = "black"
    )
    title(main = paste0(PAGE.Vector[p]), cex.main = 2)
    legend(x = "right",
           legend = groupColor3Legend,
           col = groupColor3,
           text.col = "black", cex = 2, pt.cex = 2, pch = 5
    )



    # add significant mark
    if (!is.null(pValueDT)){
      tempPvalue3 <- tempPvalue2 %>%
        filter(self_page == PAGE.Vector[p]) %>%
        ungroup() %>%
        group_by(self_dimVid, !!sym(dimV)) %>%
        arrange(self_dimVid, !!!syms(c(dimV, groupV))) %>%
        mutate(self_group = !!sym(groupV),
               trtID =  row_number(),
               radius = 0.9,
               xPos = (radius + trtID/25) * cos(self_theta),
               yPos = (radius + trtID/25) * sin(self_theta),
               label = ifelse(estimate > 0, "+", "-") ,
               markColor = groupColor[self_group]
        )


      if (nrow(tempPvalue3) > 0){
        # add p value
        text(tempPvalue3$xPos,
             tempPvalue3$yPos,
             labels = tempPvalue3$label,
             pos = 1,
             cex = 5,
             col = tempPvalue3$markColor)
      }
    }


    dev.off()

    tempP5 <- create_plot(tmpPic,
                          height = heightP,
                          width = widthP,
                          borders = 'none')

    rpt <- rpt %>%
      add_content(tempP5, page_break = T)

    allTmpPic <- c(allTmpPic, tmpPic)
    assign(x = 'allTmpPic', value = allTmpPic, envir = .GlobalEnv)

  }


  footVars2 <- footVars[footVars != '']
  rpt <- rpt %>%
    title_header(' ', tits, right = c('Page [pg] of [tpg]',' '), borders = 'bottom', blank_row = 'below') %>%
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




