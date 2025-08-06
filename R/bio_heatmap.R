#' bio_heatmap
#'
#' @description This function generates a heatmap of biomarker readouts.
#' The output is stored in 1 RTF file under REPORT/OUTPUT/subFolder under eWISE environment,
#' or current working path/subFolder under local environment.
#'
#' @note In the input dataset indt, in addition to the required parameters, there must also be the variable USUBJID.
#'
#' @param indt Required, a R dataset name which is used to create the heatmap, such as adbm.
#' @param pageV Required, a character string representing the page variable name, such as 'AGEGR1'.
#' Each value of page variable will have a heatmap on one page.
#' @param columnV Required, a character vector representing the column variable name of heatmap,
#' such as c('TRT01A', 'VIS') or 'TRT01A'.
#' @param columnOrdV A character vector representing the sorting variable name for column,
#' such as c('TRT01AN', 'VISITNUM') or 'TRT01AN'. The default value is NULL, columnV is used in this case.
#' @param columnColor Required, a list representing the colors for column variable name, such as
#' list(c('orchid3', 'darkorchid1', 'darkorchid2', 'darkorchid3', 'darkorchid4','azure3'), c("#F0F0F0", "#9A9A9A", "#545454")) or
#' list(c('orchid3', 'darkorchid1', 'darkorchid2', 'darkorchid3', 'darkorchid4','azure3')).
#' The length of each element of columnColor should be equal to the number of different values of columnV respectively,
#' @param rowV Required, a character string representing the row variable name of heatmap, such as 'PARAM'.
#' @param analysisV Required, a character string representing the variable name which is used to create heatmap,
#' such as 'AVAL', 'R2BASE' or 'FOLREISE'.
#' @param analysisType A character string representing the type of heapmap.
#' There are two values available for selection 'GM' or 'Individual', the default value is 'GM'.
#' @param tit A character string representing the title of output, the default value is 'Biomarker Heatmap'.
#' @param footn A character string representing the nth footnote displayed as footnote, n = 1, ..., 8, the default value is ''.
#' @param pgName A character string representing the program name displayed in the last footnote,
#' the default value is 'bio_heatmap'.
#' @param subFolder A character string representing sub folder name where RTF file locates, such as 'Primary', 'CMI', ...
#' the default value is NULL.
#'
#' @examples
#' \dontrun{
#'
#' adbm <- readRDS(system.file("extdata", "example.data.rds", "VacBioTLF")) %>%
#' filter(!is.na(AVAL)) %>%
#' mutate(VIS = str_split_i(VISIT, '/', 2))
#'
#' bio_heatmap(
#' indt = adbm,
#' pageV = 'AGEGR1',
#' columnV = c('TRT01A', 'VIS'),
#' columnOrdV = c('TRT01AN', 'VISITNUM'),
#' columnColor = list(c('orchid3', 'darkorchid1', 'darkorchid2', 'darkorchid3', 'darkorchid4','azure3'),
#                                        c("#F0F0F0", "#9A9A9A", "#545454")),
#' rowV = 'PARAM',
#' analysisV = 'AVAL',
#' analysisType = 'GM'
#' )
#'
#' }
#'
#' @export
bio_heatmap <- function(indt,
                        pageV,
                        columnV,
                        columnOrdV = NULL,
                        columnColor,
                        rowV,
                        analysisV,
                        analysisType = 'GM',
                        tit = 'Biomarker Heatmap',
                        foot1 = '',
                        foot2 = '',
                        foot3 = '',
                        foot4 = '',
                        foot5 = '',
                        foot6 = '',
                        foot7 = '',
                        foot8 = '',
                        pgName = 'bio_heatmap',
                        subFolder = NULL) {


  if (is.null(columnOrdV)){
    columnOrdV <- columnV
  }


  # cluster type for heatmap
  clusterType <- c('Noncluster', 'ClusterRow', 'Cluster')
  clusterRow <- c(F, T, T)
  clusterColumn <- c(F, F, T)




  temp <- indt %>%
    ungroup() %>%
    filter(!is.na(get(analysisV))) %>%
    mutate(across(c(columnV, USUBJID), ~str_remove_all(., '_')),
           self_page = get(pageV),
           self_row = get(rowV),
           self_log = log10(get(analysisV)),
           self_column = case_when(analysisType == 'GM' ~ paste(!!!rlang::syms(columnV), sep = '_'),
                                   analysisType != 'GM' ~ paste(!!!rlang::syms(columnV), USUBJID, sep = '_'))
    ) %>%
    select(starts_with('self_'), !!!rlang::syms(columnOrdV))



  # column order
  columnOrd.data <- temp %>%
    distinct(!!!rlang::syms(columnOrdV), self_column) %>%
    arrange(!!!rlang::syms(columnOrdV), self_column) %>%
    mutate(self_columnOrd = row_number()) %>%
    select(-all_of(columnOrdV))


  temp2 <- left_join(temp, columnOrd.data)


  # columnV2 is used to get annoCol for heatmap
  if (analysisType == 'GM'){

    # calculate GM
    final <- temp2 %>%
      group_by(self_page, self_column, self_columnOrd, self_row) %>%
      summarise(self_res = 10**mean(self_log)) %>%
      ungroup()


    columnV2 <- columnV

    clusterTypeCnt <- length(clusterType)

  }else{

    final <- temp2 %>%
      mutate(self_res = self_log)

    columnV2 <- c(columnV, 'USUBJID')

    clusterTypeCnt <- length(clusterType) - 1
  }




  # name the color variable and list
  annAllColors <- list()
  for (c in 1:length(columnOrdV)){

    tempColorDT <- indt %>%
      distinct(!!!rlang::syms(c(columnOrdV[c], columnV[c]))) %>%
      arrange(!!rlang::sym(columnOrdV[c]))

    names(columnColor[[c]]) <- tempColorDT[[columnV[c]]]

    annAllColors[[c]] <- columnColor[[c]]
    names(annAllColors)[[c]] <- columnV[c]
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


  # circulation vars
  PAGE.Vector <- unique(final$self_page)

  allTmpPic <- NULL
  for (p in 1:length(PAGE.Vector)){

    tempDT2 <- final %>%
      arrange(self_page, self_row, self_columnOrd) %>%
      filter(self_page == PAGE.Vector[p]) %>%
      pivot_wider(id_cols = c(self_page, self_row),
                  names_from = self_columnOrd,
                  values_from = self_res)

    tempDT3 <- tempDT2 %>%
      select(-c(self_page, self_row)) %>%
      as.data.frame()

    # assign row name for heatmap data set
    rownames(tempDT3) <- tempDT2$self_row


    # get annotate color for columns per page
    # as values for columns may vary with page variable
    annoCol <- as_tibble(names(tempDT3)) %>%
      mutate(self_columnOrd = as.numeric(value)) %>%
      left_join(columnOrd.data) %>%
      separate_wider_delim(self_column, delim = '_', names = columnV2) %>%
      select(all_of(rev(columnV)), self_columnOrd) %>%
      as.data.frame()
    # arrange(self_columnOrd)

    rownames(annoCol) <- annoCol$self_columnOrd
    annoCol2 <- annoCol %>%
      select(-self_columnOrd)



    annSubColors <- annAllColors
    for (c in 1:length(columnOrdV)){

      tempSubColor <- unlist(unique(annoCol[columnV[c]]))

      annSubColors[[c]] <- annAllColors[[c]][tempSubColor]
      names(annSubColors)[[c]] <- columnV[c]
    }




    for (c in 1:clusterTypeCnt){

      # get gap id used to separate columns
      if (c < 3 & analysisType == 'GM'){
        gapColumnDT <- annoCol %>%
          ungroup() %>%
          arrange(self_columnOrd) %>%
          mutate(self_gapID = row_number()) %>%
          group_by(!!sym(columnV[1])) %>%
          slice_tail()

        gapColumnID <- unique(gapColumnDT$self_gapID)
      }else{
        gapColumnID <- NULL
      }



      if (is.null(subFolder)){
        tmpPic <- file.path(REPO, paste0(pgName, p, c, '.png'))
      }else{
        tmpPic <- file.path(REPO, paste0(subFolder, '/', pgName, p, c, '.png'))
      }


      # standardize the data before invoking pheatmap (parameter = scale) to
      # avoid issues caused by sd = 0 for certain readouts
      tempDT4 <- t(scale(t(tempDT3))) %>%
        as.data.frame() %>%
        filter(!if_all(everything(), is.nan))


      pheatmap::pheatmap(tempDT4,
                         # scale = 'row',
                         show_colnames = F,
                         cluster_rows = clusterRow[c],
                         cluster_cols = clusterColumn[c],
                         clustering_distance_rows = 'euclidean',
                         clustering_distance_cols = 'euclidean',
                         clustering_method = 'average',
                         annotation_col = annoCol2,
                         annotation_colors = annSubColors,
                         annotation_names_col = F,
                         display_numbers = F,
                         border_color = "black",
                         gaps_col = gapColumnID,
                         fontsize = 8,
                         filename = tmpPic,
                         width = widthP,
                         height = heightP,
                         type='cairo'
      )


      # Define plot object
      tempP <- create_plot(tmpPic,
                           height = heightP,
                           width = widthP,
                           borders = 'none') %>%
        titles(paste0(str_to_title(PAGE.Vector[p]), ' with ', clusterType[c]),
               align = 'center',
               blank_row = 'none')


      rpt <- rpt %>%
        add_content(tempP, page_break = T)


      rm(tempP)
      allTmpPic <- c(allTmpPic, tmpPic)
      assign(x = 'allTmpPic', value = allTmpPic, envir = .GlobalEnv)

    }
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



