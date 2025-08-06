#' bio_box
#'
#' @description This function generates a box plot of biomarker readouts on a logarithmic scale with a base of 10.
#' At the top of the plot, the p-values of pairwise comparisons between visits are displayed,
#' and at the bottom, the numbers of each box are shown respectively.
#' The output is stored in 1 RTF file under REPORT/OUTPUT/subFolder under eWISE environment,
#' or current working path/subFolder under local environment.
#'
#' @note In the input dataset indt, in addition to the required parameters, there must also be the variable USUBJID.
#'
#' @param indt Required, a R dataset name which is used to create the box plot, such as adrc.
#' @param pageV Required, a character string representing the page variable name, such as 'PARAM'.
#' Each value of page variable will have a box plot on one page.
#' @param visCharV Required, a character string representing the character visit variable name, such as 'VISIT'.
#' @param visNumV Required, a character string representing the numeric visit variable name, such as 'VISITNUM'.
#' @param xCharV Required, a character string representing the character x variable name, such as 'TRT01A'.
#' @param xNumV Required, a character string representing the numeric x variable name, such as 'TRT01AN'.
#' @param xColor Required, a character vector representing the color for each value of xCharV,
#' so the length of xColor should be equal to the number of different values of xCharV,
#' Such as c(rep('purple', 3), 'gray') if there are 4 different values for xCharV.
#' @param yV Required, a character string representing the variable name which is used as the y axis variable, such as 'AVAL'.
#' @param groupV A character vector representing the group variable names, such as 'STUDYID' or c('STUDYID', 'AGEGR1'),
#' the default value is NULL.
#' @param tit A character string representing the title of output, the default value is 'Biomarker Box Plot'.
#' @param footn A character string representing the nth footnote displayed as footnote, n = 3, ..., 8, the default value is ''.
#' @param pgName A character string representing the program name displayed in the last footnote,
#' the default value is 'bio_box'
#' @param subFolder A character string representing sub folder name where RTF file locates, such as 'Primary', 'CMI', ...
#' the default value is NULL.
#' @param groupSize A floating number representing the font-size of group variable value (groupV) displayed at the top of box plot,
#' such as 3.5
#' @param nSize A floating number representing the font-size of each number displayed at the bottom of box plot, such as 2.5
#' @param visSize A floating number representing the font-size of visit variable value (visCharV) displayed at the top of box plot,
#' such as 2.5
#' @param pvalueSize A floating number representing the font-size of the p-values displayed at the top of box plot, such as 1.8
#'
#' @examples
#' \dontrun{
#'
#' bio_box(
#' indt = adbm,
#' pageV = 'PARAM',
#' visCharV = 'VIS',
#' visNumV = 'VISITNUM',
#' xCharV = 'TRT01A',
#' xNumV = 'TRT01AN',
#' xColor = c('gray'),
#' yV = 'AVAL',
#' )
#'
#' }
#'
#' @export
bio_box <- function(indt,
                    pageV = 'PARAM',
                    visCharV = 'VIS',
                    visNumV = 'VISITNUM',
                    xCharV = 'TRT01A',
                    xNumV = 'TRT01AN',
                    xColor = c('gray'),
                    yV = 'AVAL',
                    groupV = NULL,
                    tit = 'Biomarker Box Plot',
                    foot3 = '',
                    foot4 = '',
                    foot5 = '',
                    foot6 = '',
                    foot7 = '',
                    foot8 = '',
                    pgName = 'bio_box',
                    subFolder = NULL,
                    groupSize = NULL,
                    nSize = NULL,
                    visSize = NULL,
                    pvalueSize = NULL) {


  if (is.null(groupV)){

    temp <- indt %>%
      ungroup() %>%
      mutate(self_group = '')

  }else{

    temp <- indt %>%
      ungroup() %>%
      mutate(self_group = paste(!!!rlang::syms(groupV), sep = ' '))

  }


  temp <- temp %>%
    ungroup() %>%
    mutate(self_trtn = get(xNumV),
           self_trtc = get(xCharV),
           self_page = get(pageV),
           self_visc = get(visCharV),
           self_visn = get(visNumV),
           self_log10y = log10(get(yV))) %>%
    select(starts_with('self_'), USUBJID)


  # x axis variable order (xn) in plot
  trtDT.ord <- temp %>%
    ungroup() %>%
    distinct(self_page, self_group, self_trtn, self_trtc) %>%
    arrange(self_page, self_group, self_trtn) %>%
    group_by(self_page) %>%
    mutate(xn = row_number())

  trtDT.x.color <- trtDT.ord %>%
    ungroup() %>%
    distinct(self_trtn, self_trtc) %>%
    arrange(self_trtn)


  # in case the number of color is not equal to the number of treatment
  if (length(xColor) < nrow(trtDT.x.color) ){

    xColor <- rep(xColor, ceiling(nrow(trtDT.x.color)/length(xColor)))[1:nrow(trtDT.x.color)]

  }else if (length(xColor) > nrow(trtDT.x.color)) {

    xColor <- xColor[1:nrow(trtDT.x.color)]

  }

  trtDT.x.color$self_color <- xColor


  # add xn and color value into the analysis data set
  final <- left_join(temp, trtDT.ord) %>%
    left_join(trtDT.x.color) %>%
    mutate(fillvar = paste0(self_color, ': ', self_visc))


  # prepare the fill color for plot
  visNum <- final %>%
    ungroup() %>%
    distinct(self_visn, self_visc) %>%
    arrange(self_visn)

  visNumCnt <- nrow(visNum)

  blueP <- colorRampPalette(c('#BFEFFF', '#1C86EE'))(visNumCnt)
  yellowP <- colorRampPalette(c('#FFEC8B', '#8B814C'))(visNumCnt)
  greenP <- colorRampPalette(c('#CAFF70', '#6E8B3D'))(visNumCnt)
  purpleP <- colorRampPalette(c('darkorchid1', 'darkorchid4'))(visNumCnt)
  grayP <- colorRampPalette(c('#E0E0E0', '#545454'))(visNumCnt)

  if (visNumCnt > 12){
    otherP <-  brewer.pal(12, "Paired")
    repColorCnt <- ceiling(visNumCnt/12)
    otherP <- rep(otherP, repColorCnt)[1:visNumCnt]
  }else{
    otherP <-  brewer.pal(visNumCnt, "Paired")
  }



  visDT <- final %>%
    ungroup() %>%
    distinct(self_visn, self_visc, self_color, fillvar) %>%
    arrange(self_visn, self_color) %>%
    group_by(self_visn) %>%
    mutate(idn = cur_group_id()) %>%
    select(-self_visn)

  visDT$colorV <- sapply(1:nrow(visDT), function(i) {
    if (visDT$self_color[i] == 'green') {
      greenP[visDT$idn[i]]
    } else if (visDT$self_color[i] == 'blue'){
      blueP[visDT$idn[i]]
    } else if (visDT$self_color[i] == 'purple'){
      purpleP[visDT$idn[i]]
    } else if (visDT$self_color[i] == 'yellow'){
      yellowP[visDT$idn[i]]
    } else if (visDT$self_color[i] == 'gray'){
      grayP[visDT$idn[i]]
    }else{
      otherP[visDT$idn[i]]
    }
  })


  colorV <- visDT$colorV
  names(colorV) <- visDT$fillvar



  # LMM: p-values of pairwise comparisons between visits
  # LMM: p-values of pairwise comparisons between visits
  # LMM: p-values of pairwise comparisons between visits


  # circulation vars
  PAGE.Vector <- unique(final$self_page)

  pDT <- NULL
  for (p in 1:length(PAGE.Vector)){

    modelD <- final

    # factor in model
    modelD$self_trtn <- factor(modelD$self_trtn, levels = trtDT.x.color$self_trtn)
    modelD$self_visc <- factor(modelD$self_visc, levels = visNum$self_visc)


    modelD2 <- modelD %>%
      filter(self_page == PAGE.Vector[p])

    GROUP.Vector <- unique(modelD2$self_group)
    for (g in GROUP.Vector){

      modelD3 <- modelD2 %>%
        filter(self_group == g)

      sub.value <- modelD3 %>%
        ungroup() %>%
        distinct(USUBJID, self_log10y) %>%
        group_by(USUBJID) %>%
        summarise(cnt = n()) %>%
        ungroup() %>%
        summarise(sub.val.cnt = max(cnt))

      y.value.cnt <- length(unique(modelD3$self_log10y))
      if (y.value.cnt == 1 | sub.value$sub.val.cnt == 1){
        next
      }


      trtVar <- unique(modelD3$self_trtn)
      if (length(trtVar) == 1){

        lmm.fit <- lmer(self_log10y ~ self_visc + (1|USUBJID), modelD3)
        lmm.res <- emmeans(lmm.fit, ~ self_visc)

        modelR <- as.data.frame(pairs(lmm.res, adjust = NULL)) %>%
          mutate(contrast = str_replace_all(str_replace_all(contrast,' - ','-'),' ','-')) %>%
          separate(col = contrast, sep = "-", into = c('vis1', 'vis2'), remove = F) %>%
          mutate(trt1 = paste0('self_trtn', trtVar ), trt2 = trt1)


      }else{

        lmm.fit <- lmer(self_log10y ~ self_trtn*self_visc + (1|USUBJID),  modelD3)
        lmm.res <- emmeans(lmm.fit, ~ self_trtn*self_visc)

        modelR <- as.data.frame(pairs(lmm.res, adjust = NULL)) %>%
          mutate(contrast = str_replace_all(str_replace_all(contrast,' - ','-'),' ','-')) %>%
          separate(col = contrast, sep = "-", into = c('trt1', 'vis1', 'trt2', 'vis2'), remove = F)

      }

      modelR2 <- modelR %>%
        filter(vis1 != vis2 & trt1 == trt2) %>%
        mutate(self_page = PAGE.Vector[p],
               self_group = g,
               self_trtn = as.numeric(str_sub(trt1, 10) ),
               pValue = case_when(is.na(p.value) ~ '',
                                  p.value < 0.001 ~ '<0.001',
                                  !is.na(p.value) ~ formatC(p.value, format = 'f', digits = 3)))


      pDT <- bind_rows(pDT, modelR2)

    }
  }


  pDT <- pDT %>%
    ungroup()

  # calculate the x position for p value on each page (usually is variable PARAM)
  visPeriodCnt <- pDT %>%
    group_by(self_page) %>%
    distinct(vis1) %>%
    summarise(visPeriodCnt = n())

  visOrder1 <- pDT %>%
    distinct(self_page, vis1) %>%
    mutate(self_visc = vis1)

  visOrder2 <- pDT %>%
    distinct(self_page, vis2) %>%
    mutate(self_visc = vis2)

  visOrder3 <- bind_rows(visOrder1, visOrder2) %>%
    distinct(self_page, self_visc) %>%
    left_join(visNum) %>%
    arrange(self_page, self_visn) %>%
    group_by(self_page) %>%
    mutate(visIDn = row_number()) %>%
    select(-self_visn)


  pDT2 <- left_join(pDT, visPeriodCnt) %>%
    left_join(trtDT.ord) %>%
    left_join(visOrder3, by = join_by(self_page, vis1 == self_visc)) %>%
    rename(visIDn1 = visIDn) %>%
    left_join(visOrder3, by = join_by(self_page, vis2 == self_visc)) %>%
    rename(visIDn2 = visIDn) %>%
    mutate(xPeriod = 0.25,
           x1_position = xn - visPeriodCnt*xPeriod/2 + (visIDn1-1)*xPeriod,
           x2_position = xn - visPeriodCnt*xPeriod/2 + (visIDn2-1)*xPeriod) %>%
    ungroup() %>%
    group_by(self_page, xn) %>%
    arrange(self_page, xn, visIDn1, visIDn2) %>%
    mutate(yIDn = row_number())




  # outout to RTF
  footVars <- unlist(mget(paste0("foot", 3:8)))

  outT <- 'rtf'
  dtsource <- 'adbm'
  heightP <- 4.9-0.2*sum(footVars != '')
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


  # save data
  saveRDS(pDT, file = file.path(REPD, paste0(pgName, '.Rdata')))


  # footnotes
  foot1 <- 'Geomertic Mean is showed in red diamond.'
  foot2 <- paste0('The p value for comparison between visits is from model: log(AVAL) ~ ',
                  paste0(xCharV, '*' ,visCharV), ' + (1|USUBJID), ',
                  'and *** is for p < 0.001, ** for p < 0.01, * for p < 0.05 and ns for not significant.')
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

    tempDT2 <- final %>%
      filter(self_page == PAGE.Vector[p])


    # create dummy structure, so each treatment (x variable) has same visit position per page
    dumVisDT <- tempDT2 %>%
      distinct(self_visn, self_visc)

    dumTrtDT <- tempDT2 %>%
      distinct(self_page, self_group, self_trtc, self_color, xn) %>%
      arrange(xn)

    dumPageDT <- cross_join(dumVisDT, dumTrtDT) %>%
      mutate(fillvar = paste0(self_color, ': ', self_visc))



    pDT4 <- pDT2 %>%
      filter(self_page == PAGE.Vector[p] & pValue != '')


    # factor in input data set
    tempDT2$fillvar <- factor(tempDT2$fillvar, levels = visDT$fillvar)
    tempDT2$xn <- factor(tempDT2$xn, levels = dumTrtDT$xn)

    dumPageDT$fillvar <- factor(dumPageDT$fillvar, levels = visDT$fillvar)
    dumPageDT$xn <- factor(dumPageDT$xn, levels = dumTrtDT$xn)

    pDT4$xn <- factor(pDT4$xn, levels = dumTrtDT$xn)


    # label for x axis
    nameVar <- dumTrtDT$self_trtc
    names(nameVar) <- dumTrtDT$xn


    # Y axis
    # add more space for avoiding overlay of p value & count number with boxplot
    textYdt <- tempDT2 %>%
      ungroup() %>%
      group_by(self_page) %>%
      summarise(minY = floor(min(self_log10y*10))/10,
                maxY = ceiling(max(self_log10y*10))/10,
                minY2 = minY - 0.02*(maxY - minY),
                maxY2 = maxY + 0.22*(maxY - minY)
      )


    # keep integer y axis
    ytickV <- seq(floor(textYdt$minY2), ceiling(textYdt$maxY2), by = 1 )
    # ytickExpV <- formatC(10**ytickV, format = "f")

    digits <- c(abs(ytickV[ytickV < 0]), rep(0, length(ytickV[ytickV >= 0])))
    ytickExpV <- sprintf(paste0("%.", digits, "f"), 10**ytickV)



    minY2 <- textYdt$minY2
    maxY2 <- textYdt$maxY2


    # calculate n for each treatment and visit
    nDT <- tempDT2 %>%
      ungroup() %>%
      group_by(self_page, self_group, self_visn, self_visc, fillvar, xn, self_trtc, self_color) %>%
      summarise(n = n()) %>%
      right_join(dumPageDT) %>%
      mutate(n = ifelse(is.na(n), 0, n),
             nc = paste0('n=', n )) %>%
      left_join(textYdt)


    pvalueDT <- full_join(pDT4, textYdt) %>%
      mutate(y_position = maxY2 - (0.01 + (!is.null(groupV))*0.07)*(maxY2-minY2) - yIDn*0.04*(maxY2-minY2),
             ytxt_position = y_position + 0.02*(maxY2-minY2),
             fillvar = '',
             pMarker = case_when(p.value >= 0.05 ~ 'ns',
                                 p.value < 0.001 ~ '***',
                                 p.value < 0.01 ~ '**',
                                 p.value < 0.05 ~ '*',
             ))


    dayDT <- nDT %>%
      mutate(dayc = self_visc,
             dayYid = maxY2 - (!is.null(groupV))*0.07*(maxY2-minY2))


    colorSepLine <- nDT %>%
      ungroup() %>%
      distinct(self_page, self_group, self_color, xn, self_trtc) %>%
      arrange(self_page, self_group, xn, self_trtc, self_color) %>%
      mutate(lineId = with(rle(self_color), rep(seq_along(values), lengths))) %>%
      group_by(self_page, self_group, lineId) %>%
      slice_tail() %>%
      mutate(lineXn = as.numeric(xn) + 0.5) %>%
      arrange(lineXn)

    colorSepLine.Vactor <- sort(colorSepLine$lineXn, decreasing = T)[-1]

    pd <- position_dodge(width = 0.73)

    nLength <- nrow(nDT)
    if (is.null(groupSize)){groupSize <- 4 - 0.004*nLength}
    if (is.null(visSize)){visSize <- 4 - 0.05*nLength}
    if (is.null(pvalueSize)){pvalueSize <- 4 - 0.03*nLength}
    if (is.null(nSize)){nSize <- 4 - 0.06*nLength}


    tempDT3 <- full_join(tempDT2, dumPageDT) %>%
      left_join(textYdt) %>%
      mutate(self_log10y2 =  ifelse(!is.na(self_log10y), self_log10y, maxY * 1000))


    tempP <- ggplot(tempDT3, aes(xn, self_log10y2, fill = fillvar)) +
      geom_boxplot(outlier.size = 0.5, size = 0.3) +
      geom_vline(xintercept = colorSepLine.Vactor, linetype = "dashed", linewidth = 0.3) +
      stat_summary(fun = mean, geom = 'point', show.legend = F, color = "red",
                   size = 2, shape = 18, position = pd) +
      geom_text(data = dayDT, size = visSize,  aes(xn, dayYid, label = dayc), position = pd)

    if (!is.null(groupV)){

      groupDT <- nDT %>%
        ungroup() %>%
        distinct(self_group, maxY2, xn) %>%
        group_by(self_group, maxY2) %>%
        summarise(cnt = n()) %>%
        ungroup() %>%
        mutate(sumCnt = cumsum(cnt)) %>%
        mutate(studyXid = sumCnt - cnt/2 + 0.5,
               fillvar = '')

      tempP <- tempP +
        geom_text(data = groupDT, size = groupSize,  aes(studyXid, maxY2, label = self_group))

    }


    tempP <- tempP +
      geom_text(data = nDT, size = nSize,  aes(xn, minY2, label = nc), position = pd) +
      geom_text(data = pvalueDT, aes(x = (x1_position + x2_position)/2, y = ytxt_position, label = pMarker),
                color = '#6E6E6E', size = pvalueSize, hjust = 0.5) +
      geom_segment(data = pvalueDT, aes(x = x1_position, xend = x2_position, y = y_position, yend = y_position),
                   color = '#6E6E6E', linewidth = 0.3, linetype = 1) +
      scale_fill_manual(values = colorV) +
      scale_x_discrete(labels = nameVar) +
      scale_y_continuous(breaks = ytickV, label = ytickExpV) +
      labs(title = paste0(PAGE.Vector[p]),
           x = '',
           y = '% of positive cells (log)',
           fill = '') +
      coord_cartesian(ylim = c(minY2, maxY2)) +
      theme_classic() +
      theme(text = element_text(size = 9),
            legend.position = 'none'
      )


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

    footVars2 <- footVars[footVars != '']
    tempP5 <- create_plot(tmpPic,
                          height = heightP,
                          width = widthP,
                          borders = 'none') %>%
      title_header(' ', tits, right = c('Page [pg] of [tpg]',' '), borders = 'bottom', blank_row = 'below') %>%
      footnotes(paste0(c(foot1, foot2, footVars2, foot9), collapse = '\n'),
                borders = 'top',
                blank_row = 'none')

    rpt <- rpt %>%
      add_content(tempP5, page_break = T)

    allTmpPic <- c(allTmpPic, tmpPic)
    assign(x = 'allTmpPic', value = allTmpPic, envir = .GlobalEnv)

  }



  # Write out report
  res <- write_report(rpt)


  # remove all PNG
  print(length(allTmpPic))
  print(allTmpPic)
  file.remove(allTmpPic)

}



