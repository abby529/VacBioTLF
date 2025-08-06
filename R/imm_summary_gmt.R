#' imm_summary_gmt
#'
#' @description This function generates gmt table for immunogenicity in 1 RTF file
#' stored under REPORT/OUTPUT/subFolder under eWISE environment,
#' or current working path/subFolder under local environment.
#' @note There must be variable ATOXGD01, ..., ATOXGD09 in the input dataset indt
#'
#' @param indt A R dataset name which contains the immunogenicity data, such as adis
#' @param indt2 A R dataset name which used to calculate the big N in the header, such as adsl
#' @param groupV A character vector presents category variable names used in the table and can be found in indt,
#' such as 'PARAM' or c('AGE', 'PARAM)
#' @param groupLableV A character vector with same length of paramter groupV presents the lable of category variable,
#' default value is ''. when groupV = 'PARAM', groupLableV can be 'Parameter'
#' @param grpColSize A float vector presents the width for category variables, such as 6.5 or c(5, 6.5),
#' default value is 6.5 for each category variable
#' @param timeV A character string presents visit variable name used in the table and can be found in indt,
#' such as VISIT
#' @param timeColSize A float number presents the width for visit variable, default value is 2.5
#' @param trtCntPerPage A integer number presents number of treatment on each page, default value is 3
#' @param statColSize A float number vector of length 3 for the width of M, GMT/GMTR and (95% CI) columns,
#' the default value if c(1, 2.1, 2.3)
#' @param avalV A character string presents analysis variable name used to calculate GMT and can be found in indt,
#' the default value is 'AVAL'
#' @param foldV A character string presents fold rise variable name used to calculate GMTR and can be found in indt,
#' the default value is 'FOLDRISE'
#' @param trtnV A character string presents the numeric treatment name in indt, such as 'TRT01AN'
#' @param trtcV A character string presents the character treatment name in indt and indt2, such as 'TRT01A'.
#' And there is a 121 relationship between trtnV & trtcV
#' @param pgName A character string presents program name showing in the last footnote in the RTF file,
#' the default value is 'imm_summary_gmt'
#' @param titV A character string presents the title showing in the RTF file, the default value is ''
#' @param footxV x = 1, ...,8, a character string presents the footnote showing in the RTF file, the default value is ''
#'
#' @examples
#' \dontrun{
#'
#' imm_summary_gmt(indt = final,
#' indt2 = adsl,
#' groupV = 'PARAM',
#' groupLableV = 'Parameter',
#' grpColSize = 6.5,
#' timeV = 'VIS',
#' timeColSize = 2.5,
#' trtCntPerPage = 3,
#' statColSize = c(1, 2.1, 2.3),
#' trtnV = 'TRT01AN',
#' trtcV = 'TRT01A',
#' avalV = 'AVAL',
#' foldV = 'FOLDRISE',
#' pgName = 'imm_summary_gmt',
#' titV = 'Table 2.55 Summary of geometric means of antibody titers after vaccination - mFAS',
#' foot1V = 'M: number of participants with available data for the relevant endpoint.',
#' foot2V = 'GMT: geometric mean titer. GMTR: geometric mean titer ratio. CI = Confidence interval. NC = Not Computed.',
#' foot3V = 'The 2 sided 95% CIs is calculated based on the t-test.',
#' foot4V = 'For main cohorts, each of the unadjuvanted groups (G10, G11 or G12) are made up of the sum of participants across the main cohorts analyzed.',
#' )
#'
#'}
#'
#' @export
imm_summary_gmt <- function(indt,
                            indt2,
                            groupV = 'PARAM',
                            groupLableV = '',
                            grpColSize = 6.5,
                            timeV = 'VIS',
                            timeColSize = 2.5,
                            trtCntPerPage = 3,
                            statColSize = c(1, 2.1, 2.3),
                            avalV = 'AVAL',
                            foldV = 'FOLDRISE',
                            trtnV = 'TRT01AN',
                            trtcV = 'TRT01A',
                            pgName = 'imm_summary_gmt',
                            titV = '',
                            foot1V = '',
                            foot2V = '',
                            foot3V = '',
                            foot4V = '',
                            foot5V = '',
                            foot6V = '',
                            foot7V = '',
                            foot8V = ''
){

  tempDT <- indt %>%
    ungroup() %>%
    mutate(VISDY = case_when(str_detect(get(timeV), '(?<=D)\\d+') ~ as.numeric(str_extract(get(timeV), '(?<=D)\\d+')),
                             str_detect(get(timeV), '(?<=W)\\d+') ~ as.numeric(str_extract(get(timeV), '(?<=W)\\d+'))*7,
                             str_detect(get(timeV), '(?<=M)\\d+') ~ as.numeric(str_extract(get(timeV), '(?<=M)\\d+'))*30,
                             str_detect(get(timeV), '(?<=Y)\\d+') ~ as.numeric(str_extract(get(timeV), '(?<=Y)\\d+'))*365
    ))

  tempDT2 <- tempDT %>%
    mutate(log10V = log10(get(avalV)),
           timePoint = get(timeV),
           timePartn = 1)

  d1C <- tempDT2 %>%
    filter(VISDY == 1) %>%
    distinct(timePoint)

  tempDT3 <- tempDT2 %>%
    filter(VISDY != 1) %>%
    mutate(log10V = log10(get(foldV)),
           timePoint = paste0('Ratio ', get(timeV), '/', d1C$timePoint),
           timePartn = 2)

  tempDT4 <- bind_rows(tempDT2, tempDT3) %>%
    filter(!is.na(log10V)) %>%
    ungroup() %>%
    group_by(!!!rlang::syms(groupV), timePartn, VISDY, timePoint, !!rlang::sym(trtnV)) %>%
    summarise(m = n(), logMean = mean(log10V), sd = sd(log10V), df = m-1)



  # calculate 95% CI based on t distribution
  # calculate 95% CI based on t distribution
  alpha <- 0.05
  quanT <- qt(1-alpha/2, tempDT4$df)
  tempDT4$quanT <- quanT

  tempDT5 <- tempDT4 %>%
    mutate(logLower = logMean - quanT*sd/sqrt(m), logUpper = logMean + quanT*sd/sqrt(m))

  dumRow <- tempDT5 %>%
    ungroup() %>%
    distinct(!!!rlang::syms(groupV), timePartn, VISDY, timePoint)

  dumDT <- tempDT5 %>%
    ungroup() %>%
    distinct(!!rlang::sym(trtnV)) %>%
    cross_join(dumRow)


  tempDT6 <- left_join(dumDT, tempDT5) %>%
    mutate(mc = ifelse(is.na(m), '0', sprintf('%.0f', m)),
           gmc = ifelse(is.na(m), 'NC', as.character(signif(10**logMean, digits = 3))),
           logLowerc = ifelse(is.na(m) | m <= 5, 'NC', as.character(signif(10**logLower, digits = 3))),
           logUpperc = ifelse(is.na(m) | m <= 5, 'NC', as.character(signif(10**logUpper, digits = 3))),
           cic = paste0('(', logLowerc, ' ; ', logUpperc, ')')
    ) %>%
    pivot_wider(id_cols = c(all_of(groupV), timePartn, VISDY, timePoint),
                names_from = all_of(trtnV),
                values_from = c(mc, gmc, cic))



  # prepare for report
  statColName <- as.data.frame(names(tempDT6) ) %>%
    rename_with(~'colName', everything()) %>%
    filter(str_detect(colName, 'mc|gmc|cic')) %>%
    mutate(trtn = as.numeric(str_extract(colName, '(?<=_)\\d+')),
           colName2 = str_split_i(colName, '_', 1),
           coln = case_when(colName2 == 'mc' ~ 1,
                            colName2 == 'gmc' ~ 2,
                            colName2 == 'cic' ~ 3) ) %>%
    arrange(trtn, coln)

  statColNameV <- statColName$colName

  tb10 <- tempDT6 %>%
    select(c(all_of(groupV), timePartn, VISDY, timePoint, all_of(statColNameV))) %>%
    arrange(!!!rlang::syms(groupV), timePartn, VISDY, timePoint) %>%
    select(-c(timePartn, VISDY))




  # treatment column header & bigN
  trtDT <- indt2 %>%
    group_by(!!rlang::sym(trtnV), !!rlang::sym(trtcV)) %>%
    summarise(bigN = n()) %>%
    arrange(!!rlang::sym(trtnV)) %>%
    mutate(treatmentName = paste0(get(trtcV), '\n (N=', bigN,')')  )




  # create folder to store the RTF file
  if (!exists("REPO")){
    REPO <- here::here()
    W_STUDY <- ''
  }

  folder_path <- file.path(REPO, 'Primary')
  if (!file.exists(folder_path)) {
    dir.create(folder_path, recursive = TRUE)
  }


  # create output using reporter package
  pgname <- pgName
  outT <- 'rtf'
  fname <- paste0('Primary/', pgname, '.',outT)
  dtsource <- 'adis'

  # Create temporary path
  tmp <- file.path(REPO, fname)

  # Create the report
  rpt <-create_report(tmp, output_type = toupper(outT), font = 'Times', font_size = 9, paper_size='letter', units = 'cm') %>%
    # A4 size
    set_margins(top = 4.3, bottom = 2.5, left = 2.8, right = 1)





  tbl_tt <- create_text('')  %>%
    titles(titV,
           align = 'left',
           blank_row = 'none')


  # The last footnote
  foot9 <- paste0('Study: ', W_STUDY,' Program: ', tolower(pgname), '.R Datasets=',dtsource,
                  ' Output:/', str_replace_all(REPO,'~/wise/',''), '/', fname,
                  ' (', format(Sys.time(),'%d%b%Y %H:%M'),')')


  footVars <- unlist(mget(paste0("foot", 1:8, "V")))
  footVars2 <- footVars[footVars != '']

  tbl_ftnt <- create_text('')  %>%
    footnotes(paste0(c(footVars2, foot9), collapse = '\n'),
              borders = 'top',
              blank_row = 'none')

  pValueV <- statColNameV[-length(statColNameV)]

  # Define table
  fin <- create_table(tb10, show_cols = 'all', header_bold = T, continuous = F, border = 'top', first_row_blank = F)

  treatmentnV <- sort(unique(statColName$trtn))
  for (v in 1:length(treatmentnV)){
    fin <- fin %>%
      spanning_header(from = paste0('mc_', treatmentnV[v] ), to = paste0('cic_', treatmentnV[v] ),
                      label = trtDT$treatmentName[v],
                      bold = T, underline = F, standard_eval = T)

    if (v %% trtCntPerPage == 1 & v > trtCntPerPage){
      fin <- fin %>%
        define(paste0('mc_', treatmentnV[v]), label = 'M', width = statColSize[1], align = 'center', standard_eval = T, page_wrap = T)

    }else{
      fin <- fin %>%
        define(paste0('mc_', treatmentnV[v]), label = 'M', width = statColSize[1], align = 'center', standard_eval = T)
    }

  }

  paramCnt <- length(groupV)

  if (length(grpColSize) < paramCnt){
    grpColSize <- c(grpColSize, rep(3, paramCnt - length(grpColSize)))
  }

  for (p in 1:paramCnt){
    fin <- fin %>%
      define(groupV[p], label = groupLableV[p], width = grpColSize[p], id_var = T, align = 'left', blank_after = T, dedupe = T, standard_eval = T)
  }



  fin <- fin %>%
    define(timePoint, label = 'Time Point', width = timeColSize, id_var = T, align = 'left', blank_after = F) %>%
    define(statColNameV[str_detect(statColNameV, 'gmc_')], label = 'GMT/GMTR', width = statColSize[2], align = 'center', standard_eval = T) %>%
    define(pValueV[str_detect(pValueV, 'cic_')],  label = '(95% CI)', width = statColSize[3], align = 'center', standard_eval = T) %>%
    define(statColNameV[length(statColNameV)], label = '(95% CI)', width = statColSize[3], align = 'center', standard_eval = T)


  # Create the report
  rpt <- rpt %>%
    add_content(tbl_tt, page_break = F, blank_row = 'none') %>%
    add_content(fin, page_break = F, align = 'left' ) %>%
    add_content(tbl_ftnt, page_break = T, blank_row = 'none')


  # Write the report
  write_report(rpt)


}

