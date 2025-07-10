#' ReactoBar
#'
#' @description This function generates 2 reactogenicity barplots for Max grade over D1-3 and
#' Max grade over D1-9 in 1 RTF file stored under REPORT/OUTPUT/Primary.
#' @note There must be variable ATOXGD01, ..., ATOXGD09 in the input dataset indt
#'
#' @param indt: A R dataset name which is used to create the bar plot, such as adrc
#' @param trtn: A character string presents the numeric treatment name in indt, such as 'TRT01AN'
#' @param trtc: A character string presents the character treatment name in indt, such as 'TRT01A'.
#' And there is a 121 relationship between trtn & trtc
#' @param tits: titles for the 2 plots, a character vector of length 2.
#' default value is c('Figure 1: Barplot for reactogenicity (Max grade over D1-3)',
#' 'Figure 2: Barplot for reactogenicity (Max grade over D1-9)')
#'
#' @examples
#' \dontrun{
#'
#' ReactoBar(indt = adrc,
#' trtn ='TRT01AN',
#' trtc = 'TRT01AC',
#' tits = c('Figure 1: Barplot for reactogenicity (Max grade over D1-3)',
#' 'Figure 2: Barplot for reactogenicity (Max grade over D1-9)'))
#'
#' }
#'
#' @export
ReactoBar <- function(indt, trtn, trtc, tits = c('Figure 1: Barplot for reactogenicity (Max grade over D1-3)',
                                                 'Figure 2: Barplot for reactogenicity (Max grade over D1-9)')) {

    # prepare data
    comb <- indt %>%
        mutate(self_trtn = get(trtn), self_trtc = get(trtc)) %>%
        # cover issue: CEDECOD is missing
        mutate(CEDECOD = ifelse(CEDECOD == '', str_to_title(CETERM), CEDECOD)) %>%
        select(USUBJID, self_trtn, self_trtc, CEDECOD, CECAT, CESCAT, starts_with('ATOXGD')) %>%
        mutate(MAX3 = pmax(ATOXGD01, ATOXGD02, ATOXGD03, na.rm = T),
               MAX9 = pmax(ATOXGD01, ATOXGD02, ATOXGD03, ATOXGD04, ATOXGD05, ATOXGD06, ATOXGD07, ATOXGD08, ATOXGD09, na.rm = T),
               ordn = 1)



    # SOI, SYS
    comb2 <- comb %>%
        group_by(USUBJID, self_trtn, self_trtc, CECAT, CESCAT) %>%
        summarise(MAX3 = max(MAX3), MAX9 = max(MAX9) ) %>%
        mutate(CEDECOD = str_to_sentence(CESCAT), ordn = 2)


    #TOT
    comb3 <- comb %>%
        group_by(USUBJID, self_trtn, self_trtc, CECAT) %>%
        summarise(MAX3 = max(MAX3), MAX9 = max(MAX9) ) %>%
        mutate(CEDECOD = 'Total', ordn = 3)


    final <- bind_rows(comb, comb2, comb3)

    trtDT <- final %>%
        ungroup() %>%
        distinct(self_trtn, self_trtc) %>%
        arrange(self_trtn, self_trtc)


    final$self_trtn <- factor(final$self_trtn,
                            levels = trtDT$self_trtn,
                            labels = trtDT$self_trtc)


    chk <- final %>%
        ungroup() %>%
        distinct(self_trtn, self_trtc)

    # circulation vars
    testD <- final %>%
        distinct(ordn, CECAT, CESCAT, CEDECOD) %>%
        arrange(ordn, CECAT, CESCAT, CEDECOD)

    testV <- unique(testD$CEDECOD)


    # create folder to store the RTF file
    folder_path <- file.path(REPO, 'Primary')
    if (!file.exists(folder_path)) {
        dir.create(folder_path, recursive = TRUE)
    }


    # outout to RTF
    pgname <- 'pri_bar_reacto'
    outT <- 'rtf'
    fname = paste0('Primary/', pgname, '.',outT)
    dtsource <- 'adbm'
    heightP <- 5.5
    widthP <- 9

    titV <- tits

    invarV <- rep(c('MAX3', 'MAX9'), each = 1)
    byvarL <- rep(list(''), 2)



    # The last footnote
    foot9 <- paste0('Study: ', W_STUDY,' Program: ', tolower(pgname), '.R Datasets=',dtsource,
                    ' Output:/', str_replace_all(REPO,'~/wise/',''), '/', fname,
                    ' (', format(Sys.time(),'%d%b%Y %H:%M'),')')


    # Create temporary path
    tmp <- file.path(REPO, fname)


    # Add plot to report
    rpt <- create_report(tmp, output_type = toupper(outT), font = 'Times', font_size = 9) %>%
        set_margins(top = 1, bottom = 0.8) %>%
        options_fixed(font_size = 9)



    allTmpPic <- NULL
    for (k in 1:length(titV)){

        tit <- titV[k]
        invar <- invarV[k]
        byvar <- byvarL[[k]]


        tmpPicList <- list()
        for (i in 1: length(testV) ){

            tempDT <- final %>%
                filter(CEDECOD == testV[i] & !is.na(get(invar)))


            titLab <- paste0(invar, '.', unique(tempDT$CEDECOD))

            tempP <- ggplot(tempDT, aes(self_trtn, fill = as.factor(get(invar))) ) +
                geom_bar(position = 'fill', width = 0.5) +
                scale_fill_manual(values = c('0' = 'paleturquoise3', '1' = 'skyblue4', '2' = '#F0E68C', '3' = '#CD5555')) +
                labs(title = titLab,
                     x = 'Treatment group',
                     y = 'Percent',
                     fill = 'Grade') +
                scale_x_discrete(labels = function(x) str_wrap(x, width = 16)) +
                theme(text = element_text(size = 9),
                      axis.text.x = element_text(angle = 90, vjust = 0.5),
                      axis.ticks.x = element_blank())



            if (byvar[1] == ''){

                tempP2 <- tempP

            }else{

                tempP2 <- tempP +
                    facet_wrap(vars(!!!rlang::syms(byvar)   ))

            }

            tmpPicList[[i]] <- tempP2
        }


        fig_nrow <- floor(sqrt(length(tmpPicList)))
        fig_ncol <- ceiling(length(tmpPicList)/nrow)

        tempP3 <- ggarrange(plotlist = tmpPicList,
                            nrow = fig_nrow,
                            ncol = fig_ncol,
                            hjust = 0)


        tmpPic <- file.path(REPO, paste0('Primary/', pgname, paste0(invar,paste0(byvar, collapse = '')),'.png'))
        ggsave(tmpPic,
               plot = tempP3,
               type='cairo',
               width = 24,
               height = 16,
               units = 'in')


        # Define plot object
        tempP5 <- create_plot(tmpPic,
                              height = heightP,
                              width = widthP,
                              borders = 'none') %>%
            title_header(' ', tit, right = c('Page [pg] of [tpg]',' '), borders = 'bottom', blank_row = 'below') %>%
            footnotes(foot9, borders = 'top', blank_row = 'none')


        rpt <- rpt %>%
            add_content(tempP5, page_break = T)


        rm(tempP, tempP5)
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


