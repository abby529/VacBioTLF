#' example_data
#'
#' @description This function generates an example data for function to use
#'

# library(dplyr)
# library(tidyr)

subCnt <- 120
demog.data <- tibble(SUBJID = formatC(sample(1:200, subCnt, replace = F), width = 3, format = 'd', flag = '0')) %>%
  mutate( STUDYID = 'TEST',
          SITEID = '999',
          USUBJID = paste(STUDYID, SITEID, SUBJID, sep = '-'),
          AGE =  sample(18:80, subCnt, replace = T),
          AGEGR1N = ifelse(AGE < 60, 1, 2),
          AGEGR1 = ifelse(AGEGR1N == 1, '18-59 years', '>= 60 years'),
          COUNTRY = 'ZZZ',
          SEX = sample(c('F', 'M'), subCnt, replace = T),
          TRT01PN = rep(1:6, rep(20, 6)),
          TRT01P = paste('LNP', TRT01PN),
          TRT01AN = TRT01PN,
          TRT01A = TRT01P) %>%
  relocate(STUDYID, .before = SUBJID) %>%
  relocate(USUBJID, .after = STUDYID)


param.data <- tibble(PARAM = c("Memory B cell (%)",
                               "Memory B cell IgA+ (%)",
                               "Memory B cell IgA+specific (%)",
                               "Memory B cell IgG+ (%)",
                               "Memory B cell IgG+specific (%)",
                               "Memory B cell IgM+ (%)",
                               "Memory B cell IgM+specific (%)"),
                     PARAMCD = c("MBCELL", "MBIGA",  "MBIGAT", "MBIGG", "MBIGGT", "MBIGM", "MBIGMT")) %>%
  mutate(PARCAT1 = 'Antibody Mediated Immunity',
         PARCAT2 = 'Memory B Cell')


# reference value
set.seed(2025)
example.data.v1 <- crossing(demog.data, param.data) %>%
  rowwise() %>%
  mutate(VISIT = "VISIT 01/D1",
         VISITNUM = 10,
         AVAL = case_when(PARAMCD == "MBCELL" ~ runif(1, 1.29, 37.4),
                          PARAMCD == "MBIGA" ~ runif(1, 3.6, 68),
                          PARAMCD == "MBIGAT" ~ runif(1, 0.004, 1.14),
                          PARAMCD == "MBIGG" ~ runif(1, 15.2, 80),
                          PARAMCD == "MBIGGT" ~ runif(1, 0.007, 0.9),
                          PARAMCD == "MBIGM" ~ runif(1, 0.57, 28.9),
                          PARAMCD == "MBIGMT" ~ runif(1, 0.012, 4.43)),
         ADT = lubridate::dmy("18Apr2023") + sample(-1:3, 1))

set.seed(2025)
example.data.v2 <- example.data.v1 %>%
  rowwise() %>%
  mutate(VISIT = "VISIT 04/D9",
         VISITNUM = 40,
         AVAL = case_when(PARAMCD %in% c("MBIGAT", "MBIGGT", "MBIGMT") ~ runif(1, 1, 15)*AVAL,
                          PARAMCD != '' ~ runif(1, 1, 1.25)*AVAL),
         ADT = ADT + sample(8:10, 1))


set.seed(2025)
example.data.v3 <- example.data.v1 %>%
  mutate(VISIT = "VISIT 06/D91",
         VISITNUM = 60,
         AVAL = runif(1, 0.9, 1.2)*AVAL,
         ADT = ADT + sample(91:92, 1))



example.data.base <- example.data.v1 %>%
  rename(BASE = AVAL) %>%
  select(-starts_with('VISIT'), -ADT)

example.data <- bind_rows(example.data.v1, example.data.v2, example.data.v3) %>%
  arrange(USUBJID, PARAMCD, VISITNUM) %>%
  left_join(example.data.base) %>%
  mutate(CHG = case_when(VISITNUM != 10 ~ AVAL - BASE,
                         VISITNUM == 10 ~ NA),
         R2BASE = case_when(VISITNUM != 10 ~ AVAL/BASE,
                            VISITNUM == 10 ~ NA),
         FOLDRISE = R2BASE,
         ABLFL = case_when(VISITNUM != 10 ~ "",
                           VISITNUM == 10 ~ "Y"),
  ) %>%
  relocate(c(VISIT, VISITNUM), .after = ABLFL) %>%
  relocate(ADT, .before = ABLFL)


# save data
saveRDS(example.data, file = here::here('inst', 'extdata', 'example.data.rds'))


