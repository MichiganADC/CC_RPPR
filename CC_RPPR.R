# CC_RPPR.R

download_json <- FALSE
source("~/Box/Documents/R_helpers/config.R")
source("~/Box/Documents/R_helpers/helpers.R")

data_dump_path <- "~/Box/MADC Box Account/Data Core/DataDump/2020-03-15/"

# USEFUL LIBRARIES ----
suppressMessages( library(dplyr) )
suppressMessages( library(stringr) )
suppressMessages( library(tidyr) )
suppressMessages( library(readr) )


# _ Define Fields ----

# _ _ UDS 3 ----

# Header data
fields_u3_hd_raw <-
  c(
    "ptid"           # partic. ID
    , "redcap_event_name"  # REN
    ,"form_date"     # visit date
  )
# Form A1 Demographics
fields_u3_a1_raw <-
  c(
    "sex"            # sex
  ) %>% c(., paste0("fu_", .), paste0("tele_", .))
# Form D1 (IVP, FVP, TVP)
fields_u3_d1_raw <-
  c(
    # Diagnosis
    "normcog"    # Normal
    , "demented" # Demented
    , "amndem"   # Amnestic multidomain dementia syndrome
    , "pca"      # Posterior cortical atrophy syndrome
    , "ppasyn"   # Primary progressive aphasia (PPA) syndrome
    , "ftdsyn"   # Behavioral variant FTD (bvFTD) syndrome
    , "lbdsyn"   # Lewy body dementia syndrome
    , "namndem"  # Non-amnestic multidomain dementia syndrome
    , "mciamem"  # Amnestic MCI, single domain (aMCI SD) 
    , "mciaplus" # Amnestic MCI, multiple domains (aMCI MD)
    , "mcinon1"  # Non-amnestic MCI, single domain (naMCI SD)
    , "mcinon2"  # Non-amnestic MCI, multiple domains (naMCI MD)
    , "impnomci" # Cognitively impaired, not MCI
    # Etiology
    , "alzdis"   # Alzheimer's disease
    , "alzdisif" 
    , "lbdis"    # Lewy body disease
    , "lbdif" 
    , "msa"      # Multiple system atrophy
    , "msaif"
    , "psp"      # Progressive supranuclear palsy (PSP)
    , "pspif"
    , "cort"     # Corticobasal degeneration (CBD)
    , "cortif"
    , "ftldmo"   # FTLD with motor neuron disease
    , "ftldmoif"
    , "ftldnos"  # FTLD NOS
    , "ftldnoif"
    , "cvd"      # Vascular Brain injury (clinical or imaging evidence)
    , "cvdif" 
    , "esstrem"  # Essential tremor
    , "esstreif"
    , "downs"    # Down syndrome
    , "downsif"
    , "hunt"     # Huntington's disease
    , "huntif"
    , "prion"    # Prion disease (CJD, other)
    , "prionif"
    , "brninj"   # Traumatic brain injury
    , "brninjif"
    , "hyceph"   # Normal-pressure hydrocephalus
    , "hycephif"
    , "epilep"   # Epilepsy
    , "epilepif"
    , "neop"     # CNS neoplasm
    , "neopif" 
    , "hiv"      # Human immunodeficiency virus (HIV)
    , "hivif"
  ) %>% c(., paste0("fu_", .), paste0("tele_", .))
# Form M1 Milestone data
fields_u3_m1_raw <-
  c(
    "m1_visitmo"
    , "m1_visitday"
    , "m1_visityr"
    , "note_mlstn_type" # 0 no further contact; 1 continued contact
    , "protocol" # 1 telephone; 2 minimal; 3 in-person
    , "deceased" # 1 yes
    , "discont"  # 1 yes
  )
# Combine and collapse `fields_u3_*_raw` vectors
fields_u3_raw <- 
  c(
    fields_u3_hd_raw
    , fields_u3_a1_raw
    , fields_u3_d1_raw
    , fields_u3_m1_raw
  )
fields_u3 <- fields_u3_raw %>% paste(collapse = ",")


# _ _ UMMAP General ----

# Header data
fields_ug_head_raw <- 
  c(
    "subject_id"           # partic. ID
    , "redcap_event_name"  # REN
    , "exam_date"          # visit date
  )
# Research data
fields_ug_res_raw <-
  c(
    # "blood_drawn"          # research
    "consent_to_autopsy" # research
    # , "mri_completed"      # research
    # , "sample_given"       # research
    # , "comp_withd"         # research
  )
# Combine and collapse `fields_ms_*_raw` vectors
fields_ug_raw <-
  c(
    fields_ug_head_raw
    , fields_ug_res_raw
  )
fields_ug <- 
  fields_ug_raw %>% paste(collapse = ",")


# _ _ MiNDSet Registry ----

# Header data
fields_ms_head_raw <- 
  c(
    "subject_id"      # partic. ID
    , "date"          # RVF date
  )
# Demographic data
fields_ms_dem_raw <-
  c(
    "race_value"           # demo
  )
fields_ms_rvf_raw <-
  c(
    "owndoc"
    , "madc_website"
    , "madcnewsletter"
    , "alz_assoc"
    , "radio_announc"
    , "tv"
    , "event"
    , "health_fair"
    , "other"
    , "referred_ifother"
  )
# Combine and collapse `fields_ms_*_raw` vectors
fields_ms_raw <-
  c(
    fields_ms_head_raw
    , fields_ms_dem_raw
    , fields_ms_rvf_raw
  )
fields_ms <- fields_ms_raw %>% paste(collapse = ",")


# LOAD DATA ----

# _ JSON Retrieval ---

# # _ _ UDS 3 ---
# if (download_json) {
#   json_u3 <-
#     RCurl::postForm(
#       uri     = REDCAP_API_URI,
#       token   = REDCAP_API_TOKEN_UDS3n,
#       content = "record",
#       format  = "json",
#       type    = "flat",
#       fields  = fields_u3,
#       rawOrLabel             = "raw",
#       rawOrLabelHeaders      = "raw",
#       exportCheckboxLabel    = "false",
#       exportSurveyFields     = "false",
#       exportDataAccessGroups = "false",
#       returnFormat           = "json",
#       # .opts = list(ssl.verifypeer = TRUE, verbose = FALSE) # see note below*
#       .opts = list(ssl.verifypeer = FALSE, verbose = FALSE)
#     ) %>%
#     str_replace_all(pattern = "\r\n", replacement = " ")
# }
# # Convert JSON to tibble; convert "" values to NA
# df_u3 <- jsonlite::fromJSON(json_u3) %>% as_tibble() %>%  na_if("")


# _ _ MiNDSet Registry ---

# # Retrieve JSON object
# if (download_json) {
#   json_ms <-
#     RCurl::postForm(
#       uri     = REDCAP_API_URI,
#       token   = REDCAP_API_TOKEN_MINDSET,
#       content = "record",
#       format  = "json",
#       type    = "flat",
#       fields  = fields_ms,
#       rawOrLabel             = "label",
#       rawOrLabelHeaders      = "raw",
#       exportCheckboxLabel    = "false",
#       exportSurveyFields     = "false",
#       exportDataAccessGroups = "false",
#       returnFormat           = "json",
#       filterLogic            = '([exam_date] >= "2017-03-28")',
#       # .opts = list(ssl.verifypeer = TRUE, verbose = FALSE) # see note below*
#       .opts = list(ssl.verifypeer = FALSE, verbose = FALSE)
#     ) %>%
#     str_replace_all(pattern = "\r\n", replacement = " ")
# }
# # Convert JSON to tibble; convert "" values to NA
# df_ms <- jsonlite::fromJSON(json_ms) %>% as_tibble() %>% na_if("")

# _ Data Dump CSVs ----

# _ _ UDS 3 ----
df_u3 <- 
  read_csv(
    paste0(data_dump_path, "/UMMAPUDS3_DATA_2020-03-15_2252.csv"),
    col_types = cols(.default = col_character()),
    na = ""
  ) %>%
  select(all_of(fields_u3_raw))

# _ _ UMMAP General ----
df_ug <- 
  read_csv(
    paste0(data_dump_path, "/UMMAPGeneral_DATA_2020-03-15_2250.csv"),
    col_types = cols(.default = col_character()),
    na = ""
  ) %>%
  select(all_of(fields_ug_raw))

# _ _ MiNDSet Registry ----
df_ms <- 
  read_csv(
    paste0(data_dump_path, "/MiNDSetRegistry_DATA_2020-03-15_2248.csv"),
    col_types = cols(.default = col_character()),
    na = ""
  ) %>%
  select(all_of(fields_ms_raw))

# PROCESS DATA ----

# _ Clean Data ----

# _ _ UDS 3 ----

df_u3_cln <- df_u3 %>% 
  # deselect useless field(s)
  select(-redcap_event_name) %>%
  # remove records without visit dates
  filter(!is.na(form_date)) %>% 
  # remove double data entry (DDE) records
  filter(str_detect(ptid, pattern = "^UM\\d{8}$")) %>% 
  # # remove milestoned pts.
  # filter(!(ptid %in% df_u3_mlst)) %>% 
  # keep only distinct / non-duplicate rows
  distinct(ptid, form_date, .keep_all = TRUE) %>% 
  # coerce fields to appropriate types
  type_convert(
    col_types = cols(.default = col_integer()
                     , ptid = col_character()
                     , redcap_event_name = col_character()
                     , form_date = col_date()
    ))

# _ _ UDS 3 Milestone ----

df_u3_mlst_cln <- df_u3 %>% 
  # deselect useless field(s)
  select(-redcap_event_name) %>% 
  # select relevant fields
  select(ptid, all_of(fields_u3_m1_raw)) %>%
  # remove double data entry (DDE) records
  filter(str_detect(ptid, pattern = "^UM\\d{8}$")) %>% 
  # filter records missing `note_mlstn_type`
  filter(!is.na(note_mlstn_type)) %>% 
  # Coalesce M1 date
  mutate(
    m1_date = as.Date(paste0(m1_visityr, "-", m1_visitmo, "-", m1_visitday))
  ) %>% 
  select(-m1_visityr, -m1_visitmo, -m1_visitday) %>% 
  # keep only latest milestone form data
  group_by(ptid) %>%
  filter(m1_date == max(m1_date)) %>%
  ungroup() %>%
  # coerce fields to appropriate types
  type_convert(
    col_types = cols(.default = col_integer(),
                     ptid = col_character(),
                     m1_date = col_date())
  )

# _ _ UMMAP General ----

df_ug_cln <- df_ug %>% 
  # deselect useless field(s)
  select(-redcap_event_name) %>% 
  # remove records without visit dates
  filter(!is.na(exam_date) & exam_date >= as.Date("2017-03-15")) %>% 
  # remove non UM ID records
  filter(str_detect(subject_id, pattern = "^UM\\d{8}$")) %>% 
  # keep only distinct / non-duplicate rows
  distinct(subject_id, exam_date, .keep_all = TRUE)

# _ _ MiNDSet Registry ----

df_ms_cln <- df_ms %>% 
  # remove records without RVF dates
  filter(!is.na(date) & date >= as.Date("2015-03-15")) %>%
  # rename `race_value` field to `race`
  rename(race = race_value) %>% 
  # rename `date` to `rvf_date`
  rename(rvf_date = date) %>% 
  # keep only distinct / non-duplicate rows
  distinct(subject_id, rvf_date, .keep_all = TRUE)

# _ Mutate Data ----

# _ _ UDS 3 ----

df_u3_cln_mut <- df_u3_cln %>% 
  # coalesce IVP / FVP / TVP columns
  coalesce_ift_cols() %>% 
  # convert `sex` and `race` fields from int to string
  mutate(
    sex = case_when(
      sex == 1 ~ "Male",
      sex == 2 ~ "Female",
      TRUE ~ NA_character_
    )
  ) %>% 
  # create MADC diagnosis field
  rowwise() %>% 
  # mutate(md_sum = amndem + pca + ppasyn + ftdsyn + lbdsyn + namndem) %>% 
  mutate(md_sum = 
           sum(amndem, pca, ppasyn, ftdsyn, lbdsyn, namndem, na.rm = TRUE)) %>% 
  mutate(madc_dx = case_when(
    sum(amndem, pca, ppasyn, ftdsyn, lbdsyn, namndem, na.rm = TRUE) > 1 ~ 
      "Mixed dementia",
    # note_mlstn_type == 0 & deceased == 1 ~ "Deceased",
    # note_mlstn_type == 0 & discont  == 1 ~ "Dropped",
    normcog == 1 ~ "Normal",
    impnomci == 1 ~ "Impaired not MCI",
    demented == 0 & mciamem  == 1 ~ "MCI",
    demented == 0 & mciaplus == 1 ~ "MCI",
    demented == 0 & mcinon1  == 1 ~ "MCI",
    demented == 0 & mcinon2  == 1 ~ "MCI",
    amndem == 1 & (alzdis == 1 & alzdisif == 1) ~ "AD",
    ftdsyn == 1 ~ "FTD/PPA",
    (amndem == 1 | ppasyn == 1) & 
      ((psp == 1 & pspif == 1) |
         (cort == 1 & cortif == 1) |
         (ftldmo == 1 & ftldmoif == 1) |
         (ftldnos == 1 & ftldnoif == 1)) ~ "FTD/PPA",
    lbdsyn == 1 ~ "LBD",
    sum(amndem, pca, ppasyn, ftdsyn, lbdsyn, namndem, na.rm = TRUE) == 1 ~
      "Other dementia",
    TRUE ~ NA_character_
  )) %>% 
  ungroup() %>% 
  # simplify UDS diagnosis fields
  mutate(uds_dx_der = case_when(
    normcog  == 1 ~ "Normal",
    demented == 1 & amndem == 1 ~ 
      "Amnestic multidomain dementia syndrome",
    # "Dementia",
    demented == 1 & pca == 1 ~
      "Posterior cortical atrophy syndrome",
    # "Dementia",
    demented == 1 & ppasyn == 1 ~
      "Primary progressive aphasia syndrome",
    # "Dementia",
    demented == 1 & ftdsyn == 1 ~
      "Behavioral variant FTD syndrome",
    # "Dementia",
    demented == 1 & lbdsyn == 1 ~
      "Lewy body dementia syndrome",
    # "LBD",
    demented == 1 & namndem == 1 ~
      "Non-amnestic multidomain dementia syndrome",
    # "Dementia",
    demented == 0 & mciamem  == 1 ~ "MCI",
    demented == 0 & mciaplus == 1 ~ "MCI",
    demented == 0 & mcinon1  == 1 ~ "MCI",
    demented == 0 & mcinon2  == 1 ~ "MCI",
    demented == 0 & impnomci == 1 ~ "Impaired not MCI",
    TRUE ~ NA_character_
  )) %>% 
  # simplify UDS etiology fields
  mutate(uds_prim_etio = case_when(
    alzdis   == 1 & alzdisif == 1 ~ "AD",
    lbdis    == 1 & lbdif    == 1 ~ "LBD",
    msa      == 1 & msaif    == 1 ~ "MSA",
    psp      == 1 & pspif    == 1 ~ "PSP",
    cort     == 1 & cortif   == 1 ~ "CBD",
    ftldmo   == 1 & ftldmoif == 1 ~ "FTLD with motor neuron disease",
    ftldnos  == 1 & ftldnoif == 1 ~ "FTLD NOS",
    cvd      == 1 & cvdif    == 1 ~ "Vascular brain injury",
    esstrem  == 1 & esstreif == 1 ~ "Essential tremor",
    downs    == 1 & downsif  == 1 ~ "Down syndrome",
    hunt     == 1 & huntif   == 1 ~ "Huntington's disease",
    prion    == 1 & prionif  == 1 ~ "Prion disease",
    brninj   == 1 & brninjif == 1 ~ "Traumatic injury",
    hyceph   == 1 & hycephif == 1 ~ "Normal-pressure hydrocephalus",
    epilep   == 1 & epilepif == 1 ~ "Epilepsy",
    neop     == 1 & neopif   == 1 ~ "CNS neoplasm",
    hiv      == 1 & hivif    == 1 ~ "HIV",
    TRUE ~ NA_character_
  )) %>% 
  # calculate visit numbers
  calculate_visit_num(ptid, form_date) %>% 
  # deslect now-unnecessary fields
  select(
    -normcog,
    -demented,
    -amndem, 
    -pca, 
    -ppasyn, 
    -ftdsyn, 
    -lbdsyn, 
    -namndem,
    -mciamem,  -mciaplus, 
    -mcinon1,  -mcinon2,
    -impnomci,
    -alzdis,   -alzdisif,
    -lbdis,    -lbdif,
    -msa,      -msaif,
    -psp,      -pspif, 
    -cort,     -cortif,
    -ftldmo,   -ftldmoif, 
    -ftldnos,  -ftldnoif,
    -cvd,      -cvdif,
    -esstrem,  -esstreif,
    -downs,    -downsif,
    -hunt,     -huntif,
    -prion,    -prionif,
    -brninj,   -brninjif,
    -hyceph,   -hycephif,
    -epilep,   -epilepif,
    -neop,     -neopif,
    -hiv,      -hivif
  ) %>%
  # drop `note_mlstn_type`, `deceased`, `discont`, `protocol`
  select(-note_mlstn_type, -deceased, -discont, -protocol,
         -m1_visityr, -m1_visitmo, -m1_visitday)


# _ _ UDS 3 Milestone ----

df_u3_mlst_cln_mut <- 
  df_u3_mlst_cln %>% 
  mutate(mlst_summ = case_when(
    # note_mlstn_type == 1 & protocol == 1 ~ "Telephone Follow-up",
    # note_mlstn_type == 1 & protocol == 2 ~ "Minimal",
    # note_mlstn_type == 1 & protocol == 3 ~ "In-Person",
    note_mlstn_type == 0 & deceased == 1 ~ "Deceased",
    note_mlstn_type == 0 & discont == 1 ~ "Discontinued",
    TRUE ~ NA_character_
  )) %>% 
  select(-note_mlstn_type, -protocol, -deceased, -discont)



# _ _ UMMAP General ----

df_ug_cln_mut <- df_ug_cln %>% 
  # coerce fields to particular data types
  type_convert(
    col_types = cols(.default    = col_character()
                     , exam_date = col_date()
                     , consent_to_autopsy = col_integer()
    )) %>% 
  mutate(consent_to_autopsy = case_when(
    consent_to_autopsy == 1 ~ "Yes",
    consent_to_autopsy == 2 ~ "No",
    consent_to_autopsy == 3 ~ "Considering",
    TRUE ~ NA_character_
  ))


# _ _ MiNDSet Registry ----

df_ms_cln_mut <- df_ms_cln %>% 
  # coerce fields to particular data types
  type_convert(
    col_types = cols(.default = col_integer()
                     , subject_id  = col_character()
                     , rvf_date = col_date()
                     , referred_ifother = col_character()
    )) %>% 
  mutate(race = case_when(
    race == 1 ~ "White",
    race == 2 ~ "Black",
    race == 3 ~ "Asian",
    race == 4 ~ "Hispanic",
    race == 5 ~ "Other",
    race == 6 ~ "Unknown",
    TRUE ~ NA_character_
  ))


# _ Join Data ----

# _ _ UDS 3 + UDS 3 Milestone + UMMAP General + MiNDSet Registry ----
jndf_u3_ug_ms <- 
  df_u3_cln_mut %>%
  left_join(df_u3_mlst_cln_mut,
            by = c("ptid" = "ptid")) %>% 
  left_join(df_ug_cln_mut, 
            by = c("ptid" = "subject_id", 
                   "form_date" = "exam_date")) %>% 
  left_join(df_ms_cln_mut,
            by = c("ptid" = "subject_id")) %>% 
  mutate(madc_dx = case_when(
    !is.na(mlst_summ) ~ mlst_summ,
    TRUE ~ madc_dx
  ))


# Generate CC table for RPPR
get_four_nums <- function(dx, df) {
  total <- df %>% 
    get_visit_n(ptid, form_date, Inf) %>% 
    filter(madc_dx == dx) %>%
    nrow()
  female <- df %>% 
    get_visit_n(ptid, form_date, Inf) %>% 
    filter(madc_dx == dx) %>%
    filter(sex == "Female") %>% 
    nrow()
  minority <- df %>% 
    get_visit_n(ptid, form_date, Inf) %>% 
    filter(madc_dx == dx) %>%
    filter(race != "White") %>% 
    nrow()
  autopsy <- df %>% 
    get_visit_n(ptid, form_date, Inf) %>% 
    filter(madc_dx == dx) %>%
    filter(consent_to_autopsy == "Yes") %>%
    # filter(consent_to_autopsy == 1) %>% 
    nrow()
  tibble(madc_dx = dx,
         total = total,
         female = female,
         minority = minority,
         autopsy = autopsy)
}

uniq_madc_dx <- jndf_u3_ug_ms %>% distinct(madc_dx) %>% pull() %>% c(., NA_character_)

ummap_nums <- purrr::map_df(.x = uniq_madc_dx, 
                            .f = get_four_nums, 
                            jndf_u3_ug_ms)

ummap_nums <-
  ummap_nums[
    match(c("Normal", 
            "Impaired not MCI",
            "MCI",
            "AD",
            "FTD/PPA",
            "LBD",
            "Mixed dementia",
            "Other dementia",
            "Discontinued",
            "Deceased",
            NA_character_), 
          ummap_nums$madc_dx), 
    ]

colSums(ummap_nums[, 2:5])

bind_rows(ummap_nums, colSums(ummap_nums[, 2:5]))

readr::write_csv(bind_rows(ummap_nums, colSums(ummap_nums[, 2:5])),
                 paste0("ummap_numbers_2020-03-15.csv"), na = "Total")


# Generate OR Core table for RPPR
# write_csv(df_u3_ms, "df_u3_ms.csv", na = "")

# "madc_website"
# "madcnewsletter"
# "alz_assoc"
# "radio_announc"
# "tv"
# "event"
# "health_fair"
# "other"
# "referred_ifother"

jndf_u3_ug_ms %>%
  filter(!is.na(madc_website)) %>%
  nrow()
jndf_u3_ug_ms %>%
  filter(!is.na(madcnewsletter)) %>%
  nrow()

ref_sources <-
  c(
    "madc_website"
    , "madcnewsletter"
    , "alz_assoc"
    , "radio_announc"
    , "tv"
    , "event"
    , "health_fair"
  )

calc_ref_source_n <- function(ref_source, df) {
  ref_source <- ensym(ref_source) # string to symbol
  ref_source <- enquo(ref_source) # symbol to quosure
  df %>% 
    # filter(!!ref_source == "Source") %>% 
    filter(!is.na(!!ref_source)) %>% 
    nrow()
}

calc_ref_source_n("madc_website", jndf_u3_ug_ms)
calc_ref_source_n("madcnewsletter", jndf_u3_ug_ms)


vct_ref_sources <-
  purrr::map_int(ref_sources, calc_ref_source_n, jndf_u3_ug_ms)
names(vct_ref_sources) <- ref_sources
tibble::enframe(vct_ref_sources) %>% 
  rename(referral = name, n = value) %>% 
  write_csv("rvf_referral_2020-03-15.csv", na = "")

# purrr::reduce(purrr::map_int(ref_sources, calc_ref_source_n, df_u3_ms), `+`)
# purrr::reduce(
#   purrr::map_int(ref_sources, calc_ref_source_n, jndf_u3_ug_ms), 
#   sum
# )
sum(purrr::map_int(ref_sources, calc_ref_source_n, jndf_u3_ug_ms))

jndf_u3_ug_ms %>% 
  select(madc_website, 
         madcnewsletter, 
         alz_assoc,
         radio_announc,
         tv,
         event, 
         health_fair, 
         other, 
         referred_ifother) %>% 
  filter(is.na(madc_website)     &
           is.na(madcnewsletter) &
           is.na(alz_assoc)      &
           is.na(radio_announc)  &
           is.na(tv)             &
           is.na(event)          &
           is.na(health_fair)) %>%
  filter(!is.na(other)) %>% 
  group_by(referred_ifother) %>%
  summarize(n = n()) %>% 
  arrange(desc(n)) %>%
  replace_na(list(referred_ifother = "[empty]")) %>%
  write_csv("rvf_other_referral_2020-03-15.csv", na = "")


df_ms_cln_mut %>% 
  filter(rvf_date >= as.Date("2019-03-15")) %>% 
  nrow()

df_ms_cln_mut %>% 
  filter(rvf_date >= as.Date("2019-03-15")) %>% 
  filter(race != "White") %>% 
  nrow()










