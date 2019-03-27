# CC_RPPR.R

download_json <- FALSE
source("~/Desktop/config.R")
source("~/Desktop/helpers.R")

# USEFUL LIBRARIES ----
suppressMessages( library(dplyr) )
suppressMessages( library(stringr) )
suppressMessages( library(tidyr) )
suppressMessages( library(readr) )


# Header data
fields_u3_hd_raw <-
  c(
    'ptid'           # partic. ID
    ,'form_date'     # visit date
  )
# Form A1 Demographics
fields_u3_a1_raw <-
  c(
    "sex"            # sex
    # , "race"         # race not available for old cohort; no A1 fields!
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
    "note_mlstn_type"
    , "protocol"
    , "deceased"
    , "discont"
  )
# Combine and collapse `fields_u3_*_raw` vectors
fields_u3 <- c(fields_u3_hd_raw
               , fields_u3_a1_raw
               , fields_u3_d1_raw
               # , fields_u3_d2_raw
               , fields_u3_m1_raw
) %>% paste(collapse = ",")
rm(fields_u3_hd_raw)
rm(fields_u3_a1_raw)
rm(fields_u3_d1_raw)
rm(fields_u3_d2_raw)
rm(fields_u3_m1_raw)

# __ MiNDSet Registry ----

# Header data
fields_ms_head_raw <- 
  c(
    'subject_id'           # partic. ID
    , 'exam_date'          # visit date
  )
# Demographic data
fields_ms_dem_raw <-
  c(
    "race_value"           # demo
    # , 'county'             # demo
    # , 'zip_code'           # demo
  )
# Research data
fields_ms_res_raw <-
  c(
    # 'blood_drawn'          # research
    'consent_to_autopsy' # research
    # , 'mri_completed'      # research
    # , 'sample_given'       # research
    # , 'comp_withd'         # research
  )
# Combine and collapse `fields_ms_*_raw` vectors
fields_ms <- 
  c(
    fields_ms_head_raw
    , fields_ms_dem_raw
    , fields_ms_res_raw
  ) %>% paste(collapse = ",")
rm(fields_ms_head_raw)
rm(fields_ms_dem_raw)
# rm(fields_ms_res_raw)
# rm(fields_ms_time_raw)

# __ UDS 2

# _ Retrieve Data via REDCap API ----

# __ UDS 3

# Retrieve JSON object ----
if (download_json) {
  json_u3 <-
    RCurl::postForm(
      uri     = REDCAP_API_URI,
      token   = REDCAP_API_TOKEN_UDS3n,
      content = 'record',
      format  = 'json',
      type    = 'flat',
      fields  = fields_u3,
      rawOrLabel             = 'raw',
      rawOrLabelHeaders      = 'raw',
      exportCheckboxLabel    = 'false',
      exportSurveyFields     = 'false',
      exportDataAccessGroups = 'false',
      returnFormat           = 'json',
      # .opts = list(ssl.verifypeer = TRUE, verbose = FALSE) # see note below*
      .opts = list(ssl.verifypeer = FALSE, verbose = FALSE)
    ) %>%
    str_replace_all(pattern = "\r\n", replacement = " ")
}
# Convert JSON to tibble; convert "" values to NA
df_u3 <- jsonlite::fromJSON(json_u3) %>% as_tibble() %>%  na_if("")

# __ MiNDSet Registry ----

# Retrieve JSON object
if (download_json) {
  json_ms <-
    RCurl::postForm(
      uri     = REDCAP_API_URI,
      token   = REDCAP_API_TOKEN_MINDSET,
      content = 'record',
      format  = 'json',
      type    = 'flat',
      fields  = fields_ms,
      rawOrLabel             = 'label',
      rawOrLabelHeaders      = 'raw',
      exportCheckboxLabel    = 'false',
      exportSurveyFields     = 'false',
      exportDataAccessGroups = 'false',
      returnFormat           = 'json',
      filterLogic            = '([exam_date] >= "2017-03-28")',
      # .opts = list(ssl.verifypeer = TRUE, verbose = FALSE) # see note below*
      .opts = list(ssl.verifypeer = FALSE, verbose = FALSE)
    ) %>%
    str_replace_all(pattern = "\r\n", replacement = " ")
}
# Convert JSON to tibble; convert "" values to NA
df_ms <- jsonlite::fromJSON(json_ms) %>% as_tibble() %>% na_if("")


# PROCESS DATA ----

# _ Clean Data ----

# __ UDS 3 ----

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
  # coerce field to appropriate type
  type_convert(
    col_types = cols(.default = col_integer()
                     , ptid = col_character()
                     , redcap_event_name = col_character()
                     , form_date = col_date()
    ))

# __ MiNDSet Registry ----

df_ms_cln <- df_ms %>% 
  # deselect useless field(s)
  select(-redcap_event_name) %>% 
  # remove records without visit dates
  filter(!is.na(exam_date)) %>% 
  # remove non UM ID records
  filter(str_detect(subject_id, pattern = "^UM\\d{8}$")) %>% 
  # rename `race_value` field to `race`
  rename(race = race_value) %>% 
  # keep only distinct / non-duplicate rows
  distinct(subject_id, exam_date, .keep_all = TRUE)

# _ Mutate Data ----

# __ UDS 3 ----

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
    note_mlstn_type == 0 & deceased == 1 ~ "Deceased",
    note_mlstn_type == 0 & discont  == 1 ~ "Dropped",
    normcog == 1 ~ "Normal",
    impnomci == 1 ~ "Impaired not MCI",
    demented == 0 & mciamem  == 1 ~ "MCI",
    demented == 0 & mciaplus == 1 ~ "MCI",
    demented == 0 & mcinon1  == 1 ~ "MCI",
    demented == 0 & mcinon2  == 1 ~ "MCI",
    amndem == 1 & (alzdis == 1 & alzdisif == 1) ~ "AD",
    (amndem == 1 | ppasyn == 1) & 
      ((psp == 1 & pspif == 1) |
         (cort == 1 & cortif == 1) |
         (ftldmo == 1 & ftldmoif == 1) |
         (ftldnos == 1 & ftldnoif == 1)) ~ "FTD/PPA",
    lbdsyn == 1 ~ "LBD",
    ftdsyn == 1 ~ "FTD/PPA",
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
  select(-note_mlstn_type, -deceased, -discont, -protocol)  


# __ MiNDSet Registry ----

df_ms_cln_mut <- df_ms_cln %>% 
  # coerce fields to particular data types
  type_convert(
    col_types = cols(.default    = col_character()
                     , exam_date = col_date()
    ))


# _ Join Data ----
df_u3_ms <- left_join(x = df_u3_cln_mut, 
                      y = df_ms_cln_mut, 
                      by = c("ptid" = "subject_id", 
                             "form_date" = "exam_date"))

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
    nrow()
  tibble(madc_dx = dx,
         total = total,
         female = female,
         minority = minority,
         autopsy = autopsy)
}

uniq_madc_dx <- df_u3_ms %>% distinct(madc_dx) %>% pull()

ummap_nums <- purrr::map_df(uniq_madc_dx, get_four_nums, df_u3_ms)

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
            "Dropped",
            "Deceased"), 
          ummap_nums$madc_dx), 
    ]

colSums(ummap_nums[, 2:5])

bind_rows(ummap_nums, colSums(ummap_nums[, 2:5]))

readr::write_csv(bind_rows(ummap_nums, colSums(ummap_nums[, 2:5])),
                 "ummap_numbers.csv", na = "Total")

