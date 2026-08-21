# Run all three datasets in this order. Remove an ID here only if a partial run is needed.
STUDY_IDS <- c("UK", "US", "Pilot")

# Default: reuse the fitted .rds files and regenerate all tables and figures.
# Change to TRUE only when the Bayesian models must be refitted.
REFIT_MODELS <- FALSE

# Supplementary analyses formerly stored in Power and distributions_scripts/.
RUN_DISTRIBUTION_ANALYSIS <- TRUE
RUN_POWER_ANALYSIS <- TRUE

# Exploratory consensus/dispersion analysis requested for the manuscript.
# This module describes dispersion within each direction x relationship x norm
# cell. It does not test whether AI-condition dispersion differs from the
# matched human-condition dispersion.
RUN_DISPERSION_ANALYSIS <- TRUE
DISPERSION_BOOTSTRAP_REPS <- 5000

STUDY_CONFIG <- list(
  UK = list(
    file = "UK_Data.csv",
    display_label = "UK Full",
    output_id = "UK_Full",
    model_id = "UK"
  ),
  US = list(
    file = "US_Data.csv",
    display_label = "USA Full",
    output_id = "USA_Full",
    model_id = "US"
  ),
  Pilot = list(
    file = "Pilot_Data.csv",
    display_label = "USA Pilot",
    output_id = "USA_Pilot",
    model_id = "Pilot"
  )
)

if (!all(STUDY_IDS %in% names(STUDY_CONFIG))) {
  stop("STUDY_IDS may contain only: UK, US, Pilot")
}

run_study <- function(STUDY_ID) {
study_settings  <- STUDY_CONFIG[[STUDY_ID]]
STUDY_FILE      <- study_settings$file
STUDY_LABEL     <- study_settings$display_label
STUDY_OUTPUT_ID <- study_settings$output_id
STUDY_MODEL_ID  <- study_settings$model_id

HUMAN_AI_ID <- "human-ai"
AI_HUMAN_ID <- "ai-human"

# Save all console output from each run ----
RUN_LOG_DIR <- file.path("outputs", "run_logs")
DEMOGRAPHIC_REPORT_DIR <- file.path("outputs", "demographics")
dir.create(RUN_LOG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DEMOGRAPHIC_REPORT_DIR, recursive = TRUE, showWarnings = FALSE)

RUN_LOG_FILE <- file.path(
  RUN_LOG_DIR,
  paste0(
    "analysis_",
    STUDY_OUTPUT_ID,
    "_",
    format(Sys.time(), "%Y%m%d_%H%M%S"),
    ".txt"
  )
)

RUN_LOG_OUTPUT_SINK_LEVEL <- sink.number(type = "output")
RUN_LOG_MESSAGE_SINK_LEVEL <- sink.number(type = "message")
RUN_LOG_CONNECTION <- file(RUN_LOG_FILE, open = "wt")
sink(RUN_LOG_CONNECTION, split = TRUE)
sink(RUN_LOG_CONNECTION, type = "message")

close_run_log <- function() {
  while (sink.number(type = "message") > RUN_LOG_MESSAGE_SINK_LEVEL) {
    current_level <- sink.number(type = "message")
    try(sink(type = "message"), silent = TRUE)
    if (sink.number(type = "message") >= current_level) {
      break
    }
  }
  while (sink.number(type = "output") > RUN_LOG_OUTPUT_SINK_LEVEL) {
    current_level <- sink.number(type = "output")
    try(sink(type = "output"), silent = TRUE)
    if (sink.number(type = "output") >= current_level) {
      break
    }
  }
  if (isOpen(RUN_LOG_CONNECTION)) {
    close(RUN_LOG_CONNECTION)
  }
  invisible(NULL)
}

# Always release the output sinks, including when a study stops with an error.
on.exit(close_run_log(), add = TRUE)

cat("Run log:", normalizePath(RUN_LOG_FILE, mustWork = FALSE), "\n")

library(tidyverse)
library(stringr)
library(lme4)
library(easystats)
library(forcats)
library(boot)
library(here)
library(ggridges)
library(patchwork)
library(brms)
library(withr)

# Output folders ----
FIGURE_DIR <- here::here("outputs", "figures")
KS_EUCLID_FIGURE_DIR <- file.path(FIGURE_DIR, "ks_euclid")
NORM_RIDGE_FIGURE_DIR <- file.path(FIGURE_DIR, "norm_ridges")
SI_FIGURE_DIR <- file.path(FIGURE_DIR, "SI")
BRM_OUTPUT_DIR <- here::here("brm_median_model")

NORM_MEDIAN_DIR <- file.path(BRM_OUTPUT_DIR, "norm_median_model")
NORM_MODEL_DIR  <- file.path(NORM_MEDIAN_DIR, "models")
NORM_MEDIAN_DIFFERENCE_DIR <- file.path(NORM_MEDIAN_DIR, "median_difference")
NORM_MANUAL_CONTRAST_DIR   <- file.path(NORM_MEDIAN_DIR, "manual_contrast")

HUMANLIKENESS_MEDIAN_DIR <- file.path(BRM_OUTPUT_DIR, "humanlikeness_median_model")
HUMANLIKENESS_MODEL_DIR  <- file.path(HUMANLIKENESS_MEDIAN_DIR, "models")
HUMANLIKENESS_REPORT_DIR <- file.path(HUMANLIKENESS_MEDIAN_DIR, "reports")

DEMOGRAPHIC_REPORT_FILE <- file.path(
  DEMOGRAPHIC_REPORT_DIR,
  paste0("demographic_summary_", STUDY_OUTPUT_ID, ".txt")
)

dirs_to_make <- c(
  KS_EUCLID_FIGURE_DIR,
  NORM_RIDGE_FIGURE_DIR,
  SI_FIGURE_DIR,
  NORM_MODEL_DIR,
  NORM_MEDIAN_DIFFERENCE_DIR,
  NORM_MANUAL_CONTRAST_DIR,
  HUMANLIKENESS_MODEL_DIR,
  HUMANLIKENESS_REPORT_DIR
)

purrr::walk(dirs_to_make, ~ dir.create(.x, recursive = TRUE, showWarnings = FALSE))

# Define the processing function to ensure identical cleaning ----
## Clean function ----
clean_study_data <- function(file_name, return_all = FALSE) {
  # Loading dataset
  dataset <- read.csv(here::here("data", "raw", file_name))

  # removing top two rows
  dataset <- dataset[-c(1, 2), ]

  # including only participants who finished the survey
  dataset <- dataset %>% filter(Progress == "100")

  # implementing exclusionary criteria
  # adding new columns
  dataset <- dataset %>%
    mutate(
      exclusions = "INCLUDE",
      exclusion_reasons = ""
    )

  # marking participants who actively agreed to not read carefully
  dataset <- dataset %>%
    mutate(
      exclusions = ifelse(str_detect(ReadCarefully, "I do NOT agree to read carefully, and I will return the survey to Prolific"), "EXCLUDE", exclusions),
      exclusion_reasons = ifelse(str_detect(ReadCarefully, "I do NOT agree to read carefully, and I will return the survey to Prolific"),
        "ReadCarefully",
        exclusion_reasons
      )
    )

  # marking participants who failed CareQuiz2
  dataset <- dataset %>%
    mutate(
      exclusions = ifelse(
        !(CareQuiz == "To ensure that individuals' needs are met without feeling explicitly indebted" |
          CareQuizTake2 == "To ensure that individuals' needs are met without feeling explicitly indebted"),
        "EXCLUDE", exclusions
      ),
      exclusion_reasons = ifelse(
        !(CareQuiz == "To ensure that individuals' needs are met without feeling explicitly indebted" |
          CareQuizTake2 == "To ensure that individuals' needs are met without feeling explicitly indebted"),
        ifelse(exclusion_reasons == "", "CareQuiz, CareQuizTake2", paste(exclusion_reasons, "CareQuiz, CareQuizTake2", sep = ", ")),
        exclusion_reasons
      )
    )

  # marking participants who failed TransQuizTake2
  dataset <- dataset %>%
    mutate(
      exclusions = ifelse(
        !(TransQuiz == "To keep track of benefits given or received to ensure fair exchanges" |
          TransQuizTake2 == "To keep track of benefits given or received to ensure fair exchanges"),
        "EXCLUDE", exclusions
      ),
      exclusion_reasons = ifelse(
        !(TransQuiz == "To keep track of benefits given or received to ensure fair exchanges" |
          TransQuizTake2 == "To keep track of benefits given or received to ensure fair exchanges"),
        ifelse(exclusion_reasons == "", "TransQuiz, TransQuizTake2", paste(exclusion_reasons, "TransQuiz, TransQuizTake2", sep = ", ")),
        exclusion_reasons
      )
    )

  # marking participants who failed HierQuizTake2
  dataset <- dataset %>%
    mutate(
      exclusions = ifelse(
        !(HierQuiz == "To coordinate behavior in situations where individuals have unequal authority" |
          HierQuizTake2 == "To coordinate behavior in situations where individuals have unequal authority"),
        "EXCLUDE", exclusions
      ),
      exclusion_reasons = ifelse(
        !(HierQuiz == "To coordinate behavior in situations where individuals have unequal authority" |
          HierQuizTake2 == "To coordinate behavior in situations where individuals have unequal authority"),
        ifelse(exclusion_reasons == "", "HierQuiz, HierQuizTake2", paste(exclusion_reasons, "HierQuiz, HierQuizTake2", sep = ", ")),
        exclusion_reasons
      )
    )

  # marking participants who failed MateQuizTatke2 (please note the spelling error in this item specifically)
  dataset <- dataset %>%
    mutate(
      exclusions = ifelse(
        !(MateQuiz == "To establish or maintain a sexual or romantic connection" |
          MateQuizTatke2 == "To establish or maintain a sexual or romantic connection"),
        "EXCLUDE", exclusions
      ),
      exclusion_reasons = ifelse(
        !(MateQuiz == "To establish or maintain a sexual or romantic connection" |
          MateQuizTatke2 == "To establish or maintain a sexual or romantic connection"),
        ifelse(exclusion_reasons == "", "MateQuiz, MateQuizTatke2", paste(exclusion_reasons, "MateQuiz, MateQuizTatke2", sep = ", ")),
        exclusion_reasons
      )
    )

  # marking participants who failed HumanCheck_1 (only for participants in the "Human" condition)
  dataset <- dataset %>%
    mutate(
      exclusions = ifelse(
        Condition == "Human" &
          !(HumanCheck_1 == "They will be more intelligent than humans, be able to improve themselves, and will try to help humans" |
            HumanCheck_2 == "They will be more intelligent than humans, be able to improve themselves, and will try to help humans"),
        "EXCLUDE", exclusions
      ),
      exclusion_reasons = ifelse(
        Condition == "Human" &
          !(HumanCheck_1 == "They will be more intelligent than humans, be able to improve themselves, and will try to help humans" |
            HumanCheck_2 == "They will be more intelligent than humans, be able to improve themselves, and will try to help humans"),
        ifelse(exclusion_reasons == "", "HumanCheck_1, HumanCheck_2", paste(exclusion_reasons, "HumanCheck_1, HumanCheck_2", sep = ", ")),
        exclusion_reasons
      )
    )

  # marking participants who failed AICheck_2 (only for participants in the "AI" condition)
  dataset <- dataset %>%
    mutate(
      exclusions = ifelse(
        Condition == "AI" &
          !(AICheck_1 == "They will be more intelligent than humans, be able to improve themselves, and will try to help humans" |
            AICheck_2 == "They will be more intelligent than humans, be able to improve themselves, and will try to help humans"),
        "EXCLUDE", exclusions
      ),
      exclusion_reasons = ifelse(
        Condition == "AI" &
          !(AICheck_1 == "They will be more intelligent than humans, be able to improve themselves, and will try to help humans" |
            AICheck_2 == "They will be more intelligent than humans, be able to improve themselves, and will try to help humans"),
        ifelse(exclusion_reasons == "", "AI_Check1, AI_Check2", paste(exclusion_reasons, "AICheck_1, AICheck_2", sep = ", ")),
        exclusion_reasons
      )
    )

  # marking participants who failed Careful_Human (only for participants in the "Human" condition)
  dataset <- dataset %>%
    mutate(
      exclusions = ifelse(Condition == "Human" & Careful_Human != "I do not understand", "EXCLUDE", exclusions),
      exclusion_reasons = ifelse(Condition == "Human" & Careful_Human != "I do not understand",
        ifelse(exclusion_reasons == "", "Careful_Human", paste(exclusion_reasons, "Careful_Human", sep = ", ")),
        exclusion_reasons
      )
    )

  # marking participants who failed Careful_AI (only for participants in the "AI" condition)
  dataset <- dataset %>%
    mutate(
      exclusions = ifelse(Condition == "AI" & Careful_AI != "I do not understand", "EXCLUDE", exclusions),
      exclusion_reasons = ifelse(Condition == "AI" & Careful_AI != "I do not understand",
        ifelse(exclusion_reasons == "", "Careful_AI", paste(exclusion_reasons, "Careful_AI", sep = ", ")),
        exclusion_reasons
      )
    )

  # marking participants who failed InRealityCheck (only for participants in the "AI" condition)
  dataset <- dataset %>%
    mutate(
      exclusions = ifelse(Condition == "AI" & InRealityCheck_AI != "Think about how humans and superintelligent AI should ideally behave toward one another, regardless of whether the AI will have actual needs or can be benefitted in any way", "EXCLUDE", exclusions),
      exclusion_reasons = ifelse(Condition == "AI" & InRealityCheck_AI != "Think about how humans and superintelligent AI should ideally behave toward one another, regardless of whether the AI will have actual needs or can be benefitted in any way",
        ifelse(exclusion_reasons == "", "InRealityCheck", paste(exclusion_reasons, "InRealityCheck_AI", sep = ", ")),
        exclusion_reasons
      )
    )

  # marking participants who failed Comp_check (Correct answer = "Human-human" in Human condition, "Human and superintelligent AI" in AI condition)
  dataset <- dataset %>%
    mutate(
      exclusions = ifelse(
        (Condition == "AI" & Comp_check != "Human and superintelligent AI") |
          (Condition == "Human" & Comp_check != "Human-human"),
        "EXCLUDE", exclusions
      ),
      exclusion_reasons = ifelse(
        (Condition == "AI" & Comp_check != "Human and superintelligent AI") |
          (Condition == "Human" & Comp_check != "Human-human"),
        ifelse(exclusion_reasons == "", "Comp_check", paste(exclusion_reasons, "Comp_check", sep = ", ")),
        exclusion_reasons
      )
    )

  # marking participants who failed Bot.check.Friday and Bot.check.Friday.2
  dataset <- dataset %>%
    mutate(
      exclusions = ifelse(!(str_detect(toupper(Bot.check.Friday), "FRIDAY") | str_detect(toupper(Bot.check.Friday.2), "FRIDAY")), "EXCLUDE", exclusions),
      exclusion_reasons = ifelse(!(str_detect(toupper(Bot.check.Friday), "FRIDAY") | str_detect(toupper(Bot.check.Friday.2), "FRIDAY")),
        ifelse(exclusion_reasons == "", "Bot.check.Friday, Bot.check.Friday.2", paste(exclusion_reasons, "Bot.check.Friday, Bot.check.Friday.2", sep = ", ")),
        exclusion_reasons
      )
    )

  if (return_all) {
    return(dataset)
  }

  # excluding exclusions
  clean <- filter(dataset, exclusions == "INCLUDE")

  return(clean)
}


# Clean data
dat_clean <- clean_study_data(STUDY_FILE)
dat_clean_for_distributions <- dat_clean
dat_exclusion_all <- clean_study_data(STUDY_FILE, return_all = TRUE)

sample_size_summary <- tibble(
  metric = c("Raw N", "Exclusion N", "Final N"),
  n = c(
    nrow(dat_exclusion_all),
    sum(dat_exclusion_all$exclusions == "EXCLUDE", na.rm = TRUE),
    nrow(dat_clean)
  )
)





## Rename and recode function ----
process_dataset <- function(df) {
  # Initial Column Renaming
  df <- df %>%
    # Renaming columns to include relationship information
    # 1_ = romantic_ ; 2_ = close friends ; 3_ = workers ; 4_ = supervisor - assistant ;
    # 5_ = teacher - student; 6_ = mental health provider - patient; 8_ = customer - seller
    rename_with(~ str_replace_all(., c(
      "^X1_" = "romantic_", "^X2_" = "friends_", "^X3_" = "workers_",
      "^X4_" = "supervisor_", "^X5_" = "teacher_", "^X6_" = "therapist_", "^X8_" = "seller_"
    ))) %>%
    # _1 = romantic ; _2 = close friends ; _3 = workers ; _4 = supervisor - assistant ;
    # _5 = teacher - student; _6 = mental health provider - patient; _7 = customer - seller
    rename_with(~ str_replace_all(., c(
      "AI_Anxiety_1" = "anxiety_romantic", "AI_Anxiety_2" = "anxiety_friends",
      "AI_Anxiety_3" = "anxiety_workers", "AI_Anxiety_4" = "anxiety_supervisor",
      "AI_Anxiety_5" = "anxiety_teacher", "AI_Anxiety_6" = "anxiety_therapist",
      "AI_Anxiety_7" = "anxiety_seller"
    ))) %>%
    rename_with(~ str_replace_all(., c(
      "Believability_1$" = "realness_romantic", "Believability_2$" = "realness_friends",
      "Believability_3$" = "realness_workers", "Believability_4$" = "realness_supervisor",
      "Believability_5$" = "realness_teacher", "Believability_6$" = "realness_therapist",
      "Believability_7$" = "realness_seller"
    ))) %>%
    rename_with(~ str_replace_all(., c(
      "MindQualities_1" = "subjective experience", "MindQualities_2$" = "sentience",
      "MindQualities_3$" = "self-consciousness", "MindQualities_4$" = "well-being",
      "MindQualities_5$" = "autonomy", "MindQualities_6$" = "beneficence",
      "MindQualities_7$" = "cognition"
    )))

  # Standardization of Political Column Names
  if ("Political...Economic_1" %in% colnames(df)) {
    df <- df %>% rename(`Political - Economic_1` = `Political...Economic_1`)
  }
  if ("Political...Social_1" %in% colnames(df)) {
    df <- df %>% rename(`Political - Social_1` = `Political...Social_1`)
  }

  # Specific Likert Cleaning
  df <- df %>%
    # Clean Mind Qualities
    pivot_longer(cols = `subjective experience`:cognition, names_to = "variable", values_to = "value") %>%
    mutate(value = case_when(
      value == "Strongly disagree" ~ "1",
      value == "Somewhat disagree" ~ "2",
      value == "Neither agree nor disagree" ~ "3",
      value == "Somewhat agree" ~ "4",
      value == "Strongly agree" ~ "5",
      TRUE ~ value
    )) %>%
    pivot_wider(names_from = variable, values_from = value) %>%
    # Clean Behavioral Intent
    pivot_longer(cols = romantic_Care_Human1:seller_Transaction_AI, names_to = "variable", values_to = "value") %>%
    mutate(value = case_when(
      value == "Definitely should not\n-4\n" ~ "-4",
      value == "Neutral\n0\n" ~ "0",
      value == "Neutral\n0" ~ "0",
      value == "Definitely should\n4\n" ~ "4",
      TRUE ~ value
    )) %>%
    pivot_wider(names_from = variable, values_from = value) %>%
    # Clean Anxiety
    pivot_longer(cols = anxiety_romantic:anxiety_seller, names_to = "variable", values_to = "value") %>%
    mutate(value = case_when(
      value == "Not at all afraid\n1\n" ~ "1",
      value == "Extremely afraid\n10\n" ~ "10",
      TRUE ~ value
    )) %>%
    pivot_wider(names_from = variable, values_from = value) %>%
    # Clean Realness
    pivot_longer(cols = realness_romantic:realness_seller, names_to = "variable", values_to = "value") %>%
    mutate(value = case_when(
      value == "Not at all real\n1" ~ "1",
      value == "Completely real\n10" ~ "10",
      TRUE ~ value
    )) %>%
    pivot_wider(names_from = variable, values_from = value)

  # Convert cleaned scales back to numeric
  df <- df %>%
    mutate(across(
      `subjective experience`:`realness_seller`,
      ~ suppressWarnings(as.numeric(as.character(.)))
    )) %>%
    mutate(
      Age = as.numeric(Age),
      Gender = as.factor(Gender),
      Ethn = as.factor(Ethn),
      Relat.Status = as.factor(Relat.Status),
      Education = as.factor(Education),
      StatusLadder = as.numeric(StatusLadder),
      `Political - Economic_1` = as.numeric(`Political - Economic_1`),
      `Political - Social_1` = as.numeric(`Political - Social_1`),
      How.Religious_1 = as.numeric(How.Religious_1),
      LLM_Know = as.factor(LLM_Know),
      LLM_Past = as.factor(LLM_Past),
      LLM_Freq = as.factor(LLM_Freq)
    )

  return(df)
}

dat_clean <- process_dataset(dat_clean)



## Demographic report function ----
demographic_report <- function(df, label = "Dataset") {
  cat("Report for:", label, "\n")

  # Age
  cat("\nAge: \n")
  age_sum <- df %>% summarise(mean = mean(Age, na.rm = TRUE), SD = sd(Age, na.rm = TRUE))
  print(age_sum)
  cat("Range:", min(df$Age, na.rm = TRUE), "-", max(df$Age, na.rm = TRUE), "\n")

  # Gender
  cat("\nGender: \n")
  print(table(df$Gender))

  # Ethnicity
  cat("\nEthnicity: \n")
  ethnicity_combined <- df %>%
    mutate(Ethn = case_when(
      str_detect(Ethn, "Asian") ~ "Asian",
      str_detect(Ethn, "Black/African") ~ "Black/African",
      str_detect(Ethn, "Caucasian/White") ~ "Caucasian/White",
      str_detect(Ethn, "Hispanic/Latinx") ~ "Hispanic/Latinx",
      str_detect(Ethn, "Native American") ~ "Native American",
      str_detect(Ethn, "Pacific Islander") ~ "Pacific Islander",
      str_detect(Ethn, "Other") ~ "Other",
      str_detect(Ethn, "Prefer not to say") ~ "Prefer not to say",
      TRUE ~ NA_character_
    )) %>%
    filter(!is.na(Ethn)) %>%
    count(Ethn) %>%
    mutate(proportion = n / sum(n))
  print(ethnicity_combined)

  # Relationship Status
  cat("\nRelationship Status: \n")
  relat_summary <- df %>%
    count(Relat.Status) %>%
    mutate(proportion = n / sum(n))
  print(relat_summary)

  # Education
  cat("\nEducation: \n")
  edu_summary <- df %>%
    count(Education) %>%
    mutate(proportion = n / sum(n))
  # Apply level renaming
  if (nrow(edu_summary) >= 5) {
    levels(edu_summary$Education)[5] <- "Postgraduate/Professional degree\nor other advanced degree"
  }
  print(edu_summary)

  # Status Ladder
  cat("\nStatus Ladder: \n")
  ladder_sum <- df %>%
    summarise(mean = mean(StatusLadder, na.rm = TRUE), SD = sd(StatusLadder, na.rm = TRUE))
  print(ladder_sum)
  cat("Range:", min(df$StatusLadder, na.rm = TRUE), "-", max(df$StatusLadder, na.rm = TRUE), "\n")

  ladder_summary <- df %>%
    count(StatusLadder) %>%
    mutate(proportion = n / sum(n))
  print(ladder_summary)

  # Political Orientation
  cat("\nPolitical Orientation: \n")
  pol_long <- df %>%
    select(`Political - Economic_1`, `Political - Social_1`) %>%
    pivot_longer(cols = everything(), names_to = "type", values_to = "Score")

  pol_sum <- pol_long %>%
    group_by(type) %>%
    summarise(mean = mean(Score, na.rm = TRUE), sd = sd(Score, na.rm = TRUE))
  print(pol_sum)

  # LLM Metrics
  cat("\nLLM Metrics: \n")
  df %>%
    count(LLM_Know) %>%
    mutate(proportion = n / sum(n)) %>%
    print()
  df %>%
    count(LLM_Past) %>%
    mutate(proportion = n / sum(n)) %>%
    print()
  df %>%
    mutate(LLM_Freq = forcats::fct_na_value_to_level(na_if(LLM_Freq, ""), level = "No prior LLM use / not applicable")) %>%
    count(LLM_Freq) %>%
    mutate(proportion = n / sum(n)) %>%
    print()
}

write_demographic_report <- function(df, label, file_path, sample_size = NULL) {
  report_lines <- capture.output({
    if (!is.null(sample_size)) {
      cat("Sample size summary:\n")
      print(sample_size)
      cat("\n")
    }

    demographic_report(df, label)
  })

  cat(paste(report_lines, collapse = "\n"), "\n")
  writeLines(report_lines, file_path)
  cat("\nDemographic report saved to:", normalizePath(file_path, mustWork = FALSE), "\n")
}

write_demographic_report(
  dat_clean,
  STUDY_LABEL,
  DEMOGRAPHIC_REPORT_FILE,
  sample_size = sample_size_summary
)

glimpse(dat_clean)


# Color and label ----
fill_light <- c(
  "Care"        = "#F4C2C2",   # 浅粉
  "Hierarchy"   = "#D3D3D3",   # 浅灰
  "Mating"      = "#B0E57C",   # 浅绿
  "Transaction" = "#ADD8E6"    # 浅蓝
)

fill_dark <- c(
  "Care"        = "#C71585",
  "Hierarchy"   = "#696969",
  "Mating"      = "#228B22",
  "Transaction" = "#1E90FF"
)

fill_euclid <- c(
  "Care"        = "#F4C2C2",
  "Hierarchy"   = "#D3D3D3",
  "Mating"      = "#B0E57C",
  "Transaction" = "#ADD8E6",
  "Euclidean"   = "#262e36"
)

multicolor_theme <- theme_minimal() +
  theme(
    axis.text.x       = element_text(size = 13),
    axis.text.y       = element_text(size = 13),
    plot.title        = element_text(size = 15),
    panel.grid        = element_blank(),
    axis.ticks.y      = element_line(),
    axis.ticks.x      = element_line(),
    axis.ticks.length = unit(0.3, "lines"),
    axis.line.y       = element_line(),
    axis.line.x       = element_line()
  )

relationship_labels <- c(
  "teacher"    = "Teacher",
  "therapist"  = "Mental health provider",
  "supervisor" = "Assistant",
  "seller"     = "Seller",
  "workers"    = "Coworker",
  "friends"    = "Friend",
  "romantic"   = "Romantic\npartner"
)


# Euclidean functions ----
# compute euclidean value and use them in the K-S heatmaps
compute_euclid_summary <- function(plotting_df, dir_a, dir_b) {
  means_df <- plotting_df %>%
    filter(!is.na(value), direction %in% c(dir_a, dir_b)) %>%
    mutate(value = as.numeric(value)) %>%
    group_by(category, relationship, direction) %>%
    summarise(value = mean(value, na.rm = TRUE), .groups = "drop")

  distances <- expand_grid(
    means_df,
    means_df %>% rename(
      category.1 = category, relationship.1 = relationship,
      direction.1 = direction, value.1 = value
    )
  ) %>%
    filter(
      category == category.1,
      relationship == relationship.1,
      direction  == dir_a,
      direction.1 == dir_b
    ) %>%
    mutate(sq_dist = (value - value.1)^2)

  distances %>%
    filter(!is.na(sq_dist)) %>%
    group_by(relationship) %>%
    summarise(euclid = round(sqrt(sum(sq_dist)), 2), .groups = "drop")
}

build_euclid_augmented <- function(ks_results_df, euclid_summary) {
  euclid_col <- euclid_summary %>%
    rename(KS_Statistic = euclid) %>%
    mutate(cooperative_function = "Euclidean", P_Value = NA)

  euclid_order <- euclid_summary %>% arrange(desc(euclid)) %>% pull(relationship)

  df <- bind_rows(ks_results_df, euclid_col) %>%
    mutate(
      relationship         = factor(relationship, levels = euclid_order),
      cooperative_function = factor(cooperative_function,
                                    levels = c("Care", "Hierarchy", "Mating", "Transaction", "Euclidean"))
    )

  df <- df %>%
    group_by(is_euclid = (cooperative_function == "Euclidean")) %>%
    mutate(
      alpha_val  = scales::rescale(KS_Statistic, to = c(0.3, 1)),
      text_color = ifelse(is_euclid & alpha_val > 0.65, "white", "black"),
      font_face  = ifelse(!is_euclid & !is.na(fdr) & fdr < 0.05, "bold", "plain")
    ) %>%
    ungroup()

  return(df)
}


plot_euclid_multicolor <- function(augmented_df, title_text, subtitle_text) {
  ggplot(augmented_df, aes(x = cooperative_function, y = relationship, fill = cooperative_function)) +
    geom_tile(aes(alpha = I(alpha_val)), color = "white") +
    geom_text(aes(label = sprintf("%.2f", KS_Statistic),
                  color = I(text_color), fontface = I(font_face)), size = 5) +
    scale_fill_manual(values = fill_euclid) +
    scale_y_discrete(labels = relationship_labels) +
    labs(x = "", y = "", title = title_text, subtitle = subtitle_text) +
    multicolor_theme +
    theme(legend.position = "none")
}

format_heatmap_p_label <- function(p_value) {
  case_when(
    is.na(p_value) ~ "",
    p_value < 0.001 ~ "p < .001",
    TRUE ~ paste0("p = ", sub("^0", "", sprintf("%.3f", p_value)))
  )
}

format_heatmap_p_parenthetical <- function(p_value) {
  label <- format_heatmap_p_label(p_value)
  ifelse(label == "", "", paste0("(", label, ")"))
}

format_heatmap_p_plotmath <- function(p_value) {
  case_when(
    is.na(p_value) ~ "''",
    p_value < 0.001 ~ 'paste("(", italic(p), " < .001)")',
    TRUE ~ paste0('paste("(", italic(p), " = ', sub("^0", "", sprintf("%.3f", p_value)), ')")')
  )
}

plot_euclid_multicolor_si <- function(augmented_df, title_text, subtitle_text) {
  if (!"fdr" %in% names(augmented_df)) {
    augmented_df$fdr <- augmented_df$P_Value
  }

  plot_df <- augmented_df %>%
    mutate(
      is_euclid = cooperative_function == "Euclidean",
      is_significant = !is_euclid & !is.na(fdr) & fdr < 0.05,
      ks_label = sprintf("%.2f", KS_Statistic),
      p_label = format_heatmap_p_plotmath(fdr)
    )

  ggplot(plot_df, aes(x = cooperative_function, y = relationship, fill = cooperative_function)) +
    geom_tile(aes(alpha = I(alpha_val)), color = "white") +
    geom_text(
      data = filter(plot_df, is_euclid),
      aes(x = cooperative_function, y = relationship, label = ks_label, color = I(text_color)),
      size = 4.2,
      inherit.aes = FALSE
    ) +
    geom_text(
      data = filter(plot_df, !is_euclid),
      aes(x = cooperative_function, y = relationship, label = ks_label,
          color = I(text_color), fontface = ifelse(is_significant, "bold", "plain")),
      size = 4.2,
      vjust = -0.45,
      inherit.aes = FALSE
    ) +
    geom_text(
      data = filter(plot_df, !is_euclid),
      aes(x = cooperative_function, y = relationship, label = p_label, color = I(text_color)),
      size = 3.5,
      vjust = 1.35,
      parse = TRUE,
      inherit.aes = FALSE
    ) +
    scale_fill_manual(values = fill_euclid) +
    scale_y_discrete(labels = relationship_labels) +
    labs(x = "", y = "", title = title_text, subtitle = subtitle_text) +
    multicolor_theme +
    theme(legend.position = "none")
}




# Norm ridges plots data ----
dat_long <- dat_clean %>%
  select(ResponseId, romantic_Care_Human1:seller_Transaction_AI) %>%
  pivot_longer(
    cols      = romantic_Care_Human1:seller_Transaction_AI,
    names_to  = "variable",
    values_to = "value"
  ) %>%
  filter(!is.na(value)) %>%
  separate(variable, into = c("relationship", "category", "direction"), sep = "_")

dat_long <- dat_long %>%
  mutate(
    Condition = case_when(
      direction == "Human1" ~ "Human condition",
      direction == "Human2" ~ "Human condition",
      direction == "AI"     ~ "Superintelligent AI condition",
      direction == "Human"  ~ "Superintelligent AI condition",
      TRUE ~ direction
    ),
    role = case_when(
      direction == "Human1" ~ HUMAN_AI_ID,
      direction == "Human2" ~ AI_HUMAN_ID,
      direction == "Human"  ~ HUMAN_AI_ID,
      direction == "AI"     ~ AI_HUMAN_ID,
      TRUE ~ direction
    )
  )

dat_long

boot_median_ci <- function(x, R = 5000) {
  x <- x[!is.na(x)]
  
  b <- boot::boot(
    data = x,
    statistic = function(data, indices) {
      median(data[indices])
    },
    R = R
  )
  
  qs <- quantile(
    b$t[, 1],
    probs = c(0.025, 0.975),
    na.rm = TRUE,
    names = FALSE
  )
  
  tibble(
    median_value = median(x),
    lower_ci = qs[1],
    upper_ci = qs[2]
  )
}

ci_data_all <- dat_long %>%
  group_by(category, Condition, role, relationship) %>%
  group_nest() %>%
  arrange(category, Condition, role, relationship) %>%
  mutate(
    group_seed = 123 + row_number(),
    result = map2(
      data,
      group_seed,
      ~ with_seed(.y, boot_median_ci(.x$value, R = 5000))
    )
  ) %>%
  select(-data, -group_seed) %>%
  unnest(result)


# Human likeness data ----
{
  dat_clean_v2 <- dat_clean %>%
    mutate(across(romantic_Care_Human1:seller_Transaction_AI, as.numeric)) %>%
    pivot_longer(
      cols      = romantic_Care_Human1:seller_Transaction_AI,
      names_to  = "variable", values_to = "value"
    ) %>%
    filter(!is.na(value)) %>%
    separate(variable, into = c("relationship", "category", "direction"), sep = "_")

  dat_clean_v2 <- dat_clean_v2 %>%
    mutate(
      Condition = case_when(
        direction == "Human1" ~ "Human condition",
        direction == "Human2" ~ "Human condition",
        direction == "AI"     ~ "Superintelligent AI condition",
        direction == "Human"  ~ "Superintelligent AI condition",
        TRUE ~ direction
      ),
      role = case_when(
        direction == "Human1" ~ HUMAN_AI_ID,
        direction == "Human2" ~ AI_HUMAN_ID,
        direction == "Human"  ~ HUMAN_AI_ID,
        direction == "AI"     ~ AI_HUMAN_ID,
        TRUE ~ direction
      ),
      direction = case_when(
        direction == "Human1" ~ "Human partner 1 to human partner 2",
        direction == "Human2" ~ "Human partner 2 to human partner 1",
        direction == "Human"  ~ "Human partner to AI partner",
        direction == "AI"     ~ "AI partner to human partner",
        TRUE ~ direction
      )
    )

  dat_clean_v2$LLM_Past <- as.factor(dat_clean_v2$LLM_Past)

  dat_clean_v2 <- dat_clean_v2 %>% select(-(1:8), -(10:36))

  dat_clean_v2 <- dat_clean_v2 %>%
    mutate(
      LLM_Freq = ifelse((is.na(LLM_Freq) | LLM_Freq == "") & LLM_Past == "No",
                        "Never", as.character(LLM_Freq)),
      LLM_Freq = fct_infreq(factor(LLM_Freq))
    )

  dat_clean_v2$LLM_Freq <- factor(dat_clean_v2$LLM_Freq, levels = c(
    "Never", "About once a month", "About once every two weeks",
    "About once a week", "About once a day", "More than once a day"
  ))

  # remain only AI condition
  dat_v2_d2 <- dat_clean_v2 %>% filter(role == AI_HUMAN_ID, Condition != "Human condition")
  dat_v2_d1 <- dat_clean_v2 %>% filter(role == HUMAN_AI_ID, Condition != "Human condition")

  dat_v2_d2 <- dat_v2_d2 %>% mutate(across(`subjective experience`:realness_seller, as.numeric))
  dat_v2_d1 <- dat_v2_d1 %>% mutate(across(`subjective experience`:realness_seller, as.numeric))

  dat_v2_d2$category     <- as.factor(dat_v2_d2$category)
  dat_v2_d2$relationship <- as.factor(dat_v2_d2$relationship)
  dat_v2_d1$category     <- as.factor(dat_v2_d1$category)
  dat_v2_d1$relationship <- as.factor(dat_v2_d1$relationship)
}



# AI-human: How should social AIs treat humans? ----

## K-S Heatmaps ----
ks_AI <- function(data) {
  # Get the column names
  column_names <- names(data)

  # Extract unique X_Y combinations by removing the final component (Z)
  xy_combinations <- unique(sub("_([^_]+)$", "", column_names))

  # Initialize an empty list to store the K-S test results
  ks_results <- list()

  # Loop over each X_Y combination
  for (xy in xy_combinations) {
    # Identify the relevant Human2 and AI columns
    col_Human2 <- grep(paste0("^", xy, "_Human2$"), column_names, value = TRUE)
    col_AI <- grep(paste0("^", xy, "_AI$"), column_names, value = TRUE)

    # Only run the test if both columns exist
    if (length(col_Human2) == 1 && length(col_AI) == 1) {
      human2_data <- data[[col_Human2]]
      ai_data <- data[[col_AI]]

      # Perform the K-S test between Human2 and AI data
      ks_test <- ks.test(human2_data, ai_data)

      # Store the results
      ks_results[[length(ks_results) + 1]] <- data.frame(
        XY = xy,
        KS_Statistic = ks_test$statistic,
        P_Value = ks_test$p.value
      )
    }
  }

  # Return results
  if (length(ks_results) > 0) {
    ks_df <- do.call(rbind, ks_results)
    return(ks_df)
  } else {
    message("No valid Human2-AI column pairs found.")
    return(data.frame()) # Empty dataframe if nothing matched
  }
}


ks_d2_results <- ks_AI(dat_clean)
print(ks_d2_results)

ks_d2_results <- ks_d2_results %>%
  tidyr::separate(XY, into = c("relationship", "cooperative_function"), sep = "_", remove = FALSE)

# Add Benjamini-Hochberg p-values
ks_d2_results$fdr <- p.adjust(ks_d2_results$P_Value, method = "BH")
write.csv(
  ks_d2_results,
  file.path(SI_FIGURE_DIR, paste0("ks_results_", AI_HUMAN_ID, "_", STUDY_OUTPUT_ID, "_SI.csv")),
  row.names = FALSE
)

### Euclidean and K-S Heatmaps ----
euclid_sum_d2   <- compute_euclid_summary(dat_long, "Human2", "AI")
euclid_aug_d2   <- build_euclid_augmented(ks_d2_results, euclid_sum_d2)

p_euclid_d2 <- plot_euclid_multicolor(
  euclid_aug_d2,
  "How should social AIs treat humans?",
  paste0("Kolmogorov-Smirnov and Euclidean distances in normative expectations\nSample: ",
         STUDY_LABEL)
)
p_euclid_d2
ggsave(file.path(KS_EUCLID_FIGURE_DIR, paste0("ks_euclid_heatmap_", AI_HUMAN_ID, "_", STUDY_OUTPUT_ID, ".pdf")),
       p_euclid_d2, width = 8, height = 6.4, units = "in")

p_euclid_d2_si <- plot_euclid_multicolor_si(
  euclid_aug_d2,
  "How should social AIs treat humans?",
  paste0("Kolmogorov-Smirnov and Euclidean distances in normative expectations\nSample: ",
         STUDY_LABEL)
)
p_euclid_d2_si
ggsave(file.path(SI_FIGURE_DIR, paste0("ks_euclid_heatmap_", AI_HUMAN_ID, "_", STUDY_OUTPUT_ID, "_SI.pdf")),
       p_euclid_d2_si, width = 8, height = 6.4, units = "in")






## Norm ridges plots ----
# Filter AI-human data
dat_long_d2 <- dat_long    %>% filter(role == AI_HUMAN_ID)
ci_d2        <- ci_data_all %>% filter(role == AI_HUMAN_ID)

# Original 7-panel layout with text enlarged by exactly 10%
ORIGINAL_FONT_SCALE <- 1.10

build_norm_ridges_original_font110_legend <- function(
  data,
  title_text,
  subtitle_text = NULL
) {
  ggplot(data, aes(y = fct_relevel(category, "Transaction", "Mating", "Hierarchy", "Care"))) +
    annotate(
      "text", x = -4, y = 1, label = "Human in role", hjust = 0,
      size = 6 * ORIGINAL_FONT_SCALE, color = "black"
    ) +
    annotate("point", x = -3.45,  y = 1, color = fill_light[1], alpha = 0.4, shape = 15, size = 10) +
    annotate("point", x = -3.35,  y = 1, color = fill_light[2], alpha = 0.4, shape = 15, size = 10) +
    annotate("point", x = -3.25,  y = 1, color = fill_light[3], alpha = 0.4, shape = 15, size = 10) +
    annotate("point", x = -3.15,  y = 1, color = fill_light[4], alpha = 0.4, shape = 15, size = 10) +
    annotate(
      "text", x = -2.8, y = 1, label = "AI in role", hjust = 0,
      size = 6 * ORIGINAL_FONT_SCALE, color = "black"
    ) +
    annotate("point", x = -2.42, y = 1, color = fill_light[1], alpha = 0.9, shape = 15, size = 10) +
    annotate("point", x = -2.32, y = 1, color = fill_light[2], alpha = 0.9, shape = 15, size = 10) +
    annotate("point", x = -2.22, y = 1, color = fill_light[3], alpha = 0.9, shape = 15, size = 10) +
    annotate("point", x = -2.12, y = 1, color = fill_light[4], alpha = 0.9, shape = 15, size = 10) +
    labs(title = title_text, subtitle = subtitle_text, x = NULL, y = NULL) +
    scale_x_continuous(limits = c(-4, 4), expand = c(0, 0)) +
    coord_cartesian(clip = "off") +
    theme_minimal() +
    theme(
      plot.title.position = "panel",
      plot.title = element_text(
        size = 20 * ORIGINAL_FONT_SCALE,
        face = "bold",
        margin = margin(b = 4, t = 10)
      ),
      plot.subtitle = element_text(
        size = 16 * ORIGINAL_FONT_SCALE,
        margin = margin(b = 12)
      ),
      panel.grid = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_text(color = "transparent"),
      axis.ticks = element_blank(),
      plot.margin = margin(t = 0, r = 5, b = 0, l = 5)
    )
}

save_norm_ridges_original_font110 <- function(
  data,
  original_main_plot,
  direction_id,
  title_text
) {
  legend_plot <- build_norm_ridges_original_font110_legend(
    data = data,
    title_text = title_text,
    subtitle_text = paste0("Sample: ", STUDY_LABEL)
  )

  main_plot <- original_main_plot +
    theme(
      plot.title = element_text(
        size = 18 * ORIGINAL_FONT_SCALE,
        face = "bold"
      ),
      axis.text.x = element_text(size = 16 * ORIGINAL_FONT_SCALE),
      axis.text.y = element_text(size = 16 * ORIGINAL_FONT_SCALE),
      strip.text = element_text(
        size = 18 * ORIGINAL_FONT_SCALE,
        face = "bold"
      )
    )

  final_plot <- (legend_plot / main_plot) +
    plot_layout(heights = c(0.16, 1))

  ggsave(
    file.path(
      NORM_RIDGE_FIGURE_DIR,
      paste0("norm_ridges_", direction_id, "_", STUDY_OUTPUT_ID, ".pdf")
    ),
    final_plot,
    width = 30, height = 6, units = "in"
  )

  invisible(final_plot)
}

# norm ridges plots
p_ridges_d2 <- ggplot(dat_long_d2,
    aes(x = value, y = fct_relevel(category, "Transaction", "Mating", "Hierarchy", "Care"))) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +

  # Human-Human (top half, lighter fill)
  geom_density_ridges(
    data   = filter(dat_long_d2, Condition == "Human condition"),
    aes(fill = category, group = interaction(category, relationship)),
    scale  = 0.5, alpha = 0.4, rel_min_height = 0.01, color = NA, show.legend = FALSE
  ) +

  # Human-AI (bottom half, darker fill)
  geom_density_ridges(
    data   = filter(dat_long_d2, Condition == "Superintelligent AI condition"),
    aes(fill = category, group = interaction(category, relationship)),
    scale  = -0.5, alpha = 0.9, rel_min_height = 0.01, color = NA, show.legend = FALSE
  ) +

  # median
  geom_point(
    data = ci_d2,
    aes(x = median_value, y = category, color = category),
    size = 3, inherit.aes = FALSE,
    position = position_nudge(y = ifelse(ci_d2$Condition == "Human condition", 0.1, -0.1))
  ) +

  # errorbars
  geom_errorbarh(
    data = ci_d2,
    aes(xmin = lower_ci, xmax = upper_ci, y = category, color = category),
    height = 0.15, size = 0.6, inherit.aes = FALSE,
    position = position_nudge(y = ifelse(ci_d2$Condition == "Human condition", 0.1, -0.1))
  ) +
  facet_wrap(~relationship, ncol = 7,
    labeller = as_labeller(c(
      "teacher"   = "Teacher",    "friends"   = "Friend",
      "romantic"  = "Romantic partner",       "seller"    = "Seller",
      "supervisor"= "Assistant",  "therapist" = "Mental health provider",
      "workers"   = "Coworker"
    ))
  ) +
  labs(x = "", y = "") +
  theme_classic(base_size = 18) +
  theme(
    plot.title       = element_text(size = 18, face = "bold"),
    axis.text.x      = element_text(size = 16),
    axis.text.y      = element_text(size = 16),
    strip.text       = element_text(size = 18, face = "bold"),
    strip.background = element_rect(color = "black", fill = NA, size = 0.8),
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.ticks.x     = element_blank(),
    legend.position  = "none"
  ) +
  scale_fill_manual(values  = fill_light) +
  scale_color_manual(values = c(
    "Care"        = "#C71585",
    "Hierarchy"   = "#696969",
    "Mating"      = "#228B22",
    "Transaction" = "#3277bd"
  )) +
  scale_x_continuous(
    breaks = c(-4, 0, 4),
    labels = c("Definitely\nshould not", "Neutral", "Definitely\nshould")
  )

save_norm_ridges_original_font110(
  data = dat_long_d2,
  original_main_plot = p_ridges_d2,
  direction_id = AI_HUMAN_ID,
  title_text = "How should social AIs treat humans?"
)

## Bayesian model using median value ----
library(brms)
library(dplyr)
library(purrr)
library(tidyr)
library(posterior)
library(coda)

# raw difference of median
dat_long_d2_diff <- dat_long_d2 %>%
  group_by(category, relationship, Condition) %>%
  summarise(raw_median = median(value), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = Condition, values_from = raw_median) %>%
  mutate(diff_AI_minus_Human = `Superintelligent AI condition` - `Human condition`) %>%
  arrange(relationship, category)

write.csv(dat_long_d2_diff,
          file.path(NORM_MEDIAN_DIFFERENCE_DIR, paste0("raw_median_difference_", AI_HUMAN_ID, "_", STUDY_OUTPUT_ID, ".csv")),
          row.names = FALSE)


# set factor level
dat_long_d2$Condition <- factor(
  dat_long_d2$Condition,
  levels = c("Human condition", "Superintelligent AI condition")
)
dat_long_d2$category <- factor(
  dat_long_d2$category,
  levels = c("Care", "Hierarchy", "Mating", "Transaction")
)
dat_long_d2$relationship <- factor(
  dat_long_d2$relationship,
  levels = c("friends", "romantic", "seller", "supervisor", "teacher", "therapist", "workers")
)

levels(dat_long_d2$Condition)
levels(dat_long_d2$category)
levels(dat_long_d2$relationship)


### Model ----
if (REFIT_MODELS) {
  library(cmdstanr)
  options(brms.backend = "cmdstanr")
  cmdstan_path()
  cmdstan_version()

  model_d2 <- brm(
    formula = bf(value ~ Condition * category * relationship + (1 | ResponseId),
                 quantile = 0.5),
    data    = dat_long_d2,
    family  = asym_laplace(),
    chains  = 4, iter = 4000, warmup = 2000, cores = 4,
    control = list(max_treedepth = 15),
    seed    = 123
  )
  summary(model_d2)

  saveRDS(model_d2, file = file.path(NORM_MODEL_DIR, paste0("model_median_brms_", AI_HUMAN_ID, "_", STUDY_MODEL_ID, ".rds")))
}

### Read model ----
model_d2 <- readRDS(file.path(NORM_MODEL_DIR, paste0("model_median_brms_", AI_HUMAN_ID, "_", STUDY_MODEL_ID, ".rds")))
summary(model_d2)



### Manual posterior contrasts ----

# Helper: weighted median
weighted_median <- function(x, w) {
  ord <- order(x)
  x <- x[ord]
  w <- w[ord]
  cw <- cumsum(w) / sum(w)
  x[which(cw >= 0.5)[1]]
}

# Cell counts / weights within each category
# These weights are used for weighted marginal medians
rel_weights_d2 <- dat_long_d2 %>%
  count(category, relationship) %>%
  group_by(category) %>%
  mutate(weight = n / sum(n)) %>%
  ungroup()

rel_weights_d2


# 1. Compare Superintelligent AI vs Human within each category x relationship
compare_one_cell <- function(model, data, cat_value, rel_value) {
  
  nd <- data.frame(
    Condition = factor(
      c("Superintelligent AI condition", "Human condition"),
      levels = levels(data$Condition)
    ),
    category = factor(
      c(cat_value, cat_value),
      levels = levels(data$category)
    ),
    relationship = factor(
      c(rel_value, rel_value),
      levels = levels(data$relationship)
    ),
    ResponseId = NA
  )
  
  lin <- posterior_linpred(
    model,
    newdata = nd,
    re_formula = NA
  )
  
  # Contrast: Superintelligent AI - Human
  diff_draws <- lin[, 1] - lin[, 2]

  # using rowMeans of median to get the marginal one

  hpd <- HPDinterval(as.mcmc(diff_draws), prob = 0.95)

  data.frame(
    category = cat_value,
    relationship = rel_value,
    contrast = "Superintelligent AI - Human",
    estimate = median(diff_draws),
    ci_lower = quantile(diff_draws, 0.025),
    ci_upper = quantile(diff_draws, 0.975),
    lower_HPD = hpd[1, "lower"],
    upper_HPD = hpd[1, "upper"],
    prob_gt_0 = mean(diff_draws > 0),
    prob_lt_0 = mean(diff_draws < 0)
  )
}


# run all cell-specific contrasts
res_cellwise_d2 <- expand.grid(
  category     = levels(dat_long_d2$category),
  relationship = levels(dat_long_d2$relationship),
  stringsAsFactors = FALSE
) %>%
  split(seq_len(nrow(.))) %>%
  map_dfr(~ compare_one_cell(model_d2, dat_long_d2, .x$category, .x$relationship))


# 2. Weighted marginal median contrast within each category
compare_one_category_weighted_median <- function(model, data, cat_value, weight_df) {
  sub_w      <- weight_df %>%
    filter(category == cat_value) %>%
    arrange(match(relationship, levels(data$relationship)))
  rel_levels <- as.character(sub_w$relationship)
  w          <- sub_w$weight

  make_nd <- function(cond) data.frame(
    Condition    = factor(rep(cond, length(rel_levels)), levels = levels(data$Condition)),
    category     = factor(rep(cat_value, length(rel_levels)), levels = levels(data$category)),
    relationship = factor(rel_levels, levels = levels(data$relationship)),
    ResponseId   = NA
  )

  lin_human <- posterior_linpred(model, newdata = make_nd("Human condition"),          re_formula = NA)
  lin_ai    <- posterior_linpred(model, newdata = make_nd("Superintelligent AI condition"), re_formula = NA)

  human_wmed <- apply(lin_human, 1, weighted_median, w = w)
  ai_wmed    <- apply(lin_ai,    1, weighted_median, w = w)
  diff_draws <- ai_wmed - human_wmed
  hpd        <- HPDinterval(as.mcmc(diff_draws), prob = 0.95)

  data.frame(
    category   = cat_value, contrast = "Superintelligent AI - Human",
    estimate   = median(diff_draws),
    ci_lower   = quantile(diff_draws, 0.025), ci_upper = quantile(diff_draws, 0.975),
    lower_HPD  = hpd[1, "lower"], upper_HPD = hpd[1, "upper"],
    prob_gt_0  = mean(diff_draws > 0), prob_lt_0 = mean(diff_draws < 0)
  )
}

# Run weighted marginal median contrasts for each category
res_categorywise_d2 <- map_dfr(
  levels(dat_long_d2$category),
  ~ compare_one_category_weighted_median(model_d2, dat_long_d2, .x, rel_weights_d2)
)

# 3. tidy results
round_results <- function(df) {
  df %>% mutate(across(c(estimate, ci_lower, ci_upper, lower_HPD, upper_HPD,
                         prob_gt_0, prob_lt_0), ~ round(.x, 4)))
}
res_cellwise_d2     <- round_results(res_cellwise_d2)
res_categorywise_d2 <- round_results(res_categorywise_d2)

res_cellwise_d2
res_categorywise_d2

write.csv(res_cellwise_d2,
  file.path(NORM_MANUAL_CONTRAST_DIR, paste0("res_cellwise_", AI_HUMAN_ID, "_", STUDY_OUTPUT_ID, ".csv")),
  row.names = FALSE)
write.csv(res_categorywise_d2,
  file.path(NORM_MANUAL_CONTRAST_DIR, paste0("res_categorywise_", AI_HUMAN_ID, "_", STUDY_OUTPUT_ID, ".csv")),
  row.names = FALSE)





## Human likeness ----
# This section constructs a human-likeness outcome by reusing the fitted norm model as a human normative reference. 
# Specifically, for each category × relationship cell, the posterior median under the Human condition is extracted, 
# treating it as the model-based human benchmark for that cell. these human posterior medians are then joined back to 
# the row-level dataset and define human_likeness as the absolute distance between each observed value and the 
# corresponding human posterior median. As a result, smaller values indicate responses that are closer to the human 
# normative benchmark, whereas larger values indicate responses that are farther away. In this sense, human_likeness 
# is not a direct subjective rating of "how human-like" a response is, but rather a model-based measure of deviation 
# from the human normative standard.
# 
# human_likeness_signed is defined, defined as the signed difference
# between each observed AI value and the human posterior median. Positive values indicate AI responses above the
# human normative benchmark, whereas negative values indicate AI responses below the benchmark. This signed-difference
# analysis is intended as a supplementary analysis, while the absolute-difference model remains the primary analysis.
# 
# The final human-likeness model then uses human_likeness as the dependent variable and examines whether variation in this distance can be explained by 
# relationship, perceived realness, anxiety, other perception-related variables, LLM-use variables, and demographic covariates. 
# Compared with the earlier norm model, the focus here is no longer on comparing Human and AI across Condition × 
# category × relationship cells, but on explaining why some AI responses are closer to, or farther from, the human 
# normative benchmark. Before fitting the model, several stabilization steps are applied, including collapsing sparse 
# factor levels, standardizing continuous predictors, adding weakly informative priors, removing LLM_Past because of 
# its overlap with LLM_Freq, and using QR decomposition to reduce estimation difficulties caused by possible 
# multicollinearity among predictors. These steps are intended to improve model stability, identifiability, and 
# interpretability rather than to change the substantive question.

# get human posterior medians
get_cell_median <- function(model, data, cat_value, rel_value, condition_value) {
  nd  <- data.frame(
    Condition    = factor(condition_value, levels = levels(data$Condition)),
    category     = factor(cat_value,       levels = levels(data$category)),
    relationship = factor(rel_value,       levels = levels(data$relationship)),
    ResponseId   = NA
  )
  lin <- posterior_linpred(model, newdata = nd, re_formula = NA)
  data.frame(Condition = condition_value, category = cat_value, relationship = rel_value,
             posterior_median = median(lin[, 1]))
}


{
  # Extract Human condition posterior medians for each cell

  human_cell_medians_d2 <- expand.grid(
    category     = levels(dat_long_d2$category),
    relationship = levels(dat_long_d2$relationship),
    stringsAsFactors = FALSE
  ) %>%
    split(seq_len(nrow(.))) %>%
    map_dfr(~ get_cell_median(model_d2, dat_long_d2, .x$category, .x$relationship,
                              "Human condition"))

  glimpse(human_cell_medians_d2)

  # human_likeness_abs    = magnitude of deviation from the human normative benchmark
  # human_likeness_signed = signed deviation from the human normative benchmark
  human_median_join_d2 <- human_cell_medians_d2 %>%
    select(category, relationship, human_posterior_median = posterior_median)

  dat_v2_d2 <- dat_v2_d2 %>%
    left_join(human_median_join_d2, by = c("relationship", "category")) %>%
    mutate(
      human_likeness_signed = as.numeric(value) - human_posterior_median,
      human_likeness_abs    = abs(human_likeness_signed),
      human_likeness        = human_likeness_abs
    )

  dat_v2_d2$human_likeness
  glimpse(dat_v2_d2)

  # Match realness and anxiety to corresponding relationship
  dat_v2_d2 <- dat_v2_d2 %>%
    mutate(
      realness_matched = case_when(
        relationship == "friends"    ~ as.numeric(realness_friends),
        relationship == "romantic"   ~ as.numeric(realness_romantic),
        relationship == "workers"    ~ as.numeric(realness_workers),
        relationship == "supervisor" ~ as.numeric(realness_supervisor),
        relationship == "teacher"    ~ as.numeric(realness_teacher),
        relationship == "therapist"  ~ as.numeric(realness_therapist),
        relationship == "seller"     ~ as.numeric(realness_seller)
      ),
      anxiety_matched = case_when(
        relationship == "friends"    ~ as.numeric(anxiety_friends),
        relationship == "romantic"   ~ as.numeric(anxiety_romantic),
        relationship == "workers"    ~ as.numeric(anxiety_workers),
        relationship == "supervisor" ~ as.numeric(anxiety_supervisor),
        relationship == "teacher"    ~ as.numeric(anxiety_teacher),
        relationship == "therapist"  ~ as.numeric(anxiety_therapist),
        relationship == "seller"     ~ as.numeric(anxiety_seller)
      )
    )
}


### Model ----
dat_v2_d2 <- dat_v2_d2 %>%
  rename(
    subjective_experience = `subjective experience`,
    self_consciousness    = `self-consciousness`,
    well_being            = `well-being`,
    Political_Social      = `Political - Social_1`,
    Political_Economic    = `Political - Economic_1`
  )

options(brms.backend = "cmdstanr")

if (REFIT_MODELS) {

# Version 1: unadjusted model with all predictors
model_humanlikeness_d2_v1_rel <- brm(
  formula = bf(
    human_likeness ~ 
      relationship + # Add relationship
      realness_matched + anxiety_matched +
      subjective_experience + sentience + self_consciousness +
      well_being + autonomy + beneficence + cognition +
      LLM_Past + LLM_Freq + Gender + Age + Ethn +
      Relat.Status + Education + StatusLadder +
      Political_Social + Political_Economic + How.Religious_1 +
      (1 | ResponseId),
    quantile = 0.5
  ),
  data    = dat_v2_d2,
  family  = asym_laplace(),
  chains  = 4,
  iter    = 4000,
  warmup  = 2000,
  cores   = 4,
  control = list(max_treedepth = 10),
  seed    = 123
)

summary(model_humanlikeness_d2_v1_rel)
saveRDS(model_humanlikeness_d2_v1_rel,
        file = file.path(HUMANLIKENESS_MODEL_DIR, paste0("model_humanlikeness_", AI_HUMAN_ID, "_relationship_", STUDY_MODEL_ID, "_v1.rds")))


# Version 2: primary model for absolute-difference human likeness

# applying several stabilization steps, including collapsing sparse factor levels, 
# standardizing continuous predictors, adding weakly informative priors, removing LLM_Past because of 
# its overlap with LLM_Freq, and using QR decomposition to reduce estimation difficulties caused by possible 
# multicollinearity among predictors. These steps are intended to improve model stability, identifiability, and 
# interpretability rather than to change the substantive question.

# Note:
# min_n = 20 (lumping each option with fewer than 20 people selecting it)
# In the long format data, there are 12 lines per participant, so it was min_n: 20 * 12 = 240
lump_small <- function(x, min_n = 240, other_level = "Other_small_n") {
  forcats::fct_lump_min(factor(x), min = min_n, other_level = other_level)
}

dat_v2_d2_model <- dat_v2_d2 %>%
  mutate(
    ResponseId   = factor(ResponseId),
    
    Gender         = lump_small(Gender, min_n = 240),
    Ethn           = lump_small(Ethn, min_n = 240),
    `Relat.Status` = lump_small(`Relat.Status`, min_n = 240),
    Education      = lump_small(Education, min_n = 240),
    

    relationship = factor(
      relationship,
      levels = c("friends", "romantic", "seller", "supervisor", "teacher", "therapist", "workers")
    ),
    LLM_Freq = forcats::fct_relevel(
      factor(LLM_Freq),
      "Never"
    ),
    
    # Added: standardization of continuous predictors
    # Notes:
    # The use of prior(normal(0, 1), class = "b") as a weakly informative prior
    # generally assumes that predictors have been standardized
    # (mean = 0, SD = 1).
    # Without standardization, predictors may differ substantially in scale
    # (e.g., Age may range from 18 to 80, whereas Likert-type scales range from 1 to 7),
    # so the same normal(0, 1) prior would impose very different levels of constraint
    # across variables, undermining its interpretation as "weakly informative."
    # Here, scale() is used to perform z-score standardization, and as.numeric()
    # is applied to remove the matrix attribute returned by scale().
    
    Age              = as.numeric(scale(Age)),
    StatusLadder     = as.numeric(scale(StatusLadder)),
    How.Religious_1  = as.numeric(scale(How.Religious_1)),
    Political_Social   = as.numeric(scale(Political_Social)),
    Political_Economic = as.numeric(scale(Political_Economic)),
    
    subjective_experience = as.numeric(scale(subjective_experience)),
    sentience             = as.numeric(scale(sentience)),
    self_consciousness    = as.numeric(scale(self_consciousness)),
    well_being            = as.numeric(scale(well_being)),
    autonomy              = as.numeric(scale(autonomy)),
    beneficence           = as.numeric(scale(beneficence)),
    cognition             = as.numeric(scale(cognition)),
    
    realness_matched      = as.numeric(scale(realness_matched)),
    anxiety_matched       = as.numeric(scale(anxiety_matched))
  )

# Weakly informative priors are specified here instead of relying on the default flat priors.
# In brms, the default prior for population-level effects (class = "b") is an
# improper flat prior if no prior is explicitly set. Once prior(normal(0, 1), class = "b")
# is added, those coefficients are no longer assigned an improper flat prior, but a proper
# normal prior.
#
# Stan documentation generally recommends using informative or at least weakly informative
# priors in regression models, as they often improve estimation stability and sampling behavior.
#
# The priors used here provide conservative regularization:
# - b: normal(0, 1)       -> mildly shrinks regression coefficients toward zero,
#                            reducing weak identification and excessive drift
# - sd: exponential(1)    -> regularizes random-effect standard deviations on the positive scale
# - sigma: exponential(1) -> regularizes the residual scale on the positive scale
#
# If these priors later appear too restrictive, the prior for b can be relaxed
# to normal(0, 2).

priors_hl <- c(
  prior(normal(0, 1), class = "b"),
  prior(exponential(1), class = "sd"),
  prior(exponential(1), class = "sigma")
)


# Model

# Note 1: removing LLM_Past
# In the present results, LLM_PastYes and several LLM_Freq levels showed nearly mirrored
# large positive and negative coefficients, together with poor R-hat values. This pattern
# strongly suggests substantial information overlap between the two variables, leading to
# weak parameter identification. In this situation, it is more reasonable to retain the
# more fine-grained LLM_Freq variable and remove the coarser LLM_Past measure.
# 
# Note 2: decomp = "QR"
# The brms documentation states that decomp = "QR" can be helpful when fitting models
# with highly correlated predictors. Here, multiple perception-related variables
# are likely to be substantially correlated, so QR decomposition is enabled.
#
# Note 3: max_treedepth
# max_treedepth is not increased further here. According to the Stan diagnostics
# documentation, when max_treedepth warnings occur together with high R-hat values,
# the primary concern should be model geometry or parameter identification problems,
# rather than simply increasing max_treedepth.
# The main situation in which raising max_treedepth is justified is when R-hat is already
# acceptable (< 1.01) and ESS is sufficient, but "max treedepth exceeded" warnings still
# occur repeatedly. These warnings were not observed in the present model.

model_humanlikeness_d2_v2_rel <- brm(
  formula = bf(
    human_likeness ~
      relationship + 
      realness_matched + anxiety_matched +
      subjective_experience + sentience + self_consciousness +
      well_being + autonomy + beneficence + cognition +
      LLM_Freq + Gender + Age + Ethn +
      Relat.Status + Education + StatusLadder +
      Political_Social + Political_Economic + How.Religious_1 +
      (1 | ResponseId),
    quantile = 0.5,
    decomp   = "QR"
  ),
  data    = dat_v2_d2_model,
  family  = asym_laplace(),
  prior   = priors_hl,
  backend = "cmdstanr",
  chains  = 4,
  iter    = 4000,
  warmup  = 2000,
  cores   = 4,
  control = list(max_treedepth = 10),
  seed    = 123
)

summary(model_humanlikeness_d2_v2_rel)
saveRDS(
  model_humanlikeness_d2_v2_rel,
  file = file.path(HUMANLIKENESS_MODEL_DIR, paste0("model_humanlikeness_", AI_HUMAN_ID, "_relationship_", STUDY_MODEL_ID, "_v2.rds"))
)

# Use the relationship model as the primary model for the reporting code below.
model_humanlikeness_d2 <- model_humanlikeness_d2_v2_rel


# Supplementary signed-difference analysis
# This model asks whether the same predictors explain the direction of deviation from the
# human normative benchmark, rather than only the magnitude of deviation.
model_humanlikeness_signed_d2_v2_rel <- brm(
  formula = bf(
    human_likeness_signed ~
      relationship +
      realness_matched + anxiety_matched +
      subjective_experience + sentience + self_consciousness +
      well_being + autonomy + beneficence + cognition +
      LLM_Freq + Gender + Age + Ethn +
      Relat.Status + Education + StatusLadder +
      Political_Social + Political_Economic + How.Religious_1 +
      (1 | ResponseId),
    quantile = 0.5,
    decomp   = "QR"
  ),
  data    = dat_v2_d2_model,
  family  = asym_laplace(),
  prior   = priors_hl,
  backend = "cmdstanr",
  chains  = 4,
  iter    = 4000,
  warmup  = 2000,
  cores   = 4,
  control = list(max_treedepth = 10),
  seed    = 123
)

summary(model_humanlikeness_signed_d2_v2_rel)
saveRDS(
  model_humanlikeness_signed_d2_v2_rel,
  file = file.path(HUMANLIKENESS_MODEL_DIR, paste0("model_humanlikeness_signed_", AI_HUMAN_ID, "_relationship_", STUDY_MODEL_ID, "_v2_supplement.rds"))
)

model_humanlikeness_signed_d2 <- model_humanlikeness_signed_d2_v2_rel
} else {
  model_humanlikeness_d2 <- readRDS(
    file.path(HUMANLIKENESS_MODEL_DIR, paste0("model_humanlikeness_", AI_HUMAN_ID, "_relationship_", STUDY_MODEL_ID, "_v2.rds"))
  )
  model_humanlikeness_signed_d2 <- readRDS(
    file.path(HUMANLIKENESS_MODEL_DIR, paste0("model_humanlikeness_signed_", AI_HUMAN_ID, "_relationship_", STUDY_MODEL_ID, "_v2_supplement.rds"))
  )
}


summary(model_humanlikeness_d2)

library(dplyr)
library(brms)
library(posterior)
library(tibble)
library(purrr)
library(coda)

fixef(model_humanlikeness_d2)

fe_d2 <- as.data.frame(
  fixef(model_humanlikeness_d2,
        probs  = c(0.025, 0.975),
        robust = TRUE)
) %>%
  rownames_to_column("term") %>%
  rename(
    estimate  = Estimate,
    mad       = Est.Error, # median absolute deviation
    eti_lower = Q2.5,
    eti_upper = Q97.5
  )

fe_d2

draws_d2 <- as_draws_df(model_humanlikeness_d2)
b_cols_d2 <- grep("^b_", names(draws_d2), value = TRUE)

coef_extra_d2 <- map_dfr(b_cols_d2, function(v) {
  x <- draws_d2[[v]]

  hpd <- HPDinterval(as.mcmc(x), prob = 0.95)

  data.frame(
    term       = sub("^b_", "", v),
    hpdi_lower = unname(hpd[1, "lower"]),
    hpdi_upper = unname(hpd[1, "upper"]),
    prob_gt_0  = mean(x > 0),
    prob_lt_0  = mean(x < 0),
    pd         = max(mean(x > 0), mean(x < 0))
  )
})

coef_extra_d2


diag_d2 <- summarise_draws(
  draws_d2[, b_cols_d2],
  rhat,
  ess_bulk,
  ess_tail
) %>%
  mutate(term = sub("^b_", "", variable)) %>%
  select(term, rhat, ess_bulk, ess_tail)

diag_d2


report_table_d2 <- fe_d2 %>%
  left_join(coef_extra_d2, by = "term") %>%
  left_join(diag_d2,       by = "term") %>%
  mutate(across(where(is.numeric), ~ round(.x, 4)))

report_table_d2

write.csv(
  report_table_d2,
  file.path(HUMANLIKENESS_REPORT_DIR, paste0("humanlikeness_report_", AI_HUMAN_ID, "_relationship_", STUDY_OUTPUT_ID, "_v2.csv")),
  row.names = FALSE
)


# Supplementary report for the signed-difference model
fixef(model_humanlikeness_signed_d2)

fe_signed_d2 <- as.data.frame(
  fixef(model_humanlikeness_signed_d2,
        probs  = c(0.025, 0.975),
        robust = TRUE)
) %>%
  rownames_to_column("term") %>%
  rename(
    estimate  = Estimate,
    mad       = Est.Error,
    eti_lower = Q2.5,
    eti_upper = Q97.5
  )

draws_signed_d2 <- as_draws_df(model_humanlikeness_signed_d2)
b_cols_signed_d2 <- grep("^b_", names(draws_signed_d2), value = TRUE)

coef_extra_signed_d2 <- map_dfr(b_cols_signed_d2, function(v) {
  x <- draws_signed_d2[[v]]

  hpd <- HPDinterval(as.mcmc(x), prob = 0.95)

  data.frame(
    term       = sub("^b_", "", v),
    hpdi_lower = unname(hpd[1, "lower"]),
    hpdi_upper = unname(hpd[1, "upper"]),
    prob_gt_0  = mean(x > 0),
    prob_lt_0  = mean(x < 0),
    pd         = max(mean(x > 0), mean(x < 0))
  )
})

diag_signed_d2 <- summarise_draws(
  draws_signed_d2[, b_cols_signed_d2],
  rhat,
  ess_bulk,
  ess_tail
) %>%
  mutate(term = sub("^b_", "", variable)) %>%
  select(term, rhat, ess_bulk, ess_tail)

report_table_signed_d2 <- fe_signed_d2 %>%
  left_join(coef_extra_signed_d2, by = "term") %>%
  left_join(diag_signed_d2,       by = "term") %>%
  mutate(across(where(is.numeric), ~ round(.x, 4)))

report_table_signed_d2

write.csv(
  report_table_signed_d2,
  file.path(HUMANLIKENESS_REPORT_DIR, paste0("humanlikeness_signed_report_", AI_HUMAN_ID, "_relationship_", STUDY_OUTPUT_ID, "_v2_supplement.csv")),
  row.names = FALSE
)



# Human-AI: How should humans treat social AIs? ----

## K-S Heatmaps ----
ks_Human1_vs_Human <- function(data) {
  column_names    <- names(data)
  xy_combinations <- unique(sub("_([^_]+)$", "", column_names))
  ks_results      <- list()

  for (xy in xy_combinations) {
    col_Human1 <- grep(paste0("^", xy, "_Human1$"), column_names, value = TRUE)
    col_Human  <- grep(paste0("^", xy, "_Human$"),  column_names, value = TRUE)

    if (length(col_Human1) == 1 && length(col_Human) == 1) {
      ks_test <- ks.test(data[[col_Human1]], data[[col_Human]])
      ks_results[[length(ks_results) + 1]] <- data.frame(
        XY           = xy,
        KS_Statistic = ks_test$statistic,
        P_Value      = ks_test$p.value
      )
    }
  }

  if (length(ks_results) > 0) {
    return(do.call(rbind, ks_results))
  } else {
    message("No valid Human1-Human column pairs found.")
    return(data.frame())
  }
}

ks_d1_results <- ks_Human1_vs_Human(dat_clean)
print(ks_d1_results)

ks_d1_results <- ks_d1_results %>%
  tidyr::separate(XY, into = c("relationship", "cooperative_function"), sep = "_", remove = FALSE)

# Add Benjamini-Hochberg p-values
ks_d1_results$fdr <- p.adjust(ks_d1_results$P_Value, method = "BH")
write.csv(
  ks_d1_results,
  file.path(SI_FIGURE_DIR, paste0("ks_results_", HUMAN_AI_ID, "_", STUDY_OUTPUT_ID, "_SI.csv")),
  row.names = FALSE
)

### Euclidean and K-S heatmaps ----
euclid_sum_d1 <- compute_euclid_summary(dat_long, "Human1", "Human")
euclid_aug_d1 <- build_euclid_augmented(ks_d1_results, euclid_sum_d1)

p_euclid_d1 <- plot_euclid_multicolor(
  euclid_aug_d1,
  "How should humans treat social AIs?",
  paste0("Kolmogorov-Smirnov and Euclidean distances in normative expectations\nSample: ",
         STUDY_LABEL)
)
p_euclid_d1
ggsave(file.path(KS_EUCLID_FIGURE_DIR, paste0("ks_euclid_heatmap_", HUMAN_AI_ID, "_", STUDY_OUTPUT_ID, ".pdf")),
       p_euclid_d1, width = 8, height = 6.4, units = "in")

p_euclid_d1_si <- plot_euclid_multicolor_si(
  euclid_aug_d1,
  "How should humans treat social AIs?",
  paste0("Kolmogorov-Smirnov and Euclidean distances in normative expectations\nSample: ",
         STUDY_LABEL)
)
p_euclid_d1_si
ggsave(file.path(SI_FIGURE_DIR, paste0("ks_euclid_heatmap_", HUMAN_AI_ID, "_", STUDY_OUTPUT_ID, "_SI.pdf")),
       p_euclid_d1_si, width = 8, height = 6.4, units = "in")


## Norm ridges plots ----
# Filter human-AI data
dat_long_d1 <- dat_long    %>% filter(role == HUMAN_AI_ID)
ci_d1        <- ci_data_all %>% filter(role == HUMAN_AI_ID)


# norm ridges plots
p_ridges_d1 <- ggplot(dat_long_d1,
    aes(x = value, y = fct_relevel(category, "Transaction", "Mating", "Hierarchy", "Care"))) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +

  # Human-Human (top half, lighter fill)
  geom_density_ridges(
    data   = filter(dat_long_d1, Condition == "Human condition"),
    aes(fill = category, group = interaction(category, relationship)),
    scale  = 0.5, alpha = 0.4, rel_min_height = 0.01, color = NA, show.legend = FALSE
  ) +

  # Human-AI (bottom half, darker fill)
  geom_density_ridges(
    data   = filter(dat_long_d1, Condition == "Superintelligent AI condition"),
    aes(fill = category, group = interaction(category, relationship)),
    scale  = -0.5, alpha = 0.9, rel_min_height = 0.01, color = NA, show.legend = FALSE
  ) +

  # median
  geom_point(
    data = ci_d1,
    aes(x = median_value, y = category, color = category),
    size = 3, inherit.aes = FALSE,
    position = position_nudge(y = ifelse(ci_d1$Condition == "Human condition", 0.1, -0.1))
  ) +

  # errorbars
  geom_errorbarh(
    data = ci_d1,
    aes(xmin = lower_ci, xmax = upper_ci, y = category, color = category),
    height = 0.15, size = 0.6, inherit.aes = FALSE,
    position = position_nudge(y = ifelse(ci_d1$Condition == "Human condition", 0.1, -0.1))
  ) +
  facet_wrap(~relationship, ncol = 7,
    labeller = as_labeller(c(
      "teacher"   = "Teacher",    "friends"   = "Friend",
      "romantic"  = "Romantic partner",       "seller"    = "Seller",
      "supervisor"= "Assistant",  "therapist" = "Mental health provider",
      "workers"   = "Coworker"
    ))
  ) +
  labs(x = "", y = "") +
  theme_classic(base_size = 18) +
  theme(
    plot.title       = element_text(size = 18, face = "bold"),
    axis.text.x      = element_text(size = 16),
    axis.text.y      = element_text(size = 16),
    strip.text       = element_text(size = 18, face = "bold"),
    strip.background = element_rect(color = "black", fill = NA, size = 0.8),
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.ticks.x     = element_blank(),
    legend.position  = "none"
  ) +
  scale_fill_manual(values  = fill_light) +
  scale_color_manual(values = c(
    "Care"        = "#C71585",
    "Hierarchy"   = "#696969",
    "Mating"      = "#228B22",
    "Transaction" = "#3277bd"
  )) +
  scale_x_continuous(
    breaks = c(-4, 0, 4),
    labels = c("Definitely\nshould not", "Neutral", "Definitely\nshould")
  )

save_norm_ridges_original_font110(
  data = dat_long_d1,
  original_main_plot = p_ridges_d1,
  direction_id = HUMAN_AI_ID,
  title_text = "How should humans treat social AIs?"
)

## Bayesian model using median value ----
library(brms)
library(dplyr)
library(purrr)
library(tidyr)
library(posterior)
library(coda)

glimpse(dat_long_d1)

# raw difference of median
dat_long_d1_diff <- dat_long_d1 %>%
  group_by(category, relationship, Condition) %>%
  summarise(raw_median = median(value), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = Condition, values_from = raw_median) %>%
  mutate(diff_AI_minus_Human = `Superintelligent AI condition` - `Human condition`) %>%
  arrange(relationship, category)

write.csv(dat_long_d1_diff,
          file.path(NORM_MEDIAN_DIFFERENCE_DIR, paste0("raw_median_difference_", HUMAN_AI_ID, "_", STUDY_OUTPUT_ID, ".csv")),
          row.names = FALSE)


# set factor level
dat_long_d1$Condition <- factor(
  dat_long_d1$Condition,
  levels = c("Human condition", "Superintelligent AI condition")
)
dat_long_d1$category <- factor(
  dat_long_d1$category,
  levels = c("Care", "Hierarchy", "Mating", "Transaction")
)
dat_long_d1$relationship <- factor(
  dat_long_d1$relationship,
  levels = c("friends", "romantic", "seller", "supervisor", "teacher", "therapist", "workers")
)

levels(dat_long_d1$Condition)
levels(dat_long_d1$category)
levels(dat_long_d1$relationship)


### Model ----
if (REFIT_MODELS) {
  library(cmdstanr)
  options(brms.backend = "cmdstanr")
  cmdstan_path()
  cmdstan_version()

  model_d1 <- brm(
    formula = bf(value ~ Condition * category * relationship + (1 | ResponseId),
                 quantile = 0.5),
    data    = dat_long_d1,
    family  = asym_laplace(),
    chains  = 4, iter = 4000, warmup = 2000, cores = 4,
    control = list(max_treedepth = 15),
    seed    = 123
  )
  summary(model_d1)

  saveRDS(model_d1, file = file.path(NORM_MODEL_DIR, paste0("model_median_brms_", HUMAN_AI_ID, "_", STUDY_MODEL_ID, ".rds")))
}

### Read model ----
model_d1 <- readRDS(file.path(NORM_MODEL_DIR, paste0("model_median_brms_", HUMAN_AI_ID, "_", STUDY_MODEL_ID, ".rds")))
summary(model_d1)



### Manual posterior contrasts ----

# Helper: weighted median
weighted_median <- function(x, w) {
  ord <- order(x)
  x <- x[ord]
  w <- w[ord]
  cw <- cumsum(w) / sum(w)
  x[which(cw >= 0.5)[1]]
}

# Cell counts / weights within each category
# These weights are used for weighted marginal medians
rel_weights_d1 <- dat_long_d1 %>%
  count(category, relationship) %>%
  group_by(category) %>%
  mutate(weight = n / sum(n)) %>%
  ungroup()

rel_weights_d1


# 1. Compare Superintelligent AI vs Human within each category x relationship
compare_one_cell <- function(model, data, cat_value, rel_value) {
  
  nd <- data.frame(
    Condition = factor(
      c("Superintelligent AI condition", "Human condition"),
      levels = levels(data$Condition)
    ),
    category = factor(
      c(cat_value, cat_value),
      levels = levels(data$category)
    ),
    relationship = factor(
      c(rel_value, rel_value),
      levels = levels(data$relationship)
    ),
    ResponseId = NA
  )
  
  lin <- posterior_linpred(
    model,
    newdata = nd,
    re_formula = NA
  )
  
  # Contrast: Superintelligent AI - Human
  diff_draws <- lin[, 1] - lin[, 2]

  # using rowMeans of median to get the marginal one

  hpd <- HPDinterval(as.mcmc(diff_draws), prob = 0.95)

  data.frame(
    category = cat_value,
    relationship = rel_value,
    contrast = "Superintelligent AI - Human",
    estimate = median(diff_draws),
    ci_lower = quantile(diff_draws, 0.025),
    ci_upper = quantile(diff_draws, 0.975),
    lower_HPD = hpd[1, "lower"],
    upper_HPD = hpd[1, "upper"],
    prob_gt_0 = mean(diff_draws > 0),
    prob_lt_0 = mean(diff_draws < 0)
  )
}


# run all cell-specific contrasts
res_cellwise_d1 <- expand.grid(
  category     = levels(dat_long_d1$category),
  relationship = levels(dat_long_d1$relationship),
  stringsAsFactors = FALSE
) %>%
  split(seq_len(nrow(.))) %>%
  map_dfr(~ compare_one_cell(model_d1, dat_long_d1, .x$category, .x$relationship))


# 2. Weighted marginal median contrast within each category
compare_one_category_weighted_median <- function(model, data, cat_value, weight_df) {
  sub_w      <- weight_df %>%
    filter(category == cat_value) %>%
    arrange(match(relationship, levels(data$relationship)))
  rel_levels <- as.character(sub_w$relationship)
  w          <- sub_w$weight

  make_nd <- function(cond) data.frame(
    Condition    = factor(rep(cond, length(rel_levels)), levels = levels(data$Condition)),
    category     = factor(rep(cat_value, length(rel_levels)), levels = levels(data$category)),
    relationship = factor(rel_levels, levels = levels(data$relationship)),
    ResponseId   = NA
  )

  lin_human <- posterior_linpred(model, newdata = make_nd("Human condition"),          re_formula = NA)
  lin_ai    <- posterior_linpred(model, newdata = make_nd("Superintelligent AI condition"), re_formula = NA)

  human_wmed <- apply(lin_human, 1, weighted_median, w = w)
  ai_wmed    <- apply(lin_ai,    1, weighted_median, w = w)
  diff_draws <- ai_wmed - human_wmed
  hpd        <- HPDinterval(as.mcmc(diff_draws), prob = 0.95)

  data.frame(
    category   = cat_value, contrast = "Superintelligent AI - Human",
    estimate   = median(diff_draws),
    ci_lower   = quantile(diff_draws, 0.025), ci_upper = quantile(diff_draws, 0.975),
    lower_HPD  = hpd[1, "lower"], upper_HPD = hpd[1, "upper"],
    prob_gt_0  = mean(diff_draws > 0), prob_lt_0 = mean(diff_draws < 0)
  )
}

# Run weighted marginal median contrasts for each category
res_categorywise_d1 <- map_dfr(
  levels(dat_long_d1$category),
  ~ compare_one_category_weighted_median(model_d1, dat_long_d1, .x, rel_weights_d1)
)

# 3. tidy results
round_results <- function(df) {
  df %>% mutate(across(c(estimate, ci_lower, ci_upper, lower_HPD, upper_HPD,
                         prob_gt_0, prob_lt_0), ~ round(.x, 4)))
}
res_cellwise_d1     <- round_results(res_cellwise_d1)
res_categorywise_d1 <- round_results(res_categorywise_d1)

res_cellwise_d1
res_categorywise_d1

write.csv(res_cellwise_d1,
  file.path(NORM_MANUAL_CONTRAST_DIR, paste0("res_cellwise_", HUMAN_AI_ID, "_", STUDY_OUTPUT_ID, ".csv")),
  row.names = FALSE)
write.csv(res_categorywise_d1,
  file.path(NORM_MANUAL_CONTRAST_DIR, paste0("res_categorywise_", HUMAN_AI_ID, "_", STUDY_OUTPUT_ID, ".csv")),
  row.names = FALSE)





## Human likeness ----
# This section constructs a human-likeness outcome by reusing the fitted norm model as a human normative reference. 
# Specifically, for each category × relationship cell, the posterior median under the Human condition is extracted, 
# treating it as the model-based human benchmark for that cell. these human posterior medians are then joined back to 
# the row-level dataset and define human_likeness as the absolute distance between each observed value and the 
# corresponding human posterior median. As a result, smaller values indicate responses that are closer to the human 
# normative benchmark, whereas larger values indicate responses that are farther away. In this sense, human_likeness 
# is not a direct subjective rating of "how human-like" a response is, but rather a model-based measure of deviation 
# from the human normative standard.
#
# human_likeness_signed is defined, defined as the signed difference
# between each observed AI value and the human posterior median. Positive values indicate AI responses above the
# human normative benchmark, whereas negative values indicate AI responses below the benchmark. This signed-difference
# analysis is intended as a supplementary analysis, while the absolute-difference model remains the primary analysis.
# 
# The final human-likeness model then uses human_likeness as the dependent variable and examines whether variation in this distance can be explained by 
# relationship, perceived realness, anxiety, other perception-related variables, LLM-use variables, and demographic covariates. 
# Compared with the earlier norm model, the focus here is no longer on comparing Human and AI across Condition × 
# category × relationship cells, but on explaining why some AI responses are closer to, or farther from, the human 
# normative benchmark. Before fitting the model, several stabilization steps are applied, including collapsing sparse 
# factor levels, standardizing continuous predictors, adding weakly informative priors, removing LLM_Past because of 
# its overlap with LLM_Freq, and using QR decomposition to reduce estimation difficulties caused by possible 
# multicollinearity among predictors. These steps are intended to improve model stability, identifiability, and 
# interpretability rather than to change the substantive question.

# get human posterior medians
get_cell_median <- function(model, data, cat_value, rel_value, condition_value) {
  nd  <- data.frame(
    Condition    = factor(condition_value, levels = levels(data$Condition)),
    category     = factor(cat_value,       levels = levels(data$category)),
    relationship = factor(rel_value,       levels = levels(data$relationship)),
    ResponseId   = NA
  )
  lin <- posterior_linpred(model, newdata = nd, re_formula = NA)
  data.frame(Condition = condition_value, category = cat_value, relationship = rel_value,
             posterior_median = median(lin[, 1]))
}


{
  # Extract Human condition posterior medians for each cell

  human_cell_medians_d1 <- expand.grid(
    category     = levels(dat_long_d1$category),
    relationship = levels(dat_long_d1$relationship),
    stringsAsFactors = FALSE
  ) %>%
    split(seq_len(nrow(.))) %>%
    map_dfr(~ get_cell_median(model_d1, dat_long_d1, .x$category, .x$relationship,
                              "Human condition"))

  glimpse(human_cell_medians_d1)

  # human_likeness_abs    = magnitude of deviation from the human normative benchmark
  # human_likeness_signed = signed deviation from the human normative benchmark
  human_median_join_d1 <- human_cell_medians_d1 %>%
    select(category, relationship, human_posterior_median = posterior_median)

  dat_v2_d1 <- dat_v2_d1 %>%
    left_join(human_median_join_d1, by = c("relationship", "category")) %>%
    mutate(
      human_likeness_signed = as.numeric(value) - human_posterior_median,
      human_likeness_abs    = abs(human_likeness_signed),
      human_likeness        = human_likeness_abs
    )

  dat_v2_d1$human_likeness
  glimpse(dat_v2_d1)

  # Match realness and anxiety to corresponding relationship
  dat_v2_d1 <- dat_v2_d1 %>%
    mutate(
      realness_matched = case_when(
        relationship == "friends"    ~ as.numeric(realness_friends),
        relationship == "romantic"   ~ as.numeric(realness_romantic),
        relationship == "workers"    ~ as.numeric(realness_workers),
        relationship == "supervisor" ~ as.numeric(realness_supervisor),
        relationship == "teacher"    ~ as.numeric(realness_teacher),
        relationship == "therapist"  ~ as.numeric(realness_therapist),
        relationship == "seller"     ~ as.numeric(realness_seller)
      ),
      anxiety_matched = case_when(
        relationship == "friends"    ~ as.numeric(anxiety_friends),
        relationship == "romantic"   ~ as.numeric(anxiety_romantic),
        relationship == "workers"    ~ as.numeric(anxiety_workers),
        relationship == "supervisor" ~ as.numeric(anxiety_supervisor),
        relationship == "teacher"    ~ as.numeric(anxiety_teacher),
        relationship == "therapist"  ~ as.numeric(anxiety_therapist),
        relationship == "seller"     ~ as.numeric(anxiety_seller)
      )
    )
}


### Model ----
dat_v2_d1 <- dat_v2_d1 %>%
  rename(
    subjective_experience = `subjective experience`,
    self_consciousness    = `self-consciousness`,
    well_being            = `well-being`,
    Political_Social      = `Political - Social_1`,
    Political_Economic    = `Political - Economic_1`
  )

options(brms.backend = "cmdstanr")

if (REFIT_MODELS) {

# Version 1: unadjusted model with all predictors
model_humanlikeness_d1_v1_rel <- brm(
  formula = bf(
    human_likeness ~ 
      relationship +
      realness_matched + anxiety_matched +
      subjective_experience + sentience + self_consciousness +
      well_being + autonomy + beneficence + cognition +
      LLM_Past + LLM_Freq + Gender + Age + Ethn +
      Relat.Status + Education + StatusLadder +
      Political_Social + Political_Economic + How.Religious_1 +
      (1 | ResponseId),
    quantile = 0.5
  ),
  data    = dat_v2_d1,
  family  = asym_laplace(),
  chains  = 4,
  iter    = 4000,
  warmup  = 2000,
  cores   = 4,
  control = list(max_treedepth = 10),
  seed    = 123
)

summary(model_humanlikeness_d1_v1_rel)
saveRDS(model_humanlikeness_d1_v1_rel,
        file = file.path(HUMANLIKENESS_MODEL_DIR, paste0("model_humanlikeness_", HUMAN_AI_ID, "_relationship_", STUDY_MODEL_ID, "_v1.rds")))


# Version 2: primary model for absolute-difference human likeness
# applying several stabilization steps, including collapsing sparse factor levels, 
# standardizing continuous predictors, adding weakly informative priors, removing LLM_Past because of 
# its overlap with LLM_Freq, and using QR decomposition to reduce estimation difficulties caused by possible 
# multicollinearity among predictors. These steps are intended to improve model stability, identifiability, and 
# interpretability rather than to change the substantive question.

# Note:
# min_n = 20 (lumping each option with fewer than 20 people selecting it)
# In the long format data, there are 12 lines per participant, so it was min_n: 20 * 12 = 240
lump_small <- function(x, min_n = 240, other_level = "Other_small_n") {
  forcats::fct_lump_min(factor(x), min = min_n, other_level = other_level)
}

dat_v2_d1_model <- dat_v2_d1 %>%
  mutate(
    ResponseId   = factor(ResponseId),
    
    Gender         = lump_small(Gender, min_n = 240),
    Ethn           = lump_small(Ethn, min_n = 240),
    `Relat.Status` = lump_small(`Relat.Status`, min_n = 240),
    Education      = lump_small(Education, min_n = 240),
    
    
    relationship = factor(
      relationship,
      levels = c("friends", "romantic", "seller", "supervisor", "teacher", "therapist", "workers")
    ),
    LLM_Freq = forcats::fct_relevel(
      factor(LLM_Freq),
      "Never"
    ),
    
    # Added: standardization of continuous predictors
    # Notes:
    # The use of prior(normal(0, 1), class = "b") as a weakly informative prior
    # generally assumes that predictors have been standardized
    # (mean = 0, SD = 1).
    # Without standardization, predictors may differ substantially in scale
    # (e.g., Age may range from 18 to 80, whereas Likert-type scales range from 1 to 7),
    # so the same normal(0, 1) prior would impose very different levels of constraint
    # across variables, undermining its interpretation as "weakly informative."
    # Here, scale() is used to perform z-score standardization, and as.numeric()
    # is applied to remove the matrix attribute returned by scale().
    
    Age              = as.numeric(scale(Age)),
    StatusLadder     = as.numeric(scale(StatusLadder)),
    How.Religious_1  = as.numeric(scale(How.Religious_1)),
    Political_Social   = as.numeric(scale(Political_Social)),
    Political_Economic = as.numeric(scale(Political_Economic)),
    
    subjective_experience = as.numeric(scale(subjective_experience)),
    sentience             = as.numeric(scale(sentience)),
    self_consciousness    = as.numeric(scale(self_consciousness)),
    well_being            = as.numeric(scale(well_being)),
    autonomy              = as.numeric(scale(autonomy)),
    beneficence           = as.numeric(scale(beneficence)),
    cognition             = as.numeric(scale(cognition)),
    
    realness_matched      = as.numeric(scale(realness_matched)),
    anxiety_matched       = as.numeric(scale(anxiety_matched))
  )

# Weakly informative priors are specified here instead of relying on the default flat priors.
# In brms, the default prior for population-level effects (class = "b") is an
# improper flat prior if no prior is explicitly set. Once prior(normal(0, 1), class = "b")
# is added, those coefficients are no longer assigned an improper flat prior, but a proper
# normal prior.
#
# Stan documentation generally recommends using informative or at least weakly informative
# priors in regression models, as they often improve estimation stability and sampling behavior.
#
# The priors used here provide conservative regularization:
# - b: normal(0, 1)       -> mildly shrinks regression coefficients toward zero,
#                            reducing weak identification and excessive drift
# - sd: exponential(1)    -> regularizes random-effect standard deviations on the positive scale
# - sigma: exponential(1) -> regularizes the residual scale on the positive scale
#
# If these priors later appear too restrictive, the prior for b can be relaxed
# to normal(0, 2).

priors_hl <- c(
  prior(normal(0, 1), class = "b"),
  prior(exponential(1), class = "sd"),
  prior(exponential(1), class = "sigma")
)


# Model

# Note 1: removing LLM_Past
# In the present results, LLM_PastYes and several LLM_Freq levels showed nearly mirrored
# large positive and negative coefficients, together with poor R-hat values. This pattern
# strongly suggests substantial information overlap between the two variables, leading to
# weak parameter identification. In this situation, it is more reasonable to retain the
# more fine-grained LLM_Freq variable and remove the coarser LLM_Past measure.
# 
# Note 2: decomp = "QR"
# The brms documentation states that decomp = "QR" can be helpful when fitting models
# with highly correlated predictors. Here, multiple perception-related variables
# are likely to be substantially correlated, so QR decomposition is enabled.
#
# Note 3: max_treedepth
# max_treedepth is not increased further here. According to the Stan diagnostics
# documentation, when max_treedepth warnings occur together with high R-hat values,
# the primary concern should be model geometry or parameter identification problems,
# rather than simply increasing max_treedepth.
# The main situation in which raising max_treedepth is justified is when R-hat is already
# acceptable (< 1.01) and ESS is sufficient, but "max treedepth exceeded" warnings still
# occur repeatedly. These warnings were not observed in the present model.

model_humanlikeness_d1_v2_rel <- brm(
  formula = bf(
    human_likeness ~
      relationship + 
      realness_matched + anxiety_matched +
      subjective_experience + sentience + self_consciousness +
      well_being + autonomy + beneficence + cognition +
      LLM_Freq + Gender + Age + Ethn +
      Relat.Status + Education + StatusLadder +
      Political_Social + Political_Economic + How.Religious_1 +
      (1 | ResponseId),
    quantile = 0.5,
    decomp   = "QR"
  ),
  data    = dat_v2_d1_model,
  family  = asym_laplace(),
  prior   = priors_hl,
  backend = "cmdstanr",
  chains  = 4,
  iter    = 4000,
  warmup  = 2000,
  cores   = 4,
  control = list(max_treedepth = 10),
  seed    = 123
)

summary(model_humanlikeness_d1_v2_rel)
saveRDS(
  model_humanlikeness_d1_v2_rel,
  file = file.path(HUMANLIKENESS_MODEL_DIR, paste0("model_humanlikeness_", HUMAN_AI_ID, "_relationship_", STUDY_MODEL_ID, "_v2.rds"))
)

# Use the relationship model as the primary model for the reporting code below.
model_humanlikeness_d1 <- model_humanlikeness_d1_v2_rel


# Supplementary signed-difference analysis
# This model asks whether the same predictors explain the direction of deviation from the
# human normative benchmark, rather than only the magnitude of deviation.
model_humanlikeness_signed_d1_v2_rel <- brm(
  formula = bf(
    human_likeness_signed ~
      relationship +
      realness_matched + anxiety_matched +
      subjective_experience + sentience + self_consciousness +
      well_being + autonomy + beneficence + cognition +
      LLM_Freq + Gender + Age + Ethn +
      Relat.Status + Education + StatusLadder +
      Political_Social + Political_Economic + How.Religious_1 +
      (1 | ResponseId),
    quantile = 0.5,
    decomp   = "QR"
  ),
  data    = dat_v2_d1_model,
  family  = asym_laplace(),
  prior   = priors_hl,
  backend = "cmdstanr",
  chains  = 4,
  iter    = 4000,
  warmup  = 2000,
  cores   = 4,
  control = list(max_treedepth = 10),
  seed    = 123
)

summary(model_humanlikeness_signed_d1_v2_rel)
saveRDS(
  model_humanlikeness_signed_d1_v2_rel,
  file = file.path(HUMANLIKENESS_MODEL_DIR, paste0("model_humanlikeness_signed_", HUMAN_AI_ID, "_relationship_", STUDY_MODEL_ID, "_v2_supplement.rds"))
)

model_humanlikeness_signed_d1 <- model_humanlikeness_signed_d1_v2_rel
} else {
  model_humanlikeness_d1 <- readRDS(
    file.path(HUMANLIKENESS_MODEL_DIR, paste0("model_humanlikeness_", HUMAN_AI_ID, "_relationship_", STUDY_MODEL_ID, "_v2.rds"))
  )
  model_humanlikeness_signed_d1 <- readRDS(
    file.path(HUMANLIKENESS_MODEL_DIR, paste0("model_humanlikeness_signed_", HUMAN_AI_ID, "_relationship_", STUDY_MODEL_ID, "_v2_supplement.rds"))
  )
}


summary(model_humanlikeness_d1)

library(dplyr)
library(brms)
library(posterior)
library(tibble)
library(purrr)
library(coda)

fixef(model_humanlikeness_d1)

fe_d1 <- as.data.frame(
  fixef(model_humanlikeness_d1,
        probs  = c(0.025, 0.975),
        robust = TRUE)
) %>%
  rownames_to_column("term") %>%
  rename(
    estimate  = Estimate,
    mad       = Est.Error, # median absolute deviation
    eti_lower = Q2.5,
    eti_upper = Q97.5
  )

fe_d1



draws_d1 <- as_draws_df(model_humanlikeness_d1)
b_cols_d1 <- grep("^b_", names(draws_d1), value = TRUE)

coef_extra_d1 <- map_dfr(b_cols_d1, function(v) {
  x <- draws_d1[[v]]

  hpd <- HPDinterval(as.mcmc(x), prob = 0.95)

  data.frame(
    term       = sub("^b_", "", v),
    hpdi_lower = unname(hpd[1, "lower"]),
    hpdi_upper = unname(hpd[1, "upper"]),
    prob_gt_0  = mean(x > 0),
    prob_lt_0  = mean(x < 0),
    pd         = max(mean(x > 0), mean(x < 0))
  )
})

coef_extra_d1


diag_d1 <- summarise_draws(
  draws_d1[, b_cols_d1],
  rhat,
  ess_bulk,
  ess_tail
) %>%
  mutate(term = sub("^b_", "", variable)) %>%
  select(term, rhat, ess_bulk, ess_tail)

diag_d1


report_table_d1 <- fe_d1 %>%
  left_join(coef_extra_d1, by = "term") %>%
  left_join(diag_d1,       by = "term") %>%
  mutate(across(where(is.numeric), ~ round(.x, 4)))

report_table_d1

write.csv(
  report_table_d1,
  file.path(HUMANLIKENESS_REPORT_DIR, paste0("humanlikeness_report_", HUMAN_AI_ID, "_relationship_", STUDY_OUTPUT_ID, "_v2.csv")),
  row.names = FALSE
)


# Supplementary report for the signed-difference model
fixef(model_humanlikeness_signed_d1)

fe_signed_d1 <- as.data.frame(
  fixef(model_humanlikeness_signed_d1,
        probs  = c(0.025, 0.975),
        robust = TRUE)
) %>%
  rownames_to_column("term") %>%
  rename(
    estimate  = Estimate,
    mad       = Est.Error,
    eti_lower = Q2.5,
    eti_upper = Q97.5
  )

draws_signed_d1 <- as_draws_df(model_humanlikeness_signed_d1)
b_cols_signed_d1 <- grep("^b_", names(draws_signed_d1), value = TRUE)

coef_extra_signed_d1 <- map_dfr(b_cols_signed_d1, function(v) {
  x <- draws_signed_d1[[v]]

  hpd <- HPDinterval(as.mcmc(x), prob = 0.95)

  data.frame(
    term       = sub("^b_", "", v),
    hpdi_lower = unname(hpd[1, "lower"]),
    hpdi_upper = unname(hpd[1, "upper"]),
    prob_gt_0  = mean(x > 0),
    prob_lt_0  = mean(x < 0),
    pd         = max(mean(x > 0), mean(x < 0))
  )
})

diag_signed_d1 <- summarise_draws(
  draws_signed_d1[, b_cols_signed_d1],
  rhat,
  ess_bulk,
  ess_tail
) %>%
  mutate(term = sub("^b_", "", variable)) %>%
  select(term, rhat, ess_bulk, ess_tail)

report_table_signed_d1 <- fe_signed_d1 %>%
  left_join(coef_extra_signed_d1, by = "term") %>%
  left_join(diag_signed_d1,       by = "term") %>%
  mutate(across(where(is.numeric), ~ round(.x, 4)))

report_table_signed_d1

write.csv(
  report_table_signed_d1,
  file.path(HUMANLIKENESS_REPORT_DIR, paste0("humanlikeness_signed_report_", HUMAN_AI_ID, "_relationship_", STUDY_OUTPUT_ID, "_v2_supplement.csv")),
  row.names = FALSE
)


cat("\nRun log saved to:", normalizePath(RUN_LOG_FILE, mustWork = FALSE), "\n")
invisible(dat_clean_for_distributions)
}

# Exploratory consensus/dispersion analysis ----
#
# Manuscript question:
# Which direction x relationship x relational-norm prescriptions are tightly
# concentrated (greater apparent consensus), and which are widely dispersed?
#
# Primary measure:
# Mean absolute distance from the cell median (MAD). Here, MAD explicitly
# denotes the mean—not median—absolute distance. It is expressed in the
# original -4 to +4 response-scale units; lower values indicate greater
# concentration around the typical response.
#
# Uncertainty and repeated observations:
# Bootstrap resampling occurs at the participant level within each study and
# condition. When a participant is sampled, all of that participant's repeated
# direction x relationship x norm ratings are retained together.
#
# Scope:
# These analyses are exploratory and descriptive. They quantify dispersion
# within cells; they do not test AI-condition versus matched human-condition
# differences in dispersion. Such a comparison would require a separate
# location-scale model and is not claimed in the current manuscript text.

run_exploratory_dispersion_analysis <- function(
  cleaned_studies,
  bootstrap_reps = 5000,
  seed = 8122026
) {
  out_dir <- here::here("outputs", "dispersion")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  # This module now maintains one publication-ready export only.
  output_filename <- "dispersion_all_cells_all_samples.csv"
  existing_outputs <- list.files(
    out_dir, full.names = TRUE, all.files = TRUE, no.. = TRUE
  )
  obsolete_outputs <- existing_outputs[
    basename(existing_outputs) != output_filename
  ]
  if (length(obsolete_outputs) > 0) {
    unlink(obsolete_outputs, recursive = TRUE)
  }

  relationship_levels <- c(
    "friends", "romantic", "seller", "supervisor",
    "teacher", "therapist", "workers"
  )
  category_levels <- c("Care", "Hierarchy", "Mating", "Transaction")
  role_levels <- c("ai-human", "human-ai")

  relationship_labels <- c(
    friends = "Friend",
    romantic = "Romantic partner",
    seller = "Seller / customer",
    supervisor = "Assistant / supervisor",
    teacher = "Teacher / student",
    therapist = "Mental health provider / patient",
    workers = "Coworker"
  )

  process_norm_data <- function(df, study_id) {
    study_settings <- STUDY_CONFIG[[study_id]]

    renamed <- df %>%
      rename_with(~ str_replace_all(.x, c(
        "^X1_" = "romantic_",
        "^X2_" = "friends_",
        "^X3_" = "workers_",
        "^X4_" = "supervisor_",
        "^X5_" = "teacher_",
        "^X6_" = "therapist_",
        "^X8_" = "seller_"
      )))

    norm_pattern <- paste0(
      "^(romantic|friends|workers|supervisor|teacher|therapist|seller)",
      "_(Care|Hierarchy|Mating|Transaction)",
      "_(Human1|Human2|AI|Human)$"
    )
    norm_columns <- names(renamed)[str_detect(names(renamed), norm_pattern)]

    if (length(norm_columns) != 112) {
      stop(
        "Expected 112 relational-norm columns after renaming for ",
        study_id, "; found ", length(norm_columns), "."
      )
    }

    long <- renamed %>%
      select(ResponseId, all_of(norm_columns)) %>%
      pivot_longer(
        cols = all_of(norm_columns),
        names_to = "variable",
        values_to = "value_raw"
      ) %>%
      filter(!is.na(value_raw), as.character(value_raw) != "") %>%
      separate(
        variable,
        into = c("relationship", "category", "direction"),
        sep = "_",
        remove = TRUE
      ) %>%
      mutate(
        # parse_number handles both ordinary numeric values and the Qualtrics
        # endpoint labels such as "Definitely should not\n-4\n".
        value = suppressWarnings(readr::parse_number(as.character(value_raw))),
        Condition = case_when(
          direction %in% c("Human1", "Human2") ~ "Human condition",
          direction %in% c("AI", "Human") ~ "Superintelligent AI condition",
          TRUE ~ NA_character_
        ),
        role = case_when(
          direction %in% c("Human2", "AI") ~ "ai-human",
          direction %in% c("Human1", "Human") ~ "human-ai",
          TRUE ~ NA_character_
        ),
        study_id = study_id,
        study = study_settings$display_label,
        relationship = factor(relationship, levels = relationship_levels),
        category = factor(category, levels = category_levels),
        role = factor(role, levels = role_levels),
        cell_id = paste(role, relationship, category, sep = "___")
      ) %>%
      select(
        study_id, study, ResponseId, Condition, role,
        relationship, category, cell_id, value
      )

    invalid_values <- long %>%
      filter(is.na(value) | value < -4 | value > 4)
    if (nrow(invalid_values) > 0) {
      stop(
        "Found ", nrow(invalid_values),
        " invalid relational-norm values in ", study_id, "."
      )
    }

    duplicate_cells <- long %>%
      count(ResponseId, Condition, role, relationship, category) %>%
      filter(n > 1)
    if (nrow(duplicate_cells) > 0) {
      stop(
        "A participant contributed more than one observation to a cell in ",
        study_id, "."
      )
    }

    long
  }

  summarize_vector <- function(x) {
    x <- x[!is.na(x)]
    cell_median <- stats::median(x)

    c(
      n = length(x),
      median = cell_median,
      q1 = stats::quantile(x, 0.25, names = FALSE, type = 7),
      q3 = stats::quantile(x, 0.75, names = FALSE, type = 7),
      IQR = stats::IQR(x, type = 7),
      MAD = mean(abs(x - cell_median)),
      proportion_within_one = mean(abs(x - cell_median) <= 1)
    )
  }

  bootstrap_one_study_condition <- function(group_data, reps, group_seed) {
    group_identity <- group_data %>%
      distinct(study_id, study, Condition)
    if (nrow(group_identity) != 1) {
      stop("Bootstrap input must contain exactly one study x condition group.")
    }

    cell_lookup <- group_data %>%
      distinct(cell_id, role, relationship, category)

    wide <- group_data %>%
      select(ResponseId, cell_id, value) %>%
      tidyr::pivot_wider(names_from = cell_id, values_from = value)

    response_matrix <- as.matrix(wide %>% select(-ResponseId))
    storage.mode(response_matrix) <- "double"

    cell_lookup <- cell_lookup[
      match(colnames(response_matrix), cell_lookup$cell_id),
      ,
      drop = FALSE
    ]

    expected_cells <- length(role_levels) *
      length(relationship_levels) * length(category_levels)
    if (ncol(response_matrix) != expected_cells) {
      stop(
        "Expected ", expected_cells, " cells for ",
        group_identity$study_id, " / ", group_identity$Condition,
        "; found ", ncol(response_matrix), "."
      )
    }

    observed <- t(vapply(
      seq_len(ncol(response_matrix)),
      function(j) summarize_vector(response_matrix[, j]),
      FUN.VALUE = c(
        n = 0, median = 0, q1 = 0, q3 = 0,
        IQR = 0, MAD = 0, proportion_within_one = 0
      )
    ))

    boot_mad <- matrix(NA_real_, nrow = reps, ncol = ncol(response_matrix))
    boot_iqr <- matrix(NA_real_, nrow = reps, ncol = ncol(response_matrix))
    boot_within_one <- matrix(
      NA_real_, nrow = reps, ncol = ncol(response_matrix)
    )

    set.seed(group_seed)
    for (b in seq_len(reps)) {
      # Sampling matrix rows resamples participant IDs. Every repeated rating
      # from a sampled participant is therefore retained in the same draw.
      sampled_rows <- sample.int(
        nrow(response_matrix),
        size = nrow(response_matrix),
        replace = TRUE
      )
      sampled_matrix <- response_matrix[sampled_rows, , drop = FALSE]

      boot_stats <- t(vapply(
        seq_len(ncol(sampled_matrix)),
        function(j) summarize_vector(sampled_matrix[, j]),
        FUN.VALUE = c(
          n = 0, median = 0, q1 = 0, q3 = 0,
          IQR = 0, MAD = 0, proportion_within_one = 0
        )
      ))

      boot_mad[b, ] <- boot_stats[, "MAD"]
      boot_iqr[b, ] <- boot_stats[, "IQR"]
      boot_within_one[b, ] <- boot_stats[, "proportion_within_one"]
    }

    percentile_interval <- function(draw_matrix) {
      t(apply(
        draw_matrix,
        2,
        stats::quantile,
        probs = c(0.025, 0.975),
        na.rm = TRUE,
        names = FALSE,
        type = 7
      ))
    }

    mad_ci <- percentile_interval(boot_mad)
    iqr_ci <- percentile_interval(boot_iqr)
    within_one_ci <- percentile_interval(boot_within_one)

    result <- bind_cols(
      group_identity[rep(1, nrow(cell_lookup)), ],
      cell_lookup,
      as_tibble(observed)
    ) %>%
      mutate(
        MAD_lower = mad_ci[, 1],
        MAD_upper = mad_ci[, 2],
        IQR_lower = iqr_ci[, 1],
        IQR_upper = iqr_ci[, 2],
        proportion_within_one_lower = within_one_ci[, 1],
        proportion_within_one_upper = within_one_ci[, 2],
        bootstrap_reps = reps
      )

    result
  }

  norm_long <- purrr::imap_dfr(
    cleaned_studies,
    ~ process_norm_data(.x, .y)
  )

  # Audit how many repeated norm ratings each participant contributed. In the
  # present design, participants normally contribute 12 ratings per direction
  # (three randomly selected relationships x four norm categories).
  repeated_observation_audit <- norm_long %>%
    count(study_id, study, Condition, ResponseId, role, name = "n_ratings") %>%
    count(
      study_id, study, Condition, role, n_ratings,
      name = "n_participants"
    ) %>%
    arrange(study_id, Condition, role, n_ratings)

  bootstrap_groups <- split(
    norm_long,
    interaction(norm_long$study_id, norm_long$Condition, drop = TRUE)
  )

  message(
    "Running participant-level cluster bootstrap for dispersion (",
    bootstrap_reps, " resamples per study x condition group)..."
  )

  dispersion_summary <- purrr::imap_dfr(
    bootstrap_groups,
    ~ bootstrap_one_study_condition(
      .x,
      reps = bootstrap_reps,
      group_seed = seed + match(.y, names(bootstrap_groups))
    )
  ) %>%
    group_by(study_id, Condition, role) %>%
    mutate(
      # Rank 1 is the most concentrated cell within a direction and condition.
      MAD_rank_low_to_high = rank(MAD, ties.method = "average"),
      # Each direction x condition contains 28 cells, so this maps the rank to
      # 0 (most concentrated) through 1 (most dispersed).
      MAD_percentile = (MAD_rank_low_to_high - 1) / (dplyr::n() - 1),
      apparent_consensus_rank_high_to_low = rank(
        MAD,
        ties.method = "average"
      )
    ) %>%
    ungroup() %>%
    mutate(
      relationship_label = unname(
        relationship_labels[as.character(relationship)]
      ),
      direction_label = case_when(
        role == "ai-human" ~ "How social AIs should treat humans",
        role == "human-ai" ~ "How humans should treat social AIs"
      )
    ) %>%
    arrange(study_id, Condition, role, MAD_rank_low_to_high)

  all_cells_path <- file.path(
    out_dir,
    output_filename
  )

  # Format the only exported file as the Supporting Information table:
  # one row per question x relationship x norm x comparison role, and one
  # column per sample. Each result cell contains MAD and its 95% CI only.
  sample_metadata <- norm_long %>%
    distinct(study_id, study, ResponseId) %>%
    count(study_id, study, name = "sample_n") %>%
    mutate(
      study_order = match(study_id, STUDY_IDS),
      sample_column = paste0(study, " (N = ", sample_n, ")")
    ) %>%
    arrange(study_order)

  sample_columns <- sample_metadata$sample_column

  dispersion_export <- dispersion_summary %>%
    left_join(
      sample_metadata %>% select(study_id, study, sample_column),
      by = c("study_id", "study")
    ) %>%
    mutate(
      panel_order = match(as.character(role), role_levels),
      relationship_order = match(
        as.character(relationship), relationship_levels
      ),
      norm_order = match(as.character(category), category_levels),
      comparison_role_order = case_when(
        Condition == "Human condition" ~ 1L,
        Condition == "Superintelligent AI condition" ~ 2L
      ),
      Question = case_when(
        role == "ai-human" ~
          "(a) How should social AIs treat humans?",
        role == "human-ai" ~
          "(b) How should humans treat social AIs?"
      ),
      Relationship = relationship_label,
      Norm = as.character(category),
      Role = case_when(
        Condition == "Human condition" ~ "Human in role",
        Condition == "Superintelligent AI condition" ~ "AI in role"
      ),
      `MAD [95% CI]` = sprintf(
        "%.2f [%.2f, %.2f]", MAD, MAD_lower, MAD_upper
      )
    ) %>%
    select(
      panel_order, relationship_order, norm_order, comparison_role_order,
      Question, Relationship, Norm, Role, sample_column, `MAD [95% CI]`
    ) %>%
    tidyr::pivot_wider(
      names_from = sample_column,
      values_from = `MAD [95% CI]`
    ) %>%
    arrange(
      panel_order, relationship_order, norm_order, comparison_role_order
    ) %>%
    select(Question, Relationship, Norm, Role, all_of(sample_columns))

  readr::write_csv(dispersion_export, all_cells_path)

  ai_condition_summary <- dispersion_summary %>%
    filter(Condition == "Superintelligent AI condition")

  # Internal check of cross-sample consistency in the relative ordering of the
  # 28 AI-condition cells in each direction. These descriptive Spearman rank
  # correlations are retained in the returned R object but are not exported or
  # reported in the current manuscript or Supporting Information.
  study_ids <- names(cleaned_studies)
  study_pairs <- combn(study_ids, 2, simplify = FALSE)
  correlation_rows <- list()
  row_index <- 1

  for (role_value in role_levels) {
    for (study_pair in study_pairs) {
      first_study <- ai_condition_summary %>%
        filter(study_id == study_pair[1], role == role_value) %>%
        select(relationship, category, MAD_first = MAD)
      second_study <- ai_condition_summary %>%
        filter(study_id == study_pair[2], role == role_value) %>%
        select(relationship, category, MAD_second = MAD)

      matched_cells <- inner_join(
        first_study,
        second_study,
        by = c("relationship", "category")
      )

      correlation_rows[[row_index]] <- tibble(
        direction = role_value,
        direction_label = if_else(
          role_value == "ai-human",
          "How social AIs should treat humans",
          "How humans should treat social AIs"
        ),
        study_1 = study_pair[1],
        study_2 = study_pair[2],
        n_cells = nrow(matched_cells),
        spearman_MAD = stats::cor(
          matched_cells$MAD_first,
          matched_cells$MAD_second,
          method = "spearman",
          use = "complete.obs"
        )
      )
      row_index <- row_index + 1
    }
  }

  cross_sample_consistency <- bind_rows(correlation_rows)
  # Internal audit of cells discussed in the current manuscript. This object
  # makes it easy to check whether their qualitative patterns are similar in
  # all three samples. It is returned to the R session but is not exported.
  manuscript_cases <- tribble(
    ~case_id, ~role, ~relationship, ~category, ~manuscript_description,
    1, "ai-human", "romantic", "Mating",
    "Romantic or sexual behavior by AI romantic partners",
    2, "ai-human", "teacher", "Hierarchy",
    "Hierarchical behavior by AI teachers toward human students",
    3, "ai-human", "seller", "Mating",
    "Romantic or sexual behavior by AI sellers",
    4, "ai-human", "supervisor", "Mating",
    "Romantic or sexual behavior by AI assistants",
    5, "ai-human", "teacher", "Mating",
    "Romantic or sexual behavior by AI teachers",
    6, "ai-human", "therapist", "Mating",
    "Romantic or sexual behavior by AI mental health providers",
    7, "ai-human", "friends", "Hierarchy",
    "Hierarchical behavior by AI friends (provisional text; not an extreme)",
    8, "ai-human", "therapist", "Care",
    "Caring behavior by AI mental health providers (recommended replacement)",
    9, "human-ai", "teacher", "Hierarchy",
    "Authority exercised by human students toward AI teachers",
    10, "human-ai", "friends", "Hierarchy",
    "Authority exercised by humans toward AI friends",
    11, "human-ai", "romantic", "Mating",
    "Romantic or sexual behavior by humans toward AI romantic partners",
    12, "human-ai", "seller", "Mating",
    "Romantic or sexual behavior by humans toward AI sellers",
    13, "human-ai", "supervisor", "Mating",
    "Romantic or sexual behavior by humans toward AI assistants",
    14, "human-ai", "teacher", "Mating",
    "Romantic or sexual behavior by humans toward AI teachers",
    15, "human-ai", "therapist", "Mating",
    "Romantic or sexual behavior by humans toward AI mental health providers",
    16, "human-ai", "supervisor", "Hierarchy",
    "Authority exercised by human supervisors toward AI assistants"
  ) %>%
    mutate(
      role = factor(role, levels = role_levels),
      relationship = factor(relationship, levels = relationship_levels),
      category = factor(category, levels = category_levels)
    )

  manuscript_case_results <- ai_condition_summary %>%
    inner_join(
      manuscript_cases,
      by = c("role", "relationship", "category")
    ) %>%
    arrange(case_id, factor(study_id, levels = c("Pilot", "UK", "US"))) %>%
    select(
      case_id, manuscript_description, study_id, study,
      direction_label, relationship, category, n, median,
      q1, q3, IQR, IQR_lower, IQR_upper,
      MAD, MAD_lower, MAD_upper,
      proportion_within_one,
      proportion_within_one_lower,
      proportion_within_one_upper,
      MAD_rank_low_to_high, MAD_percentile
    )

  message("Saved exploratory dispersion outputs to: ", out_dir)
  print(cross_sample_consistency)

  invisible(list(
    all_cells = dispersion_export,
    ai_condition_cells = ai_condition_summary,
    cross_sample_consistency = cross_sample_consistency,
    manuscript_cases = manuscript_case_results,
    output_files = all_cells_path
  ))
}

# Supplementary distribution figures ----
run_distribution_analysis <- function(cleaned_studies) {
  out_dir <- here::here("outputs", "figures", "SI")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  studies <- tibble(
    study_id = names(cleaned_studies),
    study = map_chr(names(cleaned_studies), ~ STUDY_CONFIG[[.x]]$display_label)
  )

  relationship_labels <- c(
    friends = "Friend",
    romantic = "Romantic partner",
    seller = "Seller",
    supervisor = "Assistant",
    teacher = "Teacher",
    therapist = "Mental health provider",
    workers = "Coworker"
  )

  capacity_labels <- c(
    subjective_experience = "Subjective experience",
    sentience = "Sentience",
    self_consciousness = "Self-consciousness",
    well_being = "Well-being",
    autonomy = "Autonomy",
    beneficence = "Beneficence",
    cognition = "Cognition"
  )

  relationship_colors <- c(
    romantic = "#E76F51",
    friends = "#4A90D9",
    workers = "#7B8794",
    supervisor = "#E9C46A",
    teacher = "#2A9D8F",
    therapist = "#CC79A7",
    seller = "#009E73"
  )

  capacity_colors <- c(
    subjective_experience = "#264653",
    sentience = "#4A90D9",
    self_consciousness = "#E76F51",
    well_being = "#2A9D8F",
    autonomy = "#E9C46A",
    beneficence = "#CC79A7",
    cognition = "#7B8794"
  )

  theme_distribution <- function(base_size = 11) {
    theme_classic(base_size = base_size) +
      theme(
        text = element_text(family = "Arial"),
        plot.title = element_text(size = base_size + 3, face = "plain", hjust = 0.5),
        plot.margin = margin(8, 8, 20, 8),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = base_size - 1, color = "#4F5661", margin = margin(r = 8)),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = base_size - 4, color = "#4F5661"),
        axis.line.y = element_line(color = "#4F5661", linewidth = 0.3),
        axis.ticks.y = element_line(color = "#4F5661", linewidth = 0.3),
        strip.text.x = element_text(size = base_size - 3, face = "bold"),
        strip.text.y = element_text(size = base_size - 1, face = "plain", angle = 0),
        strip.background = element_rect(fill = "white", color = "#222222", linewidth = 0.7),
        panel.spacing.x = grid::unit(1.6, "lines"),
        panel.spacing.y = grid::unit(0.5, "lines"),
        legend.position = "none"
      )
  }

  process_distribution_data <- function(df) {
    df %>%
      rename_with(~ str_replace_all(., c(
        "AI_Anxiety_1" = "anxiety_romantic",
        "AI_Anxiety_2" = "anxiety_friends",
        "AI_Anxiety_3" = "anxiety_workers",
        "AI_Anxiety_4" = "anxiety_supervisor",
        "AI_Anxiety_5" = "anxiety_teacher",
        "AI_Anxiety_6" = "anxiety_therapist",
        "AI_Anxiety_7" = "anxiety_seller"
      ))) %>%
      rename_with(~ str_replace_all(., c(
        "Believability_1$" = "realness_romantic",
        "Believability_2$" = "realness_friends",
        "Believability_3$" = "realness_workers",
        "Believability_4$" = "realness_supervisor",
        "Believability_5$" = "realness_teacher",
        "Believability_6$" = "realness_therapist",
        "Believability_7$" = "realness_seller"
      ))) %>%
      rename_with(~ str_replace_all(., c(
        "MindQualities_1" = "subjective_experience",
        "MindQualities_2$" = "sentience",
        "MindQualities_3$" = "self_consciousness",
        "MindQualities_4$" = "well_being",
        "MindQualities_5$" = "autonomy",
        "MindQualities_6$" = "beneficence",
        "MindQualities_7$" = "cognition"
      ))) %>%
      mutate(across(
        all_of(names(capacity_labels)),
        ~ case_when(
          . == "Strongly disagree" ~ "1",
          . == "Somewhat disagree" ~ "2",
          . == "Neither agree nor disagree" ~ "3",
          . == "Somewhat agree" ~ "4",
          . == "Strongly agree" ~ "5",
          TRUE ~ as.character(.)
        )
      )) %>%
      mutate(across(
        starts_with("anxiety_"),
        ~ case_when(
          . == "Not at all afraid\n1\n" ~ "1",
          . == "Extremely afraid\n10\n" ~ "10",
          TRUE ~ as.character(.)
        )
      )) %>%
      mutate(across(
        starts_with("realness_"),
        ~ case_when(
          . == "Not at all real\n1" ~ "1",
          . == "Completely real\n10" ~ "10",
          TRUE ~ as.character(.)
        )
      )) %>%
      mutate(across(
        c(all_of(names(capacity_labels)), starts_with("anxiety_"), starts_with("realness_")),
        ~ suppressWarnings(as.numeric(.))
      ))
  }

  make_capacity_long <- function(df) {
    df %>%
      select(study, ResponseId, all_of(names(capacity_labels))) %>%
      pivot_longer(
        cols = all_of(names(capacity_labels)),
        names_to = "item",
        values_to = "value"
      ) %>%
      filter(!is.na(value)) %>%
      mutate(
        study = factor(study, levels = studies$study),
        item = factor(item, levels = names(capacity_labels), labels = capacity_labels)
      )
  }

  make_relationship_long <- function(df, prefix) {
    cols <- paste0(prefix, "_", names(relationship_labels))

    df %>%
      select(study, ResponseId, all_of(cols)) %>%
      pivot_longer(
        cols = all_of(cols),
        names_to = "item",
        values_to = "value"
      ) %>%
      filter(!is.na(value)) %>%
      mutate(
        item = str_remove(item, paste0("^", prefix, "_")),
        study = factor(study, levels = studies$study),
        item = factor(item, levels = names(relationship_labels), labels = relationship_labels)
      )
  }

  add_distribution_layers <- function(data, color_values) {
    medians <- data %>%
      group_by(study, item) %>%
      summarize(median_value = median(value, na.rm = TRUE), .groups = "drop")

    list(
      geom_density(
        aes(fill = item),
        adjust = 1.05,
        alpha = 0.68,
        color = NA,
        linewidth = 0,
        na.rm = TRUE
      ),
      geom_vline(
        data = medians,
        aes(xintercept = median_value, color = item),
        linetype = "dashed",
        linewidth = 0.55,
        alpha = 0.9,
        show.legend = FALSE
      ),
      scale_fill_manual(values = color_values),
      scale_color_manual(values = color_values)
    )
  }

  endpoint_label_data <- function(data, x_min, x_max, left_label, right_label) {
    bottom_study <- tail(levels(data$study), 1)

    tidyr::expand_grid(
      study = bottom_study,
      item = levels(data$item),
      endpoint = c("left", "right")
    ) %>%
      mutate(
        study = factor(study, levels = levels(data$study)),
        item = factor(item, levels = levels(data$item)),
        x = if_else(endpoint == "left", x_min, x_max),
        label = if_else(endpoint == "left", left_label, right_label),
        hjust = if_else(endpoint == "left", 0, 1)
      )
  }

  add_endpoint_labels <- function(data, x_min, x_max, left_label, right_label) {
    geom_text(
      data = endpoint_label_data(data, x_min, x_max, left_label, right_label),
      aes(x = x, y = -Inf, label = label, hjust = hjust),
      inherit.aes = FALSE,
      vjust = 1.15,
      size = 2.3,
      family = "Arial",
      color = "#4F5661",
      lineheight = 0.9
    )
  }

  audit_variable_ranges <- function(df) {
    audit_one <- function(vars, measure, expected_min, expected_max) {
      df %>%
        select(study, all_of(vars)) %>%
        pivot_longer(-study, names_to = "variable", values_to = "value") %>%
        group_by(study, variable) %>%
        summarize(
          measure = measure,
          n_nonmissing = sum(!is.na(value)),
          min = min(value, na.rm = TRUE),
          max = max(value, na.rm = TRUE),
          invalid = sum(!is.na(value) & (value < expected_min | value > expected_max)),
          .groups = "drop"
        ) %>%
        select(measure, study, variable, n_nonmissing, min, max, invalid)
    }

    bind_rows(
      audit_one(names(capacity_labels), "AI capacity evaluations", 1, 5),
      audit_one(paste0("realness_", names(relationship_labels)), "Relationship realness", 1, 10),
      audit_one(paste0("anxiety_", names(relationship_labels)), "Fearfulness of AI", 1, 10)
    ) %>%
      mutate(study = factor(study, levels = studies$study)) %>%
      arrange(measure, study, variable)
  }

  plot_distribution <- function(data, title, colors, x_min, x_max, left_label, right_label) {
    ggplot(data, aes(x = value)) +
      add_distribution_layers(data, unname(colors)) +
      facet_grid(study ~ item) +
      scale_x_continuous(
        breaks = c(x_min, x_max),
        limits = c(x_min, x_max),
        expand = expansion(mult = c(0, 0))
      ) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
      add_endpoint_labels(data, x_min, x_max, left_label, right_label) +
      coord_cartesian(clip = "off") +
      labs(title = title, x = NULL, y = NULL) +
      theme_distribution()
  }

  save_distribution_plot <- function(plot, stem, width = 13.2, height = 5.8) {
    pdf_path <- file.path(out_dir, paste0(stem, ".pdf"))
    ggsave(pdf_path, plot, width = width, height = height, units = "in", device = cairo_pdf)
    message("Saved: ", pdf_path)
    invisible(pdf_path)
  }

  all_data <- imap_dfr(cleaned_studies, function(df, study_id) {
    mutate(df, study = STUDY_CONFIG[[study_id]]$display_label)
  }) %>%
    process_distribution_data()

  message("Final analyzed sample sizes for distribution figures:")
  distribution_sample_sizes <- all_data %>%
    dplyr::mutate(study = factor(study, levels = studies$study)) %>%
    dplyr::count(study, name = "n")
  print(as.data.frame(distribution_sample_sizes), row.names = FALSE)

  capacity_long <- make_capacity_long(all_data)
  realness_long <- make_relationship_long(all_data, "realness")
  fearfulness_long <- make_relationship_long(all_data, "anxiety")

  audit_path <- file.path(out_dir, "si_distribution_variable_range_audit.csv")
  readr::write_csv(audit_variable_ranges(all_data), audit_path)
  message("Saved: ", audit_path)

  capacity_path <- save_distribution_plot(
    plot_distribution(
      capacity_long,
      "AI capacity evaluations",
      capacity_colors,
      1, 5,
      "Strongly\ndisagree",
      "Strongly\nagree"
    ),
    "fig_si_ai_capacity_distributions"
  )

  realness_path <- save_distribution_plot(
    plot_distribution(
      realness_long,
      "Relationship realness",
      relationship_colors,
      1, 10,
      "Not at all\nreal",
      "Completely\nreal"
    ),
    "fig_si_relationship_realness_distributions"
  )

  fearfulness_path <- save_distribution_plot(
    plot_distribution(
      fearfulness_long,
      "Fearfulness of AI",
      relationship_colors,
      1, 10,
      "Not at all\nafraid",
      "Extremely\nafraid"
    ),
    "fig_si_ai_fearfulness_distributions"
  )

  invisible(c(audit_path, capacity_path, realness_path, fearfulness_path))
}

# Simulation-based power analysis ----
run_simulation_power_analysis <- function() {
  suppressPackageStartupMessages(library(pbapply))

  out_dir <- here::here("outputs", "power")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  set.seed(123)

  delta <- 0.50
  sd_subject <- 0.89
  sd_resid <- 2.16
  alpha <- 0.05
  nsims <- 1000
  expected_exclusion_rate <- 0.16
  n_grid <- seq(100, 1800, by = 100)

  simulate_cellwise_power <- function(n_total) {
    message("Simulating power for N = ", n_total, "...")

    pvals <- pbapply::pbsapply(seq_len(nsims), function(i) {
      n_hai <- n_total %/% 2
      n_hhh <- n_total - n_hai
      conditions <- c(rep("Human-AI", n_hai), rep("Human-human", n_hhh))
      subject_b <- rnorm(n_total, mean = 0, sd = sd_subject)
      mu <- ifelse(conditions == "Human-AI", delta, 0)
      value <- rnorm(n_total, mean = mu + subject_b, sd = sd_resid)

      sim_data <- data.frame(
        Subject = factor(seq_len(n_total)),
        Condition = factor(conditions, levels = c("Human-human", "Human-AI")),
        value = value
      )

      fit <- tryCatch(lm(value ~ Condition, data = sim_data), error = function(e) NULL)
      if (is.null(fit)) {
        return(NA_real_)
      }

      as.numeric(summary(fit)$coefficients["ConditionHuman-AI", "Pr(>|t|)"])
    })

    tibble(
      N_total = n_total,
      n_per_condition = n_total / 2,
      Power = mean(pvals < alpha, na.rm = TRUE)
    )
  }

  power_results <- map_dfr(n_grid, simulate_cellwise_power)
  min_n_90 <- power_results %>% filter(Power >= 0.90) %>% slice(1)

  if (nrow(min_n_90) == 0) {
    stop("No tested sample size reached 90% power. Extend n_grid and rerun.")
  }

  recruit_n <- ceiling(min_n_90$N_total / (1 - expected_exclusion_rate))

  results_file <- file.path(out_dir, "simulation_power_results.csv")
  plot_file <- file.path(out_dir, "simulation_power_curve.pdf")
  summary_file <- file.path(out_dir, "simulation_power_summary.txt")
  writeup_file <- file.path(out_dir, "supplement_power_analysis_writeup.md")

  readr::write_csv(power_results, results_file)

  power_plot <- ggplot(power_results, aes(x = N_total, y = Power)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    geom_hline(yintercept = 0.90, linetype = "dashed", color = "red") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
    labs(
      title = "Power to detect a raw condition difference of 0.50",
      x = "Final analyzable total N",
      y = "Power"
    ) +
    theme_minimal(base_size = 14)

  ggsave(plot_file, power_plot, width = 7, height = 5, units = "in")

  summary_lines <- c(
    "Simulation-based power analysis",
    "",
    "Random seed: 123",
    paste0("Number of simulations per tested N: ", nsims),
    paste0("Alpha: ", alpha),
    paste0("Target raw condition difference, delta: ", delta),
    paste0("USA Pilot participant-level SD: ", sd_subject),
    paste0("USA Pilot residual SD: ", sd_resid),
    "",
    paste0(
      "Minimum tested final analyzable N reaching >= 90% power: ",
      min_n_90$N_total,
      " (power = ", sprintf("%.1f%%", 100 * min_n_90$Power), ")"
    ),
    paste0(
      "With an expected exclusion rate of ",
      sprintf("%.0f%%", 100 * expected_exclusion_rate),
      ", recruit at least ", recruit_n,
      " participants to retain ", min_n_90$N_total,
      " final analyzable cases."
    ),
    "",
    paste0("Power curve CSV: ", results_file),
    paste0("Power curve PDF: ", plot_file)
  )

  writeLines(summary_lines, summary_file)
  writeLines(summary_lines)

  writeup_lines <- c(
    "# Simulation-Based Power Analysis Write-Up",
    "",
    "## Suggested Supplement Text",
    "",
    paste0(
      "To determine the planned sample size for the two full studies, we conducted a simulation-based ",
      "power analysis using variance components estimated from the USA Pilot study. The simulation ",
      "targeted the simple Human-AI versus Human-human condition difference within a single prespecified ",
      "Relationship x Category x Direction cell. Participant-level heterogeneity had SD = ", sd_subject,
      ", residual error had SD = ", sd_resid,
      ", and the target raw mean difference was Delta = ", sprintf("%.2f", delta), "."
    ),
    "",
    paste0(
      "Using ", format(nsims, big.mark = ","),
      " simulations per tested sample size, the smallest tested final analyzable sample reaching ",
      "approximately 90% power was N = ", min_n_90$N_total,
      " (", min_n_90$n_per_condition, " per condition; estimated power = ",
      sprintf("%.1f%%", 100 * min_n_90$Power), "). With an expected exclusion rate of ",
      sprintf("%.0f%%", 100 * expected_exclusion_rate),
      ", the corresponding recruitment target was ", recruit_n, " participants."
    )
  )

  writeLines(writeup_lines, writeup_file)

  invisible(c(results_file, plot_file, summary_file, writeup_file))
}

cleaned_study_data <- setNames(vector("list", length(STUDY_IDS)), STUDY_IDS)

for (study_id in STUDY_IDS) {
  cat(
    "\n========== Running ",
    STUDY_CONFIG[[study_id]]$display_label,
    " (", study_id, ") ==========\n",
    sep = ""
  )
  cleaned_study_data[[study_id]] <- run_study(study_id)
}

if (RUN_DISPERSION_ANALYSIS) {
  cat("\n========== Running exploratory dispersion analyses ==========\n")
  dispersion_analysis_results <- run_exploratory_dispersion_analysis(
    cleaned_studies = cleaned_study_data,
    bootstrap_reps = DISPERSION_BOOTSTRAP_REPS,
    seed = 8122026
  )
}

if (RUN_DISTRIBUTION_ANALYSIS) {
  cat("\n========== Running combined distribution analyses ==========\n")
  run_distribution_analysis(cleaned_study_data)
}

if (RUN_POWER_ANALYSIS) {
  cat("\n========== Running simulation-based power analysis ==========\n")
  run_simulation_power_analysis()
}

cat("\nAll requested analyses completed.\n")
