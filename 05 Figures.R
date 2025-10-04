########################0 SETUP####################################################################
setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Scripts/CPRD-Thijs-GLP1-KF/")
source("00 Setup.R")

setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Processed data/")

for (q in c("single_episodes", "no_insulin", "ncurrtx_upto_4", "all_3_exclusions")) {
  load(paste0(today, "_hrs_sens_", q, ".Rda"))
  assign(paste0("hrs_", q), hrs)
  rm(hrs)
}

load(paste0(today, "_hrs.Rda"))
load(paste0(today, "_hrs_fg_ow.Rda"))
load(paste0(today, "_factor_hrs.Rda"))

# set default studydrug variable
studydrug_var = paste0("studydrug", main)

load(paste0(today, "_t2d_glp1_imputed_data_withweights_studydrug", main, ".Rda"))

############################1 FOREST PLOT FOR MAIN OUTCOMES################################################################

# GLP1-RA trial meta-analysis HR: HR 0·81, 95% CI 0·72–0·92 (Lancet Diabetes Endocrinol. 2025 Jan;13(1):15-28.)

hrs <- hrs %>% filter(!grepl("interaction", variable))

n.studydrug.vars <- hrs %>% .$variable %>% as.factor() %>% nlevels()

drug_reference = levels(cohort[[studydrug_var]])[1]
drug_of_interest = levels(cohort[[studydrug_var]])[2]

# Collect datasets across outcomes
# subset to outcome k and the current studydrug variable
  temp <- hrs %>%
    filter(analysis == "ow",
           outcome %in% outcomes_per_drugclass,
           variable == paste0("studydrug", main)) %>%
    mutate(HR = ifelse(is.na(string), 1, HR),
           string = ifelse(is.na(string), "1.00 (ref.)", as.character(string))
    )
  
  # ensure these cols are factors as used previously
  temp[c("outcome","contrast","variable","analysis")] <-
    lapply(temp[c("outcome","contrast","variable","analysis")], factor)
  
  
  # drug label extracted from contrast (as before)
  temp$drug <- as.factor(sub(" vs.*", "", temp$contrast))
  
  temp_k <- temp %>%
      mutate(drug = factor(drug, levels = c("SGLT2 + DPP4/SU", "SGLT2 + GLP1")))
  
  
  # --- extract the reference group's nN for this outcome (first factor level) ---
  ref_level <- levels(temp_k$drug)[1]
  
  # pull unique nN for the ref level (take first non-NA if multiple)
  ref_nN_val <- temp_k %>%
    filter(drug == ref_level) %>%
    pull(nN) %>%
    unique() %>%
    .[!is.na(.)] 
  
  # safety: if nothing found, set NA_character_
  if (is.null(ref_nN_val) || length(ref_nN_val) == 0) ref_nN_val <- NA_character_
  
  # create treated/ref columns, keep only the treated rows (drop ref rows)
  temp_k <- temp_k %>%
    filter(drug != ref_level) %>%
    mutate(treated_nN = nN,
           ref_nN = ref_nN_val) 
  
  # heading row (with column titles for treated and reference)
  heading <- tibble(
    outcome = "Outcome",
    drug = "",
    HR = NA_real_, LB = NA_real_, UB = NA_real_,
    string = "HR (95% CI)",
    treated_nN = paste0(drug_of_interest, " (n/N)"),
    ref_nN = paste0(drug_reference, " (n/N)")
  )
  
  
  


plot_df <- bind_rows(heading, temp_k)

# Define pretty labels
outcome_labels <- c(
  "ckd_egfr40" = "≥40% eGFR decline/ESKD",
  "ckd_egfr50" = "≥50% eGFR decline/ESKD",
  "mace"       = "MACE (incl. CV death)",
  "hf"         = "Hospitalisation for HF"
)

xmin = 0.10
xmax = 2.5
# Apply mapping before building plot_df
plot_df <- plot_df %>%
  mutate(outcome = recode(outcome, !!!outcome_labels),
         # also update headings if you used outcome as heading text
         drug = ifelse(drug %in% outcomes_per_drugclass, outcome, drug))

# Set y order (reverse for top-down order)
plot_df$y_order <- rev(seq_len(nrow(plot_df)))

# Forest plot with tabulated extras
forest_plot <- ggplot(plot_df, aes(y = y_order)) +
  # Forest CI + point
  geom_linerange(aes(xmin = LB, xmax = UB), na.rm = TRUE) +
  geom_point(aes(x = HR), shape = 15, size = 3, na.rm = TRUE) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  
  # Outcome headings on the left
  geom_text(aes(x = 0.10, label = outcome,
                fontface = ifelse(HR %>% is.na(), "bold", "plain")),
            hjust = 0) +
  
  # Events/subjects in middle-left
  # Treated n/N
  geom_text(aes(x = 0.2, label = treated_nN,
                fontface = ifelse(is.na(HR), "bold", "plain")),
            hjust = 0) +
  
  # Reference n/N
  geom_text(aes(x = 0.325, label = ref_nN,
                fontface = ifelse(is.na(HR), "bold", "plain")),
            hjust = 0) +
  
  
  # HR text on right-hand side
  geom_text(aes(x = 1.7, label = string,
                fontface = ifelse(HR %>% is.na(), "bold", "plain")),
            hjust = 0) +
  
  # text to indicate which drug to favour
  annotate("text", x = .75,
           y = max(plot_df$y_order), fontface = "italic",
           label = paste0("Favours ", drug_of_interest %>% substr(9, nchar(drug_of_interest)))) +
  
  annotate("text", x = 1.3,
           y = max(plot_df$y_order), fontface = "italic",
           label = paste0("Favours ", drug_reference %>% substr(9, nchar(drug_reference)))) +
  
  scale_x_continuous(trans = "log10",
                     breaks = c(0.5, 0.75, 1, 1.5),
                     limits = c(xmin, xmax)) +
  scale_y_continuous(NULL, breaks = NULL) +
  labs(x = "HR (95% CI)", y = NULL) +
  theme_classic() +
  theme(
    panel.grid = element_blank(),
    axis.line   = element_blank(),
    axis.ticks  = element_blank(),   # kill ggplot ticks
    axis.text.y = element_blank(),
    axis.title.x = element_text(
      hjust = 0.72
    )
  ) +
  
  # custom axis line
  geom_segment(aes(x = 0.49, xend = 1.52, y = 0, yend = 0),
               inherit.aes = FALSE, linewidth = 0.4, color = "black") +
  
  # custom ticks
  geom_segment(
    data = data.frame(x = c(0.5, 0.75, 1, 1.5)),
    aes(x = x, xend = x, y = 0, yend = -0.2),
    inherit.aes = FALSE,
    linewidth = 0.4,
    color = "black"
  )


setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Output/")
tiff(paste0(today, "_HR_", main, "_main_outcomes.tiff"), width=12, height=length(outcomes_per_drugclass)*0.8, units = "in", res=800)
print(forest_plot)
dev.off()


############################2 FOREST PLOT FOR SAFETY OUTCOMES################################################################
# Collect datasets across outcomes
# subset to outcome k and the current studydrug variable
temp <- hrs %>%
  filter(analysis == "ow",
         outcome %in% safety_outcomes,
         variable == paste0("studydrug", main)) %>%
  mutate(HR = ifelse(is.na(string), 1, HR),
         string = ifelse(is.na(string), "1.00 (ref.)", as.character(string))
  )

# ensure these cols are factors as used previously
temp[c("outcome","contrast","variable","analysis")] <-
  lapply(temp[c("outcome","contrast","variable","analysis")], factor)


# drug label extracted from contrast (as before)
temp$drug <- as.factor(sub(" vs.*", "", temp$contrast))

temp_k <- temp %>%
  mutate(drug = factor(drug, levels = c("SGLT2 + DPP4/SU", "SGLT2 + GLP1")))


# --- extract the reference group's nN for this outcome (first factor level) ---
ref_level <- levels(temp_k$drug)[1]

# pull unique nN for the ref level (take first non-NA if multiple)
ref_nN_val <- temp_k %>%
  filter(drug == ref_level) %>%
  pull(nN) %>%
  unique() %>%
  .[!is.na(.)] 

# safety: if nothing found, set NA_character_
if (is.null(ref_nN_val) || length(ref_nN_val) == 0) ref_nN_val <- NA_character_

# create treated/ref columns, keep only the treated rows (drop ref rows)
temp_k <- temp_k %>%
  filter(drug != ref_level) %>%
  mutate(treated_nN = nN,
         ref_nN = ref_nN_val) 

# heading row (with column titles for treated and reference)
heading <- tibble(
  outcome = "Outcome",
  drug = "",
  HR = NA_real_, LB = NA_real_, UB = NA_real_,
  string = "HR (95% CI)",
  treated_nN = paste0(drug_of_interest, " (n/N)"),
  ref_nN = paste0(drug_reference, " (n/N)")
)





plot_df <- bind_rows(heading, temp_k)

# Define pretty labels
outcome_labels <- c(
  "acutepancreatitis"      = "Acute pancreatitis",
  "retinopathy"            = "Diabetic retinopathy",
  "lowerlimbfracture"      = "Lower limb fracture"
)

xmin = 0.10
xmax = 2.10
# Apply mapping before building plot_df
plot_df <- plot_df %>%
  mutate(outcome = recode(outcome, !!!outcome_labels),
         # also update headings if you used outcome as heading text
         drug = ifelse(drug %in% safety_outcomes, outcome, drug))

# Set y order (reverse for top-down order)
plot_df$y_order <- rev(seq_len(nrow(plot_df)))

# Forest plot with tabulated extras
forest_plot <- ggplot(plot_df, aes(y = y_order)) +
  # Forest CI + point
  geom_linerange(aes(xmin = LB, xmax = UB), na.rm = TRUE) +
  geom_point(aes(x = HR), shape = 15, size = 3, na.rm = TRUE) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  
  # Outcome headings on the left
  geom_text(aes(x = 0.10, label = outcome,
                fontface = ifelse(HR %>% is.na(), "bold", "plain")),
            hjust = 0) +
  
  # Events/subjects in middle-left
  # Treated n/N
  geom_text(aes(x = 0.2, label = treated_nN,
                fontface = ifelse(is.na(HR), "bold", "plain")),
            hjust = 0) +
  
  # Reference n/N
  geom_text(aes(x = 0.325, label = ref_nN,
                fontface = ifelse(is.na(HR), "bold", "plain")),
            hjust = 0) +
  
  
  # HR text on right-hand side
  geom_text(aes(x = 1.7, label = string,
                fontface = ifelse(HR %>% is.na(), "bold", "plain")),
            hjust = 0) +
  
  # text to indicate which drug to favour
  annotate("text", x = .75,
           y = max(plot_df$y_order), fontface = "italic",
           label = paste0("Favours ", drug_of_interest %>% substr(9, nchar(drug_of_interest)))) +
  
  annotate("text", x = 1.3,
           y = max(plot_df$y_order), fontface = "italic",
           label = paste0("Favours ", drug_reference %>% substr(9, nchar(drug_reference)))) +
  
  scale_x_continuous(trans = "log10",
                     breaks = c(0.5, 0.75, 1, 1.5),
                     limits = c(xmin, xmax)) +
  scale_y_continuous(NULL, breaks = NULL) +
  labs(x = "HR (95% CI)", y = NULL) +
  theme_classic() +
  theme(
    panel.grid = element_blank(),
    axis.line   = element_blank(),
    axis.ticks  = element_blank(),   # kill ggplot ticks
    axis.text.y = element_blank(),
    axis.title.x = element_text(
      hjust = 0.72
    )
  ) +
  
  # custom axis line
  geom_segment(aes(x = 0.49, xend = 1.52, y = 0, yend = 0),
               inherit.aes = FALSE, linewidth = 0.4, color = "black") +
  
  # custom ticks
  geom_segment(
    data = data.frame(x = c(0.5, 0.75, 1, 1.5)),
    aes(x = x, xend = x, y = 0, yend = -0.2),
    inherit.aes = FALSE,
    linewidth = 0.4,
    color = "black"
  )



setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Output/")
tiff(paste0(today, "_HR_", main, "_safety_outcomes.tiff"), width=12, height=length(safety_outcomes)*0.8, units = "in", res=800)
print(forest_plot)
dev.off()


############################3 FOREST PLOT BY FACTORS################################################################

drug_reference = levels(cohort[[studydrug_var]])[1]
drug_of_interest = levels(cohort[[studydrug_var]])[2]

factors <- c("malesex", "white_ethnicity", "predrug_cvd", "predrug_heartfailure", "age_cat")

factor_hrs$condition_label <- case_when(
  factor_hrs$factor == "malesex" & factor_hrs$condition == TRUE ~ "Male",
  factor_hrs$factor == "malesex" & factor_hrs$condition == FALSE ~ "Female",
  
  factor_hrs$factor == "white_ethnicity" & factor_hrs$condition == TRUE ~ "White",
  factor_hrs$factor == "white_ethnicity" & factor_hrs$condition == FALSE ~ "Non-white",
  
  factor_hrs$factor == "predrug_cvd" & factor_hrs$condition == TRUE ~ "Present",
  factor_hrs$factor == "predrug_cvd" & factor_hrs$condition == FALSE ~ "Absent",
  
  factor_hrs$factor == "predrug_heartfailure" & factor_hrs$condition == TRUE ~ "Present",
  factor_hrs$factor == "predrug_heartfailure" & factor_hrs$condition == FALSE ~ "Absent",
  
  factor_hrs$factor == "age_cat" ~ factor_hrs$condition,
  
  TRUE ~ NA_character_
) 

# Filter to the CKD outcome only
temp <- factor_hrs %>%
  filter(analysis == "ow",
         outcome == "ckd_egfr40",
         variable == paste0("studydrug", main)) %>%
  mutate(HR = ifelse(is.na(string), 1, HR),
         string = ifelse(is.na(string), "1.00 (ref.)", as.character(string)))

temp[c("outcome","contrast","variable","analysis")] <-
  lapply(temp[c("outcome","contrast","variable","analysis")], factor)

# Extract drug labels
temp$drug <- as.factor(sub(" vs.*", "", temp$contrast))
temp <- temp %>%
  mutate(drug = factor(drug, levels = c(drug_reference, drug_of_interest)))

# Get ref n/N (same as main code)
ref_nN_val <- temp %>%
  filter(drug == drug_reference) %>%
  pull(nN) %>% unique() %>% .[!is.na(.)]
if (length(ref_nN_val) == 0) ref_nN_val <- NA_character_

# Keep only treated rows and add ref_nN
plot_df <- temp %>%
  filter(drug != drug_reference) %>%
  mutate(treated_nN = nN,
         ref_nN = ref_nN_val,
         factor_label = recode(factor,
                               "malesex" = "Sex",
                               "white_ethnicity" = "Ethnicity",
                               "predrug_cvd" = "ASCVD",
                               "predrug_heartfailure" = "Heart failure",
                               "age_cat" = "Age"
         ))

# Create the global header row 
header <- tibble(
  outcome = "Header",
  factor_label = " Header",
  condition_label = "≥40% eGFR decline/ESKD",
  HR = NA_real_, LB = NA_real_, UB = NA_real_,
  string = "HR (95% CI)",
  treated_nN = paste0(drug_of_interest, " (n/N)"),
  ref_nN = paste0(drug_reference, " (n/N)")
)

# Create subgroup heading rows (no n/N titles here) 
heading_rows <- plot_df %>%
  group_by(factor_label) %>%
  summarise() %>%
  mutate(
    HR = NA_real_, LB = NA_real_, UB = NA_real_,
    string = "",
    treated_nN = "",
    ref_nN = "",
    condition_label = factor_label,
    condition_label = paste0(
      factor_label, 
      " (p=", 
      round(factor_hrs$p_value_interaction[match(factor_label, recode(
        factor_hrs$factor,
        "malesex" = "Sex",
        "white_ethnicity" = "Ethnicity",
        "predrug_cvd" = "ASCVD",
        "predrug_heartfailure" = "Heart failure",
        "age_cat" = "Age"
      ))], 2), 
      ")"
    )
  )

# Bind together: header row first, then headings + subgroup rows 
plot_df <- bind_rows(header, heading_rows, plot_df)

custom_order <- c(" Header", "Age", "Sex", "Ethnicity", "ASCVD", "Heart failure")

# Convert factor_label to a factor with levels in the desired order
plot_df$factor_label <- factor(plot_df$factor_label, levels = custom_order)
plot_df <- plot_df %>%
  arrange(factor_label, desc(HR %>% is.na()))
# 
# # Reorder y axis
# plot_df$y_order <- rev(seq_len(nrow(plot_df)))

age_heading <- plot_df$condition_label[grepl("^Age", plot_df$condition_label)]

# Reorder Age conditions only
plot_df <-  plot_df %>%
  mutate(
    condition_label = case_when(
      grepl("^Age", condition_label) ~ factor(condition_label, levels = c(age_heading, "< 50", "50 - 60", "60 - 70", "≥ 70")),
      TRUE ~ condition_label
    )
  )

# Arrange by factor_label (groups) and HR NA (header first) 
# For Age, the condition_label factor will control the order within Age
plot_df <- plot_df %>%
  arrange(factor_label, desc(is.na(HR)), condition_label)

# Reassign y_order
plot_df$y_order <- rev(seq_len(nrow(plot_df)))

## manually sort order:
plot_df$y_order = ifelse(
  plot_df$condition_label == "< 50", 16, 
  ifelse(plot_df$condition_label == "50 - 60", 15, 
         ifelse(plot_df$condition_label == "60 - 70", 14, plot_df$y_order))
)

# plot
forest_plot <- ggplot(plot_df, aes(y = y_order)) +
  geom_linerange(aes(xmin = LB, xmax = UB), na.rm = TRUE) +
  geom_point(aes(x = HR), shape = 15, size = 3, na.rm = TRUE) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  
  # Left labels (factor headings + subgroup names)
  geom_text(aes(x = 0.1, label = condition_label,
                fontface = ifelse(is.na(HR) , "bold", "plain")),
            hjust = 0) +
  
  # n/N columns (only populated for header and subgroup rows)
  geom_text(aes(x = 0.175, label = treated_nN,
                fontface = ifelse(is.na(HR) & string != "", "bold", "plain")),
            hjust = 0) +
  geom_text(aes(x = 0.27, label = ref_nN,
                fontface = ifelse(is.na(HR) & string != "", "bold", "plain")),
            hjust = 0) +
  
  # text to indicate which drug to favour
  annotate("text", x = .65,
           y = max(plot_df$y_order), fontface = "italic",
           label = paste0("Favours ", drug_of_interest %>% substr(9, nchar(drug_of_interest)))) +
  
  annotate("text", x = 1.5,
           y = max(plot_df$y_order), fontface = "italic",
           label = paste0("Favours ", drug_reference %>% substr(9, nchar(drug_reference)))) +
  
  # HR text on right
  geom_text(aes(x = 2.3, label = string,
                fontface = ifelse(is.na(HR) & string != "", "bold", "plain")),
            hjust = 0) +
  
  scale_x_continuous(trans = "log10",
                     breaks = c(0.35, 0.5, 0.75, 1, 1.5, 2.25),
                     limits = c(0.10, 2.7)) +
  scale_y_continuous(NULL, breaks = NULL) +
  labs(x = "HR (95% CI)", y = NULL) +
  theme_classic() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.x = element_blank(),
        axis.title.x = element_text(hjust = 0.72)) +
  
  # Custom axis line + ticks
  geom_segment(aes(x = 0.35, xend = 2.25, y = 0, yend = 0),
               inherit.aes = FALSE, linewidth = 0.4, color = "black") +
  geom_segment(data = data.frame(x = c(0.35, 0.5, 0.75, 1, 1.5, 2.25)),
               aes(x = x, xend = x, y = 0, yend = -0.2),
               inherit.aes = FALSE, linewidth = 0.4, color = "black")


setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Output/")
tiff(paste0(today, "_HR_", main, "_main_outcome_by_factors.tiff"), width=14, height=6, units = "in", res=800)
print(forest_plot)
dev.off()

############################4 FOREST PLOT FOR DPP4/SU################################################################

m = 1
studydrug_var = paste0("studydrug", m)

drug_reference = levels(cohort[[studydrug_var]])[1]
drug_of_interest = levels(cohort[[studydrug_var]])[2]

# Collect datasets across outcomes
# subset to outcome k and the current studydrug variable
temp <- hrs %>%
  filter(analysis == "ow",
         outcome %in% (outcomes %>% setdiff("death")),
         variable == paste0("studydrug", m)) %>%
  mutate(HR = ifelse(is.na(string), 1, HR),
         string = ifelse(is.na(string), "1.00 (ref.)", as.character(string))
  )

# ensure these cols are factors as used previously
temp[c("outcome","contrast","variable","analysis")] <-
  lapply(temp[c("outcome","contrast","variable","analysis")], factor)


# drug label extracted from contrast (as before)
temp$drug <- as.factor(sub(" vs.*", "", temp$contrast))

temp_k <- temp %>%
  filter(!grepl(drug, "GLP1")) %>%
  mutate(drug = factor(drug, levels = c("SGLT2 + SU", "SGLT2 + DPP4")))


# --- extract the reference group's nN for this outcome (first factor level) ---
ref_level <- levels(temp_k$drug)[1]

# pull unique nN for the ref level (take first non-NA if multiple)
ref_nN_val <- temp_k %>%
  filter(drug == ref_level) %>%
  pull(nN) %>%
  unique() %>%
  .[!is.na(.)] 

# safety: if nothing found, set NA_character_
if (is.null(ref_nN_val) || length(ref_nN_val) == 0) ref_nN_val <- NA_character_

# create treated/ref columns, keep only the treated rows (drop ref rows)
temp_k <- temp_k %>%
  filter(drug != ref_level) %>%
  mutate(treated_nN = nN,
         ref_nN = ref_nN_val) 

# heading row (with column titles for treated and reference)
heading <- tibble(
  outcome = "Outcome",
  drug = "",
  HR = NA_real_, LB = NA_real_, UB = NA_real_,
  string = "HR (95% CI)",
  treated_nN = paste0(drug_of_interest, " (n/N)"),
  ref_nN = paste0(drug_reference, " (n/N)")
)





plot_df <- bind_rows(heading, temp_k)

# Define pretty labels
outcome_labels <- c(
  "ckd_egfr40" = "≥40% eGFR decline/ESKD",
  "ckd_egfr50" = "≥50% eGFR decline/ESKD",
  "mace"       = "MACE (incl. CV death)",
  "hf"         = "Hospitalisation for HF",
  "acutepancreatitis"      = "Acute pancreatitis",
  "retinopathy"            = "Diabetic retinopathy",
  "lowerlimbfracture"      = "Lower limb fracture"
)

xmin = 0.10
xmax = 2.5
# Apply mapping before building plot_df
plot_df <- plot_df %>%
  mutate(outcome = recode(outcome, !!!outcome_labels),
         # also update headings if you used outcome as heading text
         drug = ifelse(drug %in% outcomes, outcome, drug))

# Set y order (reverse for top-down order)
plot_df$y_order <- rev(seq_len(nrow(plot_df)))

# Forest plot with tabulated extras
forest_plot <- ggplot(plot_df, aes(y = y_order)) +
  # Forest CI + point
  geom_linerange(aes(xmin = LB, xmax = UB), na.rm = TRUE) +
  geom_point(aes(x = HR), shape = 15, size = 3, na.rm = TRUE) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  
  # Outcome headings on the left
  geom_text(aes(x = 0.10, label = outcome,
                fontface = ifelse(HR %>% is.na(), "bold", "plain")),
            hjust = 0) +
  
  # Events/subjects in middle-left
  # Treated n/N
  geom_text(aes(x = 0.2, label = treated_nN,
                fontface = ifelse(is.na(HR), "bold", "plain")),
            hjust = 0) +
  
  # Reference n/N
  geom_text(aes(x = 0.325, label = ref_nN,
                fontface = ifelse(is.na(HR), "bold", "plain")),
            hjust = 0) +
  
  
  # HR text on right-hand side
  geom_text(aes(x = 2.0, label = string,
                fontface = ifelse(HR %>% is.na(), "bold", "plain")),
            hjust = 0) +
  
  # text to indicate which drug to favour
  annotate("text", x = .75,
           y = max(plot_df$y_order), fontface = "italic",
           label = paste0("Favours ", drug_of_interest %>% substr(9, nchar(drug_of_interest)))) +
  
  annotate("text", x = 1.3,
           y = max(plot_df$y_order), fontface = "italic",
           label = paste0("Favours ", drug_reference %>% substr(9, nchar(drug_reference)))) +
  
  scale_x_continuous(trans = "log10",
                     breaks = c(0.5, 0.75, 1, 1.5, 2.0),
                     limits = c(xmin, xmax)) +
  scale_y_continuous(NULL, breaks = NULL) +
  labs(x = "HR (95% CI)", y = NULL) +
  theme_classic() +
  theme(
    panel.grid = element_blank(),
    axis.line   = element_blank(),
    axis.ticks  = element_blank(),   # kill ggplot ticks
    axis.text.y = element_blank(),
    axis.title.x = element_text(
      hjust = 0.72
    )
  ) +
  
  # custom axis line
  geom_segment(aes(x = 0.45, xend = 2.1, y = 0, yend = 0),
               inherit.aes = FALSE, linewidth = 0.4, color = "black") +
  
  # custom ticks
  geom_segment(
    data = data.frame(x = c(0.5, 0.75, 1, 1.5, 2.0)),
    aes(x = x, xend = x, y = 0, yend = -0.2),
    inherit.aes = FALSE,
    linewidth = 0.4,
    color = "black"
  )


setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Output/")
tiff(paste0(today, "_HR_", m, "_all_outcomes.tiff"), width=12, height=length(outcomes %>% setdiff("death"))*0.8, units = "in", res=800)
print(forest_plot)
dev.off()

############################5 FOREST PLOT BY GLP1 TYPE################################################################

m = 3
studydrug_var = paste0("studydrug", m)

drug_reference = "DPP4/SU"
drug_of_interest = "GLP1 subgroup"

# Collect datasets across outcomes
# subset to outcome k and the current studydrug variable
temp <- hrs %>%
  filter(analysis == "ow",
         outcome %in% (outcomes %>% setdiff("death")),
         variable == paste0("studydrug", m)) %>%
  mutate(HR = ifelse(is.na(string), 1, HR),
         string = ifelse(is.na(string), "1.00 (ref.)", as.character(string))
  )

# ensure these cols are factors as used previously
temp[c("outcome","contrast","variable","analysis")] <-
  lapply(temp[c("outcome","contrast","variable","analysis")], factor)


# drug label extracted from contrast (as before)
temp$drug <- as.factor(sub(" vs.*", "", temp$contrast))

temp_k <- temp %>%
  mutate(drug = factor(drug, levels = c("SGLT2 + DPP4/SU", "GLP1 with direct kidney outcome evidence", "Other GLP1")))


# --- extract the reference group's nN for this outcome (first factor level) ---
ref_level <- levels(temp_k$drug)[1]

# pull unique nN for the ref level (take first non-NA if multiple)
ref_vals <- temp %>%
  filter(drug == ref_level) %>%
  select(outcome, ref_nN = nN)

# Keep non-reference drugs, add their own nN and join correct ref_nN
temp_k <- temp %>%
  filter(drug != ref_level) %>%
  mutate(treated_nN = nN) %>%
  left_join(ref_vals, by = "outcome")


# --- 1. Create the global header row ---
header <- tibble(
  outcome = " ",
  drug = " ",
  HR = NA_real_, LB = NA_real_, UB = NA_real_,
  string = "HR (95% CI)",
  treated_nN = paste0(drug_of_interest, " (n/N)"),
  ref_nN = paste0(drug_reference, " (n/N)")
)

# --- 2. Create subgroup heading rows (no n/N titles here) ---
heading_rows <- temp_k %>%
  group_by(outcome) %>%
  summarise() %>%
  mutate(
    drug = outcome,
    HR = NA_real_, LB = NA_real_, UB = NA_real_,
    string = "",
    treated_nN = "",
    ref_nN = ""
  )

plot_df <- bind_rows(header, heading_rows, temp_k)

# Define pretty labels
outcome_labels <- c(
  "ckd_egfr40" = "≥40% eGFR decline/ESKD",
  "ckd_egfr50" = "≥50% eGFR decline/ESKD",
  "mace"       = "MACE (incl. CV death)",
  "hf"         = "Hospitalisation for HF",
  "acutepancreatitis"      = "Acute pancreatitis",
  "retinopathy"            = "Diabetic retinopathy",
  "lowerlimbfracture"      = "Lower limb fracture"
)

xmin = 0.03
xmax = 2.5

# Apply mapping before building plot_df
plot_df <- plot_df %>%
  mutate(outcome = recode(outcome, !!!outcome_labels),
         # also update headings if you used outcome as heading text
         drug = ifelse(drug %in% outcomes, outcome, drug))

# Define desired outcome order: " " first, then outcomes in  chosen order
desired_order <- c(" ", recode(outcomes, !!!outcome_labels))

plot_df <- plot_df %>%
  mutate(
    outcome = factor(outcome, levels = desired_order),
    # tag rows: 0 = header row, 1 = outcome heading, 2 = drug row
    row_type = case_when(
      outcome == " " ~ 0,
      is.na(HR) & string == "" ~ 1,
      TRUE ~ 2
    )
  ) %>%
  arrange(outcome, row_type) %>%
  select(-row_type)

# Now assign y_order
plot_df$y_order <- rev(seq_len(nrow(plot_df)))



# Forest plot with tabulated extras
forest_plot <- ggplot(plot_df, aes(y = y_order)) +
  # Forest CI + point
  geom_linerange(aes(xmin = LB, xmax = UB), na.rm = TRUE) +
  geom_point(aes(x = HR), shape = 15, size = 3, na.rm = TRUE) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  
  # Outcome headings on the left
  geom_text(aes(x = xmin, label = drug,
                fontface = ifelse(HR %>% is.na(), "bold", "plain")),
            hjust = 0) +
  
  # Events/subjects in middle-left
  # Treated n/N
  geom_text(aes(x = 0.10, label = treated_nN,
                fontface = ifelse(is.na(HR), "bold", "plain")),
            hjust = 0) +
  
  # Reference n/N
  geom_text(aes(x = 0.20, label = ref_nN,
                fontface = ifelse(is.na(HR), "bold", "plain")),
            hjust = 0) +
  
  
  # HR text on right-hand side
  geom_text(aes(x = 1.9, label = string,
                fontface = ifelse(HR %>% is.na(), "bold", "plain")),
            hjust = 0) +
  
  # text to indicate which drug to favour
  annotate("text", x = .65,
           y = max(plot_df$y_order), fontface = "italic",
           label = paste0("Favours ", drug_of_interest)) +
  
  annotate("text", x = 1.4,
           y = max(plot_df$y_order), fontface = "italic",
           label = paste0("Favours ", drug_reference)) +
  
  scale_x_continuous(trans = "log10",
                     breaks = c(0.3, 0.5, 0.75, 1, 1.5, 1.75),
                     limits = c(xmin, xmax)) +
  scale_y_continuous(NULL, breaks = NULL) +
  labs(x = "HR (95% CI)", y = NULL) +
  theme_classic() +
  theme(
    panel.grid = element_blank(),
    axis.line   = element_blank(),
    axis.ticks  = element_blank(),   # kill ggplot ticks
    axis.text.y = element_blank(),
    axis.title.x = element_text(
      hjust = 0.72
    )
  ) +
  
  # custom axis line
  geom_segment(aes(x = 0.3, xend = 1.75, y = 0, yend = 0),
               inherit.aes = FALSE, linewidth = 0.4, color = "black") +
  
  # custom ticks
  geom_segment(
    data = data.frame(x = c(0.3, 0.5, 0.75, 1, 1.5, 1.75)),
    aes(x = x, xend = x, y = 0, yend = -0.2),
    inherit.aes = FALSE,
    linewidth = 0.4,
    color = "black"
  )


setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Output/")
tiff(paste0(today, "_HR_", m, "_all_outcomes.tiff"), width=12, height=length(outcomes %>% setdiff("death"))*0.8, units = "in", res=800)
print(forest_plot)
dev.off()


############################6 FOREST PLOT FOR SENSITIVITY ANALYSES################################################################


# set default studydrug variable
studydrug_var = paste0("studydrug", main)

drug_reference = levels(cohort[[studydrug_var]])[1]
drug_of_interest = levels(cohort[[studydrug_var]])[2]

k = "ckd_egfr40"

temp <- hrs %>%
  filter(analysis == "ow",
         outcome == k,
         variable == paste0("studydrug", main)) %>%
  mutate(HR = ifelse(is.na(string), 1, HR),
         string = ifelse(is.na(string), "1.00 (ref.)", as.character(string)),
         analysis = "all"
  )

for (q in c("single_episodes", "no_insulin", "ncurrtx_upto_4", "all_3_exclusions")) {
  hrs_name = get(paste0("hrs_", q))
  temp <- temp %>% rbind(hrs_name %>%
                           filter(analysis == "ow", 
                                             outcome == k,
                                             variable == paste0("studydrug", main)) %>% mutate(analysis = q) ) 
}

temp[c("outcome","contrast","variable","analysis")] <-
  lapply(temp[c("outcome","contrast","variable","analysis")], factor)


# drug label extracted from contrast (as before)
temp$drug <- as.factor(sub(" vs.*", "", temp$contrast))

temp_k <- temp %>%
  mutate(drug = factor(drug, levels = c("SGLT2 + DPP4/SU", "SGLT2 + GLP1")))


# --- extract the reference group's nN for this outcome (first factor level) ---
ref_level <- levels(temp_k$drug)[1]

# pull unique nN for the ref level (take first non-NA if multiple)
ref_nN_val <- temp_k %>%
  filter(drug == ref_level) %>%
  pull(nN) %>%
  unique() %>%
  .[!is.na(.)] 

# safety: if nothing found, set NA_character_
if (is.null(ref_nN_val) || length(ref_nN_val) == 0) ref_nN_val <- NA_character_

# create treated/ref columns, keep only the treated rows (drop ref rows)
temp_k <- temp_k %>%
  filter(drug != ref_level) %>%
  mutate(treated_nN = nN,
         ref_nN = ref_nN_val) 

# heading row (with column titles for treated and reference)
heading <- tibble(
  outcome = k,
  analysis = "≥40% eGFR decline/ESKD",
  drug = "",
  HR = NA_real_, LB = NA_real_, UB = NA_real_,
  string = "HR (95% CI)",
  treated_nN = paste0(drug_of_interest, " (n/N)"),
  ref_nN = paste0(drug_reference, " (n/N)")
)





plot_df <- bind_rows(heading, temp_k)

# Define pretty labels
analysis_labels <- c(
  "all" = "All individuals",
  "single_episodes" = "Excluding GLP1 episodes if\npreceding DPP4/SU episode",
  "no_insulin"       = "Excluding individuals treated\nwith insulin",
  "ncurrtx_upto_4"         = "Excluding individuals treated with\n> 4 glucose-lowering treatments",
  "all_3_exclusions"       = "Excluding all of the above"
)

xmin = 0.10
xmax = 2.5
# Apply mapping before building plot_df
plot_df <- plot_df %>%
  mutate(analysis = recode(analysis, !!!analysis_labels))

# Set y order (reverse for top-down order)
plot_df$y_order <- rev(seq_len(nrow(plot_df)))

# Forest plot with tabulated extras
forest_plot <- ggplot(plot_df, aes(y = y_order)) +
  # Forest CI + point
  geom_linerange(aes(xmin = LB, xmax = UB), na.rm = TRUE) +
  geom_point(aes(x = HR), shape = 15, size = 3, na.rm = TRUE) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  
  # Outcome headings on the left
  geom_text(aes(x = 0.10, label = analysis,
                fontface = ifelse(HR %>% is.na(), "bold", "plain")),
            hjust = 0) +
  
  # Events/subjects in middle-left
  # Treated n/N
  geom_text(aes(x = 0.2, label = treated_nN,
                fontface = ifelse(is.na(HR), "bold", "plain")),
            hjust = 0) +
  
  # Reference n/N
  geom_text(aes(x = 0.325, label = ref_nN,
                fontface = ifelse(is.na(HR), "bold", "plain")),
            hjust = 0) +
  
  
  # HR text on right-hand side
  geom_text(aes(x = 1.7, label = string,
                fontface = ifelse(HR %>% is.na(), "bold", "plain")),
            hjust = 0) +
  
  # text to indicate which drug to favour
  annotate("text", x = .75,
           y = max(plot_df$y_order), fontface = "italic",
           label = paste0("Favours ", drug_of_interest %>% substr(9, nchar(drug_of_interest)))) +
  
  annotate("text", x = 1.3,
           y = max(plot_df$y_order), fontface = "italic",
           label = paste0("Favours ", drug_reference %>% substr(9, nchar(drug_reference)))) +
  
  scale_x_continuous(trans = "log10",
                     breaks = c(0.5, 0.75, 1, 1.5),
                     limits = c(xmin, xmax)) +
  scale_y_continuous(NULL, breaks = NULL) +
  labs(x = "HR (95% CI)", y = NULL) +
  theme_classic() +
  theme(
    panel.grid = element_blank(),
    axis.line   = element_blank(),
    axis.ticks  = element_blank(),   # kill ggplot ticks
    axis.text.y = element_blank(),
    axis.title.x = element_text(
      hjust = 0.72
    )
  ) +
  
  # custom axis line
  geom_segment(aes(x = 0.49, xend = 1.52, y = 0, yend = 0),
               inherit.aes = FALSE, linewidth = 0.4, color = "black") +
  
  # custom ticks
  geom_segment(
    data = data.frame(x = c(0.5, 0.75, 1, 1.5)),
    aes(x = x, xend = x, y = 0, yend = -0.2),
    inherit.aes = FALSE,
    linewidth = 0.4,
    color = "black"
  )


setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Output/")
tiff(paste0(today, "_HR_", main, "_sensitivity_analysis.tiff"), width=12, height=length(analysis_labels)*1.0, units = "in", res=800)
print(forest_plot)
dev.off()


############################7 FOREST PLOT FOR DIFFERENT ANALYSIS APPROACHES################################################################


# set default studydrug variable
studydrug_var = paste0("studydrug", main)

drug_reference = levels(cohort[[studydrug_var]])[1]
drug_of_interest = levels(cohort[[studydrug_var]])[2]

k = "ckd_egfr40"

temp <- hrs %>% rbind(hrs_fg_ow) %>%
  filter(analysis != "unadj",
         outcome == k,
         variable == paste0("studydrug", main)) %>%
  mutate(HR = ifelse(is.na(string), 1, HR),
         string = ifelse(is.na(string), "1.00 (ref.)", as.character(string))
  )

temp[c("outcome","contrast","variable","analysis")] <-
  lapply(temp[c("outcome","contrast","variable","analysis")], factor)


# drug label extracted from contrast (as before)
temp$drug <- as.factor(sub(" vs.*", "", temp$contrast))

temp_k <- temp %>%
  mutate(drug = factor(drug, levels = c("SGLT2 + DPP4/SU", "SGLT2 + GLP1")))


# --- extract the reference group's nN for this outcome (first factor level) ---
ref_level <- levels(temp_k$drug)[1]

# pull unique nN for the ref level (take first non-NA if multiple)
ref_nN_val <- temp_k %>%
  filter(drug == ref_level) %>%
  pull(nN) %>%
  trimws() %>%
  unique() %>%
  .[!is.na(.)] 

ref_nN_val <- paste("    ", ref_nN_val)

# safety: if nothing found, set NA_character_
if (is.null(ref_nN_val) || length(ref_nN_val) == 0) ref_nN_val <- NA_character_

# create treated/ref columns, keep only the treated rows (drop ref rows)
temp_k <- temp_k %>%
  filter(drug != ref_level) %>%
  mutate(treated_nN = paste("    ", trimws(nN)),
         ref_nN = ref_nN_val) 

# heading row (with column titles for treated and reference)
heading <- tibble(
  outcome = k,
  analysis = "≥40% eGFR decline/ESKD",
  drug = "",
  HR = NA_real_, LB = NA_real_, UB = NA_real_,
  string = "HR (95% CI)",
  treated_nN = paste0(drug_of_interest, " (n/N)"),
  ref_nN = paste0(drug_reference, " (n/N)")
)





plot_df <- bind_rows(heading, temp_k)

# Define pretty labels
analysis_labels <- c(
  "ow" = "Overlap-weighting",
  "iptw" = "Inverse probability of\ntreatment weighting",
  "adj"       = "Multivariable adjustment only",
  "fg_ow" = "Fine-Gray competing risk analysis\n(overlap-weighted)"
)

xmin = 0.10
xmax = 2.5
# Apply mapping before building plot_df
plot_df <- plot_df %>%
  mutate(analysis = recode(analysis, !!!analysis_labels))

# Set y order (reverse for top-down order)
plot_df$y_order <- rev(seq_len(nrow(plot_df)))

# Forest plot with tabulated extras
forest_plot <- ggplot(plot_df, aes(y = y_order)) +
  # Forest CI + point
  geom_linerange(aes(xmin = LB, xmax = UB), na.rm = TRUE) +
  geom_point(aes(x = HR), shape = 15, size = 3, na.rm = TRUE) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  
  # Outcome headings on the left
  geom_text(aes(x = 0.10, label = analysis,
                fontface = ifelse(HR %>% is.na(), "bold", "plain")),
            hjust = 0) +
  
  # Events/subjects in middle-left
  # Treated n/N
  geom_text(aes(x = 0.2, label = treated_nN,
                fontface = ifelse(is.na(HR), "bold", "plain")),
            hjust = 0) +
  
  # Reference n/N
  geom_text(aes(x = 0.325, label = ref_nN,
                fontface = ifelse(is.na(HR), "bold", "plain")),
            hjust = 0) +
  
  
  # HR text on right-hand side
  geom_text(aes(x = 1.7, label = string,
                fontface = ifelse(HR %>% is.na(), "bold", "plain")),
            hjust = 0) +
  
  # text to indicate which drug to favour
  annotate("text", x = .75,
           y = max(plot_df$y_order), fontface = "italic",
           label = paste0("Favours ", drug_of_interest %>% substr(9, nchar(drug_of_interest)))) +
  
  annotate("text", x = 1.3,
           y = max(plot_df$y_order), fontface = "italic",
           label = paste0("Favours ", drug_reference %>% substr(9, nchar(drug_reference)))) +
  
  scale_x_continuous(trans = "log10",
                     breaks = c(0.5, 0.75, 1, 1.5),
                     limits = c(xmin, xmax)) +
  scale_y_continuous(NULL, breaks = NULL) +
  labs(x = "HR (95% CI)", y = NULL) +
  theme_classic() +
  theme(
    panel.grid = element_blank(),
    axis.line   = element_blank(),
    axis.ticks  = element_blank(),   # kill ggplot ticks
    axis.text.y = element_blank(),
    axis.title.x = element_text(
      hjust = 0.72
    )
  ) +
  
  # custom axis line
  geom_segment(aes(x = 0.49, xend = 1.52, y = 0, yend = 0),
               inherit.aes = FALSE, linewidth = 0.4, color = "black") +
  
  # custom ticks
  geom_segment(
    data = data.frame(x = c(0.5, 0.75, 1, 1.5)),
    aes(x = x, xend = x, y = 0, yend = -0.2),
    inherit.aes = FALSE,
    linewidth = 0.4,
    color = "black"
  )


setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Output/")
tiff(paste0(today, "_HR_", main, "_different_analysis_approaches.tiff"), width=12, height=length(analysis_labels)*1.0, units = "in", res=800)
print(forest_plot)
dev.off()


############################8 SPLINE PLOTS################################################################

cohort <- cohort %>% 
  group_by(.imp, patid, !!sym(studydrug_var)) %>% 
  arrange(dstartdate) %>% 
  filter(!duplicated(!!sym(studydrug_var))) %>% 
  ungroup()

ddist <- cohort %>% datadist()
options(datadist = "ddist") 


numerical_covariates <- c("dstartdate_age", "prebmi", "preegfr", "uacr", "prehba1c", "ckdpc_40egfr_score")

var_labels <- c(
  dstartdate_age = "Age (years)",
  prehba1c = "HbA1c (mmol/mol)",
  prebmi = "BMI (kg/m²)",
  preegfr = "eGFR (mL/min/1.73m²)",
  uacr = "uACR (mg/mmol)",
  presbp = "Systolic blood pressure (mmHg)",
  ckdpc_40egfr_score = "CKD-PC risk score"
)



for (k in outcomes_per_drugclass) {
  print(paste0("Spline plots for outcome ", k))
  
  HR <- (hrs %>% filter(contrast == paste0(drug_of_interest, " vs ", drug_reference) & outcome == k & analysis == "ow") %>% .$HR)
  LB <- (hrs %>% filter(contrast == paste0(drug_of_interest, " vs ", drug_reference) & outcome == k & analysis == "ow") %>% .$LB)
  UB <- (hrs %>% filter(contrast == paste0(drug_of_interest, " vs ", drug_reference) & outcome == k & analysis == "ow") %>% .$UB)
  
  
  for (q in numerical_covariates) {
    
    r <- var_labels[[q]]
    
    # optimal number of knots is 3 for all variables
    
    # print(paste0("Calculating optimal number of knots for variable ", r))
    # # Define the range of knots to test
    # k_range <- 3:5
    # 
    # # Initialize empty vectors to store results
    # bic_values <- numeric(length(k_range))
    # 
    # # Loop over each value of k
    # for (i in seq_along(k_range)) {
    #   k <- k_range[i]
    #   
    #   # Fit the model with k knots
    #   model <- cph(
    #     as.formula(paste0(
    #       "Surv(ckd_egfr50_censtime_yrs, ckd_egfr50_censvar) ~ studydrug2*rcs(", q, ",", k, ") + ",
    #       paste(setdiff(covariates, q), collapse=" + ") 
    #     )),
    #     data = cohort %>% filter(.imp == n.imp), x = TRUE, y = TRUE
    #   )
    #   
    #   # Store the BIC values
    #   bic_values[i] <- BIC(model)
    # }
    # 
    # # Find the optimal k based on minimum BIC
    # optimal_k_bic <- k_range[which.min(bic_values)]
    # 
    # print(paste0("Optimal number of knots for variable ", r, ": ", optimal_k_bic))
    
    optimal_k_bic = 3
    
    print(paste0("Spline plot for variable ", r))
    
    contrast_dataframe <- data.frame()
    
    for (i in 1:n.imp) {
      
      print(paste0("Imputation number ", i))
      
      # fit model with optimal number of knots
      final_model <- cph(
        as.formula(paste0(
          "Surv(", k, "_censtime_yrs, ", k, "_censvar) ~ studydrug2*rcs(", q, ",", optimal_k_bic, ") + ",
          paste(setdiff(covariates, q), collapse=" + "))
        ),
        data = cohort %>% filter(.imp == i), x = TRUE, y = TRUE
      )
      
      
      anova(final_model)
      p_value_non_linear <- anova(final_model)[2,3] # p value for non-linear interaction term
      
      print(paste0("p-value for non-linear interaction term with ", r, ": ", p_value_non_linear))
      
      q_vals <- seq(
        quantile(cohort[[q]], 0.025, na.rm = TRUE),
        quantile(cohort[[q]], 0.975, na.rm = TRUE),
        by = 0.05
      )
      
      contrast_spline <- contrast(
        final_model, 
        setNames(list("SGLT2 + GLP1", q_vals), c("studydrug2", q)),
        setNames(list("SGLT2 + DPP4/SU",  q_vals), c("studydrug2", q))
      )
      # extract beta and SE
      contrast_dataframe <- rbind(contrast_dataframe, as.data.frame(contrast_spline[c(q,'Contrast','SE')]) %>% mutate(.imp = i))
      
    }
    
    # pool results
    contrast_spline_df <- contrast_dataframe %>% 
      group_by(!!sym(q)) %>% 
      summarise(mean.coef = mean(Contrast),
                W = mean(SE^2),
                B = var(Contrast),
                T.var = W + (1+1/n.imp)*B,
                se.coef = sqrt(T.var),
                Upper = mean.coef + se.coef*1.96,
                Lower = mean.coef - se.coef*1.96) %>%
      select(-c(W, B, T.var, se.coef)) %>%
      rename(Contrast = mean.coef)
    
    
    # plot
    setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/output/")
    
    # define scale
    x_scale <- if (q == "uacr") {
      scale_x_continuous(trans = "log1p", breaks = c(0, 1, 3, 10, 30))
    } else {
      scale_x_continuous(breaks = function(x) sort(unique(c(0, pretty(x)))))
    }
    
    # set binwidth for histogram
    range_q <- range(if (q == "uacr") {log10(contrast_spline_df[[q]])} else (contrast_spline_df[[q]]), na.rm = TRUE)
    binwidth <- (range_q[2] - range_q[1]) / 30
    
    p_spline <- ggplot(data=contrast_spline_df, aes(x=.data[[q]], y=exp(Contrast))) +
      geom_line(data=contrast_spline_df,aes(x=.data[[q]], y=exp(Contrast)), size=1) +
      xlab(r) +
      ylab("Hazard ratio") +
      x_scale +
      scale_y_log10(breaks = c(0.25, 0.50, 0.75, 1.0, 1.50, 2.0)) +
      geom_ribbon(data=contrast_spline_df, aes(x=.data[[q]], ymin=exp(Lower), ymax=exp(Upper)), alpha=0.2) +
      geom_hline(yintercept = 1, linetype = "dashed")  +
      geom_hline(aes(yintercept = HR, linetype = "hr", size="hr"), color="#56B4E9")  +
      geom_hline(aes(yintercept = LB, linetype = "hr_95", size="hr_95"), color="#56B4E9")  +
      geom_hline(aes(yintercept = UB, linetype = "hr_95", size="hr_95"), color="#56B4E9")  +
      annotate("text", x = mean(range(contrast_spline_df[[q]])), y = 0.35, 
               label = "Favours SGLT2 + GLP1", 
               size = 5, hjust = 0.5, parse = F) +
      annotate("text", x = mean(range(contrast_spline_df[[q]])), y = 1.5, 
               label = "Favours SGLT2 + DPP4/SU", 
               size = 5, hjust = 0.5, parse = F) +
      theme_bw() +
      theme(text = element_text(size = 18),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            legend.position="top",
            legend.title = element_text(size=14, face = "italic"),
            legend.text = element_text(face="italic"),
            # Add custom axis lines
            axis.line = element_line(color = "black", size = 0.5), # General axis line style
            
            # Remove top and right axes lines
            axis.line.x.top = element_blank(),    # No line on the top
            axis.line.y.right = element_blank(),) +
      scale_linetype_manual(values = c(hr = "twodash", hr_95 = "twodash"), labels = c(hr = sprintf("%.2f", HR), hr_95 = paste0("95% CI ", sprintf("%.2f", LB), "-", sprintf("%.2f", UB))), name="Overall hazard ratio") +
      scale_size_manual(values = c(hr = 1, hr_95 = 0.5), labels = c(hr = sprintf("%.2f", HR), hr_95 = paste0("95% CI ", sprintf("%.2f", LB), "-", sprintf("%.2f", UB))), name="Overall hazard ratio") +
      coord_cartesian(xlim = c(min(contrast_spline_df[[q]]), max(contrast_spline_df[[q]])), ylim = c(0.25, 2.0), expand = F)
    
    
    histogram <- ggplot(cohort %>% filter(.imp == n.imp), aes(x = .data[[q]])) +
      geom_histogram(binwidth = binwidth, fill = "grey70", color = "black") +
      xlab(NULL) +
      ylab(NULL) +
      theme_bw() +
      theme(
        text = element_text(size = 16),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(color = "black", size = 0.5),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        plot.margin = margin(0, 5, 5, 5)
      ) +
      scale_y_continuous(breaks = function(x) {
        top <- (ceiling(max(x)/ 1000) * 1000)
        mid <- ceiling(top / 2)
        c(0, top)
      }) +
      x_scale +
      coord_cartesian(xlim = c(min(contrast_spline_df[[q]]), max(contrast_spline_df[[q]])), expand = FALSE)
    
    spline_and_histogram <- p_spline / histogram + plot_layout(heights = c(6, 1)) 
    
    tiff(paste0(today, "_HR_for_", k, "_by_", q, ".tiff"), width=10, height=4, units = "in", res=600) 
    print(spline_and_histogram)
    dev.off()
    
  }
}


############################9 CUMULATIVE INCIDENCE CURVES################################################################


# fit overlap-weighted model
fit.ow <- survfit(
  Surv(ckd_egfr40_censtime_yrs, ckd_egfr40_censvar) ~ studydrug2,
  data = cohort %>% filter(.imp == n.imp),
  weights = cohort %>% filter(.imp == n.imp) %>% .$overlap2
)

# fit unweighted model to get raw numbers at risk later down
fit.unadj <- survfit(
  Surv(ckd_egfr40_censtime_yrs, ckd_egfr40_censvar) ~ studydrug2,
  data = cohort %>% filter(.imp == n.imp)
)

cols_fig = c("SGLT2 + DPP4/SU" = "#CC79A7", "SGLT2 + GLP1" = "#56B4E9")

# create cumulative risk table of unweighted data - solely to extract risk table later down
cif_unadj <- ggsurvplot(
  fit = fit.unadj,
  fun = "event",
  data = cohort[cohort$.imp == n.imp,],
  palette = unname(cols_fig),
  color = "studydrug2",
  conf.int = T,
  legend.title = "",
  legend.labs = c("SGLT2 + DPP4/SU", "SGLT2 + GLP1"),
  font.legend = 12,
  font.title = "bold",
  font.subtitle = 12,
  font.x = 12,
  font.y = 12,
  
  legend = c(0.2, 0.8),
  break.y.by = 0.01,
  ylim = c(0, 0.05),
  censor.size = 1,
  surv.scale = "percent",
  risk.table = T,
  cumevents = T,
  tables.height = 0.2,
  fontsize = 4,
  tables.y.text = F,
  tables.y.text.col = T,
  tables.theme = theme_cleantable(),
  xlab = "Follow-up time (years)",
  ylab = "Cumulative incidence",
  title = "Cumulative incidence curves of composite kidney outcome",
  subtitle = "(≥40% eGFR decline/ESKD)"
)

# cumulative incidence curve of weighted data
cif_ow <- ggsurvplot(
  fit = fit.ow,
  fun = "event",
  data = cohort[cohort$.imp == n.imp,],
  palette = unname(cols_fig),
  color = "studydrug2",
  conf.int = T,
  legend.title = "",
  legend.labs = c("SGLT2 + DPP4/SU", "SGLT2 + GLP1"),
  font.legend = 12,
  font.title = "bold",
  font.subtitle = 12,
  font.x = 12,
  font.y = 12,
  
  legend = c(0.2, 0.8),
  break.y.by = 0.01,
  ylim = c(0, 0.05),
  censor.size = 1,
  surv.scale = "percent",
  risk.table = T,
  cumevents = T,
  tables.height = 0.15,
  fontsize = 4,
  tables.y.text = F,
  tables.y.text.col = T,
  tables.theme = theme_cleantable(),
  xlab = "Follow-up time (years)",
  ylab = "Cumulative incidence",
  title = "Overlap-weighted cumulative incidence curves by treatment",
  subtitle = "Composite of ≥40% eGFR decline or ESKD"
)

# add unweighted risk table to plot
cif_ow$plot <- cif_ow$plot + coord_cartesian(xlim = c(0.12, 2.9), ylim = c(0,0.04))
cif_ow$table <- cif_unadj$table + theme(plot.title = element_text(size = 12))
cif_ow$cumevents <- cif_unadj$cumevents + theme(plot.title = element_text(size = 12))

cif_ow
setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Output/")
tiff(paste0(today, "_", k, "_cumulative_incidence_curves.tiff"), width=7.5, height=6, units = "in", res=800)
print(cif_ow)
dev.off()

