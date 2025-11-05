########################0 SETUP####################################################################
setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Scripts/CPRD-Thijs-GLP1-KF/")
source("00 Setup.R")

setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Processed data/")

for (q in c("single_episodes", "no_insulin", "ncurrtx_upto_4", "all_3_exclusions", "hes_data_only")) {
  load(paste0(today, "_hrs_sens_", q, ".Rda"))
  assign(paste0("hrs_", q), hrs)
  rm(hrs)
}

load(paste0(today, "_hrs_incl_missing_uacr.Rda"))
load(paste0(today, "_hrs.Rda"))
load(paste0(today, "_hrs_5y.Rda"))
load(paste0(today, "_hrs_pp.Rda"))
load(paste0(today, "_hrs_fg_ow.Rda"))
load(paste0(today, "_factor_hrs.Rda"))

# set default studydrug variable
studydrug_var = paste0("studydrug", main)

load(paste0(today, "_t2d_glp1_imputed_data_withweights_studydrug", main, ".Rda"))

############################1 CUMULATIVE INCIDENCE CURVES################################################################

## primary outcome
k = "ckd_egfr40"
censvar_var=paste0(k, "_censvar")
censtime_var=paste0(k, "_censtime_yrs")
studydrug_var = paste0("studydrug", main)
weights_overlap = paste0("overlap", main)

# fit overlap-weighted model
fit.ow <- survfit(
  formula = reformulate(studydrug_var, 
                        response = paste0("Surv(", censtime_var, ", ", censvar_var, ")")),
  data = cohort[cohort$.imp == n.imp,],
  weights = cohort[cohort$.imp == n.imp,][[weights_overlap]]
)

# fit unweighted model to get raw numbers at risk later down
fit.unadj <- survfit(
  formula = reformulate(studydrug_var, 
                        response = paste0("Surv(", censtime_var, ", ", censvar_var, ")")),
  data = cohort[cohort$.imp == n.imp,]
)

cols_fig = c("SGLT2i + DPP4i/SU" = "#CC79A7", "SGLT2i + GLP1-RA" = "#56B4E9")

# create cumulative risk table of unweighted data - solely to extract risk table later down
cif_unadj <- ggsurvplot(
  fit = fit.unadj,
  fun = "event",
  data = cohort[cohort$.imp == n.imp,],
  palette = unname(cols_fig),
  color = "studydrug2",
  conf.int = T,
  legend.title = "",
  legend.labs = c("SGLT2i + DPP4i/SU", "SGLT2i + GLP1-RA"),
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
  legend.labs = c("SGLT2i + DPP4i/SU", "SGLT2i + GLP1-RA"),
  font.legend = 12,
  font.title = "bold",
  font.subtitle = 12,
  font.x = 12,
  font.y = 12,
  
  legend = c(0.8, 0.15),
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
cif_ow$plot <- cif_ow$plot + coord_cartesian(xlim = c(0.12, 2.9), ylim = c(0,0.03))
cif_ow$table <- cif_unadj$table + theme(plot.title = element_text(size = 12))
cif_ow$cumevents <- cif_unadj$cumevents + theme(plot.title = element_text(size = 12))

HR <- hrs %>% filter(outcome == k, variable == studydrug_var, analysis == "ow", contrast == "SGLT2i + GLP1-RA vs SGLT2i + DPP4i/SU") %>% .$HR  
LB <- hrs %>% filter(outcome == k, variable == studydrug_var, analysis == "ow", contrast == "SGLT2i + GLP1-RA vs SGLT2i + DPP4i/SU") %>% .$LB 
UB <- hrs %>% filter(outcome == k, variable == studydrug_var, analysis == "ow", contrast == "SGLT2i + GLP1-RA vs SGLT2i + DPP4i/SU") %>% .$UB 

p_value <- 2 * (1 - pnorm(abs(log(HR) / ((log(UB) - log(LB)) / (2 * 1.96)))))

cif_text <- data.frame(x = 0.25, y = 0.02, lab = paste0("HR ", sprintf("%.2f", HR), " (95% CI ", sprintf("%.2f", LB), ", ", sprintf("%.2f", UB), ")\np=", sprintf("%.2f", p_value)))

# add hazard ratio + p value as text
cif_ow$plot$layers[[length(cif_ow$plot$layers) + 1]] <- layer(geom="text", position = "identity", stat = "identity", 
                                                          mapping = aes(x = x, y = y, label = lab, hjust = 0), data = cif_text,
                                                          params = list(size = 4.5, colour = "black"))

# remove white square behind legend
cif_ow$plot <- cif_ow$plot + theme(
  legend.background = element_rect(fill = "transparent", colour = NA),
  legend.box.background = element_rect(fill = "transparent", colour = NA)
)

setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Output/")
tiff(paste0(today, "_cumulative_incidence_curves_", k, ".tiff"), width=7.5, height=6, units = "in", res=800)
print(cif_ow)
dev.off()


# get E value for main outcome
library(EValue)
evalues.HR(HR, LB, UB, rare = T)

############################2 FOREST PLOT FOR SECONDARY OUTCOMES################################################################

studydrug_var = paste0("studydrug", main)

hrs <- hrs %>% filter(!grepl("interaction", variable))

n.studydrug.vars <- hrs %>% .$variable %>% as.factor() %>% nlevels()

drug_reference = levels(cohort[[studydrug_var]])[1]
drug_of_interest = levels(cohort[[studydrug_var]])[2]

# Collect datasets across outcomes
# subset to outcome k and the current studydrug variable
temp <- hrs %>%
  filter(analysis == "ow",
         outcome != "death" & outcome != "ckd_egfr40",
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
  mutate(drug = factor(drug, levels = c("SGLT2i + DPP4i/SU", "SGLT2i + GLP1-RA")))


#  extract the reference group's nN for this outcome (first factor level) 
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
  "mace"                   = "MACE (incl. CV death)",
  "hf"                     = "Hospitalisation for heart failure",
  "acutepancreatitis"      = "Acute pancreatitis",
  "retinopathy"            = "Incident diabetic retinopathy",
  "lowerlimbfracture"      = "Lower limb fracture"
)

xmin = 0.10
xmax = 2.5

custom_order <- c("Outcome", "mace", "hf", "acutepancreatitis", "retinopathy", "lowerlimbfracture")

# Convert outcome to a factor with levels in the desired order
plot_df$outcome <- factor(plot_df$outcome, levels = custom_order)
plot_df <- plot_df %>%
  arrange(outcome, desc(HR %>% is.na()))

# Set y order (reverse for top-down order)
plot_df$y_order <- rev(seq_len(nrow(plot_df)))

# Apply mapping before building plot_df
plot_df <- plot_df %>%
  mutate(outcome = recode(outcome, !!!outcome_labels),
         # also update headings if we used outcome as heading text
         drug = ifelse(drug %in% outcomes, outcome, drug))

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
  geom_text(aes(x = 0.185, label = treated_nN,
                fontface = ifelse(is.na(HR), "bold", "plain")),
            hjust = 0) +
  
  # Reference n/N
  geom_text(aes(x = 0.325, label = ref_nN,
                fontface = ifelse(is.na(HR), "bold", "plain")),
            hjust = 0) +
  
  
  # HR text on right-hand side
  geom_text(aes(x = 2.1, label = string,
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
    axis.ticks  = element_blank(),   
    axis.text.y = element_blank(),
    axis.title.x = element_text(
      hjust = 0.72
    )
  ) +
  
  # custom axis line
  geom_segment(aes(x = 0.49, xend = 2.04, y = 0, yend = 0),
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
tiff(paste0(today, "_HR_", main, "_secondary_outcomes.tiff"), width=12, height=length(outcomes_per_drugclass)*1.2, units = "in", res=800)
print(forest_plot)
dev.off()



############################3 FOREST PLOT BY FACTORS################################################################

drug_reference = levels(cohort[[studydrug_var]])[1]
drug_of_interest = levels(cohort[[studydrug_var]])[2]

k = "ckd_egfr40"

factors <- c("malesex", "white_ethnicity", "low_egfr", "albuminuria", "predrug_cvd", "predrug_heartfailure")

factor_hrs$condition_label <- case_when(
  factor_hrs$factor == "malesex" & factor_hrs$condition == TRUE ~ "Male",
  factor_hrs$factor == "malesex" & factor_hrs$condition == FALSE ~ "Female",
  
  factor_hrs$factor == "white_ethnicity" & factor_hrs$condition == TRUE ~ "White",
  factor_hrs$factor == "white_ethnicity" & factor_hrs$condition == FALSE ~ "Non-white",
  
  factor_hrs$factor == "low_egfr" & factor_hrs$condition == TRUE ~ "<60 mL/min per 1.73m²",
  factor_hrs$factor == "low_egfr" & factor_hrs$condition == FALSE ~ "≥60 mL/min per 1.73m²",
  
  factor_hrs$factor == "albuminuria" & factor_hrs$condition == TRUE ~ "<3 mg/mmol",
  factor_hrs$factor == "albuminuria" & factor_hrs$condition == FALSE ~ "≥3 mg/mmol",
  
  factor_hrs$factor == "predrug_cvd" & factor_hrs$condition == TRUE ~ "Present",
  factor_hrs$factor == "predrug_cvd" & factor_hrs$condition == FALSE ~ "Absent",
  
  factor_hrs$factor == "predrug_heartfailure" & factor_hrs$condition == TRUE ~ "Present",
  factor_hrs$factor == "predrug_heartfailure" & factor_hrs$condition == FALSE ~ "Absent",
  
  TRUE ~ NA_character_
) 

# Filter to the CKD outcome only
temp <- factor_hrs %>%
  filter(analysis == "ow",
         outcome == k,
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
                               "low_egfr" = "eGFR",
                               "albuminuria" = "uACR",
                               "predrug_cvd" = "Atherosclerotic CVD",
                               "predrug_heartfailure" = "Heart failure",
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
      sprintf("%.2f", round(factor_hrs$p_value_interaction[match(factor_label, recode(
        factor_hrs$factor,
        "malesex" = "Sex",
        "white_ethnicity" = "Ethnicity",
        "low_egfr" = "eGFR",
        "albuminuria" = "uACR",
        "predrug_cvd" = "Atherosclerotic CVD",
        "predrug_heartfailure" = "Heart failure",
      ))], 2)), 
      ")"
    )
  )

# Bind together: header row first, then headings + subgroup rows 
plot_df <- bind_rows(header, heading_rows, plot_df)

# custom_order <- c(" Header", "Age", "Sex", "Ethnicity", "Atherosclerotic CVD", "Heart failure")
custom_order <- c(" Header", "Sex", "Ethnicity", "eGFR", "uACR", "Atherosclerotic CVD", "Heart failure")

# Convert factor_label to a factor with levels in the desired order
plot_df$factor_label <- factor(plot_df$factor_label, levels = custom_order)
plot_df <- plot_df %>%
  arrange(factor_label, desc(HR %>% is.na()))


# Arrange by factor_label (groups) and HR NA (header first) 
plot_df <- plot_df %>%
  arrange(factor_label, desc(is.na(HR)), condition_label)

# Reassign y_order
plot_df$y_order <- rev(seq_len(nrow(plot_df)))

# plot
forest_plot <- ggplot(plot_df, aes(y = y_order)) +
  geom_linerange(aes(xmin = LB, xmax = UB), na.rm = TRUE) +
  geom_point(aes(x = HR), shape = 15, size = 3, na.rm = TRUE) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  
  # Left labels (factor headings + subgroup names)
  geom_text(aes(x = 0.09, label = condition_label,
                fontface = ifelse(is.na(HR) , "bold", "plain")),
            hjust = 0) +
  
  # n/N columns (only populated for header and subgroup rows)
  geom_text(aes(x = 0.16, label = treated_nN,
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
                     limits = c(0.09, 2.7)) +
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
tiff(paste0(today, "_HR_", main, "_main_outcome_by_factors.tiff"), width=14, height=5, units = "in", res=800)
print(forest_plot)
dev.off()

############################4 FOREST PLOT FOR DPP4i vs SU################################################################

m = 1
studydrug_var = paste0("studydrug", m)

drug_reference = levels(cohort[[studydrug_var]])[1]
drug_of_interest = levels(cohort[[studydrug_var]])[2]

# Collect datasets across outcomes
# subset to outcome k and the current studydrug variable
temp <- hrs %>%
  filter(analysis == "ow",
         outcome %in% (outcomes %>% setdiff(c("death"))),
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
  filter(!grepl(drug, "GLP1-RA")) %>%
  mutate(drug = factor(drug, levels = c("SGLT2i + SU", "SGLT2i + DPP4i")))


#  extract the reference group's nN for this outcome (first factor level) 
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

custom_order <- c("Outcome", "ckd_egfr40", "mace", "hf", "acutepancreatitis", "retinopathy", "lowerlimbfracture")

plot_df$outcome <- factor(plot_df$outcome, levels = custom_order)
plot_df <- plot_df %>%
  arrange(outcome, desc(HR %>% is.na()))

# Set y order (reverse for top-down order)
plot_df$y_order <- rev(seq_len(nrow(plot_df)))

# Define pretty labels
outcome_labels <- c(
  "ckd_egfr40"             = "≥40% eGFR decline/ESKD",
  "mace"                   = "MACE (incl. CV death)",
  "hf"                     = "Hospitalisation for HF",
  "acutepancreatitis"      = "Acute pancreatitis",
  "retinopathy"            = "Incident diabetic retinopathy",
  "lowerlimbfracture"      = "Lower limb fracture"
)

xmin = 0.10
xmax = 3.2
# Apply mapping before building plot_df
plot_df <- plot_df %>%
  mutate(outcome = recode(outcome, !!!outcome_labels),
         # also update headings if we used outcome as heading text
         drug = ifelse(drug %in% outcomes, outcome, drug))


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
  geom_text(aes(x = 0.185, label = treated_nN,
                fontface = ifelse(is.na(HR), "bold", "plain")),
            hjust = 0) +
  
  # Reference n/N
  geom_text(aes(x = 0.325, label = ref_nN,
                fontface = ifelse(is.na(HR), "bold", "plain")),
            hjust = 0) +
  
  
  # HR text on right-hand side
  geom_text(aes(x = 2.65, label = string,
                fontface = ifelse(HR %>% is.na(), "bold", "plain")),
            hjust = 0) +
  
  # text to indicate which drug to favour
  annotate("text", x = .75,
           y = max(plot_df$y_order), fontface = "italic",
           label = paste0("Favours ", drug_of_interest %>% substr(9, nchar(drug_of_interest)))) +
  
  annotate("text", x = 1.5,
           y = max(plot_df$y_order), fontface = "italic",
           label = paste0("Favours ", drug_reference %>% substr(9, nchar(drug_reference)))) +
  
  scale_x_continuous(trans = "log10",
                     breaks = c(0.5, 0.75, 1, 1.5, 2.5),
                     limits = c(xmin, xmax)) +
  scale_y_continuous(NULL, breaks = NULL) +
  labs(x = "HR (95% CI)", y = NULL) +
  theme_classic() +
  theme(
    panel.grid = element_blank(),
    axis.line   = element_blank(),
    axis.ticks  = element_blank(),   
    axis.text.y = element_blank(),
    axis.title.x = element_text(
      hjust = 0.72
    )
  ) +
  
  # custom axis line
  geom_segment(aes(x = 0.45, xend = 2.7, y = 0, yend = 0),
               inherit.aes = FALSE, linewidth = 0.4, color = "black") +
  
  # custom ticks
  geom_segment(
    data = data.frame(x = c(0.5, 0.75, 1, 1.5, 2.5)),
    aes(x = x, xend = x, y = 0, yend = -0.2),
    inherit.aes = FALSE,
    linewidth = 0.4,
    color = "black"
  )


setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Output/")
tiff(paste0(today, "_HR_", m, "_all_outcomes.tiff"), width=12, height=length(outcomes %>% setdiff("death"))*0.6, units = "in", res=800)
print(forest_plot)
dev.off()

############################5 FOREST PLOT BY GLP1-RA TYPE################################################################

for (m in 3) {
  studydrug_var = paste0("studydrug", m)
  
  drug_reference = "DPP4/SU"
  drug_of_interest = "GLP1-RA subgroup"
  
  # Collect datasets across outcomes
  # subset to outcome k and the current studydrug variable
  temp <- hrs %>%
    filter(analysis == "ow",
           outcome %in% (outcomes %>% setdiff(c("death"))),
           variable == paste0("studydrug", m)) %>%
    mutate(HR = ifelse(is.na(string), 1, HR),
           string = ifelse(is.na(string), "1.00 (ref.)", as.character(string))
    )
  
  # ensure these cols are factors as used previously
  temp[c("outcome","contrast","variable","analysis")] <-
    lapply(temp[c("outcome","contrast","variable","analysis")], factor)
  
  
  # drug label extracted from contrast (as before)
  temp$drug <- as.factor(sub(" vs.*", "", temp$contrast))
  
  if (m == 3) {
    temp_k <- temp %>%
      mutate(drug = factor(drug, levels = c("SGLT2i + DPP4i/SU",
                                            "GLP1-RA with direct kidney outcome evidence",
                                            "Other GLP1-RA")))
  } else {
    temp_k <- temp
  }
  
  #  extract the reference group's nN for this outcome (first factor level) 
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
  
  
  #  1. Create the global header row 
  header <- tibble(
    outcome = " ",
    drug = " ",
    HR = NA_real_, LB = NA_real_, UB = NA_real_,
    string = "HR (95% CI)",
    treated_nN = paste0(drug_of_interest, " (n/N)"),
    ref_nN = paste0(drug_reference, " (n/N)")
  )
  
  #  2. Create subgroup heading rows (no n/N titles here) 
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
    "ckd_egfr40"             = "≥40% eGFR decline/ESKD",
    "mace"                   = "MACE (incl. CV death)",
    "hf"                     = "Hospitalisation for HF",
    "acutepancreatitis"      = "Acute pancreatitis",
    "retinopathy"            = "Incident diabetic retinopathy",
    "lowerlimbfracture"      = "Lower limb fracture"
  )
  
  xmin = 0.02
  xmax = 2.8
  
  # Apply mapping before building plot_df
  plot_df <- plot_df %>%
    mutate(outcome = recode(outcome, !!!outcome_labels),
           # also update headings if we used outcome as heading text
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
    geom_text(aes(x = 0.085, label = treated_nN,
                  fontface = ifelse(is.na(HR), "bold", "plain")),
              hjust = 0) +
    
    # Reference n/N
    geom_text(aes(x = 0.19, label = ref_nN,
                  fontface = ifelse(is.na(HR), "bold", "plain")),
              hjust = 0) +
    
    
    # HR text on right-hand side
    geom_text(aes(x = 2.1, label = string,
                  fontface = ifelse(HR %>% is.na(), "bold", "plain")),
              hjust = 0) +
    
    # text to indicate which drug to favour
    annotate("text", x = .6,
             y = max(plot_df$y_order), fontface = "italic",
             label = paste0("Favours ", drug_of_interest)) +
    
    annotate("text", x = 1.4,
             y = max(plot_df$y_order), fontface = "italic",
             label = paste0("Favours ", drug_reference)) +
    
    scale_x_continuous(trans = "log10",
                       breaks = c(0.3, 0.5, 0.75, 1, 1.5, 2.0),
                       limits = c(xmin, xmax)) +
    scale_y_continuous(NULL, breaks = NULL) +
    labs(x = "HR (95% CI)", y = NULL) +
    theme_classic() +
    theme(
      panel.grid = element_blank(),
      axis.line   = element_blank(),
      axis.ticks  = element_blank(),   
      axis.text.y = element_blank(),
      axis.title.x = element_text(
        hjust = 0.72
      )
    ) +
    
    # custom axis line
    geom_segment(aes(x = 0.3, xend = 2.05, y = 0, yend = 0),
                 inherit.aes = FALSE, linewidth = 0.4, color = "black") +
    
    # custom ticks
    geom_segment(
      data = data.frame(x = c(0.3, 0.5, 0.75, 1, 1.5, 2.0)),
      aes(x = x, xend = x, y = 0, yend = -0.2),
      inherit.aes = FALSE,
      linewidth = 0.4,
      color = "black"
    )
  
  
  setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Output/")
  tiff(paste0(today, "_HR_", m, "_all_outcomes.tiff"), width=12, height=length(outcomes %>% setdiff("death"))*0.8, units = "in", res=800)
  print(forest_plot)
  dev.off()
}

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

for (q in c("single_episodes", "no_insulin", "ncurrtx_upto_4", "all_3_exclusions", "hes_data_only", "incl_missing_uacr")) {
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
  mutate(drug = factor(drug, levels = c("SGLT2i + DPP4i/SU", "SGLT2i + GLP1-RA")))


#  extract the reference group's nN for this outcome (first factor level) 
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

custom_order <- c("≥40% eGFR decline/ESKD", "all", "single_episodes", "no_insulin", "ncurrtx_upto_4", "all_3_exclusions", "hes_data_only", "incl_missing_uacr")

# Convert analysis to a factor with levels in the desired order
plot_df$analysis <- factor(plot_df$analysis, levels = custom_order)
plot_df <- plot_df %>%
  arrange(analysis, desc(HR %>% is.na()))

# Set y order (reverse for top-down order)
plot_df$y_order <- rev(seq_len(nrow(plot_df)))

# Define pretty labels
analysis_labels <- c(
  "all"                    = "Primary study sample",
  "single_episodes"        = "Excluding GLP1-RA initiation if\nprevious DPP4/SU initiation",
  "no_insulin"             = "Excluding individuals treated\nwith insulin",
  "ncurrtx_upto_4"         = "Excluding individuals treated with\n>4 glucose-lowering treatments",
  "all_3_exclusions"       = "Above exclusions combined",
  "hes_data_only"          = "Excluding individuals without\nlinked hospital inpatient data",
  "incl_missing_uacr"      = "Including individuals with\nmissing baseline uACR"
)

xmin = 0.08
xmax = 2.5
# Apply mapping before building plot_df
plot_df <- plot_df %>%
  mutate(analysis = recode(analysis, !!!analysis_labels))


# Forest plot with tabulated extras
forest_plot <- ggplot(plot_df, aes(y = y_order)) +
  # Forest CI + point
  geom_linerange(aes(xmin = LB, xmax = UB), na.rm = TRUE) +
  geom_point(aes(x = HR), shape = 15, size = 3, na.rm = TRUE) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  
  # Outcome headings on the left
  geom_text(aes(x = xmin, label = analysis,
                fontface = ifelse(HR %>% is.na(), "bold", "plain")),
            hjust = 0) +
  
  # Events/subjects in middle-left
  # Treated n/N
  geom_text(aes(x = 0.175, label = treated_nN,
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
    axis.ticks  = element_blank(),   
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
tiff(paste0(today, "_HR_", main, "_sensitivity_analysis.tiff"), width=12, height=length(analysis_labels)*0.8, units = "in", res=800)
print(forest_plot)
dev.off()


############################7 FOREST PLOT FOR DIFFERENT ANALYSIS APPROACHES################################################################


# set default studydrug variable
studydrug_var = paste0("studydrug", main)

drug_reference = levels(cohort[[studydrug_var]])[1]
drug_of_interest = levels(cohort[[studydrug_var]])[2]

k = "ckd_egfr40"

temp <- hrs %>% rbind(hrs_fg_ow) %>% 
  rbind(hrs_5y %>% mutate(analysis = "ow_5y")) %>%
  rbind(hrs_pp %>% mutate(analysis = "ow_pp")) %>%
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
  mutate(drug = factor(drug, levels = c("SGLT2i + DPP4i/SU", "SGLT2i + GLP1-RA")))


# extract the reference group's nN for this outcome (first factor level)
ref_level <- levels(temp_k$drug)[1]

ref_vals <- temp_k %>%
  filter(drug == ref_level) %>%
  mutate(ref_nN = paste("    ", trimws(nN))) %>%
  select(analysis, ref_nN)

# Create treated/ref columns and join by analysis
temp_k <- temp_k %>%
  filter(drug != ref_level) %>%
  mutate(treated_nN = paste("    ", trimws(nN))) %>%
  left_join(ref_vals, by = "analysis")

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

custom_order <- c("≥40% eGFR decline/ESKD", "ow", "ow_5y", "ow_pp", "adj", "iptw", "fg_ow")

# Convert analysis to a factor with levels in the desired order
plot_df$analysis <- factor(plot_df$analysis, levels = custom_order)
plot_df <- plot_df %>%
  arrange(analysis, desc(HR %>% is.na()))

# Reassign y_order
plot_df$y_order <- rev(seq_len(nrow(plot_df)))

# Define pretty labels
analysis_labels <- c(
  "ow"          = "Overlap-weighting\n(primary analysis)",
  "ow_5y"       = "Follow-up extended up to 5 years",
  "ow_pp"       = "Per-protocol analysis",
  "iptw"        = "Inverse probability of\ntreatment weighting",
  "adj"         = "Multivariable adjustment only",
  "fg_ow"       = "Fine-Gray competing risk analysis\n(overlap-weighted)"
)

xmin = 0.08
xmax = 2.5
# Apply mapping before building plot_df
plot_df <- plot_df %>%
  mutate(analysis = recode(analysis, !!!analysis_labels))

# Forest plot with tabulated extras
forest_plot <- ggplot(plot_df, aes(y = y_order)) +
  # Forest CI + point
  geom_linerange(aes(xmin = LB, xmax = UB), na.rm = TRUE) +
  geom_point(aes(x = HR), shape = 15, size = 3, na.rm = TRUE) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  
  # Outcome headings on the left
  geom_text(aes(x = xmin, label = analysis,
                fontface = ifelse(HR %>% is.na(), "bold", "plain")),
            hjust = 0) +
  
  # Events/subjects in middle-left
  # Treated n/N
  geom_text(aes(x = 0.175, label = treated_nN,
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
    axis.ticks  = element_blank(),   
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
tiff(paste0(today, "_HR_", main, "_different_analysis_approaches.tiff"), width=12, height=length(analysis_labels)*0.8, units = "in", res=800)
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
options(contrasts=c("contr.treatment", "contr.treatment"))

numerical_covariates <- c("dstartdate_age", "prebmi", "preegfr", "uacr", "prehba1c", "ckdpc_40egfr_score")

var_labels <- c(
  dstartdate_age = "Age (years)",
  prehba1c = "HbA1c (mmol/mol)",
  prebmi = "BMI (kg/m²)",
  preegfr = "eGFR (mL/min per 1.73m²)",
  uacr = "uACR (mg/mmol)",
  ckdpc_40egfr_score = "CKD-PC risk score"
)


for (k in outcomes_per_drugclass[1]) {
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
    #       "Surv(ckd_egfr40_censtime_yrs, ckd_egfr40_censvar) ~ studydrug2*rcs(", q, ",", k, ") + ",
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
    
    wald <- rep(NA, n.imp)
    
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
      
      
      a <- anova(final_model)
      
      row_id <- grep("All interactions", rownames(a), ignore.case = TRUE)[1] # total interaction
      
      # row_id <- grep("Nonlinear Interaction ", rownames(a), ignore.case = TRUE) # nonlinear interaction component only
      
      
      wald[i] <- a[row_id, "Chi-Square"]
      df_interaction <- a[row_id, "d.f."]
      
      q_vals <- seq(
        quantile(cohort[[q]], 0.025, na.rm = TRUE),
        quantile(cohort[[q]], 0.975, na.rm = TRUE),
        by = 0.05
      )
      
      contrast_spline <- contrast(
        final_model, 
        setNames(list("SGLT2i + GLP1-RA", q_vals), c("studydrug2", q)),
        setNames(list("SGLT2i + DPP4i/SU",  q_vals), c("studydrug2", q))
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
    

    m <- length(wald)
    Qbar <- mean(wald)
    B <- var(wald)
    Tstat <- Qbar + (1 + 1/m) * B
    Fstat <- Tstat / df_interaction
    R <- (1 + 1/m) * B / Qbar
    nu <- (m - 1) * (1 + 1 / R)^2
    p_value_interaction <- 1 - pf(Fstat, df1 = df_interaction, df2 = nu)
    
    print(paste0("P for interaction with ", q, ": ", p_value_interaction))
    assign(paste0("p_value_interaction_", q), p_value_interaction)
    
    # plot
    setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/output/")
    
    # define scale
    x_scale <- if (q == "uacr") {
      scale_x_continuous(trans = "log1p", breaks = c(0, 1, 3, 10, 30))
    } else {
      scale_x_continuous(breaks = function(x) sort(unique(c(0, pretty(x)))))
    }
    
    # set midpoint of scale to anchor text
    if (q == "uacr") {
      x_mid <- (exp(mean(log1p(contrast_spline_df[[q]]))) - 1) / 2.5
    } else {
      x_mid <- mean(range(contrast_spline_df[[q]]))
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
      annotate("text", x = x_mid, y = 0.35, 
               label = "Favours SGLT2i + GLP1-RA", 
               size = 5, hjust = 0.5, parse = F) +
      annotate("text", x = x_mid, y = 1.5, 
               label = "Favours SGLT2i + DPP4i/SU", 
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


