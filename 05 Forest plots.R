########################0 SETUP####################################################################
setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Scripts/CPRD-Thijs-GLP1-KF/")
source("00 Setup.R")

setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Processed data/")
load(paste0(today, "_hrs.Rda"))

############################1 FOREST PLOT FOR HR BY DRUG CLASS################################################################

# GLP1-RA trial meta-analysis HR: 0.79, 95% CI 0.66-0.95 (Circulation. 2024 Nov 26;150(22):1781-1790.)

# other GLP1-RA trial meta-analysis HR: HR 0·81, 95% CI 0·72–0·92 (Lancet Diabetes Endocrinol. 2025 Jan;13(1):15-28.)

n.studydrug.vars <- hrs %>% .$variable %>% as.factor() %>% nlevels()

for (m in 1:n.studydrug.vars) {
  for (k in outcomes_per_drugclass) {
    
    temp_hrs <- hrs %>% 
      filter(analysis!="unadj") %>%
      filter(outcome == k) %>%
      filter(variable == paste0("studydrug", m)) %>%
      mutate(string = ifelse(is.na(string), "(ref.)", as.character(string)),
             analysis_label = case_when(
               analysis == "ow" ~ "Overlap-weighted",
               analysis == "iptw" ~ "IPTW",
               analysis == "adj" ~ "Multivariable adjustment only"
             )
      )
    
    temp_hrs$drug <- as.factor(sub(" vs.*", "", temp_hrs$contrast))
    
    # Prepare labels
    labels_plot <- temp_hrs
    
    
    labels <- data.frame(matrix("", nrow = 1, ncol = length(labels_plot)))
    names(labels) <- names(labels_plot)
    labels <- labels %>% mutate(analysis = "Analysis approach",
                                string = "Hazard Ratio (95% CI)",
                                nN = "Events/subjects",
                                drug = "Treatment arm"
                                # ref_nN = if (m == 1) {"(SU)"} else {"(DPP4/SU)"}
                                )
    
    # for (d in levels(temp_hrs$analysis)) { 
    #   
    #   labels_temp <- labels
    #   labels_temp$analysis <- d
    #   labels_plot <- rbind(labels_temp, labels_plot)
    # }
    
    labels_plot <- rbind(labels, labels_plot)

    # labels_plot <- labels_plot %>%
    #   mutate(
    #     contrast = ifelse(analysis == "Overall", paste0(" ", analysis_label), albuminuria_status)
    #   )
    
    # have to coerce HR and CI to class numeric as they sometimes default to character
    
    class(temp_hrs$HR) <- class(temp_hrs$LB) <- class(temp_hrs$UB) <- "numeric"
     
    plot_expression <- ""
    
    for (n in levels(temp_hrs$drug)) {

      p_counts <- labels_plot %>% filter(drug==n) %>%
        ggplot(aes(y = factor(analysis_label, levels = rev(unique(analysis_label))))) +
        geom_text(aes(x = 1, label = nN), hjust = 1,
                  colour = ifelse(!n==unique(temp_hrs$drug)[1] & labels_plot[labels_plot$drug == n,]$nN == labels$nN,
                                  "white", "black"),
                  fontface = ifelse(n==unique(temp_hrs$drug)[1] & labels_plot[labels_plot$drug == n,]$nN == labels$nN,
                                    "bold", "plain")) +
        theme_void() +
        coord_cartesian(xlim = c(-2, 3))

      p_hr <-
        temp_hrs %>%
        filter(drug==n) %>%
        ggplot(aes(y = factor(analysis_label, levels = rev(unique(analysis_label))))) +
        scale_x_continuous(trans = "log10", breaks = c(0.5, 0.75, 1.0, 1.5, 2.25, 3.25)) +
        coord_cartesian(ylim=c(1,length(unique(labels_plot[labels_plot$drug == n,]$analysis_label)) + 1),
                        xlim=c(0.5, 3.25)) +
        theme_classic() +
        geom_point(aes(x=HR), shape=15, size=3) +
        geom_linerange(aes(xmin=LB, xmax=UB)) +
        geom_vline(xintercept = 1, linetype="dashed") +
        annotate("text", x = .65,
                 y = length(unique(temp_hrs[temp_hrs$drug == n,]$analysis_label)) + 2,
                 label = ifelse(n==unique(temp_hrs$drug)[1], "Favours treatment arm", "")) +
        annotate("text", x = 1.5,
                 y = length(unique(temp_hrs[temp_hrs$drug == n,]$analysis_label)) + 2,
                 label = ifelse(n==unique(temp_hrs$drug)[1], paste0("Favours ", temp_hrs$drug[1], collapse = ""), "")) +
        labs(x=ifelse(n==unique(temp_hrs$drug)[nlevels(as.factor(temp_hrs$drug))], "Hazard ratio", ""), y="") +
        theme(axis.line.y = element_blank(),
              axis.ticks.y= element_blank(),
              axis.text.y= element_blank(),
              axis.title.y= element_blank(),
              axis.line.x = if (!n==unique(temp_hrs$drug)[nlevels(as.factor(temp_hrs$drug))]) {element_blank()},
              axis.text.x = if (!n==unique(temp_hrs$drug)[nlevels(as.factor(temp_hrs$drug))]) {element_blank()},
              axis.ticks.x = if (!n==unique(temp_hrs$drug)[nlevels(as.factor(temp_hrs$drug))]) {element_blank()},
              plot.title = element_text(hjust = 0.5),
              plot.subtitle = element_text(hjust = 0.5))

      p_left <-
        labels_plot %>%
        filter(drug==n) %>%
        ggplot(aes(y = rev(unique((analysis_label))))) +
        geom_text(
          aes(x = 1, label = analysis_label),
          hjust = 0,
          fontface = ifelse(labels_plot[labels_plot$drug == n,]$
                              analysis_label == paste0(" ", labels_plot[labels_plot$drug == n,]$analysis_label), "bold", "plain"),
          colour = ifelse(n==unique(temp_hrs$drug)[1] & labels_plot[labels_plot$drug == n,]$
                            analysis_label == "drug" | !labels_plot[labels_plot$drug == n,]$
                            analysis_label == "drug", "black", "white")
        ) +
        theme_void() +
        coord_cartesian(xlim = c(0, 4))

      p_right <-
        labels_plot %>%
        filter(drug==n) %>%
        ggplot(aes(y = factor(string, levels = rev(unique(string))))) +
        geom_text(
          aes(x = 0, label = string),
          hjust = 0,
          fontface = ifelse(labels_plot[labels_plot$drug == n,]$string == "Hazard Ratio (95% CI)", "bold", "plain"),
          colour = ifelse(labels_plot[labels_plot$drug == n,]$string == "Hazard Ratio (95% CI)" & !n==unique(temp_hrs$drug)[1],
                          "white", "black")) +
        theme_void() +
        coord_cartesian(xlim = c(0, 4))

      assign(paste0("p_counts_", n), p_counts)
      assign(paste0("p_hr_", n), p_hr)
      assign(paste0("p_left_", n), p_left)
      assign(paste0("p_right_", n), p_right)

      plot_expression <- paste0(plot_expression, "`p_counts_", n, "` + `p_left_", n, "` + `p_hr_", n, "` + `p_right_", n, "` + ")
    }

    n.plots <- nlevels(as.factor(temp_hrs$drug))

    # layout for plots below

    i <- 1
    layout <- paste("area(t = ",(i-1)*24, ", l = 7, b = ",(i-1)*24+24,", r = 13), area(t = ",(i-1)*24, ", l = 0, b = ",(i-1)*24+24,", r = 7), area(t = ",(i-1)*24, ", l = 12, b = ",(i-1)*24+24,", r = 18), area(t = ",(i-1)*24, ", l = 19, b = ", (i-1)*24+24,", r = 24)")


    for (i in 2:n.plots) {
      layout <- paste("area(t = ",(i-1)*24, ", l = 7, b = ",(i-1)*24+24,", r = 13), area(t = ",(i-1)*24, ", l = 0, b = ",(i-1)*24+24,", r = 7), area(t = ",(i-1)*24, ", l = 12, b = ",(i-1)*24+24,", r = 18), area(t = ",(i-1)*24, ", l = 19, b = ", (i-1)*24+24,", r = 24), ", layout)
    }

    layout <- paste0("c(", layout, ")")

    layout <- eval(str2lang(layout))

    # Final plot arrangement
    
    plot_expression <- paste0(plot_expression, "plot_layout(design = layout)")
    
    setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Output/")
    tiff(paste0(today, "_HR_", m, "_", k, ".tiff"), width=18, height=5.5, units = "in", res=800) 
    eval(str2lang(plot_expression))
    dev.off()
  }
}
