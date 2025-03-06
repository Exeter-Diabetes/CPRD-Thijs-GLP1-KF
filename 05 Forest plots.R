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
      mutate(HR = ifelse(is.na(string), 1, HR),
             string = ifelse(is.na(string), "1.00 (ref.)", as.character(string)),
             analysis_label = case_when(
               analysis == "ow" ~ "Overlap-weighted",
               analysis == "iptw" ~ "IPTW",
               analysis == "adj" ~ "Multivariable adjustment only"
             )
      )
    
    temp_hrs[c("outcome", "contrast", "variable", "analysis")] <- lapply(temp_hrs[c("outcome", "contrast", "variable", "analysis")], factor)
  
    
    temp_hrs$drug <- as.factor(sub(" vs.*", "", temp_hrs$contrast))
    
    # Prepare labels
    labels_plot <- temp_hrs
    
    # set factor reference and order of levels
    
    if (m==1) {
      temp_hrs <- temp_hrs %>% mutate(
        drug = factor(drug, levels = c("SU", "DPP4", "SGLT2", "GLP1"))
      )
    }
    
    if (m==2) {  
      temp_hrs <- temp_hrs %>% mutate(
        drug = factor(drug, levels = c("DPP4/SU", "GLP1", "SGLT2", "GLP1/SGLT2"))
      )
    }
    
    if (m==3) {  
      temp_hrs <- temp_hrs %>% mutate(
        drug = factor(drug, levels = c("DPP4/SU", "GLP1", "SGLT2", "Oral semaglutide", "Subcutaneous semaglutide", "Other GLP1"))
      )
    }
    
    if (m==4) {
      temp_hrs <- temp_hrs %>% mutate(
        drug = factor(drug, levels = c("SGLT2", "DPP4/SU", "GLP1", "GLP1/SGLT2"))
      )
    }
    
    labels <- data.frame(matrix("", nrow = 1, ncol = length(labels_plot)))
    names(labels) <- names(labels_plot)
    labels <- labels %>% mutate(string = "Hazard Ratio (95% CI)",
                                nN = "Events/subjects",
                                drug = "Treatment arm"
                                )
    
    for (n in levels(labels_plot$drug)) {

      j = paste0(n, " vs ", levels(temp_hrs$drug)[1])
      
      labels_temp <- labels
      labels_temp$analysis_label <- paste0(" ", n)
      labels_temp$contrast <- j
      labels_plot <- rbind(labels_temp, labels_plot)
    }
    
    # have to coerce HR and CI to class numeric as they sometimes default to character
    
    class(labels_plot$HR) <- class(labels_plot$LB) <- class(labels_plot$UB) <- "numeric"
    
    
     
    plot_expression <- ""
    
    for (n in rev(levels(temp_hrs$drug))) {
      
      j = paste0(n, " vs ", levels(temp_hrs$drug)[1])
      
      labels_temp <- labels_plot %>% filter(contrast == j)
      
      if (n == levels(temp_hrs$drug)[1]) {
        labels_temp <- labels_temp %>% head(2)
      }

      p_counts <- labels_temp %>%
        ggplot(aes(y = factor(analysis_label, levels = rev(unique(analysis_label))))) +
        geom_text(aes(x = 1, label = nN), hjust = 1,
                  colour = ifelse(labels_temp$nN == labels$nN,
                                  "white", "black"),
                  fontface = ifelse(labels_temp$nN == labels$nN,
                                    "bold", "plain")) +
        annotate("text", x = 0.5,
                 y = length(unique(labels_temp$analysis_label)) + 1,
                 label = ifelse(j==unique(temp_hrs$contrast)[1], "Events/subjects", ""),
                 fontface = "bold") +
        theme_void() +
        coord_cartesian(xlim = c(-2, 3), ylim=c(1,ifelse(n==levels(temp_hrs$drug)[1], length(unique(labels_temp$analysis_label))+1, length(unique(labels_temp$analysis_label)))))

      p_hr <- labels_temp %>%
        ggplot(aes(y = factor(analysis_label, levels = rev(unique(analysis_label))))) +
        scale_x_continuous(trans = "log10", breaks = c(0.25, 0.5, 0.75, 1.0, 1.5, 2.0)) +
        coord_cartesian(ylim=c(1,ifelse(n==levels(temp_hrs$drug)[1], length(unique(labels_temp$analysis_label))+1, length(unique(labels_temp$analysis_label)))),
                        xlim=c(0.25, 2.0)) +
        theme_classic() +
        geom_point(aes(x=HR), shape=15, size=3) +
        geom_linerange(aes(xmin=LB, xmax=UB)) +
        geom_vline(xintercept = 1, linetype="dashed") +
        annotate("text", x = .65,
                 y = length(unique(labels_temp$analysis_label)) + 1,
                 label = ifelse(n==levels(temp_hrs$drug)[1], "Favours treatment arm", "")) +
        annotate("text", x = 1.5,
                 y = length(unique(labels_temp$analysis_label)) + 1,
                 label = ifelse(n==levels(temp_hrs$drug)[1], paste0("Favours ", temp_hrs$drug[1], collapse = ""), "")) +
        labs(x=ifelse(n==levels(temp_hrs$drug)[nlevels(as.factor(temp_hrs$drug))], "Hazard ratio", ""), y="") +
        theme(axis.line.y = element_blank(),
              axis.ticks.y= element_blank(),
              axis.text.y= element_blank(),
              axis.title.y= element_blank(),
              axis.line.x = if (n!=levels(temp_hrs$drug)[nlevels(temp_hrs$drug)]) {element_blank()},
              axis.text.x = if (n!=levels(temp_hrs$drug)[nlevels(temp_hrs$drug)]) {element_blank()},
              axis.ticks.x = if (n!=levels(temp_hrs$drug)[nlevels(temp_hrs$drug)]) {element_blank()},
              plot.title = element_text(hjust = 0.5),
              plot.subtitle = element_text(hjust = 0.5))

      p_left <- labels_temp %>%
        ggplot(aes(y = factor(analysis_label, levels = rev(unique(analysis_label))))) +
        geom_text(
          aes(x = 1, label = analysis_label),
          hjust = 0,
          fontface = ifelse(labels_temp$
                              analysis_label %in% paste0(" ", labels_temp$drug), "bold", "plain"),
          colour = ifelse(n==levels(temp_hrs$drug)[1] & !labels_temp$analysis_label %in% paste0(" ", labels_temp$drug), "white", "black")
        ) +
        theme_void() +
        coord_cartesian(xlim = c(0, 4), ylim=c(1,ifelse(n==levels(temp_hrs$drug)[1], length(unique(labels_temp$analysis_label))+1, length(unique(labels_temp$analysis_label)))))

      p_right <- labels_temp %>%
        ggplot(aes(y = factor(analysis_label, levels = rev(unique(analysis_label))))) +
        geom_text(
          aes(x = 0, label = string),
          hjust = 0,
          fontface = ifelse(labels_temp$string == "Hazard Ratio (95% CI)", "bold", "plain"),
          colour = ifelse(labels_temp$string == "Hazard Ratio (95% CI)",
                          "white", "black")) +
        annotate("text", x = 0.75,
                 y = length(unique(labels_temp$analysis_label)) + 1,
                 label = ifelse(n==levels(temp_hrs$drug)[1], "Hazard Ratio (95% CI)", ""),
                 fontface = "bold") +
        theme_void() +
        coord_cartesian(xlim = c(0, 4), ylim=c(1,ifelse(n==levels(temp_hrs$drug)[1], length(unique(labels_temp$analysis_label))+1, length(unique(labels_temp$analysis_label)))))

      assign(paste0("p_counts_", n), p_counts)
      assign(paste0("p_hr_", n), p_hr)
      assign(paste0("p_left_", n), p_left)
      assign(paste0("p_right_", n), p_right)

      plot_expression <- paste0(plot_expression, "`p_counts_", n, "` + `p_left_", n, "` + `p_hr_", n, "` + `p_right_", n, "` + ")
    }

    n.plots <- nlevels(as.factor(temp_hrs$contrast))

    # layout for plots below

    i <- 1
    
    height_all_plots <- 24
    height_first_plot <- 16
    
    layout <- paste("area(t = ",(i-1)*height_all_plots, ", l = 7, b = ",(i-1)*height_all_plots+height_first_plot,", r = 13), 
                    area(t = ",(i-1)*height_all_plots, ", l = 0, b = ",(i-1)*height_all_plots+height_first_plot,", r = 7), 
                    area(t = ",(i-1)*height_all_plots, ", l = 12, b = ",(i-1)*height_all_plots+height_first_plot,", r = 18), 
                    area(t = ",(i-1)*height_all_plots, ", l = 19, b = ", (i-1)*height_all_plots+height_first_plot,", r = 24)")


    for (i in 2:n.plots) {
      layout <- paste("area(t = ",(i-1)*height_all_plots-(24-height_first_plot), ", l = 7, b = ",(i-1)*height_all_plots+(height_first_plot),", r = 13), 
                      area(t = ",(i-1)*height_all_plots-(24-height_first_plot), ", l = 0, b = ",(i-1)*height_all_plots+(height_first_plot),", r = 7), 
                      area(t = ",(i-1)*height_all_plots-(24-height_first_plot), ", l = 12, b = ",(i-1)*height_all_plots+(height_first_plot),", r = 18), 
                      area(t = ",(i-1)*height_all_plots-(24-height_first_plot), ", l = 19, b = ", (i-1)*height_all_plots+(height_first_plot),", r = 24), ", layout)
    }

    layout <- paste0("c(", layout, ")")

    layout <- eval(str2lang(layout))

    # Final plot arrangement
    
    final_plot_expression <- paste0(plot_expression, "plot_layout(design = layout)")
    
    plot_for_saving <- eval(str2lang(final_plot_expression))
    
    setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Output/")
    tiff(paste0(today, "_HR_", m, "_", k, ".tiff"), width=18, height=5.5, units = "in", res=800) 
        print(plot_for_saving)
    dev.off()
  }
  
}


