# Trapezoidal rule
trap.rule <- function(x,y) sum(diff(x) * (y[-1] + y[-length(y)])) / 2

# Process data
source("5a. Process Shiny data.R")

# Inputs one selects on the Shiny app
adc <- "T-DXd"
pk_x_axis <- "Cmax"
Outcomes <- c("ORR", "DLT")
cui_threshold <- 2.5
weights <- c(ORR = 0.81, DLT = 0.19)


# Shiny output
# ------------
# Loop through the scenarios
for (scenario in c("Scenario 1 (only phase I)", "Scenario 2 (phases I and II)", "Scenario 3 (all phases)")) {
  for (i in seq_along(c("ORR", "PFS", "DLT"))) {
    outcomei <- c("ORR", "PFS", "DLT")[i]
    if (outcomei %in% Outcomes) {
      PK_param <- pk_x_axis
      if (adc == "T-DM1") {
        er_mod_tb <- switch(scenario,
                            `Scenario 1 (only phase I)` = tdm1_models[[1]],
                            `Scenario 2 (phases I and II)` = tdm1_models[[2]],
                            `Scenario 3 (all phases)` = tdm1_models[[3]])
      } else {
        er_mod_tb <- switch(scenario,
                            `Scenario 1 (only phase I)` = tdxd_models[[1]],
                            `Scenario 2 (phases I and II)` = tdxd_models[[2]],
                            `Scenario 3 (all phases)` = tdxd_models[[3]])
      }
      assign(paste0(outcomei, "i"), 
             unlist(er_mod_tb %>% filter(Outcome == outcomei) %>% pull(PK_param))[1], envir = .GlobalEnv) # Intercept
      assign(paste0(outcomei, "s"), 
             unlist(er_mod_tb %>% filter(Outcome == outcomei) %>% pull(PK_param))[2], envir = .GlobalEnv) # Slope
    }
  }
  
  # Link with the PK data
  pk_tb_data <- if (adc == "T-DM1") tdm1_pk_data else tdxd_pk_data
  pk_tb <- pk_tb_data %>%
    group_by(id, OCC, AMT) %>%
    reframe(AUC = max(AUC),
            Cmax = max(Cc_wghtd),
            Cmin = Cc_wghtd[time == 504], .groups = "drop") %>%
    select(id:AMT, all_of(pk_x_axis))
  colnames(pk_tb)[colnames(pk_tb) == pk_x_axis] <- "x1"
  pk_tb$x2 <- pk_tb$x1
  
  if ("ORR" %in% Outcomes) pk_tb <- pk_tb %>% mutate(pORR = inv_logit_fn(ORRi + ORRs * x1)) # # probability of ORR
  if ("PFS" %in% Outcomes) pk_tb <- pk_tb %>% mutate(pPFS = inv_logit_fn(PFSi + PFSs * x1)) # # probability of ORR
  if ("DLT" %in% Outcomes) pk_tb <- pk_tb %>% mutate(pDLT = inv_logit_fn(DLTi + DLTs * x2))
  
  # Add weights
  if ("ORR" %in% Outcomes) pk_tb$weightedpORR <- pk_tb$pORR * weights[["ORR"]]
  if ("PFS" %in% Outcomes) pk_tb$weightedpPFS <- pk_tb$pPFS * weights[["PFS"]]
  if ("DLT" %in% Outcomes) pk_tb$weightedpDLT <- (1 - pk_tb$pDLT) * weights[["DLT"]]
  
  # Compute clinical utility index based on weights
  pk_tb$pCUI <- pk_tb %>%
    ungroup() %>%
    select(contains("weighted")) %>%
    rowSums()
  
  # Add events 
  columns_to_mutate <- pk_tb %>%
    ungroup() %>%
    select(contains("p")) %>%
    select(!contains("weighted")) %>%
    colnames()
  
  # Make plot
  pk_tb_plot <- rename(pk_tb, x = x1)
  plot_label <- if (pk_x_axis == "AUC") paste0(gsub("T-{0,2}", "T-", adc), " concentrations (µg/mL/day)") else paste0(gsub("T-{0,2}", "T-", adc), " concentrations (µg/mL)")
  pk_tb_plot <- pk_tb_plot %>%
    select(OCC, x, all_of(columns_to_mutate)) %>%
    pivot_longer(!c(OCC, x)) %>%
    mutate(name = gsub("^p", "", name))
  
  # Define manual color, line types, and names
  manual_colour_names <- c("CUI", "ORR", "PFS", "DLT")
  manual_colour_values <- c("CUI" = "black", "Efficacy: ORR" = "#00BFC4", "Efficacy: PFS" = "#C77CFF", "Safety: DLT" = "#F8766D")
  manual_line_width <- c("CUI" = 1.5, "Efficacy: ORR" = 1, "Efficacy: PFS" = 1, "Safety: DLT" = 1)
  manual_colour_values <- manual_colour_values[manual_colour_names %in% c("CUI", Outcomes)]
  manual_line_width <- manual_line_width[manual_colour_names %in% c("CUI", Outcomes)]
  
  weights_plot <- tibble(Outcome = gsub("p", "", names(weights)), Weight = weights) %>%
    arrange(desc(Weight)) %>%
    mutate(Outcome = factor(Outcome, levels = unique(Outcome)))
  
  # Reorder the colours (weights plot) to be consistent
  manual_colour_plot <- tibble(Outcome = names(manual_colour_values),
                               colour = manual_colour_values) %>%
    mutate(Outcome2 = Outcome,
           Outcome = gsub(".*: ", "", Outcome))
  manual_colour_plot <- weights_plot %>%
    mutate(Outcome = as.character(Outcome)) %>%
    left_join(manual_colour_plot) %>%
    pull(colour) %>%
    as.vector()
  
  if (adc == "T-DM1") {
    # Minimum and maximum values per dose
    range_0.3 <- pk_tb_plot %>% filter(OCC == "0.3") %>% pull(x) %>% range()
    range_0.6 <- pk_tb_plot %>% filter(OCC == "0.6") %>% pull(x) %>% range()
    range_1.2 <- pk_tb_plot %>% filter(OCC == "1.2") %>% pull(x) %>% range()
    range_1.8 <- pk_tb_plot %>% filter(OCC == "1.8") %>% pull(x) %>% range()
    range_2.4 <- pk_tb_plot %>% filter(OCC == "2.4") %>% pull(x) %>% range()
    range_3.6 <- pk_tb_plot %>% filter(OCC == "3.6") %>% pull(x) %>% range()
    range_4.8 <- pk_tb_plot %>% filter(OCC == "4.8") %>% pull(x) %>% range()
    range_tb <- tibble(dose = c(0.3, 0.6, 1.2, 1.8, 2.4, 3.6, 4.8),
                       min = c(range_0.3[1], range_0.6[1], range_1.2[1], range_1.8[1], range_2.4[1], range_3.6[1], range_4.8[1]),
                       max = c(range_0.3[2], range_0.6[2], range_1.2[2], range_1.8[2], range_2.4[2], range_3.6[2], range_4.8[2]))
  } else {
    # Minimum and maximum values per dose
    range_0.8 <- pk_tb_plot %>% filter(OCC == "0.8") %>% pull(x) %>% range()
    range_1.6 <- pk_tb_plot %>% filter(OCC == "1.6") %>% pull(x) %>% range()
    range_3.2 <- pk_tb_plot %>% filter(OCC == "3.2") %>% pull(x) %>% range()
    range_5.4 <- pk_tb_plot %>% filter(OCC == "5.4") %>% pull(x) %>% range()
    range_6.4 <- pk_tb_plot %>% filter(OCC == "6.4") %>% pull(x) %>% range()
    range_7.4 <- pk_tb_plot %>% filter(OCC == "7.4") %>% pull(x) %>% range()
    range_8 <- pk_tb_plot %>% filter(OCC == "8") %>% pull(x) %>% range()
    range_tb <- tibble(dose = c(0.8, 1.6, 3.2, 5.4, 6.4, 7.4, 8),
                       min = c(range_0.8[1], range_1.6[1], range_3.2[1], range_5.4[1], range_6.4[1], range_7.4[1], range_8[1]),
                       max = c(range_0.8[2], range_1.6[2], range_3.2[2], range_5.4[2], range_6.4[2], range_7.4[2], range_8[2]))
  }
  
  # Plots
  outcome_order <- as.character(weights_plot$Outcome) %>%
    gsub("ORR", "Efficacy: ORR", .) %>%
    gsub("PFS", "Efficacy: PFS", .) %>%
    gsub("DLT", "Safety: DLT", .) %>%
    tibble() %>%
    mutate(Outcome = weights_plot$Outcome)
  outcome_order <- weights_plot %>%
    left_join(outcome_order) %>%
    pull(`.`) 
  
  pk_tb_plot$name <- pk_tb_plot$name %>%
    gsub("ORR", "Efficacy: ORR", .) %>%
    gsub("PFS", "Efficacy: PFS", .) %>%
    gsub("DLT", "Safety: DLT", .)
  
  # Get CUI AUC
  auc_cui_plot <- pk_tb_plot %>%
    filter(name == "CUI") %>%
    arrange(x)
  
  AUC_tb <- tibble(dose = unique(auc_cui_plot$OCC), AUC = NA_real_)
  for (i in seq_along(unique(auc_cui_plot$OCC))) {
    auc_cui_plot_dose <- auc_cui_plot %>% filter(OCC == AUC_tb$dose[i])
    AUC_tb$AUC[i] <- trap.rule(auc_cui_plot_dose$x, auc_cui_plot_dose$value)/(range(auc_cui_plot_dose$x)[2]-range(auc_cui_plot_dose$x)[1])
    # Divide by range of exposure so that those with high variability are not favoured  
  }
  AUC <- AUC_tb %>%
    mutate(dose = as.numeric(as.character(dose))) %>%
    arrange(dose) %>%
    mutate(improvement = (AUC - lag(AUC)) * 100/lag(AUC)) %>%
    arrange(desc(AUC)) %>%
    mutate(dose = as.character(dose),
           label = as.character(round(AUC, 3)),
           label = if_else(nchar(label) == 5, label, paste0(label, "0")),
           label = if_else(nchar(label) == 5, label, paste0(label, "0")),
           label = if_else(nchar(label) == 5, label, paste0(label, "0")),
           label2 = as.character(round(improvement, 1)),
           label2 = if_else(str_detect(label2, "\\."), label2, paste0(label2, ".0")),
           label2 = paste0(label2, "%"),
           label2 = gsub("NA%", "NA", label2),
           label = paste0(label, " (", label2, ") "),
           label = gsub("\\(NA\\)", "", label)
    )
  
  AUC_cui <- AUC %>%
    arrange(dose) %>%
    mutate(improvement = if_else(is.na(improvement), cui_threshold, improvement)) %>%
    filter(improvement >= cui_threshold) 
  AUC_cui$improvement[1] <- AUC_cui$improvement[1] - 0.01 # Add this in case the first dose has the same improvement with subsequent doses
  AUC_cui <- AUC_cui %>% 
    arrange(desc(AUC))
  
  # Let's highlight the selected dose
  AUC$highlight <- ifelse(AUC$dose == AUC_cui$dose[[1]], "highlight", "normal")
  
  # Create the bar plot
  AUC_plot <- ggplot(AUC, aes(x = dose, y = AUC, fill = highlight)) +
    geom_bar(stat = "identity", colour = "#ABABAB") +
    scale_fill_manual(values = c("highlight" = "#ABABAB", "normal" = "white")) +
    labs(
      title = "Average CUI for different doses",
      x = "Dose (mg/kg)",
      y = "Average CUI"
    ) +
    theme_minimal() +
    theme(axis.text = element_text(size = 10),
          axis.title = element_text(size = 11, face = "bold"),
          plot.title = element_text(colour = "black", size = 12, face = "bold", hjust = 0.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "none") +
    geom_text(aes(label = ifelse(highlight == "highlight", label, ""), hjust = 1), colour = "black", angle = 0, size = 3) +
    geom_text(aes(label = ifelse(highlight != "highlight", label, ""), hjust = 1), colour = "black", angle = 0, size = 3) +
    coord_flip()
  
  # Create a custom labeling function to divide values by 10
  label_divided_by_1000 <- function(x) {
    scales::label_comma()(x / 1000) # Change from ng to µg for visual clarity
  }
  
  # Main plot
  dose_levels <- if (adc == "T-DM1") c("0.3", "0.6", "1.2", "1.8", "2.4", "3.6", "4.8") else c(0.8, 1.6, 3.2, 5.4, 6.4, 7.4, 8)
  pk_tb_plot <- pk_tb_plot %>%
    mutate(OCC = factor(as.character(OCC), levels = dose_levels), 
           name = factor(name, levels = c("CUI", outcome_order))) %>%
    ggplot() +
    geom_smooth(aes(x, value, color = name, linewidth = name), 
                method = "gam", formula = y ~ s(x, bs = "cs"), se = FALSE) +
    labs(x = plot_label,
         y = "Probability of event",
         color = "Legend") +
    theme_bw() +
    theme(axis.text = element_text(size = 10),
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title = element_text(size = 11, face = "bold"),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          plot.title = element_text(colour = "black", size = 12, face = "bold", hjust = 0.5),
          legend.position = "bottom") +
    scale_color_manual(values = manual_colour_values) +
    scale_linewidth_manual(values = manual_line_width) +
    guides(color = guide_legend("Legend"), 
           linewidth = guide_legend("Legend")) +
    ggtitle(paste0("CUI plot")) +
    scale_x_continuous(labels = label_divided_by_1000) +
    coord_cartesian(ylim = c(0, 1.05)) +
    annotate("text", 
             x = filter(range_tb, dose == as.numeric(AUC_cui$dose[1]))$min, 
             y = 0.95, 
             label = paste0("\n", AUC_cui$dose[1], " mg/kg"),  
             colour = "black", angle = 90, fontface = "italic") +
    annotate("rect", 
             xmin = filter(range_tb, dose == as.numeric(AUC_cui$dose[1]))$min, 
             xmax = filter(range_tb, dose == as.numeric(AUC_cui$dose[1]))$max, 
             ymin = -Inf, ymax = Inf, alpha = 0.2, fill = "black")
  
  # Weights plot
  weights_plot <- weights_plot %>%
    ggplot(aes(x = Outcome, y = Weight, fill = Outcome)) +
    geom_bar(stat = "identity", width = 0.5) +
    theme_minimal() +
    labs(title = "ORR and DLT weights", y = "Weight", x = "Outcomes") +
    theme(axis.text = element_text(size = 10),
          axis.text.x = element_blank(),
          axis.title = element_text(size = 11, face = "bold"),
          axis.title.x = element_blank(),
          plot.title = element_text(colour = "black", size = 12, face = "bold", hjust = 0.5)) +
    scale_y_continuous(labels = scales::percent_format()) +  # Format y-axis as percentage
    scale_fill_manual(values = manual_colour_plot) +
    guides(color = guide_legend("Legend"))
  
  layout <- "
    AAAAAAA#BB
    AAAAAAA#CC
    "
  overallplot <- pk_tb_plot +
    weights_plot + AUC_plot +
    plot_layout(design = layout)
  overallplot
  x <- if (adc == "T-DM1") 0.6 else 0.6
  graph2ppt(overallplot, paste0("results/", gsub("[^[:alnum:]]", "", paste0(adc, "_", scenario))), append = FALSE, width = 16 * x, height = 9 * x)
}
