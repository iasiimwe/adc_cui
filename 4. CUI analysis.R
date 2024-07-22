# Load relevant libraries
library(tidyverse)
library(RColorBrewer)
library(patchwork)
library(mgcv) # For fitting gam model
library(ggbump)

# Trapezoidal rule
trap.rule <- function(x,y) sum(diff(x) * (y[-1] + y[-length(y)])) / 2

# Select the best CUI threshold
cui_improvements <- seq(2.5, 20, 2.5) # improvement in CUI needed to pick the next dose (at least 5%)
results_doses <- tibble()
for (ADC in c("TDM1", "TDXd")) { # Loop through the ADCs
  # Process data
  source("3. ER analysis.R")
  
  er_models <- list(er_mod_tb1, er_mod_tb2, er_mod_tb3)
  scenarios <- c("Scenario 1 (only phase I)", "Scenario 2 (phases I and II)", "Scenario 3 (all phases)")
  scenario_plots <- vector("list", length(scenarios))
  
  for (cui_improvement in cui_improvements) { 
    for (PK_metric in c("Cmax")) {
      pk_eff <- pk_safety <- pk_x_axis <- PK_metric
      for (i in seq_along(scenarios)) {
        er_mod_tb <- er_models[[i]] %>%
          filter(Outcome %in% c("ORR", "DLT")) %>%
          select(Outcome, all_of(PK_metric))
        
        # Assign weights to ORR from 0 to 100%; DLT will be 100% - assigned weight
        results_wt <- c()
        for (j in 10:90) results_wt <- c(results_wt, paste0("ORR: ", j/100, ", DLT: ", 1 - j/100))
        
        # Get the best dose per combination using average CUI
        results_tb <- tibble(Covs = results_wt,
                             best_CUI = NA_real_,
                             best_dose = NA_real_)
        for (j in seq_along(results_tb$Covs)) {
          Outcomes <- unlist(str_split(results_tb$Covs[j], ", "))
          weights <- parse_number(Outcomes)
          Outcomes <- gsub(":| |\\d|\\.", "", Outcomes)
          names(weights) <- Outcomes
          
          for (k in seq_along(c("ORR", "DLT"))) {
            outcomei <- c("ORR", "DLT")[k]
            if (outcomei %in% Outcomes) {
              PK_param <- if (outcomei %in% c("ORR_", "OS_", "PFS_")) pk_eff else pk_safety
              assign(paste0(gsub("_", "", outcomei), "i"), 
                     unlist(er_mod_tb %>% filter(Outcome == outcomei) %>% pull(PK_param))[1], envir = .GlobalEnv) # Intercept
              assign(paste0(gsub("_", "", outcomei), "s"), 
                     unlist(er_mod_tb %>% filter(Outcome == outcomei) %>% pull(PK_param))[2], envir = .GlobalEnv) # Slope
            }
          }
          
          # Link with the PK data
          pk_tb <- pk_tb_data %>%
            group_by(id, OCC, AMT) %>%
            reframe(AUC = max(AUC),
                    Cmax = max(Cc_wghtd),
                    Cmin = Cc_wghtd[time == 504], .groups = "drop") %>%
            select(id:AMT, all_of(unique(c(pk_eff, pk_safety))))
          if (pk_eff == pk_safety) {
            colnames(pk_tb)[colnames(pk_tb) == pk_eff] <- "x1"
            pk_tb$x2 <- pk_tb$x1
          } else {
            colnames(pk_tb)[colnames(pk_tb) == pk_eff] <- "x1"
            colnames(pk_tb)[colnames(pk_tb) == pk_safety] <- "x2"
          }
          
          if ("ORR" %in% Outcomes) pk_tb <- pk_tb %>% mutate(pORR = inv_logit_fn(ORRi + ORRs * x1)) # # probability of ORR
          if ("DLT" %in% Outcomes) pk_tb <- pk_tb %>% mutate(pDLT = inv_logit_fn(DLTi + DLTs * x2))
          
          # Add weights
          if ("ORR" %in% Outcomes) pk_tb$weightedpORR <- pk_tb$pORR * weights[["ORR"]]
          if ("DLT" %in% Outcomes) pk_tb$weightedpDLT <- (1 - pk_tb$pDLT) * weights[["DLT"]]
          
          # Compute clinical utility index based on weights
          pk_tb$pCUI <- pk_tb %>%
            ungroup() %>%
            select(contains("weighted")) %>%
            rowSums()
          
          # Add events 
          # set.seed(7) # For reproducibility
          columns_to_mutate <- pk_tb %>%
            ungroup() %>%
            select(contains("p")) %>%
            select(!contains("weighted")) %>%
            colnames()
          # add_events_fn <- function(x) sample(c(0, 1), 1, prob = c(1 - x, x))
          # pk_tb <- pk_tb %>%
          #   rowwise() %>%
          #   mutate_at(vars(all_of(columns_to_mutate)), add_events_fn)
          
          pk_tb_plot <- if (pk_x_axis == pk_eff) rename(pk_tb, x = x1) else rename(pk_tb, x = x2)
          pk_tb_plot <- pk_tb_plot %>%
            select(OCC, x, all_of(columns_to_mutate)) %>%
            pivot_longer(!c(OCC, x)) %>%
            mutate(name = gsub("^p", "", name))
          
          if (ADC == "T-DM1") {
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
          
          # Get CUI AUC
          # When you use method = "gam", formula = y ~ s(x, bs = "cs") in geom_smooth(), 
          # ggplot2 will fit a Generalized Additive Model to your data with a smooth curve using cubic regression splines. 
          # This can be particularly useful for capturing complex, non-linear relationships between the variables.
          auc_cui_plot <- pk_tb_plot %>%
            filter(name == "CUI") %>%
            arrange(x)
          # gam_model <- gam(value ~ s(x, bs = "cs"), data = auc_cui_plot)
          # auc_cui_plot$predicted_y <- predict(gam_model, newdata = auc_cui_plot)
          AUC_tb <- tibble(dose = unique(auc_cui_plot$OCC), AUC = NA_real_)
          for (l in seq_along(unique(auc_cui_plot$OCC))) {
            auc_cui_plot_dose <- auc_cui_plot %>% filter(OCC == AUC_tb$dose[l])
            # AUC_tb$AUC[l] <- trap.rule(auc_cui_plot_dose$x, auc_cui_plot_dose$predicted_y)/(range(auc_cui_plot_dose$x)[2]-range(auc_cui_plot_dose$x)[1])
            AUC_tb$AUC[l] <- trap.rule(auc_cui_plot_dose$x, auc_cui_plot_dose$value)/(range(auc_cui_plot_dose$x)[2]-range(auc_cui_plot_dose$x)[1])
                # Divide by range of exposure so that those with high variability are not favoured  
          }
          AUC <- AUC_tb %>% 
            arrange(dose) %>%
            mutate(improvement = (AUC - lag(AUC)) * 100/lag(AUC),
                   improvement = if_else(is.na(improvement), cui_improvement, improvement)) %>%
            filter(improvement >= cui_improvement) 
          AUC$improvement[1] <- AUC$improvement[1] - 0.01 # Add this in case the first dose has the same improvement with subsequent doses
          AUC <- AUC %>% 
            arrange(desc(AUC)) 
          results_tb$best_CUI[j] <- round(AUC$AUC[1], 3)
          results_tb$best_dose[j] <- as.numeric(as.character(AUC$dose[1]))
          
          message(paste0(ADC, ", CUI threshold ", cui_improvement, ": ", round(j * 100/length(results_tb$Covs), 3), "% complete!"))
        }
        
        # Plot for the most selected doses
        if (ADC == "TDM1") {
          doses <- results_tb %>%
            mutate(`0.3 mg/kg` = if_else(best_dose == 0.3, TRUE, FALSE),
                   `0.6 mg/kg` = if_else(best_dose == 0.6, TRUE, FALSE),
                   `1.2 mg/kg` = if_else(best_dose == 1.2, TRUE, FALSE),
                   `1.8 mg/kg` = if_else(best_dose == 1.8, TRUE, FALSE),
                   `2.4 mg/kg` = if_else(best_dose == 2.4, TRUE, FALSE),
                   `3.6 mg/kg` = if_else(best_dose == 3.6, TRUE, FALSE),
                   `4.8 mg/kg` = if_else(best_dose == 4.8, TRUE, FALSE)) %>%
            select(contains("mg"))
        } else {
          doses <- results_tb %>%
            mutate(`0.8 mg/kg` = if_else(best_dose == 0.8, TRUE, FALSE),
                   `1.6 mg/kg` = if_else(best_dose == 1.6, TRUE, FALSE),
                   `3.2 mg/kg` = if_else(best_dose == 3.2, TRUE, FALSE),
                   `5.4 mg/kg` = if_else(best_dose == 5.4, TRUE, FALSE),
                   `6.4 mg/kg` = if_else(best_dose == 6.4, TRUE, FALSE),
                   `7.4 mg/kg` = if_else(best_dose == 7.4, TRUE, FALSE),
                   `8 mg/kg` = if_else(best_dose == 8, TRUE, FALSE)) %>%
            select(contains("mg"))
        }
        
        doses <- tibble(dose = colnames(doses),
                        Percentage = colSums(doses)) %>%
          mutate(Percentage = round(Percentage * 100/nrow(results_tb), 0),
                 label = paste0(as.character(Percentage), "%"),
                 label = if_else(label == "0%", "", label)) %>%
          mutate(ADC = ADC,
                 Scenario = scenarios[i],
                 `CUI threshold` = cui_improvement)
        results_doses <- bind_rows(results_doses, doses)
      }
    }
  }
}
write_rds(results_doses, paste0("results/", "cui_thresholds.rds")) # Save this in case needed later

font_size <- 10
# x_intercept <- pull(results_doses %>%
#   filter(Scenario == "Scenario 3 (all phases)", ADC == "TDM1", dose %in% c("3.6 mg/kg")) %>%
#   group_by(dose) %>%
#   reframe(`CUI threshold` = `CUI threshold`[Percentage == max(Percentage)],
#           Percentage = max(Percentage)), `CUI threshold`)[1]
bumpplot_tdm1 <- results_doses %>%
  rename(Dose = dose) %>%
  filter(ADC == "TDM1") %>% 
  mutate(approved_dose = if_else(Dose %in% c("3.6 mg/kg"), "Yes", "No")) %>%
  ggplot(aes(x = `CUI threshold`, y = Percentage)) +
  geom_bump(aes(color = Dose, group = Dose, alpha = approved_dose, linewidth = approved_dose), show.legend = FALSE) +
  geom_point(aes(color = Dose, size = approved_dose, shape = Dose), 
             show.legend = c(alpha = FALSE, size = FALSE, linewidth = FALSE)) +
  facet_wrap(~ Scenario, nrow = 1) +
  scale_size_manual(values = c("Yes" = 2, "No" = 1)) +
  scale_alpha_manual(values = c("Yes" = 1, "No" = 0.2)) +
  scale_linewidth_manual(values = c("Yes" = 1, "No" = 0.7)) +
  scale_shape_manual(values = c(1, 1, 1, 1, 1, 16, 1)) +
  guides(color = guide_legend(title = "Dose Legend"), 
         shape = guide_legend(title = "Dose Legend")) +
  theme_bw() +
  theme(axis.title = element_text(size = font_size, face = "bold"),
        plot.title = element_text(size = font_size + 1, face = "bold"),
        axis.text = element_text(size = font_size),
        legend.position = "right",
        strip.text = element_text(size = font_size, colour = "black"))+
  labs(x = "CUI threshold (%)", y = "Percent selection", title = "A. T-DM1") 
# + geom_vline(xintercept = x_intercept, linetype = "dashed", color = "black")
# x_intercept <- pull(results_doses %>%
#                       filter(Scenario == "Scenario 3 (all phases)", ADC == "TDXd", dose %in% c("5.4 mg/kg", "6.4 mg/kg")) %>%
#                       group_by(`CUI threshold`) %>%
#                       summarise(Percentage = mean(Percentage)) %>%
#                       arrange(desc(Percentage)), `CUI threshold`)[1]
bumpplot_tdxd <- results_doses %>%
  rename(Dose = dose) %>%
  filter(ADC == "TDXd") %>% 
  mutate(approved_dose = if_else(Dose %in% c("5.4 mg/kg", "6.4 mg/kg"), "Yes", "No")) %>%
  ggplot(aes(x = `CUI threshold`, y = Percentage)) +
  geom_bump(aes(color = Dose, group = Dose, alpha = approved_dose, linewidth = approved_dose), show.legend = FALSE) +
  geom_point(aes(color = Dose, size = approved_dose, shape = Dose), 
             show.legend = c(alpha = FALSE, size = FALSE, linewidth = FALSE)) +
  facet_wrap(~ Scenario, nrow = 1) +
  scale_size_manual(values = c("Yes" = 2, "No" = 1)) +
  scale_alpha_manual(values = c("Yes" = 1, "No" = 0.2)) +
  scale_linewidth_manual(values = c("Yes" = 1, "No" = 0.7)) +
  scale_shape_manual(values = c(1, 1, 1, 16, 16, 1, 1)) +
  guides(color = guide_legend(title = "Dose Legend"), 
         shape = guide_legend(title = "Dose Legend")) +
  theme_bw() +
  theme(axis.title = element_text(size = font_size, face = "bold"),
        plot.title = element_text(size = font_size + 1, face = "bold"),
        axis.text = element_text(size = font_size),
        legend.position = "right",
        strip.text = element_text(size = font_size, colour = "black"))+
  labs(x = "CUI threshold (%)", y = "Percent selection", title = "B. T-DXd") 
# + geom_vline(xintercept = x_intercept, linetype = "dashed", color = "black")
combined_plot <- bumpplot_tdm1/bumpplot_tdxd
# Save plot
x <- 0.8
graph2ppt(combined_plot, paste0("results/", "CUI_thresholds"), append = FALSE, width = 12 * x, height = 9 * x)


# Get plots for the best scenarios using thresholds that optimize the selection of the approved doses
cutoff <- 4 # Plot the top 4 doses
cui_improvements <- c(10) # Aim for at least a 5% improvement
for (ADC in c("TDM1", "TDXd")) { # Loop through the ADCs
  # Process data
  source("3. ER analysis.R")
  
  er_models <- list(er_mod_tb1, er_mod_tb2, er_mod_tb3)
  scenarios <- c("Scenario 1 (only phase I)", "Scenario 2 (phases I and II)", "Scenario 3 (all phases)")
  scenario_plots <- vector("list", length(scenarios))
  
  for (cui_improvement in cui_improvements) { 
    for (PK_metric in c("Cmax")) {
      pk_eff <- pk_safety <- pk_x_axis <- PK_metric
      for (i in seq_along(scenarios)) {
        er_mod_tb <- er_models[[i]] %>%
          filter(Outcome %in% c("ORR", "DLT")) %>%
          select(Outcome, all_of(PK_metric))
        
        # Assign weights to ORR from 0 to 100%; DLT will be 100% - assigned weight
        results_wt <- c()
        for (j in 10:90) results_wt <- c(results_wt, paste0("ORR: ", j/100, ", DLT: ", 1 - j/100))
        
        # Get the best dose per combination using average CUI
        results_tb <- tibble(Covs = results_wt,
                             best_CUI = NA_real_,
                             best_dose = NA_real_)
        for (j in seq_along(results_tb$Covs)) {
          Outcomes <- unlist(str_split(results_tb$Covs[j], ", "))
          weights <- parse_number(Outcomes)
          Outcomes <- gsub(":| |\\d|\\.", "", Outcomes)
          names(weights) <- Outcomes
          
          for (k in seq_along(c("ORR", "DLT"))) {
            outcomei <- c("ORR", "DLT")[k]
            if (outcomei %in% Outcomes) {
              PK_param <- if (outcomei %in% c("ORR_", "OS_", "PFS_")) pk_eff else pk_safety
              assign(paste0(gsub("_", "", outcomei), "i"), 
                     unlist(er_mod_tb %>% filter(Outcome == outcomei) %>% pull(PK_param))[1], envir = .GlobalEnv) # Intercept
              assign(paste0(gsub("_", "", outcomei), "s"), 
                     unlist(er_mod_tb %>% filter(Outcome == outcomei) %>% pull(PK_param))[2], envir = .GlobalEnv) # Slope
            }
          }
          
          # Link with the PK data
          pk_tb <- pk_tb_data %>%
            group_by(id, OCC, AMT) %>%
            reframe(AUC = max(AUC),
                    Cmax = max(Cc_wghtd),
                    Cmin = Cc_wghtd[time == 504], .groups = "drop") %>%
            select(id:AMT, all_of(unique(c(pk_eff, pk_safety))))
          if (pk_eff == pk_safety) {
            colnames(pk_tb)[colnames(pk_tb) == pk_eff] <- "x1"
            pk_tb$x2 <- pk_tb$x1
          } else {
            colnames(pk_tb)[colnames(pk_tb) == pk_eff] <- "x1"
            colnames(pk_tb)[colnames(pk_tb) == pk_safety] <- "x2"
          }
          
          # Get probability
          if ("ORR" %in% Outcomes) pk_tb <- pk_tb %>% mutate(pORR = inv_logit_fn(ORRi + ORRs * x1)) # # probability of ORR
          if ("DLT" %in% Outcomes) pk_tb <- pk_tb %>% mutate(pDLT = inv_logit_fn(DLTi + DLTs * x2))

          # Add weights
          if ("ORR" %in% Outcomes) pk_tb$weightedpORR <- pk_tb$pORR * weights[["ORR"]]
          if ("DLT" %in% Outcomes) pk_tb$weightedpDLT <- (1 - pk_tb$pDLT) * weights[["DLT"]]
          
          # Compute clinical utility index based on weights
          pk_tb$pCUI <- pk_tb %>%
            ungroup() %>%
            select(contains("weighted")) %>%
            rowSums()
          
          # # Add events
          # set.seed(7) # For reproducibility
          columns_to_mutate <- pk_tb %>%
            ungroup() %>%
            select(contains("p")) %>%
            select(!contains("weighted")) %>%
            colnames()

          # Use the probabilities instead of events
          pk_tb_plot <- if (pk_x_axis == pk_eff) rename(pk_tb, x = x1) else rename(pk_tb, x = x2)
          pk_tb_plot <- pk_tb_plot %>%
            select(OCC, x, all_of(columns_to_mutate)) %>%
            pivot_longer(!c(OCC, x)) %>%
            mutate(name = gsub("^p", "", name))
          
          if (ADC == "T-DM1") {
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
          
          # Get CUI AUC
          # When you use method = "gam", formula = y ~ s(x, bs = "cs") in geom_smooth(), 
          # ggplot2 will fit a Generalized Additive Model to your data with a smooth curve using cubic regression splines. 
          # This can be particularly useful for capturing complex, non-linear relationships between the variables.
          auc_cui_plot <- pk_tb_plot %>%
            filter(name == "CUI") %>%
            arrange(x)
          # gam_model <- gam(value ~ s(x, bs = "cs"), data = auc_cui_plot)
          # auc_cui_plot$predicted_y <- predict(gam_model, newdata = auc_cui_plot)
          AUC_tb <- tibble(dose = unique(auc_cui_plot$OCC), AUC = NA_real_)
          for (l in seq_along(unique(auc_cui_plot$OCC))) {
            auc_cui_plot_dose <- auc_cui_plot %>% filter(OCC == AUC_tb$dose[l])
            # AUC_tb$AUC[l] <- trap.rule(auc_cui_plot_dose$x, auc_cui_plot_dose$predicted_y)/(range(auc_cui_plot_dose$x)[2]-range(auc_cui_plot_dose$x)[1])
            AUC_tb$AUC[l] <- trap.rule(auc_cui_plot_dose$x, auc_cui_plot_dose$value)/(range(auc_cui_plot_dose$x)[2]-range(auc_cui_plot_dose$x)[1])
            # Divide by range of exposure so that those with high variability are not favoured  
          }
          AUC <- AUC_tb %>% 
            arrange(dose) %>%
            mutate(improvement = (AUC - lag(AUC)) * 100/lag(AUC),
                   improvement = if_else(is.na(improvement), cui_improvement, improvement)) %>%
            filter(improvement >= cui_improvement) 
          AUC$improvement[1] <- AUC$improvement[1] - 0.01 # Add this in case the first dose has the same improvement with subsequent doses
          AUC <- AUC %>% 
            arrange(desc(AUC)) 
          results_tb$best_CUI[j] <- round(AUC$AUC[1], 3)
          results_tb$best_dose[j] <- as.numeric(as.character(AUC$dose[1]))
          
          message(paste0(round(j * 100/length(results_tb$Covs), 3), "% complete!"))
        }
        
        # Plot for the most selected doses
        if (ADC == "TDM1") {
          doses <- results_tb %>%
            mutate(`0.3 mg/kg` = if_else(best_dose == 0.3, TRUE, FALSE),
                   `0.6 mg/kg` = if_else(best_dose == 0.6, TRUE, FALSE),
                   `1.2 mg/kg` = if_else(best_dose == 1.2, TRUE, FALSE),
                   `1.8 mg/kg` = if_else(best_dose == 1.8, TRUE, FALSE),
                   `2.4 mg/kg` = if_else(best_dose == 2.4, TRUE, FALSE),
                   `3.6 mg/kg` = if_else(best_dose == 3.6, TRUE, FALSE),
                   `4.8 mg/kg` = if_else(best_dose == 4.8, TRUE, FALSE)) %>%
            select(contains("mg"))
        } else {
          doses <- results_tb %>%
            mutate(`0.8 mg/kg` = if_else(best_dose == 0.8, TRUE, FALSE),
                   `1.6 mg/kg` = if_else(best_dose == 1.6, TRUE, FALSE),
                   `3.2 mg/kg` = if_else(best_dose == 3.2, TRUE, FALSE),
                   `5.4 mg/kg` = if_else(best_dose == 5.4, TRUE, FALSE),
                   `6.4 mg/kg` = if_else(best_dose == 6.4, TRUE, FALSE),
                   `7.4 mg/kg` = if_else(best_dose == 7.4, TRUE, FALSE),
                   `8 mg/kg` = if_else(best_dose == 8, TRUE, FALSE)) %>%
            select(contains("mg"))
        }
        
        doses <- tibble(dose = colnames(doses),
                        Percentage = colSums(doses)) %>%
          mutate(Percentage = round(Percentage * 100/nrow(results_tb), 0),
                 label = paste0(as.character(Percentage), "%"),
                 label = if_else(label == "0%", "", label))
        bar_plot <- ggplot(doses, aes(x = Percentage, y = dose)) +
          geom_bar(stat = "identity", fill = "steelblue") +
          theme_bw() +
          scale_y_discrete(position = "right") + # Move y-axis to the right
          scale_x_reverse() + # Reverse the x-axis
          labs(x = "Percent selection of a given dose") +
          ggtitle(paste0(ifelse(i == 1, "A. ", ifelse(i == 2, "B. ", "C. ")), gsub("T", "T-", ADC), ": ", scenarios[[i]]), 
                  subtitle = "Percent selection of a given dose") +
          theme(axis.title.y = element_blank(), # Remove y-axis title
                axis.text = element_text(colour = "black"),
                axis.text.y = element_text(hjust = 1), # Right-align y-axis text
                axis.ticks.y = element_blank()) + # Remove y-axis ticks
          geom_text(aes(label = if_else(Percentage >= 10, label, "")), hjust = -0.1, color="white") +
          geom_text(aes(label = if_else(Percentage < 10, label, "")), hjust = 1.2, color="black")
        
        manual_colour_values_all <- c(brewer.pal(3, "Blues")[3], brewer.pal(4, "Reds")[4])
        order_all <- c("ORR", "DLT")
        
        doses <- doses %>%
          filter(Percentage > 0) %>% 
          arrange(desc(Percentage)) %>% # Arrange by selection frequency (most selected to least selected)
          mutate(dose_n = parse_number(dose))
        doses <- doses[1:min(nrow(doses), cutoff), ]
        
        ORR_weight <- results_tb %>%
          filter(best_dose %in% doses$dose_n) %>%
          group_by(best_dose) %>%
          mutate(n = row_number()) %>%
          separate_rows(Covs, sep = ", ") %>%
          separate(Covs, into = c("Covariate", "Percentage"), sep = ": ") %>%
          filter(Covariate == "ORR") %>%
          mutate(Percentage = as.numeric(Percentage)) %>%
          summarize(median_wt = median(Percentage) * 100, 
                    min_wt = min(Percentage) * 100,
                    max_wt = max(Percentage) * 100) %>%
          mutate(range = paste0("ORR\nweight: median ",median_wt, ", range ", min_wt, " to ", max_wt, "%)"))
        
        plot_data <- results_tb %>%
          filter(best_dose %in% doses$dose_n) %>%
          rowwise() %>%
          mutate(best_dose = paste0(doses$dose[doses$dose_n == best_dose], " (",
                                    doses$label[doses$dose_n == best_dose], " selection; ",
                                    ORR_weight$range[ORR_weight$best_dose == best_dose])) %>%
          group_by(best_dose) %>%
          mutate(n = row_number()) %>%
          separate_rows(Covs, sep = ", ") %>%
          separate(Covs, into = c("Covariate", "Percentage"), sep = ": ") %>%
          mutate(Percentage = as.numeric(Percentage) * 100,
                 Covariate = factor(Covariate, levels = order_all)) %>%
          rename(`Outcome` = Covariate)
        weights_plot <- ggplot(plot_data, aes(x = factor(n), y = Percentage, fill = `Outcome`)) +
          facet_wrap(~ best_dose, nrow = 1, scales = "free_x") + 
          geom_bar(stat = "identity", width = 1.1) +
          xlab("Different weight simulations") + 
          ylab("Percentage weight") +
          ggtitle(" ", subtitle = paste0("Top ", min(nrow(doses), cutoff), " doses")) +
          theme_minimal() +
          theme(axis.text.y = element_text(size = 10), 
                axis.text.x = element_blank(), 
                axis.title = element_text(size = 12, face = "bold"),
                axis.ticks.x = element_blank(),
                plot.title = element_text(hjust = 0, vjust = 0, size = 14, face = "bold"),
                strip.text = element_text(size = 9, margin = margin(0.05,0,0.05,0, "cm")),
                legend.title = element_text(size = 10, face = "bold"),
                legend.text = element_text(size = 10),
                panel.grid.major = element_blank(), # hide grid lines
                panel.grid.minor = element_blank()
          ) +
          scale_y_continuous(expand=expansion(mult=c(0, 0.01)))
        
        layout <- "
ABBB
"
        font_size <- 10
        scenario_plots[[i]] <- (bar_plot +
                                  theme(axis.text.y.left = element_text(size = font_size),
                                        axis.title.x = element_text(size = font_size, face = "bold"),
                                        plot.title = element_text(size = font_size + 1, face = "bold")
                                  )) +
          (weights_plot +
             theme(axis.text.y.left = element_text(size = font_size),
                   axis.title = element_text(size = font_size, face = "bold"),
                   legend.title = element_text(size = font_size - 1, face = "bold")
             )) +
          plot_layout(design = layout)
      }
      combined_plot <- (scenario_plots[[1]]/scenario_plots[[2]])/scenario_plots[[3]]
      # Save plot
      x <- if (ADC == "TDM1") 0.9 else 0.9
      graph2ppt(combined_plot, paste0("results/", ADC, "_", PK_metric, "_", gsub("\\.", "_", as.character(cui_improvement)), "_CUI"), append = FALSE, width = 16 * x, height = 12 * x)
    }
  }
}
