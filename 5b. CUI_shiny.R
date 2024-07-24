# Load relevant libraries
library(shiny)
library(shinyWidgets)
library(shinyjs)
library(tidyverse)
library(bslib)
library(patchwork)
library(mgcv) # For fitting gam model

# Trapezoidal rule
trap.rule <- function(x,y) sum(diff(x) * (y[-1] + y[-length(y)])) / 2

# Process data
source("5a. Process Shiny data.R")

# Define UI
ui <- fluidPage(
  tags$head(
    tags$style(
      HTML("
        .bold-heading {
          font-weight: bold;
        }
        .normal-heading {
          font-weight: normal;
        }
        .sidebar {
          height: calc(100vh - 100px);
          overflow-y: auto;
          position: fixed;
          width: 300px;
          background-color: #f8f9fa;
          padding: 10px;
          top: 70px;
        }
        .sidebar-left {
          left: 0;
        }
        .sidebar-right {
          right: 0;
        }
        .main-content {
          margin-left: 320px;
          margin-right: 320px;
          padding: 10px;
          margin-top: 100px;
        }
        .title-panel {
          position: fixed;
          width: 100%;
          top: 0;
          left: 0;
          background-color: white;
          padding: 10px;
          z-index: 1000;
        }
      ")
    )
  ),
  useShinyjs(),  # Initialize shinyjs
  titlePanel("Clinical Utility Index (CUI)"),
  div(class = "sidebar sidebar-left",
      h4(class = "bold-heading", "PK metrics"),
      radioButtons("adc", 
                   label = HTML('<span style="font-weight: normal;">Select ADC:</span>'),
                   choices = c("T-DM1", "T-DXd"), selected = "T-DXd"),
      radioButtons("scenario", 
                   label = HTML('<span style="font-weight: normal;">Select Scenario:</span>'),
                   choices = c("Scenario 1 (only phase I)", "Scenario 2 (phases I and II)", "Scenario 3 (all phases)"), 
                   selected = "Scenario 3 (all phases)"),
      radioButtons("pk_x_axis", 
                   label = HTML('<span style="font-weight: normal;">Select PK metric:</span>'),
                   choices = c("AUC", "Cmax", "Cmin"), selected = "Cmax"),
      h4(class = "bold-heading", "Outcomes"),
      checkboxGroupInput("outcomes",
                         label = HTML('<span style="font-weight: normal;">Select Outcomes:</span>'),
                         choices = c("ORR", "PFS", "DLT"), 
                         selected =  c("ORR", "PFS", "DLT")),
  ),
  div(class = "main-content",
      plotOutput("mainPlot")
  ),
  div(class = "sidebar sidebar-right",
      # h4(class = "bold-heading", "CUI threshold"),
      # sliderInput("cui_threshold", HTML('<span style="font-size: 1.3em;">CUI threshold</span>'),
      #               min = 0, max = 20, value = 1),
      sliderTextInput("cui_threshold", HTML('<span style="font-size: 1.3em;">CUI threshold</span>'),
                      choices = as.character(c(0, seq(2.5, 20, 2.5))), selected = "2.5", grid = TRUE),
      # h4(class = "bold-heading", "Weights"),
      radioButtons("weights_choice", HTML('<span style="font-size: 1.3em;">Weights</span>'),
                   choices = c("Principal component analysis", "User weights"), 
                   selected = "User weights"),
      uiOutput("weights_panel")
  )
)

# Define server logic
server <- function(input, output, session) {

  # Reactive expressions to get selected PK metrics and outcomes
  selected_side_bar <- reactive({
    list(
      adc = input$adc,
      scenario = input$scenario,
      pk_x_axis = input$pk_x_axis
    )
  })
  
  selected_outcomes <- reactive({
    list(
      outcomes = input$outcomes
    )
  })
  
  selected_cui_threshold <- reactive({
    list(
      cui_threshold = input$cui_threshold
    )
  })

  # Dynamic UI for weights panel based on user selection
  output$weights_panel <- renderUI({
    outcomes <- selected_outcomes()$outcomes
    n <- length(outcomes)
    
    if (input$weights_choice == "User weights") {
      # Display UI for user-defined weights
      tags$div(
        h5(class = "bold-heading", "User-defined Weights:"),
        lapply(1:n, function(i) {
          sliderInput(paste0("weight_", i), 
                      label = outcomes[i], 
                      min = 1, max = 100, value = 50)
        })
      )
    }
  })
  
    # Observe changes in weights sliders for user-defined weights
  observe({
    req(input$weights_choice == "User weights")
    outcomes <- selected_outcomes()$outcomes
    weights <- sapply(seq_along(outcomes), function(i) input[[paste0("weight_", i)]])
    # Ensure weights sum to 100
    weights <- unlist(weights)
    if (length(weights) == 0) weights <- setNames(rep(50, length(outcomes)), outcomes)

    if (any(is.na(weights))) {
      weights[is.na(weights)] <- 50  # Set default value to 50 if NA found
    }

    total_weight <- sum(weights)
    if (abs(total_weight - 100) > 2) {
      weights <- 100 * weights / total_weight
      for (i in seq_along(outcomes)) {
        updateSliderInput(session, paste0("weight_", i), value = weights[i])
      }
    }
  })
  
  # Reactive to capture user-defined weights 
  user_weights <- reactive({
    req(input$weights_choice == "User weights")
    outcomes <- selected_outcomes()$outcomes
    weights <- unlist(sapply(seq_along(outcomes), function(i) input[[paste0("weight_", i)]]))
    if (length(weights) == 0) weights <- setNames(rep(50, length(outcomes)), outcomes)
    total_weight <- sum(weights)
    if (abs(total_weight - 100) > 2) weights <- weights / total_weight
    weights <- setNames(weights / 100, outcomes)
  })
  
  # Render main plot
  output$mainPlot <- renderPlot({
    req(selected_side_bar())

    adc <- selected_side_bar()$adc # Call reactive expression to get current values
    scenario <- selected_side_bar()$scenario
    pk_x_axis <- selected_side_bar()$pk_x_axis
    Outcomes <- selected_outcomes()$outcomes
    cui_threshold <- as.numeric(selected_cui_threshold()$cui_threshold)

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
    
    # Obtain weights
    if (input$weights_choice == "Principal component analysis") {
      # Get weights using PCA
      pca_dat <- select(ungroup(pk_tb), -id, -OCC, -AMT, -x1, -x2)
      pca_result <- prcomp(pca_dat)
      weights <- (pca_result$rotation[, 1])^2 # Squared loadings used as weights
    } else {
      weights <- user_weights()
      new_names <- colnames(select(ungroup(pk_tb), -id, -OCC, -AMT, -x1, -x2))
      weights <- setNames(weights, new_names)
    }
    
    # Add weights
    if ("ORR" %in% Outcomes) pk_tb$weightedpORR <- pk_tb$pORR * weights[["pORR"]]
    if ("PFS" %in% Outcomes) pk_tb$weightedpPFS <- pk_tb$pPFS * weights[["pPFS"]]
    if ("DLT" %in% Outcomes) pk_tb$weightedpDLT <- (1 - pk_tb$pDLT) * weights[["pDLT"]]
    
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
    plot_label <- if (pk_x_axis == "AUC") paste0(gsub("T-{0,2}", "T-", adc), " concentrations (ng/mL/day)") else paste0(gsub("T-{0,2}", "T-", adc), " concentrations (ng/mL)")
    pk_tb_plot <- pk_tb_plot %>%
      select(OCC, x, all_of(columns_to_mutate)) %>%
      pivot_longer(!c(OCC, x)) %>%
      mutate(name = gsub("^p", "", name))
    
    # Define manual color, line types, and names
    manual_colour_names <- c("CUI", "ORR", "PFS", "DLT")
    manual_colour_values <- c("CUI" = "black", "Efficacy: ORR" = "#619CFF", "Efficacy: PFS" = "lightgreen", "Safety: DLT" = "#F8766D")
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
        title = "CUI AUC/exposure range\nfor different doses",
        x = "Dose (mg/kg)",
        y = "CUI AUC/exposure range"
      ) +
      theme_minimal() +
      theme(axis.text = element_text(size = 10),
            axis.title = element_text(size = 12, face = "bold"),
            plot.title = element_text(colour = "black", size = 14, face = "bold", hjust = 0.5),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.position = "none") +
      geom_text(aes(label = ifelse(highlight == "highlight", label, ""), hjust = 1), colour = "black", angle = 0, size = 3) +
      geom_text(aes(label = ifelse(highlight != "highlight", label, ""), hjust = 1), colour = "black", angle = 0, size = 3) +
      coord_flip()
    
    # Main plot
    pk_tb_plot <- pk_tb_plot %>%
      mutate(OCC = factor(as.character(OCC), levels = c("0.3", "0.6", "1.2", "1.8", "2.4", "3.6", "4.8")), 
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
            axis.title = element_text(size = 12, face = "bold"),
            panel.grid.major.x = element_blank(),
            panel.grid.minor.x = element_blank(),
            plot.title = element_text(colour = "black", size = 14, face = "bold", hjust = 0.5),
            legend.position = "bottom") +
      scale_color_manual(values = manual_colour_values) +
      scale_linewidth_manual(values = manual_line_width) +
      guides(color = guide_legend("Legend"), 
             linewidth = guide_legend("Legend")) +
      ggtitle(paste0(pk_x_axis, " as x-axis")) +
      scale_x_continuous(labels = scales::label_scientific()) +
      coord_cartesian(ylim = c(0, 1.05)) +
      annotate("text", 
               x = filter(range_tb, dose == as.numeric(AUC_cui$dose[1]))$min, 
               y = 0.95, 
               label = paste0("\n", AUC_cui$dose[1], " mg/kg"),  
               colour = "black", angle = 90, fontface = "italic")  +
      annotate("rect", 
               xmin = filter(range_tb, dose == as.numeric(AUC_cui$dose[1]))$min, 
               xmax = filter(range_tb, dose == as.numeric(AUC_cui$dose[1]))$max, 
               ymin = -Inf, ymax = Inf, alpha = 0.2, fill = "black")

    # Weights plot
    weights_plot <- weights_plot %>%
      ggplot(aes(x = Outcome, y = Weight, fill = Outcome)) +
      geom_bar(stat = "identity", width = 0.5) +
      theme_minimal() +
      labs(title = "Outcome Weights", y = "Weight", x = "Outcomes") +
      theme(axis.text = element_text(size = 10),
            axis.text.x = element_blank(),
            axis.title = element_text(size = 12, face = "bold"),
            axis.title.x = element_blank(),
            plot.title = element_text(colour = "black", size = 14, face = "bold", hjust = 0.5)) +
      # coord_cartesian(ylim = c(0, 1)) +
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
  })
}

# Run the application
shinyApp(ui = ui, server = server)
