# Import libraries and initialize Monolix
# ---------------------------------------
library(lixoftConnectors)
initializeLixoftConnectors(software = "monolix", path = "C:/Program Files/Lixoft/MonolixSuite2024R1", force = TRUE)
library(tidyverse)
library(export)
library(readxl)
library(PKNCA) # For calculating non-compartmental analysis

# Analysis date
analysis_date <- "06_07_2024" 

# Loop through the ADCs
for (ADC in c("TDM1", "TDXd")) {
  
  ## Bayesian individual dynamic predictions, https://monolix.lixoft.com/monolix-api/examples/
  # ----------------------------------------
  initializeLixoftConnectors(software = "monolix", path = "C:/Program Files/Lixoft/MonolixSuite2024R1", force = TRUE)
  loadProject(projectFile = paste0("results/", ADC, "_base.mlxtran"))
  
  # Set tasks in scenario
  scenario <- getScenario()
  scenario$tasks = c(populationParameterEstimation = TRUE,
                     conditionalDistributionSampling = TRUE,
                     conditionalModeEstimation = TRUE,
                     standardErrorEstimation = TRUE,
                     logLikelihoodEstimation = FALSE,
                     plots = FALSE)
  scenario$linearization = FALSE
  setScenario(scenario)
  
  # set pop params to estimates and fix pop params
  saveProject(paste0("results/", ADC, "_base2.mlxtran"))
  
  runPopulationParameterEstimation()
  setInitialEstimatesToLastEstimates()
  popparams <- getPopulationParameterInformation()
  popparams$method <- "FIXED"
  setPopulationParameterInformation(popparams)
  runScenario()
  
  # set MCMC settings
  setConditionalDistributionSamplingSettings(nbminiterations = 1000, ratio = 0.05, nbsimulatedparameters = 1000)
  
  # run SAEM and conditional distribution
  runPopulationParameterEstimation() # mandatory before other tasks, but nb of iterations is null as all parameters are fixed
  runConditionalDistributionSampling() # this is sampling from the posterior distribution for each individual
  
  simparams <- getSimulatedIndividualParameters()
  dict <- read_csv(paste0("results/", ADC, "_poppk_", analysis_date, ".csv"), show_col_types = FALSE) %>%
    select(id = ID, OCC) %>%
    mutate(OCC = parse_number(OCC)) %>%
    left_join(mutate(select(simparams, rep, id), id = as.numeric(as.character(id))), relationship = "many-to-many") %>%
    select(id = rep, OCC)
  simparams$id <- row.names(simparams) # replace id by rank of {rep, id}
  simparams$rep <- NULL
  write.csv(simparams, file = paste0("results/", ADC, "_table_simulated_parameters.csv"), row.names = FALSE)
  write.csv(dict, file =  paste0("results/", ADC, "_table_simulated_parameters_dict.csv"), row.names = FALSE)
  saveProject(paste0("results/", ADC, "_base2.mlxtran"))
  
  # simulate based on simulated parameters using Simulx
  initializeLixoftConnectors(software = "simulx", path = "C:/Program Files/Lixoft/MonolixSuite2024R1", force = TRUE)
  importProject(paste0("results/", ADC, "_base2.mlxtran"))
  
  ## Deal with both AUC and Cmax
  setAddLines("ddt_AUC = Cc_wghtd") # define a new element output with a regular grid for the variable "AUC" 
  defineOutputElement(name = "Cmax", element = list(data = data.frame(time = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 504)), output = "Cc_wghtd"))
  defineOutputElement(name="AUCregular", element = list(data = data.frame(time = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 504)), output = "AUC")) 
  
  # Doses
  weight_to_use <- if (ADC == "TDM1") 69.4 else 59
  dose_tb <- read_csv(paste0("results/", ADC, "_poppk_", analysis_date, ".csv"), show_col_types = FALSE) %>%
    select(OCC) %>%
    distinct() %>%
    mutate(OCC = factor(parse_number(OCC)),
           AMT = parse_number(as.character(OCC)) * weight_to_use * 1000,  # Calculate for a 70 kg subject and convert to ng
           AUC = "",
           Cmax = "",
           Cmin = "") %>% 
    arrange(OCC)
  n_round <- 3
  infusion_time <- 1.5
  for (i in 1:nrow(dose_tb)) {
    simparams <- read.csv(paste0("results/", ADC, "_table_simulated_parameters.csv")) %>%
      left_join(read.csv(paste0("results/", ADC, "_table_simulated_parameters_dict.csv"))) %>%
      distinct() %>%
      filter(OCC == dose_tb$OCC[i])
    write.csv(select(simparams, -OCC), file = "results/sim_tb.csv", quote = FALSE, row.names = FALSE)
    
    doses_tb <- data.frame(id = simparams$id,
                           time = 0,
                           amount = dose_tb$AMT[i], 
                           tInf = infusion_time)
    write.csv(doses_tb, file = "results/doses_tb.csv", quote = FALSE, row.names = FALSE)
    
    defineIndividualElement(name = "simulatedParameters", element = "results/sim_tb.csv")
    setGroupSize(group = "simulationGroup1", size = nrow(simparams))
    
    defineTreatmentElement(name="occ_doses", element=list(admID=1, data="results/doses_tb.csv"))
    setGroupElement(group = "simulationGroup1", elements = c("simulatedParameters", "Cmax", "AUCregular", "occ_doses"))
    setSharedIds(c("individual"))
    runSimulation()
    
    # AUC
    AUC_result <- getSimulationResults()$res$AUC
    AUC_res <- AUC_result %>%
      group_by(id) %>%
      summarise(AUC = max(AUC)) %>%
      pull(AUC) 
    dose_tb$AUC[i] <- paste0(round(quantile(AUC_res, probs = c(0.025, 0.5, 0.975)), n_round), collapse = " ")
    
    # Cmax
    Cmax_result <- getSimulationResults()$res$Cc_wghtd
    Cmax_res <- Cmax_result %>%
      group_by(id) %>%
      summarise(Cmax = max(Cc_wghtd)) %>%
      pull(Cmax) 
    dose_tb$Cmax[i] <- paste0(round(quantile(Cmax_res, probs = c(0.025, 0.5, 0.975)), n_round), collapse = " ")
    
    # Cmin (or c-trough, last time)
    Cmin_res <- Cmax_result %>%
      group_by(id) %>%
      filter(time > 2) %>% # Ensures time = 0 is not captured, or could have just used C504 hours
      summarise(Cmin = min(Cc_wghtd)) %>%
      pull(Cmin) 
    dose_tb$Cmin[i] <- paste0(round(quantile(Cmin_res, probs = c(0.025, 0.5, 0.975)), n_round), collapse = " ")
    
    pk_tb <- left_join(Cmax_result, AUC_result) %>%
      mutate(OCC = dose_tb$OCC[i], AMT = dose_tb$AMT[i])
    if (i == 1) pk_tb_all <- pk_tb else pk_tb_all <- bind_rows(pk_tb_all, pk_tb)
    
    message(paste0(round(i/nrow(dose_tb) * 100, 3), "% complete!"))
  }
  write_rds(pk_tb_all, paste0("results/", ADC, "_pk_tb.rds"))
  pk_tb_all %>% select(id, OCC) %>% distinct() %>% nrow() # 7000 (7 doses * 1000 replicates)
  pk_tb_all %>% select(id) %>% distinct() %>% nrow() # 1000 replicates
  
  # Combine the simulated results
  sim_results <- dose_tb %>% 
    select(-AMT) %>%
    pivot_longer(-OCC) %>%
    separate(value, into = c("lower", "value", "upper"), sep = " ") %>%
    mutate_at(vars(lower, value, upper), as.numeric) %>%
    mutate(lower = if_else(name == "AUC", lower * 24/504, lower), # Change total AUC to AUC/day
           value = if_else(name == "AUC", value * 24/504, value),
           upper = if_else(name == "AUC", upper * 24/504, upper))
  sim_results2 <- sim_results %>%
    mutate(Legend = "Simulated") %>%
    select(OCC, Legend, name, value) %>%
    mutate(OCC = factor(as.character(OCC)))
  
  # Original data 
  original_results <- read_csv(paste0("results/", ADC, "_poppk_", analysis_date, ".csv"), show_col_types = FALSE) 
  nca_tb <- tibble()
  intervals <- data.frame(start = 0, end = 504, cmin = TRUE, cmax = TRUE, auclast=TRUE)
  for (id in unique(original_results$ID)) {
    nca_data <- filter(original_results, ID == id) %>%
      select(ID:AMT, DV, OCC) %>%
      mutate(DV = parse_number(if_else(DV == ".", "0", DV)),
             AMT = parse_number(if_else(AMT == ".", "0", AMT)),
             Duration = if_else(AMT != 0, 1.5, 0))
    conc <- PKNCAconc(nca_data, DV ~ TIME | ID)
    dose <- PKNCAdose(nca_data, AMT ~ TIME | ID, route = "intravascular", duration = "Duration")
    pkdata <- PKNCAdata(conc, dose, intervals = intervals)
    results <- pk.nca(pkdata)
    summary_results <- tibble(summary(results)) %>%
      mutate(OCC = parse_number(unique(nca_data$OCC))) %>%
      select(OCC, Cmin = cmin, Cmax = cmax, AUC = auclast)
    nca_tb <- bind_rows(nca_tb, summary_results)
  }
  nca_tb <- suppressWarnings(mutate_at(nca_tb, vars(Cmin:AUC), parse_number)) %>%
    mutate(OCC = as.character(OCC),
           AUC = AUC * 24/504, # Change total AUC to AUC/day
           Legend = "Observed 2") %>%
    group_by(OCC) %>%
    summarise(Cmin = median(Cmin, na.rm = TRUE),
              Cmax = median(Cmax, na.rm = TRUE),
              AUC = median(AUC, na.rm = TRUE),
              Legend = unique(Legend))
  
  # Original data from the reported PK metrics
  sheet <- if (ADC == "TDM1") "S2" else "S4"
  PK_dat <- read_xlsx("Supplementary_Tables.xlsx", sheet = sheet, skip = 1) %>%
    select(`mg/kg dose`, contains(c("Cmin", "Cmax", "AUC_ADC"))) %>%
    mutate(Cmin_ADC_Mean = parse_number(Cmin_ADC_Mean),
           Cmax_ADC_Mean = parse_number(Cmax_ADC_Mean),
           AUC_ADC_Mean = parse_number(AUC_ADC_Mean),
           Cmin_ADC_Mean = if_else(str_detect(Cmin_ADC_Units, "ug"), Cmin_ADC_Mean * 1000, Cmin_ADC_Mean),
           Cmax_ADC_Mean = if_else(str_detect(Cmax_ADC_Units, "ug"), Cmax_ADC_Mean * 1000, Cmax_ADC_Mean),
           AUC_ADC_Mean = if_else(str_detect(AUC_ADC_Units, "ug"), AUC_ADC_Mean * 1000, AUC_ADC_Mean),
           AUC_ADC_Mean = if_else(str_detect(AUC_ADC_Units, "h_"), AUC_ADC_Mean / 24, AUC_ADC_Mean)) %>%
    rowwise() %>%
    mutate(Cmin_ADC_N = if_else(is.na(Cmin_ADC_N), 1, Cmin_ADC_N), # Times argument required
           Cmax_ADC_N = if_else(is.na(Cmax_ADC_N), 1, Cmax_ADC_N), 
           AUC_ADC_N = if_else(is.na(AUC_ADC_N), 1, AUC_ADC_N),
           Cmin = list(rep(Cmin_ADC_Mean, Cmin_ADC_N)),
           Cmax = list(rep(Cmax_ADC_Mean, Cmax_ADC_N)),
           AUC = list(rep(AUC_ADC_Mean, AUC_ADC_N))) %>%
    select(OCC = `mg/kg dose`, Cmin:AUC) %>%
    filter(!is.na(OCC)) %>%
    group_by(OCC) %>%
    summarise(Cmin = median(c(unlist(Cmin)), na.rm = TRUE),
              Cmax = median(c(unlist(Cmax)), na.rm = TRUE),
              AUC = median(c(unlist(AUC)), na.rm = TRUE)) %>%
    mutate(Legend = "Observed 1") %>%
    mutate(OCC = as.character(OCC))
  
  # Make the diagnostics plot
  diagnostics_plot_dat <- nca_tb %>%
    bind_rows(PK_dat) %>%
    mutate(OCC = factor(OCC)) %>%
    pivot_longer(Cmin:AUC) %>%
    bind_rows(sim_results2) 
  
  diagnostics_plot <- diagnostics_plot_dat %>%
    ggplot(aes(OCC, value, color = Legend)) +
    geom_point(shape = 16, size = 1, position = position_dodge(width = 0.5)) +
    geom_errorbar(data = sim_results, aes(ymin = lower, ymax = upper), 
                  width = 0.2, color = "#619CFF", position = position_nudge(x = 0.18)) +
    facet_wrap(~name, nrow = 1) +
    scale_y_log10() +
    xlab("Dose (mg/kg)") +
    ylab(paste0(gsub("T", "T-", ADC), " concentrations (ng/mL[/day])")) +
    theme_bw() +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 12, face = "bold"),
          plot.title = element_text(hjust = 0, vjust = 0, size = 14, face = "bold"),
          strip.text = element_text(size = 12, color = "black", face = "bold"))
  diagnostics_plot <- if (ADC == "TDM1") diagnostics_plot + ggtitle(paste0("A. ", gsub("T", "T-", ADC))) else diagnostics_plot + ggtitle(paste0("B. ", gsub("T", "T-", ADC)))
  x <- if (ADC == "TDM1") 0.6 else 0.6
  graph2ppt(diagnostics_plot, paste0("results/", ADC, "_PK_metrics_simulation"), append = FALSE, width = 16 * x, height = 6 * x)
}
