# Import libraries (and initialize Monolix, which will be used to fit the model)
# ----------------
library(lixoftConnectors)
initializeLixoftConnectors(software = "monolix", path = 'C:/Program Files/Lixoft/MonolixSuite2024R1', force = TRUE)
library(tidyverse)
library(readxl)
library(export)

# Add results folder
if (!file.exists("results")) dir.create("results")

# Analysis date
analysis_date <- "06_07_2024" # From: format(Sys.Date(), '%d_%m_%Y')

# Loop through the ADCs
for (ADC in c("TDM1", "TDXd")) {
  # Exploratory data analysis
  # -------------------------
  
  # Import extracted data
  sheet <- if (ADC == "TDM1") "S1" else "S3"
  dat <- read_xlsx("Supplementary_Tables.xlsx", sheet = sheet, skip = 1) %>%
    select(STUDY_ID = `Study Acronym`, CANCER = `Cancer type`, CYCLE_DAYS = `Cycle days`,
           WT = `Weight (kg)`, MG_KG = `mg/kg dose`, N = `Sample size`, ADC_PL = `ADC or PL`,
           TIME = `Time`, DV = Concentration, UNIT_TIME = `Units for time`, UNIT_DV = `Units for concentration`)
  
  # Make time and DV units consistent 
  dat <- dat %>%
    mutate(TIME = if_else(UNIT_TIME == "d", TIME * 24, TIME), # Change all time to days
           DV = if_else(UNIT_DV == "ug_mL", DV * 1000, DV)) # Change all DV to ng/mL
  
  # Impute missing weight using median  
  median_WT <- dat %>%
    filter(!is.na(WT)) %>%
    select(STUDY_ID, MG_KG, N, WT) %>%
    distinct() %>%
    rowwise() %>%
    mutate(N = list(rep(WT, N))) %>%
    pull(N) %>%
    unlist() %>%
    median() # 69.4 for T-DM1, 59 for T-DXd This approach ensures that larger studies contribute more e.g. a median of 50kg (100 participants), 
  # 55kg (1 participant) and 60kg (1 participant) is 50kg, not 55kg.
  dat <- dat %>%
    mutate(WT = if_else(is.na(WT), median_WT, WT))
  
  # Get the occasion (combinations of frequency of dosing and dose)
  dat <- dat %>%
    mutate(CYCLE_DAYS = if_else(CYCLE_DAYS == 7, "qwk", "q3wk"),
           OCC = paste0(MG_KG, "mg/kg ", CYCLE_DAYS),
           OCN = as.numeric(factor(OCC))) %>%
    select(-UNIT_TIME, -UNIT_DV, -CYCLE_DAYS)
  
  # Some EDA
  # unique(dat$STUDY_ID)
  id_levels <- if (ADC == "TDM1") c("TDM3569g", "TDM4258g", "TDM4374g", "TDM4688g", "JO22591", "BO25499", "GATSBY", 
                                    "TH3RESA", "MARIANNE", "NCT03153163", "KAMELEON_B", "KAMELEON_P", "SHR_A1201") else sort(unique(dat$STUDY_ID))
  
  # ADC plot
  nrow1 <- if (ADC == "TDM1") 3 else 2
  nrow2 <- if (ADC == "TDM1") 4 else 1
  ADC_plot <- dat %>%
    mutate(STUDY_ID = factor(STUDY_ID, levels = id_levels)) %>%
    rename(`Dosing schedule` = OCC) %>%
    filter(ADC_PL == "ADC") %>%
    ggplot(aes(x = TIME, y = DV, group = OCN, colour = `Dosing schedule`)) + 
    geom_line() + 
    geom_point(shape = 1, size = 2) +
    facet_wrap(~ STUDY_ID, nrow = nrow1) +
    xlab("Time (hours)") + 
    ylab(paste0(gsub("T", "T-", ADC), " concentrations (ng/mL)")) +
    ggtitle(gsub("T", "T-", ADC)) +
    theme_bw()  +
    theme(axis.text = element_text(size = 12, face = "bold"), 
          axis.title = element_text(size = 14, face = "bold"),
          plot.title = element_text(hjust = 0.5, vjust = 3, size = 16, face = "bold"),
          legend.position = "inside", 
          legend.position.inside = c(0.8, 0.15),
          legend.title = element_text(size = 12, face = "bold"),
          strip.text = element_text(size = 14, color = "red", face = "bold")
    ) +
    scale_y_log10() +  
    guides(colour = guide_legend(nrow = nrow2, reverse = TRUE))
  if (ADC != "TDM1") ADC_plot <- ADC_plot + theme(legend.position = "bottom") + guides(colour = guide_legend(nrow = nrow2, reverse = FALSE))
  x <- if (ADC == "TDM1") 0.8 else 0.8
  graph2ppt(ADC_plot, paste0("results/", ADC, "_ADC"), append = FALSE, width = 16 * x, height = 9 * x)
  
  # PL plot
  nrow1 <- if (ADC == "TDM1") 3 else 2
  nrow2 <- if (ADC == "TDM1") 5 else 1
  PL_plot <- dat %>%
    mutate(STUDY_ID = factor(STUDY_ID, levels = id_levels)) %>%
    rename(`Dosing schedule` = OCC) %>%
    filter(ADC_PL == "PL") %>%
    ggplot(aes(x = TIME, y = DV, group = OCN, colour = `Dosing schedule`)) + 
    geom_line() + 
    geom_point(shape = 1, size = 2) +
    facet_wrap(~ STUDY_ID, nrow = nrow1) +
    xlab("Time (hours)") + 
    ylab(paste0(gsub("T", "", ADC), " concentrations (ng/mL)")) +
    ggtitle(gsub("T", "", ADC)) +
    theme_bw()  +
    theme(axis.text = element_text(size = 12, face = "bold"), 
          axis.title = element_text(size = 14, face = "bold"),
          plot.title = element_text(hjust = 0.5, vjust = 3, size = 16, face = "bold"),
          legend.position = "inside", 
          legend.position.inside = c(0.8, 0.15),
          legend.title = element_text(size = 12, face = "bold"),
          strip.text = element_text(size = 14, color = "red", face = "bold")
    ) +
    scale_y_log10() +  
    guides(colour = guide_legend(nrow = nrow2, reverse = TRUE))
  if (ADC != "TDM1") PL_plot <- PL_plot + theme(legend.position = "bottom") + guides(colour = guide_legend(nrow = nrow2, reverse = FALSE))
  x <- if (ADC == "TDM1") 0.8 else 0.8
  graph2ppt(PL_plot, paste0("results/", ADC, "_PL"), append = FALSE, width = 16 * x, height = 9 * x)
  
  # Get Monolix datasets
  dat_amt <- dat %>%
    group_by(STUDY_ID, OCN, ADC_PL) %>%
    slice_head(n = 1) %>%
    mutate(TIME = 0,
           DV = ".") 
  dat <- dat %>%
    mutate(DV = as.character(DV)) %>%
    bind_rows(dat_amt) %>%
    mutate(ID = as.numeric(factor(STUDY_ID))) %>%
    arrange(ID, OCN, TIME) %>%
    rename(ID2 = ID) %>%
    mutate(AMT = if_else(TIME == 0, as.character(WT * MG_KG * 1000), "."),
           INF_DRN = if_else(TIME == 0, "1.5", "."),
           RATE = if_else(TIME == 0, as.character(WT * MG_KG * 1000/1.5), "."),
           ID = paste0(OCN, ID2)) %>%
    filter(ADC_PL == "ADC") %>%
    select(ID, TIME, AMT, RATE, DV, N, OCC) 
  # Save data
  write.csv(dat, paste0("results/", ADC, "_poppk_", analysis_date, ".csv"), row.names = FALSE, quote = FALSE, na = ".")
  
  
  # Run Monolix through R
  # ---------------------
  ## Install lixoftConnectors
  # install.packages("C:/Program Files/Lixoft/MonolixSuite2024R1/connectors/lixoftConnectors.tar.gz",
  #                  repos = NULL, type="source", INSTALL_opts ="--no-multiarch")
  
  # Initiate new project
  newProject(data = list(dataFile = paste0("results/", ADC, "_poppk_", analysis_date, ".csv"),
                         headerTypes = c("id", "time", "amount", "rate", "observation", "regressor", "ignore")),
             modelFile = 'base.txt')
  
  # Check you have the right settings
  str(getObservationInformation())
  str(getContinuousObservationModel())
  str(getIndividualParameterModel())
  getIndividualParameterModel()$correlationBlocks
  
  # Change error model, distribution, variability and correlations
  setErrorModel(DV = "constant")
  setIndividualParameterDistribution(Cl = "logNormal", V1 = "logNormal", Q = "logNormal", V2 = "logNormal", prop_err = "normal")
  setIndividualParameterVariability(Cl = TRUE, V1 = TRUE, Q = FALSE, V2 = FALSE, prop_err = FALSE)
  setCorrelationBlocks(id = list(c("Cl","V1")))
  
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
  
  # Give initial estimates
  getPopulationParameterInformation()
  if (ADC == "TDM1") { 
    setPopulationParameterInformation(Cl_pop = list(initialValue = 0.03, method = "MLE"), 
                                      V1_pop = list(initialValue = 2.4, method = "MLE"),
                                      Q_pop = list(initialValue = 0.02, method = "MLE"),
                                      V2_pop = list(initialValue = 1.2, method = "MLE"),
                                      prop_err_pop = list(initialValue = 0.3, method = "MLE"),
                                      corr_V1_Cl = list(initialValue = 0, method = "MLE"),
                                      omega_Cl = list(initialValue = 1, method = "MLE"),
                                      omega_V1 = list(initialValue = 1, method = "MLE"),
                                      a  = list(initialValue = 3000, method = "MLE"))
  } else {
    setPopulationParameterInformation(Cl_pop = list(initialValue = 0.02, method = "MLE"), 
                                      V1_pop = list(initialValue = 2.8, method = "MLE"),
                                      Q_pop = list(initialValue = 0.02, method = "MLE"),
                                      V2_pop = list(initialValue = 1.2, method = "MLE"),
                                      prop_err_pop = list(initialValue = 0.3, method = "MLE"),
                                      corr_V1_Cl = list(initialValue = 0, method = "MLE"),
                                      omega_Cl = list(initialValue = 1, method = "MLE"),
                                      omega_V1 = list(initialValue = 1, method = "MLE"),
                                      a  = list(initialValue = 3000, method = "MLE"))
  }
  
  # Run the estimation
  runScenario()
  
  # Get estimates
  tabestimates <- tibble(getEstimatedStandardErrors()$stochasticApproximation) %>%
    left_join(tibble(parameter = names(getEstimatedPopulationParameters()), estimate = getEstimatedPopulationParameters())) %>%
    select(param = parameter, estimate, rses = rse)
  param_order <- c("Cl_pop", "V1_pop",  "V2_pop", "Q_pop", "omega_Cl", "omega_V1", "corr_V1_Cl", "prop_err_pop", "a")
  tabestimates <- tabestimates %>%
    mutate(estimate = if_else(param %in% c("Cl_pop", "Q_pop"), estimate * 24, estimate), # Convert L/hr to L/day here
           param = factor(param, levels = param_order)) %>% 
    arrange(param) %>%
    mutate(param = as.character(param),
           param = gsub("Cl_pop", "Clearance (L/day)", param),
           param = gsub("V1_pop", "Central volume of distribution (L)", param),
           param = gsub("V2_pop", "Peripheral volume of distribution (L)", param),
           param = gsub("Q_pop", "Intercompartment clearance (L/day)", param),
           param = gsub("omega_Cl", "BSV Clearance", param),
           param = gsub("omega_V1", "BSV Central volume of distribution", param),
           param = gsub("corr_V1_Cl", "BSV Clearance ~ BSV Central volume of distribution", param),
           param = gsub("prop_err_pop", "Weighted proportional error (%)", param),
           param = gsub("^a$", "Constant error (ng/mL)", param),
           estimate = as.character(round(estimate, 3)),
           estimate = if_else(str_detect(estimate, "\\.\\d\\d\\d$"), estimate, paste0(estimate, "0")),
           estimate = if_else(str_detect(estimate, "\\.\\d\\d\\d$"), estimate, paste0(estimate, "0")),
           rses = as.character(round(rses, 1)),
           rses = if_else(str_detect(rses, "\\.\\d$"), rses, paste0(rses, ".0")),
           estimate2 = paste0(estimate, " (", rses, "%)")) %>%
    select(Parameter = param, Estimate = estimate2)
  write_csv(tabestimates, paste0("results/", ADC, "_parameter_estimates.csv"))
  
  # Plot VPC
  vpcplot <- plotVpc(obsName = "DV",
                     settings = list(outlierDots = TRUE, grid = FALSE, obs = TRUE, empirical = TRUE,
                                     theoretical = TRUE, predInterval = TRUE, linearInterpolation = TRUE,
                                     outlierDots = TRUE, outlierAreas = FALSE, legend = TRUE,
                                     level = 95, higherPercentile = 95, useCorrpred = TRUE,
                                     xBinsSettings = list(criteria = 'leastsquare', is.fixedNbBins = TRUE, nbBins = 3),
                                     ylab = paste0("Prediction-corrected concentrations (ng/mL)"), 
                                     xlab = "Time (hours)")) +
    theme_bw()  +
    theme(axis.text = element_text(size = 12), 
          axis.title = element_text(size = 14, face = "bold"),
          plot.title = element_text(hjust = 0, vjust = 0, size = 16, face = "bold"),
          legend.position = "inside", 
          legend.position.inside = c(0.8, 0.7),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank())
  vpcplot <- if (ADC == "TDM1") vpcplot + ggtitle(paste0("A. ", gsub("T", "T-", ADC))) else vpcplot + ggtitle(paste0("B. ", gsub("T", "T-", ADC))) + theme(legend.position = "none")
  x <- if (ADC == "TDM1") 0.7 else 0.7
  graph2ppt(vpcplot, paste0("results/", ADC, "_VPC"), append = FALSE, width = 13 * x, height = 10 * x)
  
  # Save project
  saveProject(paste0("results/", ADC, "_base.mlxtran"))
}
