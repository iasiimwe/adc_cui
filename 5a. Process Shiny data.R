# Import library
library(tidyverse)
library(readxl)

# Define ADC
ADCs <- c("TDM1", "TDXd")
for (z in seq_along(ADCs)) {
  ADC <- ADCs[[z]]
  # Import extracted data
  sheet <- if (ADC == "TDM1") "S2" else "S4"
  dat <- read_xlsx("Supplementary_Tables.xlsx", sheet = sheet, skip = 1) %>%
    select(`Study Phase`:Drug_delay_E) %>%
    rename(STUDY_ID = `Study Acronym`, FU = `Follow-up time`)
  
  # Get a combined DLT outcome (worst of DLT, Dose_reduction, Drug_discontinuation, Drug_delay)
  DLT_fn <- function (x, y) if(length(unique(x)) == 1 & TRUE %in% is.na(unique(x))) return(NA) else return(y[which.max(x)])
  dat_DLT <- dat %>%
    select(contains(c("DLT", "Dose_reduction", "Drug_discontinuation", "Drug_delay"))) %>%
    rowwise() %>%
    mutate(DLTN = DLT_fn(c(DLT_E, Dose_reduction_E, Drug_discontinuation_E, Drug_delay_E), c(DLT_N, Dose_reduction_N, Drug_discontinuation_N, Drug_delay_N)),
           DLTE = DLT_fn(c(DLT_E, Dose_reduction_E, Drug_discontinuation_E, Drug_delay_E), c(DLT_E, Dose_reduction_E, Drug_discontinuation_E, Drug_delay_E))) %>% 
    select(DLT_N = DLTN, DLT_E = DLTE)
  dat <- dat %>%
    select(!contains(c("DLT", "Dose_reduction", "Drug_discontinuation", "Drug_delay"))) %>%
    bind_cols(dat_DLT)
  
  # Remove rows that don't have any outcome data
  with_outcomes <- dat %>%
    select(ORR_N:DLT_E)
  rows_with_only_na <- apply(with_outcomes, 1, function(x) all(is.na(x)))
  row_indices <- which(rows_with_only_na)
  dat <- dat[!rows_with_only_na, ]
  length(unique(dat$STUDY_ID)) # 44 studies
  length(dat$STUDY_ID) # 55 studies (if a study has two arms, it is considered two different studies )
  
  # Inverse logit function
  inv_logit_fn <- function(x) exp(x)/(1+exp(x))
  
  # Find AUC to AUC_inf ratio for studies that reported both AUC and AUC_inf
  AUC_ratio <- dat %>%
    select(STUDY_ID, AUC_ADC_Mean, AUC_inf_ADC_Mean) %>%
    mutate_at(vars(-STUDY_ID), parse_number) %>%
    mutate(AUC_ratio = AUC_ADC_Mean/AUC_inf_ADC_Mean) %>%
    filter(!is.na(AUC_ratio)) %>%
    pull(AUC_ratio) %>%
    summary() %>%
    round(2)
  AUC_ratio
  # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. # For T-DM1
  # 0.95    0.97    0.98    0.98    0.99    1.00 
  AUC_ratio <- AUC_ratio["Median"]
  
  # For studies that reported AUC_inf but not AUC, estimated AUC from AUC_inf using AUC = AUC_inf * AUC_ratio
  dat <- dat %>%
    mutate(AUC_ADC_N = if_else(is.na(AUC_ADC_N) & !is.na(AUC_inf_ADC_N),
                               AUC_inf_ADC_N, AUC_ADC_N),
           AUC_ADC_Units = if_else(is.na(AUC_ADC_Units) & !is.na(AUC_inf_ADC_Units),
                                   AUC_inf_ADC_Units, AUC_ADC_Units),
           AUC_ADC_Mean = if_else(is.na(AUC_ADC_Mean) & !is.na(AUC_inf_ADC_Mean),
                                  as.character(parse_number(AUC_inf_ADC_Mean) * AUC_ratio), AUC_ADC_Mean))
  
  # For all other missing exposures, impute exposure based on that seen in other dose levels
  PK_metrics <- c("Cmin_ADC", "Cmax_ADC", "AUC_ADC")
  ER_outcomes <- c("ORR", "PFS", "DLT")
  mean_PK_tb <- tibble(Dose = sort(unique(dat$`mg/kg dose`)))
  for (i in seq_along(PK_metrics)) {
    dat2 <- dat %>%
      select(contains(PK_metrics[[i]]), `mg/kg dose`) %>%
      distinct()
    colnames(dat2) <- c("PK_N", "Mean", "Units", "Dose")
    dat2 <- dat2 %>%
      mutate(Mean = parse_number(Mean),
             Mean = if_else(str_detect(Units, "ug"), Mean * 1000, Mean), # Change all units to ng
             Mean = if_else(str_detect(Units, "d_"), Mean * 24, Mean)) %>% # Change all units to hour (* 1 day = * 24 hours)
      filter(!is.na(Mean)) %>% # Remove studies without PK data
      rowwise() %>%
      mutate(PK = list(rep(Mean, PK_N))) %>%
      group_by(Dose) %>%
      summarize(PK = median(c(unlist(PK))), .groups = "drop") # This ensures sample size is accounted for
    colnames(dat2)[2] <- paste0(PK_metrics[[i]], "_mean")
    mean_PK_tb <- mean_PK_tb %>% left_join(dat2)
  }
  dat <- dat %>%
    left_join(rename(mean_PK_tb, `mg/kg dose` = Dose)) 
  
  # Load data saved from Monolix simulations
  pk_tb_data <- tibble(read_rds(paste0("results/", ADC, "_pk_tb.rds")))
  
  # Get outcome equations
  # Scenario 1 (only phase I)
  er_mod_tb1 <- tibble(Outcome = ER_outcomes,
                       Cmin = vector("list", length(ER_outcomes)),
                       Cmax = vector("list", length(ER_outcomes)),
                       AUC = vector("list", length(ER_outcomes)))
  for (i in seq_along(PK_metrics)) {
    PK_metric <- PK_metrics[[i]]
    for (j in seq_along(ER_outcomes)) {
      ER_outcome <- ER_outcomes[[j]]
      dat2 <- dat %>% filter(`Study Phase` %in% c("I"))
      colnames(dat2)[str_detect(colnames(dat2), "months")] <- gsub("_months", "months", colnames(dat2)[str_detect(colnames(dat2), "months")])
      dat2 <- dat2 %>%
        mutate(FU = if_else(str_detect(FU, "m"), parse_number(FU) * 365.25/12, parse_number(FU))) %>%
        select(STUDY_ID, "mg/kg dose", contains(c(PK_metric, ER_outcome)))
      colnames(dat2) <- c("ID", "Dose", "PK_N", "Mean", "Units", "PK_mean", "N", "E")
      dat2 <- dat2 %>%
        mutate(Mean = parse_number(Mean),
               Mean = if_else(str_detect(Units, "ug"), Mean * 1000, Mean), # Change all units to ng
               Mean = if_else(str_detect(Units, "d_"), Mean * 24, Mean), # Change all units to hour (* 1 day = * 24 hours)
               Mean = if_else(is.na(Mean), PK_mean, Mean)) %>% # Impute the PK data
        filter(!is.na(Mean)) %>% # Remove studies without PK data (none should go)
        filter(!is.na(N)) %>% # Remove studies without outcome data
        rename(conc = Mean) %>%
        mutate(y = if_else(E/N == 0, 1E-12, E/N)) # Arbitrary low value to avoid errors with logit function
      if (nrow(dat2 > 0)) {
        wt <- sqrt(dat2$N)
        logitmod <- glm(y ~ conc, weights = wt, data = dat2, family=quasibinomial) # 'probit' gives similar estimates
        er_mod_tb1[j, i + 1] <- list(list(summary(logitmod)$coefficients[, 1]))
      }
    }
  }
  
  # Scenario 2 (phases I and II)
  er_mod_tb2 <- tibble(Outcome = ER_outcomes,
                       Cmin = vector("list", length(ER_outcomes)),
                       Cmax = vector("list", length(ER_outcomes)),
                       AUC = vector("list", length(ER_outcomes)))
  for (i in seq_along(PK_metrics)) {
    PK_metric <- PK_metrics[[i]]
    for (j in seq_along(ER_outcomes)) {
      ER_outcome <- ER_outcomes[[j]]
      dat2 <- dat %>% filter(`Study Phase` %in% c("I", "II"))
      colnames(dat2)[str_detect(colnames(dat2), "months")] <- gsub("_months", "months", colnames(dat2)[str_detect(colnames(dat2), "months")])
      dat2 <- dat2 %>%
        mutate(FU = if_else(str_detect(FU, "m"), parse_number(FU) * 365.25/12, parse_number(FU))) %>%
        select(STUDY_ID, "mg/kg dose", contains(c(PK_metric, ER_outcome)))
      colnames(dat2) <- c("ID", "Dose", "PK_N", "Mean", "Units", "PK_mean", "N", "E")
      dat2 <- dat2 %>%
        mutate(Mean = parse_number(Mean),
               Mean = if_else(str_detect(Units, "ug"), Mean * 1000, Mean), # Change all units to ng
               Mean = if_else(str_detect(Units, "d_"), Mean * 24, Mean), # Change all units to hour (* 1 day = * 24 hours)
               Mean = if_else(is.na(Mean), PK_mean, Mean)) %>% # Impute the PK data
        filter(!is.na(Mean)) %>% # Remove studies without PK data (none should go)
        filter(!is.na(N)) %>% # Remove studies without outcome data
        rename(conc = Mean) %>%
        mutate(y = if_else(E/N == 0, 1E-12, E/N)) # Arbitrary low value to avoid errors with logit function
      wt <- sqrt(dat2$N)
      logitmod <- glm(y ~ conc, weights = wt, data = dat2, family=quasibinomial) # 'probit' gives similar estimates
      er_mod_tb2[j, i + 1] <- list(list(summary(logitmod)$coefficients[, 1]))
    }
  }
  
  # Scenario 3 (all phases)
  er_mod_tb3 <- tibble(Outcome = ER_outcomes,
                       Cmin = vector("list", length(ER_outcomes)),
                       Cmax = vector("list", length(ER_outcomes)),
                       AUC = vector("list", length(ER_outcomes)))
  for (i in seq_along(PK_metrics)) {
    PK_metric <- PK_metrics[[i]]
    for (j in seq_along(ER_outcomes)) {
      ER_outcome <- ER_outcomes[[j]]
      dat2 <- dat %>%
        mutate(FU = if_else(str_detect(FU, "m"), parse_number(FU) * 365.25/12, parse_number(FU))) %>%
        select(STUDY_ID, "mg/kg dose", contains(c(PK_metric, ER_outcome)))
      colnames(dat2) <- c("ID", "Dose", "PK_N", "Mean", "Units", "PK_mean", "N", "E")
      dat2 <- dat2 %>%
        mutate(Mean = parse_number(Mean),
               Mean = if_else(str_detect(Units, "ug"), Mean * 1000, Mean), # Change all units to ng
               Mean = if_else(str_detect(Units, "d_"), Mean * 24, Mean), # Change all units to hour (* 1 day = * 24 hours)
               Mean = if_else(is.na(Mean), PK_mean, Mean)) %>% # Impute the PK data
        filter(!is.na(Mean)) %>% # Remove studies without PK data (none should go)
        filter(!is.na(N)) %>% # Remove studies without outcome data
        rename(conc = Mean) %>%
        mutate(y = if_else(E/N == 0, 1E-12, E/N)) # Arbitrary low value to avoid errors with logit function
      # wt <- dat2$N * dat2$y # Using cases favours doses that leads to higher events, we want to favour larger studies
      wt <- sqrt(dat2$N)
      logitmod <- glm(y ~ conc, weights = wt, data = dat2, family=quasibinomial) # 'probit' gives similar estimates
      er_mod_tb3[j, i + 1] <- list(list(summary(logitmod)$coefficients[, 1]))
    }
  }
  
  if (z == 1) {
    tdm1_models <- list(er_mod_tb1, er_mod_tb2, er_mod_tb3) 
    tdm1_pk_data <- pk_tb_data
    } else {
      tdxd_models <- list(er_mod_tb1, er_mod_tb2, er_mod_tb3)
      tdxd_pk_data <- pk_tb_data
    }
}
