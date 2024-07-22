# Import library
library(tidyverse)
library(export)
library(ComplexUpset)
library(patchwork)
library(readxl)
library(RColorBrewer)

# NOTE: Do not run it. Source it through "4. CUI analysis.R"

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

# ADC upset plot
dat2 <- dat %>%
  select(STUDY_ID, contains("_N")) %>%
  mutate(across(Cmin_ADC_N:DLT_N, ~ ifelse(is.na(.x), FALSE, TRUE))) 
colnames(dat2) <- gsub("_N|_ADC", "", colnames(dat2)) 

n_studies <- nrow(dat2)
covresultspossibilities <-  c("DLT", "PFS", "OS", "ORR", "Cmin", "AUC", "AUC_inf", "Cmax")
blue_colours <- brewer.pal(5, "Blues")
red_colours <- brewer.pal(3, "Reds")
green_colours <- brewer.pal(4, "Greens")

upset_plot <- upset(dat2, 
                    covresultspossibilities, 
                    sort_sets=FALSE, # Not to arrange by frequency
                    width_ratio = 0.2,
                    height_ratio = 1.5,
                    sort_intersections_by = c('cardinality'),
                    matrix=(intersection_matrix(
                      geom = geom_point(size = 3, shape='circle'),
                      segment=geom_segment(linetype='blank'),
                      outline_color=list( active='black',inactive='grey95')) +
                        theme(axis.text.y.left = element_text(size=20))),
                    queries = list(upset_query(set = 'Cmax', color = blue_colours[5], fill = blue_colours[5]),
                                   upset_query(set = 'AUC_inf', color = blue_colours[4], fill = blue_colours[4]),
                                   upset_query(set = 'AUC', color = blue_colours[3], fill = blue_colours[3]),
                                   upset_query(set = 'Cmin', color = blue_colours[2], fill = blue_colours[2]),
                                   upset_query(set = 'ORR', color = green_colours[4], fill = green_colours[4]),
                                   upset_query(set = 'OS', color = green_colours[3], fill = green_colours[3]),
                                   upset_query(set = 'PFS', color = green_colours[2], fill = green_colours[2]),
                                   upset_query(set = 'DLT', color = red_colours[3], fill = red_colours[3])),
                    set_sizes = upset_set_size(geom = geom_bar(width = 0.9),
                                               mapping = aes(y =after_stat(count)/n_studies),) +
                      geom_text(aes(label = ifelse(str_detect(as.character(round(100*after_stat(count)/n_studies, 1)), "\\."), 
                                                   as.character(round(100*after_stat(count)/n_studies, 1)), 
                                                   paste0(as.character(round(100*after_stat(count)/n_studies, 1)), ".0"))),
                                hjust = -0.1, stat = 'count', color="white", size = 3) +
                      theme(axis.text.x = element_text(angle = 30)) +
                      scale_y_reverse(labels = scales::percent_format(),limits =c(1,0),
                                      expand =expansion(mult = c(0, 0),
                                                        add = c(0, 0))) +
                      ylab('B. Percent reporting for PK metric/outcome'), 
                    wrap = FALSE, min_size = 0,
                    name = "D. PK metrics/outcomes reporting pattern",
                    guides = "collect",
                    base_annotations=list('Intersection size'= intersection_size(
                      bar_number_threshold = 1,  # show all numbers on top of bars
                      width = 0.5,   # reduce width of the bars
                      text = list(size = 2.5),
                      text_colors=c(on_background='black', on_bar='white')) + 
                        ylab(paste0('C. Analysis units with\nPK metrics/outcomes (n = ', n_studies, ')')) +
                        scale_y_continuous(expand=expansion(mult=c(0, 0.05))) + 
                        theme(panel.grid.major=element_blank(), # hide grid lines
                              panel.grid.minor=element_blank(),
                              axis.line = element_line(colour='black') # show axis lines
                        )
                    )
)

# Cross tabulate to find highest combination
PK_metrics <- c("Cmax", "AUC_inf", "AUC", "Cmin")
Outcomes <- c("DLT", "ORR", "OS", "PFS")
outcome_tb <- tibble(Outcome = Outcomes,
                     Cmax = vector("double", length(Outcomes)),
                     AUC_inf = vector("double", length(Outcomes)),
                     AUC = vector("double", length(Outcomes)),
                     Cmin = vector("double", length(Outcomes)))
for (i in seq_along(PK_metrics)) {
  PK_metric <- PK_metrics[[i]]
  for (j in seq_along(Outcomes)) {
    ER_outcome <- Outcomes[[j]]
    outcome_tb[j, i + 1] <- table(pull(dat2[PK_metric]), pull(dat2[ER_outcome]))["TRUE", "TRUE"]
  }
}
combinations_plot <- outcome_tb %>%
  pivot_longer(-Outcome, names_to = "Names", values_to = "values") %>%
  mutate(Outcome = factor(Outcome, levels = Outcomes),
         Names = factor(Names, levels = rev(PK_metrics)),
         label = as.character(values)) %>%
  ggplot(aes(Outcome, Names, fill = values)) +
  geom_tile(colour = "White") +
  geom_text(aes(label = label), colour = "white", size = 4) +
  scale_fill_gradient(low = "#FF6666", high = "#66CC66") +
  scale_x_discrete(expand = expansion(mult = c(0, 0), add = c(-0.5, -0.5))) + # position = "top" moves labels to the top
  scale_y_discrete(expand = expansion(mult = c(0, 0), add = c(-0.5, -0.5))) +
  theme(legend.position = "none",
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_text(colour = "black", angle = 90, margin = margin(t = 1)),
        axis.text.y = element_text(colour = "black"),
        plot.title = element_text(colour = "black", face = "bold", hjust = 0.5)) +
  ggtitle("A. Number of times PK metrics are reported\nwith ER outcomes")

# Types of cancers
cancer_dat <- dat %>%
  mutate(`Cancer type` = gsub("Breast \\(.*|Breast cancer \\(.*", "Breast", `Cancer type`),
         `Cancer type` = gsub("Uterine carcinosarcoma \\(previously-treated\\)", "Uterine", `Cancer type`),
         `Cancer type` = gsub("Colorectal \\(.*|Colorectal cancer$", "Colorectal", `Cancer type`),
         `Cancer type` = if_else(str_detect(`Cancer type`, "Solid|^Other|Adenocarcinoma, Gastric, Breast|Breast or gastric or gastro-oesophageal \\(refractory to standard therapy\\)"), 
                                 "Multiple", `Cancer type`),
         `Cancer type` = if_else(str_detect(`Cancer type`, "NSCLC"), "Lung", `Cancer type`),
         `Cancer type` = if_else(str_detect(`Cancer type`, "^Gastric or gastro"), "Gastric/GEJ", `Cancer type`)) %>%
  group_by(`Cancer type`) %>%
  summarize(count = n())
cancer_type_plot <- ggplot(cancer_dat, aes(x = count, y = reorder(`Cancer type`, desc(count)))) +
  geom_bar(stat = "identity") +
  labs(x = "E. Analysis units for each cancer type", y = "Cancer Type") + # title = "Number of Cancers by Type"
  geom_text(aes(label = ifelse(count >= 10, count, ""), hjust = ifelse(count >= 10, 1, 0)), colour = "white") +
  geom_text(aes(label = ifelse(count < 10, count, ""), hjust = ifelse(count < 10, 0, 1)), colour = "black") +
  theme_minimal() +
  theme(panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.x = element_line()) 

# Reporting based on the dose levels
dat2 <- dat %>%
  select(`mg/kg dose`, contains("_N")) %>%
  mutate(across(Cmin_ADC_N:DLT_N, ~ ifelse(is.na(.x), FALSE, TRUE))) %>%
  distinct()
colnames(dat2) <- gsub("_N|_ADC", "", colnames(dat2)) 
dat2 <- dat2 %>% # Reshape the data to long format
  pivot_longer(-`mg/kg dose`) %>%
  filter(value) %>% # Filter to keep only the rows where the PK metric/outcome was reported
  mutate(`mg/kg dose` = factor(`mg/kg dose`, levels = sort(unique(dat2$`mg/kg dose`))),
         name = factor(name, levels = covresultspossibilities))
manual_fill <- c("DLT" = red_colours[3], "PFS" = green_colours[2], "OS" = green_colours[3], 
                 "ORR" = green_colours[4], "Cmin" = blue_colours[2], "AUC" = blue_colours[3], 
                 "AUC_inf" = blue_colours[4], "Cmax" = blue_colours[5])
dose_report_plot <- ggplot(dat2, aes(x = `mg/kg dose`, y = name, fill = name)) +
  geom_tile(color = "black") +
  scale_y_discrete(limits = levels(dat2$name)) +  
  scale_fill_manual(values = manual_fill) +
  theme_minimal() +
  theme(axis.text.y = element_blank(), 
        axis.text.x = element_text(size = 10),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 12, face = "bold"),
        legend.position = "none",
        panel.grid.major=element_blank(), # hide grid lines
        panel.grid.minor=element_blank()
  ) +
  xlab("F. Doses (mg/kg) for which the\nPK metrics/outcomes were reported")

# Overall plot
layout <- "
ABBC
DEEF
"

font_size <- 10
overallplot <- 
  # Plot A
  (combinations_plot +
     theme(axis.text = element_text(size = font_size),
           plot.title = element_text(size = font_size + 1, face = "bold"))) +
  # Plot B
  (upset_plot[[2]] +
     theme(axis.text.y.left = element_text(size = font_size),
           axis.ticks.y.left = element_line(),
           axis.title.y.left = element_text(size = font_size + 1, face = "bold"))) +
  # Plot C
  (cancer_type_plot +
     theme(axis.text.y.left = element_text(size = font_size),
           axis.text.x = element_text(size = font_size - 2),
           axis.title.x = element_text(size = font_size + 1, face = "bold"))) +
  # Plot D
  (upset_plot[[3]] +
     theme(axis.text.x = element_text(size = font_size - 2),
           axis.title.x = element_text(size = font_size + 1, face = "bold"))) +
  # Plot E
  (upset_plot[[4]]+
     theme(axis.text.y.left = element_text(size = font_size, colour = "black"),
           axis.title.x = element_text(size = font_size + 1, face = "bold"))) +
  # Plot F
  (dose_report_plot +
     theme(axis.text.x = element_text(size = font_size - 2),
           axis.title.x = element_text(size = font_size + 1, face = "bold"))) +
  plot_layout(design = layout)

x <- if (ADC == "TDM1") 0.8 else 0.8
graph2ppt(overallplot, paste0("results/", ADC, "_upset"), append = FALSE, width = 16 * x, height = 9 * x)

## ER analyses
# ------------
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
# print(AUC_ratio)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. # For T-DM1
# 0.95    0.97    0.98    0.98    0.99    1.00 
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. # For T-DXd
# 0.87    0.92    0.93    0.93    0.94    0.96
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
ER_outcomes <- c("ORR", "OS", "PFS", "DLT")
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
er_tb1 <- tibble(Outcome = ER_outcomes,
                 Cmin = vector("list", length(ER_outcomes)),
                 Cmax = vector("list", length(ER_outcomes)),
                 AUC = vector("list", length(ER_outcomes)))
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
      er_tb1[j, i + 1] <- tibble(x = dat2$conc,
                                 Dose = dat2$Dose,
                                 Observed = dat2$y,
                                 Predicted = logitmod$fitted.values) %>%
        list() %>%
        list()
      er_mod_tb1[j, i + 1] <- list(list(summary(logitmod)$coefficients[, 1]))
    }
  }
}

# Scenario 2 (phases I and II)
er_tb2 <- tibble(Outcome = ER_outcomes,
                 Cmin = vector("list", length(ER_outcomes)),
                 Cmax = vector("list", length(ER_outcomes)),
                 AUC = vector("list", length(ER_outcomes)))
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
    er_tb2[j, i + 1] <- tibble(x = dat2$conc,
                               Dose = dat2$Dose,
                               Observed = dat2$y,
                               Predicted = logitmod$fitted.values) %>%
      list() %>%
      list()
    er_mod_tb2[j, i + 1] <- list(list(summary(logitmod)$coefficients[, 1]))
  }
}

# Scenario 3 (all phases)
er_tb3 <- tibble(Outcome = ER_outcomes,
                 Cmin = vector("list", length(ER_outcomes)),
                 Cmax = vector("list", length(ER_outcomes)),
                 AUC = vector("list", length(ER_outcomes)))
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
    er_tb3[j, i + 1] <- tibble(x = dat2$conc,
                               Dose = dat2$Dose,
                               Observed = dat2$y,
                               Predicted = logitmod$fitted.values) %>%
      list() %>%
      list()
    er_mod_tb3[j, i + 1] <- list(list(summary(logitmod)$coefficients[, 1]))
  }
}

# Make plots
# ----------
# Scenario 1 (only phase I)
for (i in seq_along(PK_metrics)) {
  PK_metric <- sub("_ADC", "", PK_metrics[[i]])
  for (j in seq_along(ER_outcomes)) {
    if(is.null(er_tb1[[j, i + 1]][[1]])) next
    er_plot <- er_tb1[[j, i + 1]][[1]] %>%
      mutate(Outcome = ER_outcomes[[j]])
    if (j == 1) er_plot_dat <- er_plot else er_plot_dat <- bind_rows(er_plot_dat, er_plot)
  }
  if (i == 1) plot_dat <- mutate(er_plot_dat, PK_metric = PK_metric) else plot_dat <- bind_rows(plot_dat, mutate(er_plot_dat, PK_metric = PK_metric))
}
scenario_1_data <- plot_dat %>%
  mutate(Scenario = paste0("Scenario 1 (only phase I, n = ", nrow(filter(dat, `Study Phase` %in% c("I"))), ")"))
# Scenario 2 (phases I and II)
for (i in seq_along(PK_metrics)) {
  PK_metric <- sub("_ADC", "", PK_metrics[[i]])
  for (j in seq_along(ER_outcomes)) {
    er_plot <- er_tb2[[j, i + 1]][[1]] %>%
      mutate(Outcome = ER_outcomes[[j]])
    if (j == 1) er_plot_dat <- er_plot else er_plot_dat <- bind_rows(er_plot_dat, er_plot)
  }
  if (i == 1) plot_dat <- mutate(er_plot_dat, PK_metric = PK_metric) else plot_dat <- bind_rows(plot_dat, mutate(er_plot_dat, PK_metric = PK_metric))
}
scenario_2_data <- plot_dat %>%
  mutate(Scenario = paste0("Scenario 2 (phases I and II, n = ", nrow(filter(dat, `Study Phase` %in% c("I", "II"))), ")"))
# Scenario 3 (all phases)
for (i in seq_along(PK_metrics)) {
  PK_metric <- sub("_ADC", "", PK_metrics[[i]])
  for (j in seq_along(ER_outcomes)) {
    er_plot <- er_tb3[[j, i + 1]][[1]] %>%
      mutate(Outcome = ER_outcomes[[j]])
    if (j == 1) er_plot_dat <- er_plot else er_plot_dat <- bind_rows(er_plot_dat, er_plot)
  }
  if (i == 1) plot_dat <- mutate(er_plot_dat, PK_metric = PK_metric) else plot_dat <- bind_rows(plot_dat, mutate(er_plot_dat, PK_metric = PK_metric))
}
scenario_3_data <- plot_dat %>%
  mutate(Scenario = paste0("Scenario 3 (all phases, n = ", nrow(dat), ")"))

# Both scenarios
er_plot <- scenario_1_data %>%
  bind_rows(scenario_2_data) %>%
  bind_rows(scenario_3_data) %>%
  ggplot(aes(x, Predicted, color = Outcome)) + # , lty = Outcome
  geom_line(linewidth = 1) +
  geom_point(aes(x, Observed, color = Outcome, shape = Outcome)) +
  scale_shape_manual(values = seq(0, length(ER_outcomes))) +
  facet_wrap(Scenario ~ PK_metric, nrow = 3, scales = "free_x") +
  xlab(paste0(gsub("T", "T-", ADC), " concentrations (ng/mL[/hr])")) + 
  ylab("Predicted outcome probability") +
  ggtitle(paste0(ifelse(ADC == "TDM1", "A. ", "B. "), gsub("T", "T-", ADC))) +
  theme_bw() +
  theme(axis.text = element_text(size = 10), 
        axis.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(hjust = 0, vjust = 0, size = 14, face = "bold"),
        strip.text = element_text(size = 12, color = "black", face = "bold"),
        legend.title = element_text(size = 10, face = "bold"),
        legend.text = element_text(size = 10)
  ) 
x <- if (ADC == "TDM1") 0.8 else 0.8
graph2ppt(er_plot, paste0("results/", ADC, "_ER_EDA"), append = FALSE, width = 16 * x, height = 12 * x)
