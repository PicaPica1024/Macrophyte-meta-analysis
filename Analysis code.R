# 01. Setup and input data -----------------------------------------------------

# 01.1. Load packages 

library(meta)
library(metafor)
library(ggplot2)
library(tidyr)
library(readxl)
library(clubSandwich)
library(dmetar)
library(MuMIn)
library(dplyr)
library(cowplot)
library(patchwork)
library(grid)
library(circlize)
library(gridExtra)
library(eulerr)
library(png)

# 01.2. Read the current analysis dataset 

project_root <- getwd()

analysis_output_dir <- file.path(project_root, "analysis outputs")

dir.create(
  analysis_output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

data_file <- file.path(project_root, "data.xlsx")
data <- read_xlsx(data_file, sheet = 1)
data <- as.data.frame(data)

# Rename two unused source columns that contain symbols not supported by some
# Windows locales. This prevents repeated encoding warnings during model-frame
# evaluation without changing any values or variables used in the analyses.
names(data)[startsWith(
  names(data),
  "Elevated.Temperature(High)"
)] <- "ElevatedTemperatureHigh"

names(data)[startsWith(
  names(data),
  "Background.Chl"
)] <- "BackgroundChlorophyll"

input_dataset_N <- length(unique(data$StudyID))
input_dataset_n <- nrow(data)

# Excel row number is retained only for tracing exclusions back to the source.
source_excel_row <- seq_len(nrow(data)) + 1L

cols_to_convert <- c(
  "C_mean", "C_sd", "C_n",
  "F1_mean", "F1_sd", "F1_n",
  "F2_mean", "F2_sd", "F2_n",
  "F1_F2_mean", "F1_F2_sd", "F1_F2_n"
)

data[cols_to_convert] <- lapply(data[cols_to_convert], as.numeric)

# Style of plots retained from the original code.
my_colours <- c("#EE7733", "#009988", "#EE3377")
my_shapes <- c(16, 17, 18)


# 01.3. Apply the prespecified analysis exclusions 

response_mean_cols <- c("C_mean", "F1_mean", "F2_mean", "F1_F2_mean")

# Non-positive response means are NOT removed here. At the user's request they
# are also passed through the multiplicative calculation so that the empirical
# comparison can show how interaction classification depends on null-model
# choice. MultiplicativePositiveResponse is retained as an audit flag because
# non-positive means violate the usual proportional-response interpretation.
data$MultiplicativePositiveResponse <- apply(
  data[response_mean_cols],
  1,
  function(x) all(is.finite(x) & x > 0)
)

# Identify invasive-species stressors from the two stressor-type columns.

factor_one_type <- tolower(
  trimws(as.character(data$FactorOneType))
)

factor_two_type <- tolower(
  trimws(as.character(data$FactorTwoType))
)

invasive_stressor <-
  factor_one_type == "invasive species" |
  factor_two_type == "invasive species"

invasive_stressor[is.na(invasive_stressor)] <- FALSE

exclusion_reason <- rep(NA_character_, nrow(data))
exclusion_reason[invasive_stressor] <-
  "Invasive-species stressor excluded from the revised analysis"

excluded <- !is.na(exclusion_reason)

# A compact audit table remains available in the R session.
exclusion_log <- data.frame(
  SourceExcelRow = source_excel_row[excluded],
  StudyID = data$StudyID[excluded],
  ObservationID = data$ObservationID[excluded],
  Title = data$Title[excluded],
  FactorPair = data$FactorPair[excluded],
  ExclusionReason = exclusion_reason[excluded],
  stringsAsFactors = FALSE
)

# Apply the exclusions after preserving their provenance in exclusion_log.
data <- data[
  !excluded,
  ,
  drop = FALSE
]

rownames(data) <- NULL

after_invasive_N <- length(unique(data$StudyID))
after_invasive_n <- nrow(data)

# 02. Calculate effect sizes under the three null models -----------------------

# The primary analysis uses the common standardization applied to all three
# null-model predictions by Morris et al. (2022). This keeps the effect-size
# calculation, small-sample correction and sampling variance identical across
# null models so that model comparisons differ only in the predicted mean.
data$Sd_p <- sqrt(
  (((data$F1_n - 1) * data$F1_sd^2) +
     ((data$F2_n - 1) * data$F2_sd^2)) /
    (data$F1_n + data$F2_n - 2)
)

data$n_p <- pmin(data$F1_n, data$F2_n)

data$j <-
  1 - 3 / (4 * (data$F1_F2_n + data$n_p - 2) - 1)

data$s <- sqrt(
  (((data$F1_F2_n - 1) * data$F1_F2_sd^2) +
     ((data$n_p - 1) * data$Sd_p^2)) /
    (data$F1_F2_n + data$n_p - 2)
)

# Calculate predicted means under the three null models. Following Morris et
# al. (2022), negative additive predictions are truncated at zero only for
# biomass responses, which cannot be negative. Relative growth rates retain
# negative predictions.
data <- data %>%
  mutate(
    AdditiveNullModelUnbounded = F1_mean + F2_mean - C_mean,
    AdditivePredictionFloored =
      grepl("biomass", tolower(as.character(ResponseVariables))) &
      AdditiveNullModelUnbounded < 0,
    AdditiveNullModel = ifelse(
      AdditivePredictionFloored,
      0,
      AdditiveNullModelUnbounded
    )
  )

data <- data %>%
  mutate(MultiplicativeNullModel = (F1_mean * F2_mean) / C_mean)

data <- data %>%
  mutate(
    DominanceNullModel = ifelse(
      abs(F1_mean - C_mean) - abs(F2_mean - C_mean) > 0,
      F1_mean,
      F2_mean
    )
  )

# Calculate small-sample-corrected Hedges' d using the same pooled standard
# deviation and variance definition for all three null models.
data <- data %>%
  mutate(
    HedgesdAdditive =
      ((F1_F2_mean - AdditiveNullModel) / s) * j
  )

data <- data %>%
  mutate(
    HedgesdAdditive = ifelse(
      AdditiveNullModel - C_mean < 0,
      -HedgesdAdditive,
      HedgesdAdditive
    )
  )

# Calculate Hedges' d for the multiplicative null model.
data <- data %>%
  mutate(
    HedgesdMultiplicative =
      ((F1_F2_mean - MultiplicativeNullModel) / s) * j
  )

data <- data %>%
  mutate(
    HedgesdMultiplicative = ifelse(
      MultiplicativeNullModel - C_mean < 0,
      -HedgesdMultiplicative,
      HedgesdMultiplicative
    )
  )

# Calculate Hedges' d for the dominance null model.
data <- data %>%
  mutate(HedgesdDominance = ((F1_F2_mean - DominanceNullModel) / s) * j)

data <- data %>%
  mutate(
    HedgesdDominance = ifelse(
      DominanceNullModel - C_mean < 0,
      -HedgesdDominance,
      HedgesdDominance
    )
  )

data$VdAdditive <-
  ((data$F1_F2_n + data$n_p) / (data$F1_F2_n * data$n_p)) +
  (((data$HedgesdAdditive)^2) /
     (2 * (data$F1_F2_n + data$n_p)))

data$VdMultiplicative <-
  ((data$F1_F2_n + data$n_p) / (data$F1_F2_n * data$n_p)) +
  (((data$HedgesdMultiplicative)^2) /
     (2 * (data$F1_F2_n + data$n_p)))

data$VdDominance <-
  ((data$F1_F2_n + data$n_p) / (data$F1_F2_n * data$n_p)) +
  (((data$HedgesdDominance)^2) /
     (2 * (data$F1_F2_n + data$n_p)))

# Calculate 95% confidence intervals for Hedges' d.
data$CIllAdditive <- data$HedgesdAdditive -
  1.96 * sqrt(data$VdAdditive)
data$CIULAdditive <- data$HedgesdAdditive +
  1.96 * sqrt(data$VdAdditive)

data$CIllMultiplicative <- data$HedgesdMultiplicative -
  1.96 * sqrt(data$VdMultiplicative)
data$CIULMultiplicative <- data$HedgesdMultiplicative +
  1.96 * sqrt(data$VdMultiplicative)

data$CIllDominance <- data$HedgesdDominance -
  1.96 * sqrt(data$VdDominance)
data$CIULDominance <- data$HedgesdDominance +
  1.96 * sqrt(data$VdDominance)

# Preserve model-specific eligibility instead of deleting non-positive
# responses before the null models are compared. The original absolute-effect
# threshold (|d| < 30) and the requirement for finite positive sampling
# variances are retained separately for each model.
data$AdditiveEligible <-
  is.finite(data$HedgesdAdditive) &
  is.finite(data$VdAdditive) &
  data$VdAdditive > 0 &
  abs(data$HedgesdAdditive) < 30

data$MultiplicativeEligible <-
  is.finite(data$HedgesdMultiplicative) &
  is.finite(data$VdMultiplicative) &
  data$VdMultiplicative > 0 &
  abs(data$HedgesdMultiplicative) < 30

data$DominanceEligible <-
  is.finite(data$HedgesdDominance) &
  is.finite(data$VdDominance) &
  data$VdDominance > 0 &
  abs(data$HedgesdDominance) < 30

# Assign interaction types under each null model.
data <- data %>%
  mutate(
    AdditiveInteraction = case_when(
      CIllAdditive <= 0 & CIULAdditive >= 0 ~ "Additive",
      CIllAdditive > 0 ~ "Positive non-additive",
      CIULAdditive < 0 ~ "Negative non-additive"
    )
  )

data <- data %>%
  mutate(
    MultiplicativeInteraction = case_when(
      CIllMultiplicative <= 0 & CIULMultiplicative >= 0 ~ "Additive",
      CIllMultiplicative > 0 ~ "Positive non-additive",
      CIULMultiplicative < 0 ~ "Negative non-additive"
    )
  )

data <- data %>%
  mutate(
    DominanceInteraction = case_when(
      CIllDominance <= 0 & CIULDominance >= 0 ~ "Additive",
      CIllDominance > 0 ~ "Positive non-additive",
      CIULDominance < 0 ~ "Negative non-additive"
    )
  )

datafull <- data

datafull <- subset(datafull, FactorPair != "not global change factor")

datafull$PooledSampleSize <-
  datafull$C_n + datafull$F1_n + datafull$F2_n + datafull$F1_F2_n

additive_floor_audit <- datafull[
  !is.na(datafull$AdditivePredictionFloored) &
    datafull$AdditivePredictionFloored,
  c(
    "StudyID", "ObservationID", "Title", "ResponseVariables",
    "C_mean", "F1_mean", "F2_mean", "F1_F2_mean",
    "AdditiveNullModelUnbounded", "AdditiveNullModel"
  ),
  drop = FALSE
]

three_model_data <- datafull[
  datafull$AdditiveEligible &
    datafull$MultiplicativeEligible &
    datafull$DominanceEligible,
  ,
  drop = FALSE
]

multiplicative_assumption_audit <- datafull[
  datafull$MultiplicativeEligible &
    !datafull$MultiplicativePositiveResponse,
  c(
    "StudyID", "ObservationID", "Title", "ResponseVariables",
    response_mean_cols, "HedgesdMultiplicative", "VdMultiplicative"
  ),
  drop = FALSE
]


# 02.1. Compare interaction classifications among null models 

library(eulerr)
library(grid)
library(gridExtra)

plot_venn <- function(type_label) {
  
  A <- three_model_data$AdditiveInteraction == type_label
  M <- three_model_data$MultiplicativeInteraction == type_label
  D <- three_model_data$DominanceInteraction == type_label
  
  fit <- euler(list(
    "Additive null model" = which(A),
    "Multiplicative null model" = which(M),
    "Dominance null model" = which(D)
  ))
  
  p <- plot(
    fit,
    fills = list(
      fill = c("#FF9999", "#99CCFF", "#99FF99"),
      alpha = 0.6
    ),
    labels = FALSE,
    edges = TRUE,
    quantities = list(cex = 1.05),
    main = paste(type_label, "results"),
    main.cex = 1
  )
  
  return(p)
}

p1 <- plot_venn("Positive non-additive")
p2 <- plot_venn("Additive")
p3 <- plot_venn("Negative non-additive")

venn_panels <- gridExtra::arrangeGrob(
  grobs = list(p1, p2, p3),
  ncol = 3
)

venn_model_legend <- grid::grobTree(
  grid::rectGrob(
    x = 0.15, y = 0.52, width = 0.026, height = 0.34,
    gp = grid::gpar(fill = "#FF9999", col = "#555555")
  ),
  grid::textGrob(
    "Additive null model", x = 0.17, y = 0.52, just = "left",
    gp = grid::gpar(cex = 0.95)
  ),
  grid::rectGrob(
    x = 0.42, y = 0.52, width = 0.026, height = 0.34,
    gp = grid::gpar(fill = "#99CCFF", col = "#555555")
  ),
  grid::textGrob(
    "Multiplicative null model", x = 0.44, y = 0.52, just = "left",
    gp = grid::gpar(cex = 0.95)
  ),
  grid::rectGrob(
    x = 0.73, y = 0.52, width = 0.026, height = 0.34,
    gp = grid::gpar(fill = "#99FF99", col = "#555555")
  ),
  grid::textGrob(
    "Dominance null model", x = 0.75, y = 0.52, just = "left",
    gp = grid::gpar(cex = 0.95)
  )
)

venn_figure <- gridExtra::arrangeGrob(
  venn_panels,
  venn_model_legend,
  ncol = 1,
  heights = c(11, 1)
)

grid::grid.newpage()
grid::grid.draw(venn_figure)

# 03. Interaction effects of individual observations and overall effects 


# 03.1. Additive null model 
additive_source <- datafull[
  datafull$AdditiveEligible,
  ,
  drop = FALSE
]

additive_data <- data.frame(
  StudyID = additive_source$StudyID,
  ObservationID = additive_source$ObservationID,
  Hedgesd = additive_source$HedgesdAdditive,
  CIll = additive_source$CIllAdditive,
  CIUL = additive_source$CIULAdditive,
  PooledSampleSize = additive_source$PooledSampleSize,
  model_type = "Additive null model"
)

# Significant positive effects: blue
# Significant negative effects: red
# Non-significant effects: grey
additive_data$color <- ifelse(
  additive_data$CIll > 0,
  "#027EB6",
  ifelse(
    additive_data$CIUL < 0,
    "#EE0220",
    "lightgrey"
  )
)

# Arrange observations according to Hedges' d
min_y_offset <- 2

additive_data$y_position <- min_y_offset +
  as.integer(
    reorder(
      additive_data$ObservationID,
      additive_data$Hedgesd
    )
  )

# Sampling variance calculated in the previous section
additive_data$V <- additive_source$VdAdditive

# Construct the provisional within-study covariance matrix
V_additive <- clubSandwich::impute_covariance_matrix(
  vi = additive_data$V,
  cluster = additive_data$StudyID,
  r = 0.6,
  smooth_vi = FALSE,
  check_PD = TRUE,
  return_list = FALSE
)

# Overall additive-model interaction effect
additive_model <- metafor::rma.mv(
  yi = Hedgesd,
  V = V_additive,
  random = ~ 1 | StudyID / ObservationID,
  data = additive_data,
  test = "t",
  method = "REML"
)

# Study-clustered CR2 inference is used for the reported overall result.
additive_overall_CR2 <- as.data.frame(
  clubSandwich::coef_test(
    additive_model,
    vcov = "CR2",
    cluster = additive_data$StudyID,
    test = "Satterthwaite"
  )
)

additive_overall_critical <- stats::qt(
  0.975,
  df = additive_overall_CR2$df_Satt[1]
)

additive_overall_effect_size <- additive_overall_CR2$beta[1]
additive_ci_lb <-
  additive_overall_effect_size -
  additive_overall_critical * additive_overall_CR2$SE[1]
additive_ci_ub <-
  additive_overall_effect_size +
  additive_overall_critical * additive_overall_CR2$SE[1]
additive_p_value <- additive_overall_CR2$p_Satt[1]

additive_I2 <- dmetar::mlm.variance.distribution(
  additive_model
)

additive_row_weights <- stats::weights(
  additive_model,
  type = "rowsum"
)

additive_study_weights <- data.frame(
  StudyID = names(
    tapply(
      additive_row_weights,
      additive_data$StudyID,
      sum
    )
  ),
  Weight = as.numeric(
    tapply(
      additive_row_weights,
      additive_data$StudyID,
      sum
    )
  ),
  NullModel = "Additive",
  row.names = NULL
)


# 03.2. Multiplicative null model 

multiplicative_source <- datafull[
  datafull$MultiplicativeEligible,
  ,
  drop = FALSE
]

multiplicative_data <- data.frame(
  StudyID = multiplicative_source$StudyID,
  ObservationID = multiplicative_source$ObservationID,
  Hedgesd = multiplicative_source$HedgesdMultiplicative,
  CIll = multiplicative_source$CIllMultiplicative,
  CIUL = multiplicative_source$CIULMultiplicative,
  PooledSampleSize = multiplicative_source$PooledSampleSize,
  model_type = "Multiplicative null model"
)

multiplicative_data$color <- ifelse(
  multiplicative_data$CIll > 0,
  "#027EB6",
  ifelse(
    multiplicative_data$CIUL < 0,
    "#EE0220",
    "lightgrey"
  )
)

multiplicative_data$y_position <- min_y_offset +
  as.integer(
    reorder(
      multiplicative_data$ObservationID,
      multiplicative_data$Hedgesd
    )
  )

multiplicative_data$V <- multiplicative_source$VdMultiplicative

V_multiplicative <- clubSandwich::impute_covariance_matrix(
  vi = multiplicative_data$V,
  cluster = multiplicative_data$StudyID,
  r = 0.6,
  smooth_vi = FALSE,
  check_PD = TRUE,
  return_list = FALSE
)

multiplicative_model <- metafor::rma.mv(
  yi = Hedgesd,
  V = V_multiplicative,
  random = ~ 1 | StudyID / ObservationID,
  data = multiplicative_data,
  test = "t",
  method = "REML"
)

multiplicative_overall_CR2 <- as.data.frame(
  clubSandwich::coef_test(
    multiplicative_model,
    vcov = "CR2",
    cluster = multiplicative_data$StudyID,
    test = "Satterthwaite"
  )
)

multiplicative_overall_critical <- stats::qt(
  0.975,
  df = multiplicative_overall_CR2$df_Satt[1]
)

multiplicative_overall_effect_size <-
  multiplicative_overall_CR2$beta[1]
multiplicative_ci_lb <-
  multiplicative_overall_effect_size -
  multiplicative_overall_critical *
    multiplicative_overall_CR2$SE[1]
multiplicative_ci_ub <-
  multiplicative_overall_effect_size +
  multiplicative_overall_critical *
    multiplicative_overall_CR2$SE[1]
multiplicative_p_value <- multiplicative_overall_CR2$p_Satt[1]

multiplicative_I2 <- dmetar::mlm.variance.distribution(
  multiplicative_model
)

multiplicative_row_weights <- stats::weights(
  multiplicative_model,
  type = "rowsum"
)

multiplicative_study_weights <- data.frame(
  StudyID = names(
    tapply(
      multiplicative_row_weights,
      multiplicative_data$StudyID,
      sum
    )
  ),
  Weight = as.numeric(
    tapply(
      multiplicative_row_weights,
      multiplicative_data$StudyID,
      sum
    )
  ),
  NullModel = "Multiplicative",
  row.names = NULL
)


# 03.3. Dominance null model 

dominance_source <- datafull[
  datafull$DominanceEligible,
  ,
  drop = FALSE
]

dominance_data <- data.frame(
  StudyID = dominance_source$StudyID,
  ObservationID = dominance_source$ObservationID,
  Hedgesd = dominance_source$HedgesdDominance,
  CIll = dominance_source$CIllDominance,
  CIUL = dominance_source$CIULDominance,
  PooledSampleSize = dominance_source$PooledSampleSize,
  model_type = "Dominance null model"
)

dominance_data$color <- ifelse(
  dominance_data$CIll > 0,
  "#027EB6",
  ifelse(
    dominance_data$CIUL < 0,
    "#EE0220",
    "lightgrey"
  )
)

dominance_data$y_position <- min_y_offset +
  as.integer(
    reorder(
      dominance_data$ObservationID,
      dominance_data$Hedgesd
    )
  )

dominance_data$V <- dominance_source$VdDominance

V_dominance <- clubSandwich::impute_covariance_matrix(
  vi = dominance_data$V,
  cluster = dominance_data$StudyID,
  r = 0.6,
  smooth_vi = FALSE,
  check_PD = TRUE,
  return_list = FALSE
)

dominance_model <- metafor::rma.mv(
  yi = Hedgesd,
  V = V_dominance,
  random = ~ 1 | StudyID / ObservationID,
  data = dominance_data,
  test = "t",
  method = "REML"
)

dominance_overall_CR2 <- as.data.frame(
  clubSandwich::coef_test(
    dominance_model,
    vcov = "CR2",
    cluster = dominance_data$StudyID,
    test = "Satterthwaite"
  )
)

dominance_overall_critical <- stats::qt(
  0.975,
  df = dominance_overall_CR2$df_Satt[1]
)

dominance_overall_effect_size <- dominance_overall_CR2$beta[1]
dominance_ci_lb <-
  dominance_overall_effect_size -
  dominance_overall_critical * dominance_overall_CR2$SE[1]
dominance_ci_ub <-
  dominance_overall_effect_size +
  dominance_overall_critical * dominance_overall_CR2$SE[1]
dominance_p_value <- dominance_overall_CR2$p_Satt[1]

dominance_I2 <- dmetar::mlm.variance.distribution(
  dominance_model
)

dominance_row_weights <- stats::weights(
  dominance_model,
  type = "rowsum"
)

dominance_study_weights <- data.frame(
  StudyID = names(
    tapply(
      dominance_row_weights,
      dominance_data$StudyID,
      sum
    )
  ),
  Weight = as.numeric(
    tapply(
      dominance_row_weights,
      dominance_data$StudyID,
      sum
    )
  ),
  NullModel = "Dominance",
  row.names = NULL
)

study_weight_results <- dplyr::bind_rows(
  additive_study_weights,
  multiplicative_study_weights,
  dominance_study_weights
)


# 03.4. Combine the results

combined_data <- rbind(
  additive_data,
  multiplicative_data,
  dominance_data
)

combined_data$model_type <- factor(
  combined_data$model_type,
  levels = c(
    "Additive null model",
    "Multiplicative null model",
    "Dominance null model"
  )
)

overall_effects <- data.frame(
  model_type = c(
    "Additive null model",
    "Multiplicative null model",
    "Dominance null model"
  ),
  overall_effect_size = c(
    additive_overall_effect_size,
    multiplicative_overall_effect_size,
    dominance_overall_effect_size
  ),
  ci.lb = c(
    additive_ci_lb,
    multiplicative_ci_lb,
    dominance_ci_lb
  ),
  ci.ub = c(
    additive_ci_ub,
    multiplicative_ci_ub,
    dominance_ci_ub
  ),
  p_value = c(
    additive_p_value,
    multiplicative_p_value,
    dominance_p_value
  )
)

overall_effects$model_type <- factor(
  overall_effects$model_type,
  levels = c(
    "Additive null model",
    "Multiplicative null model",
    "Dominance null model"
  )
)


# 03.5. Plot individual and overall interaction effects 

# The horizontal observation order follows the user-approved layout: smallest
# Hedges' d on the left and largest on the right. The x-axis labels are hidden;
# the vertical axis contains Hedges' d and every confidence interval is shown.
combined_data <- combined_data %>%
  dplyr::group_by(model_type) %>%
  dplyr::arrange(Hedgesd, .by_group = TRUE) %>%
  dplyr::mutate(
    PlotPosition = dplyr::row_number(),
    InteractionClass = dplyr::case_when(
      CIll > 0 ~ "Positive non-additive",
      CIUL < 0 ~ "Negative non-additive",
      TRUE ~ "Additive"
    )
  ) %>%
  dplyr::ungroup()

overall_effects <- overall_effects %>%
  dplyr::left_join(
    combined_data %>%
      dplyr::count(model_type, name = "ObservationCount"),
    by = "model_type"
  ) %>%
  dplyr::mutate(
    PlotPosition = (ObservationCount + 1) / 2,
    ResultLabel = sprintf(
      "d = %.3f\np = %.3f",
      overall_effect_size,
      p_value
    )
  )

combined_plot <- ggplot2::ggplot() +
  ggplot2::geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "#555555",
    linewidth = 0.45
  ) +
  ggplot2::geom_errorbar(
    data = combined_data,
    ggplot2::aes(
      x = PlotPosition,
      ymin = CIll,
      ymax = CIUL,
      color = InteractionClass
    ),
    width = 0,
    linewidth = 0.35,
    show.legend = FALSE
  ) +
  ggplot2::geom_point(
    data = combined_data,
    ggplot2::aes(
      x = PlotPosition,
      y = Hedgesd,
      color = InteractionClass
    ),
    size = 0.75,
    shape = 18
  ) +
  ggplot2::geom_errorbar(
    data = overall_effects,
    ggplot2::aes(
      x = PlotPosition,
      ymin = ci.lb,
      ymax = ci.ub,
      color = "Overall mean"
    ),
    width = 0,
    linewidth = 2.2,
    show.legend = FALSE
  ) +
  ggplot2::geom_point(
    data = overall_effects,
    ggplot2::aes(
      x = PlotPosition,
      y = overall_effect_size,
      color = "Overall mean"
    ),
    size = 2.2,
    shape = 18
  ) +
  ggplot2::geom_text(
    data = overall_effects,
    ggplot2::aes(
      x = 1,
      y = Inf,
      label = ResultLabel
    ),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 1.15,
    size = 5.0,
    color = "#333333",
    lineheight = 1.25
  ) +
  ggplot2::facet_wrap(
    ~ model_type,
    nrow = 1,
    scales = "free",
    labeller = ggplot2::as_labeller(c(
      "Additive null model" = "Additive null model",
      "Multiplicative null model" = "Multiplicative null model",
      "Dominance null model" = "Dominance null model"
    ))
  ) +
  ggplot2::scale_color_manual(
    values = c(
      "Positive non-additive" = "#027EB6",
      "Negative non-additive" = "#EE0220",
      "Additive" = "lightgrey",
      "Overall mean" = "#008A25"
    ),
    breaks = c(
      "Positive non-additive",
      "Negative non-additive",
      "Additive",
      "Overall mean"
    ),
    name = NULL
  ) +
  ggplot2::scale_x_continuous(breaks = NULL) +
  ggplot2::labs(
    x = NULL,
    y = "Interaction effect (Hedges' d)"
  ) +
  ggplot2::theme_bw(base_size = 15) +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_blank(),
    axis.ticks.x = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(size = 15, face = "bold"),
    legend.position = "bottom",
    legend.text = ggplot2::element_text(size = 12),
    plot.margin = ggplot2::margin(8, 10, 8, 8)
  )

print(combined_plot)

# 04. Mean interaction effects for each stressor pair 


# Correlation assumed among effect sizes extracted from the same article.
# Sensitivity analyses for alternative correlation values will be conducted
# separately and reported in the Supplementary Information.
sampling_correlation <- 0.6


# 04.1. Define the effect-size columns for the three null models 

null_model_columns <- data.frame(
  NullModel = c(
    "Additive",
    "Multiplicative",
    "Dominance"
  ),
  EffectColumn = c(
    "HedgesdAdditive",
    "HedgesdMultiplicative",
    "HedgesdDominance"
  ),
  VarianceColumn = c(
    "VdAdditive",
    "VdMultiplicative",
    "VdDominance"
  ),
  stringsAsFactors = FALSE
)


# 04.2. Function for multilevel heterogeneity statistics 

calculate_multilevel_heterogeneity <- function(model, V_matrix) {
  
  W_matrix <- solve(V_matrix)
  X_matrix <- as.matrix(model$X)
  
  P_matrix <- W_matrix -
    W_matrix %*%
    X_matrix %*%
    solve(
      crossprod(
        X_matrix,
        W_matrix %*% X_matrix
      )
    ) %*%
    t(X_matrix) %*%
    W_matrix
  
  typical_sampling_variance <- (
    nrow(X_matrix) - ncol(X_matrix)
  ) / sum(diag(P_matrix))
  
  variance_components <- model$sigma2
  
  total_variance <- sum(variance_components) +
    typical_sampling_variance
  
  Study_I2 <- 100 *
    variance_components[1] /
    total_variance
  
  Observation_I2 <- 100 *
    variance_components[2] /
    total_variance
  
  Total_I2 <- 100 *
    sum(variance_components) /
    total_variance
  
  data.frame(
    Tau2_Study = variance_components[1],
    Tau2_Observation = variance_components[2],
    Study_I2 = Study_I2,
    Observation_I2 = Observation_I2,
    Total_I2 = Total_I2
  )
}


# 04.3. Function for fitting one null model to one stressor pair 

fit_stressor_pair <- function(
    pair_data,
    pair_name,
    null_model,
    effect_column,
    variance_column,
    rho = 0.6
) {
  
  # Select the effect size and sampling variance
  pair_data$yi <- suppressWarnings(
    as.numeric(pair_data[[effect_column]])
  )
  
  pair_data$vi <- suppressWarnings(
    as.numeric(pair_data[[variance_column]])
  )
  
  # Retain observations with calculable effect sizes and variances
  pair_data <- pair_data[
    is.finite(pair_data$yi) &
      is.finite(pair_data$vi) &
      pair_data$vi > 0 &
      !is.na(pair_data$StudyID) &
      !is.na(pair_data$ObservationID),
    ,
    drop = FALSE
  ]
  
  # Number of articles and observations
  N_articles <- length(
    unique(pair_data$StudyID)
  )
  
  n_observations <- nrow(pair_data)
  
  # Prepare an empty result if the pair cannot be analysed
  empty_result <- data.frame(
    FactorPair = pair_name,
    NullModel = null_model,
    N = N_articles,
    n = n_observations,
    Estimate = NA_real_,
    CR2_SE = NA_real_,
    CR2_df = NA_real_,
    CI_lower = NA_real_,
    CI_upper = NA_real_,
    P_value = NA_real_,
    Interaction = NA_character_,
    Tau2_Study = NA_real_,
    Tau2_Observation = NA_real_,
    Study_I2 = NA_real_,
    Observation_I2 = NA_real_,
    Total_I2 = NA_real_,
    QE = NA_real_,
    QE_P_value = NA_real_,
    Status = NA_character_,
    stringsAsFactors = FALSE
  )
  
  # At least two independent articles are required for an inferential
  # article-level meta-analysis
  if (N_articles < 2 || n_observations < 2) {
    
    empty_result$Status <-
      "Not fitted: fewer than two independent articles"
    
    return(empty_result)
  }
  
  # Define article and observation identifiers
  pair_data$StudyID <- factor(
    pair_data$StudyID
  )
  
  pair_data$ObservationWithinStudy <- factor(
    paste(
      pair_data$StudyID,
      pair_data$ObservationID,
      sep = "_"
    )
  )
  
  # Construct the sampling variance-covariance matrix
  V_pair <- tryCatch(
    clubSandwich::impute_covariance_matrix(
      vi = pair_data$vi,
      cluster = pair_data$StudyID,
      r = rho,
      smooth_vi = FALSE,
      check_PD = TRUE,
      return_list = FALSE
    ),
    error = function(e) NULL
  )
  
  if (is.null(V_pair)) {
    
    empty_result$Status <-
      "Model not fitted: covariance matrix could not be constructed"
    
    return(empty_result)
  }
  
  # Fit the multilevel random-effects model
  pair_model <- tryCatch(
    metafor::rma.mv(
      yi = yi,
      V = V_pair,
      random = list(
        ~ 1 | StudyID,
        ~ 1 | ObservationWithinStudy
      ),
      data = pair_data,
      method = "REML",
      test = "t",
      sparse = TRUE
    ),
    error = function(e) NULL
  )
  
  if (is.null(pair_model)) {
    
    empty_result$Status <-
      "Model not fitted: multilevel model failed"
    
    return(empty_result)
  }
  
  # Study-clustered CR2 inference
  CR2_result <- tryCatch(
    clubSandwich::coef_test(
      pair_model,
      vcov = "CR2",
      cluster = pair_data$StudyID,
      test = "Satterthwaite"
    ),
    error = function(e) NULL
  )
  
  if (is.null(CR2_result)) {
    
    empty_result$Estimate <-
      as.numeric(pair_model$b[1])
    
    empty_result$Tau2_Study <-
      pair_model$sigma2[1]
    
    empty_result$Tau2_Observation <-
      pair_model$sigma2[2]
    
    empty_result$QE <-
      pair_model$QE
    
    empty_result$QE_P_value <-
      pair_model$QEp
    
    empty_result$Status <-
      "Model fitted, but CR2 inference was not estimable"
    
    return(empty_result)
  }
  
  # Extract the CR2-adjusted intercept
  effect_estimate <- as.numeric(
    CR2_result$beta[1]
  )
  
  effect_SE <- as.numeric(
    CR2_result$SE[1]
  )
  
  effect_df <- as.numeric(
    CR2_result$df_Satt[1]
  )
  
  effect_p <- as.numeric(
    CR2_result$p_Satt[1]
  )
  
  critical_value <- stats::qt(
    0.975,
    df = effect_df
  )
  
  CI_lower <- effect_estimate -
    critical_value * effect_SE
  
  CI_upper <- effect_estimate +
    critical_value * effect_SE
  
  # Classify the average interaction
  Interaction <- ifelse(
    CI_lower > 0,
    "Positive non-additive",
    ifelse(
      CI_upper < 0,
      "Negative non-additive",
      "Additive"
    )
  )
  
  # Calculate multilevel heterogeneity
  heterogeneity_result <- tryCatch(
    calculate_multilevel_heterogeneity(
      model = pair_model,
      V_matrix = V_pair
    ),
    error = function(e) {
      data.frame(
        Tau2_Study = pair_model$sigma2[1],
        Tau2_Observation = pair_model$sigma2[2],
        Study_I2 = NA_real_,
        Observation_I2 = NA_real_,
        Total_I2 = NA_real_
      )
    }
  )
  
  # Return one row for this stressor pair and null model
  data.frame(
    FactorPair = pair_name,
    NullModel = null_model,
    N = N_articles,
    n = n_observations,
    Estimate = effect_estimate,
    CR2_SE = effect_SE,
    CR2_df = effect_df,
    CI_lower = CI_lower,
    CI_upper = CI_upper,
    P_value = effect_p,
    Interaction = Interaction,
    Tau2_Study = heterogeneity_result$Tau2_Study,
    Tau2_Observation =
      heterogeneity_result$Tau2_Observation,
    Study_I2 = heterogeneity_result$Study_I2,
    Observation_I2 =
      heterogeneity_result$Observation_I2,
    Total_I2 = heterogeneity_result$Total_I2,
    QE = pair_model$QE,
    QE_P_value = pair_model$QEp,
    Status = "Successfully fitted",
    stringsAsFactors = FALSE
  )
}


# 04.4. Fit all stressor pairs under all three null models 

stressor_pair_names <- sort(
  unique(
    as.character(datafull$FactorPair)
  )
)

stressor_pair_names <- stressor_pair_names[
  !is.na(stressor_pair_names) &
    nzchar(stressor_pair_names)
]

pair_model_results <- list()
result_number <- 1L

for (pair_name in stressor_pair_names) {
  
  pair_data <- datafull[
    datafull$FactorPair == pair_name,
    ,
    drop = FALSE
  ]
  
  for (
    model_number in
    seq_len(nrow(null_model_columns))
  ) {
    
    pair_model_results[[result_number]] <-
      fit_stressor_pair(
        pair_data = pair_data,
        pair_name = pair_name,
        null_model =
          null_model_columns$NullModel[
            model_number
          ],
        effect_column =
          null_model_columns$EffectColumn[
            model_number
          ],
        variance_column =
          null_model_columns$VarianceColumn[
            model_number
          ],
        rho = sampling_correlation
      )
    
    result_number <- result_number + 1L
  }
}


# Combine all stressor-pair results
stressor_pair_results <- dplyr::bind_rows(
  pair_model_results
)

stressor_pair_results$NullModel <- factor(
  stressor_pair_results$NullModel,
  levels = c(
    "Additive",
    "Multiplicative",
    "Dominance"
  )
)

stressor_pair_results <- stressor_pair_results %>%
  dplyr::arrange(
    NullModel,
    dplyr::desc(Estimate)
  )


# 04.5. Results retained for the main stressor-pair figure 

# All combinations are calculated above. Only combinations containing at least
# 10 observations and at least two independent articles are retained for the
# main figure.
stressor_pair_results_main <- stressor_pair_results %>%
  dplyr::filter(
    n >= 10,
    N >= 2,
    Status == "Successfully fitted"
  )

print(stressor_pair_results)

# 05. Moderator analyses 

# 05.1. Prepare moderator variables 

moderator_data <- three_model_data


# Clean categorical variables
clean_category <- function(x) {
  
  x <- trimws(as.character(x))
  
  x[
    is.na(x) |
      tolower(x) %in% c(
        "",
        "na",
        "n/a",
        "nan",
        "unresolved",
        "notsinglespecies",
        "not single species"
      )
  ] <- NA_character_
  
  x
}

# Convert numeric moderator variables

moderator_data$AmbientTemperature <- suppressWarnings(
  as.numeric(
    moderator_data$`Ambient.Temperature(DegreesCelsius)`
  )
)

moderator_data$StudyDurationDays <- suppressWarnings(
  as.numeric(
    moderator_data$StudyDurationDays
  )
)


# Background nutrient variables

background_N_aq <- suppressWarnings(
  as.numeric(
    moderator_data$`Background.N(aq)(mg.L-1)`
  )
)

background_N_sed <- suppressWarnings(
  as.numeric(
    moderator_data$`Background.N(sed)(mg.g-1)`
  )
)

background_P_aq <- suppressWarnings(
  as.numeric(
    moderator_data$`Background.P(aq)(mg.L-1)`
  )
)


# Use water-column nitrogen when available.
# If it is unavailable, use sediment nitrogen.
background_N <- ifelse(
  is.finite(background_N_aq),
  background_N_aq,
  background_N_sed
)


# Classify background nutrient level:
#   TN >= 1 mg L-1 : high nutrient
#   TP > 0.1 mg L-1: high nutrient when TN is unavailable
moderator_data$BackgroundNutrientLevel <- ifelse(
  is.finite(background_N),
  ifelse(
    background_N >= 1,
    "High nutrient",
    "Low nutrient"
  ),
  ifelse(
    is.finite(background_P_aq),
    ifelse(
      background_P_aq > 0.1,
      "High nutrient",
      "Low nutrient"
    ),
    NA_character_
  )
)

moderator_data$BackgroundNutrientLevel <- factor(
  moderator_data$BackgroundNutrientLevel,
  levels = c(
    "Low nutrient",
    "High nutrient"
  )
)


# Fish and snail presence

moderator_data$FishPresence <- factor(
  clean_category(
    moderator_data$`Background.Fish.Presence(Y/N)`
  ),
  levels = c("N", "Y")
)

moderator_data$SnailPresence <- factor(
  clean_category(
    moderator_data$`Background.Snails.Presence(Y/N)`
  ),
  levels = c("N", "Y")
)


# Only single-species observations are used for plant-trait analyses
single_species_observation <- (
  tolower(
    trimws(
      as.character(
        moderator_data$SingleSpecies
      )
    )
  ) == "yes"
)


# Temperature tolerance
moderator_data$TemperatureTolerance <- clean_category(
  moderator_data$SuitableTemperatureRange
)

moderator_data$TemperatureTolerance[
  !single_species_observation
] <- NA_character_

moderator_data$TemperatureTolerance <- factor(
  moderator_data$TemperatureTolerance,
  levels = c(
    "NotWarmTolerant",
    "WarmTolerant"
  )
)


# Nutrient tolerance
moderator_data$NutrientTolerance <- clean_category(
  moderator_data$SuitableNutrientRange
)

moderator_data$NutrientTolerance[
  !single_species_observation
] <- NA_character_

moderator_data$NutrientTolerance <- factor(
  moderator_data$NutrientTolerance,
  levels = c(
    "EutrophicSensitive",
    "EutrophicTolerant"
  )
)


# Photosynthetic pathway
moderator_data$PhotosyntheticPathway <- clean_category(
  moderator_data$PhotosynthesicPathway
)

moderator_data$PhotosyntheticPathway[
  !single_species_observation
] <- NA_character_

moderator_data$PhotosyntheticPathway <- factor(
  moderator_data$PhotosyntheticPathway,
  levels = c(
    "C3",
    "NotC3"
  )
)


# Retain stressor pair for the later full model
moderator_data$FactorPair <- factor(
  moderator_data$FactorPair
)


# 05.2. Eight moderator analyses written explicitly ---------------------------
#
# Each moderator is analysed in a separate block. No fitting function or model
# loop is used, so every data filter, covariance matrix, multilevel model and
# CR2 correction can be inspected directly. The same observations are retained
# across the additive, multiplicative and dominance models within a moderator.
# Continuous moderators are standardized; their slopes therefore represent a
# one-standard-deviation increase in the original variable.

sampling_correlation <- 0.6

# 05.2.1. Background nutrient level -------------------------------------------

background_nutrient_data <- moderator_data[
    !is.na(moderator_data$BackgroundNutrientLevel) &
  trimws(as.character(moderator_data$BackgroundNutrientLevel)) != "" &
    is.finite(moderator_data$HedgesdAdditive) &
    is.finite(moderator_data$VdAdditive) &
    moderator_data$VdAdditive > 0 &
    is.finite(moderator_data$HedgesdMultiplicative) &
    is.finite(moderator_data$VdMultiplicative) &
    moderator_data$VdMultiplicative > 0 &
    is.finite(moderator_data$HedgesdDominance) &
    is.finite(moderator_data$VdDominance) &
    moderator_data$VdDominance > 0 &
    !is.na(moderator_data$StudyID) &
    !is.na(moderator_data$ObservationID),
  ,
  drop = FALSE
]

background_nutrient_data$ModeratorValue <- droplevels(
  factor(background_nutrient_data$BackgroundNutrientLevel)
)

background_nutrient_data$StudyID <- factor(
  background_nutrient_data$StudyID
)

background_nutrient_data$ObservationWithinStudy <- factor(
  paste(
    background_nutrient_data$StudyID,
    background_nutrient_data$ObservationID,
    sep = "_"
  )
)

background_nutrient_levels <- levels(
  background_nutrient_data$ModeratorValue
)

background_nutrient_N <- length(
  unique(background_nutrient_data$StudyID)
)

background_nutrient_n <- nrow(
  background_nutrient_data
)

background_nutrient_level_n_table <- as.data.frame(
  table(background_nutrient_data$ModeratorValue),
  stringsAsFactors = FALSE
)

background_nutrient_level_study_data <- unique(
  data.frame(
    Level = background_nutrient_data$ModeratorValue,
    StudyID = background_nutrient_data$StudyID
  )
)

background_nutrient_level_N_table <- aggregate(
  StudyID ~ Level,
  data = background_nutrient_level_study_data,
  FUN = length
)

background_nutrient_level_n <- as.integer(
  background_nutrient_level_n_table$Freq[
    match(
      background_nutrient_levels,
      as.character(background_nutrient_level_n_table$Var1)
    )
  ]
)

background_nutrient_level_N <- as.integer(
  background_nutrient_level_N_table$StudyID[
    match(
      background_nutrient_levels,
      as.character(background_nutrient_level_N_table$Level)
    )
  ]
)

# Additive null model
background_nutrient_additive_V <- clubSandwich::impute_covariance_matrix(
  vi = background_nutrient_data$VdAdditive,
  cluster = background_nutrient_data$StudyID,
  r = sampling_correlation,
  smooth_vi = FALSE,
  check_PD = TRUE,
  return_list = FALSE
)

background_nutrient_additive_level_model <- metafor::rma.mv(
  yi = HedgesdAdditive,
  V = background_nutrient_additive_V,
  mods = ~ 0 + ModeratorValue,
  random = list(
    ~ 1 | StudyID,
    ~ 1 | ObservationWithinStudy
  ),
  data = background_nutrient_data,
  method = "REML",
  test = "t",
  sparse = TRUE
)

background_nutrient_additive_contrast_model <- metafor::rma.mv(
  yi = HedgesdAdditive,
  V = background_nutrient_additive_V,
  mods = ~ ModeratorValue,
  random = list(
    ~ 1 | StudyID,
    ~ 1 | ObservationWithinStudy
  ),
  data = background_nutrient_data,
  method = "REML",
  test = "t",
  sparse = TRUE
)

background_nutrient_additive_level_CR2 <- as.data.frame(
  clubSandwich::coef_test(
    background_nutrient_additive_level_model,
    vcov = "CR2",
    cluster = background_nutrient_data$StudyID,
    test = "Satterthwaite"
  )
)

background_nutrient_additive_contrast_CR2 <- as.data.frame(
  clubSandwich::coef_test(
    background_nutrient_additive_contrast_model,
    vcov = "CR2",
    cluster = background_nutrient_data$StudyID,
    test = "Satterthwaite"
  )
)

background_nutrient_additive_level_critical <- stats::qt(
  0.975,
  df = background_nutrient_additive_level_CR2$df_Satt
)

background_nutrient_additive_level_lower <-
  background_nutrient_additive_level_CR2$beta -
  background_nutrient_additive_level_critical * background_nutrient_additive_level_CR2$SE

background_nutrient_additive_level_upper <-
  background_nutrient_additive_level_CR2$beta +
  background_nutrient_additive_level_critical * background_nutrient_additive_level_CR2$SE

background_nutrient_additive_contrast_CR2 <- background_nutrient_additive_contrast_CR2[
  -1,
  ,
  drop = FALSE
]

background_nutrient_additive_contrast_critical <- stats::qt(
  0.975,
  df = background_nutrient_additive_contrast_CR2$df_Satt
)

background_nutrient_additive_contrast_lower <-
  background_nutrient_additive_contrast_CR2$beta -
  background_nutrient_additive_contrast_critical * background_nutrient_additive_contrast_CR2$SE

background_nutrient_additive_contrast_upper <-
  background_nutrient_additive_contrast_CR2$beta +
  background_nutrient_additive_contrast_critical * background_nutrient_additive_contrast_CR2$SE

background_nutrient_additive_W <- solve(background_nutrient_additive_V)
background_nutrient_additive_X <- as.matrix(background_nutrient_additive_contrast_model$X)
background_nutrient_additive_P <-
  background_nutrient_additive_W -
  background_nutrient_additive_W %*% background_nutrient_additive_X %*%
  solve(crossprod(background_nutrient_additive_X, background_nutrient_additive_W %*% background_nutrient_additive_X)) %*%
  t(background_nutrient_additive_X) %*% background_nutrient_additive_W

background_nutrient_additive_typical_sampling_variance <-
  (nrow(background_nutrient_additive_X) - ncol(background_nutrient_additive_X)) /
  sum(diag(background_nutrient_additive_P))

background_nutrient_additive_total_variance <-
  sum(background_nutrient_additive_contrast_model$sigma2) +
  background_nutrient_additive_typical_sampling_variance

background_nutrient_additive_study_I2 <-
  100 * background_nutrient_additive_contrast_model$sigma2[1] /
  background_nutrient_additive_total_variance

background_nutrient_additive_observation_I2 <-
  100 * background_nutrient_additive_contrast_model$sigma2[2] /
  background_nutrient_additive_total_variance

background_nutrient_additive_total_I2 <-
  100 * sum(background_nutrient_additive_contrast_model$sigma2) /
  background_nutrient_additive_total_variance

background_nutrient_additive_level_results <- data.frame(
  Moderator = "Background nutrient level",
  NullModel = "Additive",
  Term = background_nutrient_levels,
  EffectType = "Level mean",
  ReferenceLevel = background_nutrient_levels[1],
  N = background_nutrient_N,
  n = background_nutrient_n,
  Level_N = background_nutrient_level_N,
  Level_n = background_nutrient_level_n,
  Estimate = as.numeric(background_nutrient_additive_level_CR2$beta),
  CR2_SE = as.numeric(background_nutrient_additive_level_CR2$SE),
  CR2_df = as.numeric(background_nutrient_additive_level_CR2$df_Satt),
  CI_lower = as.numeric(background_nutrient_additive_level_lower),
  CI_upper = as.numeric(background_nutrient_additive_level_upper),
  P_value = as.numeric(background_nutrient_additive_level_CR2$p_Satt),
  Moderator_QM = as.numeric(background_nutrient_additive_contrast_model$QM),
  Moderator_QM_P = as.numeric(background_nutrient_additive_contrast_model$QMp),
  Tau2_Study = background_nutrient_additive_contrast_model$sigma2[1],
  Tau2_Observation = background_nutrient_additive_contrast_model$sigma2[2],
  Study_I2 = background_nutrient_additive_study_I2,
  Observation_I2 = background_nutrient_additive_observation_I2,
  Total_I2 = background_nutrient_additive_total_I2,
  QE = as.numeric(background_nutrient_additive_contrast_model$QE),
  QE_P_value = as.numeric(background_nutrient_additive_contrast_model$QEp),
  Interpretation = ifelse(
    background_nutrient_additive_level_lower > 0,
    "Positive non-additive",
    ifelse(
      background_nutrient_additive_level_upper < 0,
      "Negative non-additive",
      "Additive"
    )
  ),
  Status = "Successfully fitted",
  stringsAsFactors = FALSE
)

background_nutrient_additive_contrast_results <- data.frame(
  Moderator = "Background nutrient level",
  NullModel = "Additive",
  Term = paste0(
    background_nutrient_levels[-1],
    " vs ",
    background_nutrient_levels[1]
  ),
  EffectType = "Contrast with reference",
  ReferenceLevel = background_nutrient_levels[1],
  N = background_nutrient_N,
  n = background_nutrient_n,
  Level_N = background_nutrient_N,
  Level_n = background_nutrient_n,
  Estimate = as.numeric(background_nutrient_additive_contrast_CR2$beta),
  CR2_SE = as.numeric(background_nutrient_additive_contrast_CR2$SE),
  CR2_df = as.numeric(background_nutrient_additive_contrast_CR2$df_Satt),
  CI_lower = as.numeric(background_nutrient_additive_contrast_lower),
  CI_upper = as.numeric(background_nutrient_additive_contrast_upper),
  P_value = as.numeric(background_nutrient_additive_contrast_CR2$p_Satt),
  Moderator_QM = as.numeric(background_nutrient_additive_contrast_model$QM),
  Moderator_QM_P = as.numeric(background_nutrient_additive_contrast_model$QMp),
  Tau2_Study = background_nutrient_additive_contrast_model$sigma2[1],
  Tau2_Observation = background_nutrient_additive_contrast_model$sigma2[2],
  Study_I2 = background_nutrient_additive_study_I2,
  Observation_I2 = background_nutrient_additive_observation_I2,
  Total_I2 = background_nutrient_additive_total_I2,
  QE = as.numeric(background_nutrient_additive_contrast_model$QE),
  QE_P_value = as.numeric(background_nutrient_additive_contrast_model$QEp),
  Interpretation = ifelse(
    background_nutrient_additive_contrast_lower > 0,
    "Higher than reference",
    ifelse(
      background_nutrient_additive_contrast_upper < 0,
      "Lower than reference",
      "No detected difference"
    )
  ),
  Status = "Successfully fitted",
  stringsAsFactors = FALSE
)

background_nutrient_additive_results <- dplyr::bind_rows(
  background_nutrient_additive_level_results,
  background_nutrient_additive_contrast_results
)


# Multiplicative null model
background_nutrient_multiplicative_V <- clubSandwich::impute_covariance_matrix(
  vi = background_nutrient_data$VdMultiplicative,
  cluster = background_nutrient_data$StudyID,
  r = sampling_correlation,
  smooth_vi = FALSE,
  check_PD = TRUE,
  return_list = FALSE
)

background_nutrient_multiplicative_level_model <- metafor::rma.mv(
  yi = HedgesdMultiplicative,
  V = background_nutrient_multiplicative_V,
  mods = ~ 0 + ModeratorValue,
  random = list(
    ~ 1 | StudyID,
    ~ 1 | ObservationWithinStudy
  ),
  data = background_nutrient_data,
  method = "REML",
  test = "t",
  sparse = TRUE
)

background_nutrient_multiplicative_contrast_model <- metafor::rma.mv(
  yi = HedgesdMultiplicative,
  V = background_nutrient_multiplicative_V,
  mods = ~ ModeratorValue,
  random = list(
    ~ 1 | StudyID,
    ~ 1 | ObservationWithinStudy
  ),
  data = background_nutrient_data,
  method = "REML",
  test = "t",
  sparse = TRUE
)

background_nutrient_multiplicative_level_CR2 <- as.data.frame(
  clubSandwich::coef_test(
    background_nutrient_multiplicative_level_model,
    vcov = "CR2",
    cluster = background_nutrient_data$StudyID,
    test = "Satterthwaite"
  )
)

background_nutrient_multiplicative_contrast_CR2 <- as.data.frame(
  clubSandwich::coef_test(
    background_nutrient_multiplicative_contrast_model,
    vcov = "CR2",
    cluster = background_nutrient_data$StudyID,
    test = "Satterthwaite"
  )
)

background_nutrient_multiplicative_level_critical <- stats::qt(
  0.975,
  df = background_nutrient_multiplicative_level_CR2$df_Satt
)

background_nutrient_multiplicative_level_lower <-
  background_nutrient_multiplicative_level_CR2$beta -
  background_nutrient_multiplicative_level_critical * background_nutrient_multiplicative_level_CR2$SE

background_nutrient_multiplicative_level_upper <-
  background_nutrient_multiplicative_level_CR2$beta +
  background_nutrient_multiplicative_level_critical * background_nutrient_multiplicative_level_CR2$SE

background_nutrient_multiplicative_contrast_CR2 <- background_nutrient_multiplicative_contrast_CR2[
  -1,
  ,
  drop = FALSE
]

background_nutrient_multiplicative_contrast_critical <- stats::qt(
  0.975,
  df = background_nutrient_multiplicative_contrast_CR2$df_Satt
)

background_nutrient_multiplicative_contrast_lower <-
  background_nutrient_multiplicative_contrast_CR2$beta -
  background_nutrient_multiplicative_contrast_critical * background_nutrient_multiplicative_contrast_CR2$SE

background_nutrient_multiplicative_contrast_upper <-
  background_nutrient_multiplicative_contrast_CR2$beta +
  background_nutrient_multiplicative_contrast_critical * background_nutrient_multiplicative_contrast_CR2$SE

background_nutrient_multiplicative_W <- solve(background_nutrient_multiplicative_V)
background_nutrient_multiplicative_X <- as.matrix(background_nutrient_multiplicative_contrast_model$X)
background_nutrient_multiplicative_P <-
  background_nutrient_multiplicative_W -
  background_nutrient_multiplicative_W %*% background_nutrient_multiplicative_X %*%
  solve(crossprod(background_nutrient_multiplicative_X, background_nutrient_multiplicative_W %*% background_nutrient_multiplicative_X)) %*%
  t(background_nutrient_multiplicative_X) %*% background_nutrient_multiplicative_W

background_nutrient_multiplicative_typical_sampling_variance <-
  (nrow(background_nutrient_multiplicative_X) - ncol(background_nutrient_multiplicative_X)) /
  sum(diag(background_nutrient_multiplicative_P))

background_nutrient_multiplicative_total_variance <-
  sum(background_nutrient_multiplicative_contrast_model$sigma2) +
  background_nutrient_multiplicative_typical_sampling_variance

background_nutrient_multiplicative_study_I2 <-
  100 * background_nutrient_multiplicative_contrast_model$sigma2[1] /
  background_nutrient_multiplicative_total_variance

background_nutrient_multiplicative_observation_I2 <-
  100 * background_nutrient_multiplicative_contrast_model$sigma2[2] /
  background_nutrient_multiplicative_total_variance

background_nutrient_multiplicative_total_I2 <-
  100 * sum(background_nutrient_multiplicative_contrast_model$sigma2) /
  background_nutrient_multiplicative_total_variance

background_nutrient_multiplicative_level_results <- data.frame(
  Moderator = "Background nutrient level",
  NullModel = "Multiplicative",
  Term = background_nutrient_levels,
  EffectType = "Level mean",
  ReferenceLevel = background_nutrient_levels[1],
  N = background_nutrient_N,
  n = background_nutrient_n,
  Level_N = background_nutrient_level_N,
  Level_n = background_nutrient_level_n,
  Estimate = as.numeric(background_nutrient_multiplicative_level_CR2$beta),
  CR2_SE = as.numeric(background_nutrient_multiplicative_level_CR2$SE),
  CR2_df = as.numeric(background_nutrient_multiplicative_level_CR2$df_Satt),
  CI_lower = as.numeric(background_nutrient_multiplicative_level_lower),
  CI_upper = as.numeric(background_nutrient_multiplicative_level_upper),
  P_value = as.numeric(background_nutrient_multiplicative_level_CR2$p_Satt),
  Moderator_QM = as.numeric(background_nutrient_multiplicative_contrast_model$QM),
  Moderator_QM_P = as.numeric(background_nutrient_multiplicative_contrast_model$QMp),
  Tau2_Study = background_nutrient_multiplicative_contrast_model$sigma2[1],
  Tau2_Observation = background_nutrient_multiplicative_contrast_model$sigma2[2],
  Study_I2 = background_nutrient_multiplicative_study_I2,
  Observation_I2 = background_nutrient_multiplicative_observation_I2,
  Total_I2 = background_nutrient_multiplicative_total_I2,
  QE = as.numeric(background_nutrient_multiplicative_contrast_model$QE),
  QE_P_value = as.numeric(background_nutrient_multiplicative_contrast_model$QEp),
  Interpretation = ifelse(
    background_nutrient_multiplicative_level_lower > 0,
    "Positive non-additive",
    ifelse(
      background_nutrient_multiplicative_level_upper < 0,
      "Negative non-additive",
      "Additive"
    )
  ),
  Status = "Successfully fitted",
  stringsAsFactors = FALSE
)

background_nutrient_multiplicative_contrast_results <- data.frame(
  Moderator = "Background nutrient level",
  NullModel = "Multiplicative",
  Term = paste0(
    background_nutrient_levels[-1],
    " vs ",
    background_nutrient_levels[1]
  ),
  EffectType = "Contrast with reference",
  ReferenceLevel = background_nutrient_levels[1],
  N = background_nutrient_N,
  n = background_nutrient_n,
  Level_N = background_nutrient_N,
  Level_n = background_nutrient_n,
  Estimate = as.numeric(background_nutrient_multiplicative_contrast_CR2$beta),
  CR2_SE = as.numeric(background_nutrient_multiplicative_contrast_CR2$SE),
  CR2_df = as.numeric(background_nutrient_multiplicative_contrast_CR2$df_Satt),
  CI_lower = as.numeric(background_nutrient_multiplicative_contrast_lower),
  CI_upper = as.numeric(background_nutrient_multiplicative_contrast_upper),
  P_value = as.numeric(background_nutrient_multiplicative_contrast_CR2$p_Satt),
  Moderator_QM = as.numeric(background_nutrient_multiplicative_contrast_model$QM),
  Moderator_QM_P = as.numeric(background_nutrient_multiplicative_contrast_model$QMp),
  Tau2_Study = background_nutrient_multiplicative_contrast_model$sigma2[1],
  Tau2_Observation = background_nutrient_multiplicative_contrast_model$sigma2[2],
  Study_I2 = background_nutrient_multiplicative_study_I2,
  Observation_I2 = background_nutrient_multiplicative_observation_I2,
  Total_I2 = background_nutrient_multiplicative_total_I2,
  QE = as.numeric(background_nutrient_multiplicative_contrast_model$QE),
  QE_P_value = as.numeric(background_nutrient_multiplicative_contrast_model$QEp),
  Interpretation = ifelse(
    background_nutrient_multiplicative_contrast_lower > 0,
    "Higher than reference",
    ifelse(
      background_nutrient_multiplicative_contrast_upper < 0,
      "Lower than reference",
      "No detected difference"
    )
  ),
  Status = "Successfully fitted",
  stringsAsFactors = FALSE
)

background_nutrient_multiplicative_results <- dplyr::bind_rows(
  background_nutrient_multiplicative_level_results,
  background_nutrient_multiplicative_contrast_results
)


# Dominance null model
background_nutrient_dominance_V <- clubSandwich::impute_covariance_matrix(
  vi = background_nutrient_data$VdDominance,
  cluster = background_nutrient_data$StudyID,
  r = sampling_correlation,
  smooth_vi = FALSE,
  check_PD = TRUE,
  return_list = FALSE
)

background_nutrient_dominance_level_model <- metafor::rma.mv(
  yi = HedgesdDominance,
  V = background_nutrient_dominance_V,
  mods = ~ 0 + ModeratorValue,
  random = list(
    ~ 1 | StudyID,
    ~ 1 | ObservationWithinStudy
  ),
  data = background_nutrient_data,
  method = "REML",
  test = "t",
  sparse = TRUE
)

background_nutrient_dominance_contrast_model <- metafor::rma.mv(
  yi = HedgesdDominance,
  V = background_nutrient_dominance_V,
  mods = ~ ModeratorValue,
  random = list(
    ~ 1 | StudyID,
    ~ 1 | ObservationWithinStudy
  ),
  data = background_nutrient_data,
  method = "REML",
  test = "t",
  sparse = TRUE
)

background_nutrient_dominance_level_CR2 <- as.data.frame(
  clubSandwich::coef_test(
    background_nutrient_dominance_level_model,
    vcov = "CR2",
    cluster = background_nutrient_data$StudyID,
    test = "Satterthwaite"
  )
)

background_nutrient_dominance_contrast_CR2 <- as.data.frame(
  clubSandwich::coef_test(
    background_nutrient_dominance_contrast_model,
    vcov = "CR2",
    cluster = background_nutrient_data$StudyID,
    test = "Satterthwaite"
  )
)

background_nutrient_dominance_level_critical <- stats::qt(
  0.975,
  df = background_nutrient_dominance_level_CR2$df_Satt
)

background_nutrient_dominance_level_lower <-
  background_nutrient_dominance_level_CR2$beta -
  background_nutrient_dominance_level_critical * background_nutrient_dominance_level_CR2$SE

background_nutrient_dominance_level_upper <-
  background_nutrient_dominance_level_CR2$beta +
  background_nutrient_dominance_level_critical * background_nutrient_dominance_level_CR2$SE

background_nutrient_dominance_contrast_CR2 <- background_nutrient_dominance_contrast_CR2[
  -1,
  ,
  drop = FALSE
]

background_nutrient_dominance_contrast_critical <- stats::qt(
  0.975,
  df = background_nutrient_dominance_contrast_CR2$df_Satt
)

background_nutrient_dominance_contrast_lower <-
  background_nutrient_dominance_contrast_CR2$beta -
  background_nutrient_dominance_contrast_critical * background_nutrient_dominance_contrast_CR2$SE

background_nutrient_dominance_contrast_upper <-
  background_nutrient_dominance_contrast_CR2$beta +
  background_nutrient_dominance_contrast_critical * background_nutrient_dominance_contrast_CR2$SE

background_nutrient_dominance_W <- solve(background_nutrient_dominance_V)
background_nutrient_dominance_X <- as.matrix(background_nutrient_dominance_contrast_model$X)
background_nutrient_dominance_P <-
  background_nutrient_dominance_W -
  background_nutrient_dominance_W %*% background_nutrient_dominance_X %*%
  solve(crossprod(background_nutrient_dominance_X, background_nutrient_dominance_W %*% background_nutrient_dominance_X)) %*%
  t(background_nutrient_dominance_X) %*% background_nutrient_dominance_W

background_nutrient_dominance_typical_sampling_variance <-
  (nrow(background_nutrient_dominance_X) - ncol(background_nutrient_dominance_X)) /
  sum(diag(background_nutrient_dominance_P))

background_nutrient_dominance_total_variance <-
  sum(background_nutrient_dominance_contrast_model$sigma2) +
  background_nutrient_dominance_typical_sampling_variance

background_nutrient_dominance_study_I2 <-
  100 * background_nutrient_dominance_contrast_model$sigma2[1] /
  background_nutrient_dominance_total_variance

background_nutrient_dominance_observation_I2 <-
  100 * background_nutrient_dominance_contrast_model$sigma2[2] /
  background_nutrient_dominance_total_variance

background_nutrient_dominance_total_I2 <-
  100 * sum(background_nutrient_dominance_contrast_model$sigma2) /
  background_nutrient_dominance_total_variance

background_nutrient_dominance_level_results <- data.frame(
  Moderator = "Background nutrient level",
  NullModel = "Dominance",
  Term = background_nutrient_levels,
  EffectType = "Level mean",
  ReferenceLevel = background_nutrient_levels[1],
  N = background_nutrient_N,
  n = background_nutrient_n,
  Level_N = background_nutrient_level_N,
  Level_n = background_nutrient_level_n,
  Estimate = as.numeric(background_nutrient_dominance_level_CR2$beta),
  CR2_SE = as.numeric(background_nutrient_dominance_level_CR2$SE),
  CR2_df = as.numeric(background_nutrient_dominance_level_CR2$df_Satt),
  CI_lower = as.numeric(background_nutrient_dominance_level_lower),
  CI_upper = as.numeric(background_nutrient_dominance_level_upper),
  P_value = as.numeric(background_nutrient_dominance_level_CR2$p_Satt),
  Moderator_QM = as.numeric(background_nutrient_dominance_contrast_model$QM),
  Moderator_QM_P = as.numeric(background_nutrient_dominance_contrast_model$QMp),
  Tau2_Study = background_nutrient_dominance_contrast_model$sigma2[1],
  Tau2_Observation = background_nutrient_dominance_contrast_model$sigma2[2],
  Study_I2 = background_nutrient_dominance_study_I2,
  Observation_I2 = background_nutrient_dominance_observation_I2,
  Total_I2 = background_nutrient_dominance_total_I2,
  QE = as.numeric(background_nutrient_dominance_contrast_model$QE),
  QE_P_value = as.numeric(background_nutrient_dominance_contrast_model$QEp),
  Interpretation = ifelse(
    background_nutrient_dominance_level_lower > 0,
    "Positive non-additive",
    ifelse(
      background_nutrient_dominance_level_upper < 0,
      "Negative non-additive",
      "Additive"
    )
  ),
  Status = "Successfully fitted",
  stringsAsFactors = FALSE
)

background_nutrient_dominance_contrast_results <- data.frame(
  Moderator = "Background nutrient level",
  NullModel = "Dominance",
  Term = paste0(
    background_nutrient_levels[-1],
    " vs ",
    background_nutrient_levels[1]
  ),
  EffectType = "Contrast with reference",
  ReferenceLevel = background_nutrient_levels[1],
  N = background_nutrient_N,
  n = background_nutrient_n,
  Level_N = background_nutrient_N,
  Level_n = background_nutrient_n,
  Estimate = as.numeric(background_nutrient_dominance_contrast_CR2$beta),
  CR2_SE = as.numeric(background_nutrient_dominance_contrast_CR2$SE),
  CR2_df = as.numeric(background_nutrient_dominance_contrast_CR2$df_Satt),
  CI_lower = as.numeric(background_nutrient_dominance_contrast_lower),
  CI_upper = as.numeric(background_nutrient_dominance_contrast_upper),
  P_value = as.numeric(background_nutrient_dominance_contrast_CR2$p_Satt),
  Moderator_QM = as.numeric(background_nutrient_dominance_contrast_model$QM),
  Moderator_QM_P = as.numeric(background_nutrient_dominance_contrast_model$QMp),
  Tau2_Study = background_nutrient_dominance_contrast_model$sigma2[1],
  Tau2_Observation = background_nutrient_dominance_contrast_model$sigma2[2],
  Study_I2 = background_nutrient_dominance_study_I2,
  Observation_I2 = background_nutrient_dominance_observation_I2,
  Total_I2 = background_nutrient_dominance_total_I2,
  QE = as.numeric(background_nutrient_dominance_contrast_model$QE),
  QE_P_value = as.numeric(background_nutrient_dominance_contrast_model$QEp),
  Interpretation = ifelse(
    background_nutrient_dominance_contrast_lower > 0,
    "Higher than reference",
    ifelse(
      background_nutrient_dominance_contrast_upper < 0,
      "Lower than reference",
      "No detected difference"
    )
  ),
  Status = "Successfully fitted",
  stringsAsFactors = FALSE
)

background_nutrient_dominance_results <- dplyr::bind_rows(
  background_nutrient_dominance_level_results,
  background_nutrient_dominance_contrast_results
)

background_nutrient_results <- dplyr::bind_rows(
  background_nutrient_additive_results,
  background_nutrient_multiplicative_results,
  background_nutrient_dominance_results
)



# 05.2.2. Fish presence --------------------------------------------------------

fish_data <- moderator_data[
    !is.na(moderator_data$FishPresence) &
  trimws(as.character(moderator_data$FishPresence)) != "" &
    is.finite(moderator_data$HedgesdAdditive) &
    is.finite(moderator_data$VdAdditive) &
    moderator_data$VdAdditive > 0 &
    is.finite(moderator_data$HedgesdMultiplicative) &
    is.finite(moderator_data$VdMultiplicative) &
    moderator_data$VdMultiplicative > 0 &
    is.finite(moderator_data$HedgesdDominance) &
    is.finite(moderator_data$VdDominance) &
    moderator_data$VdDominance > 0 &
    !is.na(moderator_data$StudyID) &
    !is.na(moderator_data$ObservationID),
  ,
  drop = FALSE
]

fish_data$ModeratorValue <- droplevels(
  factor(fish_data$FishPresence)
)

fish_data$StudyID <- factor(
  fish_data$StudyID
)

fish_data$ObservationWithinStudy <- factor(
  paste(
    fish_data$StudyID,
    fish_data$ObservationID,
    sep = "_"
  )
)

fish_levels <- levels(
  fish_data$ModeratorValue
)

fish_N <- length(
  unique(fish_data$StudyID)
)

fish_n <- nrow(
  fish_data
)

fish_level_n_table <- as.data.frame(
  table(fish_data$ModeratorValue),
  stringsAsFactors = FALSE
)

fish_level_study_data <- unique(
  data.frame(
    Level = fish_data$ModeratorValue,
    StudyID = fish_data$StudyID
  )
)

fish_level_N_table <- aggregate(
  StudyID ~ Level,
  data = fish_level_study_data,
  FUN = length
)

fish_level_n <- as.integer(
  fish_level_n_table$Freq[
    match(
      fish_levels,
      as.character(fish_level_n_table$Var1)
    )
  ]
)

fish_level_N <- as.integer(
  fish_level_N_table$StudyID[
    match(
      fish_levels,
      as.character(fish_level_N_table$Level)
    )
  ]
)

# Additive null model
fish_additive_V <- clubSandwich::impute_covariance_matrix(
  vi = fish_data$VdAdditive,
  cluster = fish_data$StudyID,
  r = sampling_correlation,
  smooth_vi = FALSE,
  check_PD = TRUE,
  return_list = FALSE
)

fish_additive_level_model <- metafor::rma.mv(
  yi = HedgesdAdditive,
  V = fish_additive_V,
  mods = ~ 0 + ModeratorValue,
  random = list(
    ~ 1 | StudyID,
    ~ 1 | ObservationWithinStudy
  ),
  data = fish_data,
  method = "REML",
  test = "t",
  sparse = TRUE
)

fish_additive_contrast_model <- metafor::rma.mv(
  yi = HedgesdAdditive,
  V = fish_additive_V,
  mods = ~ ModeratorValue,
  random = list(
    ~ 1 | StudyID,
    ~ 1 | ObservationWithinStudy
  ),
  data = fish_data,
  method = "REML",
  test = "t",
  sparse = TRUE
)

fish_additive_level_CR2 <- as.data.frame(
  clubSandwich::coef_test(
    fish_additive_level_model,
    vcov = "CR2",
    cluster = fish_data$StudyID,
    test = "Satterthwaite"
  )
)

fish_additive_contrast_CR2 <- as.data.frame(
  clubSandwich::coef_test(
    fish_additive_contrast_model,
    vcov = "CR2",
    cluster = fish_data$StudyID,
    test = "Satterthwaite"
  )
)

fish_additive_level_critical <- stats::qt(
  0.975,
  df = fish_additive_level_CR2$df_Satt
)

fish_additive_level_lower <-
  fish_additive_level_CR2$beta -
  fish_additive_level_critical * fish_additive_level_CR2$SE

fish_additive_level_upper <-
  fish_additive_level_CR2$beta +
  fish_additive_level_critical * fish_additive_level_CR2$SE

fish_additive_contrast_CR2 <- fish_additive_contrast_CR2[
  -1,
  ,
  drop = FALSE
]

fish_additive_contrast_critical <- stats::qt(
  0.975,
  df = fish_additive_contrast_CR2$df_Satt
)

fish_additive_contrast_lower <-
  fish_additive_contrast_CR2$beta -
  fish_additive_contrast_critical * fish_additive_contrast_CR2$SE

fish_additive_contrast_upper <-
  fish_additive_contrast_CR2$beta +
  fish_additive_contrast_critical * fish_additive_contrast_CR2$SE

fish_additive_W <- solve(fish_additive_V)
fish_additive_X <- as.matrix(fish_additive_contrast_model$X)
fish_additive_P <-
  fish_additive_W -
  fish_additive_W %*% fish_additive_X %*%
  solve(crossprod(fish_additive_X, fish_additive_W %*% fish_additive_X)) %*%
  t(fish_additive_X) %*% fish_additive_W

fish_additive_typical_sampling_variance <-
  (nrow(fish_additive_X) - ncol(fish_additive_X)) /
  sum(diag(fish_additive_P))

fish_additive_total_variance <-
  sum(fish_additive_contrast_model$sigma2) +
  fish_additive_typical_sampling_variance

fish_additive_study_I2 <-
  100 * fish_additive_contrast_model$sigma2[1] /
  fish_additive_total_variance

fish_additive_observation_I2 <-
  100 * fish_additive_contrast_model$sigma2[2] /
  fish_additive_total_variance

fish_additive_total_I2 <-
  100 * sum(fish_additive_contrast_model$sigma2) /
  fish_additive_total_variance

fish_additive_level_results <- data.frame(
  Moderator = "Fish presence",
  NullModel = "Additive",
  Term = fish_levels,
  EffectType = "Level mean",
  ReferenceLevel = fish_levels[1],
  N = fish_N,
  n = fish_n,
  Level_N = fish_level_N,
  Level_n = fish_level_n,
  Estimate = as.numeric(fish_additive_level_CR2$beta),
  CR2_SE = as.numeric(fish_additive_level_CR2$SE),
  CR2_df = as.numeric(fish_additive_level_CR2$df_Satt),
  CI_lower = as.numeric(fish_additive_level_lower),
  CI_upper = as.numeric(fish_additive_level_upper),
  P_value = as.numeric(fish_additive_level_CR2$p_Satt),
  Moderator_QM = as.numeric(fish_additive_contrast_model$QM),
  Moderator_QM_P = as.numeric(fish_additive_contrast_model$QMp),
  Tau2_Study = fish_additive_contrast_model$sigma2[1],
  Tau2_Observation = fish_additive_contrast_model$sigma2[2],
  Study_I2 = fish_additive_study_I2,
  Observation_I2 = fish_additive_observation_I2,
  Total_I2 = fish_additive_total_I2,
  QE = as.numeric(fish_additive_contrast_model$QE),
  QE_P_value = as.numeric(fish_additive_contrast_model$QEp),
  Interpretation = ifelse(
    fish_additive_level_lower > 0,
    "Positive non-additive",
    ifelse(
      fish_additive_level_upper < 0,
      "Negative non-additive",
      "Additive"
    )
  ),
  Status = "Successfully fitted",
  stringsAsFactors = FALSE
)

fish_additive_contrast_results <- data.frame(
  Moderator = "Fish presence",
  NullModel = "Additive",
  Term = paste0(
    fish_levels[-1],
    " vs ",
    fish_levels[1]
  ),
  EffectType = "Contrast with reference",
  ReferenceLevel = fish_levels[1],
  N = fish_N,
  n = fish_n,
  Level_N = fish_N,
  Level_n = fish_n,
  Estimate = as.numeric(fish_additive_contrast_CR2$beta),
  CR2_SE = as.numeric(fish_additive_contrast_CR2$SE),
  CR2_df = as.numeric(fish_additive_contrast_CR2$df_Satt),
  CI_lower = as.numeric(fish_additive_contrast_lower),
  CI_upper = as.numeric(fish_additive_contrast_upper),
  P_value = as.numeric(fish_additive_contrast_CR2$p_Satt),
  Moderator_QM = as.numeric(fish_additive_contrast_model$QM),
  Moderator_QM_P = as.numeric(fish_additive_contrast_model$QMp),
  Tau2_Study = fish_additive_contrast_model$sigma2[1],
  Tau2_Observation = fish_additive_contrast_model$sigma2[2],
  Study_I2 = fish_additive_study_I2,
  Observation_I2 = fish_additive_observation_I2,
  Total_I2 = fish_additive_total_I2,
  QE = as.numeric(fish_additive_contrast_model$QE),
  QE_P_value = as.numeric(fish_additive_contrast_model$QEp),
  Interpretation = ifelse(
    fish_additive_contrast_lower > 0,
    "Higher than reference",
    ifelse(
      fish_additive_contrast_upper < 0,
      "Lower than reference",
      "No detected difference"
    )
  ),
  Status = "Successfully fitted",
  stringsAsFactors = FALSE
)

fish_additive_results <- dplyr::bind_rows(
  fish_additive_level_results,
  fish_additive_contrast_results
)


# Multiplicative null model
fish_multiplicative_V <- clubSandwich::impute_covariance_matrix(
  vi = fish_data$VdMultiplicative,
  cluster = fish_data$StudyID,
  r = sampling_correlation,
  smooth_vi = FALSE,
  check_PD = TRUE,
  return_list = FALSE
)

fish_multiplicative_level_model <- metafor::rma.mv(
  yi = HedgesdMultiplicative,
  V = fish_multiplicative_V,
  mods = ~ 0 + ModeratorValue,
  random = list(
    ~ 1 | StudyID,
    ~ 1 | ObservationWithinStudy
  ),
  data = fish_data,
  method = "REML",
  test = "t",
  sparse = TRUE
)

fish_multiplicative_contrast_model <- metafor::rma.mv(
  yi = HedgesdMultiplicative,
  V = fish_multiplicative_V,
  mods = ~ ModeratorValue,
  random = list(
    ~ 1 | StudyID,
    ~ 1 | ObservationWithinStudy
  ),
  data = fish_data,
  method = "REML",
  test = "t",
  sparse = TRUE
)

fish_multiplicative_level_CR2 <- as.data.frame(
  clubSandwich::coef_test(
    fish_multiplicative_level_model,
    vcov = "CR2",
    cluster = fish_data$StudyID,
    test = "Satterthwaite"
  )
)

fish_multiplicative_contrast_CR2 <- as.data.frame(
  clubSandwich::coef_test(
    fish_multiplicative_contrast_model,
    vcov = "CR2",
    cluster = fish_data$StudyID,
    test = "Satterthwaite"
  )
)

fish_multiplicative_level_critical <- stats::qt(
  0.975,
  df = fish_multiplicative_level_CR2$df_Satt
)

fish_multiplicative_level_lower <-
  fish_multiplicative_level_CR2$beta -
  fish_multiplicative_level_critical * fish_multiplicative_level_CR2$SE

fish_multiplicative_level_upper <-
  fish_multiplicative_level_CR2$beta +
  fish_multiplicative_level_critical * fish_multiplicative_level_CR2$SE

fish_multiplicative_contrast_CR2 <- fish_multiplicative_contrast_CR2[
  -1,
  ,
  drop = FALSE
]

fish_multiplicative_contrast_critical <- stats::qt(
  0.975,
  df = fish_multiplicative_contrast_CR2$df_Satt
)

fish_multiplicative_contrast_lower <-
  fish_multiplicative_contrast_CR2$beta -
  fish_multiplicative_contrast_critical * fish_multiplicative_contrast_CR2$SE

fish_multiplicative_contrast_upper <-
  fish_multiplicative_contrast_CR2$beta +
  fish_multiplicative_contrast_critical * fish_multiplicative_contrast_CR2$SE

fish_multiplicative_W <- solve(fish_multiplicative_V)
fish_multiplicative_X <- as.matrix(fish_multiplicative_contrast_model$X)
fish_multiplicative_P <-
  fish_multiplicative_W -
  fish_multiplicative_W %*% fish_multiplicative_X %*%
  solve(crossprod(fish_multiplicative_X, fish_multiplicative_W %*% fish_multiplicative_X)) %*%
  t(fish_multiplicative_X) %*% fish_multiplicative_W

fish_multiplicative_typical_sampling_variance <-
  (nrow(fish_multiplicative_X) - ncol(fish_multiplicative_X)) /
  sum(diag(fish_multiplicative_P))

fish_multiplicative_total_variance <-
  sum(fish_multiplicative_contrast_model$sigma2) +
  fish_multiplicative_typical_sampling_variance

fish_multiplicative_study_I2 <-
  100 * fish_multiplicative_contrast_model$sigma2[1] /
  fish_multiplicative_total_variance

fish_multiplicative_observation_I2 <-
  100 * fish_multiplicative_contrast_model$sigma2[2] /
  fish_multiplicative_total_variance

fish_multiplicative_total_I2 <-
  100 * sum(fish_multiplicative_contrast_model$sigma2) /
  fish_multiplicative_total_variance

fish_multiplicative_level_results <- data.frame(
  Moderator = "Fish presence",
  NullModel = "Multiplicative",
  Term = fish_levels,
  EffectType = "Level mean",
  ReferenceLevel = fish_levels[1],
  N = fish_N,
  n = fish_n,
  Level_N = fish_level_N,
  Level_n = fish_level_n,
  Estimate = as.numeric(fish_multiplicative_level_CR2$beta),
  CR2_SE = as.numeric(fish_multiplicative_level_CR2$SE),
  CR2_df = as.numeric(fish_multiplicative_level_CR2$df_Satt),
  CI_lower = as.numeric(fish_multiplicative_level_lower),
  CI_upper = as.numeric(fish_multiplicative_level_upper),
  P_value = as.numeric(fish_multiplicative_level_CR2$p_Satt),
  Moderator_QM = as.numeric(fish_multiplicative_contrast_model$QM),
  Moderator_QM_P = as.numeric(fish_multiplicative_contrast_model$QMp),
  Tau2_Study = fish_multiplicative_contrast_model$sigma2[1],
  Tau2_Observation = fish_multiplicative_contrast_model$sigma2[2],
  Study_I2 = fish_multiplicative_study_I2,
  Observation_I2 = fish_multiplicative_observation_I2,
  Total_I2 = fish_multiplicative_total_I2,
  QE = as.numeric(fish_multiplicative_contrast_model$QE),
  QE_P_value = as.numeric(fish_multiplicative_contrast_model$QEp),
  Interpretation = ifelse(
    fish_multiplicative_level_lower > 0,
    "Positive non-additive",
    ifelse(
      fish_multiplicative_level_upper < 0,
      "Negative non-additive",
      "Additive"
    )
  ),
  Status = "Successfully fitted",
  stringsAsFactors = FALSE
)

fish_multiplicative_contrast_results <- data.frame(
  Moderator = "Fish presence",
  NullModel = "Multiplicative",
  Term = paste0(
    fish_levels[-1],
    " vs ",
    fish_levels[1]
  ),
  EffectType = "Contrast with reference",
  ReferenceLevel = fish_levels[1],
  N = fish_N,
  n = fish_n,
  Level_N = fish_N,
  Level_n = fish_n,
  Estimate = as.numeric(fish_multiplicative_contrast_CR2$beta),
  CR2_SE = as.numeric(fish_multiplicative_contrast_CR2$SE),
  CR2_df = as.numeric(fish_multiplicative_contrast_CR2$df_Satt),
  CI_lower = as.numeric(fish_multiplicative_contrast_lower),
  CI_upper = as.numeric(fish_multiplicative_contrast_upper),
  P_value = as.numeric(fish_multiplicative_contrast_CR2$p_Satt),
  Moderator_QM = as.numeric(fish_multiplicative_contrast_model$QM),
  Moderator_QM_P = as.numeric(fish_multiplicative_contrast_model$QMp),
  Tau2_Study = fish_multiplicative_contrast_model$sigma2[1],
  Tau2_Observation = fish_multiplicative_contrast_model$sigma2[2],
  Study_I2 = fish_multiplicative_study_I2,
  Observation_I2 = fish_multiplicative_observation_I2,
  Total_I2 = fish_multiplicative_total_I2,
  QE = as.numeric(fish_multiplicative_contrast_model$QE),
  QE_P_value = as.numeric(fish_multiplicative_contrast_model$QEp),
  Interpretation = ifelse(
    fish_multiplicative_contrast_lower > 0,
    "Higher than reference",
    ifelse(
      fish_multiplicative_contrast_upper < 0,
      "Lower than reference",
      "No detected difference"
    )
  ),
  Status = "Successfully fitted",
  stringsAsFactors = FALSE
)

fish_multiplicative_results <- dplyr::bind_rows(
  fish_multiplicative_level_results,
  fish_multiplicative_contrast_results
)


# Dominance null model
fish_dominance_V <- clubSandwich::impute_covariance_matrix(
  vi = fish_data$VdDominance,
  cluster = fish_data$StudyID,
  r = sampling_correlation,
  smooth_vi = FALSE,
  check_PD = TRUE,
  return_list = FALSE
)

fish_dominance_level_model <- metafor::rma.mv(
  yi = HedgesdDominance,
  V = fish_dominance_V,
  mods = ~ 0 + ModeratorValue,
  random = list(
    ~ 1 | StudyID,
    ~ 1 | ObservationWithinStudy
  ),
  data = fish_data,
  method = "REML",
  test = "t",
  sparse = TRUE
)

fish_dominance_contrast_model <- metafor::rma.mv(
  yi = HedgesdDominance,
  V = fish_dominance_V,
  mods = ~ ModeratorValue,
  random = list(
    ~ 1 | StudyID,
    ~ 1 | ObservationWithinStudy
  ),
  data = fish_data,
  method = "REML",
  test = "t",
  sparse = TRUE
)

fish_dominance_level_CR2 <- as.data.frame(
  clubSandwich::coef_test(
    fish_dominance_level_model,
    vcov = "CR2",
    cluster = fish_data$StudyID,
    test = "Satterthwaite"
  )
)

fish_dominance_contrast_CR2 <- as.data.frame(
  clubSandwich::coef_test(
    fish_dominance_contrast_model,
    vcov = "CR2",
    cluster = fish_data$StudyID,
    test = "Satterthwaite"
  )
)

fish_dominance_level_critical <- stats::qt(
  0.975,
  df = fish_dominance_level_CR2$df_Satt
)

fish_dominance_level_lower <-
  fish_dominance_level_CR2$beta -
  fish_dominance_level_critical * fish_dominance_level_CR2$SE

fish_dominance_level_upper <-
  fish_dominance_level_CR2$beta +
  fish_dominance_level_critical * fish_dominance_level_CR2$SE

fish_dominance_contrast_CR2 <- fish_dominance_contrast_CR2[
  -1,
  ,
  drop = FALSE
]

fish_dominance_contrast_critical <- stats::qt(
  0.975,
  df = fish_dominance_contrast_CR2$df_Satt
)

fish_dominance_contrast_lower <-
  fish_dominance_contrast_CR2$beta -
  fish_dominance_contrast_critical * fish_dominance_contrast_CR2$SE

fish_dominance_contrast_upper <-
  fish_dominance_contrast_CR2$beta +
  fish_dominance_contrast_critical * fish_dominance_contrast_CR2$SE

fish_dominance_W <- solve(fish_dominance_V)
fish_dominance_X <- as.matrix(fish_dominance_contrast_model$X)
fish_dominance_P <-
  fish_dominance_W -
  fish_dominance_W %*% fish_dominance_X %*%
  solve(crossprod(fish_dominance_X, fish_dominance_W %*% fish_dominance_X)) %*%
  t(fish_dominance_X) %*% fish_dominance_W

fish_dominance_typical_sampling_variance <-
  (nrow(fish_dominance_X) - ncol(fish_dominance_X)) /
  sum(diag(fish_dominance_P))

fish_dominance_total_variance <-
  sum(fish_dominance_contrast_model$sigma2) +
  fish_dominance_typical_sampling_variance

fish_dominance_study_I2 <-
  100 * fish_dominance_contrast_model$sigma2[1] /
  fish_dominance_total_variance

fish_dominance_observation_I2 <-
  100 * fish_dominance_contrast_model$sigma2[2] /
  fish_dominance_total_variance

fish_dominance_total_I2 <-
  100 * sum(fish_dominance_contrast_model$sigma2) /
  fish_dominance_total_variance

fish_dominance_level_results <- data.frame(
  Moderator = "Fish presence",
  NullModel = "Dominance",
  Term = fish_levels,
  EffectType = "Level mean",
  ReferenceLevel = fish_levels[1],
  N = fish_N,
  n = fish_n,
  Level_N = fish_level_N,
  Level_n = fish_level_n,
  Estimate = as.numeric(fish_dominance_level_CR2$beta),
  CR2_SE = as.numeric(fish_dominance_level_CR2$SE),
  CR2_df = as.numeric(fish_dominance_level_CR2$df_Satt),
  CI_lower = as.numeric(fish_dominance_level_lower),
  CI_upper = as.numeric(fish_dominance_level_upper),
  P_value = as.numeric(fish_dominance_level_CR2$p_Satt),
  Moderator_QM = as.numeric(fish_dominance_contrast_model$QM),
  Moderator_QM_P = as.numeric(fish_dominance_contrast_model$QMp),
  Tau2_Study = fish_dominance_contrast_model$sigma2[1],
  Tau2_Observation = fish_dominance_contrast_model$sigma2[2],
  Study_I2 = fish_dominance_study_I2,
  Observation_I2 = fish_dominance_observation_I2,
  Total_I2 = fish_dominance_total_I2,
  QE = as.numeric(fish_dominance_contrast_model$QE),
  QE_P_value = as.numeric(fish_dominance_contrast_model$QEp),
  Interpretation = ifelse(
    fish_dominance_level_lower > 0,
    "Positive non-additive",
    ifelse(
      fish_dominance_level_upper < 0,
      "Negative non-additive",
      "Additive"
    )
  ),
  Status = "Successfully fitted",
  stringsAsFactors = FALSE
)

fish_dominance_contrast_results <- data.frame(
  Moderator = "Fish presence",
  NullModel = "Dominance",
  Term = paste0(
    fish_levels[-1],
    " vs ",
    fish_levels[1]
  ),
  EffectType = "Contrast with reference",
  ReferenceLevel = fish_levels[1],
  N = fish_N,
  n = fish_n,
  Level_N = fish_N,
  Level_n = fish_n,
  Estimate = as.numeric(fish_dominance_contrast_CR2$beta),
  CR2_SE = as.numeric(fish_dominance_contrast_CR2$SE),
  CR2_df = as.numeric(fish_dominance_contrast_CR2$df_Satt),
  CI_lower = as.numeric(fish_dominance_contrast_lower),
  CI_upper = as.numeric(fish_dominance_contrast_upper),
  P_value = as.numeric(fish_dominance_contrast_CR2$p_Satt),
  Moderator_QM = as.numeric(fish_dominance_contrast_model$QM),
  Moderator_QM_P = as.numeric(fish_dominance_contrast_model$QMp),
  Tau2_Study = fish_dominance_contrast_model$sigma2[1],
  Tau2_Observation = fish_dominance_contrast_model$sigma2[2],
  Study_I2 = fish_dominance_study_I2,
  Observation_I2 = fish_dominance_observation_I2,
  Total_I2 = fish_dominance_total_I2,
  QE = as.numeric(fish_dominance_contrast_model$QE),
  QE_P_value = as.numeric(fish_dominance_contrast_model$QEp),
  Interpretation = ifelse(
    fish_dominance_contrast_lower > 0,
    "Higher than reference",
    ifelse(
      fish_dominance_contrast_upper < 0,
      "Lower than reference",
      "No detected difference"
    )
  ),
  Status = "Successfully fitted",
  stringsAsFactors = FALSE
)

fish_dominance_results <- dplyr::bind_rows(
  fish_dominance_level_results,
  fish_dominance_contrast_results
)

fish_results <- dplyr::bind_rows(
  fish_additive_results,
  fish_multiplicative_results,
  fish_dominance_results
)



# 05.2.3. Snail presence -------------------------------------------------------

snail_data <- moderator_data[
    !is.na(moderator_data$SnailPresence) &
  trimws(as.character(moderator_data$SnailPresence)) != "" &
    is.finite(moderator_data$HedgesdAdditive) &
    is.finite(moderator_data$VdAdditive) &
    moderator_data$VdAdditive > 0 &
    is.finite(moderator_data$HedgesdMultiplicative) &
    is.finite(moderator_data$VdMultiplicative) &
    moderator_data$VdMultiplicative > 0 &
    is.finite(moderator_data$HedgesdDominance) &
    is.finite(moderator_data$VdDominance) &
    moderator_data$VdDominance > 0 &
    !is.na(moderator_data$StudyID) &
    !is.na(moderator_data$ObservationID),
  ,
  drop = FALSE
]

snail_data$ModeratorValue <- droplevels(
  factor(snail_data$SnailPresence)
)

snail_data$StudyID <- factor(
  snail_data$StudyID
)

snail_data$ObservationWithinStudy <- factor(
  paste(
    snail_data$StudyID,
    snail_data$ObservationID,
    sep = "_"
  )
)

snail_levels <- levels(
  snail_data$ModeratorValue
)

snail_N <- length(
  unique(snail_data$StudyID)
)

snail_n <- nrow(
  snail_data
)

snail_level_n_table <- as.data.frame(
  table(snail_data$ModeratorValue),
  stringsAsFactors = FALSE
)

snail_level_study_data <- unique(
  data.frame(
    Level = snail_data$ModeratorValue,
    StudyID = snail_data$StudyID
  )
)

snail_level_N_table <- aggregate(
  StudyID ~ Level,
  data = snail_level_study_data,
  FUN = length
)

snail_level_n <- as.integer(
  snail_level_n_table$Freq[
    match(
      snail_levels,
      as.character(snail_level_n_table$Var1)
    )
  ]
)

snail_level_N <- as.integer(
  snail_level_N_table$StudyID[
    match(
      snail_levels,
      as.character(snail_level_N_table$Level)
    )
  ]
)

# Additive null model
snail_additive_V <- clubSandwich::impute_covariance_matrix(
  vi = snail_data$VdAdditive,
  cluster = snail_data$StudyID,
  r = sampling_correlation,
  smooth_vi = FALSE,
  check_PD = TRUE,
  return_list = FALSE
)

snail_additive_level_model <- metafor::rma.mv(
  yi = HedgesdAdditive,
  V = snail_additive_V,
  mods = ~ 0 + ModeratorValue,
  random = list(
    ~ 1 | StudyID,
    ~ 1 | ObservationWithinStudy
  ),
  data = snail_data,
  method = "REML",
  test = "t",
  sparse = TRUE
)

snail_additive_contrast_model <- metafor::rma.mv(
  yi = HedgesdAdditive,
  V = snail_additive_V,
  mods = ~ ModeratorValue,
  random = list(
    ~ 1 | StudyID,
    ~ 1 | ObservationWithinStudy
  ),
  data = snail_data,
  method = "REML",
  test = "t",
  sparse = TRUE
)

snail_additive_level_CR2 <- as.data.frame(
  clubSandwich::coef_test(
    snail_additive_level_model,
    vcov = "CR2",
    cluster = snail_data$StudyID,
    test = "Satterthwaite"
  )
)

snail_additive_contrast_CR2 <- as.data.frame(
  clubSandwich::coef_test(
    snail_additive_contrast_model,
    vcov = "CR2",
    cluster = snail_data$StudyID,
    test = "Satterthwaite"
  )
)

snail_additive_level_critical <- stats::qt(
  0.975,
  df = snail_additive_level_CR2$df_Satt
)

snail_additive_level_lower <-
  snail_additive_level_CR2$beta -
  snail_additive_level_critical * snail_additive_level_CR2$SE

snail_additive_level_upper <-
  snail_additive_level_CR2$beta +
  snail_additive_level_critical * snail_additive_level_CR2$SE

snail_additive_contrast_CR2 <- snail_additive_contrast_CR2[
  -1,
  ,
  drop = FALSE
]

snail_additive_contrast_critical <- stats::qt(
  0.975,
  df = snail_additive_contrast_CR2$df_Satt
)

snail_additive_contrast_lower <-
  snail_additive_contrast_CR2$beta -
  snail_additive_contrast_critical * snail_additive_contrast_CR2$SE

snail_additive_contrast_upper <-
  snail_additive_contrast_CR2$beta +
  snail_additive_contrast_critical * snail_additive_contrast_CR2$SE

snail_additive_W <- solve(snail_additive_V)
snail_additive_X <- as.matrix(snail_additive_contrast_model$X)
snail_additive_P <-
  snail_additive_W -
  snail_additive_W %*% snail_additive_X %*%
  solve(crossprod(snail_additive_X, snail_additive_W %*% snail_additive_X)) %*%
  t(snail_additive_X) %*% snail_additive_W

snail_additive_typical_sampling_variance <-
  (nrow(snail_additive_X) - ncol(snail_additive_X)) /
  sum(diag(snail_additive_P))

snail_additive_total_variance <-
  sum(snail_additive_contrast_model$sigma2) +
  snail_additive_typical_sampling_variance

snail_additive_study_I2 <-
  100 * snail_additive_contrast_model$sigma2[1] /
  snail_additive_total_variance

snail_additive_observation_I2 <-
  100 * snail_additive_contrast_model$sigma2[2] /
  snail_additive_total_variance

snail_additive_total_I2 <-
  100 * sum(snail_additive_contrast_model$sigma2) /
  snail_additive_total_variance

snail_additive_level_results <- data.frame(
  Moderator = "Snail presence",
  NullModel = "Additive",
  Term = snail_levels,
  EffectType = "Level mean",
  ReferenceLevel = snail_levels[1],
  N = snail_N,
  n = snail_n,
  Level_N = snail_level_N,
  Level_n = snail_level_n,
  Estimate = as.numeric(snail_additive_level_CR2$beta),
  CR2_SE = as.numeric(snail_additive_level_CR2$SE),
  CR2_df = as.numeric(snail_additive_level_CR2$df_Satt),
  CI_lower = as.numeric(snail_additive_level_lower),
  CI_upper = as.numeric(snail_additive_level_upper),
  P_value = as.numeric(snail_additive_level_CR2$p_Satt),
  Moderator_QM = as.numeric(snail_additive_contrast_model$QM),
  Moderator_QM_P = as.numeric(snail_additive_contrast_model$QMp),
  Tau2_Study = snail_additive_contrast_model$sigma2[1],
  Tau2_Observation = snail_additive_contrast_model$sigma2[2],
  Study_I2 = snail_additive_study_I2,
  Observation_I2 = snail_additive_observation_I2,
  Total_I2 = snail_additive_total_I2,
  QE = as.numeric(snail_additive_contrast_model$QE),
  QE_P_value = as.numeric(snail_additive_contrast_model$QEp),
  Interpretation = ifelse(
    snail_additive_level_lower > 0,
    "Positive non-additive",
    ifelse(
      snail_additive_level_upper < 0,
      "Negative non-additive",
      "Additive"
    )
  ),
  Status = "Successfully fitted",
  stringsAsFactors = FALSE
)

snail_additive_contrast_results <- data.frame(
  Moderator = "Snail presence",
  NullModel = "Additive",
  Term = paste0(
    snail_levels[-1],
    " vs ",
    snail_levels[1]
  ),
  EffectType = "Contrast with reference",
  ReferenceLevel = snail_levels[1],
  N = snail_N,
  n = snail_n,
  Level_N = snail_N,
  Level_n = snail_n,
  Estimate = as.numeric(snail_additive_contrast_CR2$beta),
  CR2_SE = as.numeric(snail_additive_contrast_CR2$SE),
  CR2_df = as.numeric(snail_additive_contrast_CR2$df_Satt),
  CI_lower = as.numeric(snail_additive_contrast_lower),
  CI_upper = as.numeric(snail_additive_contrast_upper),
  P_value = as.numeric(snail_additive_contrast_CR2$p_Satt),
  Moderator_QM = as.numeric(snail_additive_contrast_model$QM),
  Moderator_QM_P = as.numeric(snail_additive_contrast_model$QMp),
  Tau2_Study = snail_additive_contrast_model$sigma2[1],
  Tau2_Observation = snail_additive_contrast_model$sigma2[2],
  Study_I2 = snail_additive_study_I2,
  Observation_I2 = snail_additive_observation_I2,
  Total_I2 = snail_additive_total_I2,
  QE = as.numeric(snail_additive_contrast_model$QE),
  QE_P_value = as.numeric(snail_additive_contrast_model$QEp),
  Interpretation = ifelse(
    snail_additive_contrast_lower > 0,
    "Higher than reference",
    ifelse(
      snail_additive_contrast_upper < 0,
      "Lower than reference",
      "No detected difference"
    )
  ),
  Status = "Successfully fitted",
  stringsAsFactors = FALSE
)

snail_additive_results <- dplyr::bind_rows(
  snail_additive_level_results,
  snail_additive_contrast_results
)


# Multiplicative null model
snail_multiplicative_V <- clubSandwich::impute_covariance_matrix(
  vi = snail_data$VdMultiplicative,
  cluster = snail_data$StudyID,
  r = sampling_correlation,
  smooth_vi = FALSE,
  check_PD = TRUE,
  return_list = FALSE
)

snail_multiplicative_level_model <- metafor::rma.mv(
  yi = HedgesdMultiplicative,
  V = snail_multiplicative_V,
  mods = ~ 0 + ModeratorValue,
  random = list(
    ~ 1 | StudyID,
    ~ 1 | ObservationWithinStudy
  ),
  data = snail_data,
  method = "REML",
  test = "t",
  sparse = TRUE
)

snail_multiplicative_contrast_model <- metafor::rma.mv(
  yi = HedgesdMultiplicative,
  V = snail_multiplicative_V,
  mods = ~ ModeratorValue,
  random = list(
    ~ 1 | StudyID,
    ~ 1 | ObservationWithinStudy
  ),
  data = snail_data,
  method = "REML",
  test = "t",
  sparse = TRUE
)

snail_multiplicative_level_CR2 <- as.data.frame(
  clubSandwich::coef_test(
    snail_multiplicative_level_model,
    vcov = "CR2",
    cluster = snail_data$StudyID,
    test = "Satterthwaite"
  )
)

snail_multiplicative_contrast_CR2 <- as.data.frame(
  clubSandwich::coef_test(
    snail_multiplicative_contrast_model,
    vcov = "CR2",
    cluster = snail_data$StudyID,
    test = "Satterthwaite"
  )
)

snail_multiplicative_level_critical <- stats::qt(
  0.975,
  df = snail_multiplicative_level_CR2$df_Satt
)

snail_multiplicative_level_lower <-
  snail_multiplicative_level_CR2$beta -
  snail_multiplicative_level_critical * snail_multiplicative_level_CR2$SE

snail_multiplicative_level_upper <-
  snail_multiplicative_level_CR2$beta +
  snail_multiplicative_level_critical * snail_multiplicative_level_CR2$SE

snail_multiplicative_contrast_CR2 <- snail_multiplicative_contrast_CR2[
  -1,
  ,
  drop = FALSE
]

snail_multiplicative_contrast_critical <- stats::qt(
  0.975,
  df = snail_multiplicative_contrast_CR2$df_Satt
)

snail_multiplicative_contrast_lower <-
  snail_multiplicative_contrast_CR2$beta -
  snail_multiplicative_contrast_critical * snail_multiplicative_contrast_CR2$SE

snail_multiplicative_contrast_upper <-
  snail_multiplicative_contrast_CR2$beta +
  snail_multiplicative_contrast_critical * snail_multiplicative_contrast_CR2$SE

snail_multiplicative_W <- solve(snail_multiplicative_V)
snail_multiplicative_X <- as.matrix(snail_multiplicative_contrast_model$X)
snail_multiplicative_P <-
  snail_multiplicative_W -
  snail_multiplicative_W %*% snail_multiplicative_X %*%
  solve(crossprod(snail_multiplicative_X, snail_multiplicative_W %*% snail_multiplicative_X)) %*%
  t(snail_multiplicative_X) %*% snail_multiplicative_W

snail_multiplicative_typical_sampling_variance <-
  (nrow(snail_multiplicative_X) - ncol(snail_multiplicative_X)) /
  sum(diag(snail_multiplicative_P))

snail_multiplicative_total_variance <-
  sum(snail_multiplicative_contrast_model$sigma2) +
  snail_multiplicative_typical_sampling_variance

snail_multiplicative_study_I2 <-
  100 * snail_multiplicative_contrast_model$sigma2[1] /
  snail_multiplicative_total_variance

snail_multiplicative_observation_I2 <-
  100 * snail_multiplicative_contrast_model$sigma2[2] /
  snail_multiplicative_total_variance

snail_multiplicative_total_I2 <-
  100 * sum(snail_multiplicative_contrast_model$sigma2) /
  snail_multiplicative_total_variance

snail_multiplicative_level_results <- data.frame(
  Moderator = "Snail presence",
  NullModel = "Multiplicative",
  Term = snail_levels,
  EffectType = "Level mean",
  ReferenceLevel = snail_levels[1],
  N = snail_N,
  n = snail_n,
  Level_N = snail_level_N,
  Level_n = snail_level_n,
  Estimate = as.numeric(snail_multiplicative_level_CR2$beta),
  CR2_SE = as.numeric(snail_multiplicative_level_CR2$SE),
  CR2_df = as.numeric(snail_multiplicative_level_CR2$df_Satt),
  CI_lower = as.numeric(snail_multiplicative_level_lower),
  CI_upper = as.numeric(snail_multiplicative_level_upper),
  P_value = as.numeric(snail_multiplicative_level_CR2$p_Satt),
  Moderator_QM = as.numeric(snail_multiplicative_contrast_model$QM),
  Moderator_QM_P = as.numeric(snail_multiplicative_contrast_model$QMp),
  Tau2_Study = snail_multiplicative_contrast_model$sigma2[1],
  Tau2_Observation = snail_multiplicative_contrast_model$sigma2[2],
  Study_I2 = snail_multiplicative_study_I2,
  Observation_I2 = snail_multiplicative_observation_I2,
  Total_I2 = snail_multiplicative_total_I2,
  QE = as.numeric(snail_multiplicative_contrast_model$QE),
  QE_P_value = as.numeric(snail_multiplicative_contrast_model$QEp),
  Interpretation = ifelse(
    snail_multiplicative_level_lower > 0,
    "Positive non-additive",
    ifelse(
      snail_multiplicative_level_upper < 0,
      "Negative non-additive",
      "Additive"
    )
  ),
  Status = "Successfully fitted",
  stringsAsFactors = FALSE
)

snail_multiplicative_contrast_results <- data.frame(
  Moderator = "Snail presence",
  NullModel = "Multiplicative",
  Term = paste0(
    snail_levels[-1],
    " vs ",
    snail_levels[1]
  ),
  EffectType = "Contrast with reference",
  ReferenceLevel = snail_levels[1],
  N = snail_N,
  n = snail_n,
  Level_N = snail_N,
  Level_n = snail_n,
  Estimate = as.numeric(snail_multiplicative_contrast_CR2$beta),
  CR2_SE = as.numeric(snail_multiplicative_contrast_CR2$SE),
  CR2_df = as.numeric(snail_multiplicative_contrast_CR2$df_Satt),
  CI_lower = as.numeric(snail_multiplicative_contrast_lower),
  CI_upper = as.numeric(snail_multiplicative_contrast_upper),
  P_value = as.numeric(snail_multiplicative_contrast_CR2$p_Satt),
  Moderator_QM = as.numeric(snail_multiplicative_contrast_model$QM),
  Moderator_QM_P = as.numeric(snail_multiplicative_contrast_model$QMp),
  Tau2_Study = snail_multiplicative_contrast_model$sigma2[1],
  Tau2_Observation = snail_multiplicative_contrast_model$sigma2[2],
  Study_I2 = snail_multiplicative_study_I2,
  Observation_I2 = snail_multiplicative_observation_I2,
  Total_I2 = snail_multiplicative_total_I2,
  QE = as.numeric(snail_multiplicative_contrast_model$QE),
  QE_P_value = as.numeric(snail_multiplicative_contrast_model$QEp),
  Interpretation = ifelse(
    snail_multiplicative_contrast_lower > 0,
    "Higher than reference",
    ifelse(
      snail_multiplicative_contrast_upper < 0,
      "Lower than reference",
      "No detected difference"
    )
  ),
  Status = "Successfully fitted",
  stringsAsFactors = FALSE
)

snail_multiplicative_results <- dplyr::bind_rows(
  snail_multiplicative_level_results,
  snail_multiplicative_contrast_results
)


# Dominance null model
snail_dominance_V <- clubSandwich::impute_covariance_matrix(
  vi = snail_data$VdDominance,
  cluster = snail_data$StudyID,
  r = sampling_correlation,
  smooth_vi = FALSE,
  check_PD = TRUE,
  return_list = FALSE
)

snail_dominance_level_model <- metafor::rma.mv(
  yi = HedgesdDominance,
  V = snail_dominance_V,
  mods = ~ 0 + ModeratorValue,
  random = list(
    ~ 1 | StudyID,
    ~ 1 | ObservationWithinStudy
  ),
  data = snail_data,
  method = "REML",
  test = "t",
  sparse = TRUE
)

snail_dominance_contrast_model <- metafor::rma.mv(
  yi = HedgesdDominance,
  V = snail_dominance_V,
  mods = ~ ModeratorValue,
  random = list(
    ~ 1 | StudyID,
    ~ 1 | ObservationWithinStudy
  ),
  data = snail_data,
  method = "REML",
  test = "t",
  sparse = TRUE
)

snail_dominance_level_CR2 <- as.data.frame(
  clubSandwich::coef_test(
    snail_dominance_level_model,
    vcov = "CR2",
    cluster = snail_data$StudyID,
    test = "Satterthwaite"
  )
)

snail_dominance_contrast_CR2 <- as.data.frame(
  clubSandwich::coef_test(
    snail_dominance_contrast_model,
    vcov = "CR2",
    cluster = snail_data$StudyID,
    test = "Satterthwaite"
  )
)

snail_dominance_level_critical <- stats::qt(
  0.975,
  df = snail_dominance_level_CR2$df_Satt
)

snail_dominance_level_lower <-
  snail_dominance_level_CR2$beta -
  snail_dominance_level_critical * snail_dominance_level_CR2$SE

snail_dominance_level_upper <-
  snail_dominance_level_CR2$beta +
  snail_dominance_level_critical * snail_dominance_level_CR2$SE

snail_dominance_contrast_CR2 <- snail_dominance_contrast_CR2[
  -1,
  ,
  drop = FALSE
]

snail_dominance_contrast_critical <- stats::qt(
  0.975,
  df = snail_dominance_contrast_CR2$df_Satt
)

snail_dominance_contrast_lower <-
  snail_dominance_contrast_CR2$beta -
  snail_dominance_contrast_critical * snail_dominance_contrast_CR2$SE

snail_dominance_contrast_upper <-
  snail_dominance_contrast_CR2$beta +
  snail_dominance_contrast_critical * snail_dominance_contrast_CR2$SE

snail_dominance_W <- solve(snail_dominance_V)
snail_dominance_X <- as.matrix(snail_dominance_contrast_model$X)
snail_dominance_P <-
  snail_dominance_W -
  snail_dominance_W %*% snail_dominance_X %*%
  solve(crossprod(snail_dominance_X, snail_dominance_W %*% snail_dominance_X)) %*%
  t(snail_dominance_X) %*% snail_dominance_W

snail_dominance_typical_sampling_variance <-
  (nrow(snail_dominance_X) - ncol(snail_dominance_X)) /
  sum(diag(snail_dominance_P))

snail_dominance_total_variance <-
  sum(snail_dominance_contrast_model$sigma2) +
  snail_dominance_typical_sampling_variance

snail_dominance_study_I2 <-
  100 * snail_dominance_contrast_model$sigma2[1] /
  snail_dominance_total_variance

snail_dominance_observation_I2 <-
  100 * snail_dominance_contrast_model$sigma2[2] /
  snail_dominance_total_variance

snail_dominance_total_I2 <-
  100 * sum(snail_dominance_contrast_model$sigma2) /
  snail_dominance_total_variance

snail_dominance_level_results <- data.frame(
  Moderator = "Snail presence",
  NullModel = "Dominance",
  Term = snail_levels,
  EffectType = "Level mean",
  ReferenceLevel = snail_levels[1],
  N = snail_N,
  n = snail_n,
  Level_N = snail_level_N,
  Level_n = snail_level_n,
  Estimate = as.numeric(snail_dominance_level_CR2$beta),
  CR2_SE = as.numeric(snail_dominance_level_CR2$SE),
  CR2_df = as.numeric(snail_dominance_level_CR2$df_Satt),
  CI_lower = as.numeric(snail_dominance_level_lower),
  CI_upper = as.numeric(snail_dominance_level_upper),
  P_value = as.numeric(snail_dominance_level_CR2$p_Satt),
  Moderator_QM = as.numeric(snail_dominance_contrast_model$QM),
  Moderator_QM_P = as.numeric(snail_dominance_contrast_model$QMp),
  Tau2_Study = snail_dominance_contrast_model$sigma2[1],
  Tau2_Observation = snail_dominance_contrast_model$sigma2[2],
  Study_I2 = snail_dominance_study_I2,
  Observation_I2 = snail_dominance_observation_I2,
  Total_I2 = snail_dominance_total_I2,
  QE = as.numeric(snail_dominance_contrast_model$QE),
  QE_P_value = as.numeric(snail_dominance_contrast_model$QEp),
  Interpretation = ifelse(
    snail_dominance_level_lower > 0,
    "Positive non-additive",
    ifelse(
      snail_dominance_level_upper < 0,
      "Negative non-additive",
      "Additive"
    )
  ),
  Status = "Successfully fitted",
  stringsAsFactors = FALSE
)

snail_dominance_contrast_results <- data.frame(
  Moderator = "Snail presence",
  NullModel = "Dominance",
  Term = paste0(
    snail_levels[-1],
    " vs ",
    snail_levels[1]
  ),
  EffectType = "Contrast with reference",
  ReferenceLevel = snail_levels[1],
  N = snail_N,
  n = snail_n,
  Level_N = snail_N,
  Level_n = snail_n,
  Estimate = as.numeric(snail_dominance_contrast_CR2$beta),
  CR2_SE = as.numeric(snail_dominance_contrast_CR2$SE),
  CR2_df = as.numeric(snail_dominance_contrast_CR2$df_Satt),
  CI_lower = as.numeric(snail_dominance_contrast_lower),
  CI_upper = as.numeric(snail_dominance_contrast_upper),
  P_value = as.numeric(snail_dominance_contrast_CR2$p_Satt),
  Moderator_QM = as.numeric(snail_dominance_contrast_model$QM),
  Moderator_QM_P = as.numeric(snail_dominance_contrast_model$QMp),
  Tau2_Study = snail_dominance_contrast_model$sigma2[1],
  Tau2_Observation = snail_dominance_contrast_model$sigma2[2],
  Study_I2 = snail_dominance_study_I2,
  Observation_I2 = snail_dominance_observation_I2,
  Total_I2 = snail_dominance_total_I2,
  QE = as.numeric(snail_dominance_contrast_model$QE),
  QE_P_value = as.numeric(snail_dominance_contrast_model$QEp),
  Interpretation = ifelse(
    snail_dominance_contrast_lower > 0,
    "Higher than reference",
    ifelse(
      snail_dominance_contrast_upper < 0,
      "Lower than reference",
      "No detected difference"
    )
  ),
  Status = "Successfully fitted",
  stringsAsFactors = FALSE
)

snail_dominance_results <- dplyr::bind_rows(
  snail_dominance_level_results,
  snail_dominance_contrast_results
)

snail_results <- dplyr::bind_rows(
  snail_additive_results,
  snail_multiplicative_results,
  snail_dominance_results
)




# 05.2.4. Ambient temperature -------------------------------------------------

ambient_temperature_raw <- suppressWarnings(
  as.numeric(moderator_data$AmbientTemperature)
)

ambient_temperature_data <- moderator_data[
  is.finite(ambient_temperature_raw) &
    is.finite(moderator_data$HedgesdAdditive) &
    is.finite(moderator_data$VdAdditive) &
    moderator_data$VdAdditive > 0 &
    is.finite(moderator_data$HedgesdMultiplicative) &
    is.finite(moderator_data$VdMultiplicative) &
    moderator_data$VdMultiplicative > 0 &
    is.finite(moderator_data$HedgesdDominance) &
    is.finite(moderator_data$VdDominance) &
    moderator_data$VdDominance > 0 &
    !is.na(moderator_data$StudyID) &
    !is.na(moderator_data$ObservationID),
  ,
  drop = FALSE
]

ambient_temperature_data$ModeratorValue <- as.numeric(
  scale(
    suppressWarnings(
      as.numeric(ambient_temperature_data$AmbientTemperature)
    )
  )
)

ambient_temperature_data$StudyID <- factor(
  ambient_temperature_data$StudyID
)

ambient_temperature_data$ObservationWithinStudy <- factor(
  paste(
    ambient_temperature_data$StudyID,
    ambient_temperature_data$ObservationID,
    sep = "_"
  )
)

ambient_temperature_N <- length(
  unique(ambient_temperature_data$StudyID)
)

ambient_temperature_n <- nrow(
  ambient_temperature_data
)

# Additive null model
ambient_temperature_additive_V <- clubSandwich::impute_covariance_matrix(
  vi = ambient_temperature_data$VdAdditive,
  cluster = ambient_temperature_data$StudyID,
  r = sampling_correlation,
  smooth_vi = FALSE,
  check_PD = TRUE,
  return_list = FALSE
)

ambient_temperature_additive_model <- metafor::rma.mv(
  yi = HedgesdAdditive,
  V = ambient_temperature_additive_V,
  mods = ~ ModeratorValue,
  random = list(
    ~ 1 | StudyID,
    ~ 1 | ObservationWithinStudy
  ),
  data = ambient_temperature_data,
  method = "REML",
  test = "t",
  sparse = TRUE
)

ambient_temperature_additive_CR2 <- as.data.frame(
  clubSandwich::coef_test(
    ambient_temperature_additive_model,
    vcov = "CR2",
    cluster = ambient_temperature_data$StudyID,
    test = "Satterthwaite"
  )
)

ambient_temperature_additive_critical <- stats::qt(
  0.975,
  df = ambient_temperature_additive_CR2$df_Satt
)

ambient_temperature_additive_lower <-
  ambient_temperature_additive_CR2$beta -
  ambient_temperature_additive_critical * ambient_temperature_additive_CR2$SE

ambient_temperature_additive_upper <-
  ambient_temperature_additive_CR2$beta +
  ambient_temperature_additive_critical * ambient_temperature_additive_CR2$SE

ambient_temperature_additive_W <- solve(ambient_temperature_additive_V)
ambient_temperature_additive_X <- as.matrix(ambient_temperature_additive_model$X)
ambient_temperature_additive_P <-
  ambient_temperature_additive_W -
  ambient_temperature_additive_W %*% ambient_temperature_additive_X %*%
  solve(crossprod(ambient_temperature_additive_X, ambient_temperature_additive_W %*% ambient_temperature_additive_X)) %*%
  t(ambient_temperature_additive_X) %*% ambient_temperature_additive_W

ambient_temperature_additive_typical_sampling_variance <-
  (nrow(ambient_temperature_additive_X) - ncol(ambient_temperature_additive_X)) /
  sum(diag(ambient_temperature_additive_P))

ambient_temperature_additive_total_variance <-
  sum(ambient_temperature_additive_model$sigma2) +
  ambient_temperature_additive_typical_sampling_variance

ambient_temperature_additive_study_I2 <-
  100 * ambient_temperature_additive_model$sigma2[1] /
  ambient_temperature_additive_total_variance

ambient_temperature_additive_observation_I2 <-
  100 * ambient_temperature_additive_model$sigma2[2] /
  ambient_temperature_additive_total_variance

ambient_temperature_additive_total_I2 <-
  100 * sum(ambient_temperature_additive_model$sigma2) /
  ambient_temperature_additive_total_variance

ambient_temperature_additive_results <- data.frame(
  Moderator = "Ambient temperature",
  NullModel = "Additive",
  Term = c(
    "Intercept at mean",
    "Slope per SD"
  ),
  EffectType = c(
    "Intercept at mean",
    "Slope per SD"
  ),
  ReferenceLevel = NA_character_,
  N = ambient_temperature_N,
  n = ambient_temperature_n,
  Level_N = ambient_temperature_N,
  Level_n = ambient_temperature_n,
  Estimate = as.numeric(ambient_temperature_additive_CR2$beta),
  CR2_SE = as.numeric(ambient_temperature_additive_CR2$SE),
  CR2_df = as.numeric(ambient_temperature_additive_CR2$df_Satt),
  CI_lower = as.numeric(ambient_temperature_additive_lower),
  CI_upper = as.numeric(ambient_temperature_additive_upper),
  P_value = as.numeric(ambient_temperature_additive_CR2$p_Satt),
  Moderator_QM = as.numeric(ambient_temperature_additive_model$QM),
  Moderator_QM_P = as.numeric(ambient_temperature_additive_model$QMp),
  Tau2_Study = ambient_temperature_additive_model$sigma2[1],
  Tau2_Observation = ambient_temperature_additive_model$sigma2[2],
  Study_I2 = ambient_temperature_additive_study_I2,
  Observation_I2 = ambient_temperature_additive_observation_I2,
  Total_I2 = ambient_temperature_additive_total_I2,
  QE = as.numeric(ambient_temperature_additive_model$QE),
  QE_P_value = as.numeric(ambient_temperature_additive_model$QEp),
  Interpretation = c(
    ifelse(
      ambient_temperature_additive_lower[1] > 0,
      "Positive non-additive",
      ifelse(
        ambient_temperature_additive_upper[1] < 0,
        "Negative non-additive",
        "Additive"
      )
    ),
    ifelse(
      ambient_temperature_additive_lower[2] > 0,
      "Positive association",
      ifelse(
        ambient_temperature_additive_upper[2] < 0,
        "Negative association",
        "No detected association"
      )
    )
  ),
  Status = "Successfully fitted",
  stringsAsFactors = FALSE
)


# Multiplicative null model
ambient_temperature_multiplicative_V <- clubSandwich::impute_covariance_matrix(
  vi = ambient_temperature_data$VdMultiplicative,
  cluster = ambient_temperature_data$StudyID,
  r = sampling_correlation,
  smooth_vi = FALSE,
  check_PD = TRUE,
  return_list = FALSE
)

ambient_temperature_multiplicative_model <- metafor::rma.mv(
  yi = HedgesdMultiplicative,
  V = ambient_temperature_multiplicative_V,
  mods = ~ ModeratorValue,
  random = list(
    ~ 1 | StudyID,
    ~ 1 | ObservationWithinStudy
  ),
  data = ambient_temperature_data,
  method = "REML",
  test = "t",
  sparse = TRUE
)

ambient_temperature_multiplicative_CR2 <- as.data.frame(
  clubSandwich::coef_test(
    ambient_temperature_multiplicative_model,
    vcov = "CR2",
    cluster = ambient_temperature_data$StudyID,
    test = "Satterthwaite"
  )
)

ambient_temperature_multiplicative_critical <- stats::qt(
  0.975,
  df = ambient_temperature_multiplicative_CR2$df_Satt
)

ambient_temperature_multiplicative_lower <-
  ambient_temperature_multiplicative_CR2$beta -
  ambient_temperature_multiplicative_critical * ambient_temperature_multiplicative_CR2$SE

ambient_temperature_multiplicative_upper <-
  ambient_temperature_multiplicative_CR2$beta +
  ambient_temperature_multiplicative_critical * ambient_temperature_multiplicative_CR2$SE

ambient_temperature_multiplicative_W <- solve(ambient_temperature_multiplicative_V)
ambient_temperature_multiplicative_X <- as.matrix(ambient_temperature_multiplicative_model$X)
ambient_temperature_multiplicative_P <-
  ambient_temperature_multiplicative_W -
  ambient_temperature_multiplicative_W %*% ambient_temperature_multiplicative_X %*%
  solve(crossprod(ambient_temperature_multiplicative_X, ambient_temperature_multiplicative_W %*% ambient_temperature_multiplicative_X)) %*%
  t(ambient_temperature_multiplicative_X) %*% ambient_temperature_multiplicative_W

ambient_temperature_multiplicative_typical_sampling_variance <-
  (nrow(ambient_temperature_multiplicative_X) - ncol(ambient_temperature_multiplicative_X)) /
  sum(diag(ambient_temperature_multiplicative_P))

ambient_temperature_multiplicative_total_variance <-
  sum(ambient_temperature_multiplicative_model$sigma2) +
  ambient_temperature_multiplicative_typical_sampling_variance

ambient_temperature_multiplicative_study_I2 <-
  100 * ambient_temperature_multiplicative_model$sigma2[1] /
  ambient_temperature_multiplicative_total_variance

ambient_temperature_multiplicative_observation_I2 <-
  100 * ambient_temperature_multiplicative_model$sigma2[2] /
  ambient_temperature_multiplicative_total_variance

ambient_temperature_multiplicative_total_I2 <-
  100 * sum(ambient_temperature_multiplicative_model$sigma2) /
  ambient_temperature_multiplicative_total_variance

ambient_temperature_multiplicative_results <- data.frame(
  Moderator = "Ambient temperature",
  NullModel = "Multiplicative",
  Term = c(
    "Intercept at mean",
    "Slope per SD"
  ),
  EffectType = c(
    "Intercept at mean",
    "Slope per SD"
  ),
  ReferenceLevel = NA_character_,
  N = ambient_temperature_N,
  n = ambient_temperature_n,
  Level_N = ambient_temperature_N,
  Level_n = ambient_temperature_n,
  Estimate = as.numeric(ambient_temperature_multiplicative_CR2$beta),
  CR2_SE = as.numeric(ambient_temperature_multiplicative_CR2$SE),
  CR2_df = as.numeric(ambient_temperature_multiplicative_CR2$df_Satt),
  CI_lower = as.numeric(ambient_temperature_multiplicative_lower),
  CI_upper = as.numeric(ambient_temperature_multiplicative_upper),
  P_value = as.numeric(ambient_temperature_multiplicative_CR2$p_Satt),
  Moderator_QM = as.numeric(ambient_temperature_multiplicative_model$QM),
  Moderator_QM_P = as.numeric(ambient_temperature_multiplicative_model$QMp),
  Tau2_Study = ambient_temperature_multiplicative_model$sigma2[1],
  Tau2_Observation = ambient_temperature_multiplicative_model$sigma2[2],
  Study_I2 = ambient_temperature_multiplicative_study_I2,
  Observation_I2 = ambient_temperature_multiplicative_observation_I2,
  Total_I2 = ambient_temperature_multiplicative_total_I2,
  QE = as.numeric(ambient_temperature_multiplicative_model$QE),
  QE_P_value = as.numeric(ambient_temperature_multiplicative_model$QEp),
  Interpretation = c(
    ifelse(
      ambient_temperature_multiplicative_lower[1] > 0,
      "Positive non-additive",
      ifelse(
        ambient_temperature_multiplicative_upper[1] < 0,
        "Negative non-additive",
        "Additive"
      )
    ),
    ifelse(
      ambient_temperature_multiplicative_lower[2] > 0,
      "Positive association",
      ifelse(
        ambient_temperature_multiplicative_upper[2] < 0,
        "Negative association",
        "No detected association"
      )
    )
  ),
  Status = "Successfully fitted",
  stringsAsFactors = FALSE
)


# Dominance null model
ambient_temperature_dominance_V <- clubSandwich::impute_covariance_matrix(
  vi = ambient_temperature_data$VdDominance,
  cluster = ambient_temperature_data$StudyID,
  r = sampling_correlation,
  smooth_vi = FALSE,
  check_PD = TRUE,
  return_list = FALSE
)

ambient_temperature_dominance_model <- metafor::rma.mv(
  yi = HedgesdDominance,
  V = ambient_temperature_dominance_V,
  mods = ~ ModeratorValue,
  random = list(
    ~ 1 | StudyID,
    ~ 1 | ObservationWithinStudy
  ),
  data = ambient_temperature_data,
  method = "REML",
  test = "t",
  sparse = TRUE
)

ambient_temperature_dominance_CR2 <- as.data.frame(
  clubSandwich::coef_test(
    ambient_temperature_dominance_model,
    vcov = "CR2",
    cluster = ambient_temperature_data$StudyID,
    test = "Satterthwaite"
  )
)

ambient_temperature_dominance_critical <- stats::qt(
  0.975,
  df = ambient_temperature_dominance_CR2$df_Satt
)

ambient_temperature_dominance_lower <-
  ambient_temperature_dominance_CR2$beta -
  ambient_temperature_dominance_critical * ambient_temperature_dominance_CR2$SE

ambient_temperature_dominance_upper <-
  ambient_temperature_dominance_CR2$beta +
  ambient_temperature_dominance_critical * ambient_temperature_dominance_CR2$SE

ambient_temperature_dominance_W <- solve(ambient_temperature_dominance_V)
ambient_temperature_dominance_X <- as.matrix(ambient_temperature_dominance_model$X)
ambient_temperature_dominance_P <-
  ambient_temperature_dominance_W -
  ambient_temperature_dominance_W %*% ambient_temperature_dominance_X %*%
  solve(crossprod(ambient_temperature_dominance_X, ambient_temperature_dominance_W %*% ambient_temperature_dominance_X)) %*%
  t(ambient_temperature_dominance_X) %*% ambient_temperature_dominance_W

ambient_temperature_dominance_typical_sampling_variance <-
  (nrow(ambient_temperature_dominance_X) - ncol(ambient_temperature_dominance_X)) /
  sum(diag(ambient_temperature_dominance_P))

ambient_temperature_dominance_total_variance <-
  sum(ambient_temperature_dominance_model$sigma2) +
  ambient_temperature_dominance_typical_sampling_variance

ambient_temperature_dominance_study_I2 <-
  100 * ambient_temperature_dominance_model$sigma2[1] /
  ambient_temperature_dominance_total_variance

ambient_temperature_dominance_observation_I2 <-
  100 * ambient_temperature_dominance_model$sigma2[2] /
  ambient_temperature_dominance_total_variance

ambient_temperature_dominance_total_I2 <-
  100 * sum(ambient_temperature_dominance_model$sigma2) /
  ambient_temperature_dominance_total_variance

ambient_temperature_dominance_results <- data.frame(
  Moderator = "Ambient temperature",
  NullModel = "Dominance",
  Term = c(
    "Intercept at mean",
    "Slope per SD"
  ),
  EffectType = c(
    "Intercept at mean",
    "Slope per SD"
  ),
  ReferenceLevel = NA_character_,
  N = ambient_temperature_N,
  n = ambient_temperature_n,
  Level_N = ambient_temperature_N,
  Level_n = ambient_temperature_n,
  Estimate = as.numeric(ambient_temperature_dominance_CR2$beta),
  CR2_SE = as.numeric(ambient_temperature_dominance_CR2$SE),
  CR2_df = as.numeric(ambient_temperature_dominance_CR2$df_Satt),
  CI_lower = as.numeric(ambient_temperature_dominance_lower),
  CI_upper = as.numeric(ambient_temperature_dominance_upper),
  P_value = as.numeric(ambient_temperature_dominance_CR2$p_Satt),
  Moderator_QM = as.numeric(ambient_temperature_dominance_model$QM),
  Moderator_QM_P = as.numeric(ambient_temperature_dominance_model$QMp),
  Tau2_Study = ambient_temperature_dominance_model$sigma2[1],
  Tau2_Observation = ambient_temperature_dominance_model$sigma2[2],
  Study_I2 = ambient_temperature_dominance_study_I2,
  Observation_I2 = ambient_temperature_dominance_observation_I2,
  Total_I2 = ambient_temperature_dominance_total_I2,
  QE = as.numeric(ambient_temperature_dominance_model$QE),
  QE_P_value = as.numeric(ambient_temperature_dominance_model$QEp),
  Interpretation = c(
    ifelse(
      ambient_temperature_dominance_lower[1] > 0,
      "Positive non-additive",
      ifelse(
        ambient_temperature_dominance_upper[1] < 0,
        "Negative non-additive",
        "Additive"
      )
    ),
    ifelse(
      ambient_temperature_dominance_lower[2] > 0,
      "Positive association",
      ifelse(
        ambient_temperature_dominance_upper[2] < 0,
        "Negative association",
        "No detected association"
      )
    )
  ),
  Status = "Successfully fitted",
  stringsAsFactors = FALSE
)

ambient_temperature_results <- dplyr::bind_rows(
  ambient_temperature_additive_results,
  ambient_temperature_multiplicative_results,
  ambient_temperature_dominance_results
)




# 05.2.5. Study duration -------------------------------------------------------

study_duration_raw <- suppressWarnings(
  as.numeric(moderator_data$StudyDurationDays)
)

study_duration_data <- moderator_data[
  is.finite(study_duration_raw) &
    is.finite(moderator_data$HedgesdAdditive) &
    is.finite(moderator_data$VdAdditive) &
    moderator_data$VdAdditive > 0 &
    is.finite(moderator_data$HedgesdMultiplicative) &
    is.finite(moderator_data$VdMultiplicative) &
    moderator_data$VdMultiplicative > 0 &
    is.finite(moderator_data$HedgesdDominance) &
    is.finite(moderator_data$VdDominance) &
    moderator_data$VdDominance > 0 &
    !is.na(moderator_data$StudyID) &
    !is.na(moderator_data$ObservationID),
  ,
  drop = FALSE
]

study_duration_data$ModeratorValue <- as.numeric(
  scale(
    suppressWarnings(
      as.numeric(study_duration_data$StudyDurationDays)
    )
  )
)

study_duration_data$StudyID <- factor(
  study_duration_data$StudyID
)

study_duration_data$ObservationWithinStudy <- factor(
  paste(
    study_duration_data$StudyID,
    study_duration_data$ObservationID,
    sep = "_"
  )
)

study_duration_N <- length(
  unique(study_duration_data$StudyID)
)

study_duration_n <- nrow(
  study_duration_data
)

# Additive null model
study_duration_additive_V <- clubSandwich::impute_covariance_matrix(
  vi = study_duration_data$VdAdditive,
  cluster = study_duration_data$StudyID,
  r = sampling_correlation,
  smooth_vi = FALSE,
  check_PD = TRUE,
  return_list = FALSE
)

study_duration_additive_model <- metafor::rma.mv(
  yi = HedgesdAdditive,
  V = study_duration_additive_V,
  mods = ~ ModeratorValue,
  random = list(
    ~ 1 | StudyID,
    ~ 1 | ObservationWithinStudy
  ),
  data = study_duration_data,
  method = "REML",
  test = "t",
  sparse = TRUE
)

study_duration_additive_CR2 <- as.data.frame(
  clubSandwich::coef_test(
    study_duration_additive_model,
    vcov = "CR2",
    cluster = study_duration_data$StudyID,
    test = "Satterthwaite"
  )
)

study_duration_additive_critical <- stats::qt(
  0.975,
  df = study_duration_additive_CR2$df_Satt
)

study_duration_additive_lower <-
  study_duration_additive_CR2$beta -
  study_duration_additive_critical * study_duration_additive_CR2$SE

study_duration_additive_upper <-
  study_duration_additive_CR2$beta +
  study_duration_additive_critical * study_duration_additive_CR2$SE

study_duration_additive_W <- solve(study_duration_additive_V)
study_duration_additive_X <- as.matrix(study_duration_additive_model$X)
study_duration_additive_P <-
  study_duration_additive_W -
  study_duration_additive_W %*% study_duration_additive_X %*%
  solve(crossprod(study_duration_additive_X, study_duration_additive_W %*% study_duration_additive_X)) %*%
  t(study_duration_additive_X) %*% study_duration_additive_W

study_duration_additive_typical_sampling_variance <-
  (nrow(study_duration_additive_X) - ncol(study_duration_additive_X)) /
  sum(diag(study_duration_additive_P))

study_duration_additive_total_variance <-
  sum(study_duration_additive_model$sigma2) +
  study_duration_additive_typical_sampling_variance

study_duration_additive_study_I2 <-
  100 * study_duration_additive_model$sigma2[1] /
  study_duration_additive_total_variance

study_duration_additive_observation_I2 <-
  100 * study_duration_additive_model$sigma2[2] /
  study_duration_additive_total_variance

study_duration_additive_total_I2 <-
  100 * sum(study_duration_additive_model$sigma2) /
  study_duration_additive_total_variance

study_duration_additive_results <- data.frame(
  Moderator = "Study duration",
  NullModel = "Additive",
  Term = c(
    "Intercept at mean",
    "Slope per SD"
  ),
  EffectType = c(
    "Intercept at mean",
    "Slope per SD"
  ),
  ReferenceLevel = NA_character_,
  N = study_duration_N,
  n = study_duration_n,
  Level_N = study_duration_N,
  Level_n = study_duration_n,
  Estimate = as.numeric(study_duration_additive_CR2$beta),
  CR2_SE = as.numeric(study_duration_additive_CR2$SE),
  CR2_df = as.numeric(study_duration_additive_CR2$df_Satt),
  CI_lower = as.numeric(study_duration_additive_lower),
  CI_upper = as.numeric(study_duration_additive_upper),
  P_value = as.numeric(study_duration_additive_CR2$p_Satt),
  Moderator_QM = as.numeric(study_duration_additive_model$QM),
  Moderator_QM_P = as.numeric(study_duration_additive_model$QMp),
  Tau2_Study = study_duration_additive_model$sigma2[1],
  Tau2_Observation = study_duration_additive_model$sigma2[2],
  Study_I2 = study_duration_additive_study_I2,
  Observation_I2 = study_duration_additive_observation_I2,
  Total_I2 = study_duration_additive_total_I2,
  QE = as.numeric(study_duration_additive_model$QE),
  QE_P_value = as.numeric(study_duration_additive_model$QEp),
  Interpretation = c(
    ifelse(
      study_duration_additive_lower[1] > 0,
      "Positive non-additive",
      ifelse(
        study_duration_additive_upper[1] < 0,
        "Negative non-additive",
        "Additive"
      )
    ),
    ifelse(
      study_duration_additive_lower[2] > 0,
      "Positive association",
      ifelse(
        study_duration_additive_upper[2] < 0,
        "Negative association",
        "No detected association"
      )
    )
  ),
  Status = "Successfully fitted",
  stringsAsFactors = FALSE
)


# Multiplicative null model
study_duration_multiplicative_V <- clubSandwich::impute_covariance_matrix(
  vi = study_duration_data$VdMultiplicative,
  cluster = study_duration_data$StudyID,
  r = sampling_correlation,
  smooth_vi = FALSE,
  check_PD = TRUE,
  return_list = FALSE
)

study_duration_multiplicative_model <- metafor::rma.mv(
  yi = HedgesdMultiplicative,
  V = study_duration_multiplicative_V,
  mods = ~ ModeratorValue,
  random = list(
    ~ 1 | StudyID,
    ~ 1 | ObservationWithinStudy
  ),
  data = study_duration_data,
  method = "REML",
  test = "t",
  sparse = TRUE
)

study_duration_multiplicative_CR2 <- as.data.frame(
  clubSandwich::coef_test(
    study_duration_multiplicative_model,
    vcov = "CR2",
    cluster = study_duration_data$StudyID,
    test = "Satterthwaite"
  )
)

study_duration_multiplicative_critical <- stats::qt(
  0.975,
  df = study_duration_multiplicative_CR2$df_Satt
)

study_duration_multiplicative_lower <-
  study_duration_multiplicative_CR2$beta -
  study_duration_multiplicative_critical * study_duration_multiplicative_CR2$SE

study_duration_multiplicative_upper <-
  study_duration_multiplicative_CR2$beta +
  study_duration_multiplicative_critical * study_duration_multiplicative_CR2$SE

study_duration_multiplicative_W <- solve(study_duration_multiplicative_V)
study_duration_multiplicative_X <- as.matrix(study_duration_multiplicative_model$X)
study_duration_multiplicative_P <-
  study_duration_multiplicative_W -
  study_duration_multiplicative_W %*% study_duration_multiplicative_X %*%
  solve(crossprod(study_duration_multiplicative_X, study_duration_multiplicative_W %*% study_duration_multiplicative_X)) %*%
  t(study_duration_multiplicative_X) %*% study_duration_multiplicative_W

study_duration_multiplicative_typical_sampling_variance <-
  (nrow(study_duration_multiplicative_X) - ncol(study_duration_multiplicative_X)) /
  sum(diag(study_duration_multiplicative_P))

study_duration_multiplicative_total_variance <-
  sum(study_duration_multiplicative_model$sigma2) +
  study_duration_multiplicative_typical_sampling_variance

study_duration_multiplicative_study_I2 <-
  100 * study_duration_multiplicative_model$sigma2[1] /
  study_duration_multiplicative_total_variance

study_duration_multiplicative_observation_I2 <-
  100 * study_duration_multiplicative_model$sigma2[2] /
  study_duration_multiplicative_total_variance

study_duration_multiplicative_total_I2 <-
  100 * sum(study_duration_multiplicative_model$sigma2) /
  study_duration_multiplicative_total_variance

study_duration_multiplicative_results <- data.frame(
  Moderator = "Study duration",
  NullModel = "Multiplicative",
  Term = c(
    "Intercept at mean",
    "Slope per SD"
  ),
  EffectType = c(
    "Intercept at mean",
    "Slope per SD"
  ),
  ReferenceLevel = NA_character_,
  N = study_duration_N,
  n = study_duration_n,
  Level_N = study_duration_N,
  Level_n = study_duration_n,
  Estimate = as.numeric(study_duration_multiplicative_CR2$beta),
  CR2_SE = as.numeric(study_duration_multiplicative_CR2$SE),
  CR2_df = as.numeric(study_duration_multiplicative_CR2$df_Satt),
  CI_lower = as.numeric(study_duration_multiplicative_lower),
  CI_upper = as.numeric(study_duration_multiplicative_upper),
  P_value = as.numeric(study_duration_multiplicative_CR2$p_Satt),
  Moderator_QM = as.numeric(study_duration_multiplicative_model$QM),
  Moderator_QM_P = as.numeric(study_duration_multiplicative_model$QMp),
  Tau2_Study = study_duration_multiplicative_model$sigma2[1],
  Tau2_Observation = study_duration_multiplicative_model$sigma2[2],
  Study_I2 = study_duration_multiplicative_study_I2,
  Observation_I2 = study_duration_multiplicative_observation_I2,
  Total_I2 = study_duration_multiplicative_total_I2,
  QE = as.numeric(study_duration_multiplicative_model$QE),
  QE_P_value = as.numeric(study_duration_multiplicative_model$QEp),
  Interpretation = c(
    ifelse(
      study_duration_multiplicative_lower[1] > 0,
      "Positive non-additive",
      ifelse(
        study_duration_multiplicative_upper[1] < 0,
        "Negative non-additive",
        "Additive"
      )
    ),
    ifelse(
      study_duration_multiplicative_lower[2] > 0,
      "Positive association",
      ifelse(
        study_duration_multiplicative_upper[2] < 0,
        "Negative association",
        "No detected association"
      )
    )
  ),
  Status = "Successfully fitted",
  stringsAsFactors = FALSE
)


# Dominance null model
study_duration_dominance_V <- clubSandwich::impute_covariance_matrix(
  vi = study_duration_data$VdDominance,
  cluster = study_duration_data$StudyID,
  r = sampling_correlation,
  smooth_vi = FALSE,
  check_PD = TRUE,
  return_list = FALSE
)

study_duration_dominance_model <- metafor::rma.mv(
  yi = HedgesdDominance,
  V = study_duration_dominance_V,
  mods = ~ ModeratorValue,
  random = list(
    ~ 1 | StudyID,
    ~ 1 | ObservationWithinStudy
  ),
  data = study_duration_data,
  method = "REML",
  test = "t",
  sparse = TRUE
)

study_duration_dominance_CR2 <- as.data.frame(
  clubSandwich::coef_test(
    study_duration_dominance_model,
    vcov = "CR2",
    cluster = study_duration_data$StudyID,
    test = "Satterthwaite"
  )
)

study_duration_dominance_critical <- stats::qt(
  0.975,
  df = study_duration_dominance_CR2$df_Satt
)

study_duration_dominance_lower <-
  study_duration_dominance_CR2$beta -
  study_duration_dominance_critical * study_duration_dominance_CR2$SE

study_duration_dominance_upper <-
  study_duration_dominance_CR2$beta +
  study_duration_dominance_critical * study_duration_dominance_CR2$SE

study_duration_dominance_W <- solve(study_duration_dominance_V)
study_duration_dominance_X <- as.matrix(study_duration_dominance_model$X)
study_duration_dominance_P <-
  study_duration_dominance_W -
  study_duration_dominance_W %*% study_duration_dominance_X %*%
  solve(crossprod(study_duration_dominance_X, study_duration_dominance_W %*% study_duration_dominance_X)) %*%
  t(study_duration_dominance_X) %*% study_duration_dominance_W

study_duration_dominance_typical_sampling_variance <-
  (nrow(study_duration_dominance_X) - ncol(study_duration_dominance_X)) /
  sum(diag(study_duration_dominance_P))

study_duration_dominance_total_variance <-
  sum(study_duration_dominance_model$sigma2) +
  study_duration_dominance_typical_sampling_variance

study_duration_dominance_study_I2 <-
  100 * study_duration_dominance_model$sigma2[1] /
  study_duration_dominance_total_variance

study_duration_dominance_observation_I2 <-
  100 * study_duration_dominance_model$sigma2[2] /
  study_duration_dominance_total_variance

study_duration_dominance_total_I2 <-
  100 * sum(study_duration_dominance_model$sigma2) /
  study_duration_dominance_total_variance

study_duration_dominance_results <- data.frame(
  Moderator = "Study duration",
  NullModel = "Dominance",
  Term = c(
    "Intercept at mean",
    "Slope per SD"
  ),
  EffectType = c(
    "Intercept at mean",
    "Slope per SD"
  ),
  ReferenceLevel = NA_character_,
  N = study_duration_N,
  n = study_duration_n,
  Level_N = study_duration_N,
  Level_n = study_duration_n,
  Estimate = as.numeric(study_duration_dominance_CR2$beta),
  CR2_SE = as.numeric(study_duration_dominance_CR2$SE),
  CR2_df = as.numeric(study_duration_dominance_CR2$df_Satt),
  CI_lower = as.numeric(study_duration_dominance_lower),
  CI_upper = as.numeric(study_duration_dominance_upper),
  P_value = as.numeric(study_duration_dominance_CR2$p_Satt),
  Moderator_QM = as.numeric(study_duration_dominance_model$QM),
  Moderator_QM_P = as.numeric(study_duration_dominance_model$QMp),
  Tau2_Study = study_duration_dominance_model$sigma2[1],
  Tau2_Observation = study_duration_dominance_model$sigma2[2],
  Study_I2 = study_duration_dominance_study_I2,
  Observation_I2 = study_duration_dominance_observation_I2,
  Total_I2 = study_duration_dominance_total_I2,
  QE = as.numeric(study_duration_dominance_model$QE),
  QE_P_value = as.numeric(study_duration_dominance_model$QEp),
  Interpretation = c(
    ifelse(
      study_duration_dominance_lower[1] > 0,
      "Positive non-additive",
      ifelse(
        study_duration_dominance_upper[1] < 0,
        "Negative non-additive",
        "Additive"
      )
    ),
    ifelse(
      study_duration_dominance_lower[2] > 0,
      "Positive association",
      ifelse(
        study_duration_dominance_upper[2] < 0,
        "Negative association",
        "No detected association"
      )
    )
  ),
  Status = "Successfully fitted",
  stringsAsFactors = FALSE
)

study_duration_results <- dplyr::bind_rows(
  study_duration_additive_results,
  study_duration_multiplicative_results,
  study_duration_dominance_results
)




# 05.2.6. Temperature tolerance -----------------------------------------------

temperature_tolerance_data <- moderator_data[
  !is.na(moderator_data$FactorPair) &
  grepl("temperature", as.character(moderator_data$FactorPair), ignore.case = TRUE) &
  !is.na(moderator_data$TemperatureTolerance) &
  trimws(as.character(moderator_data$TemperatureTolerance)) != "" &
    is.finite(moderator_data$HedgesdAdditive) &
    is.finite(moderator_data$VdAdditive) &
    moderator_data$VdAdditive > 0 &
    is.finite(moderator_data$HedgesdMultiplicative) &
    is.finite(moderator_data$VdMultiplicative) &
    moderator_data$VdMultiplicative > 0 &
    is.finite(moderator_data$HedgesdDominance) &
    is.finite(moderator_data$VdDominance) &
    moderator_data$VdDominance > 0 &
    !is.na(moderator_data$StudyID) &
    !is.na(moderator_data$ObservationID),
  ,
  drop = FALSE
]

temperature_tolerance_data$ModeratorValue <- droplevels(
  factor(temperature_tolerance_data$TemperatureTolerance)
)

temperature_tolerance_data$StudyID <- factor(
  temperature_tolerance_data$StudyID
)

temperature_tolerance_data$ObservationWithinStudy <- factor(
  paste(
    temperature_tolerance_data$StudyID,
    temperature_tolerance_data$ObservationID,
    sep = "_"
  )
)

temperature_tolerance_levels <- levels(
  temperature_tolerance_data$ModeratorValue
)

temperature_tolerance_N <- length(
  unique(temperature_tolerance_data$StudyID)
)

temperature_tolerance_n <- nrow(
  temperature_tolerance_data
)

temperature_tolerance_level_n_table <- as.data.frame(
  table(temperature_tolerance_data$ModeratorValue),
  stringsAsFactors = FALSE
)

temperature_tolerance_level_study_data <- unique(
  data.frame(
    Level = temperature_tolerance_data$ModeratorValue,
    StudyID = temperature_tolerance_data$StudyID
  )
)

temperature_tolerance_level_N_table <- aggregate(
  StudyID ~ Level,
  data = temperature_tolerance_level_study_data,
  FUN = length
)

temperature_tolerance_level_n <- as.integer(
  temperature_tolerance_level_n_table$Freq[
    match(
      temperature_tolerance_levels,
      as.character(temperature_tolerance_level_n_table$Var1)
    )
  ]
)

temperature_tolerance_level_N <- as.integer(
  temperature_tolerance_level_N_table$StudyID[
    match(
      temperature_tolerance_levels,
      as.character(temperature_tolerance_level_N_table$Level)
    )
  ]
)

# Additive null model
temperature_tolerance_additive_V <- clubSandwich::impute_covariance_matrix(
  vi = temperature_tolerance_data$VdAdditive,
  cluster = temperature_tolerance_data$StudyID,
  r = sampling_correlation,
  smooth_vi = FALSE,
  check_PD = TRUE,
  return_list = FALSE
)

temperature_tolerance_additive_level_model <- metafor::rma.mv(
  yi = HedgesdAdditive,
  V = temperature_tolerance_additive_V,
  mods = ~ 0 + ModeratorValue,
  random = list(
    ~ 1 | StudyID,
    ~ 1 | ObservationWithinStudy
  ),
  data = temperature_tolerance_data,
  method = "REML",
  test = "t",
  sparse = TRUE
)

temperature_tolerance_additive_contrast_model <- metafor::rma.mv(
  yi = HedgesdAdditive,
  V = temperature_tolerance_additive_V,
  mods = ~ ModeratorValue,
  random = list(
    ~ 1 | StudyID,
    ~ 1 | ObservationWithinStudy
  ),
  data = temperature_tolerance_data,
  method = "REML",
  test = "t",
  sparse = TRUE
)

temperature_tolerance_additive_level_CR2 <- as.data.frame(
  clubSandwich::coef_test(
    temperature_tolerance_additive_level_model,
    vcov = "CR2",
    cluster = temperature_tolerance_data$StudyID,
    test = "Satterthwaite"
  )
)

temperature_tolerance_additive_contrast_CR2 <- as.data.frame(
  clubSandwich::coef_test(
    temperature_tolerance_additive_contrast_model,
    vcov = "CR2",
    cluster = temperature_tolerance_data$StudyID,
    test = "Satterthwaite"
  )
)

temperature_tolerance_additive_level_critical <- stats::qt(
  0.975,
  df = temperature_tolerance_additive_level_CR2$df_Satt
)

temperature_tolerance_additive_level_lower <-
  temperature_tolerance_additive_level_CR2$beta -
  temperature_tolerance_additive_level_critical * temperature_tolerance_additive_level_CR2$SE

temperature_tolerance_additive_level_upper <-
  temperature_tolerance_additive_level_CR2$beta +
  temperature_tolerance_additive_level_critical * temperature_tolerance_additive_level_CR2$SE

temperature_tolerance_additive_contrast_CR2 <- temperature_tolerance_additive_contrast_CR2[
  -1,
  ,
  drop = FALSE
]

temperature_tolerance_additive_contrast_critical <- stats::qt(
  0.975,
  df = temperature_tolerance_additive_contrast_CR2$df_Satt
)

temperature_tolerance_additive_contrast_lower <-
  temperature_tolerance_additive_contrast_CR2$beta -
  temperature_tolerance_additive_contrast_critical * temperature_tolerance_additive_contrast_CR2$SE

temperature_tolerance_additive_contrast_upper <-
  temperature_tolerance_additive_contrast_CR2$beta +
  temperature_tolerance_additive_contrast_critical * temperature_tolerance_additive_contrast_CR2$SE

temperature_tolerance_additive_W <- solve(temperature_tolerance_additive_V)
temperature_tolerance_additive_X <- as.matrix(temperature_tolerance_additive_contrast_model$X)
temperature_tolerance_additive_P <-
  temperature_tolerance_additive_W -
  temperature_tolerance_additive_W %*% temperature_tolerance_additive_X %*%
  solve(crossprod(temperature_tolerance_additive_X, temperature_tolerance_additive_W %*% temperature_tolerance_additive_X)) %*%
  t(temperature_tolerance_additive_X) %*% temperature_tolerance_additive_W

temperature_tolerance_additive_typical_sampling_variance <-
  (nrow(temperature_tolerance_additive_X) - ncol(temperature_tolerance_additive_X)) /
  sum(diag(temperature_tolerance_additive_P))

temperature_tolerance_additive_total_variance <-
  sum(temperature_tolerance_additive_contrast_model$sigma2) +
  temperature_tolerance_additive_typical_sampling_variance

temperature_tolerance_additive_study_I2 <-
  100 * temperature_tolerance_additive_contrast_model$sigma2[1] /
  temperature_tolerance_additive_total_variance

temperature_tolerance_additive_observation_I2 <-
  100 * temperature_tolerance_additive_contrast_model$sigma2[2] /
  temperature_tolerance_additive_total_variance

temperature_tolerance_additive_total_I2 <-
  100 * sum(temperature_tolerance_additive_contrast_model$sigma2) /
  temperature_tolerance_additive_total_variance

temperature_tolerance_additive_level_results <- data.frame(
  Moderator = "Temperature tolerance",
  NullModel = "Additive",
  Term = temperature_tolerance_levels,
  EffectType = "Level mean",
  ReferenceLevel = temperature_tolerance_levels[1],
  N = temperature_tolerance_N,
  n = temperature_tolerance_n,
  Level_N = temperature_tolerance_level_N,
  Level_n = temperature_tolerance_level_n,
  Estimate = as.numeric(temperature_tolerance_additive_level_CR2$beta),
  CR2_SE = as.numeric(temperature_tolerance_additive_level_CR2$SE),
  CR2_df = as.numeric(temperature_tolerance_additive_level_CR2$df_Satt),
  CI_lower = as.numeric(temperature_tolerance_additive_level_lower),
  CI_upper = as.numeric(temperature_tolerance_additive_level_upper),
  P_value = as.numeric(temperature_tolerance_additive_level_CR2$p_Satt),
  Moderator_QM = as.numeric(temperature_tolerance_additive_contrast_model$QM),
  Moderator_QM_P = as.numeric(temperature_tolerance_additive_contrast_model$QMp),
  Tau2_Study = temperature_tolerance_additive_contrast_model$sigma2[1],
  Tau2_Observation = temperature_tolerance_additive_contrast_model$sigma2[2],
  Study_I2 = temperature_tolerance_additive_study_I2,
  Observation_I2 = temperature_tolerance_additive_observation_I2,
  Total_I2 = temperature_tolerance_additive_total_I2,
  QE = as.numeric(temperature_tolerance_additive_contrast_model$QE),
  QE_P_value = as.numeric(temperature_tolerance_additive_contrast_model$QEp),
  Interpretation = ifelse(
    temperature_tolerance_additive_level_lower > 0,
    "Positive non-additive",
    ifelse(
      temperature_tolerance_additive_level_upper < 0,
      "Negative non-additive",
      "Additive"
    )
  ),
  Status = "Successfully fitted",
  stringsAsFactors = FALSE
)

temperature_tolerance_additive_contrast_results <- data.frame(
  Moderator = "Temperature tolerance",
  NullModel = "Additive",
  Term = paste0(
    temperature_tolerance_levels[-1],
    " vs ",
    temperature_tolerance_levels[1]
  ),
  EffectType = "Contrast with reference",
  ReferenceLevel = temperature_tolerance_levels[1],
  N = temperature_tolerance_N,
  n = temperature_tolerance_n,
  Level_N = temperature_tolerance_N,
  Level_n = temperature_tolerance_n,
  Estimate = as.numeric(temperature_tolerance_additive_contrast_CR2$beta),
  CR2_SE = as.numeric(temperature_tolerance_additive_contrast_CR2$SE),
  CR2_df = as.numeric(temperature_tolerance_additive_contrast_CR2$df_Satt),
  CI_lower = as.numeric(temperature_tolerance_additive_contrast_lower),
  CI_upper = as.numeric(temperature_tolerance_additive_contrast_upper),
  P_value = as.numeric(temperature_tolerance_additive_contrast_CR2$p_Satt),
  Moderator_QM = as.numeric(temperature_tolerance_additive_contrast_model$QM),
  Moderator_QM_P = as.numeric(temperature_tolerance_additive_contrast_model$QMp),
  Tau2_Study = temperature_tolerance_additive_contrast_model$sigma2[1],
  Tau2_Observation = temperature_tolerance_additive_contrast_model$sigma2[2],
  Study_I2 = temperature_tolerance_additive_study_I2,
  Observation_I2 = temperature_tolerance_additive_observation_I2,
  Total_I2 = temperature_tolerance_additive_total_I2,
  QE = as.numeric(temperature_tolerance_additive_contrast_model$QE),
  QE_P_value = as.numeric(temperature_tolerance_additive_contrast_model$QEp),
  Interpretation = ifelse(
    temperature_tolerance_additive_contrast_lower > 0,
    "Higher than reference",
    ifelse(
      temperature_tolerance_additive_contrast_upper < 0,
      "Lower than reference",
      "No detected difference"
    )
  ),
  Status = "Successfully fitted",
  stringsAsFactors = FALSE
)

temperature_tolerance_additive_results <- dplyr::bind_rows(
  temperature_tolerance_additive_level_results,
  temperature_tolerance_additive_contrast_results
)


# Multiplicative null model
temperature_tolerance_multiplicative_V <- clubSandwich::impute_covariance_matrix(
  vi = temperature_tolerance_data$VdMultiplicative,
  cluster = temperature_tolerance_data$StudyID,
  r = sampling_correlation,
  smooth_vi = FALSE,
  check_PD = TRUE,
  return_list = FALSE
)

temperature_tolerance_multiplicative_level_model <- metafor::rma.mv(
  yi = HedgesdMultiplicative,
  V = temperature_tolerance_multiplicative_V,
  mods = ~ 0 + ModeratorValue,
  random = list(
    ~ 1 | StudyID,
    ~ 1 | ObservationWithinStudy
  ),
  data = temperature_tolerance_data,
  method = "REML",
  test = "t",
  sparse = TRUE
)

temperature_tolerance_multiplicative_contrast_model <- metafor::rma.mv(
  yi = HedgesdMultiplicative,
  V = temperature_tolerance_multiplicative_V,
  mods = ~ ModeratorValue,
  random = list(
    ~ 1 | StudyID,
    ~ 1 | ObservationWithinStudy
  ),
  data = temperature_tolerance_data,
  method = "REML",
  test = "t",
  sparse = TRUE
)

temperature_tolerance_multiplicative_level_CR2 <- as.data.frame(
  clubSandwich::coef_test(
    temperature_tolerance_multiplicative_level_model,
    vcov = "CR2",
    cluster = temperature_tolerance_data$StudyID,
    test = "Satterthwaite"
  )
)

temperature_tolerance_multiplicative_contrast_CR2 <- as.data.frame(
  clubSandwich::coef_test(
    temperature_tolerance_multiplicative_contrast_model,
    vcov = "CR2",
    cluster = temperature_tolerance_data$StudyID,
    test = "Satterthwaite"
  )
)

temperature_tolerance_multiplicative_level_critical <- stats::qt(
  0.975,
  df = temperature_tolerance_multiplicative_level_CR2$df_Satt
)

temperature_tolerance_multiplicative_level_lower <-
  temperature_tolerance_multiplicative_level_CR2$beta -
  temperature_tolerance_multiplicative_level_critical * temperature_tolerance_multiplicative_level_CR2$SE

temperature_tolerance_multiplicative_level_upper <-
  temperature_tolerance_multiplicative_level_CR2$beta +
  temperature_tolerance_multiplicative_level_critical * temperature_tolerance_multiplicative_level_CR2$SE

temperature_tolerance_multiplicative_contrast_CR2 <- temperature_tolerance_multiplicative_contrast_CR2[
  -1,
  ,
  drop = FALSE
]

temperature_tolerance_multiplicative_contrast_critical <- stats::qt(
  0.975,
  df = temperature_tolerance_multiplicative_contrast_CR2$df_Satt
)

temperature_tolerance_multiplicative_contrast_lower <-
  temperature_tolerance_multiplicative_contrast_CR2$beta -
  temperature_tolerance_multiplicative_contrast_critical * temperature_tolerance_multiplicative_contrast_CR2$SE

temperature_tolerance_multiplicative_contrast_upper <-
  temperature_tolerance_multiplicative_contrast_CR2$beta +
  temperature_tolerance_multiplicative_contrast_critical * temperature_tolerance_multiplicative_contrast_CR2$SE

temperature_tolerance_multiplicative_W <- solve(temperature_tolerance_multiplicative_V)
temperature_tolerance_multiplicative_X <- as.matrix(temperature_tolerance_multiplicative_contrast_model$X)
temperature_tolerance_multiplicative_P <-
  temperature_tolerance_multiplicative_W -
  temperature_tolerance_multiplicative_W %*% temperature_tolerance_multiplicative_X %*%
  solve(crossprod(temperature_tolerance_multiplicative_X, temperature_tolerance_multiplicative_W %*% temperature_tolerance_multiplicative_X)) %*%
  t(temperature_tolerance_multiplicative_X) %*% temperature_tolerance_multiplicative_W

temperature_tolerance_multiplicative_typical_sampling_variance <-
  (nrow(temperature_tolerance_multiplicative_X) - ncol(temperature_tolerance_multiplicative_X)) /
  sum(diag(temperature_tolerance_multiplicative_P))

temperature_tolerance_multiplicative_total_variance <-
  sum(temperature_tolerance_multiplicative_contrast_model$sigma2) +
  temperature_tolerance_multiplicative_typical_sampling_variance

temperature_tolerance_multiplicative_study_I2 <-
  100 * temperature_tolerance_multiplicative_contrast_model$sigma2[1] /
  temperature_tolerance_multiplicative_total_variance

temperature_tolerance_multiplicative_observation_I2 <-
  100 * temperature_tolerance_multiplicative_contrast_model$sigma2[2] /
  temperature_tolerance_multiplicative_total_variance

temperature_tolerance_multiplicative_total_I2 <-
  100 * sum(temperature_tolerance_multiplicative_contrast_model$sigma2) /
  temperature_tolerance_multiplicative_total_variance

temperature_tolerance_multiplicative_level_results <- data.frame(
  Moderator = "Temperature tolerance",
  NullModel = "Multiplicative",
  Term = temperature_tolerance_levels,
  EffectType = "Level mean",
  ReferenceLevel = temperature_tolerance_levels[1],
  N = temperature_tolerance_N,
  n = temperature_tolerance_n,
  Level_N = temperature_tolerance_level_N,
  Level_n = temperature_tolerance_level_n,
  Estimate = as.numeric(temperature_tolerance_multiplicative_level_CR2$beta),
  CR2_SE = as.numeric(temperature_tolerance_multiplicative_level_CR2$SE),
  CR2_df = as.numeric(temperature_tolerance_multiplicative_level_CR2$df_Satt),
  CI_lower = as.numeric(temperature_tolerance_multiplicative_level_lower),
  CI_upper = as.numeric(temperature_tolerance_multiplicative_level_upper),
  P_value = as.numeric(temperature_tolerance_multiplicative_level_CR2$p_Satt),
  Moderator_QM = as.numeric(temperature_tolerance_multiplicative_contrast_model$QM),
  Moderator_QM_P = as.numeric(temperature_tolerance_multiplicative_contrast_model$QMp),
  Tau2_Study = temperature_tolerance_multiplicative_contrast_model$sigma2[1],
  Tau2_Observation = temperature_tolerance_multiplicative_contrast_model$sigma2[2],
  Study_I2 = temperature_tolerance_multiplicative_study_I2,
  Observation_I2 = temperature_tolerance_multiplicative_observation_I2,
  Total_I2 = temperature_tolerance_multiplicative_total_I2,
  QE = as.numeric(temperature_tolerance_multiplicative_contrast_model$QE),
  QE_P_value = as.numeric(temperature_tolerance_multiplicative_contrast_model$QEp),
  Interpretation = ifelse(
    temperature_tolerance_multiplicative_level_lower > 0,
    "Positive non-additive",
    ifelse(
      temperature_tolerance_multiplicative_level_upper < 0,
      "Negative non-additive",
      "Additive"
    )
  ),
  Status = "Successfully fitted",
  stringsAsFactors = FALSE
)

temperature_tolerance_multiplicative_contrast_results <- data.frame(
  Moderator = "Temperature tolerance",
  NullModel = "Multiplicative",
  Term = paste0(
    temperature_tolerance_levels[-1],
    " vs ",
    temperature_tolerance_levels[1]
  ),
  EffectType = "Contrast with reference",
  ReferenceLevel = temperature_tolerance_levels[1],
  N = temperature_tolerance_N,
  n = temperature_tolerance_n,
  Level_N = temperature_tolerance_N,
  Level_n = temperature_tolerance_n,
  Estimate = as.numeric(temperature_tolerance_multiplicative_contrast_CR2$beta),
  CR2_SE = as.numeric(temperature_tolerance_multiplicative_contrast_CR2$SE),
  CR2_df = as.numeric(temperature_tolerance_multiplicative_contrast_CR2$df_Satt),
  CI_lower = as.numeric(temperature_tolerance_multiplicative_contrast_lower),
  CI_upper = as.numeric(temperature_tolerance_multiplicative_contrast_upper),
  P_value = as.numeric(temperature_tolerance_multiplicative_contrast_CR2$p_Satt),
  Moderator_QM = as.numeric(temperature_tolerance_multiplicative_contrast_model$QM),
  Moderator_QM_P = as.numeric(temperature_tolerance_multiplicative_contrast_model$QMp),
  Tau2_Study = temperature_tolerance_multiplicative_contrast_model$sigma2[1],
  Tau2_Observation = temperature_tolerance_multiplicative_contrast_model$sigma2[2],
  Study_I2 = temperature_tolerance_multiplicative_study_I2,
  Observation_I2 = temperature_tolerance_multiplicative_observation_I2,
  Total_I2 = temperature_tolerance_multiplicative_total_I2,
  QE = as.numeric(temperature_tolerance_multiplicative_contrast_model$QE),
  QE_P_value = as.numeric(temperature_tolerance_multiplicative_contrast_model$QEp),
  Interpretation = ifelse(
    temperature_tolerance_multiplicative_contrast_lower > 0,
    "Higher than reference",
    ifelse(
      temperature_tolerance_multiplicative_contrast_upper < 0,
      "Lower than reference",
      "No detected difference"
    )
  ),
  Status = "Successfully fitted",
  stringsAsFactors = FALSE
)

temperature_tolerance_multiplicative_results <- dplyr::bind_rows(
  temperature_tolerance_multiplicative_level_results,
  temperature_tolerance_multiplicative_contrast_results
)


# Dominance null model
temperature_tolerance_dominance_V <- clubSandwich::impute_covariance_matrix(
  vi = temperature_tolerance_data$VdDominance,
  cluster = temperature_tolerance_data$StudyID,
  r = sampling_correlation,
  smooth_vi = FALSE,
  check_PD = TRUE,
  return_list = FALSE
)

temperature_tolerance_dominance_level_model <- metafor::rma.mv(
  yi = HedgesdDominance,
  V = temperature_tolerance_dominance_V,
  mods = ~ 0 + ModeratorValue,
  random = list(
    ~ 1 | StudyID,
    ~ 1 | ObservationWithinStudy
  ),
  data = temperature_tolerance_data,
  method = "REML",
  test = "t",
  sparse = TRUE
)

temperature_tolerance_dominance_contrast_model <- metafor::rma.mv(
  yi = HedgesdDominance,
  V = temperature_tolerance_dominance_V,
  mods = ~ ModeratorValue,
  random = list(
    ~ 1 | StudyID,
    ~ 1 | ObservationWithinStudy
  ),
  data = temperature_tolerance_data,
  method = "REML",
  test = "t",
  sparse = TRUE
)

temperature_tolerance_dominance_level_CR2 <- as.data.frame(
  clubSandwich::coef_test(
    temperature_tolerance_dominance_level_model,
    vcov = "CR2",
    cluster = temperature_tolerance_data$StudyID,
    test = "Satterthwaite"
  )
)

temperature_tolerance_dominance_contrast_CR2 <- as.data.frame(
  clubSandwich::coef_test(
    temperature_tolerance_dominance_contrast_model,
    vcov = "CR2",
    cluster = temperature_tolerance_data$StudyID,
    test = "Satterthwaite"
  )
)

temperature_tolerance_dominance_level_critical <- stats::qt(
  0.975,
  df = temperature_tolerance_dominance_level_CR2$df_Satt
)

temperature_tolerance_dominance_level_lower <-
  temperature_tolerance_dominance_level_CR2$beta -
  temperature_tolerance_dominance_level_critical * temperature_tolerance_dominance_level_CR2$SE

temperature_tolerance_dominance_level_upper <-
  temperature_tolerance_dominance_level_CR2$beta +
  temperature_tolerance_dominance_level_critical * temperature_tolerance_dominance_level_CR2$SE

temperature_tolerance_dominance_contrast_CR2 <- temperature_tolerance_dominance_contrast_CR2[
  -1,
  ,
  drop = FALSE
]

temperature_tolerance_dominance_contrast_critical <- stats::qt(
  0.975,
  df = temperature_tolerance_dominance_contrast_CR2$df_Satt
)

temperature_tolerance_dominance_contrast_lower <-
  temperature_tolerance_dominance_contrast_CR2$beta -
  temperature_tolerance_dominance_contrast_critical * temperature_tolerance_dominance_contrast_CR2$SE

temperature_tolerance_dominance_contrast_upper <-
  temperature_tolerance_dominance_contrast_CR2$beta +
  temperature_tolerance_dominance_contrast_critical * temperature_tolerance_dominance_contrast_CR2$SE

temperature_tolerance_dominance_W <- solve(temperature_tolerance_dominance_V)
temperature_tolerance_dominance_X <- as.matrix(temperature_tolerance_dominance_contrast_model$X)
temperature_tolerance_dominance_P <-
  temperature_tolerance_dominance_W -
  temperature_tolerance_dominance_W %*% temperature_tolerance_dominance_X %*%
  solve(crossprod(temperature_tolerance_dominance_X, temperature_tolerance_dominance_W %*% temperature_tolerance_dominance_X)) %*%
  t(temperature_tolerance_dominance_X) %*% temperature_tolerance_dominance_W

temperature_tolerance_dominance_typical_sampling_variance <-
  (nrow(temperature_tolerance_dominance_X) - ncol(temperature_tolerance_dominance_X)) /
  sum(diag(temperature_tolerance_dominance_P))

temperature_tolerance_dominance_total_variance <-
  sum(temperature_tolerance_dominance_contrast_model$sigma2) +
  temperature_tolerance_dominance_typical_sampling_variance

temperature_tolerance_dominance_study_I2 <-
  100 * temperature_tolerance_dominance_contrast_model$sigma2[1] /
  temperature_tolerance_dominance_total_variance

temperature_tolerance_dominance_observation_I2 <-
  100 * temperature_tolerance_dominance_contrast_model$sigma2[2] /
  temperature_tolerance_dominance_total_variance

temperature_tolerance_dominance_total_I2 <-
  100 * sum(temperature_tolerance_dominance_contrast_model$sigma2) /
  temperature_tolerance_dominance_total_variance

temperature_tolerance_dominance_level_results <- data.frame(
  Moderator = "Temperature tolerance",
  NullModel = "Dominance",
  Term = temperature_tolerance_levels,
  EffectType = "Level mean",
  ReferenceLevel = temperature_tolerance_levels[1],
  N = temperature_tolerance_N,
  n = temperature_tolerance_n,
  Level_N = temperature_tolerance_level_N,
  Level_n = temperature_tolerance_level_n,
  Estimate = as.numeric(temperature_tolerance_dominance_level_CR2$beta),
  CR2_SE = as.numeric(temperature_tolerance_dominance_level_CR2$SE),
  CR2_df = as.numeric(temperature_tolerance_dominance_level_CR2$df_Satt),
  CI_lower = as.numeric(temperature_tolerance_dominance_level_lower),
  CI_upper = as.numeric(temperature_tolerance_dominance_level_upper),
  P_value = as.numeric(temperature_tolerance_dominance_level_CR2$p_Satt),
  Moderator_QM = as.numeric(temperature_tolerance_dominance_contrast_model$QM),
  Moderator_QM_P = as.numeric(temperature_tolerance_dominance_contrast_model$QMp),
  Tau2_Study = temperature_tolerance_dominance_contrast_model$sigma2[1],
  Tau2_Observation = temperature_tolerance_dominance_contrast_model$sigma2[2],
  Study_I2 = temperature_tolerance_dominance_study_I2,
  Observation_I2 = temperature_tolerance_dominance_observation_I2,
  Total_I2 = temperature_tolerance_dominance_total_I2,
  QE = as.numeric(temperature_tolerance_dominance_contrast_model$QE),
  QE_P_value = as.numeric(temperature_tolerance_dominance_contrast_model$QEp),
  Interpretation = ifelse(
    temperature_tolerance_dominance_level_lower > 0,
    "Positive non-additive",
    ifelse(
      temperature_tolerance_dominance_level_upper < 0,
      "Negative non-additive",
      "Additive"
    )
  ),
  Status = "Successfully fitted",
  stringsAsFactors = FALSE
)

temperature_tolerance_dominance_contrast_results <- data.frame(
  Moderator = "Temperature tolerance",
  NullModel = "Dominance",
  Term = paste0(
    temperature_tolerance_levels[-1],
    " vs ",
    temperature_tolerance_levels[1]
  ),
  EffectType = "Contrast with reference",
  ReferenceLevel = temperature_tolerance_levels[1],
  N = temperature_tolerance_N,
  n = temperature_tolerance_n,
  Level_N = temperature_tolerance_N,
  Level_n = temperature_tolerance_n,
  Estimate = as.numeric(temperature_tolerance_dominance_contrast_CR2$beta),
  CR2_SE = as.numeric(temperature_tolerance_dominance_contrast_CR2$SE),
  CR2_df = as.numeric(temperature_tolerance_dominance_contrast_CR2$df_Satt),
  CI_lower = as.numeric(temperature_tolerance_dominance_contrast_lower),
  CI_upper = as.numeric(temperature_tolerance_dominance_contrast_upper),
  P_value = as.numeric(temperature_tolerance_dominance_contrast_CR2$p_Satt),
  Moderator_QM = as.numeric(temperature_tolerance_dominance_contrast_model$QM),
  Moderator_QM_P = as.numeric(temperature_tolerance_dominance_contrast_model$QMp),
  Tau2_Study = temperature_tolerance_dominance_contrast_model$sigma2[1],
  Tau2_Observation = temperature_tolerance_dominance_contrast_model$sigma2[2],
  Study_I2 = temperature_tolerance_dominance_study_I2,
  Observation_I2 = temperature_tolerance_dominance_observation_I2,
  Total_I2 = temperature_tolerance_dominance_total_I2,
  QE = as.numeric(temperature_tolerance_dominance_contrast_model$QE),
  QE_P_value = as.numeric(temperature_tolerance_dominance_contrast_model$QEp),
  Interpretation = ifelse(
    temperature_tolerance_dominance_contrast_lower > 0,
    "Higher than reference",
    ifelse(
      temperature_tolerance_dominance_contrast_upper < 0,
      "Lower than reference",
      "No detected difference"
    )
  ),
  Status = "Successfully fitted",
  stringsAsFactors = FALSE
)

temperature_tolerance_dominance_results <- dplyr::bind_rows(
  temperature_tolerance_dominance_level_results,
  temperature_tolerance_dominance_contrast_results
)

temperature_tolerance_results <- dplyr::bind_rows(
  temperature_tolerance_additive_results,
  temperature_tolerance_multiplicative_results,
  temperature_tolerance_dominance_results
)




# 05.2.7. Nutrient tolerance --------------------------------------------------

nutrient_tolerance_data <- moderator_data[
    !is.na(moderator_data$NutrientTolerance) &
  trimws(as.character(moderator_data$NutrientTolerance)) != "" &
    is.finite(moderator_data$HedgesdAdditive) &
    is.finite(moderator_data$VdAdditive) &
    moderator_data$VdAdditive > 0 &
    is.finite(moderator_data$HedgesdMultiplicative) &
    is.finite(moderator_data$VdMultiplicative) &
    moderator_data$VdMultiplicative > 0 &
    is.finite(moderator_data$HedgesdDominance) &
    is.finite(moderator_data$VdDominance) &
    moderator_data$VdDominance > 0 &
    !is.na(moderator_data$StudyID) &
    !is.na(moderator_data$ObservationID),
  ,
  drop = FALSE
]

nutrient_tolerance_data$ModeratorValue <- droplevels(
  factor(nutrient_tolerance_data$NutrientTolerance)
)

nutrient_tolerance_data$StudyID <- factor(
  nutrient_tolerance_data$StudyID
)

nutrient_tolerance_data$ObservationWithinStudy <- factor(
  paste(
    nutrient_tolerance_data$StudyID,
    nutrient_tolerance_data$ObservationID,
    sep = "_"
  )
)

nutrient_tolerance_levels <- levels(
  nutrient_tolerance_data$ModeratorValue
)

nutrient_tolerance_N <- length(
  unique(nutrient_tolerance_data$StudyID)
)

nutrient_tolerance_n <- nrow(
  nutrient_tolerance_data
)

nutrient_tolerance_level_n_table <- as.data.frame(
  table(nutrient_tolerance_data$ModeratorValue),
  stringsAsFactors = FALSE
)

nutrient_tolerance_level_study_data <- unique(
  data.frame(
    Level = nutrient_tolerance_data$ModeratorValue,
    StudyID = nutrient_tolerance_data$StudyID
  )
)

nutrient_tolerance_level_N_table <- aggregate(
  StudyID ~ Level,
  data = nutrient_tolerance_level_study_data,
  FUN = length
)

nutrient_tolerance_level_n <- as.integer(
  nutrient_tolerance_level_n_table$Freq[
    match(
      nutrient_tolerance_levels,
      as.character(nutrient_tolerance_level_n_table$Var1)
    )
  ]
)

nutrient_tolerance_level_N <- as.integer(
  nutrient_tolerance_level_N_table$StudyID[
    match(
      nutrient_tolerance_levels,
      as.character(nutrient_tolerance_level_N_table$Level)
    )
  ]
)

# Additive null model
nutrient_tolerance_additive_V <- clubSandwich::impute_covariance_matrix(
  vi = nutrient_tolerance_data$VdAdditive,
  cluster = nutrient_tolerance_data$StudyID,
  r = sampling_correlation,
  smooth_vi = FALSE,
  check_PD = TRUE,
  return_list = FALSE
)

nutrient_tolerance_additive_level_model <- metafor::rma.mv(
  yi = HedgesdAdditive,
  V = nutrient_tolerance_additive_V,
  mods = ~ 0 + ModeratorValue,
  random = list(
    ~ 1 | StudyID,
    ~ 1 | ObservationWithinStudy
  ),
  data = nutrient_tolerance_data,
  method = "REML",
  test = "t",
  sparse = TRUE
)

nutrient_tolerance_additive_contrast_model <- metafor::rma.mv(
  yi = HedgesdAdditive,
  V = nutrient_tolerance_additive_V,
  mods = ~ ModeratorValue,
  random = list(
    ~ 1 | StudyID,
    ~ 1 | ObservationWithinStudy
  ),
  data = nutrient_tolerance_data,
  method = "REML",
  test = "t",
  sparse = TRUE
)

nutrient_tolerance_additive_level_CR2 <- as.data.frame(
  clubSandwich::coef_test(
    nutrient_tolerance_additive_level_model,
    vcov = "CR2",
    cluster = nutrient_tolerance_data$StudyID,
    test = "Satterthwaite"
  )
)

nutrient_tolerance_additive_contrast_CR2 <- as.data.frame(
  clubSandwich::coef_test(
    nutrient_tolerance_additive_contrast_model,
    vcov = "CR2",
    cluster = nutrient_tolerance_data$StudyID,
    test = "Satterthwaite"
  )
)

nutrient_tolerance_additive_level_critical <- stats::qt(
  0.975,
  df = nutrient_tolerance_additive_level_CR2$df_Satt
)

nutrient_tolerance_additive_level_lower <-
  nutrient_tolerance_additive_level_CR2$beta -
  nutrient_tolerance_additive_level_critical * nutrient_tolerance_additive_level_CR2$SE

nutrient_tolerance_additive_level_upper <-
  nutrient_tolerance_additive_level_CR2$beta +
  nutrient_tolerance_additive_level_critical * nutrient_tolerance_additive_level_CR2$SE

nutrient_tolerance_additive_contrast_CR2 <- nutrient_tolerance_additive_contrast_CR2[
  -1,
  ,
  drop = FALSE
]

nutrient_tolerance_additive_contrast_critical <- stats::qt(
  0.975,
  df = nutrient_tolerance_additive_contrast_CR2$df_Satt
)

nutrient_tolerance_additive_contrast_lower <-
  nutrient_tolerance_additive_contrast_CR2$beta -
  nutrient_tolerance_additive_contrast_critical * nutrient_tolerance_additive_contrast_CR2$SE

nutrient_tolerance_additive_contrast_upper <-
  nutrient_tolerance_additive_contrast_CR2$beta +
  nutrient_tolerance_additive_contrast_critical * nutrient_tolerance_additive_contrast_CR2$SE

nutrient_tolerance_additive_W <- solve(nutrient_tolerance_additive_V)
nutrient_tolerance_additive_X <- as.matrix(nutrient_tolerance_additive_contrast_model$X)
nutrient_tolerance_additive_P <-
  nutrient_tolerance_additive_W -
  nutrient_tolerance_additive_W %*% nutrient_tolerance_additive_X %*%
  solve(crossprod(nutrient_tolerance_additive_X, nutrient_tolerance_additive_W %*% nutrient_tolerance_additive_X)) %*%
  t(nutrient_tolerance_additive_X) %*% nutrient_tolerance_additive_W

nutrient_tolerance_additive_typical_sampling_variance <-
  (nrow(nutrient_tolerance_additive_X) - ncol(nutrient_tolerance_additive_X)) /
  sum(diag(nutrient_tolerance_additive_P))

nutrient_tolerance_additive_total_variance <-
  sum(nutrient_tolerance_additive_contrast_model$sigma2) +
  nutrient_tolerance_additive_typical_sampling_variance

nutrient_tolerance_additive_study_I2 <-
  100 * nutrient_tolerance_additive_contrast_model$sigma2[1] /
  nutrient_tolerance_additive_total_variance

nutrient_tolerance_additive_observation_I2 <-
  100 * nutrient_tolerance_additive_contrast_model$sigma2[2] /
  nutrient_tolerance_additive_total_variance

nutrient_tolerance_additive_total_I2 <-
  100 * sum(nutrient_tolerance_additive_contrast_model$sigma2) /
  nutrient_tolerance_additive_total_variance

nutrient_tolerance_additive_level_results <- data.frame(
  Moderator = "Nutrient tolerance",
  NullModel = "Additive",
  Term = nutrient_tolerance_levels,
  EffectType = "Level mean",
  ReferenceLevel = nutrient_tolerance_levels[1],
  N = nutrient_tolerance_N,
  n = nutrient_tolerance_n,
  Level_N = nutrient_tolerance_level_N,
  Level_n = nutrient_tolerance_level_n,
  Estimate = as.numeric(nutrient_tolerance_additive_level_CR2$beta),
  CR2_SE = as.numeric(nutrient_tolerance_additive_level_CR2$SE),
  CR2_df = as.numeric(nutrient_tolerance_additive_level_CR2$df_Satt),
  CI_lower = as.numeric(nutrient_tolerance_additive_level_lower),
  CI_upper = as.numeric(nutrient_tolerance_additive_level_upper),
  P_value = as.numeric(nutrient_tolerance_additive_level_CR2$p_Satt),
  Moderator_QM = as.numeric(nutrient_tolerance_additive_contrast_model$QM),
  Moderator_QM_P = as.numeric(nutrient_tolerance_additive_contrast_model$QMp),
  Tau2_Study = nutrient_tolerance_additive_contrast_model$sigma2[1],
  Tau2_Observation = nutrient_tolerance_additive_contrast_model$sigma2[2],
  Study_I2 = nutrient_tolerance_additive_study_I2,
  Observation_I2 = nutrient_tolerance_additive_observation_I2,
  Total_I2 = nutrient_tolerance_additive_total_I2,
  QE = as.numeric(nutrient_tolerance_additive_contrast_model$QE),
  QE_P_value = as.numeric(nutrient_tolerance_additive_contrast_model$QEp),
  Interpretation = ifelse(
    nutrient_tolerance_additive_level_lower > 0,
    "Positive non-additive",
    ifelse(
      nutrient_tolerance_additive_level_upper < 0,
      "Negative non-additive",
      "Additive"
    )
  ),
  Status = "Successfully fitted",
  stringsAsFactors = FALSE
)

nutrient_tolerance_additive_contrast_results <- data.frame(
  Moderator = "Nutrient tolerance",
  NullModel = "Additive",
  Term = paste0(
    nutrient_tolerance_levels[-1],
    " vs ",
    nutrient_tolerance_levels[1]
  ),
  EffectType = "Contrast with reference",
  ReferenceLevel = nutrient_tolerance_levels[1],
  N = nutrient_tolerance_N,
  n = nutrient_tolerance_n,
  Level_N = nutrient_tolerance_N,
  Level_n = nutrient_tolerance_n,
  Estimate = as.numeric(nutrient_tolerance_additive_contrast_CR2$beta),
  CR2_SE = as.numeric(nutrient_tolerance_additive_contrast_CR2$SE),
  CR2_df = as.numeric(nutrient_tolerance_additive_contrast_CR2$df_Satt),
  CI_lower = as.numeric(nutrient_tolerance_additive_contrast_lower),
  CI_upper = as.numeric(nutrient_tolerance_additive_contrast_upper),
  P_value = as.numeric(nutrient_tolerance_additive_contrast_CR2$p_Satt),
  Moderator_QM = as.numeric(nutrient_tolerance_additive_contrast_model$QM),
  Moderator_QM_P = as.numeric(nutrient_tolerance_additive_contrast_model$QMp),
  Tau2_Study = nutrient_tolerance_additive_contrast_model$sigma2[1],
  Tau2_Observation = nutrient_tolerance_additive_contrast_model$sigma2[2],
  Study_I2 = nutrient_tolerance_additive_study_I2,
  Observation_I2 = nutrient_tolerance_additive_observation_I2,
  Total_I2 = nutrient_tolerance_additive_total_I2,
  QE = as.numeric(nutrient_tolerance_additive_contrast_model$QE),
  QE_P_value = as.numeric(nutrient_tolerance_additive_contrast_model$QEp),
  Interpretation = ifelse(
    nutrient_tolerance_additive_contrast_lower > 0,
    "Higher than reference",
    ifelse(
      nutrient_tolerance_additive_contrast_upper < 0,
      "Lower than reference",
      "No detected difference"
    )
  ),
  Status = "Successfully fitted",
  stringsAsFactors = FALSE
)

nutrient_tolerance_additive_results <- dplyr::bind_rows(
  nutrient_tolerance_additive_level_results,
  nutrient_tolerance_additive_contrast_results
)


# Multiplicative null model
nutrient_tolerance_multiplicative_V <- clubSandwich::impute_covariance_matrix(
  vi = nutrient_tolerance_data$VdMultiplicative,
  cluster = nutrient_tolerance_data$StudyID,
  r = sampling_correlation,
  smooth_vi = FALSE,
  check_PD = TRUE,
  return_list = FALSE
)

nutrient_tolerance_multiplicative_level_model <- metafor::rma.mv(
  yi = HedgesdMultiplicative,
  V = nutrient_tolerance_multiplicative_V,
  mods = ~ 0 + ModeratorValue,
  random = list(
    ~ 1 | StudyID,
    ~ 1 | ObservationWithinStudy
  ),
  data = nutrient_tolerance_data,
  method = "REML",
  test = "t",
  sparse = TRUE
)

nutrient_tolerance_multiplicative_contrast_model <- metafor::rma.mv(
  yi = HedgesdMultiplicative,
  V = nutrient_tolerance_multiplicative_V,
  mods = ~ ModeratorValue,
  random = list(
    ~ 1 | StudyID,
    ~ 1 | ObservationWithinStudy
  ),
  data = nutrient_tolerance_data,
  method = "REML",
  test = "t",
  sparse = TRUE
)

nutrient_tolerance_multiplicative_level_CR2 <- as.data.frame(
  clubSandwich::coef_test(
    nutrient_tolerance_multiplicative_level_model,
    vcov = "CR2",
    cluster = nutrient_tolerance_data$StudyID,
    test = "Satterthwaite"
  )
)

nutrient_tolerance_multiplicative_contrast_CR2 <- as.data.frame(
  clubSandwich::coef_test(
    nutrient_tolerance_multiplicative_contrast_model,
    vcov = "CR2",
    cluster = nutrient_tolerance_data$StudyID,
    test = "Satterthwaite"
  )
)

nutrient_tolerance_multiplicative_level_critical <- stats::qt(
  0.975,
  df = nutrient_tolerance_multiplicative_level_CR2$df_Satt
)

nutrient_tolerance_multiplicative_level_lower <-
  nutrient_tolerance_multiplicative_level_CR2$beta -
  nutrient_tolerance_multiplicative_level_critical * nutrient_tolerance_multiplicative_level_CR2$SE

nutrient_tolerance_multiplicative_level_upper <-
  nutrient_tolerance_multiplicative_level_CR2$beta +
  nutrient_tolerance_multiplicative_level_critical * nutrient_tolerance_multiplicative_level_CR2$SE

nutrient_tolerance_multiplicative_contrast_CR2 <- nutrient_tolerance_multiplicative_contrast_CR2[
  -1,
  ,
  drop = FALSE
]

nutrient_tolerance_multiplicative_contrast_critical <- stats::qt(
  0.975,
  df = nutrient_tolerance_multiplicative_contrast_CR2$df_Satt
)

nutrient_tolerance_multiplicative_contrast_lower <-
  nutrient_tolerance_multiplicative_contrast_CR2$beta -
  nutrient_tolerance_multiplicative_contrast_critical * nutrient_tolerance_multiplicative_contrast_CR2$SE

nutrient_tolerance_multiplicative_contrast_upper <-
  nutrient_tolerance_multiplicative_contrast_CR2$beta +
  nutrient_tolerance_multiplicative_contrast_critical * nutrient_tolerance_multiplicative_contrast_CR2$SE

nutrient_tolerance_multiplicative_W <- solve(nutrient_tolerance_multiplicative_V)
nutrient_tolerance_multiplicative_X <- as.matrix(nutrient_tolerance_multiplicative_contrast_model$X)
nutrient_tolerance_multiplicative_P <-
  nutrient_tolerance_multiplicative_W -
  nutrient_tolerance_multiplicative_W %*% nutrient_tolerance_multiplicative_X %*%
  solve(crossprod(nutrient_tolerance_multiplicative_X, nutrient_tolerance_multiplicative_W %*% nutrient_tolerance_multiplicative_X)) %*%
  t(nutrient_tolerance_multiplicative_X) %*% nutrient_tolerance_multiplicative_W

nutrient_tolerance_multiplicative_typical_sampling_variance <-
  (nrow(nutrient_tolerance_multiplicative_X) - ncol(nutrient_tolerance_multiplicative_X)) /
  sum(diag(nutrient_tolerance_multiplicative_P))

nutrient_tolerance_multiplicative_total_variance <-
  sum(nutrient_tolerance_multiplicative_contrast_model$sigma2) +
  nutrient_tolerance_multiplicative_typical_sampling_variance

nutrient_tolerance_multiplicative_study_I2 <-
  100 * nutrient_tolerance_multiplicative_contrast_model$sigma2[1] /
  nutrient_tolerance_multiplicative_total_variance

nutrient_tolerance_multiplicative_observation_I2 <-
  100 * nutrient_tolerance_multiplicative_contrast_model$sigma2[2] /
  nutrient_tolerance_multiplicative_total_variance

nutrient_tolerance_multiplicative_total_I2 <-
  100 * sum(nutrient_tolerance_multiplicative_contrast_model$sigma2) /
  nutrient_tolerance_multiplicative_total_variance

nutrient_tolerance_multiplicative_level_results <- data.frame(
  Moderator = "Nutrient tolerance",
  NullModel = "Multiplicative",
  Term = nutrient_tolerance_levels,
  EffectType = "Level mean",
  ReferenceLevel = nutrient_tolerance_levels[1],
  N = nutrient_tolerance_N,
  n = nutrient_tolerance_n,
  Level_N = nutrient_tolerance_level_N,
  Level_n = nutrient_tolerance_level_n,
  Estimate = as.numeric(nutrient_tolerance_multiplicative_level_CR2$beta),
  CR2_SE = as.numeric(nutrient_tolerance_multiplicative_level_CR2$SE),
  CR2_df = as.numeric(nutrient_tolerance_multiplicative_level_CR2$df_Satt),
  CI_lower = as.numeric(nutrient_tolerance_multiplicative_level_lower),
  CI_upper = as.numeric(nutrient_tolerance_multiplicative_level_upper),
  P_value = as.numeric(nutrient_tolerance_multiplicative_level_CR2$p_Satt),
  Moderator_QM = as.numeric(nutrient_tolerance_multiplicative_contrast_model$QM),
  Moderator_QM_P = as.numeric(nutrient_tolerance_multiplicative_contrast_model$QMp),
  Tau2_Study = nutrient_tolerance_multiplicative_contrast_model$sigma2[1],
  Tau2_Observation = nutrient_tolerance_multiplicative_contrast_model$sigma2[2],
  Study_I2 = nutrient_tolerance_multiplicative_study_I2,
  Observation_I2 = nutrient_tolerance_multiplicative_observation_I2,
  Total_I2 = nutrient_tolerance_multiplicative_total_I2,
  QE = as.numeric(nutrient_tolerance_multiplicative_contrast_model$QE),
  QE_P_value = as.numeric(nutrient_tolerance_multiplicative_contrast_model$QEp),
  Interpretation = ifelse(
    nutrient_tolerance_multiplicative_level_lower > 0,
    "Positive non-additive",
    ifelse(
      nutrient_tolerance_multiplicative_level_upper < 0,
      "Negative non-additive",
      "Additive"
    )
  ),
  Status = "Successfully fitted",
  stringsAsFactors = FALSE
)

nutrient_tolerance_multiplicative_contrast_results <- data.frame(
  Moderator = "Nutrient tolerance",
  NullModel = "Multiplicative",
  Term = paste0(
    nutrient_tolerance_levels[-1],
    " vs ",
    nutrient_tolerance_levels[1]
  ),
  EffectType = "Contrast with reference",
  ReferenceLevel = nutrient_tolerance_levels[1],
  N = nutrient_tolerance_N,
  n = nutrient_tolerance_n,
  Level_N = nutrient_tolerance_N,
  Level_n = nutrient_tolerance_n,
  Estimate = as.numeric(nutrient_tolerance_multiplicative_contrast_CR2$beta),
  CR2_SE = as.numeric(nutrient_tolerance_multiplicative_contrast_CR2$SE),
  CR2_df = as.numeric(nutrient_tolerance_multiplicative_contrast_CR2$df_Satt),
  CI_lower = as.numeric(nutrient_tolerance_multiplicative_contrast_lower),
  CI_upper = as.numeric(nutrient_tolerance_multiplicative_contrast_upper),
  P_value = as.numeric(nutrient_tolerance_multiplicative_contrast_CR2$p_Satt),
  Moderator_QM = as.numeric(nutrient_tolerance_multiplicative_contrast_model$QM),
  Moderator_QM_P = as.numeric(nutrient_tolerance_multiplicative_contrast_model$QMp),
  Tau2_Study = nutrient_tolerance_multiplicative_contrast_model$sigma2[1],
  Tau2_Observation = nutrient_tolerance_multiplicative_contrast_model$sigma2[2],
  Study_I2 = nutrient_tolerance_multiplicative_study_I2,
  Observation_I2 = nutrient_tolerance_multiplicative_observation_I2,
  Total_I2 = nutrient_tolerance_multiplicative_total_I2,
  QE = as.numeric(nutrient_tolerance_multiplicative_contrast_model$QE),
  QE_P_value = as.numeric(nutrient_tolerance_multiplicative_contrast_model$QEp),
  Interpretation = ifelse(
    nutrient_tolerance_multiplicative_contrast_lower > 0,
    "Higher than reference",
    ifelse(
      nutrient_tolerance_multiplicative_contrast_upper < 0,
      "Lower than reference",
      "No detected difference"
    )
  ),
  Status = "Successfully fitted",
  stringsAsFactors = FALSE
)

nutrient_tolerance_multiplicative_results <- dplyr::bind_rows(
  nutrient_tolerance_multiplicative_level_results,
  nutrient_tolerance_multiplicative_contrast_results
)


# Dominance null model
nutrient_tolerance_dominance_V <- clubSandwich::impute_covariance_matrix(
  vi = nutrient_tolerance_data$VdDominance,
  cluster = nutrient_tolerance_data$StudyID,
  r = sampling_correlation,
  smooth_vi = FALSE,
  check_PD = TRUE,
  return_list = FALSE
)

nutrient_tolerance_dominance_level_model <- metafor::rma.mv(
  yi = HedgesdDominance,
  V = nutrient_tolerance_dominance_V,
  mods = ~ 0 + ModeratorValue,
  random = list(
    ~ 1 | StudyID,
    ~ 1 | ObservationWithinStudy
  ),
  data = nutrient_tolerance_data,
  method = "REML",
  test = "t",
  sparse = TRUE
)

nutrient_tolerance_dominance_contrast_model <- metafor::rma.mv(
  yi = HedgesdDominance,
  V = nutrient_tolerance_dominance_V,
  mods = ~ ModeratorValue,
  random = list(
    ~ 1 | StudyID,
    ~ 1 | ObservationWithinStudy
  ),
  data = nutrient_tolerance_data,
  method = "REML",
  test = "t",
  sparse = TRUE
)

nutrient_tolerance_dominance_level_CR2 <- as.data.frame(
  clubSandwich::coef_test(
    nutrient_tolerance_dominance_level_model,
    vcov = "CR2",
    cluster = nutrient_tolerance_data$StudyID,
    test = "Satterthwaite"
  )
)

nutrient_tolerance_dominance_contrast_CR2 <- as.data.frame(
  clubSandwich::coef_test(
    nutrient_tolerance_dominance_contrast_model,
    vcov = "CR2",
    cluster = nutrient_tolerance_data$StudyID,
    test = "Satterthwaite"
  )
)

nutrient_tolerance_dominance_level_critical <- stats::qt(
  0.975,
  df = nutrient_tolerance_dominance_level_CR2$df_Satt
)

nutrient_tolerance_dominance_level_lower <-
  nutrient_tolerance_dominance_level_CR2$beta -
  nutrient_tolerance_dominance_level_critical * nutrient_tolerance_dominance_level_CR2$SE

nutrient_tolerance_dominance_level_upper <-
  nutrient_tolerance_dominance_level_CR2$beta +
  nutrient_tolerance_dominance_level_critical * nutrient_tolerance_dominance_level_CR2$SE

nutrient_tolerance_dominance_contrast_CR2 <- nutrient_tolerance_dominance_contrast_CR2[
  -1,
  ,
  drop = FALSE
]

nutrient_tolerance_dominance_contrast_critical <- stats::qt(
  0.975,
  df = nutrient_tolerance_dominance_contrast_CR2$df_Satt
)

nutrient_tolerance_dominance_contrast_lower <-
  nutrient_tolerance_dominance_contrast_CR2$beta -
  nutrient_tolerance_dominance_contrast_critical * nutrient_tolerance_dominance_contrast_CR2$SE

nutrient_tolerance_dominance_contrast_upper <-
  nutrient_tolerance_dominance_contrast_CR2$beta +
  nutrient_tolerance_dominance_contrast_critical * nutrient_tolerance_dominance_contrast_CR2$SE

nutrient_tolerance_dominance_W <- solve(nutrient_tolerance_dominance_V)
nutrient_tolerance_dominance_X <- as.matrix(nutrient_tolerance_dominance_contrast_model$X)
nutrient_tolerance_dominance_P <-
  nutrient_tolerance_dominance_W -
  nutrient_tolerance_dominance_W %*% nutrient_tolerance_dominance_X %*%
  solve(crossprod(nutrient_tolerance_dominance_X, nutrient_tolerance_dominance_W %*% nutrient_tolerance_dominance_X)) %*%
  t(nutrient_tolerance_dominance_X) %*% nutrient_tolerance_dominance_W

nutrient_tolerance_dominance_typical_sampling_variance <-
  (nrow(nutrient_tolerance_dominance_X) - ncol(nutrient_tolerance_dominance_X)) /
  sum(diag(nutrient_tolerance_dominance_P))

nutrient_tolerance_dominance_total_variance <-
  sum(nutrient_tolerance_dominance_contrast_model$sigma2) +
  nutrient_tolerance_dominance_typical_sampling_variance

nutrient_tolerance_dominance_study_I2 <-
  100 * nutrient_tolerance_dominance_contrast_model$sigma2[1] /
  nutrient_tolerance_dominance_total_variance

nutrient_tolerance_dominance_observation_I2 <-
  100 * nutrient_tolerance_dominance_contrast_model$sigma2[2] /
  nutrient_tolerance_dominance_total_variance

nutrient_tolerance_dominance_total_I2 <-
  100 * sum(nutrient_tolerance_dominance_contrast_model$sigma2) /
  nutrient_tolerance_dominance_total_variance

nutrient_tolerance_dominance_level_results <- data.frame(
  Moderator = "Nutrient tolerance",
  NullModel = "Dominance",
  Term = nutrient_tolerance_levels,
  EffectType = "Level mean",
  ReferenceLevel = nutrient_tolerance_levels[1],
  N = nutrient_tolerance_N,
  n = nutrient_tolerance_n,
  Level_N = nutrient_tolerance_level_N,
  Level_n = nutrient_tolerance_level_n,
  Estimate = as.numeric(nutrient_tolerance_dominance_level_CR2$beta),
  CR2_SE = as.numeric(nutrient_tolerance_dominance_level_CR2$SE),
  CR2_df = as.numeric(nutrient_tolerance_dominance_level_CR2$df_Satt),
  CI_lower = as.numeric(nutrient_tolerance_dominance_level_lower),
  CI_upper = as.numeric(nutrient_tolerance_dominance_level_upper),
  P_value = as.numeric(nutrient_tolerance_dominance_level_CR2$p_Satt),
  Moderator_QM = as.numeric(nutrient_tolerance_dominance_contrast_model$QM),
  Moderator_QM_P = as.numeric(nutrient_tolerance_dominance_contrast_model$QMp),
  Tau2_Study = nutrient_tolerance_dominance_contrast_model$sigma2[1],
  Tau2_Observation = nutrient_tolerance_dominance_contrast_model$sigma2[2],
  Study_I2 = nutrient_tolerance_dominance_study_I2,
  Observation_I2 = nutrient_tolerance_dominance_observation_I2,
  Total_I2 = nutrient_tolerance_dominance_total_I2,
  QE = as.numeric(nutrient_tolerance_dominance_contrast_model$QE),
  QE_P_value = as.numeric(nutrient_tolerance_dominance_contrast_model$QEp),
  Interpretation = ifelse(
    nutrient_tolerance_dominance_level_lower > 0,
    "Positive non-additive",
    ifelse(
      nutrient_tolerance_dominance_level_upper < 0,
      "Negative non-additive",
      "Additive"
    )
  ),
  Status = "Successfully fitted",
  stringsAsFactors = FALSE
)

nutrient_tolerance_dominance_contrast_results <- data.frame(
  Moderator = "Nutrient tolerance",
  NullModel = "Dominance",
  Term = paste0(
    nutrient_tolerance_levels[-1],
    " vs ",
    nutrient_tolerance_levels[1]
  ),
  EffectType = "Contrast with reference",
  ReferenceLevel = nutrient_tolerance_levels[1],
  N = nutrient_tolerance_N,
  n = nutrient_tolerance_n,
  Level_N = nutrient_tolerance_N,
  Level_n = nutrient_tolerance_n,
  Estimate = as.numeric(nutrient_tolerance_dominance_contrast_CR2$beta),
  CR2_SE = as.numeric(nutrient_tolerance_dominance_contrast_CR2$SE),
  CR2_df = as.numeric(nutrient_tolerance_dominance_contrast_CR2$df_Satt),
  CI_lower = as.numeric(nutrient_tolerance_dominance_contrast_lower),
  CI_upper = as.numeric(nutrient_tolerance_dominance_contrast_upper),
  P_value = as.numeric(nutrient_tolerance_dominance_contrast_CR2$p_Satt),
  Moderator_QM = as.numeric(nutrient_tolerance_dominance_contrast_model$QM),
  Moderator_QM_P = as.numeric(nutrient_tolerance_dominance_contrast_model$QMp),
  Tau2_Study = nutrient_tolerance_dominance_contrast_model$sigma2[1],
  Tau2_Observation = nutrient_tolerance_dominance_contrast_model$sigma2[2],
  Study_I2 = nutrient_tolerance_dominance_study_I2,
  Observation_I2 = nutrient_tolerance_dominance_observation_I2,
  Total_I2 = nutrient_tolerance_dominance_total_I2,
  QE = as.numeric(nutrient_tolerance_dominance_contrast_model$QE),
  QE_P_value = as.numeric(nutrient_tolerance_dominance_contrast_model$QEp),
  Interpretation = ifelse(
    nutrient_tolerance_dominance_contrast_lower > 0,
    "Higher than reference",
    ifelse(
      nutrient_tolerance_dominance_contrast_upper < 0,
      "Lower than reference",
      "No detected difference"
    )
  ),
  Status = "Successfully fitted",
  stringsAsFactors = FALSE
)

nutrient_tolerance_dominance_results <- dplyr::bind_rows(
  nutrient_tolerance_dominance_level_results,
  nutrient_tolerance_dominance_contrast_results
)

nutrient_tolerance_results <- dplyr::bind_rows(
  nutrient_tolerance_additive_results,
  nutrient_tolerance_multiplicative_results,
  nutrient_tolerance_dominance_results
)




# 05.2.8. Photosynthetic pathway ----------------------------------------------

photosynthetic_pathway_data <- moderator_data[
    !is.na(moderator_data$PhotosyntheticPathway) &
  trimws(as.character(moderator_data$PhotosyntheticPathway)) != "" &
    is.finite(moderator_data$HedgesdAdditive) &
    is.finite(moderator_data$VdAdditive) &
    moderator_data$VdAdditive > 0 &
    is.finite(moderator_data$HedgesdMultiplicative) &
    is.finite(moderator_data$VdMultiplicative) &
    moderator_data$VdMultiplicative > 0 &
    is.finite(moderator_data$HedgesdDominance) &
    is.finite(moderator_data$VdDominance) &
    moderator_data$VdDominance > 0 &
    !is.na(moderator_data$StudyID) &
    !is.na(moderator_data$ObservationID),
  ,
  drop = FALSE
]

photosynthetic_pathway_data$ModeratorValue <- droplevels(
  factor(photosynthetic_pathway_data$PhotosyntheticPathway)
)

photosynthetic_pathway_data$StudyID <- factor(
  photosynthetic_pathway_data$StudyID
)

photosynthetic_pathway_data$ObservationWithinStudy <- factor(
  paste(
    photosynthetic_pathway_data$StudyID,
    photosynthetic_pathway_data$ObservationID,
    sep = "_"
  )
)

photosynthetic_pathway_levels <- levels(
  photosynthetic_pathway_data$ModeratorValue
)

photosynthetic_pathway_N <- length(
  unique(photosynthetic_pathway_data$StudyID)
)

photosynthetic_pathway_n <- nrow(
  photosynthetic_pathway_data
)

photosynthetic_pathway_level_n_table <- as.data.frame(
  table(photosynthetic_pathway_data$ModeratorValue),
  stringsAsFactors = FALSE
)

photosynthetic_pathway_level_study_data <- unique(
  data.frame(
    Level = photosynthetic_pathway_data$ModeratorValue,
    StudyID = photosynthetic_pathway_data$StudyID
  )
)

photosynthetic_pathway_level_N_table <- aggregate(
  StudyID ~ Level,
  data = photosynthetic_pathway_level_study_data,
  FUN = length
)

photosynthetic_pathway_level_n <- as.integer(
  photosynthetic_pathway_level_n_table$Freq[
    match(
      photosynthetic_pathway_levels,
      as.character(photosynthetic_pathway_level_n_table$Var1)
    )
  ]
)

photosynthetic_pathway_level_N <- as.integer(
  photosynthetic_pathway_level_N_table$StudyID[
    match(
      photosynthetic_pathway_levels,
      as.character(photosynthetic_pathway_level_N_table$Level)
    )
  ]
)

# Additive null model
photosynthetic_pathway_additive_V <- clubSandwich::impute_covariance_matrix(
  vi = photosynthetic_pathway_data$VdAdditive,
  cluster = photosynthetic_pathway_data$StudyID,
  r = sampling_correlation,
  smooth_vi = FALSE,
  check_PD = TRUE,
  return_list = FALSE
)

photosynthetic_pathway_additive_level_model <- metafor::rma.mv(
  yi = HedgesdAdditive,
  V = photosynthetic_pathway_additive_V,
  mods = ~ 0 + ModeratorValue,
  random = list(
    ~ 1 | StudyID,
    ~ 1 | ObservationWithinStudy
  ),
  data = photosynthetic_pathway_data,
  method = "REML",
  test = "t",
  sparse = TRUE
)

photosynthetic_pathway_additive_contrast_model <- metafor::rma.mv(
  yi = HedgesdAdditive,
  V = photosynthetic_pathway_additive_V,
  mods = ~ ModeratorValue,
  random = list(
    ~ 1 | StudyID,
    ~ 1 | ObservationWithinStudy
  ),
  data = photosynthetic_pathway_data,
  method = "REML",
  test = "t",
  sparse = TRUE
)

photosynthetic_pathway_additive_level_CR2 <- as.data.frame(
  clubSandwich::coef_test(
    photosynthetic_pathway_additive_level_model,
    vcov = "CR2",
    cluster = photosynthetic_pathway_data$StudyID,
    test = "Satterthwaite"
  )
)

photosynthetic_pathway_additive_contrast_CR2 <- as.data.frame(
  clubSandwich::coef_test(
    photosynthetic_pathway_additive_contrast_model,
    vcov = "CR2",
    cluster = photosynthetic_pathway_data$StudyID,
    test = "Satterthwaite"
  )
)

photosynthetic_pathway_additive_level_critical <- stats::qt(
  0.975,
  df = photosynthetic_pathway_additive_level_CR2$df_Satt
)

photosynthetic_pathway_additive_level_lower <-
  photosynthetic_pathway_additive_level_CR2$beta -
  photosynthetic_pathway_additive_level_critical * photosynthetic_pathway_additive_level_CR2$SE

photosynthetic_pathway_additive_level_upper <-
  photosynthetic_pathway_additive_level_CR2$beta +
  photosynthetic_pathway_additive_level_critical * photosynthetic_pathway_additive_level_CR2$SE

photosynthetic_pathway_additive_contrast_CR2 <- photosynthetic_pathway_additive_contrast_CR2[
  -1,
  ,
  drop = FALSE
]

photosynthetic_pathway_additive_contrast_critical <- stats::qt(
  0.975,
  df = photosynthetic_pathway_additive_contrast_CR2$df_Satt
)

photosynthetic_pathway_additive_contrast_lower <-
  photosynthetic_pathway_additive_contrast_CR2$beta -
  photosynthetic_pathway_additive_contrast_critical * photosynthetic_pathway_additive_contrast_CR2$SE

photosynthetic_pathway_additive_contrast_upper <-
  photosynthetic_pathway_additive_contrast_CR2$beta +
  photosynthetic_pathway_additive_contrast_critical * photosynthetic_pathway_additive_contrast_CR2$SE

photosynthetic_pathway_additive_W <- solve(photosynthetic_pathway_additive_V)
photosynthetic_pathway_additive_X <- as.matrix(photosynthetic_pathway_additive_contrast_model$X)
photosynthetic_pathway_additive_P <-
  photosynthetic_pathway_additive_W -
  photosynthetic_pathway_additive_W %*% photosynthetic_pathway_additive_X %*%
  solve(crossprod(photosynthetic_pathway_additive_X, photosynthetic_pathway_additive_W %*% photosynthetic_pathway_additive_X)) %*%
  t(photosynthetic_pathway_additive_X) %*% photosynthetic_pathway_additive_W

photosynthetic_pathway_additive_typical_sampling_variance <-
  (nrow(photosynthetic_pathway_additive_X) - ncol(photosynthetic_pathway_additive_X)) /
  sum(diag(photosynthetic_pathway_additive_P))

photosynthetic_pathway_additive_total_variance <-
  sum(photosynthetic_pathway_additive_contrast_model$sigma2) +
  photosynthetic_pathway_additive_typical_sampling_variance

photosynthetic_pathway_additive_study_I2 <-
  100 * photosynthetic_pathway_additive_contrast_model$sigma2[1] /
  photosynthetic_pathway_additive_total_variance

photosynthetic_pathway_additive_observation_I2 <-
  100 * photosynthetic_pathway_additive_contrast_model$sigma2[2] /
  photosynthetic_pathway_additive_total_variance

photosynthetic_pathway_additive_total_I2 <-
  100 * sum(photosynthetic_pathway_additive_contrast_model$sigma2) /
  photosynthetic_pathway_additive_total_variance

photosynthetic_pathway_additive_level_results <- data.frame(
  Moderator = "Photosynthetic pathway",
  NullModel = "Additive",
  Term = photosynthetic_pathway_levels,
  EffectType = "Level mean",
  ReferenceLevel = photosynthetic_pathway_levels[1],
  N = photosynthetic_pathway_N,
  n = photosynthetic_pathway_n,
  Level_N = photosynthetic_pathway_level_N,
  Level_n = photosynthetic_pathway_level_n,
  Estimate = as.numeric(photosynthetic_pathway_additive_level_CR2$beta),
  CR2_SE = as.numeric(photosynthetic_pathway_additive_level_CR2$SE),
  CR2_df = as.numeric(photosynthetic_pathway_additive_level_CR2$df_Satt),
  CI_lower = as.numeric(photosynthetic_pathway_additive_level_lower),
  CI_upper = as.numeric(photosynthetic_pathway_additive_level_upper),
  P_value = as.numeric(photosynthetic_pathway_additive_level_CR2$p_Satt),
  Moderator_QM = as.numeric(photosynthetic_pathway_additive_contrast_model$QM),
  Moderator_QM_P = as.numeric(photosynthetic_pathway_additive_contrast_model$QMp),
  Tau2_Study = photosynthetic_pathway_additive_contrast_model$sigma2[1],
  Tau2_Observation = photosynthetic_pathway_additive_contrast_model$sigma2[2],
  Study_I2 = photosynthetic_pathway_additive_study_I2,
  Observation_I2 = photosynthetic_pathway_additive_observation_I2,
  Total_I2 = photosynthetic_pathway_additive_total_I2,
  QE = as.numeric(photosynthetic_pathway_additive_contrast_model$QE),
  QE_P_value = as.numeric(photosynthetic_pathway_additive_contrast_model$QEp),
  Interpretation = ifelse(
    photosynthetic_pathway_additive_level_lower > 0,
    "Positive non-additive",
    ifelse(
      photosynthetic_pathway_additive_level_upper < 0,
      "Negative non-additive",
      "Additive"
    )
  ),
  Status = "Successfully fitted",
  stringsAsFactors = FALSE
)

photosynthetic_pathway_additive_contrast_results <- data.frame(
  Moderator = "Photosynthetic pathway",
  NullModel = "Additive",
  Term = paste0(
    photosynthetic_pathway_levels[-1],
    " vs ",
    photosynthetic_pathway_levels[1]
  ),
  EffectType = "Contrast with reference",
  ReferenceLevel = photosynthetic_pathway_levels[1],
  N = photosynthetic_pathway_N,
  n = photosynthetic_pathway_n,
  Level_N = photosynthetic_pathway_N,
  Level_n = photosynthetic_pathway_n,
  Estimate = as.numeric(photosynthetic_pathway_additive_contrast_CR2$beta),
  CR2_SE = as.numeric(photosynthetic_pathway_additive_contrast_CR2$SE),
  CR2_df = as.numeric(photosynthetic_pathway_additive_contrast_CR2$df_Satt),
  CI_lower = as.numeric(photosynthetic_pathway_additive_contrast_lower),
  CI_upper = as.numeric(photosynthetic_pathway_additive_contrast_upper),
  P_value = as.numeric(photosynthetic_pathway_additive_contrast_CR2$p_Satt),
  Moderator_QM = as.numeric(photosynthetic_pathway_additive_contrast_model$QM),
  Moderator_QM_P = as.numeric(photosynthetic_pathway_additive_contrast_model$QMp),
  Tau2_Study = photosynthetic_pathway_additive_contrast_model$sigma2[1],
  Tau2_Observation = photosynthetic_pathway_additive_contrast_model$sigma2[2],
  Study_I2 = photosynthetic_pathway_additive_study_I2,
  Observation_I2 = photosynthetic_pathway_additive_observation_I2,
  Total_I2 = photosynthetic_pathway_additive_total_I2,
  QE = as.numeric(photosynthetic_pathway_additive_contrast_model$QE),
  QE_P_value = as.numeric(photosynthetic_pathway_additive_contrast_model$QEp),
  Interpretation = ifelse(
    photosynthetic_pathway_additive_contrast_lower > 0,
    "Higher than reference",
    ifelse(
      photosynthetic_pathway_additive_contrast_upper < 0,
      "Lower than reference",
      "No detected difference"
    )
  ),
  Status = "Successfully fitted",
  stringsAsFactors = FALSE
)

photosynthetic_pathway_additive_results <- dplyr::bind_rows(
  photosynthetic_pathway_additive_level_results,
  photosynthetic_pathway_additive_contrast_results
)


# Multiplicative null model
photosynthetic_pathway_multiplicative_V <- clubSandwich::impute_covariance_matrix(
  vi = photosynthetic_pathway_data$VdMultiplicative,
  cluster = photosynthetic_pathway_data$StudyID,
  r = sampling_correlation,
  smooth_vi = FALSE,
  check_PD = TRUE,
  return_list = FALSE
)

photosynthetic_pathway_multiplicative_level_model <- metafor::rma.mv(
  yi = HedgesdMultiplicative,
  V = photosynthetic_pathway_multiplicative_V,
  mods = ~ 0 + ModeratorValue,
  random = list(
    ~ 1 | StudyID,
    ~ 1 | ObservationWithinStudy
  ),
  data = photosynthetic_pathway_data,
  method = "REML",
  test = "t",
  sparse = TRUE
)

photosynthetic_pathway_multiplicative_contrast_model <- metafor::rma.mv(
  yi = HedgesdMultiplicative,
  V = photosynthetic_pathway_multiplicative_V,
  mods = ~ ModeratorValue,
  random = list(
    ~ 1 | StudyID,
    ~ 1 | ObservationWithinStudy
  ),
  data = photosynthetic_pathway_data,
  method = "REML",
  test = "t",
  sparse = TRUE
)

photosynthetic_pathway_multiplicative_level_CR2 <- as.data.frame(
  clubSandwich::coef_test(
    photosynthetic_pathway_multiplicative_level_model,
    vcov = "CR2",
    cluster = photosynthetic_pathway_data$StudyID,
    test = "Satterthwaite"
  )
)

photosynthetic_pathway_multiplicative_contrast_CR2 <- as.data.frame(
  clubSandwich::coef_test(
    photosynthetic_pathway_multiplicative_contrast_model,
    vcov = "CR2",
    cluster = photosynthetic_pathway_data$StudyID,
    test = "Satterthwaite"
  )
)

photosynthetic_pathway_multiplicative_level_critical <- stats::qt(
  0.975,
  df = photosynthetic_pathway_multiplicative_level_CR2$df_Satt
)

photosynthetic_pathway_multiplicative_level_lower <-
  photosynthetic_pathway_multiplicative_level_CR2$beta -
  photosynthetic_pathway_multiplicative_level_critical * photosynthetic_pathway_multiplicative_level_CR2$SE

photosynthetic_pathway_multiplicative_level_upper <-
  photosynthetic_pathway_multiplicative_level_CR2$beta +
  photosynthetic_pathway_multiplicative_level_critical * photosynthetic_pathway_multiplicative_level_CR2$SE

photosynthetic_pathway_multiplicative_contrast_CR2 <- photosynthetic_pathway_multiplicative_contrast_CR2[
  -1,
  ,
  drop = FALSE
]

photosynthetic_pathway_multiplicative_contrast_critical <- stats::qt(
  0.975,
  df = photosynthetic_pathway_multiplicative_contrast_CR2$df_Satt
)

photosynthetic_pathway_multiplicative_contrast_lower <-
  photosynthetic_pathway_multiplicative_contrast_CR2$beta -
  photosynthetic_pathway_multiplicative_contrast_critical * photosynthetic_pathway_multiplicative_contrast_CR2$SE

photosynthetic_pathway_multiplicative_contrast_upper <-
  photosynthetic_pathway_multiplicative_contrast_CR2$beta +
  photosynthetic_pathway_multiplicative_contrast_critical * photosynthetic_pathway_multiplicative_contrast_CR2$SE

photosynthetic_pathway_multiplicative_W <- solve(photosynthetic_pathway_multiplicative_V)
photosynthetic_pathway_multiplicative_X <- as.matrix(photosynthetic_pathway_multiplicative_contrast_model$X)
photosynthetic_pathway_multiplicative_P <-
  photosynthetic_pathway_multiplicative_W -
  photosynthetic_pathway_multiplicative_W %*% photosynthetic_pathway_multiplicative_X %*%
  solve(crossprod(photosynthetic_pathway_multiplicative_X, photosynthetic_pathway_multiplicative_W %*% photosynthetic_pathway_multiplicative_X)) %*%
  t(photosynthetic_pathway_multiplicative_X) %*% photosynthetic_pathway_multiplicative_W

photosynthetic_pathway_multiplicative_typical_sampling_variance <-
  (nrow(photosynthetic_pathway_multiplicative_X) - ncol(photosynthetic_pathway_multiplicative_X)) /
  sum(diag(photosynthetic_pathway_multiplicative_P))

photosynthetic_pathway_multiplicative_total_variance <-
  sum(photosynthetic_pathway_multiplicative_contrast_model$sigma2) +
  photosynthetic_pathway_multiplicative_typical_sampling_variance

photosynthetic_pathway_multiplicative_study_I2 <-
  100 * photosynthetic_pathway_multiplicative_contrast_model$sigma2[1] /
  photosynthetic_pathway_multiplicative_total_variance

photosynthetic_pathway_multiplicative_observation_I2 <-
  100 * photosynthetic_pathway_multiplicative_contrast_model$sigma2[2] /
  photosynthetic_pathway_multiplicative_total_variance

photosynthetic_pathway_multiplicative_total_I2 <-
  100 * sum(photosynthetic_pathway_multiplicative_contrast_model$sigma2) /
  photosynthetic_pathway_multiplicative_total_variance

photosynthetic_pathway_multiplicative_level_results <- data.frame(
  Moderator = "Photosynthetic pathway",
  NullModel = "Multiplicative",
  Term = photosynthetic_pathway_levels,
  EffectType = "Level mean",
  ReferenceLevel = photosynthetic_pathway_levels[1],
  N = photosynthetic_pathway_N,
  n = photosynthetic_pathway_n,
  Level_N = photosynthetic_pathway_level_N,
  Level_n = photosynthetic_pathway_level_n,
  Estimate = as.numeric(photosynthetic_pathway_multiplicative_level_CR2$beta),
  CR2_SE = as.numeric(photosynthetic_pathway_multiplicative_level_CR2$SE),
  CR2_df = as.numeric(photosynthetic_pathway_multiplicative_level_CR2$df_Satt),
  CI_lower = as.numeric(photosynthetic_pathway_multiplicative_level_lower),
  CI_upper = as.numeric(photosynthetic_pathway_multiplicative_level_upper),
  P_value = as.numeric(photosynthetic_pathway_multiplicative_level_CR2$p_Satt),
  Moderator_QM = as.numeric(photosynthetic_pathway_multiplicative_contrast_model$QM),
  Moderator_QM_P = as.numeric(photosynthetic_pathway_multiplicative_contrast_model$QMp),
  Tau2_Study = photosynthetic_pathway_multiplicative_contrast_model$sigma2[1],
  Tau2_Observation = photosynthetic_pathway_multiplicative_contrast_model$sigma2[2],
  Study_I2 = photosynthetic_pathway_multiplicative_study_I2,
  Observation_I2 = photosynthetic_pathway_multiplicative_observation_I2,
  Total_I2 = photosynthetic_pathway_multiplicative_total_I2,
  QE = as.numeric(photosynthetic_pathway_multiplicative_contrast_model$QE),
  QE_P_value = as.numeric(photosynthetic_pathway_multiplicative_contrast_model$QEp),
  Interpretation = ifelse(
    photosynthetic_pathway_multiplicative_level_lower > 0,
    "Positive non-additive",
    ifelse(
      photosynthetic_pathway_multiplicative_level_upper < 0,
      "Negative non-additive",
      "Additive"
    )
  ),
  Status = "Successfully fitted",
  stringsAsFactors = FALSE
)

photosynthetic_pathway_multiplicative_contrast_results <- data.frame(
  Moderator = "Photosynthetic pathway",
  NullModel = "Multiplicative",
  Term = paste0(
    photosynthetic_pathway_levels[-1],
    " vs ",
    photosynthetic_pathway_levels[1]
  ),
  EffectType = "Contrast with reference",
  ReferenceLevel = photosynthetic_pathway_levels[1],
  N = photosynthetic_pathway_N,
  n = photosynthetic_pathway_n,
  Level_N = photosynthetic_pathway_N,
  Level_n = photosynthetic_pathway_n,
  Estimate = as.numeric(photosynthetic_pathway_multiplicative_contrast_CR2$beta),
  CR2_SE = as.numeric(photosynthetic_pathway_multiplicative_contrast_CR2$SE),
  CR2_df = as.numeric(photosynthetic_pathway_multiplicative_contrast_CR2$df_Satt),
  CI_lower = as.numeric(photosynthetic_pathway_multiplicative_contrast_lower),
  CI_upper = as.numeric(photosynthetic_pathway_multiplicative_contrast_upper),
  P_value = as.numeric(photosynthetic_pathway_multiplicative_contrast_CR2$p_Satt),
  Moderator_QM = as.numeric(photosynthetic_pathway_multiplicative_contrast_model$QM),
  Moderator_QM_P = as.numeric(photosynthetic_pathway_multiplicative_contrast_model$QMp),
  Tau2_Study = photosynthetic_pathway_multiplicative_contrast_model$sigma2[1],
  Tau2_Observation = photosynthetic_pathway_multiplicative_contrast_model$sigma2[2],
  Study_I2 = photosynthetic_pathway_multiplicative_study_I2,
  Observation_I2 = photosynthetic_pathway_multiplicative_observation_I2,
  Total_I2 = photosynthetic_pathway_multiplicative_total_I2,
  QE = as.numeric(photosynthetic_pathway_multiplicative_contrast_model$QE),
  QE_P_value = as.numeric(photosynthetic_pathway_multiplicative_contrast_model$QEp),
  Interpretation = ifelse(
    photosynthetic_pathway_multiplicative_contrast_lower > 0,
    "Higher than reference",
    ifelse(
      photosynthetic_pathway_multiplicative_contrast_upper < 0,
      "Lower than reference",
      "No detected difference"
    )
  ),
  Status = "Successfully fitted",
  stringsAsFactors = FALSE
)

photosynthetic_pathway_multiplicative_results <- dplyr::bind_rows(
  photosynthetic_pathway_multiplicative_level_results,
  photosynthetic_pathway_multiplicative_contrast_results
)


# Dominance null model
photosynthetic_pathway_dominance_V <- clubSandwich::impute_covariance_matrix(
  vi = photosynthetic_pathway_data$VdDominance,
  cluster = photosynthetic_pathway_data$StudyID,
  r = sampling_correlation,
  smooth_vi = FALSE,
  check_PD = TRUE,
  return_list = FALSE
)

photosynthetic_pathway_dominance_level_model <- metafor::rma.mv(
  yi = HedgesdDominance,
  V = photosynthetic_pathway_dominance_V,
  mods = ~ 0 + ModeratorValue,
  random = list(
    ~ 1 | StudyID,
    ~ 1 | ObservationWithinStudy
  ),
  data = photosynthetic_pathway_data,
  method = "REML",
  test = "t",
  sparse = TRUE
)

photosynthetic_pathway_dominance_contrast_model <- metafor::rma.mv(
  yi = HedgesdDominance,
  V = photosynthetic_pathway_dominance_V,
  mods = ~ ModeratorValue,
  random = list(
    ~ 1 | StudyID,
    ~ 1 | ObservationWithinStudy
  ),
  data = photosynthetic_pathway_data,
  method = "REML",
  test = "t",
  sparse = TRUE
)

photosynthetic_pathway_dominance_level_CR2 <- as.data.frame(
  clubSandwich::coef_test(
    photosynthetic_pathway_dominance_level_model,
    vcov = "CR2",
    cluster = photosynthetic_pathway_data$StudyID,
    test = "Satterthwaite"
  )
)

photosynthetic_pathway_dominance_contrast_CR2 <- as.data.frame(
  clubSandwich::coef_test(
    photosynthetic_pathway_dominance_contrast_model,
    vcov = "CR2",
    cluster = photosynthetic_pathway_data$StudyID,
    test = "Satterthwaite"
  )
)

photosynthetic_pathway_dominance_level_critical <- stats::qt(
  0.975,
  df = photosynthetic_pathway_dominance_level_CR2$df_Satt
)

photosynthetic_pathway_dominance_level_lower <-
  photosynthetic_pathway_dominance_level_CR2$beta -
  photosynthetic_pathway_dominance_level_critical * photosynthetic_pathway_dominance_level_CR2$SE

photosynthetic_pathway_dominance_level_upper <-
  photosynthetic_pathway_dominance_level_CR2$beta +
  photosynthetic_pathway_dominance_level_critical * photosynthetic_pathway_dominance_level_CR2$SE

photosynthetic_pathway_dominance_contrast_CR2 <- photosynthetic_pathway_dominance_contrast_CR2[
  -1,
  ,
  drop = FALSE
]

photosynthetic_pathway_dominance_contrast_critical <- stats::qt(
  0.975,
  df = photosynthetic_pathway_dominance_contrast_CR2$df_Satt
)

photosynthetic_pathway_dominance_contrast_lower <-
  photosynthetic_pathway_dominance_contrast_CR2$beta -
  photosynthetic_pathway_dominance_contrast_critical * photosynthetic_pathway_dominance_contrast_CR2$SE

photosynthetic_pathway_dominance_contrast_upper <-
  photosynthetic_pathway_dominance_contrast_CR2$beta +
  photosynthetic_pathway_dominance_contrast_critical * photosynthetic_pathway_dominance_contrast_CR2$SE

photosynthetic_pathway_dominance_W <- solve(photosynthetic_pathway_dominance_V)
photosynthetic_pathway_dominance_X <- as.matrix(photosynthetic_pathway_dominance_contrast_model$X)
photosynthetic_pathway_dominance_P <-
  photosynthetic_pathway_dominance_W -
  photosynthetic_pathway_dominance_W %*% photosynthetic_pathway_dominance_X %*%
  solve(crossprod(photosynthetic_pathway_dominance_X, photosynthetic_pathway_dominance_W %*% photosynthetic_pathway_dominance_X)) %*%
  t(photosynthetic_pathway_dominance_X) %*% photosynthetic_pathway_dominance_W

photosynthetic_pathway_dominance_typical_sampling_variance <-
  (nrow(photosynthetic_pathway_dominance_X) - ncol(photosynthetic_pathway_dominance_X)) /
  sum(diag(photosynthetic_pathway_dominance_P))

photosynthetic_pathway_dominance_total_variance <-
  sum(photosynthetic_pathway_dominance_contrast_model$sigma2) +
  photosynthetic_pathway_dominance_typical_sampling_variance

photosynthetic_pathway_dominance_study_I2 <-
  100 * photosynthetic_pathway_dominance_contrast_model$sigma2[1] /
  photosynthetic_pathway_dominance_total_variance

photosynthetic_pathway_dominance_observation_I2 <-
  100 * photosynthetic_pathway_dominance_contrast_model$sigma2[2] /
  photosynthetic_pathway_dominance_total_variance

photosynthetic_pathway_dominance_total_I2 <-
  100 * sum(photosynthetic_pathway_dominance_contrast_model$sigma2) /
  photosynthetic_pathway_dominance_total_variance

photosynthetic_pathway_dominance_level_results <- data.frame(
  Moderator = "Photosynthetic pathway",
  NullModel = "Dominance",
  Term = photosynthetic_pathway_levels,
  EffectType = "Level mean",
  ReferenceLevel = photosynthetic_pathway_levels[1],
  N = photosynthetic_pathway_N,
  n = photosynthetic_pathway_n,
  Level_N = photosynthetic_pathway_level_N,
  Level_n = photosynthetic_pathway_level_n,
  Estimate = as.numeric(photosynthetic_pathway_dominance_level_CR2$beta),
  CR2_SE = as.numeric(photosynthetic_pathway_dominance_level_CR2$SE),
  CR2_df = as.numeric(photosynthetic_pathway_dominance_level_CR2$df_Satt),
  CI_lower = as.numeric(photosynthetic_pathway_dominance_level_lower),
  CI_upper = as.numeric(photosynthetic_pathway_dominance_level_upper),
  P_value = as.numeric(photosynthetic_pathway_dominance_level_CR2$p_Satt),
  Moderator_QM = as.numeric(photosynthetic_pathway_dominance_contrast_model$QM),
  Moderator_QM_P = as.numeric(photosynthetic_pathway_dominance_contrast_model$QMp),
  Tau2_Study = photosynthetic_pathway_dominance_contrast_model$sigma2[1],
  Tau2_Observation = photosynthetic_pathway_dominance_contrast_model$sigma2[2],
  Study_I2 = photosynthetic_pathway_dominance_study_I2,
  Observation_I2 = photosynthetic_pathway_dominance_observation_I2,
  Total_I2 = photosynthetic_pathway_dominance_total_I2,
  QE = as.numeric(photosynthetic_pathway_dominance_contrast_model$QE),
  QE_P_value = as.numeric(photosynthetic_pathway_dominance_contrast_model$QEp),
  Interpretation = ifelse(
    photosynthetic_pathway_dominance_level_lower > 0,
    "Positive non-additive",
    ifelse(
      photosynthetic_pathway_dominance_level_upper < 0,
      "Negative non-additive",
      "Additive"
    )
  ),
  Status = "Successfully fitted",
  stringsAsFactors = FALSE
)

photosynthetic_pathway_dominance_contrast_results <- data.frame(
  Moderator = "Photosynthetic pathway",
  NullModel = "Dominance",
  Term = paste0(
    photosynthetic_pathway_levels[-1],
    " vs ",
    photosynthetic_pathway_levels[1]
  ),
  EffectType = "Contrast with reference",
  ReferenceLevel = photosynthetic_pathway_levels[1],
  N = photosynthetic_pathway_N,
  n = photosynthetic_pathway_n,
  Level_N = photosynthetic_pathway_N,
  Level_n = photosynthetic_pathway_n,
  Estimate = as.numeric(photosynthetic_pathway_dominance_contrast_CR2$beta),
  CR2_SE = as.numeric(photosynthetic_pathway_dominance_contrast_CR2$SE),
  CR2_df = as.numeric(photosynthetic_pathway_dominance_contrast_CR2$df_Satt),
  CI_lower = as.numeric(photosynthetic_pathway_dominance_contrast_lower),
  CI_upper = as.numeric(photosynthetic_pathway_dominance_contrast_upper),
  P_value = as.numeric(photosynthetic_pathway_dominance_contrast_CR2$p_Satt),
  Moderator_QM = as.numeric(photosynthetic_pathway_dominance_contrast_model$QM),
  Moderator_QM_P = as.numeric(photosynthetic_pathway_dominance_contrast_model$QMp),
  Tau2_Study = photosynthetic_pathway_dominance_contrast_model$sigma2[1],
  Tau2_Observation = photosynthetic_pathway_dominance_contrast_model$sigma2[2],
  Study_I2 = photosynthetic_pathway_dominance_study_I2,
  Observation_I2 = photosynthetic_pathway_dominance_observation_I2,
  Total_I2 = photosynthetic_pathway_dominance_total_I2,
  QE = as.numeric(photosynthetic_pathway_dominance_contrast_model$QE),
  QE_P_value = as.numeric(photosynthetic_pathway_dominance_contrast_model$QEp),
  Interpretation = ifelse(
    photosynthetic_pathway_dominance_contrast_lower > 0,
    "Higher than reference",
    ifelse(
      photosynthetic_pathway_dominance_contrast_upper < 0,
      "Lower than reference",
      "No detected difference"
    )
  ),
  Status = "Successfully fitted",
  stringsAsFactors = FALSE
)

photosynthetic_pathway_dominance_results <- dplyr::bind_rows(
  photosynthetic_pathway_dominance_level_results,
  photosynthetic_pathway_dominance_contrast_results
)

photosynthetic_pathway_results <- dplyr::bind_rows(
  photosynthetic_pathway_additive_results,
  photosynthetic_pathway_multiplicative_results,
  photosynthetic_pathway_dominance_results
)


# 05.3. Combine and inspect moderator results ---------------------------------

moderator_results <- dplyr::bind_rows(
  background_nutrient_results,
  fish_results,
  snail_results,
  ambient_temperature_results,
  study_duration_results,
  temperature_tolerance_results,
  nutrient_tolerance_results,
  photosynthetic_pathway_results
)

moderator_results$NullModel <- factor(
  moderator_results$NullModel,
  levels = c(
    "Additive",
    "Multiplicative",
    "Dominance"
  )
)

moderator_results <- moderator_results %>%
  dplyr::arrange(
    Moderator,
    NullModel,
    EffectType,
    Term
  )

moderator_contrast_results <- moderator_results %>%
  dplyr::filter(
    EffectType %in% c(
      "Contrast with reference",
      "Slope per SD"
    )
  )

moderator_coverage <- data.frame(
  Moderator = c(
    "Background nutrient level",
    "Fish presence",
    "Snail presence",
    "Ambient temperature",
    "Study duration",
    "Temperature tolerance",
    "Nutrient tolerance",
    "Photosynthetic pathway"
  ),
  N = c(
    background_nutrient_N,
    fish_N,
    snail_N,
    ambient_temperature_N,
    study_duration_N,
    temperature_tolerance_N,
    nutrient_tolerance_N,
    photosynthetic_pathway_N
  ),
  n = c(
    background_nutrient_n,
    fish_n,
    snail_n,
    ambient_temperature_n,
    study_duration_n,
    temperature_tolerance_n,
    nutrient_tolerance_n,
    photosynthetic_pathway_n
  ),
  stringsAsFactors = FALSE
)

moderator_level_coverage <- moderator_results %>%
  dplyr::filter(
    EffectType == "Level mean",
    NullModel == "Additive"
  ) %>%
  dplyr::select(
    Moderator,
    Level = Term,
    N = Level_N,
    n = Level_n
  )

print(moderator_coverage)
print(moderator_level_coverage)
print(moderator_contrast_results)
print(moderator_results)


# 05.4. Visualize moderators with significant CR2 tests -----------------------

significant_moderator_tests <- moderator_results %>%
  dplyr::filter(
    is.finite(P_value),
    P_value < 0.05,
    EffectType %in% c(
      "Level mean",
      "Slope per SD"
    )
  )

significant_moderator_names <- unique(
  as.character(
    significant_moderator_tests$Moderator
  )
)

significant_categorical_moderators <- unique(
  as.character(
    significant_moderator_tests$Moderator[
      significant_moderator_tests$EffectType ==
        "Level mean"
    ]
  )
)

significant_categorical_model_pairs <-
  significant_moderator_tests %>%
  dplyr::filter(
    EffectType == "Level mean"
  ) %>%
  dplyr::select(
    Moderator,
    NullModel
  ) %>%
  dplyr::distinct()

significant_continuous_moderators <- unique(
  as.character(
    significant_moderator_tests$Moderator[
      significant_moderator_tests$EffectType ==
        "Slope per SD"
    ]
  )
)


# Categorical moderators: display level-specific mean interaction effects
significant_moderator_plot_data <- moderator_results %>%
  dplyr::filter(
    Moderator %in%
      significant_categorical_moderators,
    EffectType == "Level mean"
  ) %>%
  dplyr::inner_join(
    significant_categorical_model_pairs,
    by = c(
      "Moderator",
      "NullModel"
    )
  ) %>%
  dplyr::mutate(
    DisplayGroup = dplyr::case_when(
      Moderator == "Fish presence" &
        Term == "N" ~
        "Without fish",
      Moderator == "Fish presence" &
        Term == "Y" ~
        "With fish",
      Moderator == "Snail presence" &
        Term == "N" ~
        "Without snail",
      Moderator == "Snail presence" &
        Term == "Y" ~
        "With snail",
      Term == "EutrophicSensitive" ~
        "Eutrophic\nsensitive",
      Term == "EutrophicTolerant" ~
        "Eutrophic\ntolerant",
      Term == "NotWarmTolerant" ~
        "Not warm\ntolerant",
      Term == "WarmTolerant" ~
        "Warm\ntolerant",
      Term == "NotC3" ~
        "Not C3",
      TRUE ~ as.character(Term)
    ),
    GroupLabel = paste0(
      DisplayGroup,
      "\n(n = ",
      Level_n,
      ")"
    )
  )

significant_group_order <- unique(
  significant_moderator_plot_data$GroupLabel
)

significant_moderator_plot_data$GroupLabel <- factor(
  significant_moderator_plot_data$GroupLabel,
  levels = significant_group_order
)

significant_moderator_plot_data$Moderator <- factor(
  significant_moderator_plot_data$Moderator,
  levels = c(
    "Background nutrient level",
    "Fish presence",
    "Snail presence",
    "Nutrient tolerance",
    "Photosynthetic pathway",
    "Temperature tolerance"
  )
)

significant_moderator_plot <- ggplot2::ggplot(
  significant_moderator_plot_data,
  ggplot2::aes(
    x = GroupLabel,
    y = Estimate,
    color = NullModel,
    shape = NullModel
  )
) +
  ggplot2::geom_hline(
    yintercept = 0,
    color = "lightgrey",
    linewidth = 0.8,
    linetype = "dashed"
  ) +
  ggplot2::geom_point(
    position = ggplot2::position_dodge(
      width = 0.45
    ),
    size = 3
  ) +
  ggplot2::geom_errorbar(
    ggplot2::aes(
      ymin = CI_lower,
      ymax = CI_upper
    ),
    position = ggplot2::position_dodge(
      width = 0.45
    ),
    width = 0.15,
    linewidth = 0.7
  ) +
  ggplot2::facet_wrap(
    ~ Moderator,
    nrow = 2,
    ncol = 3,
    scales = "free_x",
    labeller = ggplot2::labeller(
      Moderator = ggplot2::label_wrap_gen(width = 18)
    )
  ) +
  ggplot2::scale_color_manual(
    values = c(
      "Additive" = my_colours[1],
      "Multiplicative" = my_colours[2],
      "Dominance" = my_colours[3]
    ),
    labels = c(
      "Additive" = "Additive null model",
      "Multiplicative" = "Multiplicative null model",
      "Dominance" = "Dominance null model"
    ),
    name = NULL
  ) +
  ggplot2::scale_shape_manual(
    values = c(
      "Additive" = my_shapes[1],
      "Multiplicative" = my_shapes[2],
      "Dominance" = my_shapes[3]
    ),
    labels = c(
      "Additive" = "Additive null model",
      "Multiplicative" = "Multiplicative null model",
      "Dominance" = "Dominance null model"
    ),
    name = NULL
  ) +
  ggplot2::labs(
    x = "Moderator levels",
    y = "Interaction effect (Hedges' d)"
  ) +
  ggplot2::theme_bw(
    base_size = 13
  ) +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    strip.background =
      ggplot2::element_blank(),
    strip.text =
      ggplot2::element_text(
        size = 13,
        face = "bold"
      ),
    axis.text.x =
      ggplot2::element_text(
        size = 10,
        angle = 18,
        hjust = 1,
        vjust = 1
      ),
    legend.position = "bottom",
    legend.title =
      ggplot2::element_text(
        size = 12
      ),
    legend.text =
      ggplot2::element_text(
        size = 11
      )
  )

print(significant_moderator_tests)
print(significant_moderator_plot)


# Continuous moderators are kept separate because a standardized slope is not
# directly comparable with a categorical level mean.
significant_continuous_moderator_plot_data <-
  moderator_contrast_results %>%
  dplyr::filter(
    Moderator %in%
      significant_continuous_moderators,
    EffectType == "Slope per SD",
    P_value < 0.05
  )

if (
  nrow(
    significant_continuous_moderator_plot_data
  ) > 0
) {
  
  significant_continuous_moderator_plot <-
    ggplot2::ggplot(
      significant_continuous_moderator_plot_data,
      ggplot2::aes(
        x = Moderator,
        y = Estimate,
        color = NullModel,
        shape = NullModel
      )
    ) +
    ggplot2::geom_hline(
      yintercept = 0,
      color = "lightgrey",
      linewidth = 0.8,
      linetype = "dashed"
    ) +
    ggplot2::geom_point(
      position =
        ggplot2::position_dodge(
          width = 0.45
        ),
      size = 3
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        ymin = CI_lower,
        ymax = CI_upper
      ),
      position =
        ggplot2::position_dodge(
          width = 0.45
        ),
      width = 0.15,
      linewidth = 0.7
    ) +
    ggplot2::scale_color_manual(
      values = c(
        "Additive" = my_colours[1],
        "Multiplicative" = my_colours[2],
        "Dominance" = my_colours[3]
      ),
      name = "Null model"
    ) +
    ggplot2::scale_shape_manual(
      values = c(
        "Additive" = my_shapes[1],
        "Multiplicative" = my_shapes[2],
        "Dominance" = my_shapes[3]
      ),
      name = "Null model"
    ) +
    ggplot2::labs(
      x = "Moderator",
      y = "Slope per SD"
    ) +
    ggplot2::theme_bw(
      base_size = 13
    ) +
    ggplot2::theme(
      panel.grid =
        ggplot2::element_blank(),
      legend.position = "bottom"
    )
  
  print(
    significant_continuous_moderator_plot
  )
}


# 06. Full-model selection using ML and AICc ---------------------------------

model_selection_data <- moderator_data

model_selection_data$AmbientTemperatureScaled <- as.numeric(
  scale(model_selection_data$AmbientTemperature)
)

model_selection_data$StudyDurationScaled <- as.numeric(
  scale(model_selection_data$StudyDurationDays)
)

model_selection_columns <- c(
  "StudyID",
  "ObservationID",
  "FactorPair",
  "BackgroundNutrientLevel",
  "FishPresence",
  "SnailPresence",
  "AmbientTemperatureScaled",
  "StudyDurationScaled",
  "TemperatureTolerance",
  "NutrientTolerance",
  "PhotosyntheticPathway",
  "HedgesdAdditive",
  "VdAdditive",
  "HedgesdMultiplicative",
  "VdMultiplicative",
  "HedgesdDominance",
  "VdDominance"
)

model_selection_data <- model_selection_data[
  stats::complete.cases(
    model_selection_data[
      model_selection_columns
    ]
  ),
  ,
  drop = FALSE
]

model_selection_data$StudyID <- factor(
  model_selection_data$StudyID
)

model_selection_data$ObservationID <- factor(
  model_selection_data$ObservationID
)

model_selection_data$FactorPair <- droplevels(
  factor(model_selection_data$FactorPair)
)

model_selection_data$BackgroundNutrientLevel <- droplevels(
  factor(model_selection_data$BackgroundNutrientLevel)
)

model_selection_data$FishPresence <- droplevels(
  factor(model_selection_data$FishPresence)
)

model_selection_data$SnailPresence <- droplevels(
  factor(model_selection_data$SnailPresence)
)

model_selection_data$TemperatureTolerance <- droplevels(
  factor(model_selection_data$TemperatureTolerance)
)

model_selection_data$NutrientTolerance <- droplevels(
  factor(model_selection_data$NutrientTolerance)
)

model_selection_data$PhotosyntheticPathway <- droplevels(
  factor(model_selection_data$PhotosyntheticPathway)
)

model_selection_coverage <- data.frame(
  N = length(unique(model_selection_data$StudyID)),
  n = nrow(model_selection_data)
)

model_selection_formula <-
  ~ FactorPair +
    BackgroundNutrientLevel +
    FishPresence +
    SnailPresence +
    AmbientTemperatureScaled +
    StudyDurationScaled +
    TemperatureTolerance +
    NutrientTolerance +
    PhotosyntheticPathway

V_selection_additive <- clubSandwich::impute_covariance_matrix(
  vi = model_selection_data$VdAdditive,
  cluster = model_selection_data$StudyID,
  r = 0.6,
  smooth_vi = FALSE,
  check_PD = TRUE,
  return_list = FALSE
)

V_selection_multiplicative <- clubSandwich::impute_covariance_matrix(
  vi = model_selection_data$VdMultiplicative,
  cluster = model_selection_data$StudyID,
  r = 0.6,
  smooth_vi = FALSE,
  check_PD = TRUE,
  return_list = FALSE
)

V_selection_dominance <- clubSandwich::impute_covariance_matrix(
  vi = model_selection_data$VdDominance,
  cluster = model_selection_data$StudyID,
  r = 0.6,
  smooth_vi = FALSE,
  check_PD = TRUE,
  return_list = FALSE
)

# MuMIn requires na.fail so that every candidate is fitted to the same rows.
previous_na_action <- getOption("na.action")
options(na.action = "na.fail")

eval(metafor:::.MuMIn)

full_model_additive_ML <- metafor::rma.mv(
  yi = HedgesdAdditive,
  V = V_selection_additive,
  mods = model_selection_formula,
  random = ~ 1 | StudyID / ObservationID,
  data = model_selection_data,
  test = "t",
  method = "ML"
)

additive_rankings <- MuMIn::dredge(
  full_model_additive_ML
)

additive_candidate_models <- as.data.frame(
  subset(
    additive_rankings,
    delta < 3,
    recalc.weights = FALSE
  )
)

additive_importance_raw <- as.data.frame(
  MuMIn::sw(additive_rankings)
)

additive_importance <- data.frame(
  Variable = rownames(additive_importance_raw),
  RelativeImportance = as.numeric(additive_importance_raw[[1]]),
  Model = "Additive null model",
  row.names = NULL
)

full_model_multiplicative_ML <- metafor::rma.mv(
  yi = HedgesdMultiplicative,
  V = V_selection_multiplicative,
  mods = model_selection_formula,
  random = ~ 1 | StudyID / ObservationID,
  data = model_selection_data,
  test = "t",
  method = "ML"
)

multiplicative_rankings <- MuMIn::dredge(
  full_model_multiplicative_ML
)

multiplicative_candidate_models <- as.data.frame(
  subset(
    multiplicative_rankings,
    delta < 3,
    recalc.weights = FALSE
  )
)

multiplicative_importance_raw <- as.data.frame(
  MuMIn::sw(multiplicative_rankings)
)

multiplicative_importance <- data.frame(
  Variable = rownames(multiplicative_importance_raw),
  RelativeImportance = as.numeric(multiplicative_importance_raw[[1]]),
  Model = "Multiplicative null model",
  row.names = NULL
)

full_model_dominance_ML <- metafor::rma.mv(
  yi = HedgesdDominance,
  V = V_selection_dominance,
  mods = model_selection_formula,
  random = ~ 1 | StudyID / ObservationID,
  data = model_selection_data,
  test = "t",
  method = "ML"
)

dominance_rankings <- MuMIn::dredge(
  full_model_dominance_ML
)

dominance_candidate_models <- as.data.frame(
  subset(
    dominance_rankings,
    delta < 3,
    recalc.weights = FALSE
  )
)

additive_candidate_models$NullModel <- "Additive"
multiplicative_candidate_models$NullModel <- "Multiplicative"
dominance_candidate_models$NullModel <- "Dominance"

model_selection_candidates <- rbind(
  additive_candidate_models,
  multiplicative_candidate_models,
  dominance_candidate_models
)

dominance_importance_raw <- as.data.frame(
  MuMIn::sw(dominance_rankings)
)

dominance_importance <- data.frame(
  Variable = rownames(dominance_importance_raw),
  RelativeImportance = as.numeric(dominance_importance_raw[[1]]),
  Model = "Dominance null model",
  row.names = NULL
)

options(na.action = previous_na_action)

model_importance_results <- dplyr::bind_rows(
  additive_importance,
  multiplicative_importance,
  dominance_importance
)

model_importance_results$Model <- factor(
  model_importance_results$Model,
  levels = c(
    "Dominance null model",
    "Multiplicative null model",
    "Additive null model"
  )
)

# Use the manuscript's grouped horizontal bar-chart layout. Display labels are
# separated into words while the original variable names remain in the output
# table for reproducibility.
model_importance_plot_data <- model_importance_results %>%
  dplyr::mutate(
    DisplayVariable = dplyr::recode(
      Variable,
      "StudyDurationScaled" = "Study duration days",
      "AmbientTemperatureScaled" = "Ambient temperature",
      "TemperatureTolerance" = "Temperature tolerance",
      "PhotosyntheticPathway" = "Photosynthetic pathway",
      "SnailPresence" = "Snail presence",
      "NutrientTolerance" = "Nutrient tolerance",
      "FishPresence" = "Fish presence",
      "BackgroundNutrientLevel" = "Background nutrient level",
      "FactorPair" = "Stressor pair"
    ),
    DisplayVariable = factor(
      DisplayVariable,
      levels = c(
        "Stressor pair",
        "Background nutrient level",
        "Fish presence",
        "Nutrient tolerance",
        "Snail presence",
        "Photosynthetic pathway",
        "Temperature tolerance",
        "Ambient temperature",
        "Study duration days"
      )
    )
  )

model_selection_plot <- ggplot2::ggplot(
  model_importance_plot_data,
  ggplot2::aes(
    x = RelativeImportance,
    y = DisplayVariable,
    fill = Model
  )
) +
  ggplot2::geom_col(
    position = ggplot2::position_dodge(width = 0.78),
    width = 0.72,
    colour = "#333333",
    linewidth = 0.25
  ) +
  ggplot2::scale_fill_manual(
    values = c(
      "Additive null model" = my_colours[1],
      "Multiplicative null model" = my_colours[2],
      "Dominance null model" = my_colours[3]
    ),
    breaks = c(
      "Additive null model",
      "Multiplicative null model",
      "Dominance null model"
    ),
    name = NULL
  ) +
  ggplot2::scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2),
    expand = ggplot2::expansion(mult = c(0, 0.02))
  ) +
  ggplot2::theme_minimal(
    base_size = 14
  ) +
  ggplot2::theme(
    panel.grid.major = ggplot2::element_line(
      linetype = "dashed",
      color = "grey80",
      linewidth = 0.5
    ),
    panel.grid.minor = ggplot2::element_blank(),
    panel.border = ggplot2::element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.5
    ),
    axis.text = ggplot2::element_text(size = 12),
    legend.position = "bottom",
    legend.text = ggplot2::element_text(size = 12),
    axis.title.y = ggplot2::element_blank(),
    axis.title.x = ggplot2::element_text(size = 14)
  ) +
  ggplot2::labs(
    x = "Relative importance (sum of AICc weights)",
    y = NULL
  )

print(model_selection_plot)


# 07. Publication-bias diagnostics -------------------------------------------
#
# Following the reference analysis, small-study effects are examined using
# multilevel meta-regressions against sampling variance and inverse pooled
# sample size. Inference for the slope is study-clustered CR2.

publication_variance_additive <- metafor::rma.mv(
  yi = Hedgesd,
  V = V_additive,
  mods = ~ V,
  random = ~ 1 | StudyID / ObservationID,
  data = additive_data,
  test = "t",
  method = "REML"
)

publication_variance_additive_CR2 <- as.data.frame(
  clubSandwich::coef_test(
    publication_variance_additive,
    vcov = "CR2",
    cluster = additive_data$StudyID,
    test = "Satterthwaite"
  )
)

publication_sample_additive <- metafor::rma.mv(
  yi = Hedgesd,
  V = V_additive,
  mods = ~ I(1 / PooledSampleSize),
  random = ~ 1 | StudyID / ObservationID,
  data = additive_data,
  test = "t",
  method = "REML"
)

publication_sample_additive_CR2 <- as.data.frame(
  clubSandwich::coef_test(
    publication_sample_additive,
    vcov = "CR2",
    cluster = additive_data$StudyID,
    test = "Satterthwaite"
  )
)

publication_variance_multiplicative <- metafor::rma.mv(
  yi = Hedgesd,
  V = V_multiplicative,
  mods = ~ V,
  random = ~ 1 | StudyID / ObservationID,
  data = multiplicative_data,
  test = "t",
  method = "REML"
)

publication_variance_multiplicative_CR2 <- as.data.frame(
  clubSandwich::coef_test(
    publication_variance_multiplicative,
    vcov = "CR2",
    cluster = multiplicative_data$StudyID,
    test = "Satterthwaite"
  )
)

publication_sample_multiplicative <- metafor::rma.mv(
  yi = Hedgesd,
  V = V_multiplicative,
  mods = ~ I(1 / PooledSampleSize),
  random = ~ 1 | StudyID / ObservationID,
  data = multiplicative_data,
  test = "t",
  method = "REML"
)

publication_sample_multiplicative_CR2 <- as.data.frame(
  clubSandwich::coef_test(
    publication_sample_multiplicative,
    vcov = "CR2",
    cluster = multiplicative_data$StudyID,
    test = "Satterthwaite"
  )
)

publication_variance_dominance <- metafor::rma.mv(
  yi = Hedgesd,
  V = V_dominance,
  mods = ~ V,
  random = ~ 1 | StudyID / ObservationID,
  data = dominance_data,
  test = "t",
  method = "REML"
)

publication_variance_dominance_CR2 <- as.data.frame(
  clubSandwich::coef_test(
    publication_variance_dominance,
    vcov = "CR2",
    cluster = dominance_data$StudyID,
    test = "Satterthwaite"
  )
)

publication_sample_dominance <- metafor::rma.mv(
  yi = Hedgesd,
  V = V_dominance,
  mods = ~ I(1 / PooledSampleSize),
  random = ~ 1 | StudyID / ObservationID,
  data = dominance_data,
  test = "t",
  method = "REML"
)

publication_sample_dominance_CR2 <- as.data.frame(
  clubSandwich::coef_test(
    publication_sample_dominance,
    vcov = "CR2",
    cluster = dominance_data$StudyID,
    test = "Satterthwaite"
  )
)

publication_bias_results <- dplyr::bind_rows(
  data.frame(
    NullModel = "Additive",
    Diagnostic = "Sampling variance",
    publication_variance_additive_CR2[2, , drop = FALSE]
  ),
  data.frame(
    NullModel = "Additive",
    Diagnostic = "Inverse pooled sample size",
    publication_sample_additive_CR2[2, , drop = FALSE]
  ),
  data.frame(
    NullModel = "Multiplicative",
    Diagnostic = "Sampling variance",
    publication_variance_multiplicative_CR2[2, , drop = FALSE]
  ),
  data.frame(
    NullModel = "Multiplicative",
    Diagnostic = "Inverse pooled sample size",
    publication_sample_multiplicative_CR2[2, , drop = FALSE]
  ),
  data.frame(
    NullModel = "Dominance",
    Diagnostic = "Sampling variance",
    publication_variance_dominance_CR2[2, , drop = FALSE]
  ),
  data.frame(
    NullModel = "Dominance",
    Diagnostic = "Inverse pooled sample size",
    publication_sample_dominance_CR2[2, , drop = FALSE]
  )
)

# Funnel plots retain the appearance used in the original manuscript code.
grDevices::pdf(
  file.path(analysis_output_dir, "publication_bias_funnel_plots.pdf"),
  width = 8,
  height = 7
)

metafor::funnel(
  additive_model,
  main = "Additive null model",
  xlab = "Interaction effect size",
  xlim = c(-20, 10),
  cex.lab = 1,
  refline = 0,
  back = grDevices::adjustcolor(
    "cornflowerblue",
    alpha.f = 0.1
  ),
  col = "cornflowerblue",
  lwd = 3
)

metafor::funnel(
  multiplicative_model,
  main = "Multiplicative null model",
  xlab = "Interaction effect size",
  xlim = c(-20, 10),
  cex.lab = 1,
  refline = 0,
  back = grDevices::adjustcolor(
    "cornflowerblue",
    alpha.f = 0.1
  ),
  col = "cornflowerblue",
  lwd = 3
)

metafor::funnel(
  dominance_model,
  main = "Dominance null model",
  xlab = "Interaction effect size",
  xlim = c(-20, 10),
  cex.lab = 1,
  refline = 0,
  back = grDevices::adjustcolor(
    "cornflowerblue",
    alpha.f = 0.1
  ),
  col = "cornflowerblue",
  lwd = 3
)

grDevices::dev.off()


# 08. Sensitivity analyses -----------------------------------------------------

fit_sensitivity_model <- function(input_data, rho_value) {
  V_sensitivity <- clubSandwich::impute_covariance_matrix(
    vi = input_data$V,
    cluster = input_data$StudyID,
    r = rho_value,
    smooth_vi = FALSE,
    check_PD = TRUE,
    return_list = FALSE
  )

  fitted_model <- metafor::rma.mv(
    yi = Hedgesd,
    V = V_sensitivity,
    random = ~ 1 | StudyID / ObservationID,
    data = input_data,
    test = "t",
    method = "REML"
  )

  CR2_result <- as.data.frame(
    clubSandwich::coef_test(
      fitted_model,
      vcov = "CR2",
      cluster = input_data$StudyID,
      test = "Satterthwaite"
    )
  )

  critical_value <- stats::qt(0.975, CR2_result$df_Satt[1])

  data.frame(
    Estimate = CR2_result$beta[1],
    CR2_SE = CR2_result$SE[1],
    CR2_df = CR2_result$df_Satt[1],
    CI_lower = CR2_result$beta[1] - critical_value * CR2_result$SE[1],
    CI_upper = CR2_result$beta[1] + critical_value * CR2_result$SE[1],
    P_value = CR2_result$p_Satt[1],
    N = length(unique(input_data$StudyID)),
    n = nrow(input_data)
  )
}

sensitivity_model_data <- list(
  Additive = additive_data,
  Multiplicative = multiplicative_data,
  Dominance = dominance_data
)

sensitivity_results_list <- list()
sensitivity_result_number <- 1L

for (null_model_name in names(sensitivity_model_data)) {
  current_data <- sensitivity_model_data[[null_model_name]]

  for (rho_value in c(0, 0.3, 0.6, 0.9)) {
    current_result <- tryCatch(
      fit_sensitivity_model(current_data, rho_value),
      error = function(e) NULL
    )

    if (!is.null(current_result)) {
      current_result$NullModel <- null_model_name
      current_result$Analysis <- "Sampling correlation"
      current_result$Rho <- rho_value
      current_result$OmittedStudy <- NA_character_
      sensitivity_results_list[[sensitivity_result_number]] <- current_result
      sensitivity_result_number <- sensitivity_result_number + 1L
    }
  }

  for (study_to_omit in unique(current_data$StudyID)) {
    leave_one_out_data <- current_data[
      current_data$StudyID != study_to_omit,
      ,
      drop = FALSE
    ]

    current_result <- tryCatch(
      fit_sensitivity_model(leave_one_out_data, 0.6),
      error = function(e) NULL
    )

    if (!is.null(current_result)) {
      current_result$NullModel <- null_model_name
      current_result$Analysis <- "Leave one study out"
      current_result$Rho <- 0.6
      current_result$OmittedStudy <- as.character(study_to_omit)
      sensitivity_results_list[[sensitivity_result_number]] <- current_result
      sensitivity_result_number <- sensitivity_result_number + 1L
    }
  }
}

sensitivity_results <- dplyr::bind_rows(sensitivity_results_list) %>%
  dplyr::select(
    NullModel, Analysis, Rho, OmittedStudy,
    Estimate, CR2_SE, CR2_df, CI_lower, CI_upper, P_value, N, n
  )

sensitivity_summary <- sensitivity_results %>%
  dplyr::group_by(NullModel, Analysis) %>%
  dplyr::summarise(
    Runs = dplyr::n(),
    Estimate_min = min(Estimate, na.rm = TRUE),
    Estimate_max = max(Estimate, na.rm = TRUE),
    CI_lower_min = min(CI_lower, na.rm = TRUE),
    CI_upper_max = max(CI_upper, na.rm = TRUE),
    Significant_runs = sum(CI_lower > 0 | CI_upper < 0, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    NullModel = factor(
      NullModel,
      levels = c("Additive", "Multiplicative", "Dominance")
    ),
    Analysis = factor(
      Analysis,
      levels = c("Leave one study out", "Sampling correlation")
    )
  ) %>%
  dplyr::arrange(NullModel, Analysis) %>%
  dplyr::mutate(
    NullModel = as.character(NullModel),
    Analysis = as.character(Analysis)
  )

# 09. Reproduce manuscript figures and analysis tables ------------------------

write_csv_utf8 <- function(x, filename) {
  utils::write.csv(
    x,
    file.path(analysis_output_dir, filename),
    row.names = FALSE,
    na = "",
    fileEncoding = "UTF-8"
  )
}

# Analysis-count audit used to reconcile the PRISMA/data-flow figure with R.
analysis_flow_counts <- data.frame(
  Stage = c(
    "Current merged analysis workbook",
    "After removing invasive-stressor rows",
    "After removing non-global-change pairs",
    "Additive null model after effect-size checks",
    "Multiplicative null model after effect-size checks",
    "Dominance null model after effect-size checks",
    "Common three-model comparison",
    "Strictly positive multiplicative-assumption subset"
  ),
  N = c(
    input_dataset_N,
    after_invasive_N,
    length(unique(datafull$StudyID)),
    length(unique(additive_data$StudyID)),
    length(unique(multiplicative_data$StudyID)),
    length(unique(dominance_data$StudyID)),
    length(unique(three_model_data$StudyID)),
    length(unique(datafull$StudyID[
      datafull$MultiplicativeEligible &
        datafull$MultiplicativePositiveResponse
    ]))
  ),
  n = c(
    input_dataset_n,
    after_invasive_n,
    nrow(datafull),
    nrow(additive_data),
    nrow(multiplicative_data),
    nrow(dominance_data),
    nrow(three_model_data),
    sum(
      datafull$MultiplicativeEligible &
        datafull$MultiplicativePositiveResponse
    )
  ),
  stringsAsFactors = FALSE
)

screening_counts <- data.frame(
  Stage = c(
    "Original search records identified",
    "Original reviews/books removed",
    "Original title/abstract exclusions",
    "Original full-text reports assessed",
    "Original full-text exclusions",
    "Original studies included",
    "Update records identified and screened",
    "Update title/abstract exclusions",
    "Update reports sought",
    "Update reports not retrieved",
    "Update reports assessed",
    "Update full-text exclusions",
    "Update studies included",
    "Unique studies after merging searches",
    "Studies in additive/dominance R analysis",
    "Studies in three-model R comparison"
  ),
  Count = c(
    1584, 121, 1255, 208, 153, 55,
    195, 176, 19, 1, 18, 6, 12, 88,
    length(unique(additive_data$StudyID)),
    length(unique(three_model_data$StudyID))
  ),
  stringsAsFactors = FALSE
)

# Figure 2A: stressor-pair chord diagram with study and observation counts.
stressor_abbreviations <- c(
  "Chemical pollutant" = "CP",
  "CO2 enrichment" = "CE",
  "Hydraulic stressors" = "HS",
  "Light alterations" = "LA",
  "Elevated temperature" = "ET",
  "Nutrient enrichment" = "NE"
)

factor_pair_counts <- additive_source %>%
  dplyr::count(FactorPair, name = "Frequency")

pair_parts <- strsplit(as.character(factor_pair_counts$FactorPair), " x ")
factor_pair_counts$Factor1 <- vapply(pair_parts, `[`, character(1), 1)
factor_pair_counts$Factor2 <- vapply(pair_parts, `[`, character(1), 2)
factor_pair_counts$Factor1 <- unname(stressor_abbreviations[factor_pair_counts$Factor1])
factor_pair_counts$Factor2 <- unname(stressor_abbreviations[factor_pair_counts$Factor2])
factor_pair_counts <- factor_pair_counts[
  !is.na(factor_pair_counts$Factor1) &
    !is.na(factor_pair_counts$Factor2),
  ,
  drop = FALSE
]

stressor_codes <- unique(c(
  factor_pair_counts$Factor1,
  factor_pair_counts$Factor2
))

chord_matrix <- matrix(
  0,
  nrow = length(stressor_codes),
  ncol = length(stressor_codes),
  dimnames = list(stressor_codes, stressor_codes)
)

for (row_number in seq_len(nrow(factor_pair_counts))) {
  first_code <- factor_pair_counts$Factor1[row_number]
  second_code <- factor_pair_counts$Factor2[row_number]
  frequency <- factor_pair_counts$Frequency[row_number]
  chord_matrix[first_code, second_code] <-
    chord_matrix[first_code, second_code] + frequency
  if (first_code != second_code) {
    chord_matrix[second_code, first_code] <-
      chord_matrix[second_code, first_code] + frequency
  }
}

factor_totals <- rowSums(chord_matrix) + colSums(chord_matrix)
factor_order <- names(sort(factor_totals, decreasing = TRUE))
chord_colours <- c(
  "#1a3a9fff", "#15638dff", "#318288ff",
  "#45a383ff", "#69c46eff", "#b7da5fff"
)
names(chord_colours) <- factor_order

stressor_statistics <- lapply(
  factor_order,
  function(stressor_code) {
    stressor_name <- names(stressor_abbreviations)[
      stressor_abbreviations == stressor_code
    ]
    involved <- grepl(
      paste0("(^| x )", stressor_name, "($| x )"),
      additive_source$FactorPair
    )
    data.frame(
      Code = stressor_code,
      Stressor = stressor_name,
      N = length(unique(additive_source$StudyID[involved])),
      n = sum(involved)
    )
  }
)
stressor_statistics <- dplyr::bind_rows(stressor_statistics)

figure_2a_file <- file.path(analysis_output_dir, "Figure 2A.png")

grDevices::png(
  figure_2a_file,
  width = 3200,
  height = 1900,
  res = 300
)
graphics::layout(matrix(c(1, 2), nrow = 1), widths = c(1.25, 0.75))
graphics::par(mar = c(0, 0, 0, 0))
circlize::circos.clear()
circlize::circos.par(gap.degree = 3, start.degree = 135)
circlize::chordDiagram(
  chord_matrix[factor_order, factor_order],
  grid.col = chord_colours,
  annotationTrack = c("grid", "name"),
  order = factor_order,
  transparency = 0.28,
  preAllocateTracks = list(track.height = 0.18)
)
graphics::plot.new()
graphics::legend(
  "center",
  legend = paste0(
    stressor_statistics$Code,
    " = ", stressor_statistics$Stressor,
    "  [N=", stressor_statistics$N,
    ", n=", stressor_statistics$n, "]"
  ),
  fill = chord_colours[stressor_statistics$Code],
  title = "Stressor type",
  cex = 1.18,
  bty = "n",
  title.adj = 0.5
)
circlize::circos.clear()
grDevices::dev.off()

# Figure 2B: differences in interaction classifications among null models.
figure_2b_file <- file.path(analysis_output_dir, "Figure 2B.png")

grDevices::png(
  figure_2b_file,
  width = 3600,
  height = 1400,
  res = 300
)
grid::grid.draw(venn_figure)
grDevices::dev.off()

# Combine the chord and Venn panels in the arrangement used in the manuscript.
figure_2a_grob <- grid::rasterGrob(
  png::readPNG(figure_2a_file),
  interpolate = TRUE
)

figure_2b_grob <- grid::rasterGrob(
  png::readPNG(figure_2b_file),
  interpolate = TRUE
)

figure_2a_panel <- cowplot::ggdraw() +
  cowplot::draw_grob(figure_2a_grob) +
  cowplot::draw_label(
    "(A)",
    x = 0.012,
    y = 0.985,
    hjust = 0,
    vjust = 1,
    fontface = "bold",
    size = 16
  )

figure_2b_panel <- cowplot::ggdraw() +
  cowplot::draw_grob(figure_2b_grob) +
  cowplot::draw_label(
    "(B)",
    x = 0.012,
    y = 0.985,
    hjust = 0,
    vjust = 1,
    fontface = "bold",
    size = 16
  )

figure_2 <- cowplot::plot_grid(
  figure_2a_panel,
  figure_2b_panel,
  ncol = 1,
  rel_heights = c(1.15, 1)
)

ggplot2::ggsave(
  filename = file.path(analysis_output_dir, "Figure 2.png"),
  plot = figure_2,
  width = 13.0,
  height = 11.5,
  units = "in",
  dpi = 320,
  bg = "white"
)

# Figure 3: all observations and the CR2 overall mean under each null model.
ggplot2::ggsave(
  filename = file.path(
    analysis_output_dir,
    "Figure 3.png"
  ),
  plot = combined_plot,
  width = 13.3,
  height = 7.0,
  units = "in",
  dpi = 320,
  bg = "white"
)

# Figure 4: moderators with at least one significant CR2 level estimate.
if (
  nrow(significant_moderator_plot_data) > 0 ||
    exists("significant_continuous_moderator_plot")
) {
  if (
    nrow(significant_moderator_plot_data) > 0 &&
      exists("significant_continuous_moderator_plot")
  ) {
    significant_moderators_figure <-
      significant_moderator_plot /
      significant_continuous_moderator_plot +
      patchwork::plot_layout(heights = c(1.4, 1))
  } else if (nrow(significant_moderator_plot_data) > 0) {
    significant_moderators_figure <- significant_moderator_plot
  } else {
    significant_moderators_figure <-
      significant_continuous_moderator_plot
  }

  ggplot2::ggsave(
    filename = file.path(
      analysis_output_dir,
      "Figure 4.png"
    ),
    plot = significant_moderators_figure,
    width = 11.0,
    height = 9.5,
    units = "in",
    dpi = 320,
    bg = "white"
  )
}

# Figure 5: relative importance of each moderator under the three null models.
ggplot2::ggsave(
  filename = file.path(
    analysis_output_dir,
    "Figure 5.png"
  ),
  plot = model_selection_plot,
  width = 12.8,
  height = 7.0,
  units = "in",
  dpi = 320,
  bg = "white"
)

# Prepare a concise, readable candidate-model table for slide pagination.
candidate_metadata_columns <- c(
  "df", "logLik", "AICc", "delta", "weight", "NullModel"
)
candidate_term_columns <- setdiff(
  names(model_selection_candidates),
  candidate_metadata_columns
)
candidate_term_columns <- setdiff(candidate_term_columns, "(Intercept)")

candidate_terms <- apply(
  model_selection_candidates[, candidate_term_columns, drop = FALSE],
  1,
  function(x) {
    included <- candidate_term_columns[!is.na(x)]
    if (length(included) == 0) "Intercept only" else paste(included, collapse = "; ")
  }
)

model_selection_table <- model_selection_candidates %>%
  dplyr::mutate(Terms = candidate_terms) %>%
  dplyr::group_by(NullModel) %>%
  dplyr::arrange(delta, .by_group = TRUE) %>%
  dplyr::mutate(Rank = dplyr::row_number()) %>%
  dplyr::ungroup() %>%
  dplyr::transmute(
    NullModel = paste0(NullModel, " null model"),
    Rank,
    Terms,
    df,
    AICc,
    DeltaAICc = delta,
    Weight = weight,
    N = model_selection_coverage$N[1],
    n = model_selection_coverage$n[1]
  )

plant_traits <- readxl::read_xlsx(
  file.path(project_root, "plant traits.xlsx"),
  sheet = 1
)

write_csv_utf8(analysis_flow_counts, "analysis_flow_counts.csv")
write_csv_utf8(screening_counts, "screening_counts.csv")
write_csv_utf8(additive_floor_audit, "additive_prediction_floor_audit.csv")
write_csv_utf8(multiplicative_assumption_audit, "multiplicative_nonpositive_audit.csv")
write_csv_utf8(stressor_statistics, "stressor_counts.csv")
write_csv_utf8(overall_effects, "overall_effects.csv")
write_csv_utf8(moderator_results, "moderator_results.csv")
write_csv_utf8(significant_moderator_tests, "significant_moderator_tests.csv")
write_csv_utf8(model_importance_results, "model_importance_results.csv")
write_csv_utf8(model_selection_table, "model_selection_table.csv")
write_csv_utf8(sensitivity_results, "sensitivity_results.csv")
write_csv_utf8(sensitivity_summary, "sensitivity_summary.csv")
write_csv_utf8(publication_bias_results, "publication_bias_results.csv")
write_csv_utf8(plant_traits, "plant_traits.csv")

message(
  "Final analysis outputs written to: ",
  normalizePath(analysis_output_dir, winslash = "/", mustWork = FALSE)
)
