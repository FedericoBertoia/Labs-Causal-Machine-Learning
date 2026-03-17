########################
# 0. SETUP
########################
#  Packages 
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggcorrplot)
library(fastDummies)
library(SuperLearner)
library(pROC)

#  Reproducibility 
set.seed(123)

#  Load data 
# Set this to your local path
df_path <- "C:/Users/fbertoia/OneDrive - UGent/Desktop/GitHub/Labs-Causal-Machine-Learning/Lab 2 - Super Learner/student_dropout_dataset_v3.csv"

df <- read_csv(df_path)

# Subsample for speed during this lab
df <- slice_sample(df, n = 2000)

# Initial inspection 
glimpse(df)

# Preprocessing 
df <- df %>%
  # Drop identifier and redundant columns
  select(-Student_ID, -Semester_GPA, -CGPA) %>%
  # Remove incomplete cases
  drop_na()

# Confirm no remaining missingness
stopifnot(all(colSums(is.na(df)) == 0))


########################
# 1. Data Visualization
########################

# --- Histograms for categorical variables and target variable ---
df %>%
  select(where(is.character), Dropout) %>%
  mutate(Dropout = as.character(Dropout)) %>%
  pivot_longer(cols = everything(), names_to = "variable", values_to = "value") %>%
  ggplot(aes(x = value)) +
  geom_bar(fill = "steelblue") +
  facet_wrap(~variable, scales = "free") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = NULL, y = "Count")

# --- Histograms for categorical variables by dropout ---
df %>%
  select(where(is.character), Dropout) %>%
  pivot_longer(cols = -Dropout, names_to = "variable", values_to = "value") %>%
  ggplot(aes(x = value, fill = factor(Dropout))) +
  geom_bar(position = "fill") +
  facet_wrap(~variable, scales = "free") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = NULL, y = "Dropout Rate", fill = "Dropout") +
  scale_fill_manual(values = c("steelblue", "tomato"), labels = c("No", "Yes"))

# --- Distribution of GPA ---
ggplot(df, aes(x = GPA)) +
  geom_histogram(aes(y = after_stat(density)), bins = 40,
                 fill = "steelblue", alpha = 0.7) +
  geom_density(color = "navy", linewidth = 0.8) +
  theme_minimal() +
  labs(x = "GPA", y = "Density")

# --- Distribution of GPA by Dropout  ---
df %>%
  mutate(Dropout = factor(Dropout, labels = c("No", "Yes"))) %>%
  ggplot(aes(x = GPA, fill = Dropout, color = Dropout)) +
  geom_density(alpha = 0.3, linewidth = 0.8) +
  theme_minimal() +
  scale_fill_manual(values = c("steelblue", "tomato")) +
  scale_color_manual(values = c("steelblue", "tomato")) +
  labs(x = "GPA", y = "Density", fill = "Dropout", color = "Dropout")

# --- GPA by numeric predictors  ---
df %>%
  select(where(is.numeric), -Dropout) %>%
  pivot_longer(cols = -GPA, names_to = "variable", values_to = "value") %>%
  ggplot(aes(x = value, y = GPA)) +
  geom_point(alpha = 0.2, size = 0.8, color = "steelblue") +
  geom_smooth(method = "loess", se = FALSE, color = "tomato", linewidth = 0.8) +
  facet_wrap(~variable, scales = "free_x") +
  theme_minimal() +
  labs(x = NULL, y = "GPA")

########################
# 2. Categorical Encoding
########################

# Binary: Yes/No and Male/Female -> 0/1
df <- df %>%
  mutate(
    Gender          = ifelse(Gender == "Male", 1L, 0L),
    Internet_Access = ifelse(Internet_Access == "Yes", 1L, 0L),
    Part_Time_Job   = ifelse(Part_Time_Job == "Yes", 1L, 0L),
    Scholarship     = ifelse(Scholarship == "Yes", 1L, 0L)
  )

# Ordinal: factor levels -> integer (preserving order)
df <- df %>%
  mutate(
    Semester = factor(Semester,
                      levels = c("Year 1", "Year 2", "Year 3", "Year 4"),
                      ordered = TRUE) %>% as.integer(),
    Parental_Education = factor(Parental_Education,
                                levels = c("None", "High School", "Bachelor", "Master", "PhD"),
                                ordered = TRUE) %>% as.integer()
  )

# Nominal: one-hot encode Department (drop first dummy to avoid collinearity)
df <- dummy_cols(df,
                 select_columns          = "Department",
                 remove_first_dummy      = TRUE,
                 remove_selected_columns = TRUE)

########################
# 3. Correlation Plot
########################

# --- Correlation of each predictor with Dropout ---

df %>%
  cor() %>%
  as.data.frame() %>%
  select(Dropout) %>%
  tibble::rownames_to_column("variable") %>%
  filter(variable != "Dropout") %>%
  ggplot(aes(x = reorder(variable, Dropout), y = Dropout, fill = Dropout > 0)) +
  geom_col() +
  coord_flip() +
  theme_minimal() +
  scale_fill_manual(values = c("tomato", "steelblue"), labels = c("Negative", "Positive")) +
  labs(x = NULL, y = "Correlation with Dropout", fill = "Direction")

# --- Correlation of each predictor with GPA ---

df %>%
  select(-Dropout) %>%
  cor() %>%
  as.data.frame() %>%
  select(GPA) %>%
  tibble::rownames_to_column("variable") %>%
  filter(variable != "GPA") %>%
  ggplot(aes(x = reorder(variable, GPA), y = GPA, fill = GPA > 0)) +
  geom_col() +
  coord_flip() +
  theme_minimal() +
  scale_fill_manual(values = c("tomato", "steelblue"), labels = c("Negative", "Positive")) +
  labs(x = NULL, y = "Correlation with GPA", fill = "Direction")


################################################
# 4. Super Learner on the entire data
################################################

# --- Check available learners ---
listWrappers()

#############
# DROPOUT
#############

# Define feature matrix
X_dropout <- df %>% select(-Dropout)

# Candidate learners
SL.library <- c("SL.mean", "SL.glm", "SL.ranger", "SL.gam", "SL.earth")

# Fit Super Learner on full data
sl_fit_dropout <- SuperLearner(
  Y          = df$Dropout,
  X          = X_dropout,
  SL.library = SL.library,
  family     = binomial(),
  method     = "method.NNLS",
  cvControl  = list(V = 5)
)

sl_fit_dropout

# In-sample predictions
Y_hat_dropout <- predict(sl_fit_dropout, newdata = X_dropout)$pred

# MSE
mse_dropout <- mean((df$Dropout - Y_hat_dropout)^2)
cat("In-sample MSE:", round(mse_dropout, 4), "\n")

# AUC
roc_dropout <- roc(df$Dropout, Y_hat_dropout)
cat("In-sample AUC:", round(auc(roc_dropout), 4), "\n")

# ROC curve
plot(roc_dropout,
     col         = "steelblue",
     lwd         = 3,
     legacy.axes = TRUE,
     main        = paste("ROC Curve (full data) -- AUC:", round(auc(roc_dropout), 4)),
     xlab        = "False Positive Rate",
     ylab        = "True Positive Rate")
abline(a = 1, b = -1, lty = 2, col = "tomato", lwd = 2)


#############
# GPA
#############

# Define feature matrix
X_gpa <- df %>% select(-Dropout, -GPA)

# Fit Super Learner on full data
sl_fit_gpa <- SuperLearner(
  Y          = df$GPA,
  X          = X_gpa,
  SL.library = SL.library,
  family     = gaussian(),
  method     = "method.NNLS",
  cvControl  = list(V = 5)
)

sl_fit_gpa

# In-sample predictions
Y_hat_gpa <- predict(sl_fit_gpa, newdata = X_gpa)$pred

# MSE
mse_gpa <- mean((df$GPA - Y_hat_gpa)^2)
cat("In-sample MSE:", round(mse_gpa, 4), "\n")

# R-squared
r2_gpa <- 1 - mse_gpa / var(df$GPA)
cat("In-sample R-squared:", round(r2_gpa, 4), "\n")

# Predicted vs. actual
plot(df$GPA, Y_hat_gpa,
     col  = adjustcolor("steelblue", alpha.f = 0.4),
     pch  = 16,
     xlab = "Observed GPA",
     ylab = "Predicted GPA",
     main = paste("Predicted vs. Observed GPA -- R2:", round(r2_gpa, 4)),
     ylim = c(0, max(Y_hat_gpa)))
abline(a = 0, b = 1, lty = 2, col = "tomato", lwd = 2)

################################################
# 5. Cross-Validated Risk of the Super Learner
################################################

#############
# DROPOUT
#############

cv_sl_dropout <- CV.SuperLearner(
  Y              = df$Dropout,
  X              = X_dropout,
  SL.library     = SL.library,
  family         = binomial(),
  method         = "method.NNLS",
  cvControl      = list(V = 5),
  innerCvControl = list(list(V = 5))
)

summary(cv_sl_dropout)
plot(cv_sl_dropout)

#############
# GPA
#############

cv_sl_gpa <- CV.SuperLearner(
  Y              = df$GPA,
  X              = X_gpa,
  SL.library     = SL.library,
  family         = gaussian(),
  method         = "method.NNLS",
  cvControl      = list(V = 5),
  innerCvControl = list(list(V = 5))
)

summary(cv_sl_gpa)
plot(cv_sl_gpa)



################################################
# 6. Super Learner with Cross-Fitting
################################################

# --- Create folds ---
n_folds <- 5
N       <- nrow(df)
fold_id <- sample(rep(1:n_folds, length.out = N))
folds   <- split(seq_len(N), fold_id)

#############
# DROPOUT
#############

# Container for out-of-fold predictions
Y_hat_dropout_cf <- numeric(N)

for (k in 1:n_folds) {
  
  test  <- folds[[k]]
  train <- setdiff(seq_len(N), test)
  
  sl_fit_dropout_cf <- SuperLearner(
    Y          = df$Dropout[train],
    X          = X_dropout[train, ],
    SL.library = SL.library,
    family     = binomial(),
    method     = "method.NNLS",
    cvControl  = list(V = 5)
  )
  
  # Predict on held-out fold
  Y_hat_dropout_cf[test] <- predict(sl_fit_dropout_cf, newdata = X_dropout[test, ])$pred
}

# MSE
mse_dropout_cf <- mean((df$Dropout - Y_hat_dropout_cf)^2)
cat("Cross-fitted MSE:", round(mse_dropout_cf, 4), "\n")

# AUC
roc_dropout_cf <- roc(df$Dropout, Y_hat_dropout_cf)
cat("Cross-fitted AUC:", round(auc(roc_dropout_cf), 4), "\n")

# ROC curve
plot(roc_dropout_cf,
     col         = "steelblue",
     lwd         = 3,
     legacy.axes = TRUE,
     main        = paste("ROC Curve (cross-fitted) -- AUC:", round(auc(roc_dropout_cf), 4)),
     xlab        = "False Positive Rate",
     ylab        = "True Positive Rate")
abline(a = 1, b = -1, lty = 2, col = "tomato", lwd = 2)

#############
# GPA
#############

# Container for out-of-fold predictions
Y_hat_gpa_cf <- numeric(N)

for (k in 1:n_folds) {
  
  test  <- folds[[k]]
  train <- setdiff(seq_len(N), test)
  
  sl_fit_gpa_cf <- SuperLearner(
    Y          = df$GPA[train],
    X          = X_gpa[train, ],
    SL.library = SL.library,
    family     = gaussian(),
    method     = "method.NNLS",
    cvControl  = list(V = 5)
  )
  
  # Predict on held-out fold
  Y_hat_gpa_cf[test] <- predict(sl_fit_gpa_cf, newdata = X_gpa[test, ])$pred
}

# MSE
mse_gpa_cf <- mean((df$GPA - Y_hat_gpa_cf)^2)
cat("Cross-fitted MSE:", round(mse_gpa_cf, 4), "\n")

# R-squared
r2_gpa_cf <- 1 - mse_gpa_cf / var(df$GPA)
cat("Cross-fitted R-squared:", round(r2_gpa_cf, 4), "\n")

# Predicted vs. actual
plot(df$GPA, Y_hat_gpa_cf,
     col  = adjustcolor("steelblue", alpha.f = 0.4),
     pch  = 16,
     xlab = "Observed GPA",
     ylab = "Predicted GPA",
     main = paste("Predicted vs. Observed GPA (cross-fitted) -- R2:", round(r2_gpa_cf, 4)))
abline(a = 0, b = 1, lty = 2, col = "tomato", lwd = 2)


################################################
# 7. Tuning Base Learners with create.Learner()
################################################

#####################
# Define tuned learners
#####################

# Tuned random forests
learners_ranger <- create.Learner(
  "SL.ranger",
  tune = list(
    mtry      = floor(sqrt(ncol(X_dropout)) * c(0.5, 1, 2)),
    num.trees = c(100, 500)
  )
)

# Tuned MARS
learners_earth <- create.Learner(
  "SL.earth",
  tune = list(
    degree = c(1, 2)
  )
)

# Tuned xgboost
learners_xgboost <- create.Learner(
  "SL.xgboost",
  tune = list(
    max_depth = c(2, 4, 6)
  )
)

learners_ranger

learners_earth

learners_xgboost

# Expanded library
SL.library_tuned <- c(
  "SL.mean",
  "SL.glm",
  "SL.gam",
  learners_ranger$names,
  learners_earth$names,
  learners_xgboost$names
)

cat("Candidate learners:\n")
print(SL.library_tuned)

#####################
# DROPOUT
#####################

cv_sl_tuned_dropout <- CV.SuperLearner(
  Y              = df$Dropout,
  X              = X_dropout,
  SL.library     = SL.library_tuned,
  family         = binomial(),
  method         = "method.NNLS",
  cvControl      = list(V = 5),
  innerCvControl = list(list(V = 5))
)

summary(cv_sl_tuned_dropout)
plot(cv_sl_tuned_dropout)


#####################
# GPA
#####################

cv_sl_tuned_gpa <- CV.SuperLearner(
  Y              = df$GPA,
  X              = X_gpa,
  SL.library     = SL.library_tuned,
  family         = gaussian(),
  method         = "method.NNLS",
  cvControl      = list(V = 5),
  innerCvControl = list(list(V = 5))
)

summary(cv_sl_tuned_gpa)
plot(cv_sl_tuned_gpa)


################################################
# 8. Cross-Fitting with Tuned Library
################################################

#####################
# DROPOUT
#####################

Y_hat_dropout_cf_tuned <- numeric(N)

for (k in 1:n_folds) {
  
  test  <- folds[[k]]
  train <- setdiff(seq_len(N), test)
  
  sl_fit <- SuperLearner(
    Y          = df$Dropout[train],
    X          = X_dropout[train, ],
    SL.library = SL.library_tuned,
    family     = binomial(),
    method     = "method.NNLS",
    cvControl  = list(V = 5)
  )
  
  Y_hat_dropout_cf_tuned[test] <- predict(sl_fit, newdata = X_dropout[test, ])$pred
}

mse_dropout_cf_tuned <- mean((df$Dropout - Y_hat_dropout_cf_tuned)^2)
roc_dropout_cf_tuned <- roc(df$Dropout, Y_hat_dropout_cf_tuned)

cat("Cross-fitted MSE (tuned):", round(mse_dropout_cf_tuned, 4), "\n")
cat("Cross-fitted AUC (tuned):", round(auc(roc_dropout_cf_tuned), 4), "\n")

#####################
# GPA
#####################

Y_hat_gpa_cf_tuned <- numeric(N)

for (k in 1:n_folds) {
  
  test  <- folds[[k]]
  train <- setdiff(seq_len(N), test)
  
  sl_fit <- SuperLearner(
    Y          = df$GPA[train],
    X          = X_gpa[train, ],
    SL.library = SL.library_tuned,
    family     = gaussian(),
    method     = "method.NNLS",
    cvControl  = list(V = 5)
  )
  
  Y_hat_gpa_cf_tuned[test] <- predict(sl_fit, newdata = X_gpa[test, ])$pred
}

mse_gpa_cf_tuned    <- mean((df$GPA - Y_hat_gpa_cf_tuned)^2)
r2_gpa_cf_tuned     <- 1 - mse_gpa_cf_tuned / var(df$GPA)

cat("Cross-fitted MSE (tuned):", round(mse_gpa_cf_tuned, 4), "\n")
cat("Cross-fitted R-squared (tuned):", round(r2_gpa_cf_tuned, 4), "\n")



# --- Final Comparison ---

comparison <- data.frame(
  Outcome = c("Dropout", "Dropout", "GPA", "GPA"),
  Library = c("Default", "Tuned", "Default", "Tuned"),
  MSE     = c(mse_dropout_cf,      mse_dropout_cf_tuned,
              mse_gpa_cf,          mse_gpa_cf_tuned),
  Metric2 = c(round(auc(roc_dropout_cf), 4), round(auc(roc_dropout_cf_tuned), 4),
              r2_gpa_cf,                      r2_gpa_cf_tuned),
  Metric2_name = c("AUC", "AUC", "R2", "R2")
)

print(comparison)

################################################
# BONUS: Custom Learner Solution
################################################

# Tuned glmnet: grid over alpha (0 = ridge, 1 = lasso, in-between = elastic net)
# and lambda is chosen internally by glmnet via cross-validation
learners_glmnet <- create.Learner(
  "SL.glmnet",
  tune = list(alpha = c(0, 0.5, 1))
)

learners_glmnet$names  # SL.glmnet_1, SL.glmnet_2, SL.glmnet_3