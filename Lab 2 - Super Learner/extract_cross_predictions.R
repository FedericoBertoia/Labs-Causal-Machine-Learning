library(SuperLearner)
library(tidyverse)
#  Reproducibility 
set.seed(123)

#  Load data 
# Set this to your local path
df_path <- "C:/Users/fbertoia/OneDrive - UGent/Desktop/GitHub/Labs-Causal-Machine-Learning/Lab 2 - Super Learner/student_dropout_dataset_v3.csv"

df <- read_csv(df_path)

# Subsample for speed during this lab
df <- slice_sample(df, n = 2000)

# Preprocessing 
df <- df %>%
  # Drop identifier and redundant columns
  select(-Student_ID, -Semester_GPA, -CGPA) %>%
  # Remove incomplete cases
  drop_na()

# get X and Y

X_dropout <- df %>% select(-Dropout)
Y_dropout <- df %>% pull(Dropout)

# superlearner library

SL.library <- c( "SL.glm", "SL.mean", "SL.earth") # to demonstrate, I only use learners without a stochastic element

#######################################################
####### Interesting part starts here ##################-----------------------
#######################################################
set.seed(42)

# outer_CV{ SuperLearner }  --------------------------------------------------
cv_sl_dropout <- CV.SuperLearner(
  Y              = Y_dropout,
  X              = X_dropout,
  SL.library     = SL.library,
  family         = binomial(),   # binomial() or gaussian()?
  method         = "method.NNLS",
  cvControl      = list(V = 10),
  innerCvControl = list(list(V = 5))
)

# get SL.predict object
cv_preds <- cv_sl_dropout$SL.predict 

cv_folds <- cv_sl_dropout$folds

# manual outer_CV{ SuperLearner } --------------------------------------------
set.seed(42)
folds <- CVFolds(nrow(df), id = NULL, Y = df$Dropout, cvControl = SuperLearner.CV.control(V=10)) 

# check if my folds are the same as before:
all.equal(folds, cv_folds)


# --- Create folds ---
# Hint: use n_folds greater or equal to 2. What happens if we use 1? 
n_folds <- 10
N       <- nrow(df)

# Container for out-of-fold predictions
Y_hat_dropout_cf <- numeric(N)

for (k in 1:n_folds) {
  
  test  <- folds[[k]]
  train <- setdiff(seq_len(N), test)
  
  sl_fit_dropout_cf <- SuperLearner(
    Y          = df$Dropout[train],   # Hint: use only training indices
    X          = X_dropout[train,],   # Hint: use only training indices
    SL.library = SL.library,
    family     = binomial(),   # binomial() or gaussian()?
    method     = "method.NNLS",
    cvControl  = list(V = 5)  
  )
  
  # Predict on held-out fold
  Y_hat_dropout_cf[test] <- predict(sl_fit_dropout_cf, newdata = X_dropout[test, ])$pred
}

######### Result ################################
# same?
all.equal(Y_hat_dropout_cf, cv_preds) # yes




# easily access predictions using each fold
map(cv_folds, \(index) cv_preds[index])

# or attach to df
imap(cv_folds, \(index, i) {
  
  df[index,] %>% 
    mutate(fold = i,
           pred_A = cv_preds[index], .before = 1)

  }) %>% 
  list_rbind()
