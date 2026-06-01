# Labs — Causal Machine Learning

Lab materials and solutions for the **Causal Machine Learning** course (2026 edition).
Each lab pairs a short theoretical primer with a hands-on **R** session, working through the
core toolkit of modern causal inference — from flexible nuisance estimation with the
**Super Learner** to **doubly-robust** estimation of treatment effects, **mediation**,
**instrumental variables**, and **assumption-lean inference**.

---

## 🧪 Labs

| Lab | Topic | What it covers |
|-----|-------|----------------|
| **Lab 2** | [Super Learner](#lab-2--super-learner) | Stacked ensemble learning & cross-validated prediction |
| **Lab 5** | [Average Treatment Effect](#lab-5--average-treatment-effect) | Estimating the ATE (IPW, g-computation, doubly-robust / AIPW) |
| **Lab 6** | [Mediation Analysis](#lab-6--mediation-analysis) | Decomposing direct and indirect causal effects |
| **Lab 7** | Causal Prediction | Prediction under interventions *(materials kept locally)* |
| **Lab 8** | [Assumption-Lean Inference](#lab-8--assumption-lean-inference) | Robust inference with minimal modelling assumptions |
| **Lab 9** | [Instrumental Variables](#lab-9--instrumental-variables) | Identification & estimation with instruments |
| **Lab 10** | [Solving an Old Exam](#lab-10--solving-an-old-exam) | Worked past-exam solutions & common pitfalls |

---

### Lab 2 — Super Learner
Building a **Super Learner** (stacked ensemble) and using cross-validated predictions as a
foundation for downstream causal estimators.

### Lab 5 — Average Treatment Effect
Estimating the **Average Treatment Effect** via inverse-probability weighting, g-computation,
and the **doubly-robust (AIPW)** estimator.

### Lab 6 — Mediation Analysis
Decomposing a total effect into **natural direct and indirect effects** to understand the
mechanisms through which a treatment operates.

### Lab 8 — Assumption-Lean Inference
**Assumption-lean** estimation and inference that remains valid under weak modelling
assumptions, including supporting helper functions and an appendix.

### Lab 9 — Instrumental Variables
Identification and estimation of causal effects using **instrumental variables** when
unmeasured confounding is present.

### Lab 10 — Solving an Old Exam
Fully worked solutions to a past exam, plus extra material illustrating subtleties
(e.g. when the bootstrap fails).

---

## 🗂️ Structure

Each lab folder follows a consistent layout:

```
Lab N - Topic/
├── CML_LabN_Text_2026.pdf        # assignment / problem text
├── CML_LabN_Template_2026.R      # starter code
├── CML_LabN_Solutions_2026.R     # full solutions
└── CML_LabN_Markdown_2026.Rmd    # R Markdown write-up (+ rendered .html)
```

## 🛠️ Requirements

- **R** (with R Markdown / `knitr` to render the `.Rmd` notebooks)
- Common packages used across the labs include `SuperLearner`, `tmle`/`AIPW`-style
  estimators, and standard tidyverse tooling — install per-lab as needed.

## 📝 Notes

- Previous-year materials (`last_year/`), local data folders (`data/`), and tooling
  files are intentionally git-ignored.
