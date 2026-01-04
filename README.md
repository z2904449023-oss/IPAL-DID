# Coding Samples – Stata Empirical Work  
Author: Xinyi Zhou

This repository contains three Stata projects that illustrate my experience with applied econometrics and reproducible workflows.  
The focus is on clear model specification, transparent coding, and robustness checks.  
Raw data files are **not** included due to licensing and confidentiality, but each project is fully runnable once the corresponding `.dta` file is placed in the indicated folder.

---

## 1. CFPS Smoking Behaviour and Heckman Selection (`cfps_smoking/`)

**File:** `cfps_smoking.do`  
**Data:** `Panel Data(CFPS10-16).dta`  
**Dataset:** CFPS 2010–2016 individual-level panel (China Family Panel Studies), restricted-access micro data.

**Goal**

Study the determinants of smoking behaviour and smoking intensity, and illustrate how sample selection bias can be addressed using a Heckman selection model.

**Methods**

- Linear probability model (OLS) as a baseline for the binary outcome `smoke`.
- Logit and probit models for the smoking decision with:
  - Robust standard errors
  - Marginal effects via `margins`
- Predicted probabilities from OLS/logit/probit for comparison.
- Heckman selection model for the intensive margin:
  - Outcome: number of cigarettes per day (`num_smoke`)
  - Selection equation: smoking decision (`smoke`)
  - Comparison of OLS on smokers only vs. Heckman

**How to run**

1. Place `Panel Data(CFPS10-16).dta` in the same directory as the `.do` file.
2. Open `cfps_smoking.do` in Stata.
3. Edit the `cd "YOUR_PATH_TO_CFPS_DATA"` line to your local path.
4. Run the entire script.
5. The script installs required user-written commands (e.g. `outreg2`, `asdoc`) if needed and produces regression tables and summary statistics.

---

## 2. Class Size and Achievement RDD (`class_size_rd/`)

**File:** `class_size_rd.do`  
**Data:** `lec4_grade.dta`  
**Dataset:** Angrist & Lavy–style school-level data on class size and test scores.

**Goal**

Estimate the causal effect of class size on student achievement using a regression discontinuity design at the class-size cutoff (40 students), and perform standard RD diagnostics and robustness checks.

**Methods**

- Baseline OLS regressions of average math scores (`avgmath`) on class size (`classize`) and controls (`disadv`, `enrollment`, `esquare`).
- Manual local polynomial RD around the cutoff (enrollment between 35–45):
  - Left/right local regressions and calculation of the intercept jump.
- `rdrobust` local linear RD:
  - Triangular kernel, specified bandwidth (`h=5`, `h=6`)
  - Covariate adjustment with `disadv`
- `rdplot` visualisation of the RD:
  - Different kernel choices (triangular vs. uniform)
  - Different polynomial orders (`p=1`, `p=2`)
- Density test for manipulation of the running variable using `DCdensity`.
- Placebo tests:
  - Placebo outcome (`avgverb`)
  - Placebo cutoff (c = 30)

**How to run**

1. Place `lec4_grade.dta` in the same directory as the `.do` file.
2. Open `class_size_rd.do` in Stata.
3. Edit the `cd "YOUR_PATH_TO_GRADE_DATA"` line to your local path.
4. Run the script. It will:
   - Install `rdrobust` and `DCdensity` if missing
   - Produce RD estimates and RD plots
   - Run placebo and robustness checks

---

## 3. Belt and Road Initiative (BRI) and Trade – DID (`bri_trade_did/`)

**File:** `bri_trade_did.do`  
**Data:** `bri_trade_panel.dta`  
**Dataset:** Province-level panel for China (31 provinces, 1978–2018), constructed from National Bureau of Statistics yearbooks (not included here).

**Goal**

Evaluate the impact of the Belt and Road Initiative (BRI) on provincial trade volumes using a difference-in-differences (DID) framework with fixed effects, event-study dynamics, and placebo tests.

**Variables used (in the data)**

- `province` – province identifier  
- `year` – year (1978–2018)  
- `trade` – total imports + exports  
- `primary` – primary industry value added  
- `secondary` – secondary industry value added  
- `tertiary` – tertiary industry value added  
- `treat` – indicator for BRI-affected provinces  

The `.do` file generates:

- `post` – post-BRI period indicator (year ≥ 2013)  
- `did` – DID treatment (`treat × post`)  
- Event-time dummies (`pre2`–`pre7`, `current`, `post1`–`post5`)  

**Methods**

- Two-way fixed-effects DID with clustered standard errors:
  - `trade` on `did` and sectoral controls (`primary`, `secondary`, `tertiary`)
  - Province and year fixed effects (`absorb(province year)`)
- Event-study (dynamic DID):
  - Leads and lags around the policy year (2013)
  - Coefficients interpreted as trade changes relative to the year just before BRI
  - Used to check parallel trend assumptions and dynamic policy effects
- In-time placebo:
  - Shift the policy year to a fake year (e.g. 2008)
  - Re-estimate DID to test for spurious effects
- In-space placebo:
  - Randomly reassign treatment status across provinces within each year
  - Re-estimate DID with `did_fake` to evaluate how unusual the true effect is under random assignment

**How to run**

1. Place `bri_trade_panel.dta` in the same directory as the `.do` file.
2. Open `bri_trade_did.do` in Stata.
3. Edit the `cd "YOUR_PATH_TO_BRI_DATA"` line to your local path.
4. Run the script. It will:
   - Perform the baseline DID regression
   - Estimate dynamic effects via event-study specification
   - Run in-time and in-space placebo checks

---

## Notes on data and reproducibility

- Raw `.dta` files are not included for confidentiality and licensing reasons.
- Each `.do` file is written so that, once the corresponding data file is available and placed in the right folder, the analysis can be fully reproduced.
- All non-built-in Stata commands used in these scripts (`reghdfe`, `outreg2`, `asdoc`, `rdrobust`, `DCdensity`) are installed automatically if they are not already available.

If you have any questions about the code or the underlying data construction, please feel free to contact me.
