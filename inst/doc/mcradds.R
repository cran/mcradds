## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----eval = FALSE-------------------------------------------------------------
#  browseVignettes(package = "mcradds")

## -----------------------------------------------------------------------------
library(mcradds)

## -----------------------------------------------------------------------------
data("qualData")
data("platelet")
data(creatinine, package = "mcr")
data("calcium")
data("ldlroc")
data("PDL1RP")
data("glucose")
data("adsl_sub")

## -----------------------------------------------------------------------------
size_one_prop(p1 = 0.9, p0 = 0.85, alpha = 0.05, power = 0.8)

## -----------------------------------------------------------------------------
size_ci_one_prop(p = 0.85, lr = 0.8, alpha = 0.05, method = "wilson")

## -----------------------------------------------------------------------------
size_ci_one_prop(p = 0.85, lr = 0.8, alpha = 0.05, method = "simple-asymptotic")

## -----------------------------------------------------------------------------
size_corr(r1 = 0.95, r0 = 0.9, alpha = 0.025, power = 0.8, alternative = "greater")

## -----------------------------------------------------------------------------
size_ci_corr(r = 0.9, lr = 0.85, alpha = 0.025, alternative = "greater")

## -----------------------------------------------------------------------------
adsl_sub %>%
  descfreq(
    var = "AGEGR1",
    bygroup = "TRTP",
    format = "xx (xx.x%)"
  )

## -----------------------------------------------------------------------------
adsl_sub %>%
  descfreq(
    var = c("AGEGR1", "SEX", "RACE"),
    bygroup = "TRTP",
    format = "xx (xx.x%)",
    addtot = TRUE,
    na_str = "0"
  )

## -----------------------------------------------------------------------------
adsl_sub %>%
  descvar(
    var = "AGE",
    bygroup = "TRTP"
  )

## -----------------------------------------------------------------------------
adsl_sub %>%
  descvar(
    var = c("AGE", "BMIBL", "HEIGHTBL"),
    bygroup = "TRTP",
    stats = c("N", "MEANSD", "MEDIAN", "RANGE", "IQR"),
    autodecimal = TRUE,
    addtot = TRUE
  )

## -----------------------------------------------------------------------------
head(qualData)

## -----------------------------------------------------------------------------
tb <- qualData %>%
  diagTab(
    formula = ~ CandidateN + ComparativeN,
    levels = c(1, 0)
  )
tb

## -----------------------------------------------------------------------------
dummy <- data.frame(
  id = c("1001", "1001", "1002", "1002", "1003", "1003"),
  value = c(1, 0, 0, 0, 1, 1),
  type = c("Test", "Ref", "Test", "Ref", "Test", "Ref")
) %>%
  diagTab(
    formula = type ~ value,
    bysort = "id",
    dimname = c("Test", "Ref"),
    levels = c(1, 0)
  )
dummy

## -----------------------------------------------------------------------------
# Default method is Wilson score, and digit is 4.
tb %>% getAccuracy(ref = "r")

# Alter the number of digit to 2.
tb %>% getAccuracy(ref = "r", digit = 2)

# Alter the number of digit to 2.
tb %>% getAccuracy(ref = "r", r_ci = "clopper-pearson")

## -----------------------------------------------------------------------------
# When the reference is a comparative assay, not gold standard.
tb %>% getAccuracy(ref = "nr", nr_ci = "wilson")

## -----------------------------------------------------------------------------
# Deming regression
fit <- mcreg(
  x = platelet$Comparative, y = platelet$Candidate,
  error.ratio = 1, method.reg = "Deming", method.ci = "jackknife"
)
printSummary(fit)
getCoefficients(fit)

## -----------------------------------------------------------------------------
# absolute bias.
calcBias(fit, x.levels = c(30))

# proportional bias.
calcBias(fit, x.levels = c(30), type = "proportional")

## -----------------------------------------------------------------------------
# Default difference type
blandAltman(
  x = platelet$Comparative, y = platelet$Candidate,
  type1 = 3, type2 = 5
)

# Change relative different type to 4.
blandAltman(
  x = platelet$Comparative, y = platelet$Candidate,
  type1 = 3, type2 = 4
)

## -----------------------------------------------------------------------------
# ESD approach
ba <- blandAltman(x = platelet$Comparative, y = platelet$Candidate)
out <- getOutlier(ba, method = "ESD", difference = "rel")
out$stat
out$outmat

# 4E approach
ba2 <- blandAltman(x = creatinine$serum.crea, y = creatinine$plasma.crea)
out2 <- getOutlier(ba2, method = "4E")
out2$stat
out2$outmat

## -----------------------------------------------------------------------------
x <- c(44.4, 45.9, 41.9, 53.3, 44.7, 44.1, 50.7, 45.2, 60.1)
y <- c(2.6, 3.1, 2.5, 5.0, 3.6, 4.0, 5.2, 2.8, 3.8)
pearsonTest(x, y, h0 = 0.5, alternative = "greater")

## -----------------------------------------------------------------------------
x <- c(44.4, 45.9, 41.9, 53.3, 44.7, 44.1, 50.7, 45.2, 60.1)
y <- c(2.6, 3.1, 2.5, 5.0, 3.6, 4.0, 5.2, 2.8, 3.8)
spearmanTest(x, y, h0 = 0.5, alternative = "greater")

## -----------------------------------------------------------------------------
refInterval(x = calcium$Value, RI_method = "parametric", CI_method = "parametric")

## -----------------------------------------------------------------------------
refInterval(x = calcium$Value, RI_method = "nonparametric", CI_method = "nonparametric")

## -----------------------------------------------------------------------------
refInterval(x = calcium$Value, RI_method = "robust", CI_method = "boot")

## -----------------------------------------------------------------------------
# H0 : Difference between areas = 0:
aucTest(x = ldlroc$LDL, y = ldlroc$OxLDL, response = ldlroc$Diagnosis)

## -----------------------------------------------------------------------------
# H0 : Superiority margin <= 0.1:
aucTest(
  x = ldlroc$LDL, y = ldlroc$OxLDL, response = ldlroc$Diagnosis,
  method = "superiority", h0 = 0.1
)

## -----------------------------------------------------------------------------
# H0 : Non-inferiority margin <= -0.1:
aucTest(
  x = ldlroc$LDL, y = ldlroc$OxLDL, response = ldlroc$Diagnosis,
  method = "non-inferiority", h0 = -0.1
)

## -----------------------------------------------------------------------------
reader <- PDL1RP$btw_reader
tb1 <- reader %>%
  diagTab(
    formula = Reader ~ Value,
    bysort = "Sample",
    levels = c("Positive", "Negative"),
    rep = TRUE,
    across = "Site"
  )
getAccuracy(tb1, ref = "bnr", rng.seed = 12306)

## -----------------------------------------------------------------------------
read <- PDL1RP$wtn_reader
tb2 <- read %>%
  diagTab(
    formula = Order ~ Value,
    bysort = "Sample",
    levels = c("Positive", "Negative"),
    rep = TRUE,
    across = "Sample"
  )
getAccuracy(tb2, ref = "bnr", rng.seed = 12306)

## -----------------------------------------------------------------------------
site <- PDL1RP$btw_site
tb3 <- site %>%
  diagTab(
    formula = Site ~ Value,
    bysort = "Sample",
    levels = c("Positive", "Negative"),
    rep = TRUE,
    across = "Sample"
  )
getAccuracy(tb2, ref = "bnr", rng.seed = 12306)

## -----------------------------------------------------------------------------
fit <- anovaVCA(value ~ day / run, glucose)
VCAinference(fit)

## -----------------------------------------------------------------------------
object <- blandAltman(x = platelet$Comparative, y = platelet$Candidate)

# Absolute difference plot
autoplot(object, type = "absolute")

# Relative difference plot
autoplot(object, type = "relative")

## -----------------------------------------------------------------------------
autoplot(
  object,
  type = "absolute",
  jitter = TRUE,
  fill = "lightblue",
  color = "grey",
  size = 2,
  ref.line.params = list(col = "grey"),
  loa.line.params = list(col = "grey"),
  label.digits = 2,
  label.params = list(col = "grey", size = 3, fontface = "italic"),
  x.nbreak = 6,
  main.title = "Bland-Altman Plot",
  x.title = "Mean of Test and Reference Methods",
  y.title = "Reference - Test"
)

## -----------------------------------------------------------------------------
fit <- mcreg(
  x = platelet$Comparative, y = platelet$Candidate,
  method.reg = "PaBa", method.ci = "bootstrap"
)
autoplot(fit)

## -----------------------------------------------------------------------------
autoplot(
  fit,
  identity.params = list(col = "blue", linetype = "solid"),
  reg.params = list(col = "red", linetype = "solid"),
  equal.axis = TRUE,
  legend.title = FALSE,
  legend.digits = 3,
  x.title = "Reference",
  y.title = "Test"
)

## ----sessionInfo, echo=FALSE--------------------------------------------------
sessionInfo()

