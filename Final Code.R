# Smoking Prevalence by Birth Cohort, Cigarettes per Smoker, and Lung Cancer
# Mortality, 2000-2022
#
# Statistical analysis code (R, version 4.4 or later)
#
# Input : JAMA_Replication_Dataset.xlsx, placed in the working directory.
#         Sheets used: Mortality, Control_Outcomes, CI5_Incidence.
# Output: estimates printed to the console and saved to analysis_results.csv.
# Packages: readxl, MASS. A full run takes about 1 hour on a desktop computer;
# most of the time is spent on the trend models and the sex-by-age models.
#
# Sections
#   1  Data and model set-up
#   2  Mortality trends within countries (Figure 1, eTable 7)
#   3  Primary associations (Table 2, eTables 8 and 12)
#   4  Effect modification (Table 2, eTable 9)
#   5  Birth-cohort-aligned prevalence (Table 2, Figure 2, eTables 14 and 15)
#   6  Lung cancer incidence in CI5plus registries (eTable 17)
#   7  Control causes of death (Figure 3A, eTable 18)
#   8  Histologic types (Figure 3B, eTables 19 and 20)
#   9  Sensitivity analyses (eTable 22)
#
# All models are negative binomial (NB2) regressions with a log link and the
# logarithm of population (or person-years) as an offset. Confidence intervals
# use country-clustered sandwich variances with the CR1 correction and a t
# distribution with (number of countries - 1) degrees of freedom.

library(readxl)
library(MASS)

data_file <- "JAMA_Replication_Dataset.xlsx"
ages <- c("45-49", "50-54", "55-59", "60-64", "65-69", "70-74", "75-79", "80-84")


# ---------------------------------------------------------------------------
# 1. Data and model set-up
# ---------------------------------------------------------------------------

mort_all <- as.data.frame(read_excel(data_file, sheet = "Mortality", guess_max = 30000))
controls <- as.data.frame(read_excel(data_file, sheet = "Control_Outcomes", guess_max = 30000))
ci5      <- as.data.frame(read_excel(data_file, sheet = "CI5_Incidence", guess_max = 30000))

names(mort_all)[names(mort_all) == "deaths_lung_c33_c34"] <- "deaths"
mort_all$who_icd10_list <- as.character(mort_all$who_icd10_list)
controls$who_icd10_list <- as.character(controls$who_icd10_list)

# Exposures are scaled so that 1 unit is 5 percentage points of prevalence or
# 5 cigarettes per smoker per day; centering does not affect the rate ratios.
prepare <- function(d) {
  d$prev5  <- (d$prevalence_pct - median(d$prevalence_pct)) / 5
  d$int5   <- (d$cigarettes_per_smoker_day - median(d$cigarettes_per_smoker_day)) / 5
  d$sdi_z  <- (d$sdi - mean(d$sdi)) / sd(d$sdi)
  d$year_z <- (d$outcome_year - 2011) / 10
  d$male   <- as.numeric(d$sex == "Male")
  d$age_group <- factor(d$age_group, levels = ages)
  d[order(d$iso3, d$sex, d$outcome_year, d$age_group), ]
}

d <- prepare(mort_all[mort_all$primary_sample == 1, ])
d$coh5 <- (d$cohort_aligned_prevalence_pct - median(d$cohort_aligned_prevalence_pct)) / 5

# Indicators used in the interaction models
age_no <- as.integer(d$age_group)
d$band2 <- as.numeric(age_no %in% 3:4)        # 55-64 years
d$band3 <- as.numeric(age_no %in% 5:6)        # 65-74 years
d$band4 <- as.numeric(age_no %in% 7:8)        # 75-84 years
d$age_score <- age_no - 1
d$later <- as.numeric(d$outcome_year >= 2012)
d$lower_sdi <- as.numeric(d$sdi_group == "Lower SDI")
d$stage_mid <- as.numeric(d$epidemic_stage == "Intermediate")
d$stage_adv <- as.numeric(d$epidemic_stage == "Advanced")
for (a in ages[-1]) d[[paste0("age_", a)]] <- as.numeric(d$age_group == a)
for (s in c("Female", "Male")) for (a in ages) {
  if (s == "Female" && a == ages[1]) next
  d[[paste("cell", s, a, sep = "_")]] <- as.numeric(d$sex == s & d$age_group == a)
}

# Design matrix. The fixed-effect unit is country by sex (primary analysis),
# country (prespecified design), or country by sex by age group.
design <- function(d, exposures, fe = c("country_sex", "country", "country_sex_age"),
                   interactions = NULL, extra = NULL, adjust = TRUE,
                   trends = FALSE, year_fe = FALSE) {
  fe <- match.arg(fe)
  X <- cbind("(Intercept)" = 1, as.matrix(d[, exposures, drop = FALSE]))
  for (term in interactions) {
    v <- strsplit(term, ":", fixed = TRUE)[[1]]
    X <- cbind(X, d[[v[1]]] * d[[v[2]]])
    colnames(X)[ncol(X)] <- term
  }
  if (length(extra)) X <- cbind(X, as.matrix(d[, extra, drop = FALSE]))
  if (!adjust) return(X)

  age <- model.matrix(~ age_group, d)[, -1, drop = FALSE]
  X <- cbind(X, sdi_z = d$sdi_z)
  if (fe == "country") X <- cbind(X, male = d$male)
  if (fe != "country_sex_age") X <- cbind(X, age)
  if (year_fe) {
    X <- cbind(X, model.matrix(~ factor(outcome_year), d)[, -1, drop = FALSE])
  } else {
    X <- cbind(X, year_z = d$year_z, year_z2 = d$year_z^2)
  }
  unit <- switch(fe,
                 country         = d$iso3,
                 country_sex     = paste(d$iso3, d$sex),
                 country_sex_age = paste(d$iso3, d$sex, d$age_group))
  X <- cbind(X, model.matrix(~ factor(unit))[, -1, drop = FALSE])
  if (fe != "country_sex_age") {
    male_age <- age * d$male
    colnames(male_age) <- paste0("male:", colnames(age))
    X <- cbind(X, male_age)
  }
  X <- cbind(X, "male:year_z" = d$male * d$year_z)
  if (trends) {
    slope <- model.matrix(~ factor(iso3), d)[, -1, drop = FALSE] * d$year_z
    colnames(slope) <- paste0(colnames(slope), ":year_z")
    X <- cbind(X, slope)
  }
  X
}

# Maximum likelihood fit and CR1 country-clustered variance
fit_glm <- function(X, y, offset, cluster, poisson = FALSE) {
  ctrl <- glm.control(epsilon = 1e-10, maxit = 100)
  g <- if (poisson) {
    glm(y ~ X - 1 + offset(offset), family = stats::poisson(), control = ctrl)
  } else {
    glm.nb(y ~ X - 1 + offset(offset), control = ctrl)
  }
  if (anyNA(coef(g))) stop("design matrix is not of full rank")
  beta  <- setNames(coef(g), colnames(X))
  mu    <- fitted(g)
  theta <- if (poisson) Inf else g$theta
  den   <- if (poisson) 1 else 1 + mu / theta
  bread <- chol2inv(chol(crossprod(X * (mu / den), X)))
  score <- rowsum(X * ((y - mu) / den), cluster)
  G <- nrow(score); N <- nrow(X); K <- ncol(X)
  V <- bread %*% crossprod(score) %*% bread * (G / (G - 1)) * ((N - 1) / (N - K))
  dimnames(V) <- list(colnames(X), colnames(X))
  dimnames(bread) <- dimnames(V)
  list(beta = beta, V = V, G = G, X = X, y = y, mu = mu, theta = theta,
       bread = bread, score = score, loglik = as.numeric(logLik(g)), aic = AIC(g))
}

fit_nb <- function(d, exposures, poisson = FALSE, ...) {
  X <- design(d, exposures, ...)
  fit_glm(X, d$deaths, log(d$population), d$iso3, poisson)
}

# Linear combination of coefficients on the rate ratio scale
lincom <- function(m, w) {
  est  <- sum(m$beta[names(w)] * w)
  se   <- sqrt(drop(w %*% m$V[names(w), names(w), drop = FALSE] %*% w))
  crit <- qt(0.975, m$G - 1)
  data.frame(rr = exp(est), lo = exp(est - crit * se), hi = exp(est + crit * se),
             p = 2 * pt(-abs(est / se), m$G - 1))
}
rr <- function(m, term, mult = 1) lincom(m, setNames(mult, term))

# Cluster-robust Wald test that a set of linear combinations equals 0
wald <- function(m, L) {
  if (is.character(L)) L <- diag(length(m$beta))[match(L, names(m$beta)), , drop = FALSE]
  b <- drop(L %*% m$beta)
  stat <- drop(b %*% solve(L %*% m$V %*% t(L), b)) / length(b)
  pf(stat, length(b), m$G - 1, lower.tail = FALSE)
}

# Model-standardized absolute rate difference per 100 000 person-years
# (Segi-Doll weights, ages 45-84 years), with a delta-method CI
rate_difference <- function(m, d, term, mult = 1) {
  rate  <- m$mu / d$population * 1e5
  units <- nrow(unique(d[, c("iso3", "sex", "outcome_year")]))
  w     <- d$segi_weight / units
  base  <- sum(w * rate)
  ratio <- exp(m$beta[[term]] * mult)
  grad  <- colSums(m$X * (w * rate)) * (ratio - 1)
  grad[term] <- grad[term] + base * ratio * mult
  est  <- base * (ratio - 1)
  se   <- sqrt(drop(grad %*% m$V %*% grad))
  crit <- qt(0.975, m$G - 1)
  data.frame(rd = est, lo = est - crit * se, hi = est + crit * se,
             p = 2 * pt(-abs(est / se), m$G - 1))
}

# Country-sex units without any event carry no information under
# country-by-sex fixed effects and are removed from that outcome's models.
drop_empty_units <- function(x) {
  total <- ave(x$deaths, x$iso3, x$sex, FUN = sum)
  x[total > 0, ]
}

results <- list()
result_cols <- c("estimate", "lo", "hi", "p", "p_interaction", "q_fdr", "p_test")
keep <- function(section, label, est) {
  names(est)[names(est) %in% c("rr", "apc", "rd")] <- "estimate"
  est[setdiff(result_cols, names(est))] <- NA
  est <- cbind(section = section, analysis = label, est[result_cols], row.names = NULL)
  results[[length(results) + 1]] <<- est
  invisible(est)
}
show <- function(section) {
  tab <- do.call(rbind, Filter(function(r) r$section[1] == section, results))
  tab <- tab[, c(TRUE, TRUE, colSums(!is.na(tab[result_cols])) > 0)]
  num <- vapply(tab, is.numeric, logical(1))
  tab[num] <- lapply(tab[num], signif, 4)
  cat("\n", section, "\n", sep = "")
  print(tab[, -1], row.names = FALSE)
}

cat("Countries:", length(unique(d$iso3)),
    "| country-sex units:", nrow(unique(d[, c("iso3", "sex")])),
    "| observations:", nrow(d),
    "| deaths:", sum(d$deaths),
    "| person-years:", sum(d$population), "\n")


# ---------------------------------------------------------------------------
# 2. Mortality trends within countries
# ---------------------------------------------------------------------------

# Country and age-group fixed effects, fitted separately for women and men.
# mode = "years": calendar-year effects relative to 2000; "piecewise": linear
# trend with a change of slope at the knot; "linear": single linear trend.
trend_fit <- function(x, mode, knot = 2015) {
  total <- tapply(x$deaths, x$iso3, sum)
  x <- x[x$iso3 %in% names(total)[total > 0], ]
  X <- model.matrix(~ 0 + factor(iso3), x)
  age <- droplevels(x$age_group)
  if (nlevels(age) > 1) X <- cbind(X, model.matrix(~ age)[, -1, drop = FALSE])
  yr <- x$outcome_year
  if (mode == "years") {
    yy <- model.matrix(~ factor(yr))[, -1, drop = FALSE]
    colnames(yy) <- paste0("y", sort(unique(yr))[-1])
    X <- cbind(X, yy)
  } else if (mode == "piecewise") {
    X <- cbind(X, year = yr - knot, after = pmax(yr - knot, 0))
  } else {
    X <- cbind(X, year = yr - 2011)
  }
  m <- fit_glm(X, x$deaths, log(x$population), x$iso3)
  m$countries <- length(unique(x$iso3))
  m$deaths <- sum(x$deaths)
  m
}
apc <- function(est) data.frame(apc = 100 * (est$rr - 1), lo = 100 * (est$lo - 1),
                                hi = 100 * (est$hi - 1), p = est$p)

sec <- "2. Mortality trends"
for (s in c("Female", "Male")) {
  x <- d[d$sex == s, ]
  lab <- if (s == "Female") "Women" else "Men"

  m <- trend_fit(x, "years")
  keep(sec, paste(lab, "rate ratio, 2022 vs 2000"), rr(m, "y2022"))

  groups <- c(list("All countries" = x),
              split(x, x$sdi_group),
              split(x, x$age_group))
  for (g in names(groups)) {
    mp <- trend_fit(groups[[g]], "piecewise")
    ml <- trend_fit(groups[[g]], "linear")
    keep(sec, paste(lab, g, "APC 2000-2015"), apc(rr(mp, "year")))
    keep(sec, paste(lab, g, "APC 2015-2022"), apc(lincom(mp, c(year = 1, after = 1))))
    keep(sec, paste(lab, g, "APC 2000-2022"), apc(rr(ml, "year")))
    keep(sec, paste(lab, g, "change in slope at 2015"),
         data.frame(apc = NA, lo = NA, hi = NA, p = rr(mp, "after")$p))
  }

  # Profile likelihood for the year of the change in slope
  ll <- sapply(2005:2019, function(k) trend_fit(x, "piecewise", knot = k)$loglik)
  cat(lab, ": best-fitting year of change =", (2005:2019)[which.max(ll)], "\n")
}
show(sec)


# ---------------------------------------------------------------------------
# 3. Primary associations
# ---------------------------------------------------------------------------

expo <- c("prev5", "int5")
m_prev <- fit_nb(d, "prev5")
m_main <- fit_nb(d, expo)
iqr_prev <- IQR(d$prevalence_pct)
iqr_int  <- IQR(d$cigarettes_per_smoker_day)

sec <- "3. Primary associations (rate ratios)"
keep(sec, "Prevalence, unadjusted, per 5 points", rr(fit_nb(d, "prev5", adjust = FALSE), "prev5"))
keep(sec, "Intensity, unadjusted, per 5 cigarettes", rr(fit_nb(d, "int5", adjust = FALSE), "int5"))
keep(sec, "Prevalence, adjusted, per 5 points", rr(m_main, "prev5"))
keep(sec, "Intensity, adjusted, per 5 cigarettes", rr(m_main, "int5"))
keep(sec, sprintf("Prevalence, adjusted, per IQR (%.1f points)", iqr_prev), rr(m_main, "prev5", iqr_prev / 5))
keep(sec, sprintf("Intensity, adjusted, per IQR (%.1f cigarettes)", iqr_int), rr(m_main, "int5", iqr_int / 5))
show(sec)

sec <- "3. Absolute rate differences per 100 000 person-years"
ard <- rbind(rate_difference(m_main, d, "prev5"), rate_difference(m_main, d, "int5"))
names(ard)[1] <- "rr"
keep(sec, "Prevalence, per 5 points", ard[1, ])
keep(sec, "Intensity, per 5 cigarettes", ard[2, ])
show(sec)

lr <- 2 * (m_main$loglik - m_prev$loglik)
cat("\nAdding intensity to the prevalence model: likelihood ratio statistic =",
    round(lr, 3), ", P =", signif(pchisq(lr, 1, lower.tail = FALSE), 3),
    "; AIC", round(m_prev$aic, 1), "->", round(m_main$aic, 1), "\n")


# ---------------------------------------------------------------------------
# 4. Effect modification
# ---------------------------------------------------------------------------

# One model per exposure and modifier. `levels` maps indicator columns to
# subgroup labels; the reference subgroup has all indicators equal to 0.
subgroups <- function(d, term, reference, levels, exposures = expo, extra = NULL) {
  inter <- paste(term, names(levels), sep = ":")
  m <- fit_nb(d, exposures, interactions = inter, extra = extra)
  out <- rr(m, term)
  for (i in inter) out <- rbind(out, lincom(m, setNames(c(1, 1), c(term, i))))
  cbind(subgroup = c(reference, unname(levels)), out, p_interaction = wald(m, inter))
}

modifiers <- list(
  "Sex"             = list("Women", c(male = "Men"), NULL),
  "Age group"       = list("45-54 years", c(band2 = "55-64 years", band3 = "65-74 years", band4 = "75-84 years"), NULL),
  "Calendar period" = list("2000-2011", c(later = "2012-2022"), "later"),
  "SDI group"       = list("High SDI", c(lower_sdi = "Lower SDI"), NULL),
  "Epidemic stage"  = list("Early", c(stage_mid = "Intermediate", stage_adv = "Advanced"), NULL))

sec <- "4. Effect modification"
family <- NULL
for (term in expo) for (mod in names(modifiers)) {
  s <- modifiers[[mod]]
  tab <- subgroups(d, term, s[[1]], s[[2]], extra = s[[3]])
  family <- rbind(family, cbind(exposure = term, modifier = mod, tab))
}
tests <- unique(family[, c("exposure", "modifier", "p_interaction")])
tests$q_fdr <- p.adjust(tests$p_interaction, method = "BH")   # 10 tests
family <- merge(family, tests, sort = FALSE)
for (i in seq_len(nrow(family))) {
  keep(sec, paste(family$exposure[i], family$modifier[i], family$subgroup[i], sep = " | "),
       family[i, c("rr", "lo", "hi", "p", "p_interaction", "q_fdr")])
}
show(sec)

sec <- "4. Linear trend across 5-year age groups (ratio of rate ratios)"
for (term in expo) {
  it <- paste0(term, ":age_score")
  keep(sec, term, rr(fit_nb(d, expo, interactions = it), it))
}
show(sec)


# ---------------------------------------------------------------------------
# 5. Birth-cohort-aligned prevalence
# ---------------------------------------------------------------------------

coh <- c("coh5", "int5")
m_coh   <- fit_nb(d, coh)
m_joint <- fit_nb(d, c("prev5", "coh5", "int5"))

sec <- "5. Cohort-aligned vs national prevalence"
keep(sec, "National prevalence (primary model)", rr(m_main, "prev5"))
keep(sec, "Cohort-aligned prevalence (replaces national)", rr(m_coh, "coh5"))
keep(sec, "Intensity in the cohort-aligned model", rr(m_coh, "int5"))
keep(sec, "Joint model: national prevalence", rr(m_joint, "prev5"))
keep(sec, "Joint model: cohort-aligned prevalence", rr(m_joint, "coh5"))
show(sec)
cat("AIC, national:", round(m_main$aic, 1), "| cohort-aligned:", round(m_coh$aic, 1),
    "| difference:", round(m_main$aic - m_coh$aic, 1),
    "\nCorrelation between the 2 prevalence measures:",
    round(cor(d$prevalence_pct, d$cohort_aligned_prevalence_pct), 3), "\n")

# Rate ratios by 5-year age group, women and men combined (Figure 2A)
sec <- "5. Rate ratios by age group, women and men combined"
age_levels <- setNames(ages[-1], paste0("age_", ages[-1]))
for (v in list(c("prev5", "National"), c("coh5", "Cohort-aligned"), c("int5", "Intensity"))) {
  tab <- subgroups(d, v[1], ages[1], age_levels,
                   exposures = if (v[1] == "coh5") coh else expo)
  for (i in seq_len(nrow(tab))) {
    keep(sec, paste(v[2], tab$subgroup[i]), tab[i, c("rr", "lo", "hi", "p", "p_interaction")])
  }
}
show(sec)

# Rate ratios by sex and age group (Figure 2C, eTable 15)
sec <- "5. Rate ratios by sex and age group"
cells <- grep("^cell_", names(d), value = TRUE)
for (v in list(c("prev5", "National"), c("coh5", "Cohort-aligned"))) {
  term <- v[1]
  ms <- fit_nb(d, c(term, "int5"), interactions = paste0(term, ":male"))
  keep(sec, paste(v[2], "Women, all ages"), cbind(rr(ms, term), p_test = NA))
  keep(sec, paste(v[2], "Men, all ages"),
       cbind(lincom(ms, setNames(c(1, 1), c(term, paste0(term, ":male")))),
             p_test = rr(ms, paste0(term, ":male"))$p))

  mc <- fit_nb(d, c(term, "int5"), interactions = paste(term, cells, sep = ":"))
  pick <- function(s, a) {
    w <- setNames(numeric(length(mc$beta)), names(mc$beta))
    w[term] <- 1
    cell <- paste0(term, ":cell_", s, "_", a)
    if (cell %in% names(w)) w[cell] <- 1
    w
  }
  for (s in c("Female", "Male")) for (a in ages) {
    w <- pick(s, a)
    keep(sec, paste(v[2], if (s == "Female") "Women" else "Men", a),
         cbind(lincom(mc, w[w != 0]), p_test = NA))
  }
  across <- function(s) t(sapply(ages[-1], function(a) pick(s, a) - pick(s, ages[1])))
  keep(sec, paste(v[2], "heterogeneity across age groups, women"),
       data.frame(rr = NA, lo = NA, hi = NA, p = NA, p_test = wald(mc, across("Female"))))
  keep(sec, paste(v[2], "heterogeneity across age groups, men"),
       data.frame(rr = NA, lo = NA, hi = NA, p = NA, p_test = wald(mc, across("Male"))))
  keep(sec, paste(v[2], "difference in age pattern, women vs men"),
       data.frame(rr = NA, lo = NA, hi = NA, p = NA,
                  p_test = wald(mc, across("Male") - across("Female"))))
}
show(sec)

# Absolute rate differences per 5 points of cohort-aligned prevalence by sex
# and age group: mean over the observations of each group of the fitted rate
# multiplied by (rate ratio - 1)
sec <- "5. Absolute rate differences by sex and age group (cohort-aligned)"
mc <- fit_nb(d, coh, interactions = paste("coh5", cells, sep = ":"))
fitted_rate <- mc$mu / d$population * 1e5
for (s in c("Female", "Male")) for (a in ages) {
  rows <- d$sex == s & d$age_group == a
  cell <- paste0("coh5:cell_", s, "_", a)
  terms <- intersect(c("coh5", cell), names(mc$beta))
  ratio <- exp(sum(mc$beta[terms]))
  base  <- mean(fitted_rate[rows])
  grad  <- colSums(mc$X[rows, , drop = FALSE] * fitted_rate[rows]) / sum(rows) * (ratio - 1)
  grad[terms] <- grad[terms] + base * ratio
  est <- base * (ratio - 1)
  se  <- sqrt(drop(grad %*% mc$V %*% grad))
  crit <- qt(0.975, mc$G - 1)
  keep(sec, paste(if (s == "Female") "Women" else "Men", a),
       data.frame(rr = est, lo = est - crit * se, hi = est + crit * se,
                  p = 2 * pt(-abs(est / se), mc$G - 1)))
}
show(sec)


# ---------------------------------------------------------------------------
# 6. Lung cancer incidence in CI5plus registries
# ---------------------------------------------------------------------------

# One registry series per country; the model covariates in this sheet are on
# the same scales as in the mortality analysis.
ci5$prev5 <- ci5$prevalence5
ci5$int5  <- ci5$intensity5
ci5$coh5  <- ci5$cohort5
ci5$male  <- as.numeric(ci5$sex == "Male")
ci5$outcome_year <- ci5$year
ci5$population   <- ci5$person_years
ci5$age_group <- factor(ci5$age_group, levels = ages)
ci5 <- ci5[order(ci5$iso3, ci5$sex, ci5$year, ci5$age_group), ]

incidence <- function(columns) {
  x <- ci5
  x$deaths <- rowSums(x[, columns, drop = FALSE])
  drop_empty_units(x)
}

sec <- "6. Lung cancer incidence (C33-C34)"
x <- incidence("cases_lung_c33_c34")
for (sub in list(list("43 countries", x),
                 list("countries also in the mortality sample", x[x$iso3 %in% d$iso3, ]))) {
  mn <- fit_nb(sub[[2]], expo)
  mc <- fit_nb(sub[[2]], coh)
  keep(sec, paste("National prevalence,", sub[[1]]), rr(mn, "prev5"))
  keep(sec, paste("Intensity (national model),", sub[[1]]), rr(mn, "int5"))
  keep(sec, paste("Cohort-aligned prevalence,", sub[[1]]), rr(mc, "coh5"))
  keep(sec, paste("Intensity (cohort-aligned model),", sub[[1]]), rr(mc, "int5"))
}
show(sec)


# ---------------------------------------------------------------------------
# 7. Control causes of death
# ---------------------------------------------------------------------------

# Countries and years reported on the detailed WHO cause lists
keys <- c("iso3", "sex", "outcome_year", "age_group", "who_icd10_list")
detail <- merge(d, controls, by = keys)

causes <- c("Lung cancer, detailed-list subset" = "deaths",
            "Larynx cancer (C32)"               = "deaths_larynx_c32",
            "COPD (J40-J44)"                    = "deaths_copd_j40_j44",
            "Bladder cancer (C67)"              = "deaths_bladder_c67",
            "Brain cancer (C71)"                = "deaths_brain_c71")

sec <- "7. Control causes of death"
for (cause in names(causes)) {
  x <- detail[!is.na(detail[[causes[[cause]]]]), ]
  x$deaths <- x[[causes[[cause]]]]
  x <- drop_empty_units(prepare(x))
  mn <- fit_nb(x, expo)
  mc <- fit_nb(x, coh)
  keep(sec, paste(cause, "| national prevalence"), rr(mn, "prev5"))
  keep(sec, paste(cause, "| cohort-aligned prevalence"), rr(mc, "coh5"))
  keep(sec, paste(cause, "| intensity (cohort-aligned model)"), rr(mc, "int5"))
}
show(sec)


# ---------------------------------------------------------------------------
# 8. Histologic types
# ---------------------------------------------------------------------------

types <- list(
  "Lung (C34), all types"                 = "cases_lung_c34_all",
  "Squamous cell carcinoma"               = "cases_squamous_cell",
  "Small cell carcinoma"                  = "cases_small_cell",
  "Squamous cell or small cell carcinoma" = c("cases_squamous_cell", "cases_small_cell"),
  "Large cell carcinoma"                  = "cases_large_cell",
  "Adenocarcinoma"                        = "cases_adenocarcinoma",
  "Larynx cancer (C32)"                   = "cases_larynx_c32",
  "Bladder cancer (C67)"                  = "cases_bladder_c67",
  "Brain and CNS cancer (C71-C72)"        = "cases_brain_cns_c71_c72")

sec <- "8. Incidence by histologic type and control cancers"
fits_nat <- fits_coh <- list()
for (type in names(types)) {
  x <- incidence(types[[type]])
  fits_nat[[type]] <- fit_nb(x, expo)
  fits_coh[[type]] <- fit_nb(x, coh)
  keep(sec, paste(type, "| national prevalence"), rr(fits_nat[[type]], "prev5"))
  keep(sec, paste(type, "| intensity (national model)"), rr(fits_nat[[type]], "int5"))
  keep(sec, paste(type, "| cohort-aligned prevalence"), rr(fits_coh[[type]], "coh5"))
  keep(sec, paste(type, "| intensity (cohort-aligned model)"), rr(fits_coh[[type]], "int5"))
}
show(sec)

# Ratio of the rate ratios from 2 outcome models fitted to the same registry
# series. The variance of the difference in log rate ratios uses the joint
# country-clustered covariance of the 2 sets of score contributions.
ratio_of_rr <- function(m1, m2, term) {
  common <- intersect(rownames(m1$score), rownames(m2$score))
  s1 <- m1$score[common, , drop = FALSE] %*% m1$bread
  s2 <- m2$score[common, , drop = FALSE] %*% m2$bread
  g  <- length(common)
  v  <- (sum(s1[, term]^2) + sum(s2[, term]^2) - 2 * sum(s1[, term] * s2[, term])) * g / (g - 1)
  est <- m1$beta[[term]] - m2$beta[[term]]
  crit <- qt(0.975, g - 1)
  data.frame(rr = exp(est), lo = exp(est - crit * sqrt(v)), hi = exp(est + crit * sqrt(v)),
             p = 2 * pt(-abs(est / sqrt(v)), g - 1))
}

sec <- "8. Squamous cell or small cell carcinoma vs adenocarcinoma (ratio of rate ratios)"
sq <- "Squamous cell or small cell carcinoma"
keep(sec, "National prevalence", ratio_of_rr(fits_nat[[sq]], fits_nat[["Adenocarcinoma"]], "prev5"))
keep(sec, "Intensity (national model)", ratio_of_rr(fits_nat[[sq]], fits_nat[["Adenocarcinoma"]], "int5"))
keep(sec, "Cohort-aligned prevalence", ratio_of_rr(fits_coh[[sq]], fits_coh[["Adenocarcinoma"]], "coh5"))
keep(sec, "Intensity (cohort-aligned model)", ratio_of_rr(fits_coh[[sq]], fits_coh[["Adenocarcinoma"]], "int5"))
show(sec)


# ---------------------------------------------------------------------------
# 9. Sensitivity analyses
# ---------------------------------------------------------------------------

sec <- "9. Sensitivity analyses"
sens <- function(label, data, ...) {
  m <- fit_nb(data, expo, ...)
  keep(sec, paste(label, "| prevalence"), rr(m, "prev5"))
  keep(sec, paste(label, "| intensity"), rr(m, "int5"))
}
with_deaths <- ave(d$deaths, d$iso3, d$sex, d$age_group, FUN = sum) > 0

sens("Primary model", d)
sens("Prespecified design: country fixed effects", d, fe = "country")
sens("Country-by-sex-by-age fixed effects", d[with_deaths, ], fe = "country_sex_age")
sens("Country-specific linear trends", d, trends = TRUE)
sens("Calendar-year fixed effects", d, year_fe = TRUE)
sens("Poisson pseudo-maximum likelihood", d, poisson = TRUE)
sens("Excluding 2020-2022", d[d$outcome_year <= 2019, ])
sens("ICD-10 list 104 only", d[d$who_icd10_list == "104", ])
sens("WHO population denominators only", d[d$population_source == "WHO", ])
sens("Prespecified 40-country sample", prepare(mort_all[mort_all$original_40_sample == 1, ]))
sens("All 78 eligible countries", prepare(mort_all))
show(sec)

write.csv(do.call(rbind, results), "analysis_results.csv", row.names = FALSE)
