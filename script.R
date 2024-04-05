pacman::p_load(missForest, lubridate, gdata, survminer, gtsummary, WeightIt, cobalt, smd, survival, patchwork, ggsci)
library(tidyverse)
set.cobalt.options(binary = "std")

# define ggplot2 theme...
theme_juog <- function(base_size = 10, 
                       dark_text = "#1A242F") {
  mid_text <-  monochromeR::generate_palette(
    dark_text, "go_lighter",
    n_colors = 5)[2]
  
  light_text <-  monochromeR::generate_palette(
    dark_text, "go_lighter",
    n_colors = 5)[3]
  
  theme_grey(base_size = base_size) +
    theme(text = element_text(color = mid_text, lineheight = 1.1),
          plot.title = element_text(color = dark_text, size = rel(1.2)),
          plot.subtitle = element_text(size = rel(1.1)),
          axis.text.y = element_text(color = light_text, size = rel(1)),
          axis.title.y = element_text(size = rel(1)),
          axis.text.x = element_text(color = mid_text, size = rel(1)),
          axis.title.x = element_text(size = rel(1)),
          legend.position = "top",
          legend.title = element_blank(),
          panel.grid = element_line(color = "#F3F4F5"),
          plot.caption = element_text(size = rel(1)))
}

d0 <- data %>% # import data from excel
  dplyr::mutate(id = row_number()) # data = row data

dview <- Hmisc::describe(d0) # check variables

d # cleaned data

dfull <- d %>% 
  select(fu, death_all, death_cancer, nutrfs, nutr, comp_g3, dfs, rec,
         age_rc, sex, bmi, ecog_ps_cat, smoking, utuc, nmibc, bcg, multi_tur, cis_rc, hist_rc, ct, cn,
         type, lnd, num_lnd, div, rm_other, grade, lvi, nac, ac, hgb, lognlr, alb, id, stage, stage_dum, stage_bin) %>% mutate_if(is.character, as.factor) %>% as.data.frame() # select variables

# imputation based on random forest...
set.seed(3)
dimp <- missForest(dfull, maxiter = 1, verbose = TRUE) # random forest imputation
dimp <- dimp$ximp # extract the imputed set

# define formula for PS calculation
iptw_formula <- as.formula(stage ~ 
                             age_rc + sex + bmi + ecog_ps_cat + smoking + utuc + nmibc + bcg + multi_tur + cis_rc + hist_rc + ct + cn + 
                             type + lnd + num_lnd + div + rm_other + grade + lvi + nac + ac + hgb + lognlr + alb)

# IPTW estimation 
wgt_out <- weightit(
  iptw_formula,
  data = dimp, method = "glm", stabilize = TRUE, estimand = "ATE", verbose = TRUE # evaluating average treatment effect in the treated
)

balance <- bal.tab(wgt_out, abs = TRUE, un = TRUE, thresholds = c(m = 0.1), multi.summary = TRUE) # the summary SMDs...
balance # check balance between groups

# dataset plus weights and pscores
df <- tibble(dimp, weight = wgt_out[["weights"]], pscores = wgt_out[["ps"]])

# standardized mean difference plot
label <- data.frame(
  old = c("age_rc", "sex_male", "bmi", "ecog_ps_cat_2 or more", "smoking_current", "smoking_former", "smoking_never", "utuc_yes", "nmibc_yes", "bcg_yes", 
          "multi_tur_multifocal", "cis_rc_yes", "hist_rc_Urothelial carcinoma", "hist_rc_UC with variant histology", "hist_rc_Pure non-UC", 
          "ct_cTa or cT1", "ct_cT2", "ct_cT3 or cT4", "cn_cN1", "type_Open", "type_Laparoscopic", "type_Robotic", "lnd_standard", "num_lnd", "div_Ileal conduit", "div_Orthotopic neobladder", "div_Other", 
          "grade_High grade", "rm_other_positive", "lvi_yes", "nac_yes", "ac_yes", "hgb", "lognlr", "alb"),
  new = c("Age, year (continuous)", "Sex (Male vs. Female)", "BMI (continuous)", "ECOG PS (>1 vs. 0-1)", 
          "Smoking status: Never", "Smoking status: Former", "Smoking status: Current",
          "Previous UTUC (Yes vs. No)", "Previous NMIBC (Yes vs. No)", "Previous intravesical BCG (Yes vs. No)",
          "Tumor multifocality (Yes vs. No)", "Concomitant Cis (Yes vs. No)", "Histology: UC", "Histology: UC with variant", "Histology: Pure non-UC", 
          "Clinical T stage: cTa-1", "Clinical T stage: cT2", "Clinical T stage: cT3-4", "Clinical N stage (cN1 vs. cN0)", 
          "Approach: Open", "Approach: Laparoscopic", "Approach: Robotic", 
          "Extent of LND (standard vs. extended)", "Number of lymph node removed", "Urinary diversion: Ileal conduit", "Urinary diversion: Orthotopic neobladder", "Urinary diversion: Other",
          "Tumor grade (High vs. Low)", "Soft-tissue resective margin (Yes vs. No)", "Lymphovascular invasion (Yes vs. No)",
          "Neoadjuvant chemotherapy (Yes vs. No)", "Adjuvant chemotherapy (Yes vs. No)", 
          "Baseline hemoglobin (continuous)", "Baseline NLR (log-transformed)", "Baseline albumin (continuous)")
)

# SMD plot for 3-group comparison
smdplot_three <- love.plot(wgt_out, drop.distance = TRUE, abs = FALSE, which.treat = .all,
                           colors = pal_lancet("lanonc", alpha = 1)(2), shapes = c("diamond filled", "circle filled"), sample.names = c("Unweighted population", "Weighted population"), 
                           thresholds = c(m = 0.1), size = 3.5, themes = theme_grey(), var.names = label, binary = "std") + 
  theme(plot.title = element_blank(),
        plot.subtitle = element_blank(),
        legend.title = element_blank(),
        legend.position = "top"
  ) +
  labs(
    x = "Standardized mean differences"
  )

# IPTW distributions
wdist <- ggplot(data = df,
                aes(x = weight, fill = stage)) + 
  geom_density(alpha = 0.4) + 
  theme_juog() + 
  labs(
    x = "Inverse probability of treatment weight",
    y = "Kernel density"
  ) + 
  scale_x_continuous(breaks = seq(from = 0, to = 6, by = 1)) + 
  scale_fill_lancet()

# table1 for all patients (crude data)
reset_gtsummary_theme()
tbl_theme <- list(
  "tbl_summary-str:missing_stat" = "{N_miss} ({p_miss})" # display the number of missing values and its percentage
)
set_gtsummary_theme(tbl_theme)

tbl1_crude <- dfull %>%
  dplyr::select(age_rc, sex, bmi, ecog_ps_cat, smoking, utuc, nmibc, bcg, multi_tur, hist_rc, 
                ct, cn, type, lnd, num_lnd, div, rm_other, grade, lvi, nac, ac, hgb, lognlr, alb, stage) %>%
  tbl_summary(
    by = stage,
    label = list(
      age_rc ~ "Age at radical cystectomy, year",
      sex ~ "Sex",
      bmi ~ "Body mass index at treatment initiation",
      ecog_ps_cat ~ "ECOG performance status",
      smoking ~ "Smoking status",
      utuc ~ "Previous upper-tract urothelial carcinoma",
      nmibc ~ "Previous non-muscle invasive urothelial cancer",
      bcg ~ "Prior intravesical BCG therapy",
      multi_tur ~ "Multifocality",
      hist_rc ~ "Tumor histology of radical cystectomy specimen",
      ct ~ "Clinical T stage",
      cn ~ "Clinical N stage",
      type ~ "Surgical approach",
      lnd ~ "Lymph node dissection",
      num_lnd ~ "Number of lymph node removed",
      div ~ "Urinary tract diversion",
      rm_other ~ "Soft-tissue resective margin",
      grade ~ "Tumor grade",
      lvi ~ "Lymphovascular invasion",
      nac ~ "Neoadjuvant chemotherapy",
      ac ~ "Adjuvant chemotherapy",
      hgb ~ "Baseline hemoglobin concentration",
      lognlr ~ "Neutrophil-to-lymphocyte ratio (log-transformed)",
      alb ~ "Baseline albumin"
    ),
    type = all_dichotomous() ~ "categorical",
    statistic = list(all_continuous() ~ "{mean} ({sd})",
                     all_categorical() ~ "{p}"),
    digits = list(all_continuous() ~ c(1, 1),
                  all_categorical() ~ c(1, 1)),
    missing_text = "Missing") %>%
  add_overall(statistic = list(all_continuous() ~ "{mean} ({sd})",
                               all_categorical() ~ "{n} ({p})"),
              digits = list(all_continuous() ~ c(1, 1),
                            all_categorical() ~ c(0, 1)))

# table1 for all patients (imputed data)
tbl1_imp <- df %>%
  dplyr::select(age_rc, sex, bmi, ecog_ps_cat, smoking, utuc, nmibc, bcg, multi_tur, hist_rc, 
                ct, cn, type, lnd, num_lnd, div, rm_other, grade, lvi, nac, ac, hgb, lognlr, alb, stage) %>%
  tbl_summary(
    by = stage,
    label = list(
      age_rc ~ "Age at radical cystectomy, year",
      sex ~ "Sex",
      bmi ~ "Body mass index at treatment initiation",
      ecog_ps_cat ~ "ECOG performance status",
      smoking ~ "Smoking status",
      utuc ~ "Previous upper-tract urothelial carcinoma",
      nmibc ~ "Previous non-muscle invasive urothelial cancer",
      bcg ~ "Prior intravesical BCG therapy",
      multi_tur ~ "Multifocality",
      hist_rc ~ "Tumor histology of radical cystectomy specimen",
      ct ~ "Clinical T stage",
      cn ~ "Clinical N stage",
      type ~ "Surgical approach",
      lnd ~ "Lymph node dissection",
      num_lnd ~ "Number of lymph node removed",
      div ~ "Urinary tract diversion",
      rm_other ~ "Soft-tissue resective margin",
      grade ~ "Tumor grade",
      lvi ~ "Lymphovascular invasion",
      nac ~ "Neoadjuvant chemotherapy",
      ac ~ "Adjuvant chemotherapy",
      hgb ~ "Baseline hemoglobin concentration",
      lognlr ~ "Neutrophil-to-lymphocyte ratio (log-transformed)",
      alb ~ "Baseline albumin"
    ),
    type = all_dichotomous() ~ "categorical",
    statistic = list(all_continuous() ~ "{mean} ({sd})",
                     all_categorical() ~ "{p}"),
    digits = list(all_continuous() ~ c(1, 1),
                  all_categorical() ~ c(1, 1))) %>%
  add_overall(statistic = list(all_continuous() ~ "{mean} ({sd})",
                               all_categorical() ~ "{n} ({p})"),
              digits = list(all_continuous() ~ c(1, 1),
                            all_categorical() ~ c(0, 1)))

# weighted imputed data
df_weight <- survey::svydesign(ids = ~ id, weights = ~ weight, data = df)
wtbl <- df_weight %>%
  tbl_svysummary(
    by = stage,
    include = c(age_rc, sex, bmi, ecog_ps_cat, smoking, utuc, nmibc, bcg, multi_tur, hist_rc, 
                ct, cn, type, lnd, num_lnd, div, rm_other, grade, lvi, nac, ac, hgb, lognlr, alb),
    label = list(
      age_rc ~ "Age at radical cystectomy, year",
      sex ~ "Sex",
      bmi ~ "Body mass index at treatment initiation",
      ecog_ps_cat ~ "ECOG performance status",
      smoking ~ "Smoking status",
      utuc ~ "Previous upper-tract urothelial carcinoma",
      nmibc ~ "Previous non-muscle invasive urothelial cancer",
      bcg ~ "Prior intravesical BCG therapy",
      multi_tur ~ "Multifocality",
      hist_rc ~ "Tumor histology of radical cystectomy specimen",
      ct ~ "Clinical T stage",
      cn ~ "Clinical N stage",
      type ~ "Surgical approach",
      lnd ~ "Lymph node dissection",
      num_lnd ~ "Number of lymph node removed",
      div ~ "Urinary tract diversion",
      rm_other ~ "Soft-tissue resective margin",
      grade ~ "Tumor grade",
      lvi ~ "Lymphovascular invasion",
      nac ~ "Neoadjuvant chemotherapy",
      ac ~ "Adjuvant chemotherapy",
      hgb ~ "Baseline hemoglobin concentration",
      lognlr ~ "Neutrophil-to-lymphocyte ratio (log-transformed)",
      alb ~ "Baseline albumin"
    ),
    type = all_dichotomous() ~ "categorical",
    statistic = list(all_continuous() ~ "{mean} ({sd})",
                     all_categorical() ~ "{p}"),
    digits = list(all_continuous() ~ c(1, 1),
                  all_categorical() ~ c(1, 1))
  ) 

tbl1 <- tbl_merge(tbls = list(tbl1_crude, tbl1_imp, wtbl),
                  tab_spanner = c("Unweighted population (before imputation)", "Unweighted population (after imputation)", "Weighted population"))

# survival data
theme_gtsummary_journal(journal = "jama")
dfs_fit <- survfit(Surv(dfs, rec) ~ stage, data = df, weights = weight)
os_fit <- survfit(Surv(fu, death_all) ~ stage, data = df, weights = weight)

surv <- list(dfs_fit, os_fit) %>% 
  tbl_survfit(probs = 0.5, conf.level = 0.95,
              estimate_fun = function(x) style_number(x, digits = 1)) # median DFS and OS probabilities

# kaplan-meier estimations
dfs_fit <- survfit(Surv(dfs, rec) ~ stage, data = df)
dfs_plot_crude <- ggsurvplot(dfs_fit, 
                             risk.table = TRUE, 
                             title = "Disease-free survival",
                             subtitle = "Unweighted population",
                             xlab = "Months since radical cystectomy", 
                             ylab = "DFS probability",
                             legend.labs = c("pT3-4N0", "pTanyN1", "pTanyN2-3"),
                             legend.title = "Stage 3",
                             censor = TRUE, censor.shape = "O", censor.size = 2.2,
                             palette = "lancet", size = 0.5,  break.time.by = 12,
                             ggtheme = theme_juog(), risk.table.title = "Number at risk",
                             risk.table.col = "strata",
                             tables.height = 0.12, risk.table.fontsize = 3.0,
                             tables.theme = survminer::theme_cleantable(), conf.int = TRUE) 

dfs_plot_crude$table <- dfs_plot_crude$table + 
  theme_void(base_size = 8) + 
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

os_fit <- survfit(Surv(fu, death_all) ~ stage, data = df)
os_plot_crude <- ggsurvplot(os_fit, 
                            risk.table = TRUE, 
                            title = "Overall survival",
                            subtitle = "Unweighted population",
                            xlab = "Months since radical cystectomy", 
                            ylab = "OS probability",
                            legend.labs = c("pT3-4N0", "pTanyN1", "pTanyN2-3"),
                            legend.title = "Stage 3",
                            censor = TRUE, censor.shape = "O", censor.size = 2.2,
                            palette = "lancet", size = 0.5,  break.time.by = 12,
                            ggtheme = theme_juog(), risk.table.title = "Number at risk",
                            risk.table.col = "strata",
                            tables.height = 0.12, risk.table.fontsize = 3.0,
                            tables.theme = survminer::theme_cleantable(), conf.int = TRUE) 

os_plot_crude$table <- os_plot_crude$table + 
  theme_void(base_size = 8) + 
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

dfs_fit <- survfit(Surv(dfs, rec) ~ stage, data = df, weights = weight)
dfs_plot_ipw <- ggsurvplot(dfs_fit, 
                           risk.table = TRUE, 
                           title = "Disease-free survival",
                           subtitle = "Weighted population",
                           xlab = "Months since radical cystectomy", 
                           ylab = "DFS probability",
                           legend.labs = c("pT3-4N0", "pTanyN1", "pTanyN2-3"),
                           legend.title = "Stage 3", 
                           censor = TRUE, censor.shape = "O", censor.size = 2.2,
                           palette = "lancet", size = 0.5,  break.time.by = 12,
                           ggtheme = theme_juog(), risk.table.title = "Number at risk",
                           risk.table.col = "strata",
                           tables.height = 0.12, risk.table.fontsize = 3.0,
                           tables.theme = survminer::theme_cleantable(), conf.int = TRUE) 

dfs_plot_ipw$table <- dfs_plot_ipw$table + 
  theme_void(base_size = 8) + 
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

os_fit <- survfit(Surv(fu, death_all) ~ stage, data = df, weights = weight)
os_plot_ipw <- ggsurvplot(os_fit, 
                          risk.table = TRUE, 
                          title = "Overall survival",
                          subtitle = "Weighted population",
                          xlab = "Months since radical cystectomy", 
                          ylab = "OS probability",
                          legend.labs = c("pT3-4N0", "pTanyN1", "pTanyN2-3"),
                          legend.title = "Stage 3",
                          censor = TRUE, censor.shape = "O", censor.size = 2.2,
                          palette = "lancet", size = 0.5,  break.time.by = 12,
                          ggtheme = theme_juog(), risk.table.title = "Number at risk",
                          risk.table.col = "strata",
                          tables.height = 0.12, risk.table.fontsize = 3.0,
                          tables.theme = survminer::theme_cleantable(), conf.int = TRUE) 

os_plot_ipw$table <- os_plot_ipw$table + 
  theme_void(base_size = 8) + 
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

# merge...
merge <- arrange_ggsurvplots(list(dfs_plot_crude, dfs_plot_ipw, os_plot_crude, os_plot_ipw), 
                             nrow = 2, ncol = 2, print = FALSE)

# IPTW-adjusted Cox regression model
theme_gtsummary_journal(journal = "jama")
cox_rfs <- coxph(Surv(dfs, rec) ~ stage, data = df, weights = weight) %>% 
  tbl_regression(exponentiate = TRUE) %>%
  add_n(location = "level") %>% 
  add_nevent(location = "level")
cox_os <- coxph(Surv(fu, death_all) ~ stage, data = df, weights = weight) %>% 
  tbl_regression(exponentiate = TRUE) %>%
  add_n(location = "level") %>% 
  add_nevent(location = "level")

# merge table
cox_merge <- tbl_merge(tbls = list(cox_rfs, cox_os),
                       tab_spanner = c("Disease-free survival", "Overall survival")
)

# treatment effect modifier
theme_gtsummary_journal(journal = "jama")
cox_int_os_nac <- coxph(Surv(fu, death_all) ~ nac * stage, data = df, weights = weight) %>% 
  tbl_regression(exponentiate = TRUE) %>%
  add_n(location = "level") %>% 
  add_nevent(location = "level")

cox_int_os_ac <- coxph(Surv(fu, death_all) ~ ac * stage, data = df, weights = weight) %>% 
  tbl_regression(exponentiate = TRUE) %>%
  add_n(location = "level") %>% 
  add_nevent(location = "level")

cox_int_dfs_nac <- coxph(Surv(dfs, rec) ~ nac * stage, data = df, weights = weight) %>% 
  tbl_regression(exponentiate = TRUE) %>%
  add_n(location = "level") %>% 
  add_nevent(location = "level")

cox_int_dfs_ac <- coxph(Surv(dfs, rec) ~ ac * stage, data = df, weights = weight) %>% 
  tbl_regression(exponentiate = TRUE) %>%
  add_n(location = "level") %>% 
  add_nevent(location = "level")

# subgroup analyses
theme_gtsummary_journal(journal = "jama")
# neoadjuvant chemotherapy
cox_rfs_nac <- coxph(Surv(dfs, rec) ~ stage, data = df %>% filter(nac == "yes"), weights = weight) %>% 
  tbl_regression(exponentiate = TRUE) %>%
  add_n(location = "level") %>% 
  add_nevent(location = "level")
cox_os_nac <- coxph(Surv(fu, death_all) ~ stage, data = df %>% filter(nac == "yes"), weights = weight) %>% 
  tbl_regression(exponentiate = TRUE) %>%
  add_n(location = "level") %>% 
  add_nevent(location = "level")

cox_rfs_wonac <- coxph(Surv(dfs, rec) ~ stage, data = df %>% filter(nac == "no"), weights = weight) %>% 
  tbl_regression(exponentiate = TRUE) %>%
  add_n(location = "level") %>% 
  add_nevent(location = "level")
cox_os_wonac <- coxph(Surv(fu, death_all) ~ stage, data = df %>% filter(nac == "no"), weights = weight) %>% 
  tbl_regression(exponentiate = TRUE) %>%
  add_n(location = "level") %>% 
  add_nevent(location = "level")

# adjuvant chemotherapy
cox_rfs_ac <- coxph(Surv(dfs, rec) ~ stage, data = df %>% filter(ac == "yes"), weights = weight) %>% 
  tbl_regression(exponentiate = TRUE) %>%
  add_n(location = "level") %>% 
  add_nevent(location = "level")
cox_os_ac <- coxph(Surv(fu, death_all) ~ stage, data = df %>% filter(ac == "yes"), weights = weight) %>% 
  tbl_regression(exponentiate = TRUE) %>%
  add_n(location = "level") %>% 
  add_nevent(location = "level")

cox_rfs_woac <- coxph(Surv(dfs, rec) ~ stage, data = df %>% filter(ac == "no"), weights = weight) %>% 
  tbl_regression(exponentiate = TRUE) %>%
  add_n(location = "level") %>% 
  add_nevent(location = "level")
cox_os_woac <- coxph(Surv(fu, death_all) ~ stage, data = df %>% filter(ac == "no"), weights = weight) %>% 
  tbl_regression(exponentiate = TRUE) %>%
  add_n(location = "level") %>% 
  add_nevent(location = "level")

cox_merge_nac <- tbl_merge(list(cox_rfs_nac, cox_os_nac), tab_spanner = c("Disease-free survival", "Overall survival"))
cox_merge_wonac <- tbl_merge(list(cox_rfs_wonac, cox_os_wonac), tab_spanner = c("Disease-free survival", "Overall survival"))

cox_merge_ac <- tbl_merge(list(cox_rfs_ac, cox_os_ac), tab_spanner = c("Disease-free survival", "Overall survival"))
cox_merge_woac <- tbl_merge(list(cox_rfs_woac, cox_os_woac), tab_spanner = c("Disease-free survival", "Overall survival"))

cox_nac_merge <- tbl_stack(list(cox_merge_nac, cox_merge_wonac), 
                           group_header = c("Neoadjuvant chemotherapy", "Without adjuvant chemotherapy"))

cox_ac_merge <- tbl_stack(list(cox_merge_ac, cox_merge_woac), 
                          group_header = c("Adjuvant chemotherapy", "Without adjuvant chemotherapy"))
