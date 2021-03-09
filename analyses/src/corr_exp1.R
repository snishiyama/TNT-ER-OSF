
# Set up ------------------------------------------------------------------

if (!exists('df_intr_e1_w')) {
  source(here::here('analyses/src/intrusion-rate_exp1.R'))
} 
if (!exists('df_hit_rate_e1_w')) {
  source(here::here('analyses/src/recog_exp1.R'))
} 
if (!exists('df_aro_e1_w')) {
  source(here::here('analyses/src/rating_exp1.R'))
}


# Compute correlations ----------------------------------------------------

df_indiv_ques_e1 <- 
  df_ques_e1 %>% 
  dplyr::arrange(subjID) %>%
  dplyr::ungroup() %>% 
  dplyr::mutate(
    # Positive value means beneficial consequence by memory contorl regarding the following variable 
    aro = df_aro_e1_w %>% arrange(subjID) %>% mutate(diff = b - nt) %>% pull(diff),
    val = df_val_e1_w %>% arrange(subjID) %>% mutate(diff = nt - b) %>% pull(diff)
  )

# correlation (pairs panel)
# df_indiv_ques_e1 %>%
#   dplyr::select(-subjID) %>%
#   psych::pairs.panels(lm = TRUE, ci = TRUE, ellipses = F)


# PCA and regression ------------------------------------------------------

# PCA
res_pca_e1 <- 
  df_ques_e1 %>% 
  dplyr::ungroup() %>% 
  dplyr::select(BDI = bdi, 
                RRS = rrs, 
                `STAI-S` = stais, 
                `STAI-T` = stait) %>% 
  stats::prcomp(scale = T)

smry_pca_e1 <- summary(res_pca_e1)

# Regression
df_indiv_ques_all_e1 <- 
  df_indiv_ques_e1 %>%
  dplyr::mutate(pc1 = res_pca_e1$x[,1], pc2 = res_pca_e1$x[,2])

df_indiv_lm_e1 <- 
  df_indiv_ques_all_e1 %>% 
  dplyr::select(subjID, val, aro,  pc1, pc2) %>% 
  tidyr::pivot_longer(cols = val:aro, names_to = 'vars', values_to = "value") %>% 
  dplyr::group_by(vars) %>%
  tidyr::nest() %>% 
  dplyr::mutate(
    step1 = purrr::map(
      .x = data, 
      .f = ~ lm(formula = value ~ pc1 + pc2,
                data = .x)
    ),
    step2 = purrr::map(
      .x = data, 
      .f = ~ lm(formula = value ~ pc1 + pc2 + pc1:pc2,
                data = .x)
    ),
    model_comparison = purrr::pmap(
      .l = list(step1, step2),
      .f = anova
    )
  )

# results of moderation in valence
r2_val_step1_e1 <- get_reg_stats(df_indiv_lm_e1$step1[[1]])
r2_val_step2_e1 <- get_reg_stats(df_indiv_lm_e1$step2[[1]])
model_comp_val_e1 <- get_aov_stats(df_indiv_lm_e1$model_comparison[[1]], 2)

cff_pc1_step1_val_e1 <- get_cff_stats(df_indiv_lm_e1$step1[[1]], "pc1")

# results of moderation in arousal
r2_aro_step1_e1 <- get_reg_stats(df_indiv_lm_e1$step1[[2]])
r2_aro_step2_e1 <- get_reg_stats(df_indiv_lm_e1$step2[[2]])
model_comp_aro_e1 <- get_aov_stats(df_indiv_lm_e1$model_comparison[[2]], 2)

cff_pc1_step1_aro_e1 <- get_cff_stats(df_indiv_lm_e1$step1[[2]], "pc1")

# simple slope analysis
res_simp_slope_e1 <-
  df_indiv_ques_all_e1 %>% 
  pequod::lmres(aro ~ pc1 * pc2, data = .) %>% 
  pequod::simpleSlope(pred = "pc1", mod1 = "pc2")

cff_pc1Lowpc2_step2_aro <- get_ssa_stats(res_simp_slope_e1, 1)
cff_pc1Highpc2_step2_aro <- get_ssa_stats(res_simp_slope_e1, 2)

# Code for Regression Table -----------------------------------------------

l_smry_lm_step1_e1 <- lm_to_table(df_indiv_lm_e1$step1)
l_smry_lm_step2_e1 <- lm_to_table(df_indiv_lm_e1$step2)

# purrr::map2(l_smry_lm_step1_e1, l_smry_lm_step2_e1, .f = right_join, by = c("term"))  %>%
#   purrr::map2(.y = c("val", "aro"), .f = ~mutate(.x, measure = .y)) %>%
#   dplyr::bind_rows() %>%
#   knitr::kable(format = "pandoc")

# model comparison
l_model_stats_step1_e1 <- purrr::map(df_indiv_lm_e1$step1, .f = summary) %>%
  purrr::map(.f = ~ tibble(step = 1, r_adj = .x$adj.r.squared, p = 1 - pf(.x$fstatistic["value"], .x$fstatistic["numdf"], .x$fstatistic["dendf"])))

l_model_stats_step2_e1 <- purrr::map(df_indiv_lm_e1$step2, .f = summary) %>%
  purrr::map(.f = ~ tibble(step = 2, r_adj = .x$adj.r.squared, p = 1 - pf(.x$fstatistic["value"], .x$fstatistic["numdf"], .x$fstatistic["dendf"])))

# purrr::pmap(list(l_model_stats_step1_e1, l_model_stats_step2_e1), .f = bind_rows) %>%
#   purrr::map(.f = mutate, diff = r_adj - lag(r_adj)) %>%
#   purrr::map2(.y = df_indiv_lm_e1$model_comparison,
#               .f = ~mutate(.x, p_diff = .y$`Pr(>F)`)
#   ) %>%
#   purrr::map(tidyr::pivot_longer, cols = -step) %>%
#   purrr::map(mutate,
#              type = if_else(str_detect(name, "p"), "p", "r"),
#              type2 = if_else(str_detect(name, "diff"), "diff", "adj")) %>%
#   purrr::map(select, step, type, type2, value) %>%
#   purrr::map(tidyr::pivot_wider, names_from = c(type, step), values_from = value) %>%
#   purrr::map(mutate_at, .vars = c("p_1", "p_2"), .funs = format_pval_2) %>%
#   purrr::map(.f = knitr::kable, format = "pandoc", digits = 2)
