
# Set up ------------------------------------------------------------------

if (!exists('load_data_e2')){
  source(here::here('analyses/src/setup.R'))
}
if (!exists('df_indiv_ques_e1')) {
  df_indiv_ques_e1 <- readRDS(here::here('analyses/rds/df_indiv_ques_e1.rds'))
} 
if (!exists('df_indiv_ques_e2')) {
  df_indiv_ques_e2 <- readRDS(here::here('analyses/rds/df_indiv_ques_e2.rds'))
}


# Compute correlations ----------------------------------------------------

df_indiv_ques_exps <- dplyr::bind_rows(
  df_indiv_ques_e1 %>% mutate(group = 'ds'),
  df_indiv_ques_e2 %>% mutate(group = 'ts')
)

res_pca_exps <- df_indiv_ques_exps %>% 
  dplyr::select(BDI = bdi, 
                RRS = rrs, 
                `STAI-S` = stais, 
                `STAI-T` = stait) %>% 
  stats::prcomp(scale = T)

smry_pca_exps <- summary(res_pca_exps)

df_indiv_ques_all_exps <- df_indiv_ques_exps %>% 
  dplyr::mutate(pc1 = res_pca_exps$x[,1],
                pc2 = res_pca_exps$x[,2],
                group = as_factor(group))


df_indiv_lm <- df_indiv_ques_all_exps %>% 
  tidyr::pivot_longer(cols = c(val, aro), names_to = 'vars', values_to = "value") %>% 
  dplyr::group_by(vars) %>%
  tidyr::nest() %>% 
  dplyr::mutate(
    step1 = purrr::map(
      .x = data, 
      .f = ~ lm(formula = value ~ group + pc1 + pc2,
                contrasts = list(group = -contr.sum(2)/2), 
                data = .x)
      ),
    step2 = purrr::map(
      .x = data, 
      .f = ~ lm(formula = value ~ group + pc1 + pc2 + group:pc1 + group:pc2,
                contrasts = list(group = -contr.sum(2)/2), 
                data = .x)
    ),
    model_comparison = purrr::pmap(
      .l = list(step1, step2),
      .f = anova
    )
  )

# results of moderation in valence
smry_lmfit_step1_val_both <- summary(df_indiv_lm$step1[[1]])
smry_lmfit_step2_val_both <- summary(df_indiv_lm$step2[[1]])
r2_val_step1_both <- get_reg_stats(df_indiv_lm$step1[[1]])
r2_val_step2_both <- get_reg_stats(df_indiv_lm$step2[[1]])
model_1vs2_val_both <- get_aov_stats(df_indiv_lm$model_comparison[[1]], 2)

cff_group_step1_val_both <- get_cff_stats(df_indiv_lm$step1[[1]], "group1")
cff_pc1_step1_val_both <- get_cff_stats(df_indiv_lm$step1[[1]], "pc1")

# results of moderation in arousal
smry_lmfit_step1_aro_both <- summary(df_indiv_lm$step1[[2]])
smry_lmfit_step2_aro_both <- summary(df_indiv_lm$step2[[2]])
r2_aro_step1_both <- get_reg_stats(df_indiv_lm$step1[[2]])
r2_aro_step2_both <- get_reg_stats(df_indiv_lm$step2[[2]])
model_1vs2_aro_both <- get_aov_stats(df_indiv_lm$model_comparison[[2]], 2)

# two-way interaction
cff_groupBYpc1_step2_aro <- get_cff_stats(df_indiv_lm$step2[[2]], "group1:pc1")

res_simp_slope_pc1group <- df_indiv_ques_all_exps %>% 
  dplyr::mutate(group = if_else(group == 'ds', -0.5, 0.5)) %>% 
  pequod::lmres(aro ~ group + pc1 + pc2 + group:pc1 + group:pc2, data = .) %>% 
  pequod::simpleSlope(pred = "pc1", mod1 = "group")

cff_pc1DS_step2_aro <- get_ssa_stats(res_simp_slope_pc1group, 1)
cff_pc1TS_step2_aro <- get_ssa_stats(res_simp_slope_pc1group, 2)
p.adjust(res_simp_slope_pc1group$simple_slope[, "p.value"], method = "holm") %>% 
  round(3)

# Code for Regression Table -----------------------------------------------

tbl_lm_val_step1_both <- lm_to_table(df_indiv_lm$step1[[1]])
tbl_lm_aro_step1_both <- lm_to_table(df_indiv_lm$step1[[2]])
tbl_lm_val_step2_both <- lm_to_table(df_indiv_lm$step2[[1]])
tbl_lm_aro_step2_both <- lm_to_table(df_indiv_lm$step2[[2]])


tbl_lm_val_steps_both <- right_join(tbl_lm_val_step1_both, tbl_lm_val_step2_both, by = c("term")) %>% 
  dplyr::mutate(measure = "val")
tbl_lm_aro_steps_both <- right_join(tbl_lm_aro_step1_both, tbl_lm_aro_step2_both, by = c("term")) %>% 
  dplyr::mutate(measure = "aro")

dplyr::bind_rows(tbl_lm_val_steps_both, tbl_lm_aro_steps_both) %>% 
  knitr::kable(format = "pandoc")

# model comparison
l_model_stats_1 <- purrr::map(df_indiv_lm$step1, .f = summary) %>%
  purrr::map(.f = ~ tibble(step = 1, r_adj = .x$adj.r.squared, p = 1 - pf(.x$fstatistic["value"], .x$fstatistic["numdf"], .x$fstatistic["dendf"])))

l_model_stats_2 <- purrr::map(df_indiv_lm$step2, .f = summary) %>%
  purrr::map(.f = ~ tibble(step = 2, r_adj = .x$adj.r.squared, p = 1 - pf(.x$fstatistic["value"], .x$fstatistic["numdf"], .x$fstatistic["dendf"])))


purrr::pmap(list(l_model_stats_1, l_model_stats_2), .f = bind_rows) %>%
  purrr::map(.f = mutate, diff = r_adj - lag(r_adj)) %>%
  purrr::map2(.y = df_indiv_lm$model_comparison,
              .f = ~mutate(.x, p_diff = .y$`Pr(>F)`)) %>%
  dplyr::bind_rows(.id = "dim") %>%
  tidyr::pivot_longer(cols = r_adj:p_diff) %>%
  dplyr::mutate(dim = dplyr::recode(dim, "1" = "val", "2" = "aro"),
                type = if_else(str_detect(name, "p"), "p", "r"),
                type2 = if_else(str_detect(name, "diff"), "diff", "adj")) %>%
  dplyr::select(dim, step, type, type2, value) %>%
  tidyr::pivot_wider(names_from = c(type, step), values_from = value) %>%
  dplyr::mutate(across(.cols = starts_with("p"), .fns = format_pval_2)) %>%
  knitr::kable(format = "pandoc", digits = 2)


# plots -------------------------------------------------------------------

gg_res_pca_both <- make_pca_biplot(res_pca_exps)
# ggsave(filename = here::here("analyses/figures/biplot_pca_both.png"), plot = gg_res_pca_both, width = 6, height = 4)
