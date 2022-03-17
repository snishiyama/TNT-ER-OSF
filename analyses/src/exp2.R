
# setup -------------------------------------------------------------------

if (!exists('load_data_e2')){
  source(here::here('analyses/src/setup.R'))
}

load_data_e2 <- function(phase, rm_col = T) {
  if (rm_col) {
    rm_col_e2 <- c("IAPS", "cue", "valence", "gender", "age", 'CB')
  } else {
    rm_col_e2 <- NULL
  }
  load_data(
    phase, 
    expID = "inzemi-2",
    rm_subj = c(3,5,7,11,13,16,18,204,29,30,32,33,38,41,42),
    rm_col = rm_col_e2
  )
}


# loading data ------------------------------------------------------------

df_encoding_e2 <- load_data_e2("encoding", rm_col = F)
df_cycles_e2 <- load_data_e2("cycles")
df_TNT_e2 <- load_data_e2("TNT") %>% dplyr::filter(block > 0)

df_recog_e2 <-
  map2(
    .x = c('pre-recog', 'post-recog'), 
    .y = c('pre', 'post'), 
    .f = ~ load_data_e2(.x) %>% dplyr::mutate(time = .y)
  ) %>% 
  dplyr::bind_rows() %>% 
  dplyr::mutate(TNTcond = if_else(str_detect(TNTcond, "new"), "new", TNTcond)) %>% 
  dplyr::mutate_if(is_character, as_factor)

df_rating_e2 <-
  map2(
    .x = c('pre-rating', 'post-rating'), 
    .y = c('pre', 'post'), 
    .f = ~ load_data_e2(.x) %>% dplyr::mutate(time = .y)
  ) %>% 
  dplyr::bind_rows() %>% 
  dplyr::mutate_if(is_character, as_factor)

# Questionnaires
df_rrs_e2 <- load_data_e2("rrs", rm_col = F) %>% dplyr::filter(cond == 'b') %>% calc_score("rrs")
df_stais_e2 <- load_data_e2("staiS", rm_col = F) %>% rev_value() %>% calc_score("stais")
df_stait_e2 <- load_data_e2("staiT", rm_col = F) %>% rev_value() %>% calc_score("stait")
df_bdi_e2 <- load_data_e2("bdi", rm_col = F) %>% calc_score("bdi")

df_ques_e2 <- 
  dplyr::bind_rows(df_rrs_e2, df_stais_e2, df_stait_e2, df_bdi_e2) %>% 
  tidyr::spread(key = ques, value = score)


# demographic -------------------------------------------------------------

df_demogra_e2 <- 
  df_encoding_e2 %>% 
  dplyr::select(subjID, gender, age) %>% 
  dplyr::distinct_all() %>% 
  dplyr::summarise(N = n(), male = sum(gender), Mage = mean(age), SDage = sd(age))


# intrusion rate ----------------------------------------------------------

df_recog_pre_e2 <- df_recog_e2 %>% 
  dplyr::filter(time == "pre") %>% 
  dplyr::select(subjID, itemID, correct)

# count NAs
df_intr_NA_e2 <- df_TNT_e2 %>%
  dplyr::filter(is.na(intrusion)) %>%
  dplyr::count(subjID, block, intrusion) %>%
  dplyr::arrange(desc(n))


df_intr_woex_e2 <- 
  df_TNT_e2 %>% 
  dplyr::left_join(df_recog_pre_e2, by = c("subjID", "itemID")) %>% 
  dplyr::mutate(inORnot = if_else(intrusion == 0, 0, 1)) %>% # 1 & 2 are coded 1 (i.e., intruded)
  dplyr::filter(!is.na(intrusion)) %>% 
  dplyr::select(TNTcond, itemID, subjID, block, correct, inORnot)


# mean across trials
df_intr_TNT_e2 <-  
  df_intr_woex_e2 %>% 
  dplyr::filter(correct == 1) %>% 
  dplyr::group_by(subjID,block,TNTcond) %>% 
  dplyr::summarise(int_rate = mean(inORnot), .groups = "drop")

df_intr_e2_w <- df_intr_TNT_e2 %>% 
  dplyr::filter(block %in% c(1, 10), TNTcond == 'nt') %>%
  tidyr::pivot_wider(names_from = block, values_from = int_rate, names_prefix = 'block')

ttest_intr_TNT_e2 <- t.test(df_intr_e2_w$block1, df_intr_e2_w$block10, alternative = 'greater')

smry_intr_e2 <- get_t_g_ci(df_intr_e2_w$block1, df_intr_e2_w$block10, method = 'greater')

# df_intr_TNT_e2 %>% 
#   dplyr::filter(TNTcond == 'nt') %>% 
#   dplyr::group_by(block, correct) %>% 
#   dplyr::summarise(intr = mean(int_rate)) %>% 
#   ggplot(aes(x = as.factor(block), y = intr, color = as.factor(correct), group = as.factor(correct))) +
#   geom_line()

# recognition -------------------------------------------------------------

df_hit_rate_e2 <- df_recog_e2 %>% 
  dplyr::select(subjID, time, TNTcond, itemID, correct) %>% 
  tidyr::pivot_wider(names_from = time, values_from = correct) %>%
  dplyr::filter(pre == 1, TNTcond %in% c('nt', 'b')) %>% 
  dplyr::group_by(TNTcond, subjID) %>% 
  dplyr::summarise(hit_rate = mean(post), .groups = "drop")

# t test
df_hit_rate_e2_w <- df_hit_rate_e2 %>% 
  tidyr::pivot_wider(names_from = TNTcond, values_from = hit_rate)
ttest_hit_rate_e2 <- t.test(df_hit_rate_e2_w$nt, df_hit_rate_e2_w$b, paired = T, alternative = 'less')

## create summary description
smry_hit_rate_e2 <- get_t_g_ci(df_hit_rate_e2_w$nt, df_hit_rate_e2_w$b, method = 'less')

# Number of dropped items
df_hit_dropped_e2 <- df_recog_e2 %>% 
  dplyr::select(subjID, time, TNTcond, itemID, correct) %>% 
  tidyr::spread(key = time, value = correct) %>% 
  dplyr::filter(pre == 0, TNTcond %in% c('nt','b')) %>% 
  dplyr::group_by(TNTcond, subjID) %>% 
  dplyr::tally(name = 'dropped_n') %>% 
  dplyr::summarise(mean = mean(dropped_n), sd = sd(dropped_n))

# substitutes
df_hit_sub_e2 <- df_recog_e2 %>% 
  dplyr::filter(TNTcond == 'sub') %>% 
  dplyr::select(-order, -TNTcond) %>% 
  # add TNTcond column to identify which Think or No-Think items substitutes were associated with
  dplyr::left_join(
    df_recog_e2 %>% dplyr::filter(TNTcond != 'sub', time == "pre") %>% select(TNTcond, subjID, itemID, copre = correct),
    by = c('subjID', 'itemID')
  ) %>% 
  dplyr::filter(copre == 1) %>% 
  dplyr::group_by(TNTcond, subjID) %>% 
  dplyr::summarise(hit_rate = mean(correct), .groups = "drop")

# for Table 
df_hit_sub_e2 %>%
  dplyr::group_by(TNTcond) %>%
  dplyr::summarise(mean = mean(hit_rate), sd = sd(hit_rate)) %>%
  dplyr::mutate_if(.predicate = is_double, round, digits = 2)

## ttest
df_hit_sub_e2_w <- df_hit_sub_e2 %>% 
  tidyr::pivot_wider(names_from = TNTcond, values_from = hit_rate)
t.test(df_hit_sub_e2_w$nt, df_hit_sub_e2_w$t, paired = T, alternative = 'greater')
## create summary description
smry_hit_sub_e2 <- get_t_g_ci(df_hit_sub_e2_w$nt, df_hit_sub_e2_w$t, method = 'greater')

## forced-choice recognition
df_recall_rate_e2 <- df_rating_e2 %>% 
  dplyr::select(TNTcond, itemID, subjID, correct, time) %>% 
  tidyr::pivot_wider(names_from = time, values_from = correct) %>% 
  dplyr::filter(pre == 1, TNTcond %in% c('b', 'nt', 't')) %>% 
  dplyr::group_by(subjID, TNTcond) %>% 
  dplyr::summarise(recall_rate = mean(post))

df_recall_e2_w <- df_recall_rate_e2 %>%
  tidyr::pivot_wider(names_from = TNTcond, values_from = recall_rate)
ttest_recall_rate_e2 <- t.test(df_recall_e2_w$nt, df_recall_e2_w$b, paired = T, data = ., alternative = 'less')
smry_recall_e2 <- get_t_g_ci(df_recall_e2_w$nt, df_recall_e2_w$b, method = 'less')

# affective rating --------------------------------------------------------

# Included 
# 1. scores for items that were recognized in the pre-TNT recognition test
# 2. scores for items that were rated in both pre- and post- ratings
# 3. scores from items that were recognized correctly att the end of the rating sequence.

df_recog_pre_e2 <- df_recog_e2 %>% 
  dplyr::filter(time == "pre") %>% 
  dplyr::select(subjID, TNTcond, itemID, correct_pre = correct)

df_rating_tntb_e2 <- df_rating_e2 %>% 
  dplyr::left_join(df_recog_pre_e2, by = c("subjID", "TNTcond", "itemID")) %>% 
  dplyr::filter(correct_pre == 1,
                correct == 1,
                TNTcond %in% c('t', 'nt', 'b')) %>% 
  dplyr::select(subjID, time, TNTcond, val, aro, itemID) %>% 
  tidyr::pivot_longer(cols = c(val, aro), names_to = "dim") %>% 
  tidyr::pivot_wider(names_from = time, values_from = value)

l_df_rating_subj_mean_e2 <- df_rating_tntb_e2 %>% 
  tidyr::drop_na() %>% 
  dplyr::mutate(diff = post - pre) %>% 
  dplyr::group_by(subjID, TNTcond, dim) %>% 
  dplyr::summarise(across(.cols = c(pre, post, diff), .fns = mean), .groups = "drop") %>% 
  split(.$dim)

# create summary description
## valence
df_val_e2_w <- l_df_rating_subj_mean_e2$val %>% 
  dplyr::select(subjID, TNTcond, diff) %>% 
  tidyr::pivot_wider(names_from = TNTcond, values_from = diff)
ttest_val_e2 <- t.test(df_val_e2_w$nt, df_val_e2_w$b, paired = T, alternative = 'greater')
smry_val_e2 <- get_t_g_ci(df_val_e2_w$nt, df_val_e2_w$b, method = 'greater')

## arousal
df_aro_e2_w <-  l_df_rating_subj_mean_e2$aro %>% 
  dplyr::select(subjID, TNTcond, diff) %>% 
  tidyr::pivot_wider(names_from = TNTcond, values_from = diff)
ttest_aro_e2 <- t.test(df_aro_e2_w$nt, df_aro_e2_w$b, paired = T, alternative = 'less')
smry_aro_e2 <- get_t_g_ci(df_aro_e2_w$nt, df_aro_e2_w$b, method = 'less')


# correlation -------------------------------------------------------------

df_indiv_ques_e2 <- df_ques_e2 %>% 
  dplyr::arrange(subjID) %>% 
  dplyr::mutate(
    # Positive value means beneficial consequence by memory contorl regarding the following variable 
    aro = df_aro_e2_w %>% arrange(subjID) %>% mutate(diff = b - nt) %>% pull(diff),
    val = df_val_e2_w %>% arrange(subjID) %>% mutate(diff = nt - b) %>% pull(diff)
  )

# saveRDS(df_indiv_ques_e2, here::here("analyses/rds/df_indiv_ques_e2.rds"))

# df_ques_e2 %>% 
#   dplyr::left_join(df_recall_e2_w, by = "subjID") %>% 
#   dplyr::mutate(forget = b - nt) %>% 
#   dplyr::select(bdi:stait, forget) %>% 
#   psych::pairs.panels(lm = TRUE, ci = TRUE, ellipses = F)

# correlation (paris panel)
# df_indiv_ques_e2 %>%
#   dplyr::select(-subjID) %>%
#   psych::pairs.panels(lm = TRUE, ci = TRUE, ellipses = F)

res_pca_e2 <- df_ques_e2 %>% 
  dplyr::select(BDI = bdi, RRS = rrs, `STAI-S` = stais, `STAI-T` = stait) %>% 
  stats::prcomp(scale = T)

smry_pca_e2 <- summary(res_pca_e2)

df_indiv_ques_all_e2 <- df_indiv_ques_e2 %>%
  dplyr::mutate(pc1 = res_pca_e2$x[,1], pc2 = res_pca_e2$x[,2])

# valence
lmfit_val_e2 <- lm(val ~ pc1 + pc2, data = df_indiv_ques_all_e2)
smry_lmfit_val_e2 <- summary(lmfit_val_e2)
r2_lmfit_val_e2 <- get_reg_stats(lmfit_val_e2)

# arousal
lmfit_aro_e2 <- lm(aro ~ pc1 + pc2, data = df_indiv_ques_all_e2)
smry_lmfit_aro_e2 <- summary(lmfit_aro_e2)
r2_lmfit_aro_e2 <- get_reg_stats(lmfit_aro_e2)


# tables ------------------------------------------------------------------

# old/new recognition
df_recog_e2 %>%
  dplyr::select(subjID, time, TNTcond, correct) %>%
  dplyr::mutate(time = factor(time, levels = c("pre", "post"))) %>%
  dplyr::filter(TNTcond != "filler") %>%
  dplyr::group_by(time, TNTcond, subjID) %>%
  dplyr::summarise(rate = mean(correct), .groups = "drop_last") %>%
  dplyr::summarise(mean = mean(rate), sd = sd(rate), .groups = "drop") %>%
  tidyr::pivot_longer(cols = mean:sd, names_to = "index") %>%
  tidyr::pivot_wider(names_from = TNTcond, values_from = value) %>%
  dplyr::mutate_if(.predicate = is.double, .funs = round, digits = 2) %>%
  dplyr::select(time, index, t, nt, b, sub, new)

# forced-choice recognition
df_rating_e2 %>%
  dplyr::mutate(time = factor(time, levels = c("pre", "post"))) %>%
  dplyr::filter(TNTcond != "filler") %>%
  dplyr::group_by(time, TNTcond, subjID) %>%
  dplyr::summarise(rate = mean(correct), .groups = "drop_last") %>%
  dplyr::summarise(mean = mean(rate), sd = sd(rate), .groups = "drop") %>%
  tidyr::pivot_longer(cols = mean:sd, names_to = "index") %>%
  tidyr::pivot_wider(names_from = TNTcond, values_from = "value") %>%
  dplyr::mutate_if(.predicate = is.double, .funs = round, digits = 2)

# affective rating
df_rating_tntb_e2 %>% 
  tidyr::drop_na() %>% 
  # summarize across trials
  dplyr::group_by(TNTcond, dim, subjID) %>% 
  dplyr::summarise(across(.cols = c(pre, post), .fns =  mean), 
                   .groups = "drop_last") %>% 
  # summarize across participants
  dplyr::summarise(across(.cols = c(pre, post), .fns = list("mean" = mean, "sd" = sd)),
                   .groups = "drop") %>% 
  tidyr::pivot_longer(cols = pre_mean:post_sd, names_sep = "_", names_to = c("time", "stat")) %>% 
  tidyr::pivot_wider(names_from = c(dim, TNTcond), values_from = value) %>%
  dplyr::select(time, stat, val_b, val_nt, val_t, aro_b, aro_nt, aro_t) %>%
  knitr::kable(format = "pandoc")

# regression
dplyr::bind_rows(
  list(val = lm_to_table(lmfit_val_e2), aro = lm_to_table(lmfit_aro_e2)), 
  .id = "measure") %>% 
  knitr::kable(format = "pandoc")

# for analyses on JASP
df_jasp_e2 <- l_df_rating_subj_mean_e2 %>% 
  dplyr::bind_rows(.id = "dim") %>% 
  dplyr::select(subjID, dim, TNTcond, diff) %>% 
  tidyr::pivot_wider(names_from = c(dim, TNTcond), names_prefix = "diff_", values_from = diff) %>% 
  dplyr::left_join(df_ques_e2, by = "subjID")

# write_csv(df_jasp_e2, here::here("analyses/jasp/df_jasp_e2.csv"))

# plots -------------------------------------------------------------------

gg_intr_TNT_e2 <- make_plot_intr(df_intr_TNT_e2) +
  theme(legend.position = c(0.8,0.2))

gg_hit_rate_e2 <- make_plot_hitrate(df_hit_rate_e2)

gg_aro_diff_e2 <- make_plot_rating(l_df_rating_subj_mean_e2$aro) +
  labs(y = expression(paste(Delta, "Arousal (post-pre)")))
gg_val_diff_e2 <- make_plot_rating(l_df_rating_subj_mean_e2$val) + 
  labs(y = expression(paste(Delta, "Valence (post-pre)")))

gg_pairplot_e2 <- make_pairplot(df_indiv_ques_e2)

upper_grid_e2 <- cowplot::plot_grid(gg_intr_TNT_e2, gg_hit_rate_e2, gg_val_diff_e2, gg_aro_diff_e2,
                                    labels = c("A", "B", "C", "D"))
plot_all_e2 <- cowplot::plot_grid(upper_grid_e2, ggmatrix_gtable(gg_pairplot_e2), ncol = 1, labels = c("", "E"))

# ggsave( here::here("analyses/figures/grid_plot_e2.png"), plot_all_e2, width = 10.6, height = 11.4)

gg_res_pca_e2 <- make_pca_biplot(res_pca_e2)
# ggsave(filename = here::here("analyses/figures/biplot_pca_e2.png"), plot = gg_res_pca_e2, width = 6, height = 4)

gg_corr_cued_recog <- df_ques_e2 %>% 
  dplyr::left_join(df_recall_e2_w, by = "subjID") %>%
  dplyr::left_join(df_hit_rate_e2_w %>% dplyr::mutate(`Item Recog` = b - nt, .keep = "unused"),
                   by = "subjID") %>% 
  dplyr::left_join(df_intr_TNT_e2 %>%
                     dplyr::filter(TNTcond == 'nt') %>%
                     dplyr::group_by(subjID) %>%
                     dplyr::summarise(mean_int = mean(int_rate)),
                   by = "subjID") %>%
  dplyr::mutate(`Cued Recog` = b - nt) %>%
  dplyr::select(BDI = bdi, 
                RRS = rrs, 
                `STAI-S` = stais, 
                `STAI-T` = stait, 
                `Intrusion` = mean_int,
                `Item Recog`, 
                `Cued Recog`) %>% 
  GGally::ggpairs(lower = list(continuous = wrap("smooth", alpha = 0.3)),
                  upper = list(continuous = wrap("cor", color = "black"))) +
  theme_bw(base_size = font_size, base_family = "Helvetica Neue") +
  theme(axis.text = element_text(colour = "black"),
        panel.grid = element_blank(),
        strip.background = element_blank())

# ggsave(filename = here::here("analyses/figures/pairs_plot_e2_supple.png"), plot = gg_corr_cued_recog, width = 8, height = 8)

