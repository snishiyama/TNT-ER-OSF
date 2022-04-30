
# setup -------------------------------------------------------------------

if (!exists('load_data_e1')){
  source(here::here('analyses/src/setup.R'))
}

load_data_e1 <- function(phase, col_names = T, rm_col = T) {
  if (rm_col) {
    rm_col_e1 <- c("tar", "valence", "gender", "age")
  } else {
    rm_col_e1 <- NULL
  }
  load_data(
    phase, 
    expID = "inzemi-2018", 
    col_names = col_names,
    rm_col = rm_col_e1
  )
}

# excluded ----------------------------------------------------------------

excluded <- c(
  "6_excluded/6_TNT_inzemi-2018.csv", "9_excluded/9_TNT_inzemi-2018.csv", 
  "24_excluded/24_TNT_inzemi-2018.csv", "24_excluded-2/24_TNT_inzemi-2018.csv"
)

df_tnt_ex <- purrr::map(excluded, ~here::here("exp/data/excluded", .x)) %>% 
  purrr::map_dfr(read.csv, 
                 header = F, 
                 col.names = c("subjID","gender","age","phase","block","order","itemID","TNTcond","valence","cue","tar","intrusion","RT"))


df_tnt_ex %>% 
  dplyr::group_by(subjID) %>% 
  dplyr::count(intrusion)

l_df_ex <- fs::dir_ls(here::here("exp/data/excluded"), recurse = T, regexp = "(pre|post)-rating.+csv$") %>% 
  purrr::map(read.csv,
             header = F,
             col.names = c("subjID","gender","age","phase","order","itemID","TNTcond","valence","cue","tar","val","aro"))

list(s24_1 = dplyr::bind_rows(l_df_ex[[1]], l_df_ex[[2]]),
     s24_2 = dplyr::bind_rows(l_df_ex[[3]], l_df_ex[[4]]),
     s6_1 = l_df_ex[[5]],
     s6_2 = dplyr::bind_rows(l_df_ex[[6]], l_df_ex[[7]]),
     s9 = dplyr::bind_rows(l_df_ex[[8]], l_df_ex[[9]])) %>%
  purrr::map(group_by, itemID) %>% 
  purrr::map(dplyr::summarise, skipped = any(is.na(val))) %>% 
  purrr::map(dplyr::count, skipped)

# loading data ------------------------------------------------------------

df_encoding_e1 <- load_data_e1('encoding', rm_col = F)
df_cycles_e1 <- load_data_e1('cycles', c("subjID","gender","age","phase","cycleN","order","itemID","TNTcond","valence","cue","tar","correct"))
df_TNT_e1 <- load_data_e1('TNT', c("subjID","gender","age","phase","block","order","itemID","TNTcond","valence","cue","tar","intrusion","RT"))
df_TNT_e1 <-  df_TNT_e1 %>% dplyr::filter(block != 'prac') %>% dplyr::mutate(block = as.integer(block) + 1)

df_recog_e1 <- load_data_e1('(pre|post)-recog', c("subjID","gender","age","phase","order","itemID","TNTcond","valence","tar","correct","RT"))
df_recog_e1 <- tidyr::separate(df_recog_e1, col = phase, into = c('time', 'task'))

df_rating_e1 <- load_data_e1('(pre|post)-rating', c("subjID","gender","age","phase","order","itemID","TNTcond","valence","cue","tar","val","aro"))
df_rating_e1 <- df_rating_e1 %>% tidyr::separate(col = phase, into = c('time', 'task'))

df_rrs_e1 <- load_data_e1('rrs', rm_col = F) %>% dplyr::filter(cond == 'b') %>% calc_score("rrs")
df_stais_e1 <- load_data_e1('stai1', rm_col = F) %>% rev_value() %>% calc_score("stais")
df_stait_e1 <- load_data_e1('stai2', rm_col = F) %>% rev_value() %>% calc_score("stait")
df_bdi_e1 <- load_data_e1('bdi', rm_col = F) %>% calc_score("bdi")

df_ques_e1 <- 
  dplyr::bind_rows(df_rrs_e1, df_stais_e1, df_stait_e1, df_bdi_e1) %>% 
  tidyr::spread(key = ques, value = score)


# demographic -------------------------------------------------------------

df_demogra_e1 <- 
  df_encoding_e1 %>% 
  dplyr::select(subjID, gender, age) %>% 
  dplyr::distinct_all() %>% 
  dplyr::summarise(N = n(), male = sum(gender), Mage = mean(age), SDage = sd(age))


# intrusion rate ----------------------------------------------------------

df_recog_pre_e1 <- df_recog_e1 %>% 
  dplyr::filter(time == "pre") %>% 
  dplyr::select(subjID, itemID, correct)

# count NAs
df_intr_NA_e1 <- 
  df_TNT_e1 %>%
  dplyr::filter(is.na(intrusion)) %>%
  dplyr::count(subjID, block, intrusion) %>%
  dplyr::arrange(desc(n))

df_intr_woex_e1 <- 
  df_TNT_e1 %>% 
  dplyr::left_join(df_recog_pre_e1, by = c("subjID", "itemID")) %>% 
  dplyr::mutate(inORnot = if_else(intrusion == 0, 0, 1)) %>% # 1 & 2 are coded 1 (i.e., intruded)
  dplyr::filter(!is.na(intrusion)) %>% 
  dplyr::select(TNTcond, itemID, subjID, block, correct, inORnot)

# summarise across trials
df_intr_TNT_e1 <-  
  df_intr_woex_e1 %>% 
  dplyr::filter(correct == 1) %>% 
  dplyr::group_by(subjID, block, TNTcond) %>% 
  dplyr::summarise(int_rate = mean(inORnot))

# t test
df_intr_e1_w <- df_intr_TNT_e1 %>% 
  dplyr::filter(block %in% c(1, 10), TNTcond == 'nt') %>%
  tidyr::pivot_wider(names_from = block, values_from = int_rate, names_prefix = 'block')
ttest_intr_TNT_e1 <- t.test(df_intr_e1_w$block1, df_intr_e1_w$block10, paired = T, alternative = 'greater')

## create summary description
smry_intr_e1 <- get_t_g_ci(df_intr_e1_w$block1, df_intr_e1_w$block10, method = 'greater')


# recognition -------------------------------------------------------------

# We analyzed hit rates of items that were recognized in the pre-TNT recognition test.
df_hit_rate_e1 <- df_recog_e1 %>% 
  dplyr::select(subjID, time, TNTcond, itemID, correct) %>% 
  tidyr::pivot_wider(names_from = time, values_from = correct) %>% 
  # tidyr::spread(key = time, value = correct) %>% 
  dplyr::filter(pre == 1, TNTcond %in% c('nt', 'b')) %>% 
  dplyr::group_by(TNTcond, subjID) %>% 
  dplyr::summarise(hit_rate = mean(post), .groups = "drop")

# t test
ttest_hit_rate_e1 <- df_hit_rate_e1 %>%
  dplyr::mutate(TNTcond = factor(TNTcond, levels = c("nt", "b"))) %>%
  t.test(hit_rate ~ TNTcond, paired = T, data = ., alternative = 'less')

# create summary description
df_hit_rate_e1_w <- df_hit_rate_e1 %>% spread(TNTcond, value = hit_rate)
smry_hit_rate_e1 <- get_t_g_ci(df_hit_rate_e1_w$nt, df_hit_rate_e1_w$b, method = 'less')

# Number of dropped items
df_hit_dropped_e1 <- 
  df_recog_e1 %>% 
  dplyr::select(subjID, time, TNTcond, itemID, correct) %>% 
  tidyr::pivot_wider(names_from = time, values_from = correct) %>% 
  dplyr::filter(pre == 0, TNTcond %in% c('nt','b')) %>% 
  dplyr::group_by(TNTcond, subjID) %>% 
  dplyr::tally(name = 'dropped_n') %>% 
  dplyr::summarise(mean = mean(dropped_n), sd = sd(dropped_n), .groups = "drop")

# affective rating --------------------------------------------------------

# Included 
# 1. scores for items that were recognized in the pre-TNT recognition test
# 2. scores for items that were rated in both pre- and post- rating

# skipped trials
df_rating_e1 %>% 
  dplyr::group_by(subjID, itemID) %>% 
  dplyr::summarise(skipped = any(is.na(val))) %>% 
  dplyr::count(skipped) %>% 
  dplyr::ungroup() %>% 
  tidyr::complete(subjID, nesting(skipped), fill = list(n = 0)) %>% 
  dplyr::group_by(skipped) %>% 
  dplyr::summarise(mean = mean(n), median = median(n), sd = sd(n), mad = mad(n), min(n), max(n))

df_recog_pre_e1 <- df_recog_e1 %>% 
  dplyr::filter(time == "pre") %>% 
  dplyr::select(subjID, TNTcond, itemID, correct_pre = correct)

df_rating_tntb_e1 <- df_rating_e1 %>% 
  dplyr::left_join(df_recog_pre_e1, by = c("subjID", "TNTcond", "itemID")) %>% 
  dplyr::filter(correct_pre == 1, TNTcond %in% c('t', 'nt', 'b')) %>% 
  dplyr::select(subjID, time, TNTcond, val, aro, itemID) %>% 
  tidyr::pivot_longer(cols = c(val, aro), names_to = "dim") %>% 
  tidyr::pivot_wider(names_from = time, values_from = value)

l_df_rating_subj_mean_e1 <- df_rating_tntb_e1 %>% 
  tidyr::drop_na() %>% 
  dplyr::mutate(diff = post - pre) %>% 
  dplyr::group_by(subjID, TNTcond, dim) %>% 
  dplyr::summarise(across(.cols = c(pre, post, diff), .fns = mean), .groups = "drop") %>% 
  split(.$dim)

# create summary description
## valence
df_val_e1_w <- l_df_rating_subj_mean_e1$val %>% 
  dplyr::select(subjID, TNTcond, diff) %>% 
  tidyr::pivot_wider(names_from = TNTcond, values_from = diff)
ttest_val_e1 <- t.test(df_val_e1_w$nt, df_val_e1_w$b, paired = T, alternative = 'greater')
smry_val_e1 <- get_t_g_ci(df_val_e1_w$nt, df_val_e1_w$b, method = 'greater')

## arousal
df_aro_e1_w <-  l_df_rating_subj_mean_e1$aro %>% 
  dplyr::select(subjID, TNTcond, diff) %>% 
  tidyr::pivot_wider(names_from = TNTcond, values_from = diff)
ttest_arol_e1 <- t.test(df_aro_e1_w$nt, df_aro_e1_w$b, paired = T, alternative = 'less')
smry_aro_e1 <- get_t_g_ci(df_aro_e1_w$nt, df_aro_e1_w$b, method = 'less')


# correlation -------------------------------------------------------------

df_indiv_ques_e1 <- 
  df_ques_e1 %>% 
  dplyr::arrange(subjID) %>%
  dplyr::mutate(
    # Positive value means beneficial consequence by memory contorl regarding the following variable 
    aro = df_aro_e1_w %>% arrange(subjID) %>% mutate(diff = b - nt) %>% pull(diff),
    val = df_val_e1_w %>% arrange(subjID) %>% mutate(diff = nt - b) %>% pull(diff)
  )

# saveRDS(df_indiv_ques_e1, here::here("analyses/rds/df_indiv_ques_e1.rds"))

# relationship between symptoms and intrusion proportion
# df_ques_e1 %>% 
#   dplyr::left_join(df_intr_TNT_e1 %>% dplyr::filter(TNTcond == 'nt') %>% dplyr::group_by(subjID) %>% dplyr::summarise(mean_int = mean(int_rate))) %>% 
#   dplyr::select(-c(subjID)) %>% 
#   psych::pairs.panels(lm = TRUE, ci = TRUE, ellipses = F)

# correlation (pairs panel)
df_indiv_ques_e1 %>%
  dplyr::select(-subjID) %>%
  psych::pairs.panels(lm = TRUE, ci = TRUE, ellipses = F)


# PCA
res_pca_e1 <- df_ques_e1 %>% 
  dplyr::select(BDI = bdi, RRS = rrs, `STAI-S` = stais, `STAI-T` = stait) %>% 
  stats::prcomp(scale = T)

smry_pca_e1 <- summary(res_pca_e1)

df_indiv_ques_all_e1 <- df_indiv_ques_e1 %>%
  dplyr::mutate(pc1 = res_pca_e1$x[,1], pc2 = res_pca_e1$x[,2])

# valence
lmfit_val_e1 <- lm(val ~ pc1 + pc2, data = df_indiv_ques_all_e1)
smry_lmfit_val_e1 <- summary(lmfit_val_e1)
r2_lmfit_val_e1 <- get_reg_stats(lmfit_val_e1)
cff_pc1_val_e1 <- get_cff_stats(lmfit_val_e1, "pc1")

# arousal
lmfit_aro_e1 <- lm(aro ~ pc1 + pc2, data = df_indiv_ques_all_e1)
smry_lmfit_aro_e1 <- summary(lmfit_aro_e1)
r2_lmfit_aro_e1 <- get_reg_stats(lmfit_aro_e1)
cff_pc1_aro_e1 <- get_cff_stats(lmfit_aro_e1, "pc1")

# for analyses on JASP
df_jasp_e1 <- l_df_rating_subj_mean_e1 %>% 
  dplyr::bind_rows(.id = "dim") %>% 
  dplyr::select(subjID, dim, TNTcond, diff) %>% 
  tidyr::pivot_wider(names_from = c(dim, TNTcond), names_prefix = "diff_", values_from = diff) %>% 
  dplyr::left_join(df_ques_e1, by = "subjID")

# write_csv(df_jasp_e1, here::here("analyses/jasp/df_jasp_e1.csv"))

# tables ------------------------------------------------------------------

df_recog_e1 %>%
  dplyr::select(subjID, time, TNTcond, correct) %>%
  dplyr::mutate(time = factor(time, levels = c("pre", "post")),
                TNTcond = if_else(stringr::str_detect(TNTcond, "new"), "new", TNTcond)) %>%
  dplyr::filter(TNTcond != "filler") %>%
  dplyr::group_by(time, TNTcond,  subjID) %>%
  dplyr::summarise(rate = mean(correct), .groups = "drop_last") %>%
  dplyr::summarise(mean = mean(rate), sd = sd(rate), .groups = "drop") %>%
  tidyr::pivot_longer(cols = mean:sd, names_to = "index") %>%
  tidyr::pivot_wider(names_from = TNTcond, values_from = "value") %>%
  dplyr::mutate(across(b:t, .fns = round, digits = 2)) %>% 
  dplyr::select(time, index, t, nt, b, new)

df_rating_tntb_e1 %>% 
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
  list(val = lm_to_table(lmfit_val_e1), aro = lm_to_table(lmfit_aro_e1)), 
  .id = "measure") %>% 
  knitr::kable(format = "pandoc")


# plots -------------------------------------------------------------------

gg_intr_TNT_e1 <- make_plot_intr(df_intr_TNT_e1) +
  theme(legend.position = c(0.8,0.6))

gg_hit_rate_e1 <- make_plot_hitrate(df_hit_rate_e1)

gg_aro_diff_e1 <- make_plot_rating(l_df_rating_subj_mean_e1$aro) + 
  labs(y = expression(paste(Delta, "Arousal (post-pre)")))

gg_val_diff_e1 <-  make_plot_rating(l_df_rating_subj_mean_e1$val) + 
  labs(y = expression(paste(Delta, "Valence (post-pre)")))

gg_pairplot_e1 <- make_pairplot(df_indiv_ques_e1)

upper_grid <- cowplot::plot_grid(gg_intr_TNT_e1, gg_hit_rate_e1, gg_val_diff_e1, gg_aro_diff_e1, 
                                 labels = c("A", "B", "C", "D"))

plot_all_e1 <- cowplot::plot_grid(upper_grid, ggmatrix_gtable(gg_pairplot_e1), ncol = 1, labels = c("", "E"))

# ggsave(here::here("analyses/figures/grid_plot_e1.svg"), plot_all_e1, width = 10.6, height = 11.4)

# gg_simple_slope_e1 <- 
#   sjPlot::get_model_data(
#     df_indiv_lm_e1[df_indiv_lm_e1$vars == 'aro',]$step2[[1]], 
#     type = 'pred', 
#     terms = c('pc1', 'pc2[-1, 1]')
#   ) %>%
#   ggplot() +
#   aes(x = x, y = predicted, linetype = group_col) +
#   geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.3) +
#   geom_line() +
#   labs(x = "PC 1", y = "Moderation in Arousal", linetype = "PC 2") +
#   theme_classic()

gg_res_pca_e1 <- make_pca_biplot(res_pca_e1)
# ggsave(filename = here::here("analyses/figures/biplot_pca_e1.png"), plot = gg_res_pca_e1, width = 6, height = 4)

gg_corr_intr_recog_e1 <- df_ques_e1 %>% 
  dplyr::left_join(df_hit_rate_e1_w %>% dplyr::mutate(`Item Recog` = b - nt, .keep = "unused"),
                   by = "subjID") %>% 
  dplyr::left_join(df_intr_TNT_e1 %>%
                     dplyr::filter(TNTcond == 'nt') %>%
                     dplyr::group_by(subjID) %>%
                     dplyr::summarise(mean_int = mean(int_rate)),
                   by = "subjID") %>%
  dplyr::select(BDI = bdi, 
                RRS = rrs, 
                `STAI-S` = stais, 
                `STAI-T` = stait, 
                `Intrusion` = mean_int,
                `Item Recog`) %>% 
  GGally::ggpairs(lower = list(continuous = wrap("smooth", alpha = 0.3)),
                  upper = list(continuous = wrap("cor", color = "black"))) +
  theme_bw(base_size = font_size, base_family = "Helvetica Neue") +
  theme(axis.text = element_text(colour = "black"),
        panel.grid = element_blank(),
        strip.background = element_blank())

# ggsave(filename = here::here("analyses/figures/pairs_plot_e1_supple.png"),
#        plot = gg_corr_intr_recog_e1, width = 8, height = 8)
