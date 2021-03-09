
# Setup -------------------------------------------------------------------

if (!exists('df_rating_e2')){
  source(here::here('analyses/src/setup.R'))
}


# rating ------------------------------------------------------------------

l_df_rating_e2 <- 
  df_rating_e2 %>% 
  dplyr::select(-order) %>% 
  tidyr::gather(key = 'dim', value = 'rating', val, aro) %>% 
  base::split(.$dim) %>% 
  # handle both dimension data in the same way
  purrr::map(.f = ~ .x %>% 
    dplyr::filter(TNTcond %in% c('b','nt'), correct == 1) %>% 
    tidyr::spread(key = time, value = rating) %>% 
    tidyr::drop_na() %>% 
    dplyr::mutate(diff = post - pre)
  )

l_df_rating_subj_mean_e2 <- 
  l_df_rating_e2 %>% 
  purrr::map(
    .f = ~ .x %>% 
      dplyr::group_by(subjID, TNTcond) %>% 
      dplyr::summarize_at(.vars = vars(pre, post, diff), .funs = mean)
  )

# t test for pre-post difference
l_ttest_rating_e2 <- 
  l_df_rating_subj_mean_e2 %>% 
  purrr::map2(.y = c('less', 'greater'), # less for arousal, greater for valence
    .f = ~ .x %>%
      dplyr::ungroup() %>% 
      dplyr::mutate(TNTcond = TNTcond %>% as_factor %>% fct_relevel('nt','b')) %>% 
      t.test(diff ~ TNTcond, paired = T, data = ., alternative = .y)
  )

# create summary description
# valence
df_val_e2_w <- l_df_rating_subj_mean_e2$val %>% select(subjID, TNTcond, diff) %>% spread(TNTcond, value = diff)
smry_val_e2 <- get_t_g_ci(df_val_e2_w$nt, df_val_e2_w$b, method = 'greater')
# arousal
df_aro_e2_w <- l_df_rating_subj_mean_e2$aro %>% select(subjID, TNTcond, diff) %>% spread(TNTcond, value = diff)
smry_aro_e2 <- get_t_g_ci(df_aro_e2_w$nt, df_aro_e2_w$b, method = 'less')


# for poster presentation -------------------------------------------------

# l_gg_rating_mean_e2 <- 
#   l_df_rating_subj_mean_e2 %>% 
#   purrr::map(
#     .f = ~ .x %>% 
#       dplyr::filter(TNTcond != 't') %>% 
#       dplyr::select(-diff) %>% 
#       tidyr::gather(key = time, value = mean_rating, pre, post) %>% 
#       dplyr::mutate(TNTcond = if_else(TNTcond == 'b', 'Baseline', 'No-Think'),
#                     time = as_factor(time)) %>% 
#       ggplot2::ggplot(aes(x = time, y = mean_rating, fill = time)) +
#       geom_dotplot(binwidth = 0.1, binaxis = 'y', stackdir = 'center',
#                    alpha = 0.5, dotsize = 1, binpositions = 'all', position = position_dodge()) +
#       stat_summary(fun.data = 'mean_cl_normal', position = position_dodge(width = 0.9), 
#                    show.legend = F, size = 2) +
#       scale_fill_grey(start = 0.5, end = 1) +
#       facet_grid(cols = vars(TNTcond), switch = 'x') +
#       labs(y = 'Mean Valence Rating') +
#       theme_classic(base_size = 16, base_family = 'Times New Roman') +
#       my_theme_poster(.165,.96) +
#       theme(
#         legend.position = 'none',
#         strip.background = element_rect(fill = '#D6D5D5', color = '#D6D5D5', size = 0.1),
#         strip.text = element_text(size = 42),
#         strip.placement = 'outside',
#         panel.spacing.x = unit(0, units = 'points')
#       )
#   )
# 
# l_gg_rating_diff_e2 <- 
#   l_df_rating_subj_mean_e2 %>% 
#   purrr::map(
#     .f = ~ .x %>% 
#       dplyr::filter(TNTcond != 't') %>% 
#       dplyr::mutate(TNTcond = if_else(TNTcond == 'b', 'Baseline', 'No-Think')) %>% 
#       ggplot2::ggplot(aes(x = TNTcond, y = diff)) +
#       geom_hline(yintercept = 0, linetype = 'longdash', size = 1) +
#       geom_dotplot(binwidth = 0.1, binaxis = 'y', stackdir = 'center', aes(fill = TNTcond), show.legend = F,
#                    alpha = 0.5, dotsize = 1, binpositions = 'all', position = position_dodge()) +
#       stat_summary(fun.data = 'mean_cl_normal', position = position_dodge(width = 0.9), 
#                    show.legend = F, size = 2) +
#       scale_fill_manual(values = c('Baseline' = 'blue', 'No-Think' = 'red')) +
#       theme_classic(base_size = 16, base_family = 'Times New Roman') +
#       my_theme_poster(.165,.96)
#   )
# 
# l_gg_rating_diff_e2$val <- 
#   l_gg_rating_diff_e2$val +
#   labs(y = 'Mean Valence Difference')
# 
# l_gg_rating_diff_e2$aro <- 
#   l_gg_rating_diff_e2$aro +
#   labs(y = 'Mean Arousal Difference')

# ggsave(filename = here::here('analyses/figures/mean_valence-diff_exp2.png'),
#        plot = l_gg_rating_diff_e2$val, width = 30, height = 20, units = 'cm')
# 
# ggsave(filename = here::here('analyses/figures/mean_arousal-diff_exp2.png'),
#        plot = l_gg_rating_diff_e2$aro, width = 30, height = 20, units = 'cm')
