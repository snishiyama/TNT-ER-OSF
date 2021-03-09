
# Setup -------------------------------------------------------------------

if (!exists('df_rating_e1')){
  source(here::here('analyses/src/setup.R'))
}


# Summarise scores of items rated in both phases ------------------------

l_df_rating_e1 <- 
  df_rating_e1 %>% 
  dplyr::select(-order) %>% 
  tidyr::gather(key = dim, value = value, val, aro) %>% 
  split(.$dim) %>% 
  purrr::map(
    .f = ~ .x %>% 
      # get items rated in both phases
      dplyr::filter(TNTcond %in% c('b','nt')) %>% 
      tidyr::spread(key = time, value = value) %>% 
      tidyr::drop_na() %>% 
      dplyr::mutate(diff = post - pre)
  )

l_df_rating_subj_mean_e1 <- 
  l_df_rating_e1 %>% 
  purrr::map(
    .f = ~ .x %>%
      dplyr::group_by(subjID, TNTcond) %>% 
      dplyr::summarize_at(.vars = vars(pre, post, diff), .funs = mean)
  )

# t test for pre-post difference
l_ttest_rating_e1 <- 
  l_df_rating_subj_mean_e1 %>% 
  purrr::map2(.y = c('less', 'greater'),
    .f = ~ .x %>%
      dplyr::ungroup() %>% 
      dplyr::mutate(TNTcond = TNTcond %>% as_factor %>% fct_relevel('nt','b')) %>% 
      t.test(diff ~ TNTcond, paired = T, data = ., alternative = .y)
  )

# create summary description
# valence
df_val_e1_w <- l_df_rating_subj_mean_e1$val %>% select(subjID, TNTcond, diff) %>% spread(TNTcond, value = diff)
smry_val_e1 <- get_t_g_ci(df_val_e1_w$nt, df_val_e1_w$b, method = 'greater')
# arousal
df_aro_e1_w <- l_df_rating_subj_mean_e1$aro %>% select(subjID, TNTcond, diff) %>% spread(TNTcond, value = diff)
smry_aro_e1 <- get_t_g_ci(df_aro_e1_w$nt, df_aro_e1_w$b, method = 'less')

