
# Setup -------------------------------------------------------------------

if (!exists('df_recog_e2')){
  source(here::here('analyses/src/setup.R'))
}


# Hit rate ----------------------------------------------------------------

df_hit_rate_e2 <- 
  df_recog_e2 %>% 
  dplyr::select(subjID, time, TNTcond, itemID, correct) %>% 
  tidyr::spread(key = time, value = correct) %>% 
  dplyr::filter(pre == 1, TNTcond %in% c('nt', 'b')) %>% 
  dplyr::group_by(TNTcond, subjID) %>% 
  dplyr::summarise(hit_rate = mean(post))

# t test
ttest_hit_rate_e2 <- 
  df_hit_rate_e2 %>%
  dplyr::ungroup() %>% 
  dplyr::mutate(TNTcond = TNTcond %>% fct_drop() %>% fct_relevel('nt', 'b')) %>% 
  t.test(hit_rate ~ TNTcond, paired = T, data = ., alternative = 'less')

# create summary description
df_hit_rate_e2_w <- df_hit_rate_e2 %>% spread(TNTcond, value = hit_rate)
smry_hit_rate_e2 <- get_t_g_ci(df_hit_rate_e2_w$nt, df_hit_rate_e2_w$b, method = 'less')

# Number of dropped items
df_hit_dropped_e2 <- 
  df_recog_e2 %>% 
  dplyr::select(subjID, time, TNTcond, itemID, correct) %>% 
  tidyr::spread(key = time, value = correct) %>% 
  dplyr::filter(pre == 0, TNTcond %in% c('nt','b')) %>% 
  dplyr::group_by(TNTcond, subjID) %>% 
  dplyr::tally(name = 'dropped_n') %>% 
  dplyr::summarise(mean = mean(dropped_n), sd = sd(dropped_n))


# Substitutes -------------------------------------------------------------

df_hit_sub_e2 <- 
  df_recog_e2 %>% 
  dplyr::filter(TNTcond == 'sub') %>% 
  dplyr::select(-order, -TNTcond) %>% 
  dplyr::left_join(
    df_recog_e2 %>% dplyr::filter(TNTcond != 'sub') %>% select(TNTcond, subjID, itemID),
    by = c('subjID', 'itemID')
  ) %>% 
  dplyr::group_by(TNTcond, subjID) %>% 
  dplyr::summarise(hit_rate = mean(correct))

# for Table 
# df_hit_sub_e2 %>%
#   dplyr::group_by(TNTcond) %>%
#   dplyr::summarise(mean = mean(hit_rate), sd = sd(hit_rate)) %>%
#   dplyr::mutate_if(.predicate = is_double, round, digits = 2)
  
ttest_hit_sub_e2 <-  
  df_hit_sub_e2 %>%
  dplyr::ungroup() %>% 
  dplyr::mutate(TNTcond = TNTcond %>% fct_drop() %>% fct_relevel('nt', 't')) %>% 
  t.test(hit_rate ~ TNTcond, paired = T, data = ., alternative = 'greater')

# create summary description
df_hit_sub_e2_w <- df_hit_sub_e2 %>% spread(TNTcond, value = hit_rate)
smry_hit_sub_e2 <- get_t_g_ci(df_hit_sub_e2_w$nt, df_hit_sub_e2_w$t, method = 'greater')



# forced-choice recognition -----------------------------------------------

df_recall_rate_e2 <- 
  df_rating_e2 %>% 
  dplyr::select(TNTcond, itemID, subjID, correct, time) %>% 
  tidyr::spread(key = time, value = correct) %>% 
  dplyr::filter(pre == 1, TNTcond %in% c('b', 'nt')) %>% 
  dplyr::group_by(subjID, TNTcond) %>% 
  dplyr::summarise(recall_rate = mean(post))

ttest_recall_rate_e2 <- 
  df_recall_rate_e2 %>% 
  dplyr::mutate(TNTcond = TNTcond %>% as_factor %>% fct_relevel('nt','b')) %>% 
  t.test(recall_rate ~ TNTcond, paired = T, data = ., alternative = 'less')

df_recall_e2_w <- df_recall_rate_e2 %>% spread(TNTcond, value = recall_rate)
smry_recall_e2 <- get_t_g_ci(df_recall_e2_w$nt, df_recall_e2_w$b, method = 'less')


# Tables ------------------------------------------------------------------

# old/new recognition
# df_recog_e2 %>%
#   dplyr::select(subjID, time, TNTcond, correct) %>%
#   dplyr::mutate(time = factor(time, levels = c("pre", "post"))) %>%
#   dplyr::group_by(subjID, time, TNTcond) %>%
#   dplyr::summarise(rate = mean(correct)) %>%
#   dplyr::group_by(time, TNTcond) %>%
#   dplyr::summarise(mean = mean(rate), sd = sd(rate)) %>%
#   dplyr::filter(TNTcond != "filler") %>%
#   tidyr::pivot_longer(cols = mean:sd, names_to = "index") %>%
#   tidyr::pivot_wider(names_from = TNTcond, values_from = "value") %>%
#   dplyr::mutate_if(.predicate = is.double, .funs = round, digits = 2) %>%
#   dplyr::select(time, index, t, nt, b, sub, new)

# forced-choice recognition
# df_rating_e2 %>% 
#   dplyr::mutate(time = factor(time, levels = c("pre", "post"))) %>%
#   dplyr::group_by(subjID, time, TNTcond) %>%
#   dplyr::summarise(rate = mean(correct)) %>%
#   dplyr::group_by(time, TNTcond) %>%
#   dplyr::summarise(mean = mean(rate), sd = sd(rate)) %>%
#   dplyr::filter(TNTcond != "filler") %>%
#   tidyr::pivot_longer(cols = mean:sd, names_to = "index") %>%
#   tidyr::pivot_wider(names_from = TNTcond, values_from = "value") %>%
#   dplyr::mutate_if(.predicate = is.double, .funs = round, digits = 2)

# Visualize
gg_hit_rate_e2 <- 
  df_hit_rate_e2 %>% 
  ggplot2::ggplot(aes(x = TNTcond, y = hit_rate)) +
  stat_summary(geom = 'pointrange', fun.data = 'mean_cl_normal', position = position_dodge(width = 0.95)) +
  geom_jitter(aes(color = TNTcond), width = 0.3) +
  theme_classic()

