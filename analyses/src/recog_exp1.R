
# Setup -------------------------------------------------------------------

if (!exists('df_recog_e1')){
  source(here::here('analyses/src/setup.R'))
}


# Hit rate ----------------------------------------------------------------

# We analyzed hit rates of items that were recognized in the pre-TNT recognition test.
df_hit_rate_e1 <- 
  df_recog_e1 %>% 
  dplyr::select(subjID, time, TNTcond, itemID, correct) %>% 
  tidyr::spread(key = time, value = correct) %>% 
  dplyr::filter(pre == 1, TNTcond %in% c('nt', 'b')) %>% 
  dplyr::group_by(TNTcond, subjID) %>% 
  dplyr::summarise(hit_rate = mean(post))

# t test
ttest_hit_rate_e1 <- 
  df_hit_rate_e1 %>%
  dplyr::ungroup() %>% 
  dplyr::mutate(TNTcond = TNTcond %>% fct_drop() %>% fct_relevel('nt', 'b')) %>% 
  t.test(hit_rate ~ TNTcond, paired = T, data = ., alternative = 'less')

# create summary description
df_hit_rate_e1_w <- df_hit_rate_e1 %>% spread(TNTcond, value = hit_rate)
smry_hit_rate_e1 <- get_t_g_ci(df_hit_rate_e1_w$nt, df_hit_rate_e1_w$b, method = 'less')

# Number of dropped items
df_hit_dropped_e1 <- 
  df_recog_e1 %>% 
  dplyr::select(subjID, time, TNTcond, itemID, correct) %>% 
  tidyr::spread(key = time, value = correct) %>% 
  dplyr::filter(pre == 0, TNTcond %in% c('nt','b')) %>% 
  dplyr::group_by(TNTcond, subjID) %>% 
  dplyr::tally(name = 'dropped_n') %>% 
  dplyr::summarise(mean = mean(dropped_n), sd = sd(dropped_n))


# Tables ------------------------------------------------------------------

# table
# df_recog_e1 %>%
#   dplyr::select(subjID, time, TNTcond, correct) %>%
#   dplyr::mutate(time = factor(time, levels = c("pre", "post"))) %>%
#   dplyr::group_by(subjID, time, TNTcond) %>%
#   dplyr::summarise(rate = mean(correct)) %>%
#   dplyr::mutate(TNTcond = if_else(stringr::str_detect(TNTcond, "new"), "new", TNTcond)) %>%
#   dplyr::group_by(time, TNTcond) %>%
#   dplyr::summarise(mean = mean(rate), sd = sd(rate)) %>%
#   dplyr::filter(TNTcond != "filler") %>%
#   tidyr::pivot_longer(cols = mean:sd, names_to = "index") %>%
#   tidyr::pivot_wider(names_from = TNTcond, values_from = "value") %>%
#   dplyr::mutate_at(.vars = vars(b:t), .funs = round, digits = 2) %>% 
#   dplyr::select(time, index, t, nt, b, new)


# df_recog_e1 %>% 
#   dplyr::select(subjID, time, TNTcond, itemID, correct) %>% 
#   dplyr::filter(TNTcond != "filler") %>% 
#   dplyr::mutate(TNTcond = if_else(stringr::str_detect(TNTcond, "new"), "new", TNTcond)) %>% 
#   dplyr::group_by(subjID, time, TNTcond) %>% 
#   dplyr::summarise(corr = mean(correct)) %>% 
#   ggplot2::ggplot(aes(x = TNTcond, y = corr)) +
#   stat_summary(geom = 'pointrange', fun.data = 'mean_cl_normal', position = position_dodge(width = 0.95)) +
#   geom_jitter(aes(color = TNTcond), width = 0.3) +
#   theme_classic() +
#   facet_wrap(~time)