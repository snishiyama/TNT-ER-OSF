
# Setup -------------------------------------------------------------------

if (!exists('df_TNT_e2')){
  source(here::here('analyses/src/setup.R'))
}


# Count NAs ---------------------------------------------------------------

df_intr_NA_e2 <- df_TNT_e2 %>%
  dplyr::filter(is.na(intrusion)) %>%
  dplyr::count(subjID, block, intrusion) %>%
  dplyr::arrange(desc(n))


# Count intrusions --------------------------------------------------------

df_intr_TNT_e2 <-  
  df_TNT_e2 %>% 
  dplyr::filter(!is.na(intrusion)) %>% 
  dplyr::mutate(inORnot = if_else(intrusion==0, 0, 1)) %>% # 1 & 2 are coded 1 (i.e., intruded)
  dplyr::group_by(subjID,block,TNTcond) %>% 
  dplyr::summarise(int_rate = mean(inORnot))

# t test
ttest_intr_TNT_e2 <- 
  df_intr_TNT_e2 %>% 
  dplyr::filter(block %in% c(1, 10), TNTcond == 'nt') %>%
  t.test(int_rate ~ block, paired = T, data = ., alternative = 'greater')

df_intr_e2_w <- df_intr_TNT_e2 %>% 
  dplyr::filter(block %in% c(1, 10), TNTcond == 'nt') %>%
  tidyr::pivot_wider(names_from = block, values_from = int_rate, names_prefix = 'block')

smry_intr_e2 <- get_t_g_ci(df_intr_e2_w$block1, df_intr_e2_w$block10, method = 'greater')
