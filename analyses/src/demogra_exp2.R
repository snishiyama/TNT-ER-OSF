if (!exists('df_encoding_e2')){
  source(here::here('analyses/src/setup.R'))
}

# summarize demographic information
df_demogra_e2 <- 
  df_encoding_e2 %>% 
  dplyr::select(subjID, gender, age) %>% 
  dplyr::distinct_all() %>% 
  dplyr::summarise(male = sum(gender), Mage = mean(age), SDage = sd(age))
