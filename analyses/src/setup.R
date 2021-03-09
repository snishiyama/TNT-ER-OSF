
# set up ------------------------------------------------------------------
# To install "ggbiplot" package, install.pacakges() is NOT available. Instead,
# install.pacakges(devtools)
# library(devtools)
# install_github("vqv/ggbiplot")
# see https://github.com/vqv/ggbiplot

library(tidyverse)
library(here)
library(MBESS)
library(psych)
library(ggbiplot)
library(sjPlot)
library(pequod)
library(car)
library(GGally)
library(cowplot)
library(rstatix)
knitr::opts_chunk$set(echo = F)
options(readr.num_columns = 0) # readrのnotificationをoffにする

dir_data <- here::here('exp/data')

# common functions
load_data <- function(phase, expID, rm_subj = NULL, col_names = T, rm_col = NULL) {
  temp <- 
    fs::dir_ls(
      path = dir_data, 
      regexp = str_glue("[[:digit:]]+_{phase}_{expID}.csv$")
    ) %>% 
    purrr::map_dfr(
      .f = readr::read_csv, 
      col_names = col_names, 
      col_types = list("subjID" = col_integer())
    ) %>% 
    dplyr::filter(!subjID %in% rm_subj) %>% 
    dplyr::mutate(subjID = forcats::as_factor(subjID))
  if (!is.null(rm_col)) {
    return(temp %>% dplyr::select(-!!rm_col))
  } else {
    return(temp)
  }
}

calc_score <- function(df, quesName){
  df %>% 
    dplyr::group_by(subjID) %>% 
    dplyr::summarise(score = sum(value)) %>% # meanは合成得点を算出するため。最終的に標準化するので問題なし
    dplyr::mutate(ques = quesName)
}

# for reversed items
rev_value <- function(df){
  mutate(df, value = if_else(cond == "r", 4 - value, value))
}

# for conference poster (not used for figures in the article)
# my_theme_poster <- function(legend_pos_x, legend_pos_y) {
#   list(
#     theme_classic(base_size = bsize <- 42, base_family = 'Times New Roman'),
#     theme(
#       plot.background = element_rect(fill = '#D6D5D5', color = '#D6D5D5'),
#       axis.text = element_text(color = 'black'),
#       axis.title.x = element_blank(),
#       axis.text.x = element_text(size = bsize - 2),
#       legend.position = c(legend_pos_x, legend_pos_y),
#       legend.key.size = unit(40, units = 'points'),
#       legend.title = element_blank(),
#       legend.background = element_blank()
#     )
#   )
# }

format_pval <- function(pval) {
  if (pval < .001) {
    pval_new <- "_p_ < .001"
  }
  else {
    pval_new <- pval %>% 
      round(3) %>% 
      sprintf("%.3f", .) %>% 
      stringr::str_sub(start = 2) %>% 
      str_c("_p_ =", ., sep = " ")
  }
  return(pval_new)
}

format_pval_2 <- function(pval) {
  sapply(pval, function(x) {
    if (is.na(x)) {
      return(NA_character_)
    }
    if (x < .001) {
      pval_new <- "< .001"
    }
    else {
      pval_new <- x %>% 
        round(3) %>% 
        sprintf("%.3f", .) %>% 
        stringr::str_sub(start = 2)
    }
    return(pval_new)
  })
}

# function to calculate effect size from paired t test
# retrieved from https://tjo.hatenablog.com/entry/2014/02/24/192655
calc_g_ci <- function(x, y){
  n_x <- length(x)- 1
  n_y <- length(y)- 1
  mean_diff  <- mean(x) - mean(y) # 平均値の差
  csd <- n_x * var(x) + n_y * var(y)
  csd <- csd/(n_x + n_y)
  csd <- sqrt(csd)  # Common error varianceはこれで出せる
  cd  <- mean_diff/csd  # これでCohen’s dが求まる
  return(MBESS::ci.smd(smd = cd, n.1 = n_x + 1, n.2 = n_y + 1))
}

get_t_g_ci <- function(x, y, method = c("two.sided", "less", "greater")) {
  res <- t.test(x, y, paired = T, alternative = method)
  t <- res$statistic
  dof <- res$parameter
  p_value <- res$p.value
  g_ci <- calc_g_ci(x, y)
  g <- g_ci$smd
  ci_lower <- g_ci$Lower.Conf.Limit
  ci_upper <- g_ci$Upper.Conf.Limit
  
  sprintf("_t_(%d) = %.2f, %s, _d_ = %.2f, 95%% CI [%.2f, %.2f]",
          dof, t, format_pval(p_value), g, ci_lower, ci_upper)
}


# Output functions for regression analysis --------------------------------

# variables for initial draft
id_val <- 1
id_aro <- 2

ds_lpc2 <- 1
ds_hpc2 <- 3
ts_lpc2 <- 2
ts_hpc2 <- 4

# retrieved from https://terasawat.hatenablog.jp/entry/20140322/1395490704
lm.Beta <- function(res){
  b <- res$coefficients    #非標準化回帰係数
  ysd <- sd(res$model[,1]) #結果変数のSD
  idv <- res$model[,-1]
  if( class(idv) == "data.frame"){ 
    N <- ncol(idv) 
  }else{ 
    N <- 1
  }
  res.beta <- b[-1]  
  for(j in 1 : N ){
    if( N > 1){
      xxx <- idv[,j] 
    }else{
      xxx <- idv
    }
    if( class(xxx) == "numeric"){ #もし数値じゃなければ以下をかます
      # for(i in 2:length(levels(xxx)) ){
      #   dummy <- as.numeric(xxx == levels(xxx)[i])   # 1/0 の変数化
      #   lab <- paste(colnames(idv)[j],levels(xxx)[i],sep="") #対象の変数
      #   beta <-  b[lab] * ( sd(dummy) / ysd ) #ベータの計算方法：2変数のSDの商に回帰係数をかける
      #   res.beta[lab] <- beta 
      # }
      
      
    # }else{ #数値だったら簡単にβを計算できる
      lab <- paste(colnames(idv)[j],sep="") #対象の変数
      beta <- b[lab] * ( sd(xxx) / ysd )
      res.beta[lab] <- beta 
    }
  }
  res.beta
}

get_reg_stats <- function(lm_obj){
  res <- lm_obj %>% summary
  adj_r <- res$adj.r.squared
  dfa <- res$fstatistic[["numdf"]]
  dfe <- res$fstatistic[["dendf"]]
  f_val <- res$fstatistic[["value"]]
  p_val <- 1 - pf(f_val, dfa, dfe)
  
  sprintf("$R^2_{adj}$ = %.2f, _F_(%d, %d) = %.2f, %s",
          adj_r, dfa, dfe, f_val, format_pval(p_val))
}

get_cff_stats <- function(lm_obj, v_name){
  res <- lm_obj %>% summary
  beta <- res$coefficients[v_name, "Estimate"]
  p_val <- res$coefficients[v_name, "Pr(>|t|)"]
  se <- res$coefficients[v_name, "Std. Error"]
  b_std <- lm.Beta(lm_obj)[v_name]
  ci_l <- beta - 1.96 * se
  ci_u <- beta + 1.96 * se
  
  sprintf("_B_ = %.2f, 95%% CI [%.2f, %.2f], $\\beta$ = %.2f, %s",
          beta, ci_l, ci_u, b_std, format_pval(p_val))
}

get_aov_stats <- function(aov_obj, model_num){
  df1 <- aov_obj$Df[model_num]
  df2 <- aov_obj$Res.Df[model_num]
  f_val <- aov_obj$`F`[model_num]
  p_val <- aov_obj$`Pr(>F)`[model_num]
  
  sprintf("_F_(%d, %d) = %.2f, %s",
          df1, df2, f_val, format_pval(p_val))
}

get_ssa_stats <- function(simple_slope_obj, rownum){
  stats_ss <- simple_slope_obj$simple_slope
  slope <- stats_ss[rownum,'simple slope']
  se <- stats_ss[rownum,'standard error']
  ci_u <- slope + 1.96 * se
  ci_l <- slope - 1.96 * se
  p_val <- stats_ss[rownum,'p.value']
  sprintf("B = %.2f, 95%% CI [%.2f, %.2f], %s",
          slope, ci_l, ci_u, format_pval(p_val))
}

lm_to_table <- function(lm_obj) {
  l_ivs <- purrr::map(lm_obj, .f = ~ .x$coefficients %>% names())
  l_lm_Beta <- 
    lm_obj %>% 
    purrr::map(.f = lm.Beta) %>%
    purrr::map2(.y = l_ivs, 
                .f = ~ tibble(term = tail(.y, length(.y) - 1), beta = .x))

  l_tab_like_df <-
    lm_obj %>%
    purrr::map(.f = broom::tidy) %>%
    purrr::map(.f = mutate, ci = sprintf("[%.2f,%.2f]", estimate - 1.96 * std.error, estimate + 1.96 * std.error)) %>%
    purrr::map2(.y = l_lm_Beta, .f = ~left_join(x = .x, y = .y, by = c("term"))) %>%
    purrr::map(.f = select, term, estimate, ci, beta, p.value) %>%
    purrr::map(.f = mutate_at, .vars = c("estimate", "beta"), .funs = round, digits = 2) %>%
    purrr::map(.f = mutate_at, .vars = "p.value", .funs = format_pval_2)
  return(l_tab_like_df)
}

# import exp1 data --------------------------------------------------------

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

# import exp2 data -------------------------------------------------------------

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
