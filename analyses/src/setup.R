
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
    dplyr::summarise(score = sum(value), .groups = "drop") %>% # meanは合成得点を算出するため。最終的に標準化するので問題なし
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

get_t_g_ci <- function(x, y, method = c("two.sided", "less", "greater")) {
  res <- t.test(x, y, paired = T, alternative = method)
  t <- res$statistic
  dof <- res$parameter
  p_value <- res$p.value
  dz_ci <- MBESS::ci.sm(ncp = t, N = dof + 1)
  dz <- dz_ci$Standardized.Mean
  ci_lower <- dz_ci$Lower.Conf.Limit.Standardized.Mean # g_ci$Lower.Conf.Limit
  ci_upper <- dz_ci$Upper.Conf.Limit.Standardized.Mean# g_ci$Upper.Conf.Limit
  
  
  sprintf("_t_(%d) = %.2f, %s, _dz_ = %.2f, 95%% CI [%.2f, %.2f]",
          dof, t, format_pval(p_value), dz, ci_lower, ci_upper)
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


# plot functions ----------------------------------------------------------

font_size <- 12

make_plot_intr <- function(df_intr) {
  df_intr %>% 
    dplyr::mutate(TNTcond = forcats::fct_relevel(TNTcond, "t","nt")) %>% 
    ggplot2::ggplot()+
    aes(x=block, y=int_rate, group=TNTcond) +
    stat_summary(geom = "line", fun.y = "mean", aes(linetype = TNTcond)) +
    stat_summary(geom = "errorbar", fun.data = "mean_se", width=0.3, color="black") +
    stat_summary(geom = "point", fun.y = "mean", size=3, mapping = aes(shape = TNTcond)) +
    scale_shape_manual(values = c("nt"=16, "t"=17), labels = c("nt"="No-Think", "t"="Think"), name="TNT status") +
    scale_linetype_manual(values = c("nt"="dashed", "t"="solid"), labels = c("nt"="No-Think", "t"="Think"), name="TNT status") +
    scale_x_continuous(breaks=seq(1,10,by=1)) +
    scale_y_continuous(breaks=seq(0,1,by=0.1), limits=c(0,1)) +
    labs(x="Block", y="Intrusion Proportion") +
    theme_bw(base_size = font_size, base_family = "Helvetica Neue") +
    theme(panel.grid = element_blank(),
          axis.text = element_text(color = "black"),
          legend.title = element_blank(),
          legend.background = element_blank())
}

make_plot_intr_supp <- function(df_intr_woex_1, df_intr_woex_2) {
  list(`Experiment 1` = df_intr_woex_1, `Experiment 2` = df_intr_woex_2) %>% 
    dplyr::bind_rows(.id = "expt") %>% 
    dplyr::filter(TNTcond == 'nt') %>% 
    dplyr::mutate(block = as.factor(block),
                  correct = factor(as.character(correct), levels = c('1', '0'), labels = c("Recognized", "Non-Recognized"))) %>% 
    dplyr::group_by(subjID, block, TNTcond, correct, expt) %>% 
    dplyr::summarise(int_rate = mean(inORnot), .groups = "drop") %>% 
    ggplot(aes(x = block, y = int_rate, group = correct)) +
    stat_summary(geom = "line", fun.y = "mean", aes(linetype = correct)) +
    stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0.3) +
    stat_summary(geom = "point", fun.y = "mean", size = 3, aes(shape = correct)) +
    scale_shape_manual(values = c("Recognized" = 16, "Non-Recognized" = 4)) +
    labs(x="Block", y="Intrusion Proportion") +
    theme_bw(base_size = font_size, base_family = "Helvetica Neue") +
    theme(panel.grid = element_blank(),
          axis.text = element_text(color = "black"),
          legend.title = element_blank(),
          legend.background = element_blank(),
          legend.position = c(0.7, 0.2)) +
    facet_grid(cols = vars(expt))
}

make_plot_hitrate <- function(df_hitrate) {
  df_hitrate %>% 
    dplyr::ungroup() %>% 
    dplyr::mutate(TNTcond = if_else(TNTcond == "b", "Baseline", "No-Think")) %>% 
    ggplot2::ggplot() +
    aes(x = TNTcond, y = hit_rate) +
    stat_summary(geom = "errorbar", fun.data = 'mean_se', position = position_dodge(width = 0.9), 
                 show.legend = F, width = 0.07) +
    stat_summary(geom = "point", fun = "mean", position = position_dodge(width = 0.9), 
                 show.legend = F, fill = "grey95", color = "black", shape = 16, size = 3) +
    geom_jitter(width = 0.2, color = "black", alpha = 0.3, height = 0, size = 2) +
    labs(y = "Hit rate") +
    theme_bw(base_size = font_size, base_family = "Helvetica Neue") +
    theme(axis.text = element_text(color = "black"),
          axis.title.x = element_blank(),
          panel.grid = element_blank(),
          legend.position = "none")
}

make_plot_rating <- function(df_rating) {
  df_rating %>% 
    dplyr::mutate(TNTcond = if_else(TNTcond == 'b', 'Baseline', 'No-Think')) %>% 
    ggplot2::ggplot(aes(x = TNTcond, y = diff)) +
    geom_hline(yintercept = 0, linetype = 'longdash') +
    stat_summary(geom = "errorbar", fun.data = 'mean_se', position = position_dodge(width = 0.9), 
                 show.legend = F, width = 0.07)+
    stat_summary(geom = "point", fun = "mean", position = position_dodge(width = 0.9), 
                 show.legend = F, fill = "grey95", color = "black", shape = 16, size = 3) +
    geom_jitter(width = 0.2, color = "black", alpha = 0.3, height = 0, size = 2) +
    theme_bw(base_size = font_size, base_family = 'Helvetica Neue') +
    theme(
      axis.text = element_text(color = 'black'),
      axis.title.x = element_blank(),
      panel.grid = element_blank(),
      legend.title = element_blank(),
      legend.background = element_blank()
    )
}

make_pca_biplot <- function(pca_obj) {
  ggbiplot(
    pca_obj,
    scale = 1,
    ellipse = TRUE, 
    circle = F,
    alpha = 0.5
  )  +
    coord_cartesian(xlim = c(-2, 2)) +
    theme_bw(base_size = font_size, base_family = "Helvetica Neue") +
    theme(panel.grid = element_blank(),
          axis.text = element_text(color = "black"),
          plot.margin = margin(t = 10, b = 10, r = 10, l = 10))
}

make_pairplot <- function(df_indiv_ques) {
  df_indiv_ques %>% 
    dplyr::select(BDI = bdi, 
                  RRS = rrs, 
                  `STAI-S` = stais, 
                  `STAI-T` = stait, 
                  Valence = val,
                  Arosal = aro) %>%
    GGally::ggpairs(lower = list(continuous = wrap("smooth", alpha = 0.3)),
                    upper = list(continuous = wrap("cor", color = "black"))) +
    theme_bw(base_size = font_size, base_family = "Helvetica Neue") +
    theme(axis.text = element_text(colour = "black"),
          panel.grid = element_blank(),
          strip.background = element_blank())
}
