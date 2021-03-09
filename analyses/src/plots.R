
# setup -------------------------------------------------------------------

if (!exists("df_indiv_lm")) {
  source(here::here("analyses/src/compare_exp1-2.R"))
}

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

# Experiment 1 ------------------------------------------------------------

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

ggsave(here::here("analyses/figures/grid_plot_e1.png"), plot_all_e1, width = 10.6, height = 11.4)

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
ggsave(filename = here::here("analyses/figures/biplot_pca_e1.png"), plot = gg_res_pca_e1, width = 6, height = 4)

# Experiment 2 ------------------------------------------------------------

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

ggsave( here::here("analyses/figures/grid_plot_e2.png"), plot_all_e2, width = 10.6, height = 11.4)

gg_res_pca_e2 <- make_pca_biplot(res_pca_e2)
ggsave(filename = here::here("analyses/figures/biplot_pca_e2.png"), plot = gg_res_pca_e2, width = 6, height = 4)

# Integrated Analyses -----------------------------------------------------

gg_res_pca_both <- make_pca_biplot(res_pca_ques)
ggsave(filename = here::here("analyses/figures/biplot_pca_both.png"), plot = gg_res_pca_both, width = 6, height = 4)
