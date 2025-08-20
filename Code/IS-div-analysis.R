#-------------- modeling diversity ~ interaction strength (IS)-----------------#
#---------------- analyze results of pairwise non-linear model-----------------#

# platform  "x86_64-w64-mingw32"
# R version 4.3.3 (2024-02-29 ucrt) "Angel Food Cake"

rm(list=ls())

# load required packages
library(ggplot2)
library(tidyverse)
library(scales)
library(cowplot) # for get_legend
library(gridExtra)


# load required data from Yu et al 2024
load("stan_dat_1.RData")
alpha_matrix <- stan_dat$alpha

# load estimated parameter output form pairwise non-linear model (model.10)
sum_10 <- read.csv("pars_pairwise_nl.csv", row.names = 1)

# extract net HOI effects gamma
par_gamma_10 <- sum_10[grep("^gamma", rownames(sum_10)), "mean"]
gamma_mat_10 <- matrix(par_gamma_10, nrow = 8, ncol = 8, byrow = T)

# create index matrix
index_matrix <- matrix(1:64, nrow = 8, ncol = 8, byrow = TRUE)

# extract the gamma we actually have data to estimate (37 in total)
index_gamma <- sort(index_matrix[stan_dat$div_idx])  

# prepare data.frame for later plotting
df_gamma <- data.frame(alpha=c(t(alpha_matrix)), gamma=par_gamma_10)
# interaction classification: pos, neg, intra
group_id <- sapply(df_gamma$alpha, function(x)if (x>0) {return ("pos")} else {return ("neg")})
group_id[diag(index_matrix)] <- "intra"
df_gamma$group <- group_id

# calculate pairwise interactions across diversity level based on net HOI effects
div <- c(1,2,4,8)
gamma_mat <- gamma_mat_10

res_mat <- div |> 
  lapply(function(x) stan_dat$alpha + gamma_mat*log2(x/stan_dat$div_ave)) |> # log (non_linear)
  lapply(function(x) as.vector(t(x))) |>
  (function(x) do.call(cbind, x))() |> 
  (function(x) as.data.frame(x))() 

colnames(res_mat) <- c("mono", "two", "four", "eight")
res_mat$group <- df_gamma$group


# create null_mat for simulation
# intra-specific interactions from mono-culture 
# and inter-specific interactions from 2 species mixture (replicate it across div)
null_mat <- res_mat[ ,c("mono","two")]
null_mat$two[res_mat$group == "intra"] <- null_mat$mono[res_mat$group == "intra"]
null_mat$four <- null_mat$two
null_mat$eight <- null_mat$two
null_mat$group <- res_mat$group

group_id_null <- sapply(null_mat$two,function(x)if (x>0) {return ("pos")} else {return ("neg")})
group_id_null[diag(index_matrix)] <- "intra"
df_gamma$group_null <- group_id_null
res_mat$group_null <- group_id_null
null_mat$group_null <- group_id_null

# extract the gammas for which we have info to estimate
df_gamma <- df_gamma[index_gamma, ] 
res_mat_real <- res_mat[index_gamma, ] 


# make histogram of gamma
p1 <- ggplot(df_gamma, aes(x=gamma, fill=group_null)) +
  geom_vline(xintercept = 0, col="grey", linetype = "longdash", linewidth = 1) +
  geom_density( color="#e9ecef", alpha=0.6, position = "identity") + # identity is for the overlapping effect
  xlim(-1.5, 2.5) +
  scale_fill_manual(values=c("#A6A15E", "#B6495C",  "#6D9DC5"),
                    labels=c("Intra-specific", "Negative inter-specific", "Positive inter-specific")) +
  labs(fill="", title = "Distribution of Net HOI Effects",
       x="Net HOI effects", y="Density", tag = "A" )+  #expression(gamma[s(i)] [s(j)])
  theme_classic() +
  theme(legend.position = "none",  #"inside"
        plot.title = element_text(hjust = 0.5, size = 9),
        plot.tag = element_text(face = "bold", size = 10),
        plot.tag.position = c(0.05, 1),
        axis.title.x = element_text(size = 9),
        axis.title.y = element_text(size = 9),
        axis.text.x = element_text(size=8),
        axis.text.y = element_text(size=8)  ) 


# covert to long data format for plotting
df_IS_div <- res_mat_real %>%
  rownames_to_column("row_id") %>%
  pivot_longer(
    cols = c(mono, two, four, eight),  
    names_to = "div",
    values_to = "IS"
  ) %>%
  mutate(
    div = case_when(
      div == "mono"   ~ 1,
      div == "two"    ~ 2,
      div == "four"   ~ 4,
      div == "eight"  ~ 8,
      TRUE ~ NA_real_
    )
  ) %>%
  arrange(row_id, div)


# plot the div~IS relationship for pairwise models
p2 <- ggplot(df_IS_div, aes(x=log2(div), y=IS, colour=group_null)) +
  #geom_point(alpha=0.85, size=4, shape = 19, stroke = 0.5 )+
  geom_line( aes(group = row_id), colour="#F0EAD6", linewidth = 2.4, linetype = "solid" ) +
  geom_line( aes(group = row_id), alpha = 0.7, linewidth = 1.4, linetype = "solid" ) +
  scale_color_manual(values=c("#A6A15E", "#B6495C",  "#6D9DC5"),
                     labels=c("Intra-specific", "Negative inter-specific", "Positive inter-specific")) +
  labs(x = expression(log[2]("div")), y = "Pairwise Interaction",
       title = "Pairwise Interaction across Diversity", tag = "B" ) +
  theme_classic()+
  theme(legend.position = "none", # remove legend
        plot.title = element_text(hjust = 0.5, size = 9),
        plot.tag = element_text(face = "bold", size = 10),
        plot.tag.position = c(0.05, 1),
        axis.title.x = element_text(size = 9),
        axis.title.y = element_text(size = 9),
        axis.text.x = element_text(size=8),
        axis.text.y = element_text(size=8))


# get common legend for both plots
legend <- get_legend( p1 + 
                        theme(legend.position = "right",# defaut legend position
                              legend.text = element_text(size = 8),
                              legend.margin = margin(t = -5, unit = "mm")) + 
                        guides(fill = guide_legend(nrow = 1)) )  



# pool each IS by group together and fit a linear model through each group
p3 <- ggplot(df_IS_div, aes(x=log2(div), y=IS, colour=group_null)) +
  geom_jitter(data=filter(df_IS_div, group_null=="intra"), aes(x=log2(div), y=IS), colour="#A6A15E", width = 0.3, alpha = 0.5, size = 3) +
  geom_jitter(data=filter(df_IS_div, group_null=="neg" & div!=1), aes(x=log2(div), y=IS), colour="#B6495C", width = 0.3, alpha = 0.5, size = 3) +
  geom_jitter(data=filter(df_IS_div, group_null=="pos" & div!=1 ), aes(x=log2(div), y=IS), colour = "#6D9DC5", width = 0.3, alpha = 0.5, size = 3) +
  geom_hline(yintercept = 0, colour="darkgrey") +
  geom_smooth( aes(group = group_null), method = "lm",se = F,  colour="#F0EAD6", linewidth = 2.4, linetype = "solid" ) +
  geom_smooth(method = "lm", linewidth = 1.4, alpha=0.2)+
  scale_color_manual(values=c("#A6A15E", "#B6495C",  "#6D9DC5"),
                     labels=c("Intra-specific", "Negative Inter-specific", "Positive Inter-specific")) +
  labs(x = expression(log[2]("div")), y = "Pairwise Interaction",
       title = "Averaged Trend across diversity", tag = "C" ) +
  theme_classic() +
  theme(legend.position = "none",
        #legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 9),
        plot.tag = element_text(face = "bold", size = 10),
        plot.tag.position = c(0.05, 1),
        axis.title.x = element_text(size = 9),
        axis.title.y = element_text(size = 9),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8))


# show the slopes for each group
summary(lm(IS ~ log2(div), data = filter(df_IS_div, group_null=="intra")))
summary(lm(IS ~ log2(div), data = filter(df_IS_div, group_null=="pos" & div !=1)))
summary(lm(IS ~ log2(div), data = filter(df_IS_div, group_null=="neg" & div !=1)))


fig_2 <-grid.arrange(arrangeGrob(p1, p2, p3, ncol = 3), 
                     legend,  nrow = 2,  
                     heights = c(8, 1)) 


ggsave(filename = "fig_2.png",  plot = fig_2, device = "png",  
       width = 180,  height = 75,  units = "mm",  dpi = 600,  
       bg = "white" )

