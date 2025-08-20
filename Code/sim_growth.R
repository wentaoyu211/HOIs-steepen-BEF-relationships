#------------------------------------------------------------------------------#
#----------simulating tree growth over years and calculate biomass-------------#
#--------------- use random planting schemes to create replicates -------------#
#------------------------------------------------------------------------------#

# data required:
#                res_mat: interaction matrix with net HOI effects
#                null_mat: interaction matrix without net HOI effects
# IS-div-analysis.R need to be run first to have res_mat and null_mat available


# platform  "x86_64-w64-mingw32"
# R version 4.3.3 (2024-02-29 ucrt) "Angel Food Cake"


# load required packages
library(dplyr)
library(ggplot2)
library(cowplot)

# load the empirical data 
load("stan_dat_1.RData")

# define required parameters for simulation
b <- stan_dat$b
theta <- stan_dat$theta
beta <- stan_dat$beta
alpha_matrix <- stan_dat$alpha
year <- 7


# create matrices store the species ID (compositions) across diversity gradient 
# this is the real composition and need to be maintained due to broken-stick design
sp_comp <- list(1,2,3,4,5,6,7,8,c(1,5), c(2,6), c(3,8), c(4,7), c(1,4,5,7), c(2,3,6,8), c(1:8))
# number of replicates of same plot composition (only affected by starting biomass)
rep <- 3
# number of trees per plot
nind <- 6*6

# in BEF China 4 and 8 mixture is 12 by 12 instead of 6 by 6 in mono and 2 mixture
# here I quadruple the replicates in these mixtures
sp_comp_rep <- unlist(
  lapply(seq_along(sp_comp), function(i) rep(list(sp_comp[[i]]), ifelse(i <= 12, rep, rep*4))),
  recursive = FALSE
)


# load required functions for simulation
source("functions.R")

# run simulations in parallel
library(parallel)
library(doParallel)
print(detectCores())
cl <- makeCluster(8)
print(cl)  # check if it's a valid connection
registerDoParallel(cl)

# number of replicated BEF experiments 
set.seed(100)   # 25.8
n_sim = 1000 # simulation take a bit time

system.time(
  output_list <- foreach(i=1:n_sim) %dopar% {
    output <- sim_bef(sp_comp_rep, nind, res_mat[,1:4], null_mat[,1:4])
    output$BEF_id <- i
  })

stopCluster(cl)

output_list <- vector("list", length = n_sim)
system.time(
  for (i in 1:n_sim){
    output <- sim_bef(sp_comp_rep, nind, res_mat[,1:4], null_mat[,1:4])
    output$BEF_id <- i
    output_list[[i]] <- output
  }
)

output_df <- do.call(rbind, output_list)
write.csv(output_df, file = "sim_output.csv") 



#----------------- analyze and visualize simulation results -------------------#
# if it's two long, the simulation output can be loaded directly
# output_df <- read.csv("sim_output.csv")

# linear regression of BEF relationship
model_HOI <- lm(std_y7 ~ log2(div), data = filter(output_df, group=="HOI") )
model_null <- lm(std_y7 ~ log2(div), data = filter(output_df, group=="no HOI") )

# the increase of BEF relationships by HOIs effects
(model_HOI$coefficient[2] - model_null$coefficient[2])/min(model_HOI$coefficient[2], model_null$coefficient[2])



# prepare required data for plotting the BEF relationship with and without net HOIs
df_HOI <- filter(output_df, group=="HOI") |> (function(x) split(x, f=x$BEF_id))()
# result is a 71(newdata) by 1000 matrix
BEF_HOI <- sapply(df_HOI, function(x) predict(lm(y7 ~ log2(div), data=x), 
                                              newdata=data.frame(div=seq(1,8,0.1))), simplify = T)
Q_HOI <- as.data.frame(t(apply(BEF_HOI, 1, function(x) quantile(x, probs=c(0.05, 0.95)))))
Q_HOI$div <- seq(1, 8, 0.1)
colnames(Q_HOI)[1:2] <- c("lower", "upper")
Q_HOI$y7 <- predict(model_HOI, newdata=data.frame(div=seq(1,8,0.1)))


df_null <- filter(output_df, group=="no HOI") |> (function(x) split(x, f=x$BEF_id))()
BEF_null <- sapply(df_null, function(x) predict(lm(y7 ~ log2(div), data=x), 
                                                newdata=data.frame(div=seq(1,8,0.1))), simplify = T)
Q_null <- as.data.frame(t(apply(BEF_null, 1, function(x) quantile(x, probs=c(0.05, 0.95)))) )
Q_null$div <- seq(1, 8, 0.1)
colnames(Q_null)[1:2] <- c("lower", "upper")
Q_null$y7 <- predict(model_null, newdata=data.frame(div=seq(1,8,0.1)))


p_sim <- ggplot() +
  geom_ribbon( data=Q_HOI, aes(x=log2(div),y=y7, ymin = lower, ymax= upper), alpha=0.1, fill="#6C8A9E") +
  geom_ribbon( data=Q_null, aes(x=log2(div),y=y7, ymin = lower, ymax= upper ), alpha=0.1, fill="#C4C27D") +
  
  geom_smooth(data=output_df, aes(log2(div), y7, color = factor(group)),
              method = "lm", se = FALSE, linewidth = 1.6) +
  
  geom_line(data=Q_HOI, aes(x=log2(div),y=lower), alpha=0.5, colour="#6C8A9E", linewidth=0.6, linetype="dashed") +
  geom_line(data=Q_HOI, aes(x=log2(div),y=upper), alpha=0.5, colour="#6C8A9E", linewidth=0.6, linetype="dashed") +
  
  geom_line(data=Q_null, aes(x=log2(div),y=lower), alpha=0.5, colour="#C4C27D", linewidth=0.6, linetype="dashed") +
  geom_line(data=Q_null, aes(x=log2(div),y=upper), alpha=0.5, colour="#C4C27D", linewidth=0.6, linetype="dashed") +
  
  scale_color_manual(values = c("#6C8A9E","#C4C27D"),
                     name = "",  # Change legend title
                     labels = c("with net HOI effects", "no net HOI effects")) + 
  labs(x = expression(log[2]("div")), y =  "Productivity (g)",
       title = "BEF relationships", tag = "A" ) +
  theme_classic() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 12),
        plot.tag = element_text(face = "bold", size = 12),
        plot.tag.position = c(0.1, 1),
        axis.title.x = element_text(size = 10, margin = margin(b = 5, unit = "pt")),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        plot.margin = margin(r = -5,t=5,b=3,l = 1, unit = "pt"))


ggsave(filename = "fig_3_new.png",  plot = p_sim, device = "png",  
       width = 75,  height = 100,  units = "mm",  dpi = 600,  
       bg = "white" )




#------------calculate BEF using empirical data (7th year growth) -------------#
# load the 7th year tree growth data from focused plots
year_7_gr <- read.csv("year_7_growth.csv")

# calculated empirical BEF relationships at the 7th year
emp_bef <- lm(gr_std ~ log2(plot_richness), data = year_7_gr)$coefficient[2]
emp_se <- summary(lm(gr_std ~ log2(plot_richness), data = year_7_gr))$coefficients[2, "Std. Error"]

# standardized BEF slope with net HOI effects
HOI_bef <- sapply(df_HOI, function(x) lm(std_y7 ~ log2(div), data = x)$coefficient[2])

# standardized BEF slope without net HOI effects
null_bef <- sapply(df_null, function(x) lm(std_y7 ~ log2(div), data = x)$coefficient[2])

# create data for plotting
df_slope <- data.frame(mean_slope = c(mean(null_bef), mean(HOI_bef), emp_bef),
                       CI_lower = c(quantile(null_bef, probs = 0.05),quantile(HOI_bef, probs = 0.05),(emp_bef-emp_se)),
                       CI_upper = c(quantile(null_bef, probs = 0.95),quantile(HOI_bef, probs = 0.95),(emp_bef+emp_se)),
                       group = c("no HOI", "HOI", "Empirical"))


p_slope <- ggplot(df_slope, aes(x = group, y = mean_slope, fill = group, color = group)) +
  geom_errorbar(
    aes(ymin = CI_lower, ymax = CI_upper), width = 0.15, alpha = 0.7,
    linewidth = 1, show.legend = F) +
  geom_point(size = 4, shape = 21, alpha = 0.9) +
  scale_fill_manual(values = c("no HOI" = "#C4C27D", "HOI" = "#6C8A9E", "Empirical" = "#874F64")) +
  scale_color_manual(values = c("no HOI" = "#C4C27D", "HOI" = "#6C8A9E", "Empirical" = "#874F64")) +
  labs(x = "", y = "Effect size", tag = "B", 
       title = "Standardized BEF relationships") +
  coord_flip() +
  theme_classic() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 12),
        plot.tag = element_text(face = "bold", size = 12),
        plot.tag.position = c(0.1, 1),
        axis.title.x = element_text(size = 10, margin = margin(b = 5, unit = "pt")),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        plot.margin = margin(t=5,b=4, r=3,l = -7, unit = "pt"))



# create standalone legend 
legend_data <- data.frame(
  group = rep(c("with net HOI effects", "no net HOI effects", "Empirical"), each = 10),
  x = rep(1:3, each = 10),
  y = c(rnorm(10, 1,0.01), rnorm(10, 2, 0.02), rnorm(10, 3, 0.03)) )


legend_plot <- ggplot(legend_data, aes(x, y, color = group)) +
  geom_line(linewidth=1.6) +
  scale_color_manual(values = c("#6C8A9E", "#C4C27D", "#874F64"), name = "",
                     labels = c("with net HOI effects", "no net HOI effects", "Empirical")) +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.text = element_text(size = 10, margin = margin(r = 10, unit = "pt")),
    legend.key.height = unit(10, "mm"),
    legend.key.width = unit(10, "mm")
  )

# extract legend grob
plot_legend <- get_legend( legend_plot + 
                        theme(legend.position = "right",# defaut legend position
                              legend.text = element_text(size = 10),
                              legend.margin = margin(t = -1, unit = "mm")) + 
                        guides(fill = guide_legend(nrow = 1)) )


# combine plots and legend
fig_3 <- plot_grid(
  plot_grid(p_sim, p_slope, nrow = 1, align = 'h'),
  plot_legend,
  nrow = 2,
  rel_heights = c(1, 0.1)  # Adjust legend height as needed
)


ggsave(filename = "fig_4_growth.png",  plot = fig_3, device = "png",  
       width = 150,  height = 100,  units = "mm",  dpi = 600,  
       bg = "white" )

