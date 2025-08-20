#------------------------------------------------------------------------------#
#------ quantify interactions per plot for empirical plot arrangement ---------#
#----- structure equations to quantify effects of interaction on biomass-------#
#------------------------------------------------------------------------------#

# IS-div-analysis.R need to be run first to have res_mat and null_mat available

#--------------------- extract the plot map for each plot ---------------------#
# 6x6 of mono-cultures and 2-species mixtures, 12x12 of 4, 8 species mixtures 
# load configuration of each plot
BEF_A_8 <- read.csv("BEF_plot_map.csv", row.names = 1)

#plot_count <- BEF_A_8 %>% count(plotID_2) # should be 74 plots
plot_count <- BEF_A_8 %>% group_by(plotID_2, plot_richness) %>%  summarise(n = n(), .groups = "drop")
plot_count %>% count(n)

# replace species name with numbers to correspond to interaction matrix
sp_level <- c( "Castanea henryi" ,"Castanopsis sclerophylla", "Choerospondias axillaris",
               "Liquidambar formosana" , "Nyssa sinensis", "Quercus serrata",         
               "Sapindus mukorossi", "Sapium sebiferum")

match("Sapium sebiferum", sp_level )

# split the data into each plot
plot_map_list <- split(BEF_A_8, f=BEF_A_8$plotID_2)

# split big plots(12 by 12 to 6 by 6)
# create the info for correctly subsetting the big plots 12X12
mat_idx_12 <- matrix(1:12^2, nrow=12, byrow = T)

# define sub-matrix indices
sub_matrices <- list( mat_idx_12[1:6, 1:6],  mat_idx_12[1:6, 7:12], 
                      mat_idx_12[7:12, 1:6], mat_idx_12[7:12, 7:12] )

# extract positions for each sub-matrix
sub_rows <- lapply(sub_matrices, function(mat) {mat_idx_12[mat]} )

# split plots to ensure same number of trees and store the plot map
plot_list <- list() # length should be 116
index <- 1  # Track position in plot_list
# correct the coordinates of tree for each plot
for (i in 1:length(plot_map_list)) {
  output <- plot_structure(plot_map_list[[i]])  
  
  if (is.data.frame(output)) {
    plot_list[[index]] <- output  # store the single data frame
    index <- index + 1
  } else if (is.list(output)) {
    for (df in output) {
      plot_list[[index]] <- df  # store each data frame separately
      index <- index + 1
    }
  }
}

# check the misplanting
print(unique(sapply(plot_list, function(x) length(unique(x$sp))))) 
misplant_id <- which(! sapply(plot_list, function(x) length(unique(x$sp))) %in% c(1,2,4,8))

# remove weird plot (too many misplanting)
misplant_id <- misplant_id[ !misplant_id %in% c(60, 90, 96)]
# manually checked all the 13 misplanted plot, most plots only one tree was misplanted
# however plot_list[[60]] should have 4 species, but only 3 species present in reality should be excluded
# plot_list[[96]] has 10 species, Koelreuteria bipinnata has 3 trees, Quercus fabri has 1

# correct the misplanting by randomly replacing misplanted species with specie in that plot
plot_list[[90]][plot_list[[90]]$sp=="Quercus acutissima" |
                  plot_list[[90]]$sp=="Sapium sebiferum", "sp"] <- "Castanea henryi"

plot_list[[96]][plot_list[[96]]$sp=="Koelreuteria bipinnata" |
                  plot_list[[96]]$sp=="Quercus fabri",  "sp"] <- "Quercus serrata"

for (i in 1:length(misplant_id)) {
  idx = misplant_id[i]
  plot_list[[idx]] = correct_misplant(plot_list[[idx]])
}

# check again 3 species plot:plot_list[[60]] should be exclude in the end from SEM
print(unique(sapply(plot_list, function(x) length(unique(x$sp))))) 


# number of individual per plot
nind <- 36 

# extract the interaction matrix for each plot
plot_IS_div <- vector("list", length=length(plot_list)) 
plot_IS_null <- vector("list", length=length(plot_list)) 
empi_HOI <- vector("list", length=length(plot_list)) 
empi_null <- vector("list", length=length(plot_list))

# extract the IS group: pos, neg corresponding to interaction matrix for each plot
plot_IS_group <- vector("list", length=length(plot_list))

set.seed(1230) 
for (i in 1:length(plot_list)){
  
  dat <- plot_list[[i]]
  dat$sp_id <- sapply(dat$sp, function(x) match(x, sp_level) )
  mat <-  matrix(dat$sp_id, nrow = sqrt(nrow(dat)), byrow = T)
  coord <- as.matrix(dat[ ,c("cell_x", "cell_y")])
  
  # random starting biomass
  start_bm <- matrix(rlnorm(nrow(dat), meanlog = 3.5, sdlog = 1), nrow = sqrt(nrow(dat)))
  
  # get the start biomass for focal tree for each plot (by column)
  plot_start_bm <- c(start_bm)
  # get the neighbor biomass for each plot (by column)
  plot_nb_bm <- do.call(rbind, apply(coord, MARGIN = 1, FUN = neighbour_info, mat= start_bm, simplify = F))
  
  # get the neighbor id for each plot
  df_nb_id <- do.call(rbind, apply(coord, MARGIN = 1, FUN = neighbour_info, mat= mat, simplify = F))
  
  # get the interaction strength based on plot map WITH interaction modification
  # if misplanted, all IS are 0 for the following three lists
  plot_IS_div[[i]] <- search_IS(focal_id = dat$sp_id, df_nb_id, res_mat[ ,1:4])
  
  # interaction strength based on plot map WITHOUT interaction modification
  plot_IS_null[[i]] <- search_IS(focal_id = dat$sp_id, df_nb_id, null_mat[ ,1:4])
  
  # assign each interaction into pos and neg group according to null_mat or alpha_matrix?!!
  plot_IS_group[[i]] <- interaction_group_2(sp_id = dat$sp_id, nb_id = df_nb_id, ref_IS = null_mat$two  )  
  
  # simulate tree growth for net HOI and without net HOI
  empi_HOI[[i]] <- sim_growth(focal_index = coord,
                              focal_id = dat$sp_id,
                              focal_bm = plot_start_bm,
                              div = length(unique(dat$sp_id)),
                              df_IS = plot_IS_div[[i]],
                              nb_bm = plot_nb_bm)
  empi_HOI[[i]]$plot_id <- i
  
}


# calculate average biomass per plot and create data.frame for plotting
# remove plot 60 as it has only 1 alive tree 
empi_plot_HOI <- do.call(rbind, empi_HOI[-60]) %>% group_by(plot_id) %>% summarise(y7 = mean(y7), div=first(div))



#------------------------------------------------------------------------------#
#------------------------ SEM based on empirical arrangement ------------------#
# load packages for SEM
library(lavaan)
library(lavaanPlot)

# calculate summary statistics
summary_vals <- vector("list", length(plot_IS_div))
stopifnot(length(plot_list) == length(plot_IS_div))

for (i in 1:length(plot_IS_div)){
  summary_vals[[i]] <- plot_IS_summary( plot_IS_div[[i]])  #plot_IS_group[[i]]
}


# combine interaction metric with biomass
dat_sem <- cbind(empi_plot_HOI, do.call(rbind, summary_vals)[-60, ])
colnames(dat_sem)[4:ncol(dat_sem)] <- c("pos", "neg", "pos_sum", "neg_sum" , "pos_fr", "neg_fr")
dat_sem$log2div <- log2(dat_sem$div)
dat_sem$y7_log <-log(dat_sem$y7)


# just check the NAs, mono-culture don't have pos, neg, make sense
# this should not hold (sometimes some plots don't have pos or neg)
na_indices <- sapply(summary_vals, function(x) any(is.na(x)))  
identical(which(na_indices), which(dat_sem$div==1))

# change NA to 0
dat_sem <- dat_sem %>% replace(is.na(.), 0)

# SEM model
model_test <- 'y7_log ~ pos + neg + log2div
               pos ~ log2div
               neg ~ log2div'

fit_test <- sem(model_test, data=dat_sem) 

summary(fit_test, rsq=T, standardized=T, fit.measures=T)

lavaanPlot(model = fit_test, coefs = TRUE, stand=TRUE, sig = 0.05)


