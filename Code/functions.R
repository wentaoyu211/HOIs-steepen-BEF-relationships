#----------------------------required functions -------------------------------#

# 1.----------------------------------------------------------------------------
# function to extract neighbour identity and neigbors' biomass depending on input mat
# coord(i,j) indicates the location of the focal tree in the matrix
# mat is the spatical scheme matrix n by n
# applied periodic boundary condition to avoid edge effect, each tree has 8 neighbours
neighbour_info = function(coord, mat){ 
  
  i = coord[1] 
  j = coord[2] 
  
  # x, y are the locations of the neighboring tree in the grid
  x = c(i-1, i-1, i-1,  i,   i,  i+1, i+1, i+1) 
  y = c(j-1,  j,  j+1, j-1, j+1, j-1,  j,  j+1)
  
  # Apply periodic boundary conditions
  x = ((x - 1) %% nrow(mat)) + 1  # Wrap around row indices
  y = ((y - 1) %% ncol(mat)) + 1  # Wrap around column indices
  
  nb_id <- mat[cbind(x, y)]
  
  # check if it's 8 neighbors
  stopifnot(length(nb_id) == 8)
  
  return(nb_id)
}


# 2.-------------fill in interaction between focal and neighbors--------------#
search_IS <- function(focal_id, df_nb_id, IS_div_mat){
    div = length(unique(focal_id)) # div

    interaction_matrix <- matrix(IS_div_mat[, log2(div)+1], nrow = 8, byrow = T)
    df_nb_IS <- matrix(interaction_matrix[cbind(rep(focal_id, each = 8), as.vector(t(df_nb_id)))], 
                            nrow = nrow(df_nb_id), ncol =ncol(df_nb_id), byrow = TRUE)
    
    return(df_nb_IS)
  }


# 3. ---------------------------------------------------------------------------
# tree growth for each plot
sim_growth <- function(focal_index, focal_id, focal_bm, div, df_IS, nb_bm){
  
  dat <- data.frame(div=rep(div, nind) )

for (j in 1:year){
  # calculating tree growth across years
  interactions <- rowSums(df_IS*nb_bm^b)
  
  bm_next <-  beta[focal_id]*focal_bm^theta + interactions
  
  # if there is negative value it denotes dead tree
  # but need to make sure they stay dead
  bm_next[bm_next <= 0] <- 0.0
  #bm_next[focal_bm > bm_next] <- 0 # when there is negative growth, tree is dead
  bm_next[focal_bm == 0] <- 0.0
  
  dat[[paste0("y", j)]] <- bm_next
  
  focal_bm <- bm_next
  
  nb_bm <- do.call(rbind, apply(focal_index, MARGIN = 1, FUN = neighbour_info, 
                                   mat= matrix(bm_next, nrow = sqrt(nrow(dat)), byrow = T), simplify = F))
}
  return(dat)
}


# 4. ---------------------------------------------------------------------------
# function to simulate tree growth for one virtual BEF experiment
# with net HOI versus without HOI effects
# for mixture depends on both starting biomass and spatial planting
# record each year 
sim_bef <- function(sp_comp, nind, IS_div_mat, IS_null_mat){
  
  # created index for focal tree way always by column
  focal_index <- as.matrix(expand.grid(1:sqrt(nind), 1:sqrt(nind)))
  
  # create empty list to store output of each plot
  dat_HOI <- vector("list", length=length(sp_comp))
  dat_null <- vector("list", length=length(sp_comp))
  
  # loop through diversity gradient/plot 
  for (i in 1:length(sp_comp)){
    div <- length(sp_comp[[i]])
    # generate plot configuration and draw initial biomass
    plot_map <- sample(rep(sp_comp[[i]], nind), size = nind) |> (\(x)matrix(x, ncol=sqrt(nind)))()
    start_bm <- matrix(rlnorm(nind, meanlog = 3.5, sdlog = 1), nrow = sqrt(nind))
    #start_bm <- matrix(rep(30, nind), nrow = sqrt(nind))
    
    # focal tree id (by column) and biomass
    focal_id <- c(plot_map)
    focal_bm <- c(start_bm)
    stopifnot(length(focal_id)==length(focal_bm))
    
    # extract id of tree neighbors and their biomass
    df_nb_id <- do.call(rbind, apply(focal_index, MARGIN = 1, FUN = neighbour_info, mat= plot_map, simplify = F))
    df_nb_bm <- do.call(rbind, apply(focal_index, MARGIN = 1, FUN = neighbour_info, mat= start_bm, simplify = F))
    
    # extract interaction between focal tree and neighbor over the years
    # for scenarios with HOI versus without HOI
    df_IS_HOI <- search_IS(focal_id, df_nb_id, IS_div_mat)
    df_IS_null <- search_IS(focal_id, df_nb_id, IS_null_mat)
    
    dat_HOI[[i]] <- sim_growth(focal_index, focal_id, focal_bm,
                               div, df_IS_HOI, df_nb_bm)
    dat_HOI[[i]]$plot_id <- i
    
    dat_null[[i]] <- sim_growth(focal_index, focal_id, focal_bm,
                                div, df_IS_null, df_nb_bm)
    dat_null[[i]]$plot_id <- i
  }
  # standardize y7 biomass
  HOI <- do.call(rbind, dat_HOI)
  HOI$std_y7 <- scale(HOI$y7)
  no_HOI <- do.call(rbind, dat_null)
  no_HOI$std_y7 <- scale(no_HOI$y7)
  output <- rbind(HOI, no_HOI)
  output$group <- rep(c("HOI", "no HOI"), each= nrow(output)/2)
  return(output)
}


# 5.---------------------------------------------------------------------------#
# need a function to straighten the plot configuration due to the misplanting
# plot_map: a data frame for each plot contain tree coordinates and other info
plot_structure <- function(plot_map){
  stopifnot(c("cell_x", "cell_y") %in% colnames(plot_map)) # check the data
  # some cell_x, cell_y are not in order due to misplanting, reorder them
  plot_df <- plot_map[order(plot_map[ ,"cell_x"], plot_map[ ,"cell_y"]), ]
  plot_df_list <- vector("list", length = 144/36)
  # and rescale the cell_x, cell_y 
  if (nrow(plot_map) == 36) { plot_df$cell_x = plot_df$cell_x -6
  plot_df$cell_y = plot_df$cell_y -6  
  return (plot_df)}
  if (nrow(plot_map) == 144){ plot_df$cell_x = plot_df$cell_x -3
  plot_df$cell_y = plot_df$cell_y -3  
  plot_df_list <- lapply(sub_rows, function(x) { df <- plot_df[x, ] 
  df[, c("cell_x", "cell_y")] <- expand.grid(1:6, 1:6) |> (\(x) x[, c(2,1)])()       
  return (df)} )
  }
}

# 6.---------------------------------------------------------------------------#
# function to replace misplanted trees with specie in that plot
correct_misplant <- function(plot_map){
  sp_count = plot_map %>% count(sp)
  wrong_sp <- sp_count[sp_count$n == 1, "sp"]
  replace_sp <- sp_count[sp_count$n == sort(sp_count$n)[2], "sp"]
  if (length(replace_sp)==2 | length(replace_sp)==3){
    plot_map[plot_map$sp==wrong_sp, "sp"] = replace_sp[1]
  } else {plot_map[plot_map$sp==wrong_sp, "sp"] = replace_sp}
  return(plot_map)
}
# 17 18 25(double) 30(double) 72 73(double) 74 75 90(different double) 92 95


# 7. ---------------------------------------------------------------------------
# function to create a matrix to identify  pos, neg IS groups
# use null_mat, intra-specific interaction from mono culture
# inter-specific interaction from div 2 to define IS groups
# sp_id: identity of focal tree
# nb_id, data.frame of identities of neighbouring trees
# ref_IS: reference interaction strength (IS) as mentioned

interaction_group_2 <- function(sp_id, nb_id, ref_IS){
  
  ref_mat <- matrix(ref_IS, nrow = 8, byrow = T)
  
  pos_idx <- which(ref_mat > 0 , arr.ind = TRUE)
  neg_idx <- which(ref_mat < 0 , arr.ind = TRUE)
  
  inter_group <- matrix(0, nrow = nrow(nb_id), ncol = ncol(nb_id))
  for (i in 1:nrow(nb_id)) {
    focal_species <- sp_id[i]  # Get focal tree species
      idx_function <- function(nb, focal_id){
        if (any(apply(pos_idx, 1, function(row) all(row == c(focal_id, nb))))) {return ("pos")}
        else if (any(apply(neg_idx, 1, function(row) all(row == c(focal_id, nb))))) {return ("neg")}
        else {return (as.numeric(0)) }  
      }
      # identify the group based species identity
      inter_group[i, ] <- sapply(nb_id[i, ], FUN = idx_function, focal_id = focal_species )
    
  }
  return(inter_group)
}


# 8. ---------------------------------------------------------------------------
# function to calculate mean of positive and negative IS for each plot (2 IS groups)
plot_IS2_summary <- function(plot_IS, IS_group){
  
  pos = mean(plot_IS[which(IS_group == "pos", arr.ind=T)])
  pos_sum = sum(plot_IS[which(IS_group == "pos", arr.ind=T)])
  pos_fr = length(which(IS_group == "pos"))/length(which(IS_group != "0" ))
  
  neg = mean(plot_IS[which(IS_group == "neg", arr.ind=T)])
  neg_sum = sum(plot_IS[which(IS_group == "neg", arr.ind=T)])
  neg_fr = length(which(IS_group == "neg"))/length(which(IS_group != "0" ))
  
  
  return(c(pos, neg,
           pos_sum, neg_sum, 
           pos_fr, neg_fr))
}

#------calculate mean of positive and negative IS for each plot just based 
#------ on positive and negative interaction values
plot_IS_summary <- function(plot_IS){
  
  stopifnot(any(plot_IS!=0)) # check if there is 0
  
  pos = ifelse(length(which(plot_IS>0))>0, mean(plot_IS[plot_IS > 0]), 0) # prevent NA
  pos_sum = sum(plot_IS[plot_IS > 0])
  pos_fr = length(which(plot_IS > 0))/(nind*8)
  
  neg = ifelse(length(which(plot_IS<0))>0, mean(plot_IS[plot_IS < 0]), 0)
  neg_sum = sum(plot_IS[plot_IS < 0])
  neg_fr = length(which(plot_IS < 0))/(nind*8)
  
  return(c(pos, neg,
           pos_sum, neg_sum, 
           pos_fr, neg_fr))
}




