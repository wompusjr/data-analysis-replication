#necessary libraries
library(tidyverse)
library(skimr)
library(purrr)
library(igraph)
library(sna)
library(corrplot)
library(doParallel)
library(foreach)
library(tibble)
library(reshape2)
library(viridis)
library(ecp)
library(fossil)
#Load data
d <- read_csv("data/data_all copy.csv")
ceramic <- read_csv("data/Ceramic_type_master copy.csv")
coords <- read_csv("data/attr_all copy.csv")
#BR similarity matrix function #adapted from someone elses function
BRsim <- function(x, correction = TRUE, rescale = TRUE, ncores = parallel::detectCores() - 1) {
  cl <- parallel::makeCluster(ncores)
  doParallel::registerDoParallel(cl)
  x <- as.matrix(x)
  rd <- nrow(x)
  rownames_x <- rownames(x)
  norm_x <- sweep(x, 1, rowSums(x), "/")
  sim_matrix <- matrix(NA, nrow = rd, ncol = rd)
  rownames(sim_matrix) <- rownames_x
  colnames(sim_matrix) <- rownames_x
  row_pairs <- expand.grid(s1 = 1:rd, s2 = 1:rd)
  results <- foreach(i = 1:nrow(row_pairs), .combine = rbind) %dopar% {
    s1 <- row_pairs[i, 1]
    s2 <- row_pairs[i, 2]
    diff <- abs(norm_x[s1, ] - norm_x[s2, ])
    sim_score <- 1 - (sum(diff) / 2)
    if (correction) {
      zero_a <- sum(x[s1, ] == 0)
      zero_b <- sum(x[s2, ] == 0)
      joint_zero <- sum(x[s1, ] == 0 & x[s2, ] == 0)
      divisor <- if (zero_a == zero_b) 1 else max(zero_a, zero_b) - joint_zero + 0.5
      sim_score <- sim_score / divisor
    }
    sim_score <- if (rescale) sim_score else sim_score * 200
    c(s1, s2, round(sim_score, 3))
  }
  for (i in 1:nrow(results)) {
    s1 <- results[i, 1]
    s2 <- results[i, 2]
    sim <- results[i, 3]
    sim_matrix[s1, s2] <- sim
  }
  sim_matrix[lower.tri(sim_matrix)] <- t(sim_matrix)[lower.tri(sim_matrix)]
  parallel::stopCluster(cl)
  return(sim_matrix)
}
#Clean and join to main dataframe
dj <- d |>
  left_join(ceramic |> mutate(Type = SWSN_Type) |> select(Type, Ware = SWSN_Ware), by = "Type") |>
  left_join(coords |> mutate(Site = SWSN_Site) |> select(Site, EASTING, NORTHING), by = "Site")
#Filter by connection to cyberSW, at least 30 painted sherds, and between 1000 and 1450
d2 <- dj |>
  filter(cyberSW == 1, painted == 1, Begin >= 1000, End <= 1450)
temp <- d2 |>
  group_by(Site) |>
  summarise(True_Count = sum(Count)) |>
  mutate(remove = if_else(True_Count >= 30, "1", "0"))
dj2 <- d2 |>
  left_join(temp |> select(Site, remove), by = "Site") |>
  filter(remove == 1) |>
  select(-remove)
#Binning temporal intervals
df <- dj2 |>
  rowwise() |>
  mutate(bin_starts = list(seq(floor(Begin / 50) * 50, floor(End / 50) * 50, by = 50))) |>
  unnest(bin_starts) |>
  mutate(bin_end = bin_starts + 49,
         overlap_start = pmax(Begin, bin_starts),
         overlap_end = pmin(End, bin_end),
         overlap_years = overlap_end - overlap_start + 1) |>
  filter(overlap_years > 10) |>
  mutate(interval = paste0(bin_starts, "-", bin_end)) |>
  select(-overlap_start, -overlap_end, -overlap_years) |>
  ungroup()
##descriptive analysis
length(unique(df$Site)) #1766
length(unique(d$Type)) #over 1156 gets culled down to
length(unique(df$Type)) #295
###these align closely with what was described in the article; however, the sites are off by 24. 
###im not entirely certain, but it could relate to the room-size variable that wasn't on the dataset
#Building data for 1400–1449 interval
d1400 <- df |>
  filter(interval == "1400-1449") |>
  group_by(Site, Ware) |>
  summarise(wareCount = sum(Count), .groups = "drop") |>
  pivot_wider(names_from = Ware, values_from = wareCount, values_fill = 0) |>
  column_to_rownames("Site")
#Run similarity test on the subsetted data
sim1400 <- BRsim(d1400)
#create the distance matrix
coord1400 <- df |> filter(interval == "1400-1449") |> select(Site, EASTING, NORTHING) |> distinct()
coordmatrix <- as.matrix(dist(coord1400[, c("EASTING", "NORTHING")]))
rownames(coordmatrix) <- coord1400$Site
colnames(coordmatrix) <- coord1400$Site
# Filter matrices for matching site names - make sure everything is honky-dory
common_sites <- intersect(rownames(sim1400), rownames(coordmatrix))
sim1400 <- sim1400[common_sites, common_sites]
coordmatrix <- coordmatrix[common_sites, common_sites]

#Louvain community detection
#Setting up result storage
res <- matrix(NA, nrow(sim1400), 100)
#Loop through 100 different distance thresholds to generate community partitions
for (i in 1:100) {
  # Convert coordmatrix to binary matrix based on percentile threshold
  distmatrix <- event2dichot(coordmatrix, method = 'quantile', thresh = 1 - (i / 100))
  #creating the community matrix
  community <- distmatrix * sim1400
  #replace N/A's w/ 0
  community[is.na(community)] <- 0
  community[is.nan(community)] <- 0
  # Create weighted graph
  G <- graph_from_adjacency_matrix(community, mode = "undirected", weighted = TRUE, diag = TRUE)
  # Community detection with Louvain algorithm
  set.seed(63246)
  res[, i] <- cluster_louvain(G, weights = E(G)$weight)$membership
}
#Calculate Adjusted Rand Index
#setting up cores
registerDoParallel(cores = parallel::detectCores() - 1)
#creating a rand index with :fossil:
rand_index <- foreach(i = 1:100, .combine = rbind) %dopar% {
  sapply(1:100, function(n) {
    fossil::adj.rand.index(res[, i], res[, n])
  })
}
rand_index[is.na(rand_index)] <- 0 #sorting out NAs
#getting protoypical scales using :ecp:
scale1 <- e.agglo(rand_I, alpha = 1, penalty = function(cp, Xts)0) # calculate across all distances
scale1a <- scale1$estimates - 1 #isolating the estimates
scale2 <- e.agglo(rand_I[1:25, 1:25], alpha = 1.5, penalty = function(cp, Xts)0) # calculate within the lowest quartile
scale2a <- scale2$estimates - 1 #isolating the estimates
scale2a <- scale2a[-length(scale2a)]
# combine and remove duplicates
scale <- sort(unique(c(scale1a, scale2a)))
scale[1] <- 1
scale <- unique(scale)
#calculate prototypical distances for each partition and create data frame
proto.typical <- NULL
for (i in 1:(length(scale) - 1)) {
  lv <-  scale[i]:scale[i + 1]
  proto.typical[i] <-  lv[which.max(rowSums(rand_I[lv, lv]))]
}
pscale <- as.data.frame(proto.typical)
colnames(pscale) <- c('X')
view(pscale)
#Plot ARI heatmap
#creating a dataframe for ggplot
ari_df <- as.data.frame(rand_I) |>
  rownames_to_column(var = "x") |>
  pivot_longer(-x, names_to = "y", values_to = "value") |>
  mutate(x = as.integer(x), y = as.integer(y))
#making the plot
ari <- ggplot(ari_df, aes(x = x, y = y, fill = value)) +
  geom_tile() +
  scale_fill_viridis(name = "ARI", option = "D", limits = c(0, 1)) +
  geom_vline(xintercept = scale + 0.5, col = "white") +
  geom_hline(yintercept = scale + 0.5, col = "white") +
  geom_point(data = pscale, mapping = aes(x = X, y = X), size = 4, col = 'red', inherit.aes = FALSE) +
  scale_x_continuous(expand = c(0, 0), breaks = seq(0, 100, 25)) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 100, 25)) +
  coord_equal() +
  labs(x = "Percent of Distance Considered", y = "Percent of Distance Considered", title = "1400–1449 CE") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_line(color = "white", size = 1),
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 16, face = "bold")
  )
#showing the plot
par(mfrow=c(2,2))
print("images/fig4")
print(ari) #this displays that during the final 50year span, there were 6 key prototypical spans