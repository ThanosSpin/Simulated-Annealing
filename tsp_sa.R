# Run after tsp.R script

# CONFIGURATION -------------------------------------------------------------------------------------------------------
# the same points with heuristic solution

# DISTANCES -----------------------------------------------------------------------------------------------------------

# Get the distances using the Google Maps API. However, for production purposes this would be done with a local OSRM
# instance.

# Re order the addresses data frame to initial state
addresses_sa <- addresses[order(addresses$label), ]

# Calculate the distances through Google's API
distances_sa <- gmapsdistance(origin = addresses_sa$address,
                              destination = addresses_sa$address,
                                           key = (""),
                                           combinations = "all",
                                           mode = "driving")$Distance

# Change the column names and row names of distances data frame
names(distances_sa)[which(names(distances_sa) %in% c("or"))] <- c("x/x")
colnames(distances_sa)[-which(names(distances_sa) %in% c("x/x"))] <- addresses_sa$label
rownames(distances_sa) <- addresses_sa$label

# Remove the first column (longlat)
distances_sa <- distances_sa[, -1]

# Scale to km.
distances_mat <- (as.matrix(distances_sa))/ 1000

# Change the order of the points
order_mat <- (x = 1:ncol(distances_mat))

# compute the total distance given the ordering of the cities in the path
distance <- function(x) {
  total_distance <- 0
  
  for (i in 1:(length(x) - 1)) {
    total_distance <- total_distance + distances_mat[x[i], x[i+1]]
  }
  
  total_distance <- total_distance + distances_mat[x[length(order_mat)], x[1]]
  return (total_distance)
}

distance(order_mat)

E <- distance(order_mat)

# swap two randomly picked cities in the ordering
swap <- function(y) {
  indices <- sample(1:length(y), 2)
  temp = y[indices[1]]
  y[indices[1]] = y[indices[2]]
  y[indices[2]] = temp
  return(y)
}

# returns the probability with which a solution should be accepted
accptFun <- function(delE, t) {
  return(exp(-delE/t))
} 

# Initilization
t = 1e5
iter <- 1
loop <- 1
Error <- c(E)

# Start the SA optimization process
start_time_sa <- Sys.time()
startTime <- proc.time()
while(proc.time() - startTime < 12) {
  newOrder <- swap(y = order_mat)
  newE <- distance(x = newOrder)
  
  rand <- runif(1, min = 0, max = 1)
  acpt <- accptFun(newE - E, t)
  
  if( (newE - E) < 0 | rand <= acpt) {
    delE <- newE - E
    
    cat("Iteration: ", iter, ", Time: ", proc.time() - startTime, "Temperature: ", t, "\n")
    cat("Order: ", order_mat, ", New Order: ", newOrder, "\n")
    cat("Total Distance: ", E, "Del E: ", delE, ", New Total Distance: ", newE, "\n")
    cat("\n\n")
    
    E <- newE
    order_mat <- newOrder
    t <- t * 0.95
    iter <- iter + 1
    Error[iter] = E
  }
  loop <- loop + 1
}
end_time_sa <- Sys.time()
end_time_sa - start_time_sa

# Plot distance in each iteration
plot(1:iter, Error, xlab = "Accepted Routes", ylab = "Total Distance (km)", ylim = c(0,500), type = "l")


# Sort addresses.
#
addresses_sa <- addresses_sa[order_mat, ]

# BUILD ROUTE ---------------------------------------------------------------------------------------------------------

latlong_sa <- c()

for (i in 1:nrow(addresses_sa)) {
  latlong_sa <- c(latlong_sa, list(paste(addresses_sa$lat[i], addresses_sa$lon[i], sep = ",")))
}

route_sa <- lapply(seq(nrow(addresses_sa) - 1), function(x) {
  Sys.sleep(1)
  googleway::google_directions(latlong_sa[[x]], latlong_sa[[x + 1]], mode = "driving") 
})

# ---------------------------------------------------------------------------------------------------------------------

latlong_from_to_sa <- c()

for (i in 1:length(route_sa)) {
  latlong_from_to_sa <-c(latlong_from_to_sa, route_sa[[i]]$routes[1])
}

# Keep the pairs lang-long 
latlong_from_to_sa <- gsub("[^0-9\\.,]", "", (as.character(latlong_from_to_sa)))

# Split the pairs of lang-long 
split_latlong_sa <- c()

for (i in 1:length(route_sa)) {
  split_latlong_sa <-c(split_latlong_sa, strsplit(latlong_from_to_sa[i], ","))
}

# Lang - long df (SA)
latlong_sa_df <- c()

for (i in 1:length(split_latlong_sa)) {
  latlong_sa_df <-  rbind.data.frame(latlong_sa_df, matrix(strsplit(split_latlong_sa[[i]], " "), 2, 2, byrow = T))
}

# Change the variable names
names(latlong_sa_df) <- c("Lat", "Long")

# Polylines (SA) 
pl_sa <- c()
for (i in 1:length(split_latlong_sa)) {
  pl_sa <-  c(pl_sa, direction_polyline(route_sa[[i]]))
}

# data frame pl_sa
df_sa_pl <- data.frame(polyline = pl_sa)

# Map with best route of TSP (SA optimization)
map_sa <-  google_map(location = c(latlong_sa_df$Lat[1], latlong_sa_df$Long[1]), zoom = 12) %>%
  add_polylines(data = df_sa_pl, polyline = "polyline", stroke_weight = 9)
map_sa
