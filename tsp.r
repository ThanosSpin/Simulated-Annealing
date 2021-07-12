
# Clear R Environment
rm(list = ls())

# CONFIGURATION -------------------------------------------------------------------------------------------------------

ADDRESSES = c(
  "Vizantinon Aftokratoron 13, Nea Ionia, Greece",
  "Kolokotroni 3, Nea Ionia",
  "Kolokotroni 12, Nea Ionia",
  "Kolokotroni 58, Nea Ionia",
  "Marinou Antipa 88, Iraklio",
  "Leoforos Kimis 52, Nea Ionia",
  "Leoforos Kimis 112, Nea Ionia",
  "Leof. Galatsiou 54, Athina 111 41",
  "Leof. Galatsiou 132, Athina 111 41",
  "Leof. Galatsiou 33, Athina 111 41",
  "Leof. Irakliou 54, Athina 111 41",
  "Leof. Irakliou 125, Athina 111 41",
  "Andrea Papandreou 35, Marousi 151 22",
  "Leof. Kifisias 37A, Marousi 151 23",
  "Patision 135, Athina 112 51",
  "3is Septemvriou 65, Athina 104 33",
  "Tatoiou 14, Metamorfosi 144 51",
  "Agia Olga Periferiako Geniko Nosokomio Neas Ionias Konstantopoulio, Agias Olgas, Nea Ionia",
  "Leoforos Mesogeion 12, Athens",
  "Leof. Mesogeion 264, Cholargos 155 62",
  "Katechaki 12, Athina 115 25",
  "Leof. Alexandras 104",
  "Leof. Andrea Siggrou 115, Athina 117 45",
  "Leof. Andrea Siggrou 95, Athina 117 45",
  "Leof. Andrea Siggrou 100, Athina 117 41",
  "Leof. Andrea Siggrou 90, Athina 117 41",
  "Pl. Sintagmatos, Athina 105 63",
  "Leoforos Vasilisis Amalias 10, Athina 105 57",
  "Kallimarmaro, Leoforos Vasileos Konstantinou, Athens",
  "Vasilissis Sofias 46, Athina 115 28",
  "National Technical Univeristy of Athens, Zografou, Greece",
  "National and Kapodistrian University of Athens, Athens, Greece",
  "Leof. Kifisias 99, Marousi 151 24",
  "Leof. Kifisias 29, Marousi 151 23",
  "Leof. Kifisias 129, Marousi 151 24",
  "Leof. Kifisias 220, Kifisia 145 62",
   "Leof. Kifisias 39, Marousi 151 23",
  "Argiroupoleos 62, Argiroupoli 164 51",
  "Leof. Kiprou 94, Vironas 162 32",
  "Leof. Chrisostomou Smirnis 90, Vironas 162 32",
  "Leof. Kiprou 102-104, Vironas 162 32",
  "Voriou Ipirou 2, Vironas 162 32",
  "Agias Sofias 43, Vironas 162 32",
   "Pefkon 1, Iraklio 141 22",
  "Pefkon 54, Iraklio, Greece",
  "Pefkon 79, Iraklio, Greece",
  "Patission 129, Athens, Greece",
  "Patission 39, Athens, Greece",
  "Patission 84, Athens, Greece",
  "Patission 124, Athens, Greece"
)

# LIBRARIES -----------------------------------------------------------------------------------------------------------

library(dplyr)
library(purrr)
library(ggmap)
library(gmapsdistance)
library(TSP)
library(googleway)

# GEOCODE -------------------------------------------------------------------------------------------------------------

#Setting the API key:
set.api.key("")
set_key("")

# Empty vector
reply <- c()

# Functio for getting addresses from Google's API
geocode_df <- function(address) {
for (i in 1:length(address)) {
    reply <- c(reply, google_geocode(address[i], key = ""))
}
reply
}
  
# Empty vector
addresses <- c()  

# Get the addresses 
addresses <- geocode_df(ADDRESSES)

address_format <- c()

for (i in seq(1, length(addresses), by = 2)) {
  address_format <-  rbind(address_format, tibble(
    address = addresses[[i]]$formatted_address,
    lon = addresses[[i]]$geometry$location$lng,
    lat = addresses[[i]]$geometry$location$lat
  )
  )
}


addresses <- address_format %>%
  mutate(latlon = sprintf("%f+%f", lat, lon)) %>%
  mutate(label = LETTERS[1:nrow(.)]) %>%
  select(label, everything())
for (i in 1:24) {addresses$label[26 + i] <-  paste("A", addresses$label[i], sep = "")}

# DISTANCES -----------------------------------------------------------------------------------------------------------

# Get the distances using the Google Maps API. However, for production purposes this would be done with a local OSRM
# instance.


addresses$address <- gsub(" ", "+", addresses$address)
 
distances <- gmapsdistance(origin = addresses$address,
                           destination = addresses$address,
                           key = ("pass"),
                           combinations = "all",
                           mode = "driving")$Distance[, -1]

# Scale to km.
distances <- as.matrix(distances) / 1000

colnames(distances) <- addresses$label
# colnames(distances)[27:30] <- c("AA", "AB", "AC", "AD")
#[which(addresses$label %in% c("A", "B", "C", "D", "E"))]
rownames(distances) <- addresses$label
# rownames(distances)[27:30] <- c("AA", "AB", "AC", "AD")
#[which(addresses$label %in% c("A", "B", "C", "D" , "E"))]

# Convert to distance matrix.
distances <- as.dist(distances)


# TRAVELLING SALESMAN -------------------------------------------------------------------------------------------------

tsp <- TSP(distances)

methods <- c(
  "nearest_insertion",
  "farthest_insertion",
  "cheapest_insertion",
  "arbitrary_insertion",
  "nn",
  "repetitive_nn",
  "two_opt"
)

tours <- methods %>% map(function(method) {
  solve_TSP(tsp, method)
})

start_time <- Sys.time()
tour_km <- c()
tours <- c()
x <- solve_TSP(tsp)
tours <- cbind(matrix(as.integer(x), 1, length(as.integer(x))),  tour_length(x))
for (i in 1:3500) {
  tour <- solve_TSP(tsp)
  x1 <- cbind(t(matrix(as.integer(tour))), tour_length(tour))
  tours <- rbind(tours, c(x1))
  tour_km[i] <- tour_length(tour)
}
end_time <- Sys.time()
end_time - start_time
#
# Order of locations in tour.
# Pick one tour randomly with minimum distance
tour_order <- as.integer(tours[sample(which(tour_km == min(tour_km)), 1), 1:length(addresses$label)])

# Find similar paths in tours data set
which(apply(tours[, 1:length(addresses$label)], 1, function(x) return(all(x == c(tour_order)))))
tours[which(apply(tours[, 1:length(addresses$label)], 1, function(x) return(all(x == c(tour_order)))))[1], ]
tours[which(apply(tours[, 1:length(addresses$label)], 1, function(x) return(all(x == c(tour_order)))))[2], ]
#
# Sort addresses.
#
addresses <- addresses[tour_order,]

# BUILD ROUTE ---------------------------------------------------------------------------------------------------------

latlong <- c()

for (i in 1:nrow(addresses)) {
  latlong <- c(latlong, list(paste(addresses$lat[i], addresses$lon[i], sep = ",")))
  }

route <- lapply(seq(nrow(addresses) - 1), function(x) {
  Sys.sleep(1)
  googleway::google_directions(latlong[[x]], latlong[[x + 1]], mode = "driving") 
})

# ---------------------------------------------------------------------------------------------------------------------

latlong_from_to <- c()

for (i in 1:length(route)) {
  latlong_from_to <-c(latlong_from_to, route[[i]]$routes[1])
  }

# Keep the pairs lang-long 
latlong_from_to <- gsub("[^0-9\\.,]", "", (as.character(latlong_from_to)))

# Split the pairs of lang-long 
split_latlong <- c()

for (i in 1:length(route)) {
  split_latlong <-c(split_latlong, strsplit(latlong_from_to[i], ","))
  }

# Lang - long df
latlong_df <- c()

for (i in 1:length(split_latlong)) {
  latlong_df <-  rbind.data.frame(latlong_df, matrix(strsplit(split_latlong[[i]], " "), 2, 2, byrow = T))
  }

# Change the variable names
names(latlong_df) <- c("Lat", "Long")

# Polylines 
pl <- c()
for (i in 1:length(split_latlong)) {
  pl <-  c(pl, direction_polyline(route[[i]]))
  }

# data frame pl
df_pl <- data.frame(polyline = pl)

# Map with best route of TSP
map <-  google_map(location = c(latlong_df$Lat[1], latlong_df$Long[1]), zoom = 12) %>%
    add_polylines(data = df_pl, polyline = "polyline", stroke_weight = 9)
map
