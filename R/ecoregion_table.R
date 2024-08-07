ecoregion.code.l1 <- data.frame(
  code = c(1:25),
  name = c(
    "Arctic Cordillera",
    "Tundra",
    "Taiga",
    "Hudson Plains",
    "Northern Forests",
    "Northwestern Forested Mountains",
    "Marine West Coast Forests",
    "Eastern Temperate Forests",
    "Great Plains",
    "Noth american Deserts",
    "Mediterranean California",
    "Southern Semi-Arid Highlands",
    "Temperate Sierras",
    "Tropical Dry Forests",
    "Tropical Humid Forests",
    "Caribbean Islands",
    "Northern Andes",
    "Central Andes",
    "Southern Andes",
    "Amazonian-Orinocan Lowland",
    "Eastern Highlands",
    "Gran Chaco",
    "Pampas",
    "Monte-Patagonian",
    "Hawaii"
  )
)

ecoregion.ls <- list(
  ecoregion.1 = c("1", "2"),
  ecoregion.2 = c("3", "4"),
  ecoregion.3 = c("5"),
  ecoregion.4 = c("6", "7"),
  ecoregion.5 = c("8"),
  ecoregion.6 = c("9"),
  ecoregion.7 = c("11"),
  ecoregion.8 = c("10", "12", "13"),
  ecoregion.9 = c("14", "15", "16"),
  ecoregion.10 = c("17", "18", "19"),
  ecoregion.11 = c("20", "21"),
  ecoregion.12 = c("22", "23", "24"),
  ecoregion.13 = c("25")
)

ecoregion.text.ls <- c(
  "Tundra",
  "Taiga-Hudson plains",
  "Northern forest",
  "Northwestern forests",
  "Eastern forests",
  "Great plains",
  "Mediterranean California",
  "Desert-Semiarid highliands",
  "Tropical forests",
  "Andes",
  "Amazonian-Orinocan",
  "Pampas-Patagonian",
  "Hawaii"
)

ecoregion.color.table <-
  data.frame(
    ecoregion = c(as.character(c(1:25))),
    plot.group = c(1, 1,
                   2, 2,
                   3,
                   4, 4,
                   5,
                   6,
                   8, 7, 8, 8,
                   9, 9, 9,
                   10, 10, 10,
                   11, 11,
                   12, 12, 12,
                   13),
    n.site = NA,
    n.site2 = NA,
    # plot.order = c(14, 4, 1, 5, 2, 3,
    #                8, 9, 6, 7,
    #                10, 13, 11, 17, 12, 15, 16),
    plot.text = c(
      "Arctic cordillera",
      "Tundra",
      "Taiga",
      "Hudson plains",
      "Northern forests",
      "NW forests",
      "W coast forests",
      "E temperate forests",
      "Great plains",
      "N American deserts",
      "Mediterranean California",
      "S Semi-arid highlands",
      "Temperate sierras",
      "Tropical dry Forests",
      "Tropical humid Forests",
      "Caribbean islands",
      "N Andes",
      "C Andes",
      "S Andes",
      "Amazonian-Orinocan",
      "Eastern Highlands",
      "Gran Chaco",
      "Pampas",
      "Monte-Patagonian",
      "Hawaii"
    ),
    col = c(
      rgb(0, 0, 255, 255, maxColorValue = 255),
      rgb(0, 50, 255, 255, maxColorValue = 255),
      rgb(0, 120, 255, 255, maxColorValue = 255),
      rgb(0, 170, 255, 255, maxColorValue = 255),
      rgb(200, 255, 80, 255, maxColorValue = 255),
      rgb(100, 255, 50, 255, maxColorValue = 255),
      rgb(50, 255, 25, 255, maxColorValue = 255),
      rgb(0, 255, 0, 255, maxColorValue = 255),
      rgb(255, 220, 100, 255, maxColorValue = 255),
      rgb(255, 0, 0, 255, maxColorValue = 255),
      rgb(255, 150, 150, 255, maxColorValue = 255),
      rgb(245, 165, 0, 255, maxColorValue = 255),
      rgb(255, 255, 0, 255, maxColorValue = 255),
      rgb(150, 255, 150, 255, maxColorValue = 255),
      rgb(100, 200, 100, 255, maxColorValue = 255),
      rgb(50, 200, 50, 255, maxColorValue = 255),
      rgb(170, 175, 255, 255, maxColorValue = 255),
      rgb(90, 120, 220, 255, maxColorValue = 255),
      rgb(75, 80, 180, 255, maxColorValue = 255),
      rgb(255, 0, 255, 255, maxColorValue = 255),
      rgb(200, 0, 200, 255, maxColorValue = 255),
      rgb(0, 255, 255, 255, maxColorValue = 255),
      rgb(55, 200, 255, 255, maxColorValue = 255),
      rgb(0, 125, 125, 255, maxColorValue = 255),
      rgb(169, 169, 169, 255, maxColorValue = 255)
    ),
    col2 = c(
      rgb(0, 0, 255, 150, maxColorValue = 255),
      rgb(0, 50, 255, 150, maxColorValue = 255),
      rgb(0, 120, 255, 150, maxColorValue = 255),
      rgb(0, 170, 255, 150, maxColorValue = 255),
      rgb(200, 255, 80, 150, maxColorValue = 255),
      rgb(100, 255, 50, 150, maxColorValue = 255),
      rgb(50, 255, 25, 150, maxColorValue = 255),
      rgb(0, 255, 0, 150, maxColorValue = 255),
      rgb(255, 220, 100, 150, maxColorValue = 255),
      rgb(255, 0, 0, 150, maxColorValue = 255),
      rgb(255, 150, 150, 150, maxColorValue = 255),
      rgb(245, 165, 0, 150, maxColorValue = 255),
      rgb(255, 255, 0, 150, maxColorValue = 255),
      rgb(150, 255, 150, 150, maxColorValue = 255),
      rgb(100, 200, 100, 150, maxColorValue = 255),
      rgb(50, 200, 50, 150, maxColorValue = 255),
      rgb(170, 175, 255, 150, maxColorValue = 255),
      rgb(90, 120, 220, 150, maxColorValue = 255),
      rgb(75, 80, 180, 150, maxColorValue = 255),
      rgb(255, 0, 255, 150, maxColorValue = 255),
      rgb(200, 0, 200, 150, maxColorValue = 255),
      rgb(0, 255, 255, 150, maxColorValue = 255),
      rgb(55, 200, 255, 150, maxColorValue = 255),
      rgb(0, 125, 125, 150, maxColorValue = 255),
      rgb(169, 169, 169, 150, maxColorValue = 255)
    ),
    stringsAsFactors = F
  )

ecoregion.color.table2 <-
  list(
    ecoregion = list(
      ecoregion.1 = c("1", "2"),
      ecoregion.2 = c("3", "4"),
      ecoregion.3 = c("5"),
      ecoregion.4 = c("6", "7"),
      ecoregion.5 = c("8"),
      ecoregion.6 = c("9"),
      ecoregion.7 = c("11"),
      ecoregion.8 = c("10", "12", "13"),
      ecoregion.9 = c("14", "15", "16"),
      ecoregion.10 = c("17", "18", "19"),
      ecoregion.11 = c("20", "21"),
      ecoregion.12 = c("22", "23", "24"),
      ecoregion.13 = c("25")
    ),
    n.site = rep(NA, 13),
    n.site2 = rep(NA, 13),
    plot.text = c(
      "Tundra",
      "Taiga-Hudson plains",
      "Northern forest",
      "Northwestern forests",
      "Eastern forests",
      "Great plains",
      "Mediterranean California",
      "Desert-Semiarid highliands",
      "Tropical forests",
      "Andes",
      "Amazonian-Orinocan",
      "Pampas-Patagonian",
      "Hawaii"
    ),
    col = c(
      rgb(0, 0, 255, 255, maxColorValue = 255),
      rgb(0, 170, 255, 255, maxColorValue = 255),
      rgb(200, 255, 80, 255, maxColorValue = 255),
      rgb(150, 255, 50, 255, maxColorValue = 255),
      rgb(0, 255, 150, 255, maxColorValue = 255),
      rgb(255, 220, 100, 255, maxColorValue = 255),
      rgb(255, 150, 150, 255, maxColorValue = 255),
      rgb(255, 0, 0, 255, maxColorValue = 255),
      rgb(100, 200, 100, 255, maxColorValue = 255),
      rgb(170, 175, 255, 255, maxColorValue = 255),
      rgb(255, 0, 255, 255, maxColorValue = 255),
      rgb(0, 125, 125, 255, maxColorValue = 255),
      rgb(169, 169, 169, 255, maxColorValue = 255)
    ),
    col2 = c(
      rgb(0, 0, 255, 150, maxColorValue = 255),
      rgb(0, 170, 255, 150, maxColorValue = 255),
      rgb(200, 255, 80, 150, maxColorValue = 255),
      rgb(150, 255, 50, 150, maxColorValue = 255),
      rgb(0, 255, 150, 150, maxColorValue = 255),
      rgb(255, 220, 100, 150, maxColorValue = 255),
      rgb(255, 150, 150, 150, maxColorValue = 255),
      rgb(255, 0, 0, 150, maxColorValue = 255),
      rgb(100, 200, 100, 150, maxColorValue = 255),
      rgb(170, 175, 255, 150, maxColorValue = 255),
      rgb(255, 0, 255, 150, maxColorValue = 255),
      rgb(0, 125, 125, 150, maxColorValue = 255),
      rgb(169, 169, 169, 150, maxColorValue = 255)
    )
  )
