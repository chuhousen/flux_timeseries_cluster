## following a full list IGBP groups
igbp.ls <- list(
  igbp.b = c("EBF", "DBF"),
  igbp.m = c("MF"),
  igbp.n = c("ENF", "DNF"),
  igbp.s = c("SAV", "WSA"),
  igbp.r = c("OSH", "CSH"),
  igbp.g = c("GRA"),
  igbp.c = c("CRO", "CVM"),
  igbp.w = c("WET", "SNO", "WAT"),
  igbp.u = c("URB", "BSV")
)

igbp.text.ls <- c(
  "Broadleaf forest",
  "Mixed forest",
  "Needleleaf forest",
  "Savanna",
  "Shrubland",
  "Grassland",
  "Cropland",
  "Wetland / Ice / Water",
  "Urban / Barren"
)

igbp.color.table <-
  data.frame(
    igbp = c(
      "WAT",
      "ENF",
      "EBF",
      "DNF",
      "DBF",
      "MF",
      "CSH",
      "OSH",
      "WSA",
      "SAV",
      "GRA",
      "WET",
      "CRO",
      "URB",
      "CVM",
      "SNO",
      "BSV"
    ),
    plot.group = c(8, 3, 1, 3, 1, 2,
                   5, 5, 4, 4,
                   6, 8, 7, 9, 7, 8, 9),
    n.site = NA,
    n.site2 = NA,
    plot.order = c(14, 4, 1, 5, 2, 3,
                   8, 9, 6, 7,
                   10, 13, 11, 17, 12, 15, 16),
    plot.text = c(
      "Water",
      "Evergreen needleleaf forest",
      "Evergreen broadleaf forest",
      "Decidous needleleaf forest",
      "Decidous broadleaf forest",
      "Mixed forest",
      "Closed shrubland",
      "Open shrubland",
      "Woody savanna",
      "Savanna",
      "Grassland",
      "Wetland",
      "Cropland",
      "Urban",
      "Cropland-vegetation mosaics",
      "Snow/ice",
      "Barren land"
    ),
    col = c(
      rgb(55, 8, 255, 255, maxColorValue = 255),
      rgb(0, 141, 2, 255, maxColorValue = 255),
      rgb(0, 249, 10, 255, maxColorValue = 255),
      rgb(158, 208, 3, 255, maxColorValue = 255),
      rgb(206, 254, 212, 255, maxColorValue = 255),
      rgb(49, 162, 120, 255, maxColorValue = 255),
      rgb(168, 73, 116, 255, maxColorValue = 255),
      rgb(255, 211, 160, 255, maxColorValue = 255),
      rgb(202, 254, 143, 255, maxColorValue = 255),
      rgb(247, 214, 0, 255, maxColorValue = 255),
      rgb(255, 168, 0, 255, maxColorValue = 255),
      rgb(0, 140, 148, 255, maxColorValue = 255),
      rgb(250, 254, 3, 255, maxColorValue = 255),
      rgb(255, 42, 0, 255, maxColorValue = 255),
      rgb(138, 144, 0, 255, maxColorValue = 255),
      rgb(200, 200, 200, 255, maxColorValue = 255),
      rgb(170, 170, 170, 255, maxColorValue = 255)
    ),
    col2 = c(
      rgb(55, 8, 255, 150, maxColorValue = 255),
      rgb(0, 141, 2, 150, maxColorValue = 255),
      rgb(0, 249, 10, 150, maxColorValue = 255),
      rgb(158, 208, 3, 150, maxColorValue = 255),
      rgb(206, 254, 212, 150, maxColorValue = 255),
      rgb(49, 162, 120, 150, maxColorValue = 255),
      rgb(168, 73, 116, 150, maxColorValue = 255),
      rgb(255, 211, 160, 150, maxColorValue = 255),
      rgb(202, 254, 143, 150, maxColorValue = 255),
      rgb(247, 214, 0, 150, maxColorValue = 255),
      rgb(255, 168, 0, 150, maxColorValue = 255),
      rgb(0, 140, 148, 150, maxColorValue = 255),
      rgb(250, 254, 3, 150, maxColorValue = 255),
      rgb(255, 42, 0, 150, maxColorValue = 255),
      rgb(138, 144, 0, 150, maxColorValue = 255),
      rgb(200, 200, 200, 150, maxColorValue = 255),
      rgb(170, 170, 170, 150, maxColorValue = 255)
    ),
    stringsAsFactors = F
  )

## consolidated version 
igbp.color.table2 <-
  list(
    igbp = list(igbp.b = c("EBF", "DBF"),
                igbp.m = c("MF"),
                igbp.n = c("ENF", "DNF"),
                igbp.s = c("SAV", "WSA"),
                igbp.r = c("OSH", "CSH"),
                igbp.g = c("GRA"),
                igbp.c = c("CRO", "CVM"),
                igbp.w = c("WET", "SNO", "WAT"),
                igbp.u = c("URB", "BSV")),
    n.site = rep(NA, 9),
    n.site2 = rep(NA, 9),
    plot.text = c(
      "Broadleaf forest",
      "Mixed forest",
      "Needleleaf forest",
      "Savanna",
      "Shrubland",
      "Grassland",
      "Cropland",
      "Wetland / Ice / Water",
      "Urban / Barren"
    ),
    col = c(
      rgb(0, 249, 10, 255, maxColorValue = 255),
      rgb(49, 162, 120, 255, maxColorValue = 255),
      rgb(0, 141, 2, 255, maxColorValue = 255),
      rgb(247, 214, 0, 255, maxColorValue = 255),
      rgb(49, 162, 120, 255, maxColorValue = 255),
      rgb(255, 168, 0, 255, maxColorValue = 255),
      rgb(250, 254, 3, 255, maxColorValue = 255),
      rgb(55, 8, 255, 255, maxColorValue = 255),
      rgb(170, 170, 170, 255, maxColorValue = 255)
    ),
    col2 = c(
      rgb(0, 249, 10, 150, maxColorValue = 255),
      rgb(49, 162, 120, 150, maxColorValue = 255),
      rgb(0, 141, 2, 150, maxColorValue = 255),
      rgb(247, 214, 0, 150, maxColorValue = 255),
      rgb(49, 162, 120, 150, maxColorValue = 255),
      rgb(255, 168, 0, 150, maxColorValue = 255),
      rgb(250, 254, 3, 150, maxColorValue = 255),
      rgb(55, 8, 255, 150, maxColorValue = 255),
      rgb(170, 170, 170, 150, maxColorValue = 255)
    )
  )
