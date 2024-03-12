
  library(rgbif)
  library(ggplot2)
  library(ggspatial)
  library(spThin)
  library(sdmpredictors)
  library(ecospat)
  library(patchwork)
  library(raster)
  library(viridis)
  library(grid)
  library(terra)
  library(ENMTools)
  library(dplyr)
  library(grid)
  library(gridExtra)
  library(rgeos)
  library(rgdal)
  library(sf)
  library(ade4)
  library(factoextra)
  library(gridExtra)
  library(maps)
  library(ggpubr)
  library(magick)
  library(igraph)
  library(ENMTML)
  library(png)
  library(maptools)
  library(rlang)
  library(FactoMineR)
  
  
  
  # read occ records - present
  occ_cur <- read.table("occ/present.txt", header = TRUE)
  occ_cur$age <- "Present"
  
  # read occ records - mid Holocene
  occ_pal <- read.table("occ/paleo.txt", header = TRUE)
  occ_pal$age <- "Holocene"
  
  # Plot occ and background
  data("wrld_simpl")
  occ <- rbind.data.frame(occ_cur, occ_pal)
  
  
  v <- vect(cbind(occ$x, occ$y), crs="+proj=longlat")
  cv <- convHull(v)
  #### convert to spatialpolygondatafeam
  pl <- sf::st_as_sf(cv)
  
  coordinates(occ)=~x+y
  proj4string(occ)<- CRS("+proj=longlat +datum=WGS84")
  #### buffer occurrences enmtools
  sp_buf2 <- background.shape.buffer(occ, 1700000)
  #### transform buffer in SpatialPolygonsDataFrame
  sp_buf2 <- sp::spTransform(sp_buf2,CRS("+proj=longlat +ellps=WGS84"))
  #### # transform buffer in SpatialPolygonsDataFrame
  sp_buf2 <- SpatialPolygonsDataFrame(sp_buf2, data=as.data.frame("sp_buf2"))
  
  pl <- as_Spatial(pl$geometry)
  
  WorldData <- map_data('world')
  
  poli <- data.frame(x = c(-115, -115, -110, -110),
                     y = c(-10, -15, -15, -10), id = "poli")
  
  buf <- data.frame(x = c(-115, -115, -110, -110),
                    y = c(-18, -23, -23, -18), id = "buf")
  occ1 <- rbind.data.frame(occ_cur, occ_pal)
  
  
  ##0##
  
  map_bkg_pol <- ggplot() +
    
    geom_polygon(pl, mapping = aes(x=long, y=lat, group = group),
                 fill = "gray", alpha = 0,  col = "black", linewidth = 2.5,
                 linetype = 1) +
    geom_polygon(sp_buf2, mapping = aes(x=long, y=lat, group = group),
                 fill = "gray", alpha = 0,  col = "black",  linewidth = 2.5,
                 linetype = 14) +
    
    geom_polygon(WorldData, mapping = aes(x=long, y=lat, group = group),
                 fill = "gray", alpha = 1,  col = "black") +
    
    geom_point(occ1, mapping = aes(x = x, y = y, group = 1, color = as.factor(age)),
               size = 10) +
    
    scale_colour_manual(values = c("blue", "red")) +
    labs(color = "Occurrences") +
    
    geom_text() +
    
    
    geom_polygon(data = poli, aes(x = x, y = y, group = id),
                 stat = "identity", position = "identity",
                 rule = "evenodd",
                 fill ="transparent", colour = "black",
                 linewidth = 2) +
    annotate("text", x = -98, y =  -13, label = "Polygon background",
             size = 14) +
    
    geom_polygon(data = buf, aes(x = x, y = y, group = id),
                 stat = "identity", position = "identity",
                 rule = "evenodd",
                 fill ="transparent", colour = "black",
                 linewidth = 2, linetype = 2) +
    annotate("text", x = -99, y =  -20, label = "Buffer background",
             size = 14) +
    
    
    coord_sf(xlim = c(-120, -18), ylim = c(-35, 50), expand = FALSE) +
    
    
    theme_classic() +
    
    theme(legend.position = c(0.85,0.80),
          legend.key.size = unit(2, 'cm'), #change legend key size
          legend.key.height = unit(2, 'cm'), #change legend key height
          legend.key.width = unit(2, 'cm'), #change legend key width
          legend.title = element_text(size = 50), #change legend title font size
          legend.text = element_text(size = 45)) +
    labs(x = "Longitude", y = "Latitude") +
    theme(axis.text=element_text(size = 45),
          axis.title=element_text(size = 45))
  map_bkg_pol <- map_bkg_pol +
    annotation_scale(location = "bl", width_hint = 0.2) +
    
    ggspatial::annotation_north_arrow(
      location = "bl", which_north = "grid",
      pad_x = unit(0.2, "in"), pad_y = unit(0.3, "in"),
      style = ggspatial::north_arrow_nautical(
        fill = c("grey40", "white"),
        line_col = "grey20",
        text_family = "ArcherPro Book"))
  
  dir.create("results")
  dir.create("results/plots")
  
    ggsave(
    filename = "occ_bkg.png",
    plot = map_bkg_pol,
    device = "png",
    path = "results/plots",
    scale = 1,
    width = 9000,
    height = 7000,
    units = "px",
    dpi = 300,
    limitsize = TRUE,
    bg = NULL
  )
  
  
  ### extract values from occurrences points and studu area
  
  # load environmental layers
  # present layers
  env_pre_pat <- stack(list.files("var/present", pattern = "tif",
                                  full.names = TRUE, include.dirs = FALSE, all.files = FALSE,
                                  ignore.case = TRUE, recursive = FALSE))
  # mid holocene layers
  env_pal_pat <- stack(list.files("var/paleo", pattern = "tif",
                                  full.names = TRUE, include.dirs = FALSE, all.files = FALSE,
                                  ignore.case = TRUE, recursive = FALSE))
  #
  # extract values from occurrences records
  #### present data
  occ_cur_ext1 <- cbind(occ_cur[, 2:3], raster::extract(env_pre_pat, occ_cur[, 2:3]))
  occ_cur_ext <- cbind(occ_cur_ext1, occ_cur[1])
  occ_cur_ext <- na.omit(occ_cur_ext)
  occ_cur_ext$age <- "Present"
 
  # 
  #### Paleo data
  occ_pal_ext1 <- cbind(occ_pal[, 2:3], raster::extract(env_pal_pat, occ_pal[, 2:3]))
  occ_pal_ext <- cbind(occ_pal_ext1, occ_pal[1])
  occ_pal_ext <- na.omit(occ_pal_ext)
  occ_pal_ext$age <- "Holocene"
  
  ### union occurrence values
  occ_data_cur <- rbind.data.frame(occ_cur_ext, occ_pal_ext)
  occ_data_val <- occ_data_cur
  posi <- 12
  occ_data_val <- occ_data_val %>%
    select(-posi)
  #
  #
  # extract values from polygon area
  
  #### present data
  #### transform dataframe in spatial data
  v <- vect(cbind(occ_cur$x, occ_cur$y), crs="+proj=longlat")
  cv <- convHull(v)
  
  #### convert to spatialpolygondatafeam
  pol <- as(cv, "Spatial")
  
  #### crop layer
  env_cur_crop <- mask(env_pre_pat, pol)
  
  # extract values from M area
  M_ext_cur <- rasterToPoints(env_cur_crop, progress="text")
  M_ext_cur <- data.frame(M_ext_cur)
  M_ext_cur$age <- "Present"
  
  
  #### Paleo data
  #### transform dataframe in spatial data
  v2 <- vect(cbind(occ_pal$x, occ_pal$y), crs="+proj=longlat")
  cv2 <- convHull(v2)
  #### convert to spatialpolygondatafeam
  pl2 <- sf::st_as_sf(cv2)
  pol2 <- as(cv2, "Spatial")
  
  #### crop layer
  env_pal_crop <- mask(env_pal_pat, pl2)
  
  # extract values from M area
  M_ext_pal <- rasterToPoints(env_pal_crop, progress="text")
  M_ext_pal <- data.frame(M_ext_pal)
  M_ext_pal$age <- "Holocene"
  
  ### union background values
  bkg_data_cur <- rbind.data.frame(M_ext_cur, M_ext_pal)
  bkg_data_pol <- bkg_data_cur
  #
  #
  
  # extract values from buffer area
  #### present data
  # create a background area

  
  # read occ records - present
  occ_cur <- read.table("occ/present.txt", header = TRUE)
  occ_cur$age <- "Present"
  # read occ records - mid Holocene
  occ_pal <- read.table("occ/paleo.txt", header = TRUE)
  occ_pal$age <- "Holocene"
  
  
  # load environmental layers
  # present layers
  env_pre_pat <- stack(list.files("var/present", pattern = "tif",
                                  full.names = TRUE, include.dirs = FALSE, all.files = FALSE,
                                  ignore.case = TRUE, recursive = FALSE))
  # mid holocene layers
  env_pal_pat <- stack(list.files("var/paleo", pattern = "tif",
                                  full.names = TRUE, include.dirs = FALSE, all.files = FALSE,
                                  ignore.case = TRUE, recursive = FALSE))
  #
  # extract values from occurrences records
  #### present data
  occ_cur_ext1 <- cbind(occ_cur[, 2:3], raster::extract(env_pre_pat, occ_cur[, 2:3]))
  occ_cur_ext <- cbind(occ_cur_ext1, occ_cur[1])
  occ_cur_ext <- na.omit(occ_cur_ext)
  occ_cur_ext$age <- "Present"
  
  # 
  #### Paleo data
  occ_pal_ext1 <- cbind(occ_pal[, 2:3], raster::extract(env_pal_pat, occ_pal[, 2:3]))
  occ_pal_ext <- cbind(occ_pal_ext1, occ_pal[1])
  occ_pal_ext <- na.omit(occ_pal_ext)
  occ_pal_ext$age <- "Holocene"
  
  ### union occurrence values
  occ_data_cur <- rbind.data.frame(occ_cur_ext, occ_pal_ext)
  occ_data_val <- occ_data_cur
  posi <- 12
  occ_data_val <- occ_data_val %>%
    select(-posi)
  
  
  coordinates(occ_cur)=~x+y
  proj4string(occ_cur)<- CRS("+proj=longlat +datum=WGS84")
  #### buffer occurrences enmtools
  sp_buf <- background.shape.buffer(occ_cur, 1700000)
  #### transform buffer in SpatialPolygonsDataFrame
  sp_buf <- sp::spTransform(sp_buf,CRS("+proj=longlat +ellps=WGS84"))
  #### # transform buffer in SpatialPolygonsDataFrame
  sp_buf <- SpatialPolygonsDataFrame(sp_buf,data=as.data.frame("sp_buf"))
  
  #### crop layer
  env_cur_crop <- mask(env_pre_pat, sp_buf)
  
  # extract values from M area
  buf_cur <- rasterToPoints(env_cur_crop, progress="text")
  buf_cur <- data.frame(buf_cur)
  buf_cur$age <- "Present"
  
  
  #### paleo data
  #### paleo data
  # create a background area
    #### buffer occurrences enmtools
  coordinates(occ_pal)=~x+y
  proj4string(occ_pal)<- CRS("+proj=longlat +datum=WGS84")
  
  sp_buf2 <- background.shape.buffer(occ_pal, 1700000)
  #### transform buffer in SpatialPolygonsDataFrame
  sp_buf2 <- sp::spTransform(sp_buf2,CRS("+proj=longlat +ellps=WGS84"))
  #### # transform buffer in SpatialPolygonsDataFrame
  sp_buf2 <- SpatialPolygonsDataFrame(sp_buf2, data=as.data.frame("sp_buf2"))
  #### crop layer
  
  env_pal_crop <- mask(env_pal_pat, sp_buf2)
  
  
  # extract values from M area
  buf_pal <- rasterToPoints(env_pal_crop, progress="text")
  buf_pal <- data.frame(buf_pal)
  buf_pal$age <- "Holocene"
  
  ### union buffer values
  buf_data_pal <- rbind.data.frame(buf_cur, buf_pal)
  bkg_data_buf <- buf_data_pal
  #
  #
  #
  #
  #
  #
  ### shapiro test
  dir.create("results/shapiro")
  # present occurrences values
  write.table(occ_cur_ext[3:11],
              "results/shapiro/occ_cur_val.txt",
              sep = "\t", col.names = TRUE, row.names = TRUE, append = TRUE, quote = TRUE)
  
  
  write.table(occ_pal_ext[3:11],
              "results/shapiro/occ_pal_val.txt",
              sep = "\t", col.names = TRUE, row.names = TRUE, append = TRUE, quote = TRUE)
  
  ## Shapiro test
  env_shap <- occ_cur_ext[3:11]
  
  results <- list()
  
  for (i in colnames(env_shap)) {
    results_test <- shapiro.test(env_shap[[i]])
    results[[i]] <- results_test}
  
  shap <- matrix(results)
  write.csv(
    shap,
    "results/shapiro/occ_cur_shap.txt",
    sep = "\t", col.names = TRUE, row.names = TRUE, append = TRUE,
    quote = TRUE)
  
  
  # paleo occurrences values
  write.table(occ_pal_ext[3:11],
              "results/shapiro/occ_pal_val.txt",
              sep = "\t", col.names = TRUE, row.names = TRUE,
              append = TRUE, quote = TRUE)
  
  
  ## Shapiro test
  env_shap2 <- occ_pal_ext[3:11]
  
  results2 <- list()
  
  for (i in colnames(env_shap2)) {
    results_test2 <- shapiro.test(env_shap2[[i]])
    results2[[i]] <- results_test2}
  
  shap2 <- matrix(results2)      
  
  write.csv(
    shap2,
    "results/shapiro/occ_pal_shap.txt", sep = "\t", col.names = TRUE, row.names = TRUE, append = TRUE,
    quote = TRUE)
  #
  #
  ### t test
  ## occurrences values
  # salinity
  sss_mean_t_test <- data_frame(x = c(occ_pal_ext$sss_mean, occ_cur_ext$sss_mean),
                                y = c(occ_pal_ext$age, occ_cur_ext$age))
  sss_mean_boxplot <- ggboxplot(sss_mean_t_test, x = "y", y = "x",
                                color = "y", palette = c("blue", "red"),
                                add = "jitter", title = "Sea surface salinity mean",
                                xlab = "", ylab = "Salinity") +
    stat_compare_means(method = "t.test") +
    theme(plot.title = element_text(size = 10))
  #
  sss_min_t_test <- data_frame(x = c(occ_pal_ext$sss_min, occ_cur_ext$sss_min),
                               y = c(occ_pal_ext$age, occ_cur_ext$age))
  sss_min_boxplot <- ggboxplot(sss_min_t_test, x = "y", y = "x",
                               color = "y", palette = c("blue", "red"),
                               add = "jitter", title = "Sea surface salinity minimun",
                               xlab = "", ylab = "") +
    stat_compare_means(method = "t.test") +
    theme(plot.title = element_text(size = 10))
  #
  sss_max_t_test <- data_frame(x = c(occ_pal_ext$sss_max, occ_cur_ext$sss_max),
                               y = c(occ_pal_ext$age, occ_cur_ext$age))
  sss_max_boxplot <- ggboxplot(sss_max_t_test, x = "y", y = "x",
                               color = "y", palette = c("blue", "red"),
                               add = "jitter", title = "Sea surface salinity maximun",
                               xlab = "", ylab = "") +
    stat_compare_means(method = "t.test") +
    theme(plot.title = element_text(size = 10))
  #
  sss_range_t_test <- data_frame(x = c(occ_pal_ext$sss_range, occ_cur_ext$sss_range),
                                 y = c(occ_pal_ext$age, occ_cur_ext$age))
  sss_range_boxplot <- ggboxplot(sss_range_t_test, x = "y", y = "x",
                                 color = "y", palette = c("blue", "red"),
                                 add = "jitter", title = "Sea surface salinity range",
                                 xlab = "", ylab = "") +
    stat_compare_means(method = "t.test") +
    theme(plot.title = element_text(size = 10))
  ### temperature
  sst_mean_t_test <- data_frame(x = c(occ_pal_ext$sst_mean, occ_cur_ext$sst_mean),
                                y = c(occ_pal_ext$age, occ_cur_ext$age))
  sst_mean_boxplot <- ggboxplot(sst_mean_t_test, x = "y", y = "x",
                                color = "y", palette = c("blue", "red"),
                                add = "jitter", title = "Sea surface temperature mean",
                                xlab = "", ylab = "Temperature") +
    stat_compare_means(method = "t.test") +
    theme(plot.title = element_text(size = 10))
  #
  sst_min_t_test <- data_frame(x = c(occ_pal_ext$sst_mean, occ_cur_ext$sst_mean),
                               y = c(occ_pal_ext$age, occ_cur_ext$age))
  sst_min_boxplot <- ggboxplot(sst_min_t_test, x = "y", y = "x",
                               color = "y", palette = c("blue", "red"),
                               add = "jitter", title = "Sea surface temperature minimun",
                               xlab = "", ylab = "") +
    stat_compare_means(method = "t.test") +
    theme(plot.title = element_text(size = 10))
  #
  sst_max_t_test <- data_frame(x = c(occ_pal_ext$sst_max, occ_cur_ext$sst_max),
                               y = c(occ_pal_ext$age, occ_cur_ext$age))
  sst_max_boxplot <- ggboxplot(sst_max_t_test, x = "y", y = "x",
                               color = "y", palette = c("blue", "red"),
                               add = "jitter", title = "Sea surface temperature maximun",
                               xlab = "", ylab = "") +
    stat_compare_means(method = "t.test") +
    theme(plot.title = element_text(size = 10))
  #
  sst_range_t_test <- data_frame(x = c(occ_pal_ext$sst_range, occ_cur_ext$sst_range),
                                 y = c(occ_pal_ext$age, occ_cur_ext$age))
  sst_range_boxplot <- ggboxplot(sst_range_t_test, x = "y", y = "x",
                                 color = "y", palette = c("blue", "red"),
                                 add = "jitter", title = "Sea surface temperature range",
                                 xlab = "", ylab = "") +
    stat_compare_means(method = "t.test") +
    theme(plot.title = element_text(size = 10))
  #
  bat_t_test <- data_frame(x = c(occ_pal_ext$bat, occ_cur_ext$bat),
                                 y = c(occ_pal_ext$age, occ_cur_ext$age))
  bat_boxplot <- ggboxplot(bat_t_test, x = "y", y = "x",
                                 color = "y", palette = c("blue", "red"),
                                 add = "jitter", title = "Sea bathymetry",
                                 xlab = "", ylab = "") +
    stat_compare_means(method = "t.test") +
    theme(plot.title = element_text(size = 10))
  
  
  png("results/plots/shapiro.png",
      units="px", width= 7000, height = 9000, res = 800)
  (bat_boxplot + sss_mean_boxplot + sss_min_boxplot) /
    (sss_max_boxplot + sss_range_boxplot + sst_mean_boxplot) /
    (sst_min_boxplot + sst_max_boxplot + sst_range_boxplot)
  
  dev.off()
  #
  #
  #
  ## ECOSPAT MPC
  
  # list of scenaries
  cen.lst <- levels(occ_data_val[, 12])
  cen.lst <- c("Present", "Holocene")
   
  sp.nbocc2 <- c()
  
  for (i in 1:length(cen.lst)) {
    sp.nbocc2 <- c(sp.nbocc2, length(which(occ_data_val[, 12] == cen.lst[i])))
  }
  
  # calculate the nb of occurences per species
  nb.sp2 <- length(cen.lst) # nb of species
  
  # try with all and then try only worldclim Variables
  Xvar <- c(3:11)
  nvar <- length(Xvar)
  # number of interation for the tests of equivalency and similarity
  iterations <- 100
  # resolution of the gridding of the climate space
  R <- 100
  
  #################################### PCA-ENVIRONMENT ##################################
  ### poygon data
  
  data_ext <- rbind(occ_data_val[, Xvar], bkg_data_pol[, Xvar])
  data_ext_pol <- na.omit(data_ext)
  data_ext_pol <- data_ext
  # tirei o row.w = w2
  #
  ### PCA polygon
  pca.cal2 <- dudi.pca(data_ext_pol, center = TRUE, scale = TRUE, scannf = FALSE, nf = 2)
  
  dir.create("results/PCA")
  # write PCA results
  write.csv(pca.cal2$co, "results/PCA/var_contribution_pol.csv")
  
  # cumulative.variance.percent
  cumul_pol1 <- get_eig(pca.cal2)
  cumul_pol <- cumul_pol1$cumulative.variance.percent
  write.csv(cumul_pol, "results/PCA/cumulative_variance_pol.csv")
  
  #
  #
  #
  ###### ecospat ######
  
  sp.combn2 <- combn(1:2, 2)
  
  for (i in 1:ncol(sp.combn2)) {
    row.pal <- which(occ_data_val[, 12] == cen.lst[sp.combn2[2, i]]) # rows in data corresponding to paleo
    row.cur <- which(occ_data_val[, 12] == cen.lst[sp.combn2[1, i]]) # rows in data corresponding to current
    name.pal <- cen.lst[sp.combn2[2, i]]
    name.cur <- cen.lst[sp.combn2[1, i]]
    # predict the scores on the axes
    scores.clim2 <- pca.cal2$li[(nrow(occ_data_val) + 1):nrow(data_ext_pol), ] # scores for global climate
    scores.pal2 <- pca.cal2$li[row.pal, ] # scores for paleo
    scores.cur2 <- pca.cal2$li[row.cur, ] # scores for current
  }
  
  #  calculation of occurence density and test of niche equivalency and similarity
  # with the default kernel method
  za1 <- ecospat.grid.clim.dyn(scores.clim2, scores.clim2, scores.pal2)
  za2 <- ecospat.grid.clim.dyn(scores.clim2, scores.clim2, scores.cur2)
  
  
  dir.create("results/overlap")
  # niche overlap test
  overl_pol <- ecospat.niche.overlap (za1, za2, cor = FALSE)
  write.csv(overl_pol, "results/overlap/overl_pol.csv")
  #
  #
  #
  #
  #
  ## ECOSPAT BUFFER
  # list of scenaries
  cen.lst <- levels(bkg_data_buf[, 12])
  cen.lst <- c("Present", "Holocene")
  
  sp.nbocc3 <- c()
  
  for (i in 1:length(cen.lst)) {
    sp.nbocc3 <- c(sp.nbocc3, length(which(bkg_data_buf[, 12] == cen.lst[i])))
  }
  
  # calculate the nb of occurences per species
  nb.sp3 <- length(cen.lst) # nb of species
  
  # try with all and then try only worldclim Variables
  Xvar2 <- c(3:11)
  nvar <- length(Xvar2)
  # number of interation for the tests of equivalency and similarity
  iterations <- 100
  # resolution of the gridding of the climate space
  R <- 100
  
  #################################### PCA-ENVIRONMENT ##################################
  ### buffer data
  data_ext_buf <- rbind(occ_data_val[, Xvar2], bkg_data_buf[, Xvar2])
  data_ext_buf <- na.omit(data_ext_buf)
  ### PCA buffer
  pca_buf <- dudi.pca(data_ext_buf, center = TRUE, scale = TRUE, scannf = FALSE, nf = 2)
  pca_buf$li[2] <- pca_buf$li[2] * -1
  #     pca_buf$li$Axis2 <- pca_buf$li$Axis2 * -1
  
  # write PCA results
  write.csv(pca_buf$co, "results/PCA/var_contribution_buf.csv")
  # cumulative.variance.percent
  cumul_buf1 <- get_eig(pca_buf)
  cumul_buf <- cumul_buf1$cumulative.variance.percent
  write.csv(cumul_buf, "results/PCA/cumulative_variance_buf.csv")
  
  
  ###### selection of species ######
  sp.combn3 <- combn(1:2, 2)
  
  for (i in 1:ncol(sp.combn3)) {
    row.pal3 <- which(occ_data_val[, 12] == cen.lst[sp.combn3[2, i]]) # rows in data corresponding to paleo
    row.cur3 <- which(occ_data_val[, 12] == cen.lst[sp.combn3[1, i]]) # rows in data corresponding to current
    name.pal3 <- cen.lst[sp.combn3[2, i]]
    name.cur3 <- cen.lst[sp.combn3[1, i]]
    # predict the scores on the axes
    scores.clim3 <- pca_buf$li[(nrow(occ_data_val) + 1):nrow(bkg_data_buf), ] # scores for global climate
    scores.pal3 <- pca_buf$li[row.pal3, ] # scores for paleo
    scores.cur3 <- pca_buf$li[row.cur3, ] # scores for current
  }
  
  #  calculation of occurence density and test of niche equivalency and similarity
  # with the default kernel method
  za3 <- ecospat.grid.clim.dyn(scores.clim3,
                               scores.clim3, scores.pal3, R = 100)
  za4 <- ecospat.grid.clim.dyn(scores.clim3,
                               scores.clim3, scores.cur3, R = 100)
  
  #
  
  # niche overlap test
  overl_buf <- ecospat.niche.overlap (za3, za4, cor = FALSE)
  write.csv(overl_buf, "results/overlap/overl_buf.csv")
  
  ###############
  # PLOT
  
  ## PCA
  xlab_pol =  round(cumul_pol[1], 2)
  ylab_pol <- round(cumul_pol[2]-cumul_pol[1], 2)
  
  
  pca_pol_plot <- fviz_pca_biplot(pca.cal2, label ="var", col.ind="cos2", col.var = "black",
                                  gradient.cols = c("red", "green", "blue"),
                                  geom = c("point", "text"),
                                  ggtheme = theme_minimal(),
                                  labelsize = 3,
                                  title = "Polygon background",
                                  xlab = paste0("PC1 = ", xlab_pol, "%"),
                                  ylab = paste0("PC1 = ", ylab_pol, "%"),
                                  legend = "none") +
    annotate("text", x = -12, y = -22, label = "(a)", size = 6)
  
  
  
  xlab_buf =  round(cumul_buf[1], 2)
  ylab_buf <- round(cumul_buf[2]-cumul_buf[1], 2)
  pca_buf_plot <- fviz_pca_biplot(pca_buf, label ="var", col.ind="cos2",
                                  col.var = "black",                                
                                  gradient.cols = c("red", "green", "blue"),
                                  geom = c("point", "text"),
                                  title = "Buffer background",
                                  labelsize = 3,
                                  ggtheme = theme_minimal(),
                                  xlab = paste0("PC1 = ", xlab_buf, "%"),
                                  ylab = paste0("PC1 = ", ylab_buf, "%"),
                                  legend = "none") +
    
    annotate("text", x = -12, y = -22, label = "(b)", size = 6)
  #
  
  pca_pol_buf <- grid.arrange(pca_pol_plot, pca_buf_plot, ncol =2)
  ### save plot PCA
  
  ggsave(
    filename = "pca_pol_buf.png",
    plot = pca_pol_buf,
    device = "png",
    path = "results/plots",
    scale = 1,
    width = 4000,
    height = 2000,
    units = "px",
    dpi = 300,
    limitsize = TRUE,
    bg = NULL
  )
  
  ### #### kernel density
  #### density polygon present
  za1_lyr <- raster(za1$z)
  za1_plot <- ggplot(za1_lyr) +
    geom_raster(
      mapping = aes(x = x, y = y, fill = lyr.1),
      stat = "identity",
      position = "identity",
      hjust = 0.5,
      vjust = 0.5,
      na.rm = FALSE) +
    scale_fill_viridis(name = "", option="turbo",
                       begin = 0.0) +
    ggtitle("Holocene (polygon background)") +
    labs(x = "PC1", y = "PC2") +
    theme_minimal() +
    annotate("text", x = 15, y = -27, label = "(c)", size = 6,
             color = "white") +
    theme(plot.title = element_text(size = 10)) +
    theme(legend.position = "none") +
    ylim(c(-20,3))
  
  #### density polygon holocene
  za2_lyr <- raster(za2$z)
  za2_plot <- ggplot(za2_lyr) +
    geom_raster(
      mapping = aes(x = x, y = y, fill = lyr.1),
      stat = "identity",
      position = "identity",
      hjust = 0.5,
      vjust = 0.5,
      na.rm = FALSE) +
    scale_fill_viridis(name = "", option="turbo",
                       begin = 0.0) +
    
    ggtitle("Present (polygon background)") +
    labs(x = "PC1", y = "PC2") +
    annotate("text", x = 15, y = -27, label = "(d)", size = 6,
             color = "white") +
    theme_minimal()  +
    theme(plot.title = element_text(size = 10)) +
    theme(legend.position = "none") +
    ylim(c(-20,3))
  
  #### density buffer present
  za3_lyr <- raster(za3$z)
  za3_plot <- ggplot(za3_lyr) +
    geom_raster(
      mapping = aes(x = x, y = y, fill = lyr.1),
      stat = "identity",
      position = "identity",
      hjust = 0.5,
      vjust = 0.5,
      na.rm = FALSE) +
    scale_fill_viridis(name = "", option="turbo",
                       begin = 0.0) +
    ggtitle("Holocene (buffer background)") +
    labs(x = "PC1", y = "PC2") +
    xlim(c(-2,10)) +

    annotate("text", x = 8, y = -18, label = "(e)", size = 6,
             color = "white") +
    theme_minimal()  +
    theme(plot.title = element_text(size = 10)) +
    theme(legend.position = "none")
    
  #### density buffer holocene
  za4_lyr <- raster(za4$z)
  za4_plot <- ggplot(za4_lyr) +
    geom_raster(
      mapping = aes(x = x, y = y, fill = lyr.1),
      stat = "identity",
      position = "identity",
      hjust = 0.5,
      vjust = 0.5,
      na.rm = FALSE) +
    scale_fill_viridis(name = "", option="turbo",
                       begin = 0.0) +
    ggtitle("Present (buffer background)") +
    labs(x = "PC1", y = "PC2") +
    xlim(c(-2,10)) +
   
    annotate("text", x = 8, y = -18, label = "(f)", size = 6,
             color = "white") +
    
    theme_minimal() +
    theme(plot.title = element_text(size = 10)) +
    theme(legend.position = "none")
  
  
  png("results/plots/kernel.png",
      units="px", width= 9000, height = 6000, res = 800)
  
  (za1_plot | za2_plot | za3_plot | za4_plot)
  dev.off()
  
  
  png("results/plots/pca_kernel.png",
      units="px", width= 9000, height = 6000, res = 800)
  
  (pca_pol_plot | pca_buf_plot) /
    (za1_plot | za2_plot | za3_plot | za4_plot)
  dev.off()
  
  
  #  Niche index
  #
  niche_index <- ecospat.niche.dyn.index(za1, za2, intersection = NA)
  ext_pol <- niche_index$dynamic.index.w[1]
  ext_pol <- round(ext_pol[1], 2)
  
  sta_pol <- niche_index$dynamic.index.w[2]
  sta_pol <- round(sta_pol, 2)
  
  unf_pol <- niche_index$dynamic.index.w[3]
  unf_pol <- round(unf_pol, 2) 
  
  
  
  niche_index_buf <- ecospat.niche.dyn.index (za3, za4, intersection = NA)
  ext_buf <- niche_index_buf$dynamic.index.w[1]
  ext_buf <- round(ext_buf, 2)
  
  sta_buf <- niche_index_buf$dynamic.index.w[2]
  sta_buf <- round(sta_buf, 2)
  
  unf_buf <- niche_index_buf$dynamic.index.w[3]
  unf_buf <- round(unf_buf, 2)
  

    ### Plot niche dynamic
  
  ### Plot niche dynamic
  #### Overlap polygon
  over_pol_res = round(overl_pol[[1]], 2)
   png(
    filename = "results/plots/over_pol.png",
    width = 5000, height = 5000, units = "px", res = 800)
  
  #### Niche overlap polygon
  
  ecospat.plot.niche.dyn(za1, za2, quant = 0.9,
                         title = "Polygon background",
                         name.axis1 = "PC1",
                         name.axis2 = "PC2",
                         interest = 1, col.exp = "orange",
                         col.stab = "blue", col.unf = "red",
                         colZ1 = "brown", colZ2 = "black",
                         transparency = 40)
  
  
  ecospat.shift.centroids(scores.cur2, scores.pal2, scores.clim2, scores.clim2)
  
  text(paste0("Overlap = ", over_pol_res), x = 14, y = 1, cex = 1)
  text(paste0("Expansion = ", ext_pol), x = 14, y = -0.4, cex = 1, col = "orange")
  text(paste0("Loss = ", unf_pol), x = 14, y = -2, cex = 1, col = "red")
  text(paste0("Stability = ", sta_pol), x = 14, y = -3.6, cex = 1, col = "blue")
  text("(g)", x = -4, y = -20, cex = 1.3, col = "black")
  dev.off()
  
  
  ##### Overlap buffer
  overl_buf_res =  round(overl_buf[[1]], 2)
 
   png(
    filename = "results/plots/over_buf.png",
    width = 5000, height = 5000, units = "px", res = 800)
   #### Niche overlap polygon
   
  ecospat.plot.niche.dyn(za3, za4, quant = 0.9,
                         title = "Buffer background",
                         name.axis1 = "PC1",
                         name.axis2 = "PC2",
                         interest = 1, col.exp = "orange", col.stab = "blue",
                         col.unf = "red", colZ1 = "brown", colZ2 = "black",
                         transparency = 40,
                         xlim = c(-5, 15))
  
 
  ecospat.shift.centroids(scores.cur3, scores.pal3, scores.clim3, scores.clim3)
  
  text(paste0("Overlap = ", overl_buf_res), x = 1, y = -12, cex = 1, col = "black")
  text(paste0("Expansion = ", ext_buf), x = 1, y = -13.4, cex = 1, col = "orange")
  text(paste0("Loss = ", unf_buf), x = 1, y = -14.8, cex = 1, col = "red")
  text(paste0("Stability = ", sta_buf), x = 1, y = -16, cex = 1, col = "blue")
  text("(h)", x = -1, y = -19, cex = 1.3, col = "black")
  dev.off()
  
  
  # Carregar as duas imagens (substitua os caminhos pelos caminhos das suas imagens)
  image_path1 <- "results/plots/over_pol.png"
  image_path2 <- "results/plots/over_buf.png"
  # Carregar as imagens usando a função image_read()
  img1 <- image_read(image_path1)
  img2 <- image_read(image_path2)
  # Unir as imagens verticalmente usando a função image_append()
  imagem_unida <- image_append(c(img1, img2), stack = FALSE)
  
  # Salvar a imagem resultante (opcional)
  output_path <- "results/plots/over_union.png"
  image_write(imagem_unida, path = output_path)
  
  
  # Carregar as duas imagens (substitua os caminhos pelos caminhos das suas imagens)
  image_path3 <- "results/plots/pca_kernel.png"
  image_path4 <- "results/plots/over_union.png"
  # Carregar as imagens usando a função image_read()
  img3 <- image_read(image_path3)
  img4 <- image_read(image_path4)
  # Unir as imagens verticalmente usando a função image_append()
  imagem_unida2 <- image_append(c(img3, img4), stack = TRUE)
  
  # Salvar a imagem resultante (opcional)
  output_path <- "results/plots/pca_kernel_over.png"
  image_write(imagem_unida2, path = output_path)

  ###########################################################
  ###########################################################
  ###########################################################
  ###########################################################
  ###########################################################
  ###########################################################
  ###########################################################
  ###########################################################
  ###########################################################
  ###########################################################
  ###########################################################
  ###########################################################
  ###########################################################
  
  ### Pre modeling process
  
  #### Model present polygon background
  
  library(ENMTML)
  library(raster)
  library(kuenm)
  
  dir.create("ENMTML/models")
  dir.create("ENMTML/models/mod1_pres_pol")
  dir.create("ENMTML/models/mod1_pres_pol/env")
  dir.create("ENMTML/models/mod1_pres_pol/proj")
  dir.create("ENMTML/models/mod1_pres_pol/proj/world")
  
  dir.create("ENMTML/models/mod2_pres_buf")
  dir.create("ENMTML/models/mod2_pres_buf/env")
  dir.create("ENMTML/models/mod2_pres_buf/proj")
  dir.create("ENMTML/models/mod2_pres_buf/proj/world")
  
  
  dir.create("ENMTML/models/mod3_pal_pol")
  dir.create("ENMTML/models/mod3_pal_pol/env")
  dir.create("ENMTML/models/mod3_pal_pol/proj")
  dir.create("ENMTML/models/mod3_pal_pol/proj/world")
  
  dir.create("ENMTML/models/mod4_pal_buf")
  dir.create("ENMTML/models/mod4_pal_buf/env")
  dir.create("ENMTML/models/mod4_pal_buf/proj")
  dir.create("ENMTML/models/mod4_pal_buf/proj/world")
  
  
  ##### Present layers
  ##### Criate a PCA layers
  env_pres <- stack(list.files("var/present",
                          pattern = "tif", full.names = TRUE,
                          include.dirs = FALSE, all.files = FALSE,
                          ignore.case = TRUE, recursive = FALSE))
  
  kuenm_rpca(env_pres, in.format = "tif", var.scale = TRUE,
             write.result = TRUE, out.format = "GTiff",
             out.dir = "var/ENMTML/present/pca", project = FALSE, n.pcs = 9)  
  ##### load pca envs
  env_pca_pres <- stack(list.files("var/ENMTML/present/pca/initial",
                   pattern = "tif", full.names = TRUE,
                   include.dirs = FALSE, all.files = FALSE,
                   ignore.case = TRUE, recursive = FALSE))

  
  env_pal <- stack(list.files("var/paleo",
                               pattern = "tif", full.names = TRUE,
                               include.dirs = FALSE, all.files = FALSE,
                               ignore.case = TRUE, recursive = FALSE))
  kuenm_rpca(env_pal, in.format = "tif", var.scale = TRUE,
             write.result = TRUE, out.format = "GTiff",
             out.dir = "var/ENMTML/paleo/pca", project = FALSE, n.pcs = 9)
  
  env_pca_pal <- stack(list.files("var/ENMTML/paleo/pca/initial",
                                   pattern = "tif", full.names = TRUE,
                                   include.dirs = FALSE, all.files = FALSE,
                                   ignore.case = TRUE, recursive = FALSE))
  
  ##### Crop layers
  # read occ records - present
  occ_cur <- read.table("occ/present.txt", header = TRUE)
  occ_cur$age <- "Present"
  
  # read occ records - mid Holocene
  occ_pal <- read.table("occ/paleo.txt", header = TRUE)
  occ_pal$age <- "Holocene"
  occ <- rbind.data.frame(occ_cur, occ_pal)
  
  ##### polygon background
  library(terra)
  v <- vect(cbind(occ$x, occ$y), crs="+proj=longlat")
  cv <- convHull(v)
  #### convert to spatialpolygondatafeam
  pl <- sf::st_as_sf(cv)
  env_pres_pol <- mask(env_pca_pres, pl)
  env_pal_pol <- mask(env_pca_pal, pl)
  
  ##### buffer background
  library(ENMTools)
  coordinates(occ)=~x+y
  proj4string(occ)<- CRS("+proj=longlat +datum=WGS84")
  #### buffer occurrences enmtools
  sp_buf <- background.shape.buffer(occ, 1700000)
  #### transform buffer in SpatialPolygonsDataFrame
  sp_buf <- sp::spTransform(sp_buf,CRS("+proj=longlat +ellps=WGS84"))
  #### # transform buffer in SpatialPolygonsDataFrame
  sp_buf <- SpatialPolygonsDataFrame(sp_buf, data=as.data.frame("sp_buf"))
  
  ##### crop pca envs
  env_pres_buf <- mask(env_pca_pres, sp_buf)
  env_pal_buf <- mask(env_pca_pal, sp_buf)
  ##### write crop pca present polygon
  writeRaster(env_pres_pol[[1]], "ENMTML/models/mod1_pres_pol/env/PC1.tif",
              overwrite = TRUE)
  writeRaster(env_pres_pol[[2]], "ENMTML/models/mod1_pres_pol/env/PC2.tif",
              overwrite = TRUE)
  writeRaster(env_pca_pres[[1]], "ENMTML/models/mod1_pres_pol/proj/world/PC1.tif",
              overwrite = TRUE)
  writeRaster(env_pres_pol[[2]], "ENMTML/models/mod1_pres_pol/proj/world/PC2.tif",
              overwrite = TRUE)
  ##### write crop pca present buffer
  writeRaster(env_pres_buf[[1]], "ENMTML/models/mod2_pres_buf/env/PC1.tif",
              overwrite = TRUE)
  writeRaster(env_pres_buf[[2]], "ENMTML/models/mod2_pres_buf/env/PC2.tif",
              overwrite = TRUE)
  writeRaster(env_pca_pres[[1]], "ENMTML/models/mod2_pres_buf/proj/world/PC1.tif",
              overwrite = TRUE)
  writeRaster(env_pca_pres[[2]], "ENMTML/models/mod2_pres_buf/proj/world/PC2.tif",
              overwrite = TRUE)
  
  
  ##### write crop pca Holocene polygon
  writeRaster(env_pal_pol[[1]], "ENMTML/models/mod3_pal_pol/env/PC1.tif",
              overwrite = TRUE)
  writeRaster(env_pal_pol[[2]], "ENMTML/models/mod3_pal_pol/env/PC2.tif",
              overwrite = TRUE)
  writeRaster(env_pca_pal[[1]], "ENMTML/models/mod3_pal_pol/proj/world/PC1.tif",
              overwrite = TRUE)
  writeRaster(env_pca_pal[[2]], "ENMTML/models/mod3_pal_pol/proj/world/PC2.tif",
              overwrite = TRUE)
  
  ##### write crop pca Holocene buffer
  writeRaster(env_pal_buf[[1]], "ENMTML/models/mod4_pal_buf/env/PC1.tif",
              overwrite = TRUE)
  writeRaster(env_pal_buf[[2]], "ENMTML/models/mod4_pal_buf/env/PC2.tif",
              overwrite = TRUE)
  writeRaster(env_pca_pal[[1]], "ENMTML/models/mod4_pal_buf/proj/world/PC1.tif",
              overwrite = TRUE)
  writeRaster(env_pca_pal[[2]], "ENMTML/models/mod4_pal_buf/proj/world/PC2.tif",
              overwrite = TRUE)

  
  
  ##### Filter occurrences
  occ_thined <- thin(loc.data = occ_cur, 
                     lat.col = "y", long.col = "x",
                     spec.col = "sp", 
                     thin.par = 10,
                     reps = 1,
                     max.files = 1,
                     out.dir = "ENMTML/models/mod1_pres_pol",
                     out.base = "present",
                     locs.thinned.list.return = TRUE,
                     write.files = TRUE,
                     write.log.file = FALSE,
                     log.file = "present.txt",
                     verbose = TRUE)
  
  occ_cur <- read.csv("ENMTML/models/mod1_pres_pol/present_thin1.csv", header = TRUE)
  
  write.table(occ_cur, "ENMTML/models/mod1_pres_pol/occ.txt",
              sep = "\t", row.names = FALSE)
  write.table(occ_cur, "ENMTML/models/mod2_pres_buf/occ.txt",
              sep = "\t", row.names = FALSE)
  
  ##### paleo occurrences
  occ_thined <- thin(loc.data = occ_pal, 
                     lat.col = "y", long.col = "x",
                     spec.col = "sp", 
                     thin.par = 10,
                     reps = 1,
                     max.files = 1,
                     out.dir = "ENMTML/models/mod3_pal_pol",
                     out.base = "paleo",
                     locs.thinned.list.return = TRUE,
                     write.files = TRUE,
                     write.log.file = FALSE,
                     log.file = "paleo.txt",
                     verbose = TRUE)
  
  occ_pal <- read.csv("ENMTML/models/mod3_pal_pol/paleo_thin1.csv", header = TRUE)
  
  write.table(occ_pal, "ENMTML/models/mod3_pal_pol/occ.txt",
              sep = "\t", row.names = FALSE)
  write.table(occ_pal, "ENMTML/models/mod4_pal_buf/occ.txt",
              sep = "\t", row.names = FALSE)
  
  
  ### Modelling process
  
  library(ENMTML)
  library(terra)
  # PRESENT
  # polygon background
  #### mod1_pres_pol
  ENMTML(
    pred_dir = "ENMTML/models/mod1_pres_pol/env",
    proj_dir = "ENMTML/models/mod1_pres_pol/proj",
    result_dir = "ENMTML/models/mod1_pres_pol/results",
    occ_file = "ENMTML/models/mod1_pres_pol/occ.txt",
    sp = 'sp',
    x = 'x',
    y = 'y',
    min_occ = 2,
    thin_occ = c(method='USER-DEFINED', distance='10'),
    eval_occ = NULL,
    colin_var = NULL,
    imp_var = TRUE,
    sp_accessible_area = NULL,
    pseudoabs_method=c(method='GEO_ENV_KM_CONST', width='50'),
    pres_abs_ratio = 1,
    part = c(method='BOOT', replicates='100', proportion='0.7'),
    save_part = FALSE,
    save_final = TRUE,
    algorithm = "MXD",
    thr=c(type='MAX_TSS'),
    msdm=NULL,
    ensemble = NULL,
    extrapolation=TRUE,
    cores = 1)
    #
  
  #
  #### mod2_pres_buf
  ENMTML(
  pred_dir = "ENMTML/models/mod2_pres_buf/env",
  proj_dir = "ENMTML/models/mod2_pres_buf/proj",
  result_dir = "ENMTML/models/mod2_pres_buf/results",
  occ_file = "ENMTML/models/mod2_pres_buf/occ.txt",
  sp = 'sp',
  x = 'x',
  y = 'y',
  min_occ = 2,
  thin_occ = c(method='USER-DEFINED', distance='10'),
  eval_occ = NULL,
  colin_var = NULL,
  imp_var = TRUE,
  sp_accessible_area = NULL,
  pseudoabs_method=c(method='GEO_ENV_KM_CONST', width='50'),
  pres_abs_ratio = 1,
  part = c(method='BOOT', replicates='100', proportion='0.7'),
  save_part = FALSE,
  save_final = TRUE,
  algorithm =  "MXD",
  thr=c(type='MAX_TSS'),
  msdm=NULL,
  ensemble = NULL,
  extrapolation=TRUE,
  cores = 1)
  #
  
  #
  #
  # PALEO
  # polygon background
  #### mod3_pal_pol
  
  ENMTML(
    pred_dir = "ENMTML/models/mod3_pal_pol/env",
    proj_dir = "ENMTML/models/mod3_pal_pol/proj",
    result_dir = "ENMTML/models/mod3_pal_pol/results",
    occ_file = "ENMTML/models/mod3_pal_pol/occ.txt",
    sp = 'sp',
    x = 'x',
    y = 'y',
    min_occ = 2,
    thin_occ = c(method='USER-DEFINED', distance='10'),
    eval_occ = NULL,
    colin_var = NULL,
    imp_var = TRUE,
    sp_accessible_area = NULL,
    pseudoabs_method=c(method='GEO_ENV_KM_CONST', width='50'),
    pres_abs_ratio = 1,
    part = c(method='BOOT', replicates='100', proportion='0.7'),
    save_part = FALSE,
    save_final = TRUE,
    algorithm = "MXD",
    thr=c(type='MAX_TSS'),
    msdm=NULL,
    ensemble = NULL,
    extrapolation=TRUE,
    cores = 1)
  #
  # buff background
  #### mod4_pal_buf
  
  ENMTML(
    pred_dir = "ENMTML/models/mod4_pal_buf/env",
    proj_dir = "ENMTML/models/mod4_pal_buf/proj",
    result_dir = "ENMTML/models/mod4_pal_buf/results",
    occ_file = "ENMTML/models/mod4_pal_buf/occ.txt",
    sp = 'sp',
    x = 'x',
    y = 'y',
    min_occ = 2,
    thin_occ = c(method='USER-DEFINED', distance='10'),
    eval_occ = NULL,
    colin_var = NULL,
    imp_var = TRUE,
    sp_accessible_area = NULL,
    pseudoabs_method=c(method='GEO_ENV_KM_CONST', width='50'),
    pres_abs_ratio = 1,
    part = c(method='BOOT', replicates='100', proportion='0.7'),
    save_part = FALSE,
    save_final = TRUE,
    algorithm = "MXD" ,
    thr=c(type='MAX_TSS'),
    msdm=NULL,
    ensemble = NULL,
    extrapolation=TRUE,
    cores = 1)  
  
  ##################################################################
  ##################################################################
  ##################################################################
  ##################################################################
  ##################################################################
  ################
  ################
  
  #### Results (AUC)
  mod1_resul <- read.table("ENMTML/models/mod1_pres_pol/results/Evaluation_Table.txt", header = TRUE)
  mod1_AUC <- round(mod1_resul[1,5], 2)
  
  mod2_resul <- read.table("ENMTML/models/mod2_pres_buf/results/Evaluation_Table.txt", header = TRUE)
  mod2_AUC <- round(mod2_resul[1,5], 2)
  
  mod3_resul <- read.table("ENMTML/models/mod3_pal_pol/results/Evaluation_Table.txt", header = TRUE)
  mod3_AUC <- round(mod3_resul[1,5], 2)
  
  mod4_resul <- read.table("ENMTML/models/mod4_pal_buf/results/Evaluation_Table.txt", header = TRUE)
  mod4_AUC <- round(mod4_resul[1,5], 2)
  #### Geographical niche dynamic
  
  library(raster)
  mod1_pol_pres <- raster("ENMTML/models/mod1_pres_pol/results/Projection/world/MXD/p_astreoides.tif")
  
  mod2_buf_pres <- raster("ENMTML/models/mod2_pres_buf/results/Projection/world/MXD/p_astreoides.tif")
  
  mod3_pal_pol <- raster("ENMTML/models/mod3_pal_pol/results/Projection/world/MXD/p_astreoides.tif")
  
  mod4_pal_buf <- raster("ENMTML/models/mod4_pal_buf/results/Projection/world/MXD/p_astreoides.tif")
  
  ##### 
  
  # read occ records - present
  occ_cur <- read.table("occ/present.txt", header = TRUE)
  occ_cur$age <- "Present"
  # read occ records - mid Holocene
  
  occ_pal <- read.table("occ/paleo.txt", header = TRUE)
  occ_pal$age <- "Holocene"
  
  occ <- rbind.data.frame(occ_cur, occ_pal)
  
  # Plot occ and background
  WorldData <- map_data('world')
  
  v <- vect(cbind(occ$x, occ$y), crs="+proj=longlat")
  cv <- convHull(v)
  pol <- sf::st_as_sf(cv)
  pol <- as(cv, "Spatial")
  proj4string(pol)<- CRS("+proj=longlat +datum=WGS84")
  pol <- SpatialPolygonsDataFrame(pol,data=as.data.frame("pol"))
  
  
  occ_cur <- read.table("occ/present.txt", header = TRUE)
  occ_cur$age <- "Present"
  
  occ_pal <- read.table("occ/paleo.txt", header = TRUE)
  occ_pal$age <- "Holocene"
  occ <- rbind.data.frame(occ_cur, occ_pal)
  
  coordinates(occ)=~x+y
  proj4string(occ)<- CRS("+proj=longlat +datum=WGS84")
  #### buffer occurrences enmtools
  sp_buf2 <- background.shape.buffer(occ, 1700000)
  #### transform buffer in SpatialPolygonsDataFrame
  sp_buf2 <- sp::spTransform(sp_buf2,CRS("+proj=longlat +ellps=WGS84"))
  #### # transform buffer in SpatialPolygonsDataFrame
  sp_buf2 <- SpatialPolygonsDataFrame(sp_buf2, data=as.data.frame("sp_buf2"))
  
  
  #### Plot continuous models
  
  mod1_plot <- ggplot() +
    geom_raster(mod1_pol_pres, mapping = aes(x = x, y = y, fill = p_astreoides),
                stat = "identity",
                position = "identity",
                hjust = 0.5,
                vjust = 0.5,
                na.rm = FALSE) +
    scale_fill_viridis(name = "Adequability", option="turbo",
                       begin = 0.1) +
    
    geom_polygon(pol, mapping = aes(x=long, y=lat, group = 1),
                 fill = "transparent", alpha = 0.001,  col = "black", size = 1) +
    
    geom_polygon(WorldData, mapping = aes(x=long, y=lat, group = group),
                 fill = "gray", alpha = 1,  col = "gray1",
                 linewidth = .2) +
    labs(x = "Longitude", y = "Latitude") +
    
    ggtitle("Polygon background (Present)") +
    
    annotate("text", x = c(-105, -90), y =  c(-5, -5),
             label = c("AUC = ", mod1_AUC), size = 4,
             color = "white") +
    annotate("text", x = -23, y =  -33,
             label = "b)", size = 6,
             color = "black") +
    
    theme_classic() +
    
    theme(axis.text=element_text(size = 8),
          axis.title=element_text(size = 8)) +
    theme(legend.position = c(0.85, 0.8),
          legend.key.size = unit(0.1, 'cm'), #change legend key size
          legend.key.height = unit(0.15, 'cm'), #change legend key height
          legend.key.width = unit(0.3, 'cm'), #change legend key width
          legend.spacing.y = unit(0.3, 'cm'),
          legend.text = element_text(size = 7),
          legend.title = element_text(size = 7)) +
    theme(legend.key = element_rect(fill = "blue")) +
  coord_sf(xlim = c(-113, -22), ylim = c(48, -33)) 
    
    mod1_plot1 <- mod1_plot +
    ggspatial::annotation_north_arrow(
      location = "bl", which_north = "grid",
      pad_x = unit(0.2, "in"), pad_y = unit(0.3, "in"),
      style = ggspatial::north_arrow_nautical(
        fill = c("grey40", "white"),
        line_col = "grey20",
        text_family = "ArcherPro Book"))

  ##### mod 2
  mod2_plot <- ggplot() +
    geom_raster(mod2_buf_pres, mapping = aes(x = x, y = y, fill = p_astreoides),
                stat = "identity",
                position = "identity",
                hjust = 0.5,
                vjust = 0.5,
                na.rm = FALSE) +
    scale_fill_viridis(name = "Adequability", option="turbo",
                       begin = 0.1) +
    
    geom_polygon(sp_buf2, mapping = aes(x=long, y=lat, group = group),
                 fill = "gray", alpha = 0,  col = "black",  linewidth = 1.5,
                 linetype = 14) +
    
    geom_polygon(WorldData, mapping = aes(x=long, y=lat, group = group),
                 fill = "gray", alpha = 1,  col = "gray1",
                 linewidth = .2) +
    labs(x = "Longitude", y = "Latitude") +
    
    ggtitle("Buffer background (Present)") +
    
    annotate("text", x = c(-108, -93), y =  c(-8, -8),
             label = c("AUC = ", mod2_AUC), size = 4,
             color = "white") +
    annotate("text", x = -23, y =  -33,
             label = "b)", size = 6,
             color = "black") +
    
    theme_classic() +
    
    theme(axis.text=element_text(size = 8),
          axis.title=element_text(size = 8)) +
    theme(legend.position = c(0.85, 0.8),
          legend.key.size = unit(0.1, 'cm'), #change legend key size
          legend.key.height = unit(0.15, 'cm'), #change legend key height
          legend.key.width = unit(0.3, 'cm'), #change legend key width
          legend.spacing.y = unit(0.3, 'cm'),
          legend.text = element_text(size = 7),
          legend.title = element_text(size = 7)) +
    theme(legend.key = element_rect(fill = "blue")) +
    coord_sf(xlim = c(-113, -22), ylim = c(48, -33))
  
  mod2_plot1 <- mod2_plot +
    annotation_scale(location = "bl", width_hint = 0.2) +
    ggspatial::annotation_north_arrow(
      location = "bl", which_north = "grid",
      pad_x = unit(0.2, "in"), pad_y = unit(0.3, "in"),
      style = ggspatial::north_arrow_nautical(
        fill = c("grey40", "white"),
        line_col = "grey20",
        text_family = "ArcherPro Book"))
  
  ##### mod3
  mod3_plot <- ggplot() +
    geom_raster(mod3_pal_pol, mapping = aes(x = x, y = y, fill = p_astreoides),
                stat = "identity",
                position = "identity",
                hjust = 0.5,
                vjust = 0.5,
                na.rm = FALSE) +
    scale_fill_viridis(name = "Adequability", option="turbo",
                       begin = 0.1) +
    
    geom_polygon(pol, mapping = aes(x=long, y=lat, group = 1),
                 fill = "transparent", alpha = 0.001,  col = "black", size = 1) +
    
    geom_polygon(WorldData, mapping = aes(x=long, y=lat, group = group),
                 fill = "gray", alpha = 1,  col = "gray1",
                 linewidth = .2) +
    labs(x = "Longitude", y = "Latitude") +
    
    ggtitle("Polygon background (Holocene)") +
    
    annotate("text", x = c(-105, -90), y =  c(-5, -5),
             label = c("AUC = ", mod1_AUC), size = 4,
             color = "white") +
    annotate("text", x = -23, y =  -33,
             label = "b)", size = 6,
             color = "black") +
    
    theme_classic() +
    
    theme(axis.text=element_text(size = 8),
          axis.title=element_text(size = 8)) +
    theme(legend.position = c(0.85, 0.8),
          legend.key.size = unit(0.1, 'cm'), #change legend key size
          legend.key.height = unit(0.15, 'cm'), #change legend key height
          legend.key.width = unit(0.3, 'cm'), #change legend key width
          legend.spacing.y = unit(0.3, 'cm'),
          legend.text = element_text(size = 7),
          legend.title = element_text(size = 7)) +
    theme(legend.key = element_rect(fill = "blue")) +
    coord_sf(xlim = c(-113, -22), ylim = c(48, -33))
  
  mod3_plot1 <- mod3_plot +
    annotation_scale(location = "bl", width_hint = 0.2) +
    ggspatial::annotation_north_arrow(
      location = "bl", which_north = "grid",
      pad_x = unit(0.2, "in"), pad_y = unit(0.3, "in"),
      style = ggspatial::north_arrow_nautical(
        fill = c("grey40", "white"),
        line_col = "grey20",
        text_family = "ArcherPro Book"))
  
  # mod4
  mod4_plot <- ggplot() +
    geom_raster(mod4_pal_buf, mapping = aes(x = x, y = y, fill = p_astreoides),
                stat = "identity",
                position = "identity",
                hjust = 0.5,
                vjust = 0.5,
                na.rm = FALSE) +
    scale_fill_viridis(name = "Adequability", option="turbo",
                       begin = 0.1) +
    
    geom_polygon(sp_buf2, mapping = aes(x=long, y=lat, group = group),
                 fill = "gray", alpha = 0,  col = "black",  linewidth = 1.5,
                 linetype = 14) +
    
    geom_polygon(WorldData, mapping = aes(x=long, y=lat, group = group),
                 fill = "gray", alpha = 1,  col = "gray1",
                 linewidth = .2) +
    labs(x = "Longitude", y = "Latitude") +
    
    ggtitle("Polygon background (Holocene)") +
    
    annotate("text", x = c(-108, -93), y =  c(-8, -8),
             label = c("AUC = ", mod1_AUC), size = 4,
             color = "white") +
    annotate("text", x = -23, y =  -33,
             label = "b)", size = 6,
             color = "black") +
    
    theme_classic() +
    
    theme(axis.text=element_text(size = 8),
          axis.title=element_text(size = 8)) +
    theme(legend.position = c(0.85, 0.8),
          legend.key.size = unit(0.1, 'cm'), #change legend key size
          legend.key.height = unit(0.15, 'cm'), #change legend key height
          legend.key.width = unit(0.3, 'cm'), #change legend key width
          legend.spacing.y = unit(0.3, 'cm'),
          legend.text = element_text(size = 7),
          legend.title = element_text(size = 7)) +
    theme(legend.key = element_rect(fill = "blue")) +
    coord_sf(xlim = c(-113, -22), ylim = c(48, -33))
  
  mod4_plot1 <- mod4_plot +
    annotation_scale(location = "bl", width_hint = 0.2) +
    ggspatial::annotation_north_arrow(
      location = "bl", which_north = "grid",
      pad_x = unit(0.2, "in"), pad_y = unit(0.3, "in"),
      style = ggspatial::north_arrow_nautical(
        fill = c("grey40", "white"),
        line_col = "grey20",
        text_family = "ArcherPro Book"))
  
  
  png("results/plots/cont_models.png",
      units="px", width= 9000, height = 6000, res = 800)
  
   (mod3_plot1 + mod1_plot1) /
                      (mod4_plot1 + mod2_plot1)
  dev.off()
  ################################################################
  
  ### Binary models
  
  library(raster)
  mod1_bin <- raster("ENMTML/models/mod1_pres_pol/results/Projection/world/MXD/MAX_TSS/p_astreoides.tif")
  mod2_bin <- raster("ENMTML/models/mod2_pres_buf/results/Projection/world/MXD/MAX_TSS/p_astreoides.tif")
  mod3_bin <- raster("ENMTML/models/mod3_pal_pol/results/Projection/world/MXD/MAX_TSS/p_astreoides.tif")
  mod4_bin <- raster("ENMTML/models/mod4_pal_buf/results/Projection/world/MXD/MAX_TSS/p_astreoides.tif")

  #### crop binary models
  mod1_crop <- mask(mod1_bin, pol)
  mod2_crop <- mask(mod2_bin, sp_buf2)
  mod3_crop <- mask(mod3_bin, pol)
  mod4_crop <- mask(mod4_bin, sp_buf2)

  
  #### Binary map to presence map
  mod1_bin_01 <- rasterToPolygons(mod1_crop, fun=function(x){x == 1},
                                  n=4, na.rm=TRUE, digits=12, dissolve=TRUE)
  str(mod1_bin_01, max.level = 3)
  mod_pres_pol_pl1 <- mod1_bin_01@polygons
  mod_pres_pol_pl <- SpatialPolygons(mod_pres_pol_pl1, proj4string = mod1_bin_01@proj4string)
  
  #####  mod 2
  mod2_bin_01 <- rasterToPolygons(mod2_crop, fun=function(x){x == 1},
                                  n=4, na.rm=TRUE, digits=12, dissolve=TRUE)
  str(mod2_bin_01, max.level = 3)
  mod2_pres_buf_pl1  <- mod2_bin_01@polygons
  mod2_pres_buf_pl <- SpatialPolygons(mod2_pres_buf_pl1, proj4string = mod2_bin_01@proj4string)
  
  #####  mod3
  mod3_bin_01 <- rasterToPolygons(mod3_crop, fun=function(x){x == 1},
                                  n=4, na.rm=TRUE, digits=12, dissolve=TRUE)
  str(mod3_bin_01, max.level = 3)
  mod_pal_pol_pl1  <- mod3_bin_01@polygons
  mod_pal_pol_pl <- SpatialPolygons(mod_pal_pol_pl1, proj4string = mod3_bin_01@proj4string)

  ##### mod4
  mod4_bin_01 <- rasterToPolygons(mod4_crop, fun=function(x){x == 1},
                                  n=4, na.rm=TRUE, digits=12, dissolve=TRUE)
  str(mod4_bin_01, max.level = 3)
  mod4_pal_buf_pl1  <- mod4_bin_01@polygons
  mod4_pal_buf_pl <- SpatialPolygons(mod4_pal_buf_pl1, proj4string = mod4_bin_01@proj4string)

  
  ### total area of models
  mod_pres_pol_area <- gArea(mod_pres_pol_pl)
  mod_pal_pol_area <- gArea(mod_pal_pol_pl)
  mod2_pres_buf_area <- gArea(mod2_pres_buf_pl)
  mod4_pal_buf_area <- gArea(mod4_pal_buf_pl)  
  
  
  ### areas POLYGON
  ### Intersection (stability)  polygon paleo - present
  ##### Total areas
  area_pres_pol <- gArea(mod_pres_pol_pl)
  area_pal_pol <- gArea(mod_pal_pol_pl)
  
  ##### percentage area
  inters_pol <- gIntersection(mod_pal_pol_pl, mod_pres_pol_pl, byid = TRUE)
  inters_pol_area <- gArea(inters_pol)
  inters_pol_perc <- inters_pol_area / area_pres_pol * 100
  
  ##### Stability spatialpolygon
  mod_pres_pol_pl_st <- sf::st_as_sf(mod_pres_pol_pl)
  mod_pal_pol_pl_st <- sf::st_as_sf(mod_pal_pol_pl)
  stability <- st_intersection(mod_pal_pol_pl_st, mod_pres_pol_pl_st)
  stability_pol <- sf::as_Spatial(st_geometry(stability), IDs =  as.character(1:nrow(stability)))
  
  ### Difference (expansion) polygon paleo  * present
  exp_area_pol <- gDifference(mod_pres_pol_pl, mod_pal_pol_pl, byid = TRUE)
  exp_area_pol1 <- gArea(exp_area_pol)
  exp_area_pol2 <- gArea(mod_pres_pol_pl)
  # Calculando a porcentagem de expansion
  exp_area_pol_perc <- (exp_area_pol1 / exp_area_pol2) * 100
  

  ### Difference (lost) polygon paleo - present
  ## Lost area
  los_area_pol <- gDifference(mod_pal_pol_pl, mod_pres_pol_pl, byid = TRUE)
  los_area_pol1 <- gArea(los_area_pol)
  los_area_pol2 <- gArea(mod_pal_pol_pl)
  
  los_area_pol_perc <- los_area_pol1 / los_area_pol2 * 100
  
  
  #### Niche overlap
  mod1_over <- mask(mod1_pol_pres, pol)
  mod1_spat <- as(mod1_over, "SpatRaster")
  mod2_over <- mask(mod2_buf_pres, sp_buf2)
  mod2_spat <- as(mod2_over, "SpatRaster")
  mod3_over <- mask(mod3_pal_pol, pol)
  mod3_spat <- as(mod3_over, "SpatRaster")
  mod4_over <- mask(mod4_pal_buf, sp_buf2)
  mod4_spat <- as(mod4_over, "SpatRaster")
  
  over_pol <- raster.overlap(mod1_spat, mod3_spat, verbose = FALSE)
  over_buf <- raster.overlap(mod2_spat, mod4_spat, verbose = FALSE)
  
  ############
  
  plot_niche_din_pol <- ggplot() +
    geom_polygon(pol, mapping = aes(x=long, y=lat, group = 1),
                 fill = "transparent", alpha = 0.001,  col = "black", size = 1) +
    
    geom_polygon(mod_pal_pol_pl, mapping = aes(long, lat, group = group),
                 fill = "transparent", color = "transparent",
                 size = 1) +
    geom_polygon(mod_pres_pol_pl, mapping = aes(long, lat, group = group),
                 fill = "transparent", color = "transparent",
                 size = 1) +
    geom_polygon(stability_pol, mapping = aes(long, lat, group = group),
                 fill = "blue", alpha = 0.6) +
    geom_polygon(los_area_pol, mapping = aes(long, lat, group = group),
                 fill = "red", alpha = 0.6) +
    geom_polygon(exp_area_pol, mapping = aes(long, lat, group = group),
                 fill = "orange", alpha = 0.6) +
    
    
    geom_polygon(WorldData, mapping = aes(x=long, y=lat, group = group),
                 fill = "gray", alpha = 1,  col = "gray1",
                 linewidth = .2) +
    
    annotate("text", x = c(-38, -23.8), y =  c(41, 41),
             label = c("Overlap = ", round(over_pol$D, 2)), size = 4.5,
             color = "black") +
    
    annotate("text", x = c(-41, -23), y =  c(36, 36),
             label = c("Expansion = ", round(exp_area_pol_perc, 2)), size = 4.5,
             color = "orange") +
    annotate("text", x = c(-37, -25.1), y =  c(32, 32),
             label = c("Loss = ", round(los_area_pol_perc, 2)), size = 4.5,
             color = "red") +
    
    annotate("text", x = c(-41, -24), y =  c(27, 27),
             label = c("Stability = ", round(inters_pol_perc, 2)), size = 4.5,
             color = "blue") +
    
    annotate("text", x = -23, y =  -33,
             label = "c)", size = 6,
             color = "black") +
    theme_minimal() +
    labs(x = "Latitude", y = "Longitude") +

    coord_sf(xlim = c(-113, -22), ylim = c(48, -33)) +
    theme_classic()
  
  plot_niche_din_pol1 <- plot_niche_din_pol +
    
    annotation_scale(crs = 3995) +
    
    ggspatial::annotation_north_arrow(
      location = "bl", which_north = "grid",
      pad_x = unit(0.2, "in"), pad_y = unit(0.3, "in"),
      style = ggspatial::north_arrow_nautical(
        fill = c("grey40", "white"),
        line_col = "grey20",
        text_family = "ArcherPro Book")) 

  
  #######################
  ###### Buffer niche dynamic
  
  ### Intersection (stability)  buffer paleo - present
  ##### Total areas
  area_pres_buf <- gArea(mod2_pres_buf_pl)
  area_pal_buf <- gArea(mod4_pal_buf_pl)
  
  ##### percentage area
  inters_buf <- gIntersection(mod4_pal_buf_pl, mod2_pres_buf_pl, byid = TRUE)
  inters_buf_area <- gArea(inters_buf)
  inters_buf_perc <- inters_buf_area / area_pres_buf * 100
  
  ##### Stability spatialpolygon
  mod_pres_buf_pl_st <- sf::st_as_sf(mod2_pres_buf_pl)
  mod_pal_buf_pl_st <- sf::st_as_sf(mod4_pal_buf_pl)
  stability_buf <- st_intersection(mod_pal_buf_pl_st, mod_pres_buf_pl_st)
  stability_buf <- sf::as_Spatial(st_geometry(stability_buf), IDs =  as.character(1:nrow(stability_buf)))
  
  
  ### Difference (expansion) buffer paleo  * present
  exp_area_buf <- gDifference(mod2_pres_buf_pl, mod4_pal_buf_pl, byid = TRUE)
  exp_area_buf1 <- gArea(exp_area_buf)
  exp_area_buf2 <- gArea(mod2_pres_buf_pl)
  # Calculating the expansion percentage
  exp_area_buf_perc <- (exp_area_buf1 / exp_area_buf2) * 100
  
  
  ### Difference (lost) buffer paleo - present
  ## Lost area
  los_area_buf <- gDifference(mod4_pal_buf_pl, mod2_pres_buf_pl, byid = TRUE)
  los_area_buf1 <- gArea(los_area_buf)
  los_area_buf2 <- gArea(mod2_pres_buf_pl)
  los_area_buf_perc <- los_area_buf1 / los_area_buf2 * 100
  
  ##########
  plot_niche_din_buf <- ggplot() +
    
    geom_polygon(sp_buf2, mapping = aes(x=long, y=lat, group = group),
                 fill = "gray", alpha = 0,  col = "black",  linewidth = 1.5,
                 linetype = 14) +
    
    geom_polygon(stability_buf, mapping = aes(long, lat, group = group),
                 fill = "blue", alpha = 0.6) +
    geom_polygon(los_area_buf, mapping = aes(long, lat, group = group),
                 fill = "red", alpha = 0.6) +
    geom_polygon(exp_area_buf, mapping = aes(long, lat, group = group),
                 fill = "orange", alpha = 0.6) +
    
    geom_polygon(WorldData, mapping = aes(x=long, y=lat, group = group),
                 fill = "gray", alpha = 1,  col = "gray1",
                 linewidth = .2) +
    coord_sf(xlim = c(-113, -22), ylim = c(48, -33)) +
    
    annotate("text", x = c(-38, -23.8), y =  c(41, 41),
             label = c("Overlap = ", round(over_buf$D, 2)), size = 4.5,
             color = "black") +
    
    annotate("text", x = c(-41, -23), y =  c(36, 36),
             label = c("Expansion = ", round(exp_area_buf_perc, 2)), size = 4.5,
             color = "orange") +
    annotate("text", x = c(-37, -25.1), y =  c(32, 32),
             label = c("Loss = ", round(los_area_buf_perc, 2)), size = 4.5,
             color = "red") +
    
    annotate("text", x = c(-41, -24), y =  c(27, 27),
             label = c("Stability = ", round(inters_buf_perc, 2)), size = 4.5,
             color = "blue") +
    
    annotate("text", x = -23, y =  -33,
             label = "c)", size = 6,
             color = "black") +
    theme_classic()
  
  plot_niche_din_buf1 <- plot_niche_din_buf +
    
    annotation_scale(crs = 3995) +
    
    ggspatial::annotation_north_arrow(
      location = "bl", which_north = "grid",
      pad_x = unit(0.2, "in"), pad_y = unit(0.3, "in"),
      style = ggspatial::north_arrow_nautical(
        fill = c("grey40", "white"),
        line_col = "grey20",
        text_family = "ArcherPro Book")) 
  
  
  png("results/plots/mod_overlap_ok.png",
      units="px", width= 12000, height = 9000, res = 800)
  (mod3_plot + mod1_plot | plot_niche_din_pol1) /
    (mod4_plot + mod2_plot | plot_niche_din_buf1)
  dev.off()
  