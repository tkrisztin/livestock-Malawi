rm(list=ls(all=TRUE)) 

### PACKAGES
require(tidyverse)
require(readxl)
require(stringr)
require(forcats)
require(scales)
require(RColorBrewer)
# Spatial packages
require(ggmap)
require(raster)
require(rasterVis)
require(sp)
require(sf)
require(fmesher)
require(mapproj)
require(proj4)
require(gdalUtilities)
# Additional packages
require(fields)
require(splancs)
library(gridExtra)
require(exactextractr)

if (!require(INLA) | !require(inlabru)) {
  install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
} else {
  require(INLA)
  require(inlabru)
}

### Projections
crs <- "+init=EPSG:4326" # WSG84
crs_m <- "+init=EPSG:32736" # Zone 26 for Malawi

RUN_YEAR = "2019"
# RUN_YEAR = "2010"


### LOAD DATA AND MAPS
adm <- readRDS("./data/2010/maps/gaul/adm_2010_MWI.rds")
adm <- spTransform(adm, CRS(crs_m)) 
adm_withIsland = adm
adm_islands = adm[adm$adm_gaul == "LIKOMA",]
adm = adm[adm$adm_gaul != "LIKOMA",]

# Livestock data
lvst_raw = read_rds(paste0("./data/MWI/Household_survey/mwi_lvst_geo_",RUN_YEAR,".rds"))
lvst_raw <- lvst_raw %>%
  group_by(type, ea_id, lat, lon) %>%
  summarize(quantity = sum(quantity, na.rm = T))

lvst_dis = read_rds(paste0("./data/MWI/Household_survey/mwi_lvst_dis_",RUN_YEAR,".rds")) %>% 
  dplyr::select(-year) %>% 
  pivot_wider(names_from = type,values_from = quantity,values_fill = 0)

# travel_time
tt <- raster(file.path("./data/2010/maps/travel_time/travel_time_30sec_2010_MWI.tif"))
names(tt) <- "tt"
tt <- projectRaster(from = tt, crs = crs_m)
## Create an upper bound for travel time (3 hours)
ttemp = getValues(tt); ttemp[ttemp > 180] = 180; 
tt = setValues(tt,log1p(ttemp))
# Aridity Index
ai <- raster(file.path(".", "data/2010/maps/ai/ai_30sec_2010_MWI.tif"))
ai =  setValues(ai,log1p(getValues(ai))) 
names(ai) <- "ai"
ai <- projectRaster(ai, crs = crs_m)
# elevation
el <- raster(file.path(".", "data/2010/maps/elevation/elevation_30sec_2010_MWI.tif"))
el =  setValues(el,log1p(getValues(el)))
names(el) <- "el"
el <- projectRaster(el, crs = crs_m)
# pop  density
pop <- raster(paste0("./data/2010/maps/pop/pop_",RUN_YEAR,"_raster_30sec.tif"))
names(pop) <- "pop"
# input 0 values for population density where missing
pop_val = getValues(pop)
pop_val[which(is.na(pop_val) & !is.na(getValues(el)))] = 0
pop = setValues(pop,log1p(pop_val))
# road buffer
road1  <- raster("./data/2010/maps/roads/road_GRIP_raster_30sec_1km.tif")
road2p5  <- raster("./data/2010/maps/roads/road_GRIP_raster_30sec_2.5km.tif")
road5  <- raster("./data/2010/maps/roads/road_GRIP_raster_30sec_5km.tif")
names(road1) <- "road1"
names(road2p5) <- "road2p5"
names(road5) <- "road5"
# protected areas
protected  <- raster(paste0("./data/2010/maps/wdpa/wdpa_",RUN_YEAR,"_raster_30sec.tif"))
names(protected) <- "protected"
# tsetse flies
tsetse  <- raster("./data/2010/maps/tsetse/tsetse_raster_30sec.tif")
names(tsetse) <- "tsetse"
# grassland share
grass  <- raster(paste0("./data/2010/maps/landcover/grassland_",RUN_YEAR,"_raster_30sec.tif"))
names(grass) <- "grass"
# cropland share
crop  <- raster(paste0("./data/2010/maps/landcover/cropland_",RUN_YEAR,"_raster_30sec.tif"))
names(crop) <- "crop"
# forest share
forest  <- raster(paste0("./data/2010/maps/landcover/forest_",RUN_YEAR,"_raster_30sec.tif"))
names(forest) <- "forest"
# urban share
urban  <- raster(paste0("./data/2010/maps/landcover/settlement_",RUN_YEAR,"_raster_30sec.tif"))
names(urban) <- "urban"


### attribute pixels to regions ### 
pxl_ids = getValues(el)
pxl_ids[!is.na(pxl_ids)] = 1:sum(!is.na(pxl_ids))
pxls = setValues(el,pxl_ids)
pxls_list = exact_extract(pxls,adm)
adm_pxl = matrix(0,max(pxl_ids,na.rm = TRUE),nrow(adm))
for (rr in 1:length(pxls_list)) {
  non_na_ids = pxls_list[[rr]]$value
  adm_pxl[non_na_ids[!is.na(non_na_ids)],rr] = pxls_list[[rr]]$coverage_fraction[!is.na(non_na_ids)]
}
adm_pxl2 = apply(adm_pxl,c(1),
                 function(x) {if (all(x == 0)) {return(0)} else {return(which(x == max(x))[1])}})
pxl_ids[!is.na(pxl_ids)] = adm_pxl2
pxls_adm = setValues(el,pxl_ids)



### PREPARE LVST
# Spread data and add 0 for NA
lvst <- lvst_raw %>%
  spread(type, quantity, fill = 0) %>%
  mutate(x = lon, y = lat)
# Make spatial, reproject and rescale
coordinates(lvst) <- c("x", "y")
proj4string(lvst) <- crs
lvst <- spTransform(lvst, CRS(crs_m))

# Create file with locations of lvst observations and divide by 1000 so everything is in km and computation goes faster
lvst_df <- as.data.frame(lvst) %>%
  mutate(x = coords.x1/1000, y=coords.x2/1000)

coords <- dplyr::select(lvst_df, x, y) %>% 
  as.matrix


### PREPARE OUTPUT AND INPUT DATA
# Stack data
input_r <- stack(tt, ai, el, pop,road2p5,protected,tsetse,
                 grass,crop,forest,urban,road1,road5)

expl_vars = c("tt","ai","el","pop","road2p5","protected","tsetse",
              "grass","crop","forest","urban")
expl_vars_name = c("Travel time","Aridity","Elevation","Population",
                   "Road buffer",
                   "Protected areas","Tsetse species",
                   "Grassland","Cropland","Forest","Settlement")
expl_vars_plots = list()

for (ii in 1:length(expl_vars)) {
  expl_vars_plots[[ii]] = ggplot() +
    geom_raster(data = raster(input_r,layer = ii) %>%
                  as.data.frame(xy = TRUE) %>% 
                  rename_with(~"value", .cols = expl_vars[ii]), aes(x = x, y = y, fill = value)) +
    geom_sf(data = st_as_sf(adm), fill = NA, color = "grey80") +
    ggtitle("",subtitle = expl_vars_name[ii]) + 
    scale_fill_viridis_c(option = "A",na.value = NA,direction = -1) + labs(fill = expl_vars_name[ii]) + 
    theme_bw() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          #legend.position = "bottom",
          legend.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank())
}

scl = 1
pdf(paste0("./expl_vars_",RUN_YEAR,".pdf"), height = 16.5, width = 11.75, onefile = TRUE)
grid.arrange(expl_vars_plots[[1]] , expl_vars_plots[[2]] ,
               expl_vars_plots[[3]] , expl_vars_plots[[4]] ,
               expl_vars_plots[[5]], expl_vars_plots[[6]], 
             expl_vars_plots[[7]], expl_vars_plots[[8]],
             expl_vars_plots[[9]],expl_vars_plots[[10]],expl_vars_plots[[11]], 
             ncol = 4,padding = 0)
dev.off()

# Extract covariates for DHS
data <- as.data.frame(raster::extract(input_r, lvst)) %>%
  bind_cols(lvst_df) 
isOK_data = complete.cases(data) #& !data$protected
data = data[isOK_data,]
coords = coords[isOK_data,]

# Union  adm level 
border = as(st_union(st_as_sf(adm)),"Spatial")
island_border = as(st_union(st_as_sf(adm_islands)),"Spatial")
border2 = spTransform(border, CRS(crs_m))

# Create boundary and divide by 1000 to get kilometers
bdry = fm_as_segm(border)
bdry$loc[,1] <- sapply(bdry$loc[,1], function (x) x/1000) 
bdry$loc[,2] <- sapply(bdry$loc[,2], function (x) x/1000)

# Create mesh
# mesh_b <- inla.mesh.2d(boundary = bdry, 
#                        max.edge=c(8, 60), # maximum triange edge length for inner and outer boundary (in SAME SCALE UNIT as projection, here km!!)
#                        cutoff = 15) # Minimum allowed distance between points. Helps avoiding very small triangles at the border!
mesh_b <- fm_mesh_2d_inla(boundary = bdry,
                      max.edge=c(8, 60), # maximum triange edge length for inner and outer boundary (in SAME SCALE UNIT as projection, here km!!)
                      cutoff = 15) # Minimum allowed distance between points. Helps avoiding very small triangles at the border!
mesh_b$n # Number of nodes (not)
plot(mesh_b, main = "")
axis(1)
axis(2)
points(coords, pch = 19, col =2)


### MAPPING BETWEEN MESH AND CONTINOUS SPACE
# Projector matrix for coordinates of DSH values
A <- inla.spde.make.A(mesh = mesh_b, loc = coords)

## Build for prediction
# build grid
# Create grid
x_res <- round(dim(pop)[2]/1)
y_res <- round(dim(pop)[1]/1)
box = bbox(pop)
seq.x.grid = seq(from = box[1,1],to = box[1,2],
                 length=x_res)/1000
seq.y.grid = seq(from = box[2,1], to = box[2,2],
                 length=y_res)/1000
pred_grid = expand.grid(seq.x.grid,seq.y.grid)
AA = SpatialPoints(pred_grid*1000,proj4string = CRS(crs_m))

predict_matrix = matrix(NA,y_res,x_res)

xy.in = !is.na(over(AA,border)) 
isIsland = !is.na(over(AA,island_border))

data.pred <- as.data.frame(raster::extract(input_r, AA))
isOK = complete.cases(data.pred)
isProtected = !is.na(data.pred$protected) & data.pred$protected == 1
data.pred = data.pred[xy.in & isOK & !isProtected,]

pxls_adm_pred = raster::extract(pxls_adm, AA)[xy.in & isOK & !isProtected]

xy.in[(xy.in & !isOK) | (xy.in & isProtected )] = FALSE

A.pred <- inla.spde.make.A(mesh = mesh_b,loc = as.matrix(pred_grid[xy.in,]))

### CREATE SPDE MODEL
pcprec <- list(prior='pcprec', param=c(1, 0.01))
# SPDE model definition
spde <- inla.spde2.matern(mesh = mesh_b,
                          theta.prior.mean = c(0,0),
                          theta.prior.prec = c(0.1,0.1))
s.index = inla.spde.make.index(name = "spatial.field",n.spde = spde$n.spde)


lvst_type = c("catt","pigs","smru","duck","poul")
lvst_names = c("Cattle", "Pigs", "Small ruminants", "Ducks","Poultry")

### Generate random hold out samples per livestock type, 
###  but preserve proportion of zeros
nr_test_obs = round(nrow(data)*.3)
NR_TEST = 5
set.seed(123)
test_index_data = list()
for (curr_lvst in lvst_type) {
  ind_zero = which(data[,curr_lvst] == 0)
  ind_not_zero = which(data[,curr_lvst] != 0)
  
  nr_zero_samples = round(length(ind_zero) / nrow(data) * nr_test_obs)
  nr_not_zero_samples = nr_test_obs  - nr_zero_samples
  
  sample_dat = sapply(1:NR_TEST,
                      function(x) {
                        c(sample(ind_zero,size = nr_zero_samples,replace = FALSE),
                          sample(ind_not_zero,size = nr_not_zero_samples,replace = FALSE))})
  test_index_data[[curr_lvst]] = sample_dat
}

### INLA model
mod_full <- lvst ~ -1 + intercept + tt + ai + pop + el + road2p5  + tsetse +
  crop + grass + forest + urban +
  f(spatial.field, model = spde)
# mod_full <- lvst ~ -1 + intercept + tt + I(tt^2) + ai + pop + el + road2p5  + tsetse +
#   crop + grass + forest + urban +
#   f(spatial.field, model = spde)


# Setup for the est. loop
#curr.family = "zeroinflatednbinomial1"
curr.family = "nbinomial"

out_spec = paste0("y",RUN_YEAR,"_",curr.family)
# out_spec = paste0(curr.family,"_ttSqrLogPop")

dat.rmse = data.frame()
output_list = list()

### Loop through livestock and test cases

for (curr_lvst in lvst_type) {
  cat(curr_lvst,"\n")
  
  ### ORGANIZE DATA
  # Stack data
  stack.est <- inla.stack(
    data = list(lvst = data[,curr_lvst]),
    A = list(A,1,1,1,1,1,1,1,1,1,1,1,1) ,
    effects = list(
      c(s.index, list(intercept = 1)),
      list(tt = data$tt),
      list(ai = data$ai),
      list(pop = data$pop),
      list(el = data$el),
      list(road1 = data$road1),
      list(road2p5 = data$road2p5),
      list(road5 = data$road5),
      list(tsetse = data$tsetse),
      list(crop = data$crop),
      list(grass = data$grass),
      list(forest = data$forest),
      list(urban = data$urban)
    ),
    tag = "full_est")

  stack.pred.response = inla.stack(
    data = list(lvst =NA),
    A = list(A.pred,1,1,1,1,1,1,1,1,1,1,1,1) ,
    effects = list(
      c(s.index, list(intercept = 1)),
      list(tt = data.pred$tt),
      list(ai = data.pred$ai),
      list(pop = data.pred$pop),
      list(el = data.pred$el),
      list(road1 = data.pred$road1),
      list(road2p5 = data.pred$road2p5),
      list(road5 = data.pred$road5),
      list(tsetse = data.pred$tsetse),
      list(crop = data.pred$crop),
      list(grass = data.pred$grass),
      list(forest = data.pred$forest),
      list(urban = data.pred$urban)
    ),
    tag = "pred")

  join_full_stack = inla.stack(stack.est,stack.pred.response)

  # 1. Do overall prediction with full sample
  full.output <- inla(mod_full,
                       data = inla.stack.data(join_full_stack),
                       family = curr.family,
                       control.predictor = list(A=inla.stack.A(join_full_stack), link = 1,compute = TRUE),
                       #control.fixed = list(mean = 0,prec = .0001),
                       control.compute = list(dic = TRUE,cpo = TRUE,config = TRUE,
                                              return.marginals.predictor=TRUE),
                       verbose = FALSE)
  # In-sample fit
  index.full_est = inla.stack.index(join_full_stack,tag="full_est")$data
  index.pred = inla.stack.index(join_full_stack,tag="pred")$data
  post.mean.full_est.fitted = full.output$summary.fitted.values[index.full_est,"mean"]
  obs.full = data[,curr_lvst]
  RMSE.full = sqrt(mean((obs.full - post.mean.full_est.fitted)^2))

  # Plot projected results
  glw <- raster(paste0("./data/2010/maps/glw/",curr_lvst,"_DA.tif"))
  extent(glw) = extent(el)
  names(glw) = "value"
  glw_reg <- raster(paste0("./data/2010/maps/glw/",curr_lvst,"_AW.tif"))
  extent(glw_reg) = extent(el)
  names(glw_reg) = "value"
  max_value = max(getValues(glw),na.rm = TRUE) * 2

  sum_integral = function(x,total,max_val) {
    check_val = sum(x) / total
    while (abs(1 - check_val) > 10^-6) {
      x = x / check_val
      if (any(x > max_val)) {
        x[x > max_val] = max_val
      }
      if (any(x < 0)) {
        x[x< 0] = 0
      }
      check_val = sum(x) / total
    }
    return(x)
  }
  sum_integral2 = function(x,totals,mapping,max_val,
                           return_nr_cut = FALSE) {
    nr_cut = 0
    for (j in 1:length(totals)) {
      curr_x = x[mapping == j]
      curr_total = totals[j]
      if (curr_total == 0) {
        x[mapping == j] = 0
      } else {
        check_val = sum(curr_x) / curr_total
        while (abs(1 - check_val) > 10^-6) {
          curr_x = curr_x / check_val
          if (any(curr_x > max_val)) {
            nr_cut = nr_cut + sum(curr_x > max_val)
            curr_x[curr_x > max_val] = max_val
          }
          if (any(curr_x < 0)) {
            curr_x[curr_x< 0] = 0
          }
          check_val = sum(curr_x) / curr_total
        }
        x[mapping == j] = curr_x
      }
    }
    if (return_nr_cut) {
      return(list(nr_cut = nr_cut,x = x))
    } else {
      return(x)
    }
  }

  post.mean.pred.fitted = full.output$summary.fitted.values[index.pred,"mean"]
  
  post.res_list =
    sum_integral2(
      post.mean.pred.fitted,
      totals = unlist(lvst_dis[match(adm$ADM2_NAME,lvst_dis$district),curr_lvst]),
      mapping = pxls_adm_pred,max_val = max_value,
      return_nr_cut = TRUE)
  cat(paste0("Max. value: ",round(max_value,2),
             " Posterior max:",round(max(post.res_list$x),2),
             " Number cut: ",post.res_list$nr_cut,", ",
             round(post.res_list$nr_cut / length(post.mean.pred.fitted),2),"% \n"))
   

  post.0.025_cut = sum_integral2(
    full.output$summary.fitted.values[index.pred,"0.025quant"],
    total = unlist(lvst_dis[match(adm$ADM2_NAME,lvst_dis$district),curr_lvst]),
                            mapping = pxls_adm_pred,max_val = max_value,
    return_nr_cut = TRUE)$nr_cut
  cat(paste0("  2.5th quantile - Number cut: ",post.0.025_cut,", ",
             round(post.0.025_cut / length(post.mean.pred.fitted),2),"% \n"))
  post.0.975_cut = sum_integral2(
    full.output$summary.fitted.values[index.pred,"0.975quant"],
    total = unlist(lvst_dis[match(adm$ADM2_NAME,lvst_dis$district),curr_lvst]),
    mapping = pxls_adm_pred,max_val = max_value,
    return_nr_cut = TRUE)$nr_cut
  cat(paste0("  97.5th quantile - Number cut: ",post.0.975_cut,", ",
             round(post.0.975_cut / length(post.mean.pred.fitted),2),"% \n"))
  
  cut_nr_values = c(post.res_list$nr_cut,post.0.025_cut,post.0.975_cut)
  names(cut_nr_values) = c("mean","q2.5","q97.5")
  
  # 2. Loop over train / test cases
  curr_test_nr = 1
  curr_rmses = c()
  curr_rmses_intercept = c()
  curr_tests = list()

  for (curr_test_nr in 1:NR_TEST) {
    cat(curr_test_nr," of ",NR_TEST,"\n")
    curr_test_ind = test_index_data[[curr_lvst]][,curr_test_nr]
    data.train = data[-curr_test_ind,]
    data.test = data[curr_test_ind,]
    
    A.train = inla.spde.make.A(mesh = mesh_b, loc = coords[-curr_test_ind,])
    A.test = inla.spde.make.A(mesh = mesh_b, loc = coords[curr_test_ind,])
    
    # Training data
    stack.train <- inla.stack(
      data = list(lvst = data.train[,curr_lvst]),
      A = list(A.train,1,1,1,1) ,
      effects = list(
        c(s.index, list(intercept = 1)),
        list(tt = data.train$tt),
        list(ai = data.train$ai),
        list(pop = data.train$pop),
        list(el = data.train$el)
      ),
      tag = "train")
    
    stack.train_intercept <- inla.stack(
      data = list(lvst = data.train[,curr_lvst]),
      A = list(A.train) ,
      effects = list(
        c(s.index, list(intercept = 1))
      ),
      tag = "train")
    
    # Validation data
    stack.test  = inla.stack(
      data = list(lvst =NA),
      A = list(A.test,1,1,1,1) ,
      effects = list(
        c(s.index, list(intercept = 1)),
        list(tt = data.test$tt),
        list(ai = data.test$ai),
        list(pop = data.test$pop),
        list(el = data.test$el)
      ),
      tag = "test")
    stack.test_intercept  = inla.stack(
      data = list(lvst =NA),
      A = list(A.test) ,
      effects = list(
        c(s.index, list(intercept = 1))
      ),
      tag = "test")
    
    join_test_stack = inla.stack(stack.train,stack.pred.response,stack.test)
    join_test_stack_intercept = inla.stack(stack.train_intercept,stack.pred.response,stack.test_intercept)
    
    test.output <- inla(mod_full,
                        data = inla.stack.data(join_test_stack),
                        family = curr.family,
                        control.predictor = list(A=inla.stack.A(join_test_stack), link = 1,compute = TRUE),
                        #control.fixed = list(mean = 0,prec = 1),
                        control.compute = list(dic = TRUE,cpo = TRUE,config = TRUE,
                                               return.marginals.predictor=TRUE),
                        verbose = FALSE)
    
    test.output_intercept <- inla(mod_full, 
                        data = inla.stack.data(join_test_stack_intercept),
                        family = curr.family, 
                        control.predictor = list(A=inla.stack.A(join_test_stack_intercept), link = 1,compute = TRUE),
                        #control.fixed = list(mean = 0,prec = 1),
                        control.compute = list(dic = TRUE,cpo = TRUE,config = TRUE,
                                               return.marginals.predictor=TRUE),
                        verbose = FALSE)
    
    
    index.pred.test = inla.stack.index(join_test_stack,tag="test")$data
    post.mean.test.fitted = test.output$summary.fitted.values[index.pred.test,"mean"]
    obs.test = data[curr_test_ind,curr_lvst]
    min_val = min(mean(obs.test) - sd(obs.test),0)
    max_val = mean(obs.test) + sd(obs.test)
    post.mean.test.fitted[post.mean.test.fitted<min_val] = min_val
    post.mean.test.fitted[post.mean.test.fitted>max_val] = max_val
    RMSE.test = sqrt(mean((obs.test - post.mean.test.fitted)^2))
    
    index.pred.test_intercept = inla.stack.index(join_test_stack_intercept,tag="test")$data
    post.mean.test.fitted_intercept = test.output_intercept$summary.fitted.values[index.pred.test_intercept,"mean"]
    obs.test = data[curr_test_ind,curr_lvst]
    RMSE.test_intercept = sqrt(mean((obs.test - post.mean.test.fitted_intercept)^2))
    
    curr_rmses = c(curr_rmses,RMSE.test)
    curr_rmses_intercept = c(curr_rmses_intercept,RMSE.test_intercept)
  }
  
  output_list[[curr_lvst]] = list(
                                  rmses = curr_rmses,
                                  curr_rmses_intercept = curr_rmses_intercept,
                                  nr_cut = cut_nr_values)
  
}

rmse_summary = data.frame()
for (curr_lvst in lvst_type) {
  rmse_summary = rmse_summary %>%
    bind_rows(
      data.frame(lvst = curr_lvst,
                 outSampleRMSE_mean = mean(output_list[[curr_lvst]]$rmses),
                 outSampleRMSE_ratio = mean(output_list[[curr_lvst]]$rmses / output_list[[curr_lvst]]$curr_rmses_intercept),
                 outSampleRMSE_sd = sd(output_list[[curr_lvst]]$rmses))
    )
}
rmse_summary_intercept = data.frame()
for (curr_lvst in lvst_type) {
  rmse_summary_intercept = rmse_summary_intercept %>%
    bind_rows(
      data.frame(lvst = curr_lvst,
                 outSampleRMSE_mean = mean(output_list[[curr_lvst]]$curr_rmses_intercept),
                 outSampleRMSE_sd = sd(output_list[[curr_lvst]]$curr_rmses_intercept))
    )
}
require(openxlsx)
fname = paste0("Review_tests_",Sys.Date(),"_",NR_TEST)
wb <- createWorkbook()
addWorksheet(wb,"RMSEs")
writeData(wb, sheet = 1,rmse_summary)
saveWorkbook(wb, paste0(fname,".xlsx"), overwrite = TRUE)


filename = paste0("./output_",out_spec,"_",Sys.Date())

save(data,coords,output_list,lvst_type,lvst_names,
     file = paste0(filename,"_reviewSummary.RData"))





