
########################################################################################################
########  Geostatistical downscaling (prediction and simulation) with SoilGrids data  ##################
########################################################################################################


##### 1. Load and process SoilGrids maps (more info at https://www.isric.org/explore/soilgrids) #####
##### Maps can be downloaded from the site https://soilgrids.org/ ###################################
##### or at the link: https://files.isric.org/soilgrids/data/recent/ ################################

### In the working directory there must be:
### - A directory 'SoilGrids' containing SoilGrids' maps of the variable in exam covering the area of interest
### - A raster (.tif or .tiff) file of the Digital Elevation Model for the area in exam.

# Install all packages used (if required).
packages = c('raster','gstat','compositions', 'dissever', 'fields', 'soiltexture', 'viridis','psych','latex2exp')
for (pck in packages) {
  if ( ! (pck %in% installed.packages()[,'Package']) ) install.packages(pck) 
}

library(raster) # To handle raster data.

# Load the Digital Elevation Model (DEM). The resolution, extent and coordinate system will be used as reference.
dem=raster('DEM.tif')

# Plot
plot(dem)

# Coordinates Reference System (CRS) of the dem used as reference.
crs_ref=crs(dem) 
# Reference extent.
extent_ref=extent(dem)  

# Load particle-size fractions (psf) into a stack (raster with multiple data associated to each pixel) 
psf = stack( c('SoilGrids/CLYPPT_M_sl1_250m.tiff',    # clay
               'SoilGrids/SLTPPT_M_sl1_250m.tiff',    # silt
               'SoilGrids/SNDPPT_M_sl1_250m.tiff' ) ) # sand

### The .tiff files are searched in the directory 'SoilGrids'.
### The '1' in 'sl1' in the names of the raster files indicates that the values considered afre those of the topsoil
### (0 cm from the top). Number 2 correspond to 15 cm, and so on up to number 7 following the standard depths 
### (0-15-30-60-100-200cm).

#To load psf at different depths
# psf5cm = stack( c('SoilGrids/CLYPPT_M_sl2_250m.tiff',  # clay
#                'SoilGrids/SLTPPT_M_sl2_250m.tiff',     # silt
#                'SoilGrids/SNDPPT_M_sl2_250m.tiff' ) )  # sand

# Change the names of the layers to CLAY, SILT and SAND.
names(psf) = c('CLAY','SILT','SAND') 

# Change the coordinate system of the stack of SoilGrids maps to match the reference one.
psf=projectRaster(psf, crs = crs_ref) 

# Crop the maps to match the extent of the DEM (SoilGrids maps are originally broader)
psf=crop(psf, extent_ref) #, filename='psf.tiff')  ( The processed data can be saved )

# Create a Data Frame (df) containing the values of each coarse 'pixel'.
psf_df = as.data.frame(psf) 

# Each row must sum to 100 (%), but there is a residual (res) error due to rounding
res=rep(100,nrow(psf_df))-apply(psf_df,1,sum,na.rm=TRUE) 

# Apply (uniform) correction to the data frame and the rasters. 
for(i in 1:ncol(psf_df)) {
  psf_df[,i] = psf_df[,i] + res/ncol(psf_df)
  values(psf[[i]]) = psf_df[,i] 
}

# Load Absolute Depth to Bedrock (ADB), Censored Depth to Bedrock (CDB, up to 200 cm), and 
# Probability of occurrence of the R Horizon (i.e., the bedrock) within the threshold of 2 m.

# NOT USED
# ADB=raster('SoilGrids/BDTICM_M_250m.tiff')
# ADB=projectRaster(ADB, crs = crs_ref)
# ADB=crop(ADB, extent_ref)#, filename='ADB.tiff')
# # Convert to meters.
# ADB = ADB/100 

CDB=raster('SoilGrids/BDRICM_M_250m.tiff')
CDB=projectRaster(CDB, crs = crs_ref)
CDB=crop(CDB, extent_ref) # , filename='CDB.tiff')

PRH=raster('SoilGrids/BDRLOG_M_250m.tiff')
PRH=projectRaster(PRH, crs = crs_ref)
PRH=crop(PRH, extent_ref) # , filename='PRH.tiff')

##### Plots #####

# viridis colors (for pretty plots)
library(viridis)

### Functions ###

# Plot the histogram (with optional fitted density) of raster values.
plot_raster_hist = function(rstr, length = 40, colors = terrain.colors(40), title = "", 
                            density_col = "blue", xlim = NA, fit.density = TRUE) {
  
  ### Arguments:
  # rstr = raster map
  # length = number of intervals 
  # density_col = color of the denisty line
  # fit.density = boolean, if true fits a density and plots it
  
  vals = values(rstr)
  min_val <-min(vals, na.rm = TRUE)
  max_val <-max(vals, na.rm = TRUE)
  if(is.na(xlim)){xlim = c(min_val,max_val)}
  
  hist(vals, breaks=seq(min_val,max_val,length = length), 
       col=colors, xlab="",
       main=title,
       xlim=xlim, 
       freq=FALSE) 
  if(fit.density) { lines(density(vals, na.rm=TRUE), col = density_col, lwd = 2) }
}

# Plot a map and the histogram of its values on the side. 
quartz(width = 10, height = 5) ### For windows it is not 'quartz()' but 'x11()' !!! (consider using 'replace all').
par(mfrow=c(1,2),mai = c(1,1,1,1)) # For multiple plots in the same figure. mai = margins
plot(psf$CLAY, main = "Map of clay %")
plot_raster_hist(psf$CLAY, title = "Histogram of clay %") # Function defined above

# Plot two rasters side by side.
quartz(width = 10, height = 5)
par(mfrow=c(1,2),mai = c(1,1,1,1))
plot(CDB, main = "Map of CDB (cm)", col = viridis(20))
plot(PRH, main = "Map of PRH (%)", col = viridis(20))


##### 2. Analysis of soil texture  ###############################################################
##### https://cran.r-project.org/web/packages/soiltexture/vignettes/soiltexture_vignette.pdf #####

library(soiltexture)

# Plot psf on the soil texture triangle with the different soil types of the USDA classification.
geo = TT.geo.get()
quartz(width = 6, height = 6)
TT.plot(class.sys = 'USDA.TT')
TT.plot(class.sys = 'USDA.TT', tri.data = psf_df[-which(is.na(psf_df$CLAY)),], geo = geo, #grid.show = FALSE, 
        col = "blue", cex=.5, lwd=.5) #, add = FALSE)

# Get classes (one-hot encoding) for each observation.
soil_type_one_hot= TT.points.in.classes(class.sys = 'USDA.TT', tri.data = psf_df[-which(is.na(psf_df$CLAY)),])
classes = colnames(soil_type_one_hot)
# Get numeric encoding.
soil_type = apply(soil_type_one_hot,1,FUN = function(x) as.vector(x) %*% 1:ncol(soil_type_one_hot) )

# To get the name of the class of the pixel in position i
#classes[soil_type[i]]

# Map of soil types (numerical code)
soil_type_map = psf$CLAY
values(soil_type_map)[-which(is.na(psf_df$CLAY))] = soil_type
plot(soil_type_map)

# Get modal soil type over the area.
Mode = function(x) {
  ux = unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}
classes[Mode(soil_type)]


##### 3. Downscaling and simulation of psf through Isometric Log-Ratio Area-to-Point Regression Kriging (ILR_ATPRK) ###

### Functions ###################

variogram_deconvolution = function(coarse_raster, vg_type = "Sph", nblocks = 4, maxiter = 100, cutoff = "default",
                                   tol1 = 1e-2, tol2 = 1e-6) { 
  
  ##### Variogram deconvolution procedure (Goovaerts, 2008) #########
  ##### This is a special version for regular grids (raster data) ###
  
  require(gstat)
  require(fields)
  
  ########### Arguments: ###################################################################################
  
  # coarse_raster = coarse resolution raster. 
  # vg_type = type of variogram to be fitted.
  # nblocks = the regularized variogram is computed using R function vgmArea, but only at short distances,
  #           nblocks is the number of adjacent blocks that are considered near enough to require vgmArea.
  # maxiter = maximum number of iterations.
  # cutoff = cutoff. If = "default" it is set equal to half the raster extent.
  # tol1 = tolerance for Di/D0, if the value is lower the iterations stop.
  # tol2 = tolerance for abs(D_i-D_opt)/D_opt.
  
  ##########################################################################################################
  
  # Borders and extent of the raster map.
  xmin = coarse_raster@extent@xmin
  xmax = coarse_raster@extent@xmax
  ymin = coarse_raster@extent@ymin
  ymax = coarse_raster@extent@ymax
  coarse_res = res(coarse_raster)
  x_extent = xmax - xmin
  y_extent = ymax - ymin
  
  if (cutoff == "default") {
    cutoff = min(c(x_extent,y_extent))/2 # cutoff equal to half the raster extent by default.
  }
  
  # Change the name to a generic "z".
  names(coarse_raster) = 'z' 
  
  # Convert the raster to a list of square polygons corresponding to the coarse pixels (=blocks).
  poly_coarse = rasterToPolygons(coarse_raster, fun=NULL, n=16, na.rm=TRUE, digits=12, dissolve=FALSE)
  
  # 1. Compute empirical variogram on areal data and fit a model.
  gamma_v_hat = variogram(z ~ 1, data = poly_coarse, cutoff=cutoff, width = min(coarse_res)/2)
  gamma_v_exp = fit.variogram(gamma_v_hat, vgm(mean(tail(gamma_v_hat$gamma),10), vg_type, 
                                               cutoff/2, min(gamma_v_hat$gamma)))
  
  # 2. Initial variogram.
  gamma_0 = gamma_v_exp
  
  # 3. Variogram regularization using 'vgmArea' for low distances, and Journel approximation for high distances.
  # This procedure is specific for regular grids and allows to significantly speed up the algorithm,
  # avoiding the repeated computation of block covariances, or the computation of block covariances at great
  # distances, which can be approximated using the Journel formula (Journel, 1978).
  # help(vgmArea)
  new_xmax = xmin + nblocks*coarse_res[1]
  new_ymin = ymax - nblocks*coarse_res[2]
  inrange = crop(coarse_raster, extent(c(xmin,new_xmax,new_ymin,ymax)))
  poly_ref = polygons(rasterToPolygons(inrange, fun=NULL, n=16, na.rm=TRUE, digits=12, dissolve=FALSE))
  # Now poly_ref contains a square of blocks, given the regularity of the problem we can drop the upper part
  # of the square since covariance only depends on the distance of the blocks.
  polygon_indexes = upper.tri(matrix(1:(length(poly_ref)),nrow = nblocks, ncol = nblocks), diag = TRUE)
  polygon_indexes = as.vector(polygon_indexes)
  poly_ref = poly_ref[polygon_indexes]
  # Compute gamma^(v,v_h) for small lags. 
  gamma_A_0 = vgmArea(x = poly_ref[1],y = poly_ref, vgm = gamma_0, covariance = FALSE)
  gamma_vv_0 = gamma_A_0[1] # gamma(v,v), unique since the grid is regular.
  coords_ref = as.matrix(coordinates(poly_ref))
  short_dist =  as.vector( rdist(t(coords_ref[1,]),coords_ref[2:nrow(coords_ref),]) ) # Distances of neighboring blocks.
  drop_duplicates = !duplicated(short_dist)
  # Drop distance duplicates (the block covariance only depend on distance).
  short_dist = short_dist[drop_duplicates]
  ordered = order(short_dist)
  short_dist = short_dist[ordered]
  gamma_v_0 = as.vector( gamma_A_0[2:length(gamma_A_0)] )
  gamma_v_0 = gamma_v_0[drop_duplicates]
  gamma_v_0 = gamma_v_0[ordered] # Values of regularized variogram for small distances
  # Add lags (greater distances).
  dist_tail = seq(short_dist[length(short_dist)]+max(coarse_res),cutoff,max(coarse_res))
  ref_dist = c(short_dist, dist_tail)
  gamma_v_0_tail = variogramLine(gamma_0,dist_vector=dist_tail)$gamma
  gamma_v_0 = c(gamma_v_0,gamma_v_0_tail)
  gamma_v_0 = gamma_v_0 - gamma_vv_0 # Regularize
  
  # 4. Quantify deviation.
  D_0 = ( 1/length(ref_dist) ) * 
    sum( abs(gamma_v_0-variogramLine(gamma_v_exp,dist_vector=ref_dist)$gamma)/
           variogramLine(gamma_v_exp,dist_vector=ref_dist)$gamma )
  
  # 5. Define initial optimal variograms.
  gamma_opt = gamma_0
  gamma_v_opt = gamma_v_0
  D_opt = D_0
  rescaling_flag = 0
  
  # Start loop.
  for (i in 1:maxiter){
    print(paste0("iter: ", i))
    
    # 6. Compute experimental values for the new point support semivariogram through a rescaling of the optimal 
    #    point support model. 
    if (!rescaling_flag){
      w = 1 + ( 1/(gamma_v_exp$psill[2]* (i^(1/2)) ) )*
        (variogramLine(gamma_v_exp, dist_vector=ref_dist)$gamma - gamma_v_opt)
    }
    else {
      w = 1+(w-1)/2
      rescaling_flag = 0
    }
    
    # Empirical values to which a variogram model must be fitted.
    gamma_hat_i_values = variogramLine(gamma_opt,dist_vector=ref_dist)$gamma * w
    gamma_hat_i = gamma_v_hat[1:length(ref_dist),]
    gamma_hat_i$np=rep(1,length(ref_dist))
    gamma_hat_i$dist = ref_dist 
    gamma_hat_i$gamma = gamma_hat_i_values
    gamma_hat_i$dir.hor[is.na(gamma_hat_i$dir.hor)] = 0
    gamma_hat_i$dir.ver[is.na(gamma_hat_i$dir.ver)] = 0
    gamma_hat_i$id[is.na(gamma_hat_i$id)] = gamma_hat_i$id[1]
    
    # 7. Fit a new model.
    gamma_i = fit.variogram(object = gamma_hat_i, vgm(mean(tail(gamma_hat_i$gamma),2), vg_type, 
                                                      cutoff, min(gamma_hat_i$gamma)))
    
    # 8. Regularize gamma_i.
    gamma_A_i = vgmArea(x = poly_ref[1],y = poly_ref, vgm = gamma_i, covariance = FALSE)
    gamma_vv_i = gamma_A_i[1]
    gamma_v_i = as.vector( gamma_A_i[2:length(gamma_A_i)] )
    gamma_v_i = gamma_v_i[drop_duplicates]
    gamma_v_i = gamma_v_i[ordered]
    
    gamma_v_i_tail = variogramLine(gamma_i,dist_vector=dist_tail)$gamma
    gamma_v_i = c(gamma_v_i,gamma_v_i_tail)
    gamma_v_i = gamma_v_i - gamma_vv_i
    
    # 9. Compute D_i.
    D_i= ( 1/length(ref_dist) ) * 
      sum( abs(gamma_v_i-variogramLine(gamma_v_exp,dist_vector=ref_dist)$gamma)/
             variogramLine(gamma_v_exp,dist_vector=ref_dist)$gamma )
    
    # 10. Stopping criteria. 
    if ( (D_i/D_0 < tol1) || abs(D_i-D_opt)/D_opt < tol2) break
    if (D_i < D_opt) {
      gamma_opt = gamma_i
      gamma_v_opt = gamma_v_i
      D_opt = D_i
    }
    else {
      rescaling_flag = 1
    }
    
  }
  return(gamma_opt)
}

ATPlm = function(coarse_raster, covariates_stack){
  
  ### Area-to-Point Linear Regression.
  ### Remark: the target coarse resolution raster and the covariate maps must have the same coordinate system
  ### and must overlap.
  ### The function returns an object containing a fine resolution map (matching the resolution of the covariates)
  ### of the fitted target variable and an lm-type object containig the regression results.
  
  df = as.data.frame(coarse_raster[[1]])
  xnames = names(covariates_stack) # Covariate names.
  yname = names(coarse_raster)[1] # Name of target variable to be used in the formula.
  nx = length(xnames)
  # Upscale covariate maps for the regression.
  for (i in 1:nx){
    df[xnames[i]] = values( resample(covariates_stack[[i]],coarse_raster, method = "bilinear") )
  }
  f <- as.formula( paste(yname, "~ .") )
  LM = lm(formula = f, data = df) # Linear regression
  predicted = predict(LM, as.data.frame(covariates_stack))
  map = covariates_stack[[1]]
  values(map) = predicted
  out = list(map = map, LM = LM)
  return(out)
}

ATPK = function(coarse_raster, fine_raster, variogram, npoints = 8, nsim = 0, beta = NA, nmax = 10,
                frac = 1.0, noise_sd = 0.0){
  
  ### Area-to-Point kriging for the downscaling of coarse raster data.
  ### Converts the coarse pixels into spatial polygons and uses the ATPK method of Kyriakydis (2004).
  ### The function also performs Block Sequential Gaussian Simulation (BSGS).
  ### Returns the downscaled map or a stack of the simulations.
  
  require(gstat)
  
  ### Arguments: 
  # coarse_raster = coarse resolution raster to be downscaled.
  # fine_raster = fine resolution raster covering the extent of the coarse raster. It will be used to define 
  #               the target resolution for the downscaling. The values in this raster are not used.
  # frac = fraction of coarse data to use when performing block sequential simulation (increases the variability)
  # npoints = number of oints used to approximate the blocks (coarse pixels). Can be 4,8 or 16.
  # For the other parameters consult the documentation of 'krige' function - help(krige) 
  
  # If simulations are required, cosnider only the specified fraction of areal data.
  if (nsim > 0){
    drop = 1.0-frac
    ndrops = floor(drop*length(coarse_raster))
    values(coarse_raster)[sample(1:length(coarse_raster),ndrops,replace = FALSE)] = NA
  }
  
  BlockMap = rasterToPolygons(coarse_raster, fun=NULL, n=npoints, na.rm=TRUE, digits=12, dissolve=FALSE)
  names(BlockMap) = "Z"
  
  # Add optional noise to the block data
  if(noise_sd){
    print("Adding an uniform noise to the block data")
    noise = rnorm(nrow(BlockMap),mean = 0, sd = noise_sd)
    BlockMap@data = BlockMap@data + noise
  }
  
  d_new <- SpatialPoints(coordinates(fine_raster)[!is.na(values(fine_raster)),], proj4string = crs(BlockMap))
  if(is.na(beta)) {
    downscale = krige(Z~1, BlockMap, newdata = d_new, model = variogram, nmax = nmax, nsim=nsim)
  } else {
    downscale = krige(Z~1, BlockMap, newdata = d_new, model = variogram, nmax = nmax, nsim=nsim, beta = beta)
  }
  if (nsim == 0) {
    downscaled_map = fine_raster
    values(downscaled_map)[!is.na(values(fine_raster))] = downscale$var1.pred
    return(downscaled_map)
  } else {
    sim_data = downscale@data 
    sims = stack(replicate(nsim,fine_raster))
    for (i in 1:nsim) values(sims[[i]])[!is.na(values(fine_raster))] = sim_data[,i]
    return(sims)
  }
}  

#################################

### Compositional data analysis (Aitchison, 1986).

library(compositions) # For compositional data analysis.

### ILR transformation (Egozcue et al., 2003).
help(ilr)

# Apply closure operation to the psf dataset (acomp turns a df into a compositional object).
psf_comp = acomp(psf_df[,c(1,3,2)]) 
### Remark: silt and sand were switched to match the order of the thesis!

# Apply Isometric Log-Ratio to the data (with default basis), and put the results in a data frame.
ilr_df = data.frame( ilr(x = psf_comp) )
colnames(ilr_df) = c('ILR1','ILR2')

# When applying ilr, NAs are turned into 0s, set to NA the original NAs.  
ilr_df[is.na(psf_df$CLAY),] = NA 

# Create ilr rasters.
ILR1 = psf[[1]] # Copy clay map 
ILR2 = ILR1
values(ILR1) = ilr_df$ILR1 # Change values 
values(ILR2) = ilr_df$ILR2

library(latex2exp) # For formulas using LaTeX syntax in plot titles.
quartz(width=10,height = 5)
par(mfrow=c(1,2))
plot(ILR1, main = TeX("$\\bar{ilr}_1$")) # The overline indicates the coarse resolution.
plot(ILR2, main = TeX("$\\bar{ilr}_2$"))

### ATPRK (Area-to-Point Regression Kriging).

# Regression step using Digital Elevation. 
# THIS STEP IS NOT NECESSARY. 
# Upscale covariate map using function 'resample' with method bilinear (upscales by averaging).
coarse_dem = resample(dem, ILR1, method="bilinear")
plot(coarse_dem)

# Add column 'DEM' to the data frame with the ilrs.
ilr_df['DEM'] = values(coarse_dem)

library(psych) # For pretty scatterplots.
quartz(width=6, height=6)
pairs.panels(ilr_df[c('DEM','ILR1','ILR2')], 
             method = "pearson", # correlation method
             hist.col = "#00AFBB",
             density = TRUE,  # show density plots
             ellipses = FALSE, # show correlation ellipses
             labels= c(TeX("$\\bar{ilr}_1$"),TeX("$\\bar{ilr}_2$"),TeX("$\\bar{DEM}$"))
)

# Perform Area-to_Point Linear Regression using dem and dem^2 as covariates.
lm_dem1 = ATPlm(ILR1,stack(c(dem,dem^2)))
lm_dem2 = ATPlm(ILR2,stack(c(dem,dem^2)))

plot(lm_dem2$map)
plot(values(ILR2)[!is.na(values(ILR2))], lm_dem2$LM$fitted.values, col = "blue", cex=.5,
     xlab = "observed", ylab = "fitted", 
     main = TeX("$\\mathbf{E}\\[ILR_2\\] = \\beta_0^{(2)}+ \\beta_1^{(2)} \\cdot DEM + \\beta_2^{(2)} \\cdot DEM^2$") )
abline(a=0,b=1, col = "red", lwd = 2)

### Area to point kriging and simulation.

library(gstat) # Kriging, ATPK, variogram estimation, (Block) Sequential Gaussian Simulation.

# Compute residuals.
RES1 = ILR1 - resample(lm_dem1$map, ILR1, method="bilinear")
RES2 = ILR2 - resample(lm_dem2$map, ILR2, method="bilinear")

# Compute empirical variograms (as a reference to choose the variogram model).
emp_vg1 = variogram(layer ~ 1, data = as(RES1, 'SpatialPointsDataFrame')) 
emp_vg2 = variogram(layer ~ 1, data = as(RES2, 'SpatialPointsDataFrame')) 

# v1 = fit.variogram(emp_vg1, vgm(model = "Sph"))
# plot(emp_vg1, v1, xlab=TeX("$h$"),ylab=TeX("$\\gamma(h)$")) # Spherical
# v2 = fit.variogram(emp_vg2, vgm(model = "Exp"))
# plot(emp_vg2, v2) # Exponential

# Perform variogram deconvolution.
vd1 = variogram_deconvolution(RES1, vg_type = "Sph", nblocks = 4, maxiter = 100) 
vd2 = variogram_deconvolution(RES2, vg_type = "Exp", nblocks = 4, maxiter = 100) 

########################## Plot the variograms. ################################################
quartz(width=8, height=6)
plot(emp_vg1$dist, emp_vg1$gamma, ylim = c(0,0.012), col = 'blue', lwd=2,
     main = TeX("Variogram deconvolution of $\\gamma_1$"), ylab = "semivariance", xlab = "distance")
lines(seq(0,4500,45),variogramLine(v1,dist_vector=seq(0,4500,45))$gamma, type='l', col = 'green', lwd = 2)
lines(seq(0,4500,45),variogramLine(vd1,dist_vector=seq(0,4500,45))$gamma, type='l', col = 'red', lwd = 2)
grid(col = 'cornsilk2')
legend(x = 3000,0.006, legend=c('Empirical variogram','Initial fit', 'Deconvoluted variogram'), 
       col = c('blue','green','red'), lty = c(2,1,1), cex = .8)

quartz(width=8, height=6)
plot(emp_vg2$dist, emp_vg2$gamma, ylim = c(0,0.003), col = 'blue', lwd=2,
     main = TeX("Variogram deconvolution of $\\gamma_2$"), ylab = "semivariance", xlab = "distance")
lines(seq(0,4500,45),variogramLine(v2,dist_vector=seq(0,4500,45))$gamma, type='l', col = 'green', lwd = 2)
lines(seq(0,4500,45),variogramLine(vd2,dist_vector=seq(0,4500,45))$gamma, type='l', col = 'red', lwd = 2)
grid(col = 'cornsilk2')
legend(x = 3000,0.001, legend=c('Empirical variogram','Initial fit', 'Deconvoluted variogram'), 
       col = c('blue','green','red'), lty = c(2,1,1), cex = .8)
################################################################################################

### ATPK of the residuals.

dem2 = aggregate(dem, fact=4) # Reduce the resolution
# Sequential simulation can take several minutes when simulating very high resolution maps, so the resolution
# is reduced by at least a 2 factor (this is done by fuction aggregate).

# Create fine resolution mask with null values over the lake.
mask = resample(psf[[1]], dem2, method='ngb')
values(mask)[!is.na(values(mask))] = 1

# Downscaling.
down1 = ATPK(RES1, mask, variogram = vd1, nmax = 8, beta = 0)
down2 = ATPK(RES2, mask, variogram = vd2, nmax = 8, beta = 0)

# Simulation.
sim1 = ATPK(RES1, mask, variogram = vd1, nmax = 8, beta = 0, nsim = 1, frac = .5, noise_sd = 0.1)
sim2 = ATPK(RES2, mask, variogram = vd2, nmax = 8, beta = 0, nsim = 1, frac = .5, noise_sd = 0.1)

# Add the baseline (regression results).
baseline1 = aggregate(lm_dem1$map, fact = 4)
baseline2 = aggregate(lm_dem2$map, fact = 4)
for (i in 1:length(sim1@layers)) {
  sim1[[i]] = sim1[[i]] + baseline1
  sim2[[i]] = sim2[[i]] + baseline2
}

# Back transform the results and create stack of rasters
# This time particle size fractions are not stacked, but simulations of each are
clay_sim_stack = stack()
sand_sim_stack = stack()
silt_sim_stack = stack()
for (i in 1:length(sim1@layers)) {
  fine_ilr = as.data.frame(stack(c(sim1[[i]],sim2[[i]])))
  psf_sim_df = ilrInv(fine_ilr)
  # colnames(psf_sim_df) = c("CLAY", "SAND", "SILT")
  clay_sim_i = mask
  sand_sim_i = mask
  silt_sim_i = mask
  values(clay_sim_i) = psf_sim_df[,1]
  values(sand_sim_i) = psf_sim_df[,2]
  values(silt_sim_i) = psf_sim_df[,3]
  clay_sim_stack = stack(clay_sim_stack,clay_sim_i)
  sand_sim_stack = stack(sand_sim_stack,sand_sim_i)
  silt_sim_stack = stack(silt_sim_stack,silt_sim_i)
}

# Plot.
quartz()
par(mfrow=c(1,2))
plot(clay_sim_stack[[1]]*mask, main = "Simulated clay %")
plot_raster_hist(clay_sim_stack[[1]], title = "Histogram of simulated clay %")


# Save the simulations.
for (i in 1:length(sim1@layers)){
  aux_clay = clay_sim_stack[[i]]
  aux_sand = sand_sim_stack[[i]]
  aux_silt = silt_sim_stack[[i]]
  save( aux_clay, file=paste0('clay_sim_',i,'.tiff') )
  save( aux_sand, file=paste0('sand_sim_',i,'.tiff') )
  save( aux_silt, file=paste0('silt_sim_',i,'.tiff') )
}

#########################################################################################################

##### 4. Downscaling and simulation of soil thickness with dissever #####

# Recall: 
# CDB = Censored depth to bedrock up to 200 cm.
# PRH = probability of occurrence of the R-horizon (bedrock) within the first 2 meters.

plot(CDB)
# Set to NA the censored values (previously set = 200).
values(CDB)[which(values(CDB)>199)] = NA
NAs = which(is.na(values(CDB)))

# Absolute slope in radians (discrete gradient angle).
slope = terrain(dem, opt='slope', unit='radians', neighbors=8) # , filename='slope.tif')

help(terrain)
plot(slope, col = viridis(40))

# Upscale DEM and Slope (JUST for graphic comparison!)
coarse_slope = resample(slope, PRH, method = "bilinear")
coarse_dem = resample(dem, PRH, method = "bilinear")
soil_thickness_df = data.frame(DEM = values(coarse_dem), Slope = values(coarse_slope),
                               CDB = values(CDB), PRH = values(PRH) )
pairs.panels(soil_thickness_df, 
             method = "pearson", # correlation method
             hist.col = "#00AFBB",
             density = TRUE,  # show density plots
             ellipses = FALSE # show correlation ellipses
             )


### Dissever algorithm (Malone and McBratney, 2012).

library(dissever)
help(dissever)

# Downscale PRH using dissever with cubist method.
dPRH = dissever(coarse = PRH, fine = stack(c(dem,slope)), method= "cubist", p = .1, min_iter=1, max_iter = 4)
plot(dPRH$map, col = viridis(100))
plot(dem)

#dPRH = dissever(coarse = PRH, fine = slope, method= "lm", p = 0.1, min_iter=1, max_iter = 4)
PRHo = PRH/100
values(PRHo) = log(values(PRHo)/(1-values(PRHo)))
dPRHo = ATPlm(PRHo,stack(c(dem,slope)))
values(dPRHo$map) = 1/(1+exp(-values(dPRHo$map)))
plot(dPRHo$map, col = viridis(100))

values(dPRHo$map)[values(dPRHo$map)<0.37]=NA

quartz()
par(mfrow=c(1,2))
plot(PRH, col = viridis(20))
plot(dPRH$map, col=viridis(40))


# Analysis of the performance.
quartz(width=10, height = 5)
par(mfrow=c(1,2))
plot(dPRH)
plot(dPRH, type='perf')

preds = extractPrediction( list(dPRH$fit) )
plotObsVsPred(preds[sample(1:nrow(preds),3000,replace = FALSE),])

# Set probability threshold t (if PRH > t, then soil thickness < 2m).
t = 35 #%
# Map of thick classes, if thick_class(x) = 1, at location x the soil is considered thinner than 2m.
thick_class = dPRH$map
values(thick_class)[values(thick_class) <= t] = 0
values(thick_class)[values(thick_class) > t] = 1
plot(thick_class)

### Downscaling and simulation of CDB on the locations where thick_class=1 (PRH>t).

# Downscale CDB using ATPlm.
CDBlm = ATPlm(CDB,stack(c(dem,slope)))
plot(CDBlm$map)

CDB_down = CDBlm$map
values(CDBlm$map) = (values(CDB))

cdb_lm_map = resample(CDBlm$map, thick_class, method = 'ngb') 
# cdb_lm_map is the map of the predicted censored depth to bedrock at fine scale.
# Since the extent of the downscaled map does not match exactly the extent of the map thick_class,
# a correction is applied to match the extent and resolution of the two maps. The method used is nearest neighbours.
# Lower case for CDB to indicate that it is at a fine scale.
# Censored depth to bedrock at fine scale is indicated 'cdb'.
# Only locations with PRH>t will be considered, hence the map will be indicated with a 't' (for 'thin') at the end.
# The model assumption for cdb is that the values at each location are equal to an average value given by the linear 
# model, plus an areal residual, plus an additional white noise with zero mean. 

cdbt = resample(CDB, thick_class, method = 'ngb') # Each location has a value equal to the value of the block it belongs to. 
values(cdbt)[values(thick_class) == 0] = NA # Drop block for wich PRH < t at each location.
CDBt = resample(cdbt, CDB, method = 'bilinear') 
values(cdb_lm_map)[values(thick_class) == 0] = NA 

# Compute areal residuals.
CDBres = CDBt - resample(cdb_lm_map, CDBt, method='bilinear')
plot(CDBres)
cdb_res = resample(CDBres, cdbt, method = 'ngb')
# Finally compute point residuals (epsilon)
cdb_epsilon = cdbt-cdb_res-cdb_lm_map

quartz(width = 10, height = 4)
par(mfrow=c(1,3))
plot(cdb_lm_map, col = viridis(40))
plot(cdb_res, col = viridis(40))
plot(cdb_epsilon, col = viridis(40))

plot_raster_hist(cdb_epsilon)

# Simulation of the point residuals (white noise)
SD = sd(values(cdb_epsilon), na.rm = TRUE)
cdb_sim = cdb_lm_map + cdb_res
notNAs = which(!is.na(values(cdb_sim)))
L = length(notNAs)
epsilon = rnorm(L, mean = 0, sd = SD)
values(cdb_sim)[notNAs] = values(cdb_sim)[notNAs] + epsilon
values(cdb_sim)[values(cdb_sim) >= 200] = NA
plot(cdb_sim)
save(cdb_sim, 'Censored Depth to Bedrock simulation.tif')
