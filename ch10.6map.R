#	ch10.6 maps in R

library(sp)  # classes for spatial data
library(raster)  # grids, rasters
library(rasterVis)  # raster visualisation
library(maptools)
library(rgeos)
library(dismo)
# and their dependencies

mymap<-gmap("United States")
plot(mymap)
# capture area in Harms2009
drawExtent()
class       : Extent 
xmin        : -12945681 
xmax        : -9625552 
ymin        : 2907865 
ymax        : 6250534 
e<-drawExtent()
mymap<-gmap(e, type="satellite")
plot(mymap, interpolate=TRUE)


# try USGS data
# data is available from https://gapanalysis.usgs.gov/gaplandcover/data/download/
# put unzipped files in directory called compressed 

setwd( "compressed/gaplndcov_il")
r<-raster('gaplndcov_il/')
r
#class       : RasterLayer 
#dimensions  : 15718, 9143, 143709674  (nrow, ncol, ncell)
#resolution  : 30, 30  (x, y)
#extent      : 687436.9, 961726.9, 1666924, 2138464  (xmin, xmax, ymin, ymax)
#coord. ref. : +proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs 
#data source : /Users/jclarke/Documents/predictive/Book/ch10/VegOut_Jennifer_Clarke/USGS_land_cover/gaplndcov_il/gaplndcov_il 
#names       : gaplndcov_il
#values      : 38, 584  (min, max)
#attributes  :
# Var.1  ID  COUNT VALUE_1  COUNT_1     RED   GREEN    BLUE CL
# from:  38 308785      38 88928526 0.65882 0.74902 0.17647 01
# to  : 584 443412     584 24785929 0.78824 0.30196 0.25882 08
#                   NVC_CLASS   SC         NVC_SUBCL     FRM ...
#           Forest & Woodland 01.C  Temperate Forest 01.C.01 ...
# Developed & Other Human Use 08.A Developed & Urban 08.A.01 ...
plot(r)

setwd("compressed")
for(i in 2:16){
	setwd("compressed")
	dirs<-dir()[i]
	setwd(paste(getwd(),"/",dir()[i],sep=""))
	assign(paste('r',(i-1),sep=""),raster(paste(dirs,"/",sep="")))
}
allr<-merge(r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,r11,r12,r13,r14,r15)
allr
#class       : RasterLayer 
#dimensions  : 91118, 75055, 6838861490  (nrow, ncol, ncell)
#resolution  : 30, 30  (x, y)
#extent      : -1497973, 753676.9, 310894.4, 3044434  (xmin, xmax, ymin, ymax)
#coord. ref. : +proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs 

dim(attr(r, "data")@attributes[[1]])
[1] 54 20
 unique(colTab[as.character(valTab)])
 [1] "#A8BF2D" "#26AD00" "#6E7300" "#368200" "#266E00" "#267357" "#005782" "#666196"
 [9] "#349C79" "#C99964" "#F6A35A" "#C99675" "#6FD6E3" "#8CAFE3" "#C2C9BF" "#8C8F91"
[17] "#FEFEC1" "#A1459C" "#872E26" "#002DC2" "#E01242" "#C94D42"
unique(attr(r, "data")@attributes[[1]][,17])
 [1] Southeastern North American Ruderal Forest & Plantation 
 [2] Central Oak-Hardwood & Pine Forest                      
 [3] Eastern North American Ruderal Forest & Plantation      
 [4] Northern Mesic Hardwood & Conifer Forest                
 [5] Central Mesophytic Hardwood Forest                      
 [6] Northern & Eastern Pine - Oak Forest, Woodland & Barrens
 [7] Northern & Central Floodplain Forest & Scrub            
 [8] Northern & Central Swamp Forest                         
 [9]                                                         
[10] Southern Floodplain Hardwood Forest                     
[11] Great Plains Tallgrass Prairie & Shrubland              
[12] Northern & Central Alvar & Glade                        
[13] Eastern North American Coastal Grassland & Shrubland    
[14] Eastern North American Wet Meadow & Marsh               
[15] Great Plains Wet Meadow, Wet Prairie & Marsh            
[16] Eastern North American Cliff & Rock Vegetation          
[17] Barren                                                  
[18] Herbaceous Agricultural Vegetation                      
[19] Introduced & Semi Natural Vegetation                    
[20] Recently Disturbed or Modified                          
[21] Open Water                                              
[22] Quarries, Mines, Gravel Pits and Oil Wells              
[23] Developed & Urban                                       

# add station locations
# get shape file from Tsegaye so have correct .prj
setwd("/Users/jclarke/Documents/predictive/Book/ch10/Jennifer/")
opt<-readOGR("vegOut_stns.shp", "vegOut_stns")
OGR data source with driver: ESRI Shapefile 
Source: "vegOut_stns.shp", layer: "vegOut_stns"
with 1536 features
It has 9 fields
opt@proj4string
CRS arguments:
 +proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +a=6370997 +b=6370997
+units=m +no_defs 
plot(opt)

CRS.new <- CRS("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
opt.proj<-spTransform(opt,CRS.new) 

par(xpd = FALSE,mar=c(0.5, 2.1, 2.1, 9.1))
plot(opt.proj, pch=16, col="red")	# col needs to be replaced with border
plot(allr, col=attr(r,"legend")@colortable, legend=FALSE, xaxt="n",yaxt="n",add=T)
par(xpd = TRUE)
legend( par()$usr[2], 2000000, legend=unique(attr(r,"data")@attributes[[1]][,17])[c(1:8,10:23)], fill=col, cex=0.35)
plot(opt.proj, pch=16, col="red", add=T)
# save as region3.pdf

vegout<-read.csv('vegout1345test.csv',header=T)
d<-data.frame(lon=vegout$Long, lat=vegout$Lat)
d.uniq<-unique(d)
coordinates(d.uniq) <- c("lon", "lat")
proj4string(d.uniq) <- CRS("+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +a=6370997 +b=6370997+units=m +no_defs")
#proj4string(d.uniq) <- CRS("+init=epsg:4326")	# WGS84
#proj4string(d.uniq) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
d.proj <- spTransform(d.uniq, CRS.new)
points(d.proj, pch=16, col="red")
# with this added save as region3_wd.pdf

# updated map with all states ....

for(i in 2:23){
	setwd("compressed")
	dirs<-dir()[i]
	setwd(paste(getwd(),"/",dir()[i],sep=""))
	assign(paste('r',(i-1),sep=""),raster(paste(dirs,"/",sep="")))
}
allr2<-merge(r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,r11,r12,r13,r14,r15,r16,r17,r18,r19,r20,r21,r22)
allr2
class       : RasterLayer 
dimensions  : 95546, 103827, 9920254542  (nrow, ncol, ncell)
resolution  : 30, 30  (x, y)
extent      : -2361133, 753676.9, 310894.4, 3177274  (xmin, xmax, ymin, ymax)
coord. ref. : +proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs 
data source : /private/var/folders/1h/3_101pxs5zgd6v5prb9fk54w0000gq/T/R_raster_jclarke/r_tmp_2016-01-07_104635_12898_20353.grd 
names       : layer 
values      : 0, 584  (min, max)

par(xpd = FALSE,mar=c(0.5, 2.1, 2.1, 9.1))
plot(opt.proj, pch=16)
plot(allr2, col=attr(r,"legend")@colortable, legend=FALSE, xaxt="n",yaxt="n",add=T)
par(xpd = TRUE)
legend( par()$usr[2], 2000000, legend=unique(attr(r,"data")@attributes[[1]][,17])[c(1:8,10:23)], fill=col, cex=0.35)
plot(opt.proj, pch=16, border="red", add=T, lwd=2.0)
# save as region4.pdf
points(d.proj, pch=16, col="cyan", cex=0.3)
# save as region5.pdf

## remake in bw for book
library(RColorBrewer)
myColours<-colorRampPalette(brewer.pal(9,"YlOrRd"), space = "Lab")
par(xpd = FALSE,mar=c(0.5, 2.1, 2.1, 9.1))
plot(opt.proj, pch=16)
plot(allr2, col=myColours(22), legend=FALSE, xaxt="n",yaxt="n",add=T)
par(xpd = TRUE)
legend( par()$usr[2]-100000, 2000000, legend=unique(attr(r,"data")@attributes[[1]][,17])[c(1:8,10:23)],fill=myColours(22), cex=0.4)
plot(opt.proj, pch=16, border="white", add=T, lwd=2.0)
# save as region5_bw.pdf

