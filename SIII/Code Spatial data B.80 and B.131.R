library(raster)
library(gstat)
library(viridisLite)
library(automap)


#...........................BUILDING 80............................


#.............Preparation..............
# Get data
DataB80<- read.csv("Spatial data B.80.csv")
coordinates(DataB80) <- ~X+Y

# grid to interpolate values
r <- raster(extent(DataB80))
res(r) <- 0.05
r.grd <- as(r, "SpatialGrid")

# shape of the sampled area
shp <- shapefile("MaskB80.shp")


#.................B.80 Spatial analysis Starch Grains (variograms)..............


#1) Triticeae Starch variogram

var1 <- variogram(TriticeaeST~1, DataB80, cutoff=2.5, width=0.5)
fve1 <- fit.variogram(var1, vgm(1,c("Exp", "Gau", "Sph"), 1.2, .04))
tritiB80<- plot(var1, fve1, main= "Triticeae Starch", cex.main=1)
plot(tritiB80)


#2) Panicoideae starch variogram

var2 <- variogram(PanicoideaeST~1, DataB80, cutoff=4.5, width=0.6)
fve2<- fit.variogram(var2, vgm(180, "Exp", 1, 0.01))
panistB80<- plot(var2,fve2, main= "Panicoideae Starch", cex.main=1)
plot(panistB80)

#3) Faboideae starch variogram

var3 <- variogram(FaboideaeST~1, DataB80, cutoff=5, width=0.8)
fve3 <- fit.variogram(var3, vgm(0.25, "Exp", .4, 0.1))
faboB80<-plot(var3, fve3, main= "Faboideae Starch", cex.main=1)
plot(faboB80)

#4) USO starch variogram

var4<- variogram(USO~1, DataB80, cutoff=3, width=0.5)
fve4<- fit.variogram(var4, vgm(3.4, "Exp", 0.4, 1))
USO80<-plot(var4, fve4, main= "USO Starch", cex.main=1)
plot(USO80)

#5) USO Group A starch variogram

var5 <- variogram(USOA~1, DataB80, cutoff=2.5, width=0.5)
fve5<- fit.variogram(var5, vgm(.50, "Exp", .5, .01))
USOA<-plot(var5, fve5, main= "USO Group A", cex.main=1)
plot(USOA)


#6) USO Group B starch variogram

var6 <- variogram(USOB~1, DataB80, cutoff=2, width=0.3)
fve6<- fit.variogram(var6, vgm(0.34, "Exp", .50, 0));plot(var6,fve6)
USOB<-plot(var6, fve6, main= "USO Group B", cex.main=1)
plot(USOB)

#7) USO Group C starch variogram

var7 <- variogram(USOC~1, DataB80, cutoff=5, width= 1)
fve7<- fit.variogram(var7, vgm(1.35, "Exp", 0.50, 0))
USOC<-plot(var7, fve7, main= "USO Group C", cex.main=1)
plot(USOC)

#.................B.80 Spatial analysis Starch Grains (Kriging)..............



#1) Krigign Triticeae Starch

krg1<- krige(TriticeaeST~1, DataB80, r.grd, fve1)
pred1<- mask(raster(krg1), shp)
var1<- mask(raster(krg1, "var1.var"), shp)
plot(pred1, main= "Triticeae Starch", cex.main= 1,col= turbo (100), axes= F)
plot(var1, main= "Triticeae Starch variance", cex.main= 1,col= mako (100), axes= F)


#2) Kriging Panicoideae Starch

krg2<- krige(PanicoideaeST~1, DataB80, r.grd, fve2)
pred2<- mask(raster(krg2), shp)
var2<- mask(raster(krg2, "var1.var"), shp)
plot(pred2, main= "Panicoideae Starch", cex.main= 1,col= turbo (100), axes= F)
plot(var2, main= "Panicoideae Starch variance", cex.main= 1,col= mako (100), axes= F)


#3) Kriging Faboideae starch

krg3<- krige(FaboideaeST~1, DataB80, r.grd, fve3)
pred3<- mask(raster(krg3), shp)
var3<- mask(raster(krg3, "var1.var"), shp)
plot(pred3, main= "Faboideae Starch", cex.main= 1,col= turbo (100), axes= F)
plot(var3, main= "Faboideae Starch variance", cex.main= 1,col= mako (100), axes= F)


#4) Kriging USO

krg4<- krige(USO~1, DataB80, r.grd, fve4)
pred4<- mask(raster(krg4), shp)
var4<- mask(raster(krg4, "var1.var"), shp)
plot(pred4, main= "USO Starch", cex.main= 1,col= turbo (100), axes= F)
plot(var4, main= "USO Starch variance", cex.main= 1,col= mako (100), axes= F)


#5) Kriging USO Group A

krg5<- krige(USOA~1, DataB80, r.grd, fve5)
pred5<- mask(raster(krg5), shp)
var5<- mask(raster(krg5, "var1.var"), shp)
plot(pred5, main= "USO Group A", cex.main= 1,col= turbo (100), axes= F)
plot(var5, main= "USO Group A variance", cex.main= 1,col= mako (100), axes= F)


#6) Kriging USO Group B

krg6<- krige(USOB~1, DataB80, r.grd, fve6)
pred6<- mask(raster(krg6), shp)
var6<- mask(raster(krg6, "var1.var"), shp)
plot(pred6, main= "USO Group B", cex.main= 1,col= turbo (100), axes= F)
plot(var6, main= "USO Group B variance", cex.main= 1,col= mako (100), axes= F)




#7) Kriging USO Group C

krg7<- krige(USOC~1, DataB80, r.grd, fve7)
pred7<- mask(raster(krg7), shp)
var7<- mask(raster(krg7, "var1.var"), shp)
plot(pred7, main= "USO Group C", cex.main= 1,col= turbo (100), axes= F)
plot(var7, main= "USO Group C variance", cex.main= 1,col= mako (100), axes= F)

#Export kriging prediction as rasters.


writeRaster(pred1, "Triticeae Starch 80.asc")
writeRaster(pred2, "Faboideae Starch 80.asc")
writeRaster(pred3, "Panicoideae Starch 80.asc")
writeRaster(pred4, "USO Starch 80.asc")
writeRaster(pred5, "USO Group A.asc")
writeRaster(pred6, "USO Group B.asc")
writeRaster(pred7, "USO Group C.asc")



#.................B.80 Spatial analysis Phytoltihs (Variograms)..............



#8) Infloresence phytoltiths variogram

var8<- variogram(Grass.inflorescence~1, DataB80, cutoff=4, width= 0.4)
fve8<- fit.variogram(var8, vgm(650, "Exp", 0.50, 300))
Inflo80<- plot(var8, fve8, main= "Inflorescence phytoliths", cex.main=1)
plot(Inflo80)


#9) Grass.Culm.leaves phytoltiths variogram

var9<- variogram(Grass.Culm.Leaves~1, DataB80, cutoff=5, width= 0.7)
fve9<- fit.variogram(var9, vgm(900, "Exp", 1.30, 150))
Grassculm80<- plot(var9, fve9, main= "Grass Culm.leaves phytoliths", cex.main=1)
plot(Grassculm80)


#10) Panicoideae phytoltiths variogram

var10<- variogram(Panicoideae.grass~1, DataB80, cutoff=5.5)
fve10<- fit.variogram(var10, vgm(14, "Exp", 2, 0))
PanicoPh80<-plot(var10, fve10, main= "Panicoideae phytoliths", cex.main=1)
plot(PanicoPh80)


#11) Arecaceae phytoltiths variogram

var11<- variogram(Arecaceae~1, DataB80, cutoff=4.5)
fve11<- fit.variogram(var11, vgm(0.20, "Exp", .30, 0))
Palm80<-plot(var11, fve11, main= "Arecaceae phytoliths", cex.main=1)
plot(Palm80)


#.................B.80 Spatial analysis Phytoltihs (Kriging)..............



#8) Kriging Inflorescence

krg8<- krige(Grass.inflorescence~1, DataB80, r.grd, fve8)
pred8<- mask(raster(krg8), shp)
var8<- mask(raster(krg8, "var1.var"), shp)
plot(pred8, main= "Grass.inflorescence", cex.main= 1,col= turbo (100), axes= F)
plot(var8, main= "Grass.inflorescence phytoliths variance", cex.main= 1,col= mako (100), axes= F)


#9) kriging Culms.leaves

krg9 <- krige(Grass.Culm.Leaves~1, DataB80, r.grd, fve9)
pred9<- mask(raster(krg9), shp)
var9<- mask(raster(krg9, "var1.var"), shp)
plot(pred9, main= "Grass culm.leaves", cex.main= 1,col= turbo (100), axes= F)
plot(var9, main= "Grass culm.leaves phytoliths variance", cex.main= 1,col= mako (100), axes= F)


#10) Kriging Panicoideae Phyto

krg10 <- krige(Panicoideae.grass~1, DataB80, r.grd, fve10)
pred10<- mask(raster(krg10), shp)
var10 <- mask(raster(krg10, "var1.var"), shp)
plot(pred10, main= "Panicoideae Phytoliths", cex.main= 1,col= turbo (100), axes= F)
plot(var10, main= "Panicoideae Phytoliths variance", cex.main= 1,col= mako (100), axes= F)


#11) Kriging Arecaceae 

krg11<- krige(Arecaceae~1, DataB80, r.grd, fve11)
pred11<- mask(raster(krg11), shp)
var11<- mask(raster(krg11, "var1.var"), shp)
plot(pred11, main= "Arecaceae", cex.main= 1,col= turbo (100), axes= F)
plot(var11, main= "Arecaceae variance", cex.main= 1,col= mako (100), axes= F)


#..............Export kriging prediction as rasters............

writeRaster(pred8, "Inflorescence 80.asc")
writeRaster(pred9, "Grass.culm 80.asc")
writeRaster(pred10, "Panicoideae Phyto 80.asc")
writeRaster(pred11, "Arecaceae 80.asc")






#...........................BUILDING 131............................



#.............Preparation..............

# Get data
DataB131<- read.csv("Spatial data B.131.csv")
coordinates(DataB131) <- ~X+Y

# grid to interpolate values to
r131 <- raster(extent(DataB131))
res(r131) <- 0.05
r.grd131 <- as(r131, "SpatialGrid")

# shape of the sampled area
shp131 <- shapefile("MaskB131.shp")


#.................B.131 Spatial analysis Starch Grains (variograms)..............



#12) Triticeae starch variogram

var12 <- variogram(TriticeaeST~1, DataB131, cutoff=2.5, width=0.4)
fve12<- fit.variogram(var12, vgm(90, "Exp", .5, 0))
Triti131<- plot(var12, fve12, main= "Triticeae starch", cex.main=1)
plot(Triti131)



#13) Panicoideae starch variogram

var13<- variogram(PanicoideaeST~1, DataB131, cutoff=2.5, width=0.3)
fve13<- fit.variogram(var13, vgm(.70, "Exp", .30, 0))
Panist131<-plot(var13, fve13, main= "Panicoideae starch", cex.main=1)
plot(Panist131)


#14) USO starch variogram

var14<- variogram(USO~1, DataB131, cutoff=2.5, width= 0.6)
fve14<- fit.variogram(var14, vgm(10, "Exp", .34, 3))
USO131<-plot(var14, fve14, main= "USO starch", cex.main=1)
plot(USO131)


#.................B.131 Spatial analysis Starch Grains (Kriging/IDW)..............



#12) Kriging Triticeae

krg12<- krige(TriticeaeST~1, DataB131, r.grd131,fve12)
pred12<- mask(raster(krg12), shp131)
var12<- mask(raster(krg12, "var1.var"), shp131)
plot(pred12, main= "Triticeae starch", cex.main= 1,col= turbo (100), axes= F)
plot(var12, main= "Triticeae starch variance", cex.main= 1,col= mako (100), axes= F)


#13) Kriging Panicoideae

krg13<- krige(PanicoideaeST~1, DataB131, r.grd131,fve13)
pred13<- mask(raster(krg13), shp131)
var13<- mask(raster(krg13, "var1.var"), shp131)
plot(pred13, main= "Panicoideae starch", cex.main= 1,col= turbo (100), axes= F)
plot(var13, main= "Panicoideae starch variance", cex.main= 1,col= mako (100), axes= F)


#14) Kriging USO

krg14<- krige(USO~1, DataB131, r.grd131,fve14)
pred14<- mask(raster(krg14), shp131)
var14<- mask(raster(krg14, "var1.var"), shp131)
plot(pred14, main= "USO starch", cex.main= 1,col= turbo (100), axes= F)
plot(var14, main= "USO starch variance", cex.main= 1,col= mako (100), axes= F)


#15) IDW Faboideae starch

gs15<- gstat(formula=FaboideaeST~1, location= DataB131,nmax = 5)
idw15<- interpolate(r131,gs15)
idw15<- mask(idw15, shp131)
plot(idw15, main= "Faboideae starch", cex.main= 1,col= turbo (100),axes= F)



#......Export kriging prediction as rasters.....

writeRaster(pred12, "Triticeae Starch 131.asc")
writeRaster(pred13, "Panicoideae Starch 131.asc")
writeRaster(pred14, "USO Starch 131.asc")
writeRaster(idw15, "Faboideae Starch 131.asc")




#.................B131 Spatial analysis Phytoltihs (Variograms)..............



#16) Panicoideae phytolith variogram

var16<- variogram(Panicoideae.grass~1, DataB131, cutoff=2.5, width=0.2)
fve16<- fit.variogram(var16, vgm(2.8, "Exp", .28,0))
Paniph131<-plot(var16, fve16, main= "Panicoideae phytoliths", cex.main=1)
plot(Paniph131)


#17) Arecaceae phytolith variogram

var17<- variogram(Arecaceae~1, DataB131, cutoff=3, width= 0.4)
fve17<- fit.variogram(var17, vgm(0.14, "Exp", 4.10, 0.09))
Arec131<-plot(var17, fve17, main= "Arecaceae phytoliths", cex.main=1)
plot(Arec131)

#.................B.131 Spatial analysis Phytoliths (Kriging/IDW)..............


#16) Kriging Panicoideae phytoltihs

krg16<- krige(Panicoideae.grass~1, DataB131, r.grd131,fve16)
pred16<- mask(raster(krg16), shp131)
var16<- mask(raster(krg16, "var1.var"), shp131)
plot(pred16, main= "Panicodeae phytoliths", cex.main= 1,col= turbo (100), axes= F)
plot(var16, main= "Panicodeae phytoliths variance", cex.main= 1,col= mako (100), axes= F)


#17) Kriging Arecaceae

krg17<- krige(Arecaceae~1, DataB131, r.grd131,fve17)
pred17<- mask(raster(krg17), shp131)
var17<- mask(raster(krg17, "var1.var"), shp131)
plot(pred17, main= "Arecaceae phytoliths", cex.main= 1,col= turbo (100), axes= F)
plot(var17, main= "Arecaceae phytoliths variance", cex.main= 1,col= mako (100), axes= F)



#18) IDW Grass inflorescence

Inf18<- gstat(formula=Grass.inflorescence~1, location=DataB131,nmax = 5)
idw18<- interpolate(r131,Inf18)
idw18<- mask(idw18, shp131)
plot(idw18,main= "Inflorescence", cex.main= 1,col= turbo (100),axes= F)


#19) IDW Grass culm/leaves  

gs19<- gstat(formula=Grass.Culm.Leaves~1, location=DataB131,nmax = 5)
idw19<- interpolate(r131,gs19)
idw19<- mask(idw19, shp131)
plot(idw19,main= "Grass culm.leaves", cex.main= 1,col= turbo (100),axes= F)



#..........Export kriging prediction as rasters.....


writeRaster(pred16, "Panicoideae Phyto 131.asc")
writeRaster(pred17, "Arecaceae 131.asc")
writeRaster(idw18, "Inflorescence 131.asc")
writeRaster(idw19, "Grass.culm 131.asc")