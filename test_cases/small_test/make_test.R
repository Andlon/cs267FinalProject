#Take a slice of herten (horizontal so fairly isotropic and 2D)
height <- sort(unique(coordinates(domain)[,3]))[round(length(z_inc)/2)] #about halfway
slice <- domain$lnK[coordinates(domain)[,3]==height]
  
#View slice
test_grid <- expand.grid(x=unique(coordinates(domain)[,1]),
                         y=unique(coordinates(domain)[,2]))
coordinates(test_grid) <- ~x+y
test_field <- data.frame(lnK=slice)
coordinates(test_field) <- coordinates(test_grid)
spplot(test_field,col.regions=heat.colors(5))
#image.plot(x_inc,y_inc,t(matrix(slice,ncol=length(xax),byrow=TRUE)))

#Check variogram is coherent
plot(variogram(lnK~1,data=test_field))

#Print slice to gstat format
#no header; x y v
ofile <- "C:/cs267FinalProject/small_test/test_data.txt"
write.table(as.data.frame(test_field),ofile,col.names=FALSE,row.names=FALSE)
