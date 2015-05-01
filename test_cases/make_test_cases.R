mainDir <- ("/Users/heathersavoy/GitHub/cs267FinalProject/test_cases/")

#Test cases for parallel variogram estimator
library(gstat)
library(sp)
library(gdata)
#1. Specify size
gst <- c("data(small_test): 'test_data.txt', v=1,x=2,y=3;",
         "method: semivariogram;",
         "set fraction=1.0;"
         "set output='test_output.txt';")
bases <- seq(5,50,5)
for(i in 1:length(bases)){
  base <- bases[i]
  width <- base*10
  n <- width*width
  I <- width/10
  mod <- vgm(1,"Exp",I)
  # unconditional simulation
  xy <- expand.grid(1:width, 1:width)
  names(xy) <- c("x","y")
  g.dummy <- gstat(formula = z~1, locations = ~x+y, dummy = TRUE, beta = 0,
                   model = mod, nmax=20)
  yy <- predict(g.dummy, newdata = xy, nsim = 1)
  # show realisation:
  gridded(yy) = ~x+y
  spplot(yy)
  vg <- variogram(sim1~1,data=yy,cutoff=sqrt(2*width^2))
  # print field to input file
  subDir <- paste("n",n,sep="")
  dir.create(file.path(mainDir, subDir), showWarnings = FALSE)
  setwd(file.path(mainDir, subDir))
  write.fwf(vg[,1:3],"correct.txt",
            colnames=FALSE,rownames=FALSE)
  ofile <- "test_data.txt"
  write.fwf(as.data.frame(yy),ofile,
              colnames=FALSE,rownames=FALSE)
  write.table(gst,"test.gst",col.names=FALSE,row.names=FALSE,quote=FALSE)
}

