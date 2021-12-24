# load packages

devtools::install_github("geomorphR/geomorph", ref = "Stable", build_vignettes = TRUE)
library(here)
library(geomorph)
library(tidyverse)
library(wesanderson)

# read GM data
source('readmulti.csv.R')
setwd("./data")
filelist <- list.files(pattern = ".csv")
coords <- readmulti.csv(filelist)
setwd("../")

# read qualitative data
qdata <- read.csv("qdata.csv", 
                  header = TRUE, 
                  row.names = 1)
qdata <- qdata[match(dimnames(coords)[[3]],
                     rownames(qdata)),]

# gpa
Y.gpa <- gpagen(coords, 
                PrinAxes = TRUE, 
                ProcD = TRUE, 
                Proj = TRUE, 
                print.progress = FALSE)

# output + consensus configuration coords
Y.gpa

# combmorph data frame
gdf <- geomorph.data.frame(shape = Y.gpa$coords, 
                           size = Y.gpa$Csize,
                           geo = qdata$geo,
                           time = qdata$time,
                           comb = qdata$comb)

# render 3d gpa plot
#plot(Y.gpa)

# gpa plot
knitr::include_graphics('images/gpa3d.png')

# add centroid size to qdata
qdata$csz <- Y.gpa$Csize

# principal components analysis
pca<-gm.prcomp(Y.gpa$coords)
summary(pca)

# size.geo
fit.size <- procD.lm(size ~ geo * time,
                     data = gdf,
                     print.progress = FALSE,
                     iter = 9999)

# shape.geo
fit.shape <- procD.lm(shape ~ geo * time,
                      data = gdf,
                      print.progress = FALSE,
                      iter = 9999)

# attributes for boxplot
csz <- qdata$csz # centroid size
comb <-  qdata$comb # comb + time

# boxplot of Caddo bottle centroid size by comb (north/south)
csz.comb <- ggplot(qdata, aes(x = comb, y = csz, color = comb)) + 
  geom_boxplot() +
  geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 0.3) +
  scale_colour_manual(values = wes_palette("Moonrise2")) +
  theme(legend.position = "none") +
  labs(x = 'North/South + Formative/Early or Late/Historic?', y = 'Centroid Size')

# render plot
csz.comb

# set plot parameters
pch.gps.comb <- c(15:18)[as.factor(comb)]
col.gps.comb <- wes_palette("Moonrise2")[as.factor(comb)]
col.hull <- c("#C27D38","#798E87","#CCC591","#29211F")

# plot pca by comb
pc.plot <- plot(pca, 
                asp = 1,
                pch = pch.gps.comb,
                col = col.gps.comb)
shapeHulls(pc.plot, 
           groups = comb,
           group.cols = col.hull)

picknplot.shape(pc.plot)

# trajectory analysis::shape
TA <- trajectory.analysis(fit.shape, 
                          groups = qdata$geo, 
                          traj.pts = qdata$time, 
                          print.progress = FALSE)

# magnitude difference
summary(TA, attribute = "MD")

# plot
TP <- plot(TA, 
           pch = as.numeric(qdata$geo), 
           bg = as.numeric(qdata$time),
           cex = 0.9,
           col = "gray")
add.trajectories(TP, traj.pch = c(15, 19),
                 start.bg = 1, 
                 end.bg = 2)


# subset landmark coordinates to produce mean shapes
new.coords <- coords.subset(A = Y.gpa$coords,
                            group = qdata$comb)
names(new.coords)

# group shape means
mean <- lapply(new.coords, mshape)

## plot mean shape north
plot(mean$`north, formative early`)
plot(mean$`north, late historic`)
plot(mean$`south, formative early`)
plot(mean$`south, late historic`)

## comparison plot
plotRefToTarget(mean$`north, formative early`,
                mean$`south, formative early`,
                method = "points",
                mag = 1,
                useRefPts = TRUE)

plotRefToTarget(mean$`north, late historic`,
                mean$`south, late historic`,
                method = "points",
                mag = 1,
                useRefPts = TRUE)
