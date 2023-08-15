
# Load packages and source code
devtools::install_github("RGLab/flowCore")

4
install.packages("RcppArmadillo")
library("flowCore") # for putting data into a usable form
library("ggplot2") # for plotting
library("ggcyto") # for plotting

library("tidyverse") # for data manipulation
library("knitr") 
library("dplyr") # for data manipulation

library("flowViz") # for plotting
library("Phenoflow") # for fingerprinting
library("flowAI") # for denoising

# Set a fixed seed to ensure reproducible analysis
set.seed(777)

# load data
fs <- read.flowSet(path = "C:/Users/mirte/Documents/UCSD/Data/Flow Cytometry/All of Nikki's samples in one place", pattern = ".fcs",alter.names = T)

# retrieve and alter sample names
pData(fs)

pData(fs)$well <- gsub(".*_(.*)_.*.fcs","\\1",sampleNames(fs)) # extract well from name and add new 'well' column
pData(fs) # check successful

# change channel names to something more appropriate
colnames(fs)

colnames(fs)[colnames(fs)=="FITC.A"] <- "SYBRGreen.A" # this is the FL that stains DNA in both algae and bacteria
colnames(fs)[colnames(fs)=="Pacific.Blue.A"] <- "425-475nm.A"
colnames(fs)[colnames(fs)=="PerCP.Cy5.5.A"] <- "670-735nm.A" # this picks up chlorophyll autofluourescense strongly, and mitochondrial autofluourescense weakly

colnames(fs)[colnames(fs)=="FITC.W"] <- "SYBRGreen.W"
colnames(fs)[colnames(fs)=="Pacific.Blue.W"] <- "425-475nm.W"
colnames(fs)[colnames(fs)=="PerCP.Cy5.5.W"] <- "670-735nm.W"

colnames(fs)[colnames(fs)=="FITC.H"] <- "SYBRGreen.H"
colnames(fs)[colnames(fs)=="Pacific.Blue.H"] <- "425-475nm.H"
colnames(fs)[colnames(fs)=="PerCP.Cy5.5.H"] <- "670-735nm.H"

colnames(fs)


### Extract metadata from sample names
metadata <- data.frame(do.call(rbind, lapply(strsplit(flowCore::sampleNames(fs),"_"), rbind)))
colnames(metadata) <- c("Sample Replicate", "Run Replicate", "Sample", "ID within run")

# find non-samples which will mess with some analyses
ind <- as.data.frame(pData(fs))
ind <- cbind(ind, seq(1, length(metadata$Sample), 1))
colnames(ind) <- c("File", "Sample", "Unique")

not.samples <- as.vector(c(1:5, 54:61, 90:92))

## Make a second set of metadata for those that cannot include samples
metadata2 <- metadata %>%
  slice(-not.samples)

# Do same for data as for metadata
attributes(fs)

sfs <- fs[-not.samples]

attributes(sfs)


### block change ####


# Select phenotypic features of interest and transform parameters
fs_transformed <- transform(fs,`SYBRGreen.A`=asinh(`SYBRGreen.A`), `SYBRGreen.H`=asinh(`SYBRGreen.H`), `SYBRGreen.W`=asinh(`SYBRGreen.W`), 
                            `SSC.A`=asinh(`SSC.A`), `SSC.W`=asinh(`SSC.W`), `SSC.H`=asinh(`SSC.H`), 
                            `425-475nm.A`=asinh(`425-475nm.A`), `425-475nm.H`=asinh(`425-475nm.H`), `425-475nm.W`=asinh(`425-475nm.W`), 
                            `670-735nm.A` = asinh(`670-735nm.A`),`670-735nm.H` = asinh(`670-735nm.H`),`670-735nm.W` = asinh(`670-735nm.W`), 
                            `FSC.A`=asinh(`FSC.A`), `FSC.H`=asinh(`FSC.H`), `FSC.W`=asinh(`FSC.W`),
                            `PE.A`=asinh(`PE.A`), `PE.H`=asinh(`PE.H`), `PE.W`=asinh(`PE.W`),
                            `PE.Cy7.A`=asinh(`PE.Cy7.A`), `PE.Cy7.H`=asinh(`PE.Cy7.H`), `PE.Cy7.W`=asinh(`PE.Cy7.W`),
                            `APC.A`=asinh(`APC.A`), `APC.H`=asinh(`APC.H`), `APC.W`=asinh(`APC.W`),
                            `APC.Cy7.A`=asinh(`APC.Cy7.A`), `APC.Cy7.H`=asinh(`APC.Cy7.H`), `APC.Cy7.W`=asinh(`APC.Cy7.W`),
                            `AmCyan.A`=asinh(`AmCyan.A`), `AmCyan.H`=asinh(`AmCyan.H`), `AmCyan.W`=asinh(`AmCyan.W`))

sfs_transformed <- transform(sfs,`SYBRGreen.A`=asinh(`SYBRGreen.A`), `SYBRGreen.H`=asinh(`SYBRGreen.H`), `SYBRGreen.W`=asinh(`SYBRGreen.W`),
                             `SSC.A`=asinh(`SSC.A`), `SSC.W`=asinh(`SSC.W`), `SSC.H`=asinh(`SSC.H`), 
                             `425-475nm.A`=asinh(`425-475nm.A`), `425-475nm.H`=asinh(`425-475nm.H`), `425-475nm.W`=asinh(`425-475nm.W`), 
                             `670-735nm.A` = asinh(`670-735nm.A`),`670-735nm.H` = asinh(`670-735nm.H`),`670-735nm.W` = asinh(`670-735nm.W`),
                             `FSC.A`=asinh(`FSC.A`), `FSC.H`=asinh(`FSC.H`), `FSC.W`=asinh(`FSC.W`),
                             `PE.A`=asinh(`PE.A`), `PE.H`=asinh(`PE.H`), `PE.W`=asinh(`PE.W`),
                             `PE.Cy7.A`=asinh(`PE.Cy7.A`), `PE.Cy7.H`=asinh(`PE.Cy7.H`), `PE.Cy7.W`=asinh(`PE.Cy7.W`),
                             `APC.A`=asinh(`APC.A`), `APC.H`=asinh(`APC.H`), `APC.W`=asinh(`APC.W`),
                             `APC.Cy7.A`=asinh(`APC.Cy7.A`), `APC.Cy7.H`=asinh(`APC.Cy7.H`), `APC.Cy7.W`=asinh(`APC.Cy7.W`),
                             `AmCyan.A`=asinh(`AmCyan.A`), `AmCyan.H`=asinh(`AmCyan.H`), `AmCyan.W`=asinh(`AmCyan.W`))

# make a vector of the transformed column names
param <- c("SYBRGreen.A","425-475nm.A","SSC.A","FSC.A", "670-735nm.A", "PE.A", "PE.Cy7.A", "APC.A", "APC.Cy7.A", "AmCyan.A",
           "SYBRGreen.H","425-475nm.H","SSC.H","FSC.H", "670-735nm.H", "PE.H", "PE.Cy7.H", "APC.H", "APC.Cy7.H", "AmCyan.H",
           "SYBRGreen.W","425-475nm.W","SSC.W","FSC.W", "670-735nm.W", "PE.W", "PE.Cy7.W", "APC.W", "APC.Cy7.W", "AmCyan.W")


#Create a PolygonGate for denoising the dataset
### Create a PolygonGate for denoising the dataset

#                     bottom left; top left, top middle left, top middle right 
mB <- as.matrix(rbind(c(3, 5), c(5.25, 5), c (6.25, 8.5), c(6.25,10),
                      c(7, 12), c(4,10.5), c(3,9)))
#                     top right, bottom right, bottom middle
colnames(mB) <- c("425-475nm.A", "SYBRGreen.A")
polyGate <- polygonGate(.gate=mB, filterId = "Bacterial Cells")


xyplot(`425-475nm.A` ~ `SYBRGreen.A`, data=fs_transformed[11],
       scales=list(y=list(limits=c(2,14)),
                   x=list(limits=c(4,14))),
       axis = axis.default, nbin=200, filter=polyGate,
       checkName = FALSE, par.strip.text=list(col="white", font=2, cex=2), smooth=FALSE)
#in the above plot try at a minimum: 11, 35, 66, 78, 114, 129, 161

xyplot(`425-475nm.A` ~ `SYBRGreen.A`, data=fs_transformed[1],
       scales=list(y=list(limits=c(2,14)),
                   x=list(limits=c(4,14))),
       axis = axis.default, nbin=200, filter=polyGate,
       checkName = FALSE, par.strip.text=list(col="white", font=2, cex=2), smooth=FALSE)

xyplot(`425-475nm.A` ~ `SYBRGreen.A`, data=fs_transformed[2],
       scales=list(y=list(limits=c(2,14)),
                   x=list(limits=c(4,14))),
       axis = axis.default, nbin=200, filter=polyGate,
       checkName = FALSE, par.strip.text=list(col="white", font=2, cex=2), smooth=FALSE)

xyplot(`425-475nm.A` ~ `SYBRGreen.A`, data=fs_transformed[55],
       scales=list(y=list(limits=c(2,14)),
                   x=list(limits=c(4,14))),
       axis = axis.default, nbin=200, filter=polyGate,
       checkName = FALSE, par.strip.text=list(col="white", font=2, cex=2), smooth=FALSE)

xyplot(`425-475nm.A` ~ `SYBRGreen.A`, data=fs_transformed[3],
       scales=list(y=list(limits=c(2,14)),
                   x=list(limits=c(4,14))),
       axis = axis.default, nbin=200, filter=polyGate,
       checkName = FALSE, par.strip.text=list(col="white", font=2, cex=2), smooth=FALSE)

xyplot(`425-475nm.A` ~ `SYBRGreen.A`, data=fs_transformed[90],
       scales=list(y=list(limits=c(2,14)),
                   x=list(limits=c(4,14))),
       axis = axis.default, nbin=200, filter=polyGate,
       checkName = FALSE, par.strip.text=list(col="white", font=2, cex=2), smooth=FALSE)

xyplot(`425-475nm.A` ~ `SYBRGreen.A`, data=fs_transformed[4],
       scales=list(y=list(limits=c(2,14)),
                   x=list(limits=c(4,14))),
       axis = axis.default, nbin=200, filter=polyGate,
       checkName = FALSE, par.strip.text=list(col="white", font=2, cex=2), smooth=FALSE)


# Can come back and try improve these further, once you get the full job working.


### Isolate only the cellular information based on the polyGate1
fs_subset <- Subset(fs_transformed, polyGate)
sfs_subset <- Subset(sfs_transformed, polyGate)

### Check gate

xyplot(`670-735nm.A` ~ `SYBRGreen.A`, data=fs_subset[100],
       scales=list(y=list(limits=c(0,14)),
                   x=list(limits=c(0,13))),
       axis = axis.default, nbin=200, filter=mB, checkName = FALSE, par.strip.text=list(col="white", font=2, cex=2), smooth=FALSE)

xyplot(`670-735nm.A` ~ `SYBRGreen.A`, data=fs_subset[100],
       scales=list(y=list(limits=c(0,14)),
                   x=list(limits=c(0,15))),
       axis = axis.default, nbin=200, checkName = FALSE, par.strip.text=list(col="white", font=2, cex=2), smooth=FALSE)

xyplot(`SSC.A` ~ `FSC.A`, data=fs_transformed[100],
       scales=list(y=list(limits=c(0,14)),
                   x=list(limits=c(0,15))),
       axis = axis.default, nbin=200, checkName = FALSE, par.strip.text=list(col="white", font=2, cex=2), smooth=FALSE)









