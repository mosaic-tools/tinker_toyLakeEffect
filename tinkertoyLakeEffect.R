# script to use tinker_toy model and Savannah river to explore influence network structure (location of lakes, number of lakes, etc.) on biogeochemistry

# source tinker_toy functions
source("tinkerToy_analytical.R")

# load stream-only network (Savannah river out of NHD)
NHDout=read.csv("savannah_river_v2.csv",stringsAsFactors=FALSE)

# scaling relationship between discharge and width from Raymond et al. 2012 "scaling the gas transfer velocity and hydraulic geometry on streams and small rivers"
a=12.88
b=0.42

# scaling relationship between discharge and depth from Lauren
c=0.408
d=0.294

# reachNumber arbitrary ordering from NHD output
# reachID from NHD (COMID)
# length from NHD (LENGTHKM)
# width from hydraulic geometry
# depth from (discharge/velocity)/(width from hydraulic geometry)
# localArea from NHD (AreaSqKM)
# localPET standardized across catchment and 
#     * google says 35-year average annual precip is 43.28 in or 1.099 m
#     * ET = 1 m/yr best guess from maps on google
#     * daily PET is (1.099-1)/365 = 0.0003 m/d
#     Annual Estimates of Recharge, Quick-Flow Runoff, and Evapotranspiration for the Contiguous U.S. Using Empirical Regression Equations
#     Reitz et al. 2017
# localCin standardized across catchment and arbitrary
# connect derived from information provided by John Gardner (fromCOMID)
# order from NHD for inserting lakes

#  convert NHD output to network dataframe for  use in tinker_toy; note  we add order and lake fields for use with model experiments
network=data.frame(reachNumber=1:nrow(NHDout), 
                  reachID=NHDout$COMID, 
                  length=NHDout$LENGTHKM*1000, 
                  width=a*(NHDout$QE_MA*0.0283)^b, 
                  depth=c*(NHDout$QE_MA*0.0283)^d, 
                  localArea=NHDout$AreaSqKM*1e6, 
                  localPET=0.005, 
                  localCin=10, 
                  d=0.1, 
                  connect=NHDout$fromCOMID,
                  order=NHDout$StreamOrde,
                  lake=0)

# removing some artifacts from how John created the contributing reaches
network$connect=gsub("c(","",network$connect,fixed=TRUE)
network$connect=gsub(")","",network$connect,fixed=TRUE)
network$connect=gsub(" ","",network$connect,fixed=TRUE)
network$connect=gsub(":",",",network$connect,fixed=TRUE)
network$connect[is.na(network$connect)]=0

# fill in NHD reaches with  a discharge (QE_MA) of 0
meanWidth_firstorder=mean(network$width[network$order==1 & network$width>0])
meanDepth_firstorder=mean(network$depth[network$order==1 & network$depth>0])
network$width[network$width==0]=meanWidth_firstorder
network$depth[network$depth==0]=meanDepth_firstorder

# sort network
#sorted=sortNetwork(network=network,counterMax=1000,verbose=TRUE)
#write.csv(sorted,"SORTED_savannah_river_v2.csv",row.names=FALSE)
sorted=read.csv("SORTED_savannah_river_v2.csv",header=TRUE,stringsAsFactors=FALSE)

# solve equilibrium for base network
baseProcess=solveNetwork(network=sorted)

#*******************************#
#   Running model experiments   #
#*******************************#

# add 40 lakes to 6th order streams
lakes40=addRandomLake(network=sorted,order=rep(6,40),width=10000,depth=10,lake_d=0.005,verbose=TRUE)

# solve equilibrium for base network
lakes40Process=solveNetwork(network=lakes40)

