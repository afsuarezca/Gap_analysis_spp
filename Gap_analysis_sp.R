#__________________________________#
# Code written by Felipe Suarez    #
# Version :  01-09-2020            #
#__________________________________#

#clear workspace
rm(list=ls())

library(SpaDES)
library(raster)
library(xfun)
library(rgdal)
library(ReIns)

#Read rasters species data####

setwd("~/rasters_aves/aves-shp/aves_filter_elev_eco_prj")

myfiles = list.files(pattern="*.tif") 

list_raster<-list()
for (i in 1:length(myfiles)){
  print(file_ext(myfiles[[i]]))
  if(file_ext(myfiles[[i]]) == "xml"|file_ext(myfiles[[i]]) == "dbf"
     |file_ext(myfiles[[i]]) == "ovr"|file_ext(myfiles[[i]]) == "cpg"
     |file_ext(myfiles[[i]]) == "tfw"|file_ext(myfiles[[i]]) == "lock"){
    list_raster[[i]]<-NULL
  }else{
    list_raster[[i]]<-myfiles[[i]]
  }
}

list_raster<-plyr::compact(list_raster)


#Read RUNAP data multiple years (1970,1990,2000,2010,2020)
setwd("~/RUNAP/RUNAP_rasters")

myfiles = list.files(pattern="*.tif") 

list_raster_PNN<-list()
for (i in 1:length(myfiles)){
  print(file_ext(myfiles[[i]]))
  if(file_ext(myfiles[[i]]) == "xml"|file_ext(myfiles[[i]]) == "dbf"
     |file_ext(myfiles[[i]]) == "ovr"|file_ext(myfiles[[i]]) == "cpg"
     |file_ext(myfiles[[i]]) == "tfw"|file_ext(myfiles[[i]]) == "lock"){
    list_raster_PNN[[i]]<-NULL
  }else{
    list_raster_PNN[[i]]<-myfiles[[i]]
  }
}

list_raster_PNN<-plyr::compact(list_raster_PNN)

#FOR 1990 YOU HAVE TO TRANSFORM THE 128 DATA TO NA

#Function to calculate area to protect per sp####

TempDF = data.frame(x=c(1000,2.5e5), y=c(100,10))
m = diff(TempDF$y)/diff(TempDF$x)
b = ((TempDF$y[[1]]*TempDF$x[[2]])-(TempDF$y[[2]]*TempDF$x[[1]]))/(TempDF$x[[2]]-TempDF$x[[1]])
lm.func<-function(x,m,b){
  y = m*x + b
  return(y)
}
x <- seq(1000, 250000, 1000)
plot(x, lm.func(x,m,b), xlab="x", ylab="PDF", type="l")

sp_stats<-list()#list to save stats per species
PNN_allsp_stats<-list()# list to save stats for all years

for (h in c(1)) {#list of PNN (RUNAP) rasters
  setwd("~/RUNAP/RUNAP_rasters")
  PNN<-raster(list_raster_PNN[[h]])
  #PNN[which(PNN[]>126)]<-NA#ONLY FOR 1990 DATA!
  for (j in 1:length(list_raster)){#list_pol is the list with all the species
    setwd("~/rasters_aves/aves-shp/aves_filter_elev_eco_prj")
    sp_focus<-raster(list_raster[[j]])#read each sp raster
    if (is.null(intersect(extent(PNN), sp_focus))){#check if raster and shapefile (each sp) intersect
      dfsp<-as.data.frame(matrix(nrow = 1,ncol = 3))
      colnames(dfsp)<-c("PNN","freq")
      dfsp$species<-names(sp_focus)
      dfsp$year_PNN<-strsplit(names(PNN),"_")[[1]][[2]]
    }else{
      df_sp<-as.data.frame(rasterToPoints(sp_focus))#extract values raster (dataframe with x,y,value)
      df_sp<-subset(df_sp,df_sp[3]==1)#subset df_sp with species occurrence (= 1)
      PNN_data<-raster::extract(PNN,df_sp[1:2])#extract data of PAs
      df_sp$PNN<-PNN_data#add column with PAs_data
      if(nrow(df_sp)>0){
        df_sp<-as.data.frame(df_sp)#convert to data frame
        #df_sp<-subset(df_sp,!is.na(df_sp[4]))#remove NAs
        df_sp<-plyr::count(df_sp[4])#count number of pixels inside and outside PAs
        df_sp$species<-names(sp_focus)#add column with sp name
        df_sp$year_PNN<-strsplit(names(PNN),"_")[[1]][[2]] #add year of PAs network
      }
    }
    #summarize data inside and outside PAs 
    df_sp[which(is.na(df_sp$PNN)),"in_PNN"]<-0 
    df_sp[which(!is.na(df_sp$PNN)),"in_PNN"]<-1
    count_pnn<-aggregate(df_sp$freq,list(PNN=df_sp$in_PNN),"sum")
    total_area<-sum(count_pnn$x)
    total_area_km2<-total_area*((300/1000)^2)#transform to area
    total_protected<-count_pnn[which(count_pnn$PNN == "1"),2]
    
    #calculate the proportion of the sp distribution area that is protected
    if(length(total_protected)>0){
      df_sp$per_protected<-(total_protected/total_area)*100
    }else{
      df_sp$per_protected<-0
    }
    #calculate representativeness target
    if(total_area_km2<=1000){
      df_sp$needed_to_protect<-100
    }else{
      if(total_area_km2>=250000){
        df_sp$needed_to_protect<-10 
      }else{
        df_sp$needed_to_protect<-lm.func(total_area_km2,m,b)
      }
    }
    df_sp$total_area_km2<-total_area_km2
    df_sp$per_achieved<-df_sp$per_protected/df_sp$needed_to_protect
    if(df_sp$per_achieved[[1]]>0.9){
      df_sp$achieved<-1
    }else{
      df_sp$achieved<-0
    }
    print(paste(j))
    sp_stats[[j]]<-df_sp[1,]
    gc()
  }
  PNN_allsp_stats[[h]]<-do.call(rbind,sp_stats)
}

write.csv(do.call(rbind,PNN_allsp_stats),"statistics_gap_analysis_sp.csv")