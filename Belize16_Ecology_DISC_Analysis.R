          ### Belize 2016 Population-Level DISC Analysis: Ecological Context ###



## NOTE!!! to search for publication plots ctrl + f "PUBLICATION PLOT"
# to search for supmat plots ctrl + f "SUPMAT PLOT"


# second pass ideas:
#
# - a model to see what predicts behaviour_out
# - kinematics variation instead of means
# - bearing variable
# - treat drift direction as a circular variable for circular-circular correlations




## Set Up Shop ------------------------

# set working directory 

# setwd("D:/Belize_DISC_2016/DATA")
setwd("C:/Users/mattf/Desktop/RSMAS/ecological_context")

# load relevant libraries

library("discr")
library("plyr")
library("stringr")
library("ggplot2")
library("lubridate")
library("circular")
library("mapproj")
library("maptools")
library("ggmap")
library("rgeos")
library("grid")
library("maps")
library("geosphere")
source("./sun_position.R")
source("./map_aids.R")

# set paths to discA & discB

pathA <- str_c(getwd(),"/","DISC_A")
pathB <- str_c(getwd(),"/","DISC_B")




## Assemble a statistics file from deployment folders -------------

statsAf <- disc_assemble("stats", deploy.dir = str_c(pathA,"/deployments"))
statsBf <- disc_assemble("stats", deploy.dir = str_c(pathB,"/deployments"))


# discA, for circular stats
statsAc <- subset(statsAf,rotation=="rotated")

# discA, for speed stats
statsAs <- subset(statsAf,rotation=="raw")

# merge them
statsA <- statsAc
statsA[,c("speed.n","speed.mean","speed.sd","speed.median","speed.mad")] <- 
  statsAs[,c("speed.n","speed.mean","speed.sd","speed.median","speed.mad")]

# discB, for circular stats
statsBc <- subset(statsBf,rotation=="rotated")

# discB, for speed stats
statsBs <- subset(statsBf,rotation=="raw")

# merge them
statsB <- statsBc
statsB[,c("speed.n","speed.mean","speed.sd","speed.median","speed.mad")] <- 
  statsBs[,c("speed.n","speed.mean","speed.sd","speed.median","speed.mad")]


# and combine the A and B data into one file

stats <- rbind(statsA,statsB)

# and save this as a .csv file for later use
setwd("C:/Users/Matt/Desktop/belize_to_do/belize16_pop_analysis/ecological_context")
write.csv(stats, "belize16_stats_file.csv")




## Determine which deployments are artifacts --------------

setwd("E:/Belize_DISC_2016/DATA")

pathA <- str_c(getwd(),"/","DISC_A")
pathB <- str_c(getwd(),"/","DISC_B")

statsAf <- disc_assemble("stats", deploy.dir = str_c(pathA,"/deployments"))
statsBf <- disc_assemble("stats", deploy.dir = str_c(pathB,"/deployments"))

stats <- rbind(statsAf,statsBf)

artifacts <- ddply(stats, ~deploy_id, function(x) {
  rw <- x[x$rotation=="raw",c("r")]
  rr <- x[x$rotation=="rotated",c("r")]
  art <- FALSE
  if(rw-rr > (.05*rr)){ art= TRUE }
  return(data.frame(rw,rr,art))
})


tracksAf <- disc_assemble("rotated_larvae_tracks", deploy.dir = str_c(pathA,"/deployments"))
tracksBf <- disc_assemble("rotated_larvae_tracks", deploy.dir = str_c(pathB,"/deployments"))

tracks <- rbind(tracksAf,tracksBf)

t <- tracks[!is.na(tracks$cameraHeading) & tracks$rotation=="raw",]
t$comp = conversion.circular(t$cameraHeading, units="radians", rotation="counter", zero=0, modulo="asis")

rotation = ddply(t, ~deploy_id, function(x){
  tmp = data.frame(a=x$comp, b=c(x$comp[2:nrow(x)], circular(NA)))
  ranges = apply(tmp, 1, range.circular)
  c(range=sum(ranges, na.rm=T))
}, .progress="text")

# rotation2 = ddply(t, ~deploy_id, function(x) c(range=range.circular(x$comp)))

avgN = round(mean(ddply(t, ~deploy_id, nrow)$V1))
avgR = mean(stats[stats$rotation=="raw", "r"])

dispersion = seq(0,2,by=0.1)
correspondingRs = c()
for (disp in dispersion) {
  currentRs = c()
  for (i in 1:100) {
    angles = rvonmises(avgN, circular(0), disp, control.circular=list(units="degree"))
    currentRs = c(currentRs, rho.circular(angles))
  }
  correspondingRs = c(correspondingRs, mean(currentRs))
}
disp = approx(x=correspondingRs, y=dispersion, xout=avgR)$y

rotAngles = seq(20, 290, by=5)
rotEffect = data.frame()
for (rotAngle in rotAngles) {
  cat(rotAngle)
  for (i in 1:100) {
    if (i%%10 == 0) {
      cat(".")
    }
    
    angles = rvonmises(avgN, circular(0), disp, control.circular=list(units="degree"))
    
    anglesRot = angles+seq(0, rotAngle, length.out=length(angles))
    
    diffR = rho.circular(angles) - rho.circular(anglesRot)
    
    # test of differences
    Utest = watson.two.test(angles, anglesRot)$statistic
    wwTest = watson.williams.test(list(angles, anglesRot))$p.value
    
    # test of significance after rotation
    rayTest = rayleigh.test(anglesRot)$p.value
    
    rotEffect = rbind(rotEffect, data.frame(rotAngle, diffR, Utest, wwTest, rayTest))
  }
  cat("\n")
}

rotEffect = ddply(rotEffect, ~rotAngle, summarise,
                   Utest = mean(Utest),
                   wwTest = mean(wwTest),
                   rayTest = mean(rayTest))

rotEffect$Usignif = rotEffect$Utest > 0.187
rotEffect$wwSignif = rotEffect$wwTest < 0.05
rotEffect$raySignif = rotEffect$rayTest < 0.05

limitAngle1 = min(rotEffect$rotAngle[rotEffect$Usignif])
limitAngle2 = min(rotEffect$rotAngle[!rotEffect$wwSignif])
limitAngle3 = min(rotEffect$rotAngle[!rotEffect$raySignif])

limitAngle <- min(limitAngle1,limitAngle2,limitAngle3)

artifacts$rotationRange = rotation$range
artifacts$enoughRotation = rotation$range > limitAngle


setwd("C:/Users/Matt/Desktop/belize_to_do/belize16_pop_analysis/ecological_context")
write.csv(artifacts, "belize16_artifacts_file.csv")




## Extract max speed & mean of last quartile of speeds --------------

tracksAf <- disc_assemble("rotated_larvae_tracks", deploy.dir = str_c(pathA,"/deployments"))
tracksBf <- disc_assemble("rotated_larvae_tracks", deploy.dir = str_c(pathB,"/deployments"))


tracks <- rbind(tracksAf,tracksBf)

t <- tracks[!is.na(tracks$cameraHeading) & tracks$rotation=="raw",]

t<- do.call(data.frame,lapply(t, function(x) replace(x, is.infinite(x),NA)))


speedStats <- ddply(t,~deploy_id, summarise,
                    meanSpeed <- mean(speed,na.rm=TRUE),
                    quartileMean <- mean(speed[speed >= data.frame(quantile(speed,na.rm=TRUE))[4,1]],na.rm=TRUE),
                    maxSpeed <- max(speed,na.rm=TRUE)
)

              
colnames(speedStats) <- c("deploy_id","meanSpeed","meanThirdQuartile","maxSpeed")







## Extract sensor data summaries from the deployment folders --------

setwd("E:/Belize_DISC_2016/DATA")

pathA <- str_c(getwd(),"/","DISC_A")
pathB <- str_c(getwd(),"/","DISC_B")

hoboAf <- disc_assemble("hobo", deploy.dir = str_c(pathA,"/deployments"))
hoboBf <- disc_assemble("hobo", deploy.dir = str_c(pathB,"/deployments"))

hoboA <- ddply(hoboAf, ~deploy_id, summarise,
               meanTemp = mean(temp),
               maxTemp = max(temp),
               minTemp = min(temp),
               meanLight = mean(light),
               maxLight = max(light),
               minLight = min(light))

hoboB <- ddply(hoboBf, ~deploy_id, summarise,
               meanTemp = mean(temp),
               maxTemp = max(temp),
               minTemp = min(temp),
               meanLight = mean(light),
               maxLight = max(light),
               minLight = min(light))

hobo <- rbind(hoboA,hoboB)

setwd("C:/Users/Matt/Desktop/belize_to_do/belize16_pop_analysis/ecological_context")

write.csv(hobo, "belize16_hobo_file.csv")





## Load and Clean all the Data Files ---------------

setwd("C:/Users/mattf/Desktop/RSMAS/ecological_context")

# read log and stats files, join them
log <- read.csv("deployment_log2.csv")
log <- within(log, rm("X","X.1"))

stats <- read.csv("belize16_stats_file.csv")
stats <- within(stats, rm("X"))

d <- join(log,stats,by="deploy_id",type="right") # keep only deployments with stats data

# read the hobo data in
hobo <- read.csv("belize16_hobo_file.csv")

d <- join(d,hobo,by="deploy_id",type="left") # keep all deployments, match hobo where we have it

# correct a couple rows with celsius temperatures

d[!is.na(d$meanTemp) & d$meanTemp < 50,c("meanTemp")] <- d[!is.na(d$meanTemp) & d$meanTemp < 50,c("meanTemp")] * (9/5) + 32

# read the artifacts data
artifacts <- read.csv("belize16_artifacts_file.csv")

d <- join(d,artifacts,by="deploy_id",type="left") # keep all deployments, match artifacts where we have it



# change any data types if necessary, e.g.
str(d)
d$deploy_id <- as.factor(d$deploy_id)
d$depth..m. <- as.factor(d$depth..m.)
d$fish_age..dph. <- as.numeric(d$fish_age..dph.)
d$bin_id <- as.factor(d$bin_id)
d$batch_id <- as.factor(d$batch_id)
d$wind_dir <- as.circular(as.numeric(as.character(d$wind_dir)), units="degrees")
d$drift_direction <- as.circular(as.numeric(as.character(d$drift_direction)), units="degrees")
d$drift_distance <- as.numeric(d$drift_distance)
d$In_Behavior <- as.factor(d$In_Behavior)
d$mean <- as.circular(d$mean, units="degrees")
d$se.mean <- as.circular(d$se.mean, units="degrees")
d$dir.mean <- as.circular(d$dir.mean, units="degrees")
d$dir.se.mean <- as.circular(d$dir.se.mean, units="degrees")
d$turn.abs.mean <- as.circular(d$turn.abs.mean, units="degrees")
d$moonAlt <- as.circular(as.numeric(as.character(d$moonAlt)), units="degrees")
d$moonAz <- as.circular(as.numeric(as.character(d$moonAz)), units="degrees")
d$moonHor <- as.logical(d$moonHor)

# add a discrete dateTime variable
# see ?strptime
timeIn <- as.POSIXct(str_c(d$date_start," ", d$time_start,sep=""),format="%m/%d/%Y %H:%M:%S")
d$dateTime <- timeIn + minutes(10) # we asssume deployments are ~20 minutes and we use the midpoint of the deployment

# fix the gps variables

lat1 <- str_split_fixed(as.character(d$lat_start),pattern="'",n=2)
d$lat_start <- (as.numeric(lat1[,1]) + as.numeric(lat1[,2])/60)

lat2 <- str_split_fixed(as.character(d$lat_stop),pattern="'",n=2)
d$lat_stop <- (as.numeric(lat2[,1]) + as.numeric(lat2[,2])/60)

lon1 <- str_split_fixed(as.character(d$lon_start),pattern="'",n=2)
d$lon_start <- -(as.numeric(lon1[,1]) + as.numeric(lon1[,2])/60)

lon2 <- str_split_fixed(as.character(d$lon_stop),pattern="'",n=2)
d$lon_stop <- -(as.numeric(lon2[,1]) + as.numeric(lon2[,2])/60)
  





## Plot a map of the deployment locations and transect -----------------------

# convert total drift to cm/s
d$current_strength <- d$drift_distance * 100 / 25 / 60

# get rid of inf values
d <- do.call(data.frame,lapply(d, function(x) replace(x, is.infinite(x),NA)))


# get google sat of the region
test <- get_map(location=c(lon=-88.0766,lat=16.82),zoom=14,source="google",maptype="satellite")


# compute mean location, rreverse signs for lat/lon if needed)
d$lat <- apply(d[,c("lat_start","lat_stop")], 1, mean, na.rm=T)
d$lon <- apply(d[,c("lon_start","lon_stop")], 1, mean, na.rm=T)

# load transect waypoints
wp <- read.csv("transect_waypoints.csv")


d2 <- subset(d, fish_age..dph. >17 & fish_age..dph. <25)


#### PUBLICATION PLOT - MAP
# plot location of all ecology-context exp points
ggmap(test) + geom_point(aes(x=lon_start,y=lat_start,color=signif),size=2,data=d2) + # plot basemap and deploys
  scale_color_manual(values=c("hotpink","white")) +
  geom_point(aes(x=lon,y=lat),data=wp,size=2,color="yellow") +
  geom_path(aes(x=lon,y=lat),data=wp,size=1,color="yellow") +
  scale_x_continuous(limits = c(-88.085, -88.0655),breaks = c(-88.084,-88.0785,-88.072,-88.0655), expand = c(0, 0)) +
  scale_y_continuous(limits = c(16.8055, 16.837),breaks = c(16.806,16.816,16.8265,16.837), expand = c(0, 0)) +
  xlab("Longitude") + ylab("Latitude") +
  theme(axis.text=element_text(size=12,vjust=.5)) +
  theme(axis.title = element_text(size=18,face="bold",vjust=.8)) + 
  theme(legend.position="none") #+
  # scaleBar(lon=-88.076,lat=16.807,distanceLon=.25,distanceLat=0.03,distanceLegend=.07,dist.unit="km",legend.colour="white",
  #          orientation=TRUE,arrow.length=0.12,arrow.distance=0.1,arrow.North.size=5)
  # 


## Plot a map of current / orientation / net vector for DNF only --------------

# just flood
dnf <- subset(d2, tide=="flood" & location=="near" & depth..m.=="18")


# computing the net vector
dnf$Cy <- as.numeric(dnf$current_strength * sin(dnf$drift_direction*pi/180))
dnf$Cx <- as.numeric(dnf$current_strength * cos(dnf$drift_direction*pi/180))

dnf$Sy <- as.numeric(2 * sin(dnf$mean*pi/180))
dnf$Sx <- as.numeric(2 * cos(dnf$mean*pi/180))

dnf$Ty <- dnf$Sy + dnf$Cy # y component of displacement vector
dnf$Tx <- dnf$Sx + dnf$Cx # x component of displacement vector

dnf$Tr <- sqrt(dnf$Ty^2 + dnf$Tx^2) # length of displacement vector
dnf$Tt <- atan2(dnf$Ty,dnf$Tx)*180/pi # angle of displacement vector (in degrees)

if (dnf$Tt <0) { dnf$Tt <- 360 + dnf$Tt}

# use this to scale lines on the plot

# note: the swim speeds and current strengths are in cm/s, but we need them
# to be in lat/lon per 20 minutes... sooo...
# cm/s * 60 s/min * 25 minutes * 1 m/100cm = meters
# then, use estimate that 111,111 m in y is 1 degree of latitude
# and then 111,111*cos(lat) is 1 degree of longitude, where cos(lat)=0.957
# but this 0.957 ends up being trivial, so its:
#  *60*25/100/111,111 = *0.001
# 

# note - these Cy and Cx vectors are reversed bc of R angle conventions
dnf$xe2 <- 0.001*dnf$Cy + dnf$lon
dnf$ye2 <- 0.001*dnf$Cx + dnf$lat

dnf$ax2 <- 0.001*dnf$Sy + dnf$lon
dnf$ay2 <- 0.001*dnf$Sx + dnf$lat

dnf$ttx <- 0.001*dnf$Ty + dnf$lon
dnf$tty <- 0.001*dnf$Tx + dnf$lat


ggmap(test) + geom_point(aes(x=lon,y=lat),size=3,data=dnf) + # plot basemap and deploys
  facet_wrap(~tide+location+depth..m.) + # split by the treatments
  geom_segment(aes(x=lon,y=lat,xend=xe2,yend=ye2),alpha=.7,color="green",size=1,arrow = arrow(length = unit(0.1,"cm")), data=dnf) +
  geom_segment(aes(x=lon,y=lat,xend=ax2,yend=ay2),alpha=.7,color="brown",size=1,arrow = arrow(length = unit(0.1,"cm")), data=dnf) +
  geom_segment(aes(x=lon,y=lat,xend=ttx,yend=tty),alpha=.7,color="orange",size=1,arrow = arrow(length = unit(0.1,"cm")), data=dnf) +
  theme(strip.background = element_blank(),strip.text.x = element_blank()) +
  xlab("Longitude") + ylab("Latitude") +
  scale_x_continuous(limits = c(-88.09, -88.06), expand = c(0, 0)) +
  scale_y_continuous(limits = c(16.8, 16.845), expand = c(0, 0))



## Stats for the length of current vectors versus net vectors for all signif data --------------

d2s <- subset(d2, signif=="TRUE")
d2s <- subset(d2s, location=="near" & depth..m.==18)

# computing the net vector
d2s$Cy <- as.numeric(d2s$current_strength * sin(d2s$drift_direction*pi/180))
d2s$Cx <- as.numeric(d2s$current_strength * cos(d2s$drift_direction*pi/180))

d2s$Sy <- as.numeric(2 * sin(d2s$mean*pi/180))
d2s$Sx <- as.numeric(2 * cos(d2s$mean*pi/180))

d2s$Ty <- d2s$Sy + d2s$Cy # y component of displacement vector
d2s$Tx <- d2s$Sx + d2s$Cx # x component of displacement vector

d2s$Tr <- sqrt(d2s$Ty^2 + d2s$Tx^2) # length of displacement vector
d2s$Tt <- atan2(d2s$Ty,d2s$Tx)*180/pi # angle of displacement vector (in degrees)

# shapiro test reveals that current_strength and Tr are not normally distributed
# use wilcoxon signed rank test (non-parametric for paired data) instead

d2s <- subset(d2s, current_strength >= 1)

wilcox.test(d2s$current_strength, d2s$Tr,paired=TRUE)





## Proportion of orienting individuals ----------------


# what are we working with?
summaryStats <- ddply(d,.(depth..m.,tide,location),function(x){
  N <-nrow(x)
})

colnames(summaryStats) <- c("Depth (m)", "Tide", "Location", "N")
summaryStats <- summaryStats[c(3,1,2,4)]

propTable <- count(d$signif)
artTable <- count(d$art)

twoTable <- count(d, vars=c("signif","art"))


# let artifacts remove significancy
d[d$art==TRUE,c("signif")] <- FALSE


# And split by environmental context
d2 <- subset(d, fish_age..dph. >17 & fish_age..dph. <25)

propStatsEco <- ddply(d2,.(depth..m.,tide,location),function(x){
  propTableEco <- count(x,vars=c("signif"))
})

logitFit <- glm(signif ~ tide+depth..m.+location, data=d2, family="binomial")
summary(logitFit)


# try a big logit model 


logitFitFull <- glm(signif ~ tide+depth..m.+location+fish_age..dph.+sky+dateTime, 
                    data=d, family="binomial")
summary(logitFitFull)


# proportion orientating by ontogeny

dneT <- subset(d, depth..m.=="18" & location=="near" & tide=="ebb")

logitFitOnt <- glm(signif ~ fish_age..dph., 
                    data=dneT, family="binomial")
summary(logitFitOnt)
# significant finding - they don't like rainy days


# something going on with the rain...







## Basic population-level orientation analysis ---------------

# remove not directional individuals
d <- d[d$signif == "TRUE",]

 # define the basic function
circStats <- function(x) { # define descriptive statistics function
  n = nrow(x)
  bearing = mean.circular(x$mean)
  sdbearing = sd.circular(x$mean)
  meanR = mean(x$r, na.rm=T)
  sdR = sd(x$r, na.rm=T)
  ray = rayleigh.test(na.omit(x$mean))
  r = ray$statistic
  p = ray$p.value
  signif = p < 0.05
  return(data.frame(n, bearing,sdbearing, meanR, sdR, r, p, signif))
}

# and draw a circle to use for plots
circleFun <- function(center = c(0,0),diameter = 6, npoints = 100){  
  r = diameter / 2                                                 # note the radius of the circle = 3
  tt <- seq(0,2*pi,length.out = npoints)
  xx <- center[1] + r * cos(tt)
  yy <- center[2] + r * sin(tt)
  return(data.frame(x = xx, y = yy))
}
# now we use the function
circle = circleFun() 
circle = car2pol(circle) # converting to polar coordinates
circle$theta = (circle$theta * 180)/pi # converting to degrees
circle[100,c("theta")] = 359.99999 # and changing the last point so the circle draws completely
circle2 = circleFun(diameter=20)
circle2 = car2pol(circle2)
circle2$theta = (circle2$theta * 180)/pi
circle2[100,c("theta")] = 359.99999 




# bin angles for plot
bin <- 5 # binning parameter 
d$theta <- as.numeric(round_any(d$mean, bin))
d$theta[d$theta==360] <- 0
# create a data.frame with count
counts <- count(d, vars="theta")
  
# repeat each point the appropriate number of times
db <- adply(counts, 1, function(x) {
  data.frame(theta=x$theta, count=1:x$freq)
}, .expand=F)
  
# make the scale prettier
db$count <- 10 + db$count


circStatsAll <- circStats(d) # run the statistics


p <- ggplot(db) + polar() + # set the data and the polar layout
  # draw circles on the plot for the r-axis
  geom_line(data=circle, aes(theta,rho), color="gray") +
  geom_line(data=circle2, aes(theta,rho), color="gray") +
  # add the data
  geom_point(aes(x=theta, y=count), size=1) + # add the bearings as stacked dots
  geom_segment(aes(x=bearing, y=0, xend=bearing, yend=r*10, color=signif), size = .6, data=circStatsAll) + # add the mean lines, w/ type determined by significance
  scale_color_manual(values=c("gray", "black")) + # set black line to be significant
  # groom the plot a bit
  scale_x_continuous("", limits=c(0,360), breaks=seq(0,360-1,by=45), labels=c("N","","","","S","","","")) + # set the X axis (angle axis)
  theme_bw() + labs(color="Statistically Significant") + # remove the gray background and adjust legend label
  theme(legend.title = element_text(size=8)) + theme(legend.text = element_text(size=6)) + # adjust legend text sizes
  theme(panel.grid.minor = element_line(size = .2, color = "white")) + theme(panel.grid.major = element_line(size=.5, color = "white")) + # remove x-axis grid, make y-axis grid black
  theme(plot.title = element_text(face="bold", size=18)) + theme(axis.title.y = element_text(size=8, face="bold")) + # format the plot title and r-axis label
  theme(axis.text.x = element_text(size=14, vjust=.05,face="bold")) + theme(axis.text.y = element_blank(), axis.ticks = element_blank()) + # adjust axis text
  theme(strip.text.x = element_text(size=6, face="bold"), strip.background = element_rect(colour="black", fill="#CCCCFF"))  # format facet labels

p




## Population-level orientation by environmental context --------------------

d2 <- subset(d, fish_age..dph. >17 & fish_age..dph. <25)
d2$expcode <- 0
d2[d2$depth..m.=="18" & d2$location=="near" & d2$tide=="ebb",c("expcode")] <- 1
d2[d2$depth..m.=="9" & d2$location=="near" & d2$tide=="ebb",c("expcode")] <- 2
d2[d2$depth..m.=="18" & d2$location=="far" & d2$tide=="ebb",c("expcode")] <- 3
d2[d2$depth..m.=="18" & d2$location=="near" & d2$tide=="flood",c("expcode")] <- 4
   
   

db2 <- ddply(d2,.(expcode), function(x, bin) {  # 'splitting' by method here b/c theres only 1 method and its easier than changing code to remove ddply
  x$mean <- as.numeric(round_any(x$mean, bin))
  x$mean[x$mean==360] <- 0
  
  # create a data.frame with count
  counts <- count(x, "mean")
  
  # repeat each point the appropriate number of times
  d <- adply(counts, 1, function(x) {
    data.frame(theta=x$mean, count=1:x$freq)
  }, .expand=F)
  
  # make the scale prettier
  d$count <- 10 + d$count
  
  return(d)
}, bin=bin)

circStatsEco <- ddply(d2,.(expcode),circStats) # run the statistics
circStatsEco[circStatsEco$bearing<0,c("bearing")] <- circStatsEco[circStatsEco$bearing<0,c("bearing")] + 360

# PUBLICATION PLOT ---- FIGURE 2 MULTIPLANEL CARDINAL ORIENTATION
p2 <- ggplot(db2) + polar() + # set the data and the polar layout
  # draw circles on the plot for the r-axis
  geom_line(data=circle, aes(theta,rho),size=1, color="gray") +
  geom_line(data=circle2, aes(theta,rho),size=1, color="gray") +
  # add the data
  geom_point(aes(x=theta, y=count), size=2) + # add the bearings as stacked dots
  geom_segment(aes(x=bearing, y=0, xend=bearing, yend=r*10),color="black", size = 1, data=circStatsEco) + # add the mean lines, w/ type determined by significance
  facet_wrap(~expcode) +
  # groom the plot a bit
  scale_x_continuous("", limits=c(0,360), breaks=seq(0,360-1,by=45), labels=c("N","","","","S","","","")) + # set the X axis (angle axis)
  theme_bw() + labs(color="Statistically Significant") + # remove the gray background and adjust legend label
  theme(legend.title = element_text(size=8)) + theme(legend.text = element_text(size=6)) + # adjust legend text sizes
  theme(panel.grid.minor = element_line(size = .2, color = "white")) + theme(panel.grid.major = element_line(size=.5, color = "white")) + # remove x-axis grid, make y-axis grid black
  theme(plot.title = element_text(face="bold", size=18)) + theme(axis.title.y = element_text(size=16, face="bold")) + # format the plot title and r-axis label
  theme(axis.text.x = element_text(size=14, face="bold", vjust=.05)) + theme(axis.text.y = element_text(size=16)) + # adjust axis text
  scale_y_continuous(breaks=c(0,3,6,10),labels=c(0,0.3,0.6,1)) +
  # theme(strip.text.x = element_text(size=6, face="bold"), strip.background = element_rect(colour="black", fill="#CCCCFF"))  # format facet labels
  theme(strip.background = element_blank(), strip.text.x = element_blank()) # or remove them

p2




## PUBLICATION PLOT - Shallow,Near,Ebb Orientation

sne <- subset(d2,expcode==2)

dbsne <- ddply(sne,.(expcode), function(x, bin) {  # 'splitting' by method here b/c theres only 1 method and its easier than changing code to remove ddply
  x$mean <- as.numeric(round_any(x$mean, bin))
  x$mean[x$mean==360] <- 0
  
  # create a data.frame with count
  counts <- count(x, "mean")
  
  # repeat each point the appropriate number of times
  d <- adply(counts, 1, function(x) {
    data.frame(theta=x$mean, count=1:x$freq)
  }, .expand=F)
  
  # make the scale prettier
  d$count <- 10 + d$count
  
  return(d)
}, bin=bin)

circStatssne <- circStats(sne) # run the statistics
circStatssne[circStatssne$bearing < 0,c("bearing")] <- 360 + circStatssne[circStatssne$bearing < 0,c("bearing")]

psne <- ggplot(dbsne) + polar() + # set the data and the polar layout
  # draw circles on the plot for the r-axis
  geom_line(data=circle, aes(theta,rho), color="gray",size=1.5) +
  geom_line(data=circle2, aes(theta,rho), color="gray",size=1.5) +
  # add the data
  geom_point(aes(x=theta, y=count), size=4) + 
  geom_segment(aes(x=bearing, y=0, xend=bearing, yend=r*10), size = 1, data=circStatssne) + 
  scale_color_manual(values=c("gray", "black")) + 
  # groom the plot
  theme_bw() + 
  theme(panel.grid.minor = element_line(size = .2, color = "white")) + 
  theme(panel.grid.major = element_line(size=.5, color = "white")) +
  theme(axis.title.y = element_text(size=18, face="bold")) + 
  theme(axis.text.y = element_text(size=16)) +
  scale_y_continuous(breaks=c(0,3,6,10),labels=c(0,0.3,0.6,1)) +
  theme(axis.text.x = element_text(face="bold",size=22, vjust=.05)) + 
  scale_x_continuous("", limits=c(0,360), breaks=seq(0,360-1,by=45), labels=c("N","","","","S","","","")) 

psne


## PUBLICATION PLOT - Deep,Far,Ebb Orientaiton

dfe <- subset(d2,expcode==3)

dbdfe <- ddply(dfe,.(expcode), function(x, bin) {  # 'splitting' by method here b/c theres only 1 method and its easier than changing code to remove ddply
  x$mean <- as.numeric(round_any(x$mean, bin))
  x$mean[x$mean==360] <- 0
  
  # create a data.frame with count
  counts <- count(x, "mean")
  
  # repeat each point the appropriate number of times
  d <- adply(counts, 1, function(x) {
    data.frame(theta=x$mean, count=1:x$freq)
  }, .expand=F)
  
  # make the scale prettier
  d$count <- 10 + d$count
  
  return(d)
}, bin=bin)

circStatsdfe <- circStats(dfe) # run the statistics
circStatsdfe[circStatsdfe$bearing < 0,c("bearing")] <- 360 + circStatsdfe[circStatsdfe$bearing < 0,c("bearing")]

pdfe <- ggplot(dbdfe) + polar() + # set the data and the polar layout
  # draw circles on the plot for the r-axis
  geom_line(data=circle, aes(theta,rho), color="gray",size=1.5) +
  geom_line(data=circle2, aes(theta,rho), color="gray",size=1.5) +
  # add the data
  geom_point(aes(x=theta, y=count), size=4) + 
  geom_segment(aes(x=bearing, y=0, xend=bearing, yend=r*10), size = 1, data=circStatsdfe) + 
  scale_color_manual(values=c("gray", "black")) + 
  # groom the plot
  theme_bw() + 
  theme(panel.grid.minor = element_line(size = .2, color = "white")) + 
  theme(panel.grid.major = element_line(size=.5, color = "white")) +
  theme(axis.title.y = element_text(size=18, face="bold")) + 
  theme(axis.text.y = element_text(size=16)) +
  scale_y_continuous(breaks=c(0,3,6,10),labels=c(0,0.3,0.6,1)) +
  theme(axis.text.x = element_text(face="bold",size=22, vjust=.05)) + 
  scale_x_continuous("", limits=c(0,360), breaks=seq(0,360-1,by=45), labels=c("N","","","","S","","","")) 

pdfe


## PUBLICATOIN PLOT - DFE versus Reference Turning Angles

dfedne <- subset(d2,expcode==3 | expcode==1)

dbdfedne <- ddply(dfedne,.(expcode), function(x, bin) {  # 'splitting' by method here b/c theres only 1 method and its easier than changing code to remove ddply
  x$turn <- as.numeric(round_any(x$turn.abs.mean, bin))
  x$turn[x$turn==360] <- 0
  
  # create a data.frame with count
  counts <- count(x, "turn")
  
  # repeat each point the appropriate number of times
  d <- adply(counts, 1, function(x) {
    data.frame(theta=x$turn, count=1:x$freq)
  }, .expand=F)
  
  # make the scale prettier
  d$count <- 10 + d$count
  
  return(d)
}, bin=bin)

turnStatsdfedne <- ddply(dfedne,~expcode,circStatsTurn) # run the statistics

pdfedne <- ggplot(dbdfedne) + polar() + # set the data and the polar layout
  # draw circles on the plot for the r-axis
  geom_line(data=circle, aes(theta,rho), color="gray",size=1.5) +
  geom_line(data=circle2, aes(theta,rho), color="gray",size=1.5) +
  # add the data
  geom_point(aes(x=theta, y=count), size=2) + 
  geom_segment(aes(x=meanTurn, y=0, xend=meanTurn, yend=r*10), size = 1, data=turnStatsdfedne) + 
  scale_color_manual(values=c("gray", "black")) + 
  facet_wrap(~expcode) +
  # add a comparison line
  geom_segment(aes(x=55.75840,y=0,xend=55.75840,yend=10),size=1,linetype=2) +
  # groom the plot
  theme_bw() + 
  theme(strip.background = element_blank(), strip.text.x = element_blank()) +
  theme(panel.grid.minor = element_line(size = .2, color = "white")) + 
  theme(panel.grid.major = element_line(size=.5, color = "white")) +
  theme(axis.title.y = element_text(size=18, face="bold")) + 
  theme(axis.text.y = element_text(size=16)) +
  scale_y_continuous(breaks=c(0,3,6,10),labels=c(0,0.3,0.6,1)) +
  theme(axis.text.x = element_text(face="bold",size=18, vjust=.05)) + 
  scale_x_continuous("", limits=c(0,360), breaks=seq(0,360-1,by=45), labels=c("0° Turn","","","","180° Turn","","","")) 

pdfedne



## PUBLICATION PLOT - DNE versus DNF Orientation

dnednf <- subset(d2,expcode==1 | expcode==4)

dbdnednf <- ddply(dnednf,.(expcode), function(x, bin) {  # 'splitting' by method here b/c theres only 1 method and its easier than changing code to remove ddply
  x$mean <- as.numeric(round_any(x$mean, bin))
  x$mean[x$mean==360] <- 0
  
  # create a data.frame with count
  counts <- count(x, "mean")
  
  # repeat each point the appropriate number of times
  d <- adply(counts, 1, function(x) {
    data.frame(theta=x$mean, count=1:x$freq)
  }, .expand=F)
  
  # make the scale prettier
  d$count <- 10 + d$count
  
  return(d)
}, bin=bin)

circStatsdnednf <- ddply(dnednf,~expcode,circStats)

pdnednf <- ggplot(dbdnednf) + polar() + # set the data and the polar layout
  # draw circles on the plot for the r-axis
  geom_line(data=circle, aes(theta,rho), color="gray",size=1.5) +
  geom_line(data=circle2, aes(theta,rho), color="gray",size=1.5) +
  # add the data
  geom_point(aes(x=theta, y=count), size=3) + 
  geom_segment(aes(x=bearing, y=0, xend=bearing, yend=r*10), size = 1, data=circStatsdnednf) + 
  scale_color_manual(values=c("gray", "black")) + 
  facet_wrap(~expcode) +
 # groom the plot
  theme_bw() + 
  theme(strip.background = element_blank(), strip.text.x = element_blank()) +
  theme(panel.grid.minor = element_line(size = .2, color = "white")) + 
  theme(panel.grid.major = element_line(size=.5, color = "white")) +
  theme(axis.title.y = element_text(size=18, face="bold")) + 
  theme(axis.text.y = element_text(size=16)) +
  scale_y_continuous(breaks=c(0,3,6,10),labels=c(0,0.3,0.6,1)) +
  theme(axis.text.x = element_text(face="bold",size=22, vjust=.05)) + 
  scale_x_continuous("", limits=c(0,360), breaks=seq(0,360-1,by=45), labels=c("N","","","","S","","","")) 

pdnednf







## Analysis for each condition independently --------------

dne <- d2[d2$depth..m.=="18" & d2$location=="near" & d2$tide=="ebb",]

dbdne <- ddply(dne,.(expcode), function(x, bin) {  # 'splitting' by method here b/c theres only 1 method and its easier than changing code to remove ddply
  x$mean <- as.numeric(round_any(x$mean, bin))
  x$mean[x$mean==360] <- 0
  
  # create a data.frame with count
  counts <- count(x, "mean")
  
  # repeat each point the appropriate number of times
  d <- adply(counts, 1, function(x) {
    data.frame(theta=x$mean, count=1:x$freq)
  }, .expand=F)
  
  # make the scale prettier
  d$count <- 10 + d$count
  
  return(d)
}, bin=bin)

circStatsdne <- circStats(dne) # run the statistics


pdne <- ggplot(dbdne) + polar() + # set the data and the polar layout
  # draw circles on the plot for the r-axis
  geom_line(data=circle, aes(theta,rho), color="gray") +
  geom_line(data=circle2, aes(theta,rho), color="gray") +
  # add the data
  geom_point(aes(x=theta, y=count), size=1) + # add the bearings as stacked dots
  geom_segment(aes(x=bearing, y=0, xend=bearing, yend=r*10, color=signif), size = .6, data=circStatsdne) + # add the mean lines, w/ type determined by significance
  scale_color_manual(values=c("gray", "black")) + # set black line to be signif
  facet_wrap(~expcode) +
  # groom the plot a bit
  scale_x_continuous("", limits=c(0,360), breaks=seq(0,360-1,by=45), labels=c("N","","","","S","","","")) + # set the X axis (angle axis)
  theme_bw() + labs(color="Statistically Significant") + # remove the gray background and adjust legend label
  theme(legend.title = element_text(size=8)) + theme(legend.text = element_text(size=6)) + # adjust legend text sizes
  theme(panel.grid.minor = element_line(size = .2, color = "white")) + theme(panel.grid.major = element_line(size=.5, color = "white")) + # remove x-axis grid, make y-axis grid black
  theme(plot.title = element_text(face="bold", size=18)) + theme(axis.title.y = element_text(size=8, face="bold")) + # format the plot title and r-axis label
  theme(axis.text.x = element_text(size=10, vjust=.05,face="bold")) + theme(axis.text.y = element_blank(), axis.ticks = element_blank()) + # adjust axis text
  # theme(strip.text.x = element_text(size=6, face="bold"), strip.background = element_rect(colour="black", fill="#CCCCFF"))  # format facet labels
  theme(strip.background = element_blank(), strip.text.x = element_blank()) # or remove them

pdne


sne <- d2[d2$depth..m.=="9" & d2$location=="near" & d2$tide=="ebb",]

dbsne <- ddply(sne,.(expcode), function(x, bin) {  # 'splitting' by method here b/c theres only 1 method and its easier than changing code to remove ddply
  x$mean <- as.numeric(round_any(x$mean, bin))
  x$mean[x$mean==360] <- 0
  
  # create a data.frame with count
  counts <- count(x, "mean")
  
  # repeat each point the appropriate number of times
  d <- adply(counts, 1, function(x) {
    data.frame(theta=x$mean, count=1:x$freq)
  }, .expand=F)
  
  # make the scale prettier
  d$count <- 10 + d$count
  
  return(d)
}, bin=bin)

circStatssne <- circStats(sne) # run the statistics


psne <- ggplot(dbsne) + polar() + # set the data and the polar layout
  # draw circles on the plot for the r-axis
  geom_line(data=circle, aes(theta,rho), color="gray") +
  geom_line(data=circle2, aes(theta,rho), color="gray") +
  # add the data
  geom_point(aes(x=theta, y=count), size=1) + # add the bearings as stacked dots
  geom_segment(aes(x=bearing, y=0, xend=bearing, yend=r*10, color=signif), size = .6, data=circStatssne) + # add the mean lines, w/ type determined by significance
  scale_color_manual(values=c("gray", "black")) + # set black line to be signif
  facet_wrap(~expcode) +
  # groom the plot a bit
  scale_x_continuous("", limits=c(0,360), breaks=seq(0,360-1,by=45), labels=c("N","","","","S","","","")) + # set the X axis (angle axis)
  theme_bw() + labs(color="Statistically Significant") + # remove the gray background and adjust legend label
  theme(legend.title = element_text(size=8)) + theme(legend.text = element_text(size=6)) + # adjust legend text sizes
  theme(panel.grid.minor = element_line(size = .2, color = "white")) + theme(panel.grid.major = element_line(size=.5, color = "white")) + # remove x-axis grid, make y-axis grid black
  theme(plot.title = element_text(face="bold", size=18)) + theme(axis.title.y = element_text(size=8, face="bold")) + # format the plot title and r-axis label
  theme(axis.text.x = element_text(size=10, vjust=.05,face="bold")) + theme(axis.text.y = element_blank(), axis.ticks = element_blank()) + # adjust axis text
  # theme(strip.text.x = element_text(size=6, face="bold"), strip.background = element_rect(colour="black", fill="#CCCCFF"))  # format facet labels
  theme(strip.background = element_blank(), strip.text.x = element_blank()) # or remove them

psne


dfe <- d2[d2$depth..m.=="18" & d2$location=="far" & d2$tide=="ebb",]

dbdfe <- ddply(dfe,.(expcode), function(x, bin) {  # 'splitting' by method here b/c theres only 1 method and its easier than changing code to remove ddply
  x$mean <- as.numeric(round_any(x$mean, bin))
  x$mean[x$mean==360] <- 0
  
  # create a data.frame with count
  counts <- count(x, "mean")
  
  # repeat each point the appropriate number of times
  d <- adply(counts, 1, function(x) {
    data.frame(theta=x$mean, count=1:x$freq)
  }, .expand=F)
  
  # make the scale prettier
  d$count <- 10 + d$count
  
  return(d)
}, bin=bin)

circStatsdfe <- circStats(dfe) # run the statistics


pdfe <- ggplot(dbdfe) + polar() + # set the data and the polar layout
  # draw circles on the plot for the r-axis
  geom_line(data=circle, aes(theta,rho), color="gray") +
  geom_line(data=circle2, aes(theta,rho), color="gray") +
  # add the data
  geom_point(aes(x=theta, y=count), size=1) + # add the bearings as stacked dots
  geom_segment(aes(x=bearing, y=0, xend=bearing, yend=r*10, color=signif), size = .6, data=circStatsdfe) + # add the mean lines, w/ type determined by significance
  scale_color_manual(values=c("gray", "black")) + # set black line to be signif
  facet_wrap(~expcode) +
  # groom the plot a bit
  scale_x_continuous("", limits=c(0,360), breaks=seq(0,360-1,by=45), labels=c("N","","","","S","","","")) + # set the X axis (angle axis)
  theme_bw() + labs(color="Statistically Significant") + # remove the gray background and adjust legend label
  theme(legend.title = element_text(size=8)) + theme(legend.text = element_text(size=6)) + # adjust legend text sizes
  theme(panel.grid.minor = element_line(size = .2, color = "white")) + theme(panel.grid.major = element_line(size=.5, color = "white")) + # remove x-axis grid, make y-axis grid black
  theme(plot.title = element_text(face="bold", size=18)) + theme(axis.title.y = element_text(size=8, face="bold")) + # format the plot title and r-axis label
  theme(axis.text.x = element_text(size=10, vjust=.05,face="bold")) + theme(axis.text.y = element_blank(), axis.ticks = element_blank()) + # adjust axis text
  # theme(strip.text.x = element_text(size=6, face="bold"), strip.background = element_rect(colour="black", fill="#CCCCFF"))  # format facet labels
  theme(strip.background = element_blank(), strip.text.x = element_blank()) # or remove them

pdfe



dnf <- d2[d2$depth..m.=="18" & d2$location=="near" & d2$tide=="flood",]

dbdnf <- ddply(dnf,.(expcode), function(x, bin) {  # 'splitting' by method here b/c theres only 1 method and its easier than changing code to remove ddply
  x$mean <- as.numeric(round_any(x$mean, bin))
  x$mean[x$mean==360] <- 0
  
  # create a data.frame with count
  counts <- count(x, "mean")
  
  # repeat each point the appropriate number of times
  d <- adply(counts, 1, function(x) {
    data.frame(theta=x$mean, count=1:x$freq)
  }, .expand=F)
  
  # make the scale prettier
  d$count <- 10 + d$count
  
  return(d)
}, bin=bin)

circStatsdnf <- circStats(dnf) # run the statistics


pdnf <- ggplot(dbdnf) + polar() + # set the data and the polar layout
  # draw circles on the plot for the r-axis
  geom_line(data=circle, aes(theta,rho), color="gray") +
  geom_line(data=circle2, aes(theta,rho), color="gray") +
  # add the data
  geom_point(aes(x=theta, y=count), size=1) + # add the bearings as stacked dots
  geom_segment(aes(x=bearing, y=0, xend=bearing, yend=r*10, color=signif), size = .6, data=circStatsdnf) + # add the mean lines, w/ type determined by significance
  scale_color_manual(values=c("gray", "black")) + # set black line to be signif
  facet_wrap(~expcode) +
  # groom the plot a bit
  scale_x_continuous("", limits=c(0,360), breaks=seq(0,360-1,by=45), labels=c("N","","","","S","","","")) + # set the X axis (angle axis)
  theme_bw() + labs(color="Statistically Significant") + # remove the gray background and adjust legend label
  theme(legend.title = element_text(size=8)) + theme(legend.text = element_text(size=6)) + # adjust legend text sizes
  theme(panel.grid.minor = element_line(size = .2, color = "white")) + theme(panel.grid.major = element_line(size=.5, color = "white")) + # remove x-axis grid, make y-axis grid black
  theme(plot.title = element_text(face="bold", size=18)) + theme(axis.title.y = element_text(size=8, face="bold")) + # format the plot title and r-axis label
  theme(axis.text.x = element_text(size=10, vjust=.05,face="bold")) + theme(axis.text.y = element_blank(), axis.ticks = element_blank()) + # adjust axis text
  # theme(strip.text.x = element_text(size=6, face="bold"), strip.background = element_rect(colour="black", fill="#CCCCFF"))  # format facet labels
  theme(strip.background = element_blank(), strip.text.x = element_blank()) # or remove them

pdnf





## Population-level orientation with age ---------------

cla <- lm.circular(y=d$mean,x=d$fish_age..dph.,type="c-l",init=0)
cla$p.values
# close to signif


# only reference group
dneT <- subset(d, depth..m.=="18" & location=="near" & tide=="ebb")
cladneT <- lm.circular(y=dneT$mean,x=dneT$fish_age..dph.,type="c-l",init=0)
# terrible





## Orientation by specific age groups -------------

circStatsAge <- ddply(d,.(fish_age..dph.),circStats) # run the statistics
circStatsAge$age <- circStatsAge$fish_age..dph.
for (i in 1:nrow(circStatsAge)){  
  if (circStatsAge[i,c("bearing")]<=0){
    circStatsAge[i,c("bearing")] <- 360 + circStatsAge[i,c("bearing")]
  }
}

# bin angles for plot
dneT <- subset(d, depth..m.=="18" & location=="near" & tide=="ebb")
t <- dneT
t$theta <- as.numeric(round_any(d$mean, bin))
t$theta[t$theta==360] <- 0
# create a data.frame with count
counts <- count(t, vars=c("theta","fish_age..dph."))

# repeat each point the appropriate number of times
tb <- adply(counts, 1, function(x) {
  data.frame(theta=x$theta,age=x$fish_age..dph., count=1:x$freq)
}, .expand=F)

# make the scale prettier
tb$count <- 10 + tb$count



poa <- ggplot(tb) + polar() + # set the data and the polar layout
  # draw circles on the plot for the r-axis
  geom_line(data=circle, aes(theta,rho), color="gray") +
  geom_line(data=circle2, aes(theta,rho), color="gray") +
  # add the data
  geom_point(aes(x=theta, y=count), size=1) + # add the bearings as stacked dots
  geom_segment(aes(x=bearing, y=0, xend=bearing, yend=r*10, color=signif), size = .6, data=circStatsAge) + # add the mean lines, w/ type determined by significance
  scale_color_manual(values=c("gray", "black")) + # set black line to be significant
  facet_wrap(~age) +
  # groom the plot a bit
  scale_x_continuous("", limits=c(0,360), breaks=seq(0,360-1,by=45), labels=c("N","","","","S","","","")) + # set the X axis (angle axis)
  theme_bw() + labs(color="Statistically Significant") + # remove the gray background and adjust legend label
  theme(legend.title = element_text(size=8)) + theme(legend.text = element_text(size=6)) + # adjust legend text sizes
  theme(panel.grid.minor = element_line(size = .2, color = "white")) + theme(panel.grid.major = element_line(size=.5, color = "white")) + # remove x-axis grid, make y-axis grid black
  theme(plot.title = element_text(face="bold", size=18)) + theme(axis.title.y = element_text(size=8, face="bold")) + # format the plot title and r-axis label
  theme(axis.text.x = element_text(size=5, vjust=.05)) + theme(axis.text.y = element_blank(), axis.ticks = element_blank()) + # adjust axis text
  theme(strip.text.x = element_text(size=6, face="bold"), strip.background = element_rect(colour="black", fill="#CCCCFF"))  # format facet labels

poa


# try larger groups
tg <- dneT
tg$ageG <- 0
tg[tg$fish_age..dph. <= 10,c("ageG")] <- "10&Under"
tg[tg$fish_age..dph. > 20,c("ageG")] <- "22&Over"
tg[tg$ageG == 0, c("ageG")] <- "12to20"
tg$ageG <- as.factor(tg$ageG)


circStatsAgeG <- ddply(tg,.(ageG),circStats) # run the statistics
circStatsAgeG$age <- circStatsAgeG$ageG
for (i in 1:nrow(circStatsAgeG)){  
  if (circStatsAgeG[i,c("bearing")]<=0){
    circStatsAgeG[i,c("bearing")] <- 360 + circStatsAgeG[i,c("bearing")]
  }
}

# bin angles for plot
tg$theta <- as.numeric(round_any(tg$mean, bin))
tg$theta[tg$theta==360] <- 0
# create a data.frame with count
counts <- count(tg, vars=c("theta","ageG"))

# repeat each point the appropriate number of times
tbg <- adply(counts, 1, function(x) {
  data.frame(theta=x$theta,age=x$ageG, count=1:x$freq)
}, .expand=F)

# make the scale prettier
tbg$count <- 10 + tbg$count



pga <- ggplot(tbg) + polar() + # set the data and the polar layout
  # draw circles on the plot for the r-axis
  geom_line(data=circle, aes(theta,rho), color="gray") +
  geom_line(data=circle2, aes(theta,rho), color="gray") +
  # add the data
  geom_point(aes(x=theta, y=count), size=1) + # add the bearings as stacked dots
  geom_segment(aes(x=bearing, y=0, xend=bearing, yend=r*10, color=signif), size = .6, data=circStatsAgeG) + # add the mean lines, w/ type determined by significance
  scale_color_manual(values=c("gray", "black")) + # set black line to be significant
  facet_wrap(~age) +
  # groom the plot a bit
  scale_x_continuous("", limits=c(0,360), breaks=seq(0,360-1,by=45), labels=c("N","","","","S","","","")) + # set the X axis (angle axis)
  theme_bw() + labs(color="Statistically Significant") + # remove the gray background and adjust legend label
  theme(legend.title = element_text(size=8)) + theme(legend.text = element_text(size=6)) + # adjust legend text sizes
  theme(panel.grid.minor = element_line(size = .2, color = "white")) + theme(panel.grid.major = element_line(size=.5, color = "white")) + # remove x-axis grid, make y-axis grid black
  theme(plot.title = element_text(face="bold", size=18)) + theme(axis.title.y = element_text(size=8, face="bold")) + # format the plot title and r-axis label
  theme(axis.text.x = element_text(size=5, vjust=.05)) + theme(axis.text.y = element_blank(), axis.ticks = element_blank()) + # adjust axis text
  theme(strip.text.x = element_text(size=6, face="bold"), strip.background = element_rect(colour="black", fill="#CCCCFF"))  # format facet labels

pga





## Precision of orientation by environmental context -----------------

# SUPMAT PLOT

rmodel <- lm(r ~ tide+depth..m.+location, data=d2)
summary(rmodel)

pr <- ggplot(d2) +
  geom_boxplot(aes(tide,r,fill=location,linetype=depth..m.)) +
  theme_bw() + 
  ylab("Rho-Value") + xlab("Tide")
pr





## Precision of orientation and age ---------------

dneT <- do.call(data.frame,lapply(dneT, function(x) replace(x, is.infinite(x),NA)))

rmodel2 <- lm(r ~ fish_age..dph.+speed.mean, data=dneT)
summary(rmodel2)
# significant finding - r value decreases with age

ggplot(dneTs,aes(y=r,x=fish_age..dph.,color=speed.mean)) + geom_point() + geom_smooth(method="lm") + 
  scale_y_continuous(limits=c(0,1)) + scale_colour_gradient(limits=c(0, 2), low="blue", high="red")





## Swimming speed by environmental context ----------------

# SUPMAT PLOT

d22 <- do.call(data.frame,lapply(d2, function(x) replace(x, is.infinite(x),NA)))

smodel <- lm(speed.mean ~ tide+depth..m.+location, data=d22)
summary(smodel)

ps <- ggplot(d22) +
  geom_boxplot(aes(tide,speed.mean,fill=location,linetype=depth..m.)) +
  scale_y_continuous(limits=c(0,1.25)) +
  theme_bw() + ylab("Activity Level (cm/s)") + xlab("Tide")
ps



## look at effect of ontogeny/temp on speed
dneT <- do.call(data.frame,lapply(dneT, function(x) replace(x, is.infinite(x),NA)))

tempsmodel <- lm(speed.mean ~ meanTemp + fish_age..dph.,data=dneTs)
summary(tempsmodel)

ggplot(dneT,aes(y=speed.mean,x=meanTemp,color=fish_age..dph.)) + geom_point() + geom_smooth(method="lm") + 
  scale_y_continuous(limits=c(0,1)) + scale_colour_gradient(limits=c(2, 30), low="blue", high="red")




## Swimming speed by age ------------------

dneT <- do.call(data.frame,lapply(dneT, function(x) replace(x, is.infinite(x),NA)))
dneTs <- dneT[dneT$speed.mean <=2,]

smodel2 <- lm(speed.mean ~ fish_age..dph., data=dneTs)
summary(smodel2)
# significant finding - speed increases with age

ggplot(dneTs,aes(y=speed.mean,x=fish_age..dph.)) + geom_point() + 
  geom_smooth(method="lm",color="black",linetype=3) +
  theme_bw() + xlab("Fish Age (DPH)") + ylab("Activity Level (cm/s)")


smodel3 <- lm(speed.mean ~ fish_age..dph. + drift_distance, data=dneTs)
summary(smodel3)


## Relation between r value and speed --------------

dneTs <- dneT[dneT$speed.mean <=2,]

rsmodel <- lm(r ~ speed.mean, data=dneTs)
summary(rsmodel)
# significant finding - increasing speed decreases r value

ggplot(dneTs,aes(y=r,x=speed.mean)) + geom_point() + geom_smooth(method=lm) #+ scale_x_continuous(limits=c(0,1.3))





## Precision by age (and speed as controlling factor) -----------

rsamodel <- lm(r ~ fish_age..dph.*speed.mean, data=dneTs)
summary(rsamodel)
# significant finding - r value decreases with age

ggplot(dneTs,aes(y=r,x=speed.mean,color=fish_age..dph.)) + geom_point() + geom_smooth(method="lm") + 
  scale_y_continuous(limits=c(0,1)) + scale_colour_gradient(limits=c(2, 30), low="blue", high="red")





## Turning angle by environmental context --------------

# define function to look at turning angles instead of position
circStatsTurn <- function(x) { # define descriptive statistiics function
  n = nrow(x)
  meanTurn = mean.circular(x$turn.abs.mean,na.rm=TRUE)
  sdTurn = sd.circular(x$turn.abs.mean,na.rm=TRUE)
  meanR = mean(x$r, na.rm=T)
  sdR = sd(x$r, na.rm=T)
  ray = rayleigh.test(na.omit(x$turn.abs.mean))
  r = ray$statistic
  p = ray$p.value
  signif = p < 0.05
  return(data.frame(n, meanTurn,sdTurn, meanR, sdR, r, p, signif))
}

# run it
circStatsA <- ddply(d2,.(depth..m.,tide,location),circStatsTurn)

# bin the turn angles
db3 <- ddply(d2,.(depth..m.,tide,location), function(x, bin) {  # 'splitting' by method here b/c theres only 1 method and its easier than changing code to remove ddply
  x$mean <- as.numeric(round_any(x$turn.abs.mean, bin))
  x$mean[x$mean==360] <- 0
  
  # create a data.frame with count
  counts <- count(x, "mean")
  
  # repeat each point the appropriate number of times
  d <- adply(counts, 1, function(x) {
    data.frame(theta=x$mean, count=1:x$freq)
  }, .expand=F)
  
  # make the scale prettier
  d$count <- 10 + d$count
  
  return(d)
}, bin=bin)

# plot them by condition
pa <- ggplot(db3) + polar() + # set the data and the polar layout
  # draw circles on the plot for the r-axis
  geom_line(data=circle, aes(theta,rho), color="gray") +
  geom_line(data=circle2, aes(theta,rho), color="gray") +
  # add the data
  geom_point(aes(x=theta, y=count), size=1) + # add the bearings as stacked dots
  geom_segment(aes(x=meanTurn, y=0, xend=meanTurn, yend=r*10, color=signif), size = .6, data=circStatsA) + # add the mean lines, w/ type determined by significance
  scale_color_manual(values=c("gray", "black")) + # set black line to be significant
  facet_wrap(~tide+location+depth..m.) +
  # groom the plot a bit
  scale_x_continuous("", limits=c(0,360), breaks=seq(0,360-1,by=45), labels=c("N","","","","S","","","")) + # set the X axis (angle axis)
  theme_bw() + labs(color="Statistically Significant") + # remove the gray background and adjust legend label
  theme(legend.title = element_text(size=8)) + theme(legend.text = element_text(size=6)) + # adjust legend text sizes
  theme(panel.grid.minor = element_line(size = .2, color = "white")) + theme(panel.grid.major = element_line(size=.5, color = "white")) + # remove x-axis grid, make y-axis grid black
  theme(plot.title = element_text(face="bold", size=18)) + theme(axis.title.y = element_text(size=8, face="bold")) + # format the plot title and r-axis label
  theme(axis.text.x = element_text(size=5, vjust=.05)) + theme(axis.text.y = element_blank(), axis.ticks = element_blank()) + # adjust axis text
  theme(strip.text.x = element_text(size=6, face="bold"), strip.background = element_rect(colour="black", fill="#CCCCFF"))  # format facet labels

pa


# pairwise comparison of angles

dne <- subset(d2, depth..m.=="18" & location=="near" & tide=="ebb")
dnf <- subset(d2, depth..m.=="18" & location=="near" & tide=="flood")
dfe <- subset(d2, depth..m.=="18" & location=="far" & tide=="ebb")
sne <- subset(d2, depth..m.=="9" & location=="near" & tide=="ebb")

watson.two.test(dne$turn.abs.mean, dnf$turn.abs.mean)
watson.two.test(dne$turn.abs.mean, dfe$turn.abs.mean) 
# ^^^ significant finding!! larger turning angles farther from the reef
watson.two.test(dne$turn.abs.mean, sne$turn.abs.mean)
watson.two.test(dnf$turn.abs.mean, dfe$turn.abs.mean)
watson.two.test(dnf$turn.abs.mean, sne$turn.abs.mean)
watson.two.test(dfe$turn.abs.mean, sne$turn.abs.mean)



# plot for the significant pairwise comparison


da <- subset(d2, expcode == 1 | expcode == 3)


da2 <- ddply(da,.(expcode), function(x, bin) {  # 'splitting' by method here b/c theres only 1 method and its easier than changing code to remove ddply
  x$mean <- as.numeric(round_any(x$turn.abs.mean, bin))
  x$mean[x$mean==360] <- 0
  
  # create a data.frame with count
  counts <- count(x, "mean")
  
  # repeat each point the appropriate number of times
  d <- adply(counts, 1, function(x) {
    data.frame(theta=x$mean, count=1:x$freq)
  }, .expand=F)
  
  # make the scale prettier
  d$count <- 10 + d$count
  
  return(d)
}, bin=bin)

circStatsda <- ddply(da,.(expcode),circStatsTurn) # run the statistics

pa2 <- ggplot(da2) + polar() + # set the data and the polar layout
  # draw circles on the plot for the r-axis
  geom_line(data=circle, aes(theta,rho), color="gray") +
  geom_line(data=circle2, aes(theta,rho), color="gray") +
  # add the data
  geom_point(aes(x=theta, y=count), size=1) + # add the bearings as stacked dots
  geom_segment(aes(x=meanTurn, y=0, xend=meanTurn, yend=r*10, color=signif), size = .6, data=circStatsda) + # add the mean lines, w/ type determined by significance
  scale_color_manual(values=c("gray", "black")) + # set black line to be signif
  facet_grid(expcode~.) +
  # groom the plot a bit
  scale_x_continuous("", limits=c(0,360), breaks=seq(0,360-1,by=45), labels=c("Straight","","","","Reverse","","","")) + # set the X axis (angle axis)
  theme_bw() + labs(color="Statistically Significant") + # remove the gray background and adjust legend label
  theme(legend.title = element_text(size=8)) + theme(legend.text = element_text(size=6)) + # adjust legend text sizes
  theme(panel.grid.minor = element_line(size = .2, color = "white")) + theme(panel.grid.major = element_line(size=.5, color = "white")) + # remove x-axis grid, make y-axis grid black
  theme(plot.title = element_text(face="bold", size=18)) + theme(axis.title.y = element_text(size=8, face="bold")) + # format the plot title and r-axis label
  theme(axis.text.x = element_text(size=8, vjust=.05,face="bold")) + theme(axis.text.y = element_blank(), axis.ticks = element_blank()) + # adjust axis text
  # theme(strip.text.x = element_text(size=6, face="bold"), strip.background = element_rect(colour="black", fill="#CCCCFF"))  # format facet labels
  theme(strip.background = element_blank(), strip.text.x = element_blank()) # or remove them

pa2






## Turning angles and age --------------

d = d[!is.na(d$turn.abs.mean),] 

dneT <- subset(d, depth..m.=="18" & location=="near" & tide=="ebb")

clt <- lm.circular(y=dneT$turn.abs.mean,x=dneT$fish_age..dph.,type="c-l",init=0)
clt$p.values
# significant finding - turning angles correlated with age

circStatsdt <- circStatsTurn(dneT) # run the statistics
dneT$meanTurn <- dneT$turn.abs.mean
dneT$count <- 11

pda <- ggplot(dneT) + polar() + # set the data and the polar layout
  # draw circles on the plot for the r-axis
  geom_line(data=circle, aes(theta,rho), color="gray") +
  geom_line(data=circle2, aes(theta,rho), color="gray") +
  # add the data
  geom_point(aes(x=meanTurn, y=count,fill=fish_age..dph.),pch=21, size=3) + # add the bearings as stacked dots
  scale_fill_gradient(limits=c(2, 30), low="blue", high="red") +
  geom_segment(aes(x=meanTurn, y=0, xend=meanTurn, yend=r*10, color=signif), size = .6, data=circStatsdt) + # add the mean lines, w/ type determined by significance
  scale_color_manual(values=c("gray", "black")) + # set black line to be significant
  # groom the plot a bit
  scale_x_continuous("", limits=c(0,360), breaks=seq(0,360-1,by=45), labels=c("N","","","","S","","","")) + # set the X axis (angle axis)
  theme_bw() + labs(color="Statistically Significant") + # remove the gray background and adjust legend label
  theme(legend.title = element_text(size=8)) + theme(legend.text = element_text(size=6)) + # adjust legend text sizes
  theme(panel.grid.minor = element_line(size = .2, color = "white")) + theme(panel.grid.major = element_line(size=.5, color = "white")) + # remove x-axis grid, make y-axis grid black
  theme(plot.title = element_text(face="bold", size=18)) + theme(axis.title.y = element_text(size=8, face="bold")) + # format the plot title and r-axis label
  theme(axis.text.x = element_text(size=5, vjust=.05)) + theme(axis.text.y = element_blank(), axis.ticks = element_blank()) + # adjust axis text
  theme(strip.text.x = element_text(size=6, face="bold"), strip.background = element_rect(colour="black", fill="#CCCCFF"))  # format facet labels

pda


tg <- dneT
tg$ageG <- 0
tg[tg$fish_age..dph. <= 10,c("ageG")] <- "10&Under"
tg[tg$fish_age..dph. > 20,c("ageG")] <- "22&Over"
tg[tg$ageG == 0, c("ageG")] <- "12to20"
tg$ageG <- as.factor(tg$ageG)


circStatsTurn2 <- ddply(tg,.(ageG),circStatsTurn) # run the statistics
circStatsTurn2$age <- circStatsTurn2$ageG
for (i in 1:nrow(circStatsTurn2)){  
  if (circStatsTurn2[i,c("bearing")]<=0){
    circStatsTurn2[i,c("bearing")] <- 360 + circStatsTurn2[i,c("bearing")]
  }
}

# bin angles for plot
tg$theta <- as.numeric(round_any(tg$turn.abs.mean, bin))
tg$theta[tg$theta==360] <- 0
# create a data.frame with count
counts <- count(tg, vars=c("theta","ageG"))

# repeat each point the appropriate number of times
tbg <- adply(counts, 1, function(x) {
  data.frame(theta=x$theta,age=x$ageG, count=1:x$freq)
}, .expand=F)

# make the scale prettier
tbg$count <- 10 + tbg$count



pta <- ggplot(tbg) + polar() + # set the data and the polar layout
  # draw circles on the plot for the r-axis
  geom_line(data=circle, aes(theta,rho), color="gray") +
  geom_line(data=circle2, aes(theta,rho), color="gray") +
  # add the data
  geom_point(aes(x=theta, y=count), size=1) + # add the bearings as stacked dots
  geom_segment(aes(x=meanTurn, y=0, xend=meanTurn, yend=r*10, color=signif), size = .6, data=circStatsTurn2) + # add the mean lines, w/ type determined by significance
  scale_color_manual(values=c("gray", "black")) + # set black line to be significant
  facet_wrap(~age) +
  # groom the plot a bit
  scale_x_continuous("", limits=c(0,360), breaks=seq(0,360-1,by=45), labels=c("N","","","","S","","","")) + # set the X axis (angle axis)
  theme_bw() + labs(color="Statistically Significant") + # remove the gray background and adjust legend label
  theme(legend.title = element_text(size=8)) + theme(legend.text = element_text(size=6)) + # adjust legend text sizes
  theme(panel.grid.minor = element_line(size = .2, color = "white")) + theme(panel.grid.major = element_line(size=.5, color = "white")) + # remove x-axis grid, make y-axis grid black
  theme(plot.title = element_text(face="bold", size=18)) + theme(axis.title.y = element_text(size=8, face="bold")) + # format the plot title and r-axis label
  theme(axis.text.x = element_text(size=5, vjust=.05)) + theme(axis.text.y = element_blank(), axis.ticks = element_blank()) + # adjust axis text
  theme(strip.text.x = element_text(size=6, face="bold"), strip.background = element_rect(colour="black", fill="#CCCCFF"))  # format facet labels

pta


tg10 <- subset(tg, ageG=="10&Under")
tg20 <- subset(tg, ageG=="12to20")
tg30 <- subset(tg, ageG=="22&Over")

watson.two.test(tg10$turn.abs.mean, tg20$turn.abs.mean) 
watson.two.test(tg10$turn.abs.mean, tg30$turn.abs.mean) 
watson.two.test(tg20$turn.abs.mean, tg30$turn.abs.mean)



tmodel <- lm(as.numeric(turn.abs.mean) ~ fish_age..dph., data=dneT)
summary(tmodel)
    
ggplot(dneT,aes(y=as.numeric(turn.abs.mean),x=fish_age..dph.)) +
  geom_point() + geom_smooth(method=lm,color="black",linetype=3) +
  xlab("Fish Age (DPH)") + ylab("Average Turning Angle (abs value)") + theme_bw()




## Orientation to the current direction ----------
d2 <- subset(d, fish_age..dph. >17 & fish_age..dph. <25)

# correct angles so headings are relative to current direction
d2$current_bearing <- d2$mean - d2$drift_direction
d2 = d2[!is.na(d2$current_bearing),] # remove trials with no current data

# correct for negative angle values produced
for (i in 1:nrow(d2)){  
  if (d2[i,c("current_bearing")]<=0){
    d2[i,c("current_bearing")] <- 360 + d2[i,c("current_bearing")]
  }
}

# define function to look at current bearings
circStatsC <- function(x) {
  n = nrow(x)
  currentBearing = mean.circular(x$current_bearing)
  sdbearing = sd.circular(x$current_bearing)
  meanR = mean(x$r, na.rm=T)
  sdR = sd(x$r, na.rm=T)
  ray = rayleigh.test(na.omit(x$current_bearing))
  r = ray$statistic
  p = ray$p.value
  signif = p < 0.05
  return(data.frame(n, currentBearing,sdbearing, meanR, sdR, r, p, signif))
}



dbC <- ddply(d2,.(depth..m.,tide,location), function(x, bin) {  # 'splitting' by method here b/c theres only 1 method and its easier than changing code to remove ddply
  x$mean <- as.numeric(round_any(x$current_bearing, bin))
  x$mean[x$mean==360] <- 0
  
  # create a data.frame with count
  counts <- count(x, "mean")
  
  # repeat each point the appropriate number of times
  d <- adply(counts, 1, function(x) {
    data.frame(theta=x$mean, count=1:x$freq)
  }, .expand=F)
  
  # make the scale prettier
  d$count <- 10 + d$count
  
  return(d)
}, bin=bin)

circStatsCur <- ddply(d2,.(depth..m.,tide,location),circStatsC) 
# significant finding! swim against current in the flood tide
for (i in 1:nrow(circStatsCur)){  
  if (circStatsCur[i,c("currentBearing")]<=0){
    circStatsCur[i,c("currentBearing")] <- 360 + circStatsCur[i,c("currentBearing")]
  }
}

pC <- ggplot(dbC) + polar() + # set the data and the polar layout
  # draw circles on the plot for the r-axis
  geom_line(data=circle, aes(theta,rho), color="gray") +
  geom_line(data=circle2, aes(theta,rho), color="gray") +
  # add the data
  geom_point(aes(x=theta, y=count), size=1) + # add the bearings as stacked dots
  geom_segment(aes(x=currentBearing, y=0, xend=currentBearing, yend=r*10, color=signif), size = .6, data=circStatsCur) + # add the mean lines, w/ type determined by significance
  scale_color_manual(values=c("gray", "black")) + # set black line to be significant
  facet_wrap(~tide+location+depth..m.) +
  # groom the plot a bit
  scale_x_continuous("", limits=c(0,360), breaks=seq(0,360-1,by=45), labels=c("N","","","","S","","","")) + # set the X axis (angle axis)
  theme_bw() + labs(color="Statistically Significant") + # remove the gray background and adjust legend label
  theme(legend.title = element_text(size=8)) + theme(legend.text = element_text(size=6)) + # adjust legend text sizes
  theme(panel.grid.minor = element_line(size = .2, color = "white")) + theme(panel.grid.major = element_line(size=.5, color = "white")) + # remove x-axis grid, make y-axis grid black
  theme(plot.title = element_text(face="bold", size=18)) + theme(axis.title.y = element_text(size=8, face="bold")) + # format the plot title and r-axis label
  theme(axis.text.x = element_text(size=5, vjust=.05)) + theme(axis.text.y = element_blank(), axis.ticks = element_blank()) + # adjust axis text
  theme(strip.text.x = element_text(size=6, face="bold"), strip.background = element_rect(colour="black", fill="#CCCCFF"))  # format facet labels

pC


# for far / shallow conditions only - SUPMAT PLOT
dbC2 <- subset(dbC, depth..m. == "9" | location == "far")
circStatsCur2 <- subset(circStatsCur, depth..m. == "9" | location == "far")



pC2 <- ggplot(dbC2) + polar() + # set the data and the polar layout
  # draw circles on the plot for the r-axis
  geom_line(data=circle, aes(theta,rho), color="gray") +
  geom_line(data=circle2, aes(theta,rho), color="gray") +
  # add the data
  geom_point(aes(x=theta, y=count), size=1) + # add the bearings as stacked dots
  geom_segment(aes(x=currentBearing, y=0, xend=currentBearing, yend=r*10, color=signif), size = .6, data=circStatsCur2) + # add the mean lines, w/ type determined by significance
  scale_color_manual(values=c("gray", "black")) + # set black line to be significant
  facet_wrap(~tide+location+depth..m.) +
  # groom the plot a bit
  theme_bw() + ylab("Rho-Value") +
  theme(panel.grid.minor = element_line(size = .2, color = "white")) + 
  theme(panel.grid.major = element_line(size=.5, color = "white")) +
  theme(axis.title.y = element_text(size=16, face="bold")) + 
  theme(axis.text.y = element_text(size=16)) +
  scale_y_continuous(breaks=c(0,3,6,10),labels=c(0,0.3,0.6,1)) +
  theme(axis.text.x = element_text(face="bold",size=12, vjust=.05)) + 
  scale_x_continuous("", limits=c(0,360), breaks=seq(0,360-1,by=45), labels=c("With Current","","","","Against Current","","","")) + 
  theme(strip.background = element_blank(),strip.text.x = element_blank())

  
pC2




# for flood condition only

dnf <- d2[d2$depth..m.=="18" & d2$location=="near" & d2$tide=="flood",]

dbcdnf <- ddply(dnf,.(expcode), function(x, bin) {  # 'splitting' by method here b/c theres only 1 method and its easier than changing code to remove ddply
  x$mean <- as.numeric(round_any(x$current_bearing, bin))
  x$mean[x$mean==360] <- 0
  
  # create a data.frame with count
  counts <- count(x, "mean")
  
  # repeat each point the appropriate number of times
  d <- adply(counts, 1, function(x) {
    data.frame(theta=x$mean, count=1:x$freq)
  }, .expand=F)
  
  # make the scale prettier
  d$count <- 10 + d$count
  
  return(d)
}, bin=bin)

circStatsCurdnf <- circStatsC(dnf) # run the statistics


pcdnf <- ggplot(dbcdnf) + polar() + # set the data and the polar layout
  # draw circles on the plot for the r-axis
  geom_line(data=circle, aes(theta,rho), color="gray") +
  geom_line(data=circle2, aes(theta,rho), color="gray") +
  # add the data
  geom_point(aes(x=theta, y=count), size=1) + # add the bearings as stacked dots
  geom_segment(aes(x=currentBearing, y=0, xend=currentBearing, yend=r*10, color=signif), size = .6, data=circStatsCurdnf) + # add the mean lines, w/ type determined by significance
  scale_color_manual(values=c("gray", "black")) + # set black line to be signif
  # groom the plot a bit
  scale_x_continuous("", limits=c(0,360), breaks=seq(0,360-1,by=45), labels=c("With Current","","","","Against Current","","","")) + # set the X axis (angle axis)
  theme_bw() + labs(color="Statistically Significant") + # remove the gray background and adjust legend label
  theme(legend.title = element_text(size=8)) + theme(legend.text = element_text(size=6)) + # adjust legend text sizes
  theme(panel.grid.minor = element_line(size = .2, color = "white")) + theme(panel.grid.major = element_line(size=.5, color = "white")) + # remove x-axis grid, make y-axis grid black
  theme(plot.title = element_text(face="bold", size=18)) + theme(axis.title.y = element_text(size=8, face="bold")) + # format the plot title and r-axis label
  theme(axis.text.x = element_text(size=10, vjust=.05,face="bold")) + theme(axis.text.y = element_blank(), axis.ticks = element_blank()) + # adjust axis text
  # theme(strip.text.x = element_text(size=6, face="bold"), strip.background = element_rect(colour="black", fill="#CCCCFF"))  # format facet labels
  theme(strip.background = element_blank(), strip.text.x = element_blank()) # or remove them

pcdnf




# combine flood and ebb (not shallow or far)
dnf <- d2[d2$depth..m.=="18" & d2$location=="near" & d2$tide=="flood",]
dne <- d2[d2$depth..m.=="18" & d2$location=="near" & d2$tide=="ebb",]

dbz <- rbind(dnf,dne)
dbz$expcode <- 1

dbcdbz <- ddply(dbz,.(expcode), function(x, bin) {  # 'splitting' by method here b/c theres only 1 method and its easier than changing code to remove ddply
  x$mean <- as.numeric(round_any(x$current_bearing, bin))
  x$mean[x$mean==360] <- 0
  
  # create a data.frame with count
  counts <- count(x, "mean")
  
  # repeat each point the appropriate number of times
  d <- adply(counts, 1, function(x) {
    data.frame(theta=x$mean, count=1:x$freq)
  }, .expand=F)
  
  # make the scale prettier
  d$count <- 10 + d$count
  
  return(d)
}, bin=bin)

circStatsCurdbz <- circStatsC(dbz) # run the statistics
for (i in 1:nrow(circStatsCurdbz)){  
  if (circStatsCurdbz[i,c("currentBearing")]<=0){
    circStatsCurdbz[i,c("currentBearing")] <- 360 + circStatsCurdbz[i,c("currentBearing")]
  }
}

pcdbz <- ggplot(dbcdbz) + polar() + # set the data and the polar layout
  # draw circles on the plot for the r-axis
  geom_line(data=circle, aes(theta,rho), color="gray") +
  geom_line(data=circle2, aes(theta,rho), color="gray") +
  # add the data
  geom_point(aes(x=theta, y=count), size=1) + # add the bearings as stacked dots
  geom_segment(aes(x=currentBearing, y=0, xend=currentBearing, yend=r*10, color=signif), size = .6, data=circStatsCurdbz) + # add the mean lines, w/ type determined by significance
  scale_color_manual(values=c("gray", "black")) + # set black line to be signif
  # groom the plot a bit
  scale_x_continuous("", limits=c(0,360), breaks=seq(0,360-1,by=45), labels=c("With Current","","","","Against Current","","","")) + # set the X axis (angle axis)
  theme_bw() + labs(color="Statistically Significant") + # remove the gray background and adjust legend label
  theme(legend.title = element_text(size=8)) + theme(legend.text = element_text(size=6)) + # adjust legend text sizes
  theme(panel.grid.minor = element_line(size = .2, color = "white")) + theme(panel.grid.major = element_line(size=.5, color = "white")) + # remove x-axis grid, make y-axis grid black
  theme(plot.title = element_text(face="bold", size=18)) + theme(axis.title.y = element_text(size=8, face="bold")) + # format the plot title and r-axis label
  theme(axis.text.x = element_text(size=10, vjust=.05,face="bold")) + theme(axis.text.y = element_blank(), axis.ticks = element_blank()) + # adjust axis text
  # theme(strip.text.x = element_text(size=6, face="bold"), strip.background = element_rect(colour="black", fill="#CCCCFF"))  # format facet labels
  theme(strip.background = element_blank(), strip.text.x = element_blank()) # or remove them

pcdbz


## PUBLICATION PLOT - Orientation to Current in DNE/DNF


pcurO <- ggplot(dbcdbz) + polar() + # set the data and the polar layout
  # draw circles on the plot for the r-axis
  geom_line(data=circle, aes(theta,rho), color="gray",size=1.5) +
  geom_line(data=circle2, aes(theta,rho), color="gray",size=1.5) +
  # add the data
  geom_point(aes(x=theta, y=count), size=4) + 
  geom_segment(aes(x=currentBearing, y=0, xend=currentBearing, yend=r*10), size = 1, data=circStatsCurdbz) + 
  scale_color_manual(values=c("gray", "black")) + 
  # groom the plot
  theme_bw() + 
  theme(panel.grid.minor = element_line(size = .2, color = "white")) + 
  theme(panel.grid.major = element_line(size=.5, color = "white")) +
  theme(axis.title.y = element_text(size=18, face="bold")) + 
  theme(axis.text.y = element_text(size=16)) +
  scale_y_continuous(breaks=c(0,3,6,10),labels=c(0,0.3,0.6,1)) +
  theme(axis.text.x = element_text(face="bold",size=22, vjust=.05)) + 
  scale_x_continuous("", limits=c(0,360), breaks=seq(0,360-1,by=45), labels=c("With Current","","","","Against Current","","","")) 

pcurO



## Current Patterns by Environmental Condition ---------------

# define function to look at current direction
circStatsCr <- function(x) {
  n = nrow(x)
  currentDirection = mean.circular(x$drift_direction,na.rm=TRUE)
  sdbearing = sd.circular(x$drift_direction)
  meanR = mean(x$r, na.rm=T)
  sdR = sd(x$r, na.rm=T)
  ray = rayleigh.test(na.omit(x$drift_direction))
  r = ray$statistic
  p = ray$p.value
  signif = p < 0.05
  return(data.frame(n, currentDirection,sdbearing, meanR, sdR, r, p, signif))
}


# create binned data
dbCr <- ddply(d2,.(depth..m.,tide,location), function(x, bin) {  
  x$mean <- as.numeric(round_any(x$drift_direction, bin))
  x$mean[x$mean==360] <- 0
  
  # create a data.frame with count
  counts <- count(x, "mean")
  
  # repeat each point the appropriate number of times
  d <- adply(counts, 1, function(x) {
    data.frame(theta=x$mean, count=1:x$freq)
  }, .expand=F)
  
  # make the scale prettier
  d$count <- 10 + d$count
  
  return(d)
}, bin=bin)


# run stats and correct for negative angles
circStatsCurr <- ddply(d2,.(depth..m.,tide,location),circStatsCr) 
# signifcant! swim against current in the flood tide
for (i in 1:nrow(circStatsCurr)){  
  if (circStatsCurr[i,c("currentDirection")]<=0){
    circStatsCurr[i,c("currentDirection")] <- 360 + circStatsCurr[i,c("currentDirection")]
  }
}


# plot
pCr <- ggplot(dbCr) + polar() + # set the data and the polar layout
  # draw circles on the plot for the r-axis
  geom_line(data=circle, aes(theta,rho), color="gray") +
  geom_line(data=circle2, aes(theta,rho), color="gray") +
  # add the data
  geom_point(aes(x=theta, y=count), size=.7) + # add the bearings as stacked dots
  geom_segment(aes(x=currentDirection, y=0, xend=currentDirection, yend=r*10, color=signif), size = .6, data=circStatsCurr) + # add the mean lines, w/ type determined by significance
  scale_color_manual(values=c("gray", "black")) + # set black line to be significant
  facet_wrap(~tide+location+depth..m.) +
  # groom the plot a bit
  scale_x_continuous("", limits=c(0,360), breaks=seq(0,360-1,by=45), labels=c("N","","","","S","","","")) + # set the X axis (angle axis)
  theme_bw() + labs(color="Statistically Significant") + # remove the gray background and adjust legend label
  theme(legend.title = element_text(size=8)) + theme(legend.text = element_text(size=6)) + # adjust legend text sizes
  theme(panel.grid.minor = element_line(size = .2, color = "white")) + theme(panel.grid.major = element_line(size=.5, color = "white")) + # remove x-axis grid, make y-axis grid black
  theme(plot.title = element_text(face="bold", size=18)) + theme(axis.title.y = element_text(size=8, face="bold")) + # format the plot title and r-axis label
  theme(axis.text.x = element_text(size=5, vjust=.05)) + theme(axis.text.y = element_blank(), axis.ticks = element_blank()) + # adjust axis text
  theme(strip.text.x = element_text(size=6, face="bold"), strip.background = element_rect(colour="black", fill="#CCCCFF"))  # format facet labels

pCr
# significant finding - all current N-S except far is just north


# for only the flood/reference conditions - PUBLICATION PLOT
dnf <- d2[d2$depth..m.=="18" & d2$location=="near" & d2$tide=="flood",]
dne <- d2[d2$depth..m.=="18" & d2$location=="near" & d2$tide=="ebb",]

dbz <- rbind(dnf,dne)
dbz$expcode <- 0

dbcrdbz <- ddply(dbz,.(tide), function(x, bin) {  # 'splitting' by method here b/c theres only 1 method and its easier than changing code to remove ddply
  x$mean <- as.numeric(round_any(x$drift_direction, bin))
  x$mean[x$mean==360] <- 0
  
  # create a data.frame with count
  counts <- count(x, "mean")
  
  # repeat each point the appropriate number of times
  d <- adply(counts, 1, function(x) {
    data.frame(theta=x$mean, count=1:x$freq)
  }, .expand=F)
  
  # make the scale prettier
  d$count <- 10 + d$count
  
  return(d)
}, bin=bin)

circStatsCurrdnf <- ddply(dbz, ~tide, circStatsCr) # run the statistics
circStatsCurrdnf$currentDirection <- 360 + circStatsCurrdnf$currentDirection

pcrdnf <- ggplot(dbcrdbz) + polar() + # set the data and the polar layout
  # draw circles on the plot for the r-axis
  geom_line(data=circle, aes(theta,rho), color="gray") +
  geom_line(data=circle2, aes(theta,rho), color="gray") +
  # add the data
  geom_point(aes(x=theta, y=count), size=1) + # add the bearings as stacked dots
  geom_segment(aes(x=currentDirection, y=0, xend=currentDirection, yend=r*10, color=signif), size = .6, data=circStatsCurrdnf) + # add the mean lines, w/ type determined by significance
  scale_color_manual(values=c("gray", "black")) + # set black line to be signif
  # groom the plot a bit
  facet_wrap(~tide) +
  scale_x_continuous("", limits=c(0,360), breaks=seq(0,360-1,by=45), labels=c("N","","","","S","","","")) + # set the X axis (angle axis)
  theme_bw() + labs(color="Statistically Significant") + # remove the gray background and adjust legend label
  theme(legend.title = element_text(size=8)) + theme(legend.text = element_text(size=6)) + # adjust legend text sizes
  theme(panel.grid.minor = element_line(size = .2, color = "white")) + theme(panel.grid.major = element_line(size=.5, color = "white")) + # remove x-axis grid, make y-axis grid black
  theme(plot.title = element_text(face="bold", size=18)) + theme(axis.title.y = element_text(size=8, face="bold")) + # format the plot title and r-axis label
  theme(axis.text.x = element_text(size=10, vjust=.05,face="bold")) + theme(axis.text.y = element_blank(), axis.ticks = element_blank())  # adjust axis text
  # theme(strip.text.x = element_text(size=6, face="bold"), strip.background = element_rect(colour="black", fill="#CCCCFF"))  # format facet labels
  # theme(strip.background = element_blank(), strip.text.x = element_blank()) # or remove them

pcrdnf


## PUBLICATION PLOT - Current Directions in DNE / DNF

pcurC <- ggplot(dbcrdbz) + polar() + # set the data and the polar layout
  # draw circles on the plot for the r-axis
  geom_line(data=circle, aes(theta,rho), color="gray",size=1.5) +
  geom_line(data=circle2, aes(theta,rho), color="gray",size=1.5) +
  # add the data
  geom_point(aes(x=theta, y=count), size=3) + 
  geom_segment(aes(x=currentDirection, y=0, xend=currentDirection, yend=r*10), size = 1, data=circStatsCurrdnf) + 
  scale_color_manual(values=c("gray", "black")) + 
  facet_wrap(~tide) +
  # groom the plot
  theme_bw() + 
  theme(strip.background = element_blank(), strip.text.x = element_blank()) +
  theme(panel.grid.minor = element_line(size = .2, color = "white")) + 
  theme(panel.grid.major = element_line(size=.5, color = "white")) +
  theme(axis.title.y = element_text(size=18, face="bold")) + 
  theme(axis.text.y = element_text(size=16)) +
  scale_y_continuous(breaks=c(0,3,6,10),labels=c(0,0.3,0.6,1)) +
  theme(axis.text.x = element_text(face="bold",size=22, vjust=.05)) + 
  scale_x_continuous("", limits=c(0,360), breaks=seq(0,360-1,by=45), labels=c("N","","","","S","","","")) 

pcurC


# lm for current speed
dmodel <- lm(drift_distance ~ tide+depth..m.+location, data=d2)
summary(dmodel)

pd <- ggplot(d2) +
  geom_boxplot(aes(tide,drift_distance,fill=location,linetype=depth..m.))
pd



dmodel <- lm(drift_distance ~ tide, data=d2)
summary(dmodel)


dnf <- d2[d2$depth..m.=="18" & d2$location=="near" & d2$tide=="flood",]
dne <- d2[d2$depth..m.=="18" & d2$location=="near" & d2$tide=="ebb",]

mean(dnf$drift_distance)
sd(dnf$drift_distance)
mean(dne$drift_distance)
sd(dne$drift_distance)


# current speed in deep / shallow waters

dns <- d2[d2$depth..m.=="9",]
dnd <- d2[d2$depth..m.=="18",]

dnsS <- mean(dns$current_strength)
dndS <- mean(dnd$current_strength)


## Current orientation through ontogeny --------------

dneT <- subset(d, depth..m.=="18" & location=="near" & tide=="ebb")

dneT$current_bearing <- dneT$mean - dneT$drift_direction
dneT = dneT[!is.na(dneT$current_bearing),] # remove trials with no current data

# correct for negative angle values produced
for (i in 1:nrow(dneT)){  
  if (dneT[i,c("current_bearing")]<=0){
    dneT[i,c("current_bearing")] <- 360 + dneT[i,c("current_bearing")]
  }
}


cltc <- lm.circular(y=dneT$current_bearing,x=dneT$fish_age..dph.,type="c-l",init=0)


dbTC <- ddply(dneT,.(tide), function(x, bin) {  # 'splitting' by method here b/c theres only 1 method and its easier than changing code to remove ddply
  x$mean <- as.numeric(round_any(x$current_bearing, bin))
  x$mean[x$mean==360] <- 0
  
  # create a data.frame with count
  counts <- count(x, "mean")
  
  # repeat each point the appropriate number of times
  d <- adply(counts, 1, function(x) {
    data.frame(theta=x$mean, count=1:x$freq)
  }, .expand=F)
  
  # make the scale prettier
  d$count <- 10 + d$count
  
  return(d)
}, bin=bin)

circStatsTCur <- circStatsC(dneT)
# significant finding! swim against current in the flood tide
for (i in 1:nrow(circStatsTCur)){  
  if (circStatsTCur[i,c("currentBearing")]<=0){
    circStatsTCur[i,c("currentBearing")] <- 360 + circStatsTCur[i,c("currentBearing")]
  }
}

pTC <- ggplot(dbTC) + polar() + # set the data and the polar layout
  # draw circles on the plot for the r-axis
  geom_line(data=circle, aes(theta,rho), color="gray") +
  geom_line(data=circle2, aes(theta,rho), color="gray") +
  # add the data
  geom_point(aes(x=theta, y=count), size=1) + # add the bearings as stacked dots
  geom_segment(aes(x=currentBearing, y=0, xend=currentBearing, yend=r*10, color=signif), size = .6, data=circStatsTCur) + # add the mean lines, w/ type determined by significance
  scale_color_manual(values=c("gray", "black")) + # set black line to be significant
  # groom the plot a bit
  scale_x_continuous("", limits=c(0,360), breaks=seq(0,360-1,by=45), labels=c("With Current","","","","Against Current","","","")) + # set the X axis (angle axis)
  theme_bw() + labs(color="Statistically Significant") + # remove the gray background and adjust legend label
  theme(legend.title = element_text(size=8)) + theme(legend.text = element_text(size=6)) + # adjust legend text sizes
  theme(panel.grid.minor = element_line(size = .2, color = "white")) + theme(panel.grid.major = element_line(size=.5, color = "white")) + # remove x-axis grid, make y-axis grid black
  theme(plot.title = element_text(face="bold", size=18)) + theme(axis.title.y = element_text(size=8, face="bold")) + # format the plot title and r-axis label
  theme(axis.text.x = element_text(size=8, vjust=.05)) + theme(axis.text.y = element_blank(), axis.ticks = element_blank()) + # adjust axis text
  theme(strip.text.x = element_text(size=6, face="bold"), strip.background = element_rect(colour="black", fill="#CCCCFF"))  # format facet labels

pTC


# now try with specific age gropus

circStatsAgeC <- ddply(dneT,.(fish_age..dph.),circStatsC) # run the statistics
circStatsAgeC$age <- circStatsAgeC$fish_age..dph.
for (i in 1:nrow(circStatsAgeC)){  
  if (circStatsAgeC[i,c("currentBearing")]<=0){
    circStatsAgeC[i,c("currentBearing")] <- 360 + circStatsAgeC[i,c("currentBearing")]
  }
}

# bin angles for plot
t <- dneT
t$theta <- as.numeric(round_any(t$current_bearing, bin))
t$theta[t$theta==360] <- 0
# create a data.frame with count
counts <- count(t, vars=c("theta","fish_age..dph."))

# repeat each point the appropriate number of times
tbC <- adply(counts, 1, function(x) {
  data.frame(theta=x$theta,age=x$fish_age..dph., count=1:x$freq)
}, .expand=F)

# make the scale prettier
tbC$count <- 10 + tbC$count



poaC <- ggplot(tbC) + polar() + # set the data and the polar layout
  # draw circles on the plot for the r-axis
  geom_line(data=circle, aes(theta,rho), color="gray") +
  geom_line(data=circle2, aes(theta,rho), color="gray") +
  # add the data
  geom_point(aes(x=theta, y=count), size=1) + # add the bearings as stacked dots
  geom_segment(aes(x=currentBearing, y=0, xend=currentBearing, yend=r*10, color=signif), size = .6, data=circStatsAgeC) + # add the mean lines, w/ type determined by significance
  scale_color_manual(values=c("gray", "black")) + # set black line to be significant
  facet_wrap(~age) +
  # groom the plot a bit
  scale_x_continuous("", limits=c(0,360), breaks=seq(0,360-1,by=45), labels=c("With Current","","","","Against Current","","","")) + # set the X axis (angle axis)
  theme_bw() + labs(color="Statistically Significant") + # remove the gray background and adjust legend label
  theme(legend.title = element_text(size=8)) + theme(legend.text = element_text(size=6)) + # adjust legend text sizes
  theme(panel.grid.minor = element_line(size = .2, color = "white")) + theme(panel.grid.major = element_line(size=.5, color = "white")) + # remove x-axis grid, make y-axis grid black
  theme(plot.title = element_text(face="bold", size=18)) + theme(axis.title.y = element_text(size=8, face="bold")) + # format the plot title and r-axis label
  theme(axis.text.x = element_text(size=8, vjust=.05)) + theme(axis.text.y = element_blank(), axis.ticks = element_blank()) + # adjust axis text
  theme(strip.text.x = element_text(size=6, face="bold"), strip.background = element_rect(colour="black", fill="#CCCCFF"))  # format facet labels

poaC


# try larger groups
tg <- dneT
tg$ageG <- 0
tg[tg$fish_age..dph. <= 10,c("ageG")] <- "10&Under"
tg[tg$fish_age..dph. > 20,c("ageG")] <- "22&Over"
tg[tg$ageG == 0, c("ageG")] <- "12to20"
tg$ageG <- as.factor(tg$ageG)


circStatsAgeGC <- ddply(tg,.(ageG),circStatsC) # run the statistics
circStatsAgeGC$age <- circStatsAgeGC$ageG
for (i in 1:nrow(circStatsAgeGC)){  
  if (circStatsAgeGC[i,c("currentBearing")]<=0){
    circStatsAgeGC[i,c("currentBearing")] <- 360 + circStatsAgeGC[i,c("currentBearing")]
  }
}

# bin angles for plot
tg$theta <- as.numeric(round_any(tg$current_bearing, bin))
tg$theta[tg$theta==360] <- 0
# create a data.frame with count
counts <- count(tg, vars=c("theta","ageG"))

# repeat each point the appropriate number of times
tbg <- adply(counts, 1, function(x) {
  data.frame(theta=x$theta,age=x$ageG, count=1:x$freq)
}, .expand=F)

# make the scale prettier
tbg$count <- 10 + tbg$count



pgaC <- ggplot(tbg) + polar() + # set the data and the polar layout
  # draw circles on the plot for the r-axis
  geom_line(data=circle, aes(theta,rho), color="gray") +
  geom_line(data=circle2, aes(theta,rho), color="gray") +
  # add the data
  geom_point(aes(x=theta, y=count), size=1) + # add the bearings as stacked dots
  geom_segment(aes(x=currentBearing, y=0, xend=currentBearing, yend=r*10, color=signif), size = .6, data=circStatsAgeGC) + # add the mean lines, w/ type determined by significance
  scale_color_manual(values=c("gray", "black")) + # set black line to be significant
  facet_wrap(~age) +
  # groom the plot a bit
  scale_x_continuous("", limits=c(0,360), breaks=seq(0,360-1,by=45), labels=c("With Current","","","","Against Current","","","")) + # set the X axis (angle axis)
  theme_bw() + labs(color="Statistically Significant") + # remove the gray background and adjust legend label
  theme(legend.title = element_text(size=8)) + theme(legend.text = element_text(size=6)) + # adjust legend text sizes
  theme(panel.grid.minor = element_line(size = .2, color = "white")) + theme(panel.grid.major = element_line(size=.5, color = "white")) + # remove x-axis grid, make y-axis grid black
  theme(plot.title = element_text(face="bold", size=18)) + theme(axis.title.y = element_text(size=8, face="bold")) + # format the plot title and r-axis label
  theme(axis.text.x = element_text(size=5, vjust=.05)) + theme(axis.text.y = element_blank(), axis.ticks = element_blank()) + # adjust axis text
  theme(strip.text.x = element_text(size=6, face="bold"), strip.background = element_rect(colour="black", fill="#CCCCFF"))  # format facet labels

pgaC








## Orientation to the wind direction ----------- 

#SUPMAT PLOT

d2 = d2[!is.na(d2$wind_dir),] # remove trials with no wind data
# correct angles so headings are relative to current direction
d2$wind_bearing <- d2$mean - as.numeric(d2$wind_dir)

# correct for negative angle values produced
for (i in 1:nrow(d2)){  
  if (d2[i,c("wind_bearing")]<=0){
    d2[i,c("wind_bearing")] <- 360 + d2[i,c("wind_bearing")]
  }
}

# define function to look at current bearings
circStatsW <- function(x) {
  n = nrow(x)
  windBearing = mean.circular(x$wind_bearing)
  sdbearing = sd.circular(x$wind_bearing)
  meanR = mean(x$r, na.rm=T)
  sdR = sd(x$r, na.rm=T)
  ray = rayleigh.test(na.omit(x$wind_bearing))
  r = ray$statistic
  p = ray$p.value
  signif = p < 0.05
  return(data.frame(n, windBearing,sdbearing, meanR, sdR, r, p, signif))
}



dbW <- ddply(d2,.(depth..m.,tide,location), function(x, bin) {  # 'splitting' by method here b/c theres only 1 method and its easier than changing code to remove ddply
  x$mean <- as.numeric(round_any(x$wind_bearing, bin))
  x$mean[x$mean==360] <- 0
  
  # create a data.frame with count
  counts <- count(x, "mean")
  
  # repeat each point the appropriate number of times
  d <- adply(counts, 1, function(x) {
    data.frame(theta=x$mean, count=1:x$freq)
  }, .expand=F)
  
  # make the scale prettier
  d$count <- 10 + d$count
  
  return(d)
}, bin=bin)

circStatsWind <- ddply(d2,.(depth..m.,tide,location),circStatsW) # run the statistics
circStatsWind[circStatsWind$windBearing <= 0,c("windBearing")] <- circStatsWind[circStatsWind$windBearing <= 0,c("windBearing")] + 360

pW <- ggplot(dbW) + polar() + # set the data and the polar layout
  # draw circles on the plot for the r-axis
  geom_line(data=circle, aes(theta,rho), color="gray") +
  geom_line(data=circle2, aes(theta,rho), color="gray") +
  # add the data
  geom_point(aes(x=theta, y=count), size=1) + # add the bearings as stacked dots
  geom_segment(aes(x=windBearing, y=0, xend=windBearing, yend=r*10, color=signif), size = .6, data=circStatsWind) + # add the mean lines, w/ type determined by significance
  scale_color_manual(values=c("gray", "black")) + # set black line to be significant
  facet_wrap(~tide+location+depth..m.) +
  # groom the plot a bit
  theme_bw() + ylab("Rho-Value") +
  theme(panel.grid.minor = element_line(size = .2, color = "white")) + 
  theme(panel.grid.major = element_line(size=.5, color = "white")) +
  theme(axis.title.y = element_text(size=16, face="bold")) + 
  theme(axis.text.y = element_text(size=16)) +
  scale_y_continuous(breaks=c(0,3,6,10),labels=c(0,0.3,0.6,1)) +
  theme(axis.text.x = element_text(face="bold",size=12, vjust=.05)) + 
  scale_x_continuous("", limits=c(0,360), breaks=seq(0,360-1,by=45), labels=c("With Wind","","","","Against Wind","","","")) + 
  theme(strip.background = element_blank(),strip.text.x = element_blank())

pW






## Wind directions by environmental condition --------------


# define function to look at current bearings
circStatsW2 <- function(x) {
  n = nrow(x)
  windDir = mean.circular(x$wind_dir)
  sdbearing = sd.circular(x$wind_dir)
  meanR = mean(x$r, na.rm=T)
  sdR = sd(x$r, na.rm=T)
  ray = rayleigh.test(na.omit(x$wind_dir))
  r = ray$statistic
  p = ray$p.value
  signif = p < 0.05
  return(data.frame(n, windDir,sdbearing, meanR, sdR, r, p, signif))
}



dbW2 <- ddply(d2,.(depth..m.,tide,location), function(x, bin) {  # 'splitting' by method here b/c theres only 1 method and its easier than changing code to remove ddply
  x$mean <- as.numeric(round_any(x$wind_dir, bin))
  x$mean[x$mean==360] <- 0
  
  # create a data.frame with count
  counts <- count(x, "mean")
  
  # repeat each point the appropriate number of times
  d <- adply(counts, 1, function(x) {
    data.frame(theta=x$mean, count=1:x$freq)
  }, .expand=F)
  
  # make the scale prettier
  d$count <- 10 + d$count
  
  return(d)
}, bin=bin)

circStatsWindDir <- ddply(d2,.(depth..m.,tide,location),circStatsW2) # run the statistics


pW2 <- ggplot(dbW2) + polar() + # set the data and the polar layout
  # draw circles on the plot for the r-axis
  geom_line(data=circle, aes(theta,rho), color="gray") +
  geom_line(data=circle2, aes(theta,rho), color="gray") +
  # add the data
  geom_point(aes(x=theta, y=count), size=1) + # add the bearings as stacked dots
  geom_segment(aes(x=windDir, y=0, xend=windDir, yend=r*10, color=signif), size = .6, data=circStatsWindDir) + # add the mean lines, w/ type determined by significance
  scale_color_manual(values=c("gray", "black")) + # set black line to be significant
  facet_wrap(~tide+location+depth..m.) +
  # groom the plot a bit
  scale_x_continuous("", limits=c(0,360), breaks=seq(0,360-1,by=45), labels=c("N","","","","S","","","")) + # set the X axis (angle axis)
  theme_bw() + labs(color="Statistically Significant") + # remove the gray background and adjust legend label
  theme(legend.title = element_text(size=8)) + theme(legend.text = element_text(size=6)) + # adjust legend text sizes
  theme(panel.grid.minor = element_line(size = .2, color = "white")) + theme(panel.grid.major = element_line(size=.5, color = "white")) + # remove x-axis grid, make y-axis grid black
  theme(plot.title = element_text(face="bold", size=18)) + theme(axis.title.y = element_text(size=8, face="bold")) + # format the plot title and r-axis label
  theme(axis.text.x = element_text(size=5, vjust=.05)) + theme(axis.text.y = element_blank(), axis.ticks = element_blank()) + # adjust axis text
  theme(strip.text.x = element_text(size=6, face="bold"), strip.background = element_rect(colour="black", fill="#CCCCFF"))  # format facet labels

pW2





## Wind orientation through ontogeny ---------


dneT <- subset(d, depth..m.=="18" & location=="near" & tide=="ebb")

dneT$wind_bearing <- dneT$mean - dneT$wind_dir
dneT = dneT[!is.na(dneT$wind_bearing),] # remove trials with no wind data

# correct for negative angle values produced
for (i in 1:nrow(dneT)){  
  if (dneT[i,c("wind_bearing")]<=0){
    dneT[i,c("wind_bearing")] <- 360 + dneT[i,c("wind_bearing")]
  }
}


cltw <- lm.circular(y=dneT$wind_bearing,x=dneT$fish_age..dph.,type="c-l",init=0)


dbTW <- ddply(dneT,.(tide), function(x, bin) {  # 'splitting' by method here b/c theres only 1 method and its easier than changing code to remove ddply
  x$mean <- as.numeric(round_any(x$wind_bearing, bin))
  x$mean[x$mean==360] <- 0
  
  # create a data.frame with count
  counts <- count(x, "mean")
  
  # repeat each point the appropriate number of times
  d <- adply(counts, 1, function(x) {
    data.frame(theta=x$mean, count=1:x$freq)
  }, .expand=F)
  
  # make the scale prettier
  d$count <- 10 + d$count
  
  return(d)
}, bin=bin)

circStatsTW <- circStatsW(dneT)

for (i in 1:nrow(circStatsTW)){  
  if (circStatsTW[i,c("windBearing")]<=0){
    circStatsTW[i,c("windBearing")] <- 360 + circStatsTW[i,c("windBearing")]
  }
}

pTW <- ggplot(dbTW) + polar() + # set the data and the polar layout
  # draw circles on the plot for the r-axis
  geom_line(data=circle, aes(theta,rho), color="gray") +
  geom_line(data=circle2, aes(theta,rho), color="gray") +
  # add the data
  geom_point(aes(x=theta, y=count), size=1) + # add the bearings as stacked dots
  geom_segment(aes(x=windBearing, y=0, xend=windBearing, yend=r*10, color=signif), size = .6, data=circStatsTW) + # add the mean lines, w/ type determined by significance
  scale_color_manual(values=c("gray", "black")) + # set black line to be significant
  # groom the plot a bit
  scale_x_continuous("", limits=c(0,360), breaks=seq(0,360-1,by=45), labels=c("Against wind","","","","With wind","","","")) + # set the X axis (angle axis)
  theme_bw() + labs(color="Statistically Significant") + # remove the gray background and adjust legend label
  theme(legend.title = element_text(size=8)) + theme(legend.text = element_text(size=6)) + # adjust legend text sizes
  theme(panel.grid.minor = element_line(size = .2, color = "white")) + theme(panel.grid.major = element_line(size=.5, color = "white")) + # remove x-axis grid, make y-axis grid black
  theme(plot.title = element_text(face="bold", size=18)) + theme(axis.title.y = element_text(size=8, face="bold")) + # format the plot title and r-axis label
  theme(axis.text.x = element_text(size=8, vjust=.05)) + theme(axis.text.y = element_blank(), axis.ticks = element_blank()) + # adjust axis text
  theme(strip.text.x = element_text(size=6, face="bold"), strip.background = element_rect(colour="black", fill="#CCCCFF"))  # format facet labels

pTW


# now try with specific age gropus

circStatsAgeW <- ddply(dneT,.(fish_age..dph.),circStatsW) # run the statistics
circStatsAgeW$age <- circStatsAgeW$fish_age..dph.
for (i in 1:nrow(circStatsAgeW)){  
  if (circStatsAgeW[i,c("windBearing")]<=0){
    circStatsAgeW[i,c("windBearing")] <- 360 + circStatsAgeW[i,c("windBearing")]
  }
}

# bin angles for plot
t <- dneT
t$theta <- as.numeric(round_any(t$wind_bearing, bin))
t$theta[t$theta==360] <- 0
# create a data.frame with count
counts <- count(t, vars=c("theta","fish_age..dph."))

# repeat each point the appropriate number of times
tbW <- adply(counts, 1, function(x) {
  data.frame(theta=x$theta,age=x$fish_age..dph., count=1:x$freq)
}, .expand=F)

# make the scale prettier
tbW$count <- 10 + tbW$count



poaW <- ggplot(tbW) + polar() + # set the data and the polar layout
  # draw circles on the plot for the r-axis
  geom_line(data=circle, aes(theta,rho), color="gray") +
  geom_line(data=circle2, aes(theta,rho), color="gray") +
  # add the data
  geom_point(aes(x=theta, y=count), size=1) + # add the bearings as stacked dots
  geom_segment(aes(x=windBearing, y=0, xend=windBearing, yend=r*10, color=signif), size = .6, data=circStatsAgeW) + # add the mean lines, w/ type determined by significance
  scale_color_manual(values=c("gray", "black")) + # set black line to be significant
  facet_wrap(~age) +
  # groom the plot a bit
  scale_x_continuous("", limits=c(0,360), breaks=seq(0,360-1,by=45), labels=c("Against wind","","","","With wind","","","")) + # set the X axis (angle axis)
  theme_bw() + labs(color="Statistically Significant") + # remove the gray background and adjust legend label
  theme(legend.title = element_text(size=8)) + theme(legend.text = element_text(size=6)) + # adjust legend text sizes
  theme(panel.grid.minor = element_line(size = .2, color = "white")) + theme(panel.grid.major = element_line(size=.5, color = "white")) + # remove x-axis grid, make y-axis grid black
  theme(plot.title = element_text(face="bold", size=18)) + theme(axis.title.y = element_text(size=8, face="bold")) + # format the plot title and r-axis label
  theme(axis.text.x = element_text(size=8, vjust=.05)) + theme(axis.text.y = element_blank(), axis.ticks = element_blank()) + # adjust axis text
  theme(strip.text.x = element_text(size=6, face="bold"), strip.background = element_rect(colour="black", fill="#CCCCFF"))  # format facet labels

poaW


# try larger groups
tg <- dneT
tg$ageG <- 0
tg[tg$fish_age..dph. <= 10,c("ageG")] <- "10&Under"
tg[tg$fish_age..dph. > 20,c("ageG")] <- "22&Over"
tg[tg$ageG == 0, c("ageG")] <- "12to20"
tg$ageG <- as.factor(tg$ageG)


circStatsAgeGW <- ddply(tg,.(ageG),circStatsW) # run the statistics
circStatsAgeGW$age <- circStatsAgeGW$ageG
for (i in 1:nrow(circStatsAgeGC)){  
  if (circStatsAgeGW[i,c("windBearing")]<=0){
    circStatsAgeGW[i,c("windBearing")] <- 360 + circStatsAgeGW[i,c("windBearing")]
  }
}

# bin angles for plot
tg$theta <- as.numeric(round_any(tg$wind_bearing, bin))
tg$theta[tg$theta==360] <- 0
# create a data.frame with count
counts <- count(tg, vars=c("theta","ageG"))

# repeat each point the appropriate number of times
tbg <- adply(counts, 1, function(x) {
  data.frame(theta=x$theta,age=x$ageG, count=1:x$freq)
}, .expand=F)

# make the scale prettier
tbg$count <- 10 + tbg$count



pgaW <- ggplot(tbg) + polar() + # set the data and the polar layout
  # draw circles on the plot for the r-axis
  geom_line(data=circle, aes(theta,rho), color="gray") +
  geom_line(data=circle2, aes(theta,rho), color="gray") +
  # add the data
  geom_point(aes(x=theta, y=count), size=1) + # add the bearings as stacked dots
  geom_segment(aes(x=windBearing, y=0, xend=windBearing, yend=r*10, color=signif), size = .6, data=circStatsAgeGW) + # add the mean lines, w/ type determined by significance
  scale_color_manual(values=c("gray", "black")) + # set black line to be significant
  facet_wrap(~age) +
  # groom the plot a bit
  scale_x_continuous("", limits=c(0,360), breaks=seq(0,360-1,by=45), labels=c("Against wind","","","","With wind","","","")) + # set the X axis (angle axis)
  theme_bw() + labs(color="Statistically Significant") + # remove the gray background and adjust legend label
  theme(legend.title = element_text(size=8)) + theme(legend.text = element_text(size=6)) + # adjust legend text sizes
  theme(panel.grid.minor = element_line(size = .2, color = "white")) + theme(panel.grid.major = element_line(size=.5, color = "white")) + # remove x-axis grid, make y-axis grid black
  theme(plot.title = element_text(face="bold", size=18)) + theme(axis.title.y = element_text(size=8, face="bold")) + # format the plot title and r-axis label
  theme(axis.text.x = element_text(size=5, vjust=.05)) + theme(axis.text.y = element_blank(), axis.ticks = element_blank()) + # adjust axis text
  theme(strip.text.x = element_text(size=6, face="bold"), strip.background = element_rect(colour="black", fill="#CCCCFF"))  # format facet labels

pgaW





## Precision of orientation throughout the day -----------

d$roughlat <- 16.815663
d$roughlon <- -88.082242


# sun.position requires time in GMT. Belize is 5 hours back from GMT)
d$dateTimeGMT = d$dateTime + 6*3600
sunPos = sun.position(d$dateTimeGMT, d$roughlat, d$roughlon)
# add sun azimuth/zenith to the data
d = cbind(d, sunPos)
d$azimuth = circular(d$azimuth, units="degrees", template="geographics", modulo="2pi")
d$zenith = circular(d$zenith, units="degrees", template="geographics", modulo="2pi")


zenithRplot <- ggplot(d) + 
  # add the deployments 
  geom_point(aes(x=zenith, y=r), size=2) +
  # add a smooth line and bounds 
  geom_smooth(aes(x=zenith, y=r), method=loess, colour="black") +
  # make it pretty
  theme_bw() + ylab("Rayleigh's R Value") + xlab("Zenith Angle (degrees)") +
  theme(panel.grid.minor = element_line(size = .2, color = "white")) + theme(panel.grid.major = element_line(size=.5, color = "white")) +
  theme(axis.text.x = element_text(size=12), axis.title.x = element_text(size=12, face="bold",vjust=.01)) +
  theme(axis.text.y = element_text(size=12), axis.title.y = element_text(size=12, face="bold",vjust=.3)) # +
#   facet_grid(.~year) # optional, to split view by year, also uncomment the "+" sign in line above ^
zenithRplot


d$time <- as.POSIXct(str_c("2016-01-01"," ",d$time_start,sep=""))

timePlot <- ggplot(d) + 
  # add the deployments 
  geom_point(aes(x=time, y=r), size=2) +
  # add a smooth line and bounds 
  geom_smooth(aes(x=time, y=r), method=lm, colour="black") +
  geom_smooth(aes(x=time, y=r), method=loess, colour="blue") +
  # make it pretty
  theme_bw() + ylab("Rayleigh's R Value") + xlab("Time of Day") +
  theme(panel.grid.minor = element_line(size = .2, color = "white")) + theme(panel.grid.major = element_line(size=.5, color = "white")) +
  theme(axis.text.x = element_text(size=12), axis.title.x = element_text(size=12, face="bold",vjust=.01)) +
  theme(axis.text.y = element_text(size=12), axis.title.y = element_text(size=12, face="bold",vjust=.3)) # +
#   facet_grid(.~year) # optional, to split view by year, also uncomment the "+" sign in line above ^
timePlot


summary(lm(r ~ poly(time,2), data=d))
# significant finding - r is related to the time of the day

## Sun Compass -------------- SUPMAT PLOT

circStatsAz <- function(x) {
  n = nrow(x)
  azBearing = mean.circular(x$azmean)
  sdbearing = sd.circular(x$azmean)
  meanR = mean(x$r, na.rm=T)
  sdR = sd(x$r, na.rm=T)
  ray = rayleigh.test(na.omit(x$azmean))
  r = ray$statistic
  p = ray$p.value
  signif = p < 0.05
  return(data.frame(n, azBearing,sdbearing, meanR, sdR, r, p, signif))
}


d$azmean <- d$mean - d$azimuth
for (i in 1:nrow(d)){
  if (d[i,c("azmean")]<=0){
    d[i,c("azmean")] <- 360 + d[i,c("azmean")]
  }
}


circStatsSunA <- circStatsAz (d)



d2 <- subset(d, fish_age..dph. >17 & fish_age..dph. <25)

d2$azmean <- d2$mean - d2$azimuth
for (i in 1:nrow(d2)){
  if (d2[i,c("azmean")]<=0){
    d2[i,c("azmean")] <- 360 + d2[i,c("azmean")]
  }
}

circStatsSun <- ddply(d2, .(tide,depth..m.,location), circStatsAz)
for (i in 1:nrow(circStatsSun)){
  if (circStatsSun[i,c("azBearing")]<=0){
    circStatsSun[i,c("azBearing")] <- 360 + circStatsSun[i,c("azBearing")]
  }
}


# bin angles for each treatment 
BinnedDataAz <- ddply(d2, .(tide,depth..m.,location), function(x, bin) {
  x$azmean <- as.numeric(round_any(x$azmean, bin))
  x$azmean[x$azmean==360] <- 0
  
  # create a data.frame with count
  counts <- count(x, "azmean")
  
  # repeat each point the appropriate number of times to create the histogram
  d <- adply(counts, 1, function(x) {
    data.frame(theta=x$azmean, count=1:x$freq)
  }, .expand=F)
  
  # make the scale prettier
  d$count <- 10 + d$count
  
  return(d)
}, bin=bin)



pAz <- ggplot(BinnedDataAz) + polar() + # set the data and the polar layout
  # draw circles on the plot for the r-axis
  geom_line(data=circle, aes(theta,rho), color="gray") +
  geom_line(data=circle2, aes(theta,rho), color="gray") +
  # add the data
  geom_point(aes(x=theta, y=count), size=1) + # add the bearings as stacked dots
  geom_segment(aes(x=azBearing, y=0, xend=azBearing, yend=r*10, color=signif), size = .6, data=circStatsSun) + # add the mean lines, w/ type determined by significance
  scale_color_manual(values=c("gray", "black")) + # set black line to be significant
  facet_wrap(~tide+location+depth..m.) + # split the plot by the treatment variable
  # groom the plot a bit
  theme_bw() + ylab("Rho-Value") +
  theme(panel.grid.minor = element_line(size = .2, color = "white")) + 
  theme(panel.grid.major = element_line(size=.5, color = "white")) +
  theme(axis.title.y = element_text(size=16, face="bold")) + 
  theme(axis.text.y = element_text(size=16)) +
  scale_y_continuous(breaks=c(0,3,6,10),labels=c(0,0.3,0.6,1)) +
  theme(axis.text.x = element_text(face="bold",size=12, vjust=.05)) + 
  scale_x_continuous("", limits=c(0,360), breaks=seq(0,360-1,by=45), labels=c("Towards Azimuth","","","","Away from Azimuth","","","")) + 
  theme(strip.background = element_blank(),strip.text.x = element_blank())

pAz




# sun compass (type 2)

d$azmean2 <- d$azimuth - d$mean
for (i in 1:nrow(d)){
  if (d[i,c("azmean2")]<=0){
    d[i,c("azmean2")] <- 360 + d[i,c("azmean2")]
  }
}

timePlot2 <- ggplot(d) + 
  # add the deployments 
  geom_point(aes(x=time, y=as.numeric(azmean2)), size=2) +
  # add a smooth line and bounds 
  geom_smooth(aes(x=time, y=as.numeric(azmean2)), method=loess, colour="black") +
  # make it pretty
  theme_bw() + ylab("Azimuth - Bearing") + xlab("Time of Day") +
  theme(panel.grid.minor = element_line(size = .2, color = "white")) + theme(panel.grid.major = element_line(size=.5, color = "white")) +
  theme(axis.text.x = element_text(size=12), axis.title.x = element_text(size=12, face="bold",vjust=.01)) +
  theme(axis.text.y = element_text(size=12), axis.title.y = element_text(size=12, face="bold",vjust=.3)) # +
#   facet_grid(.~year) # optional, to split view by year, also uncomment the "+" sign in line above ^
timePlot2


# type 2 for BEST orienters and ONLY DNE group
# make sure you've selected best orienters only first

dne <- subset(d, depth..m.=="18" & location=="near" & tide=="ebb")

timePlot3 <- ggplot(dne) + 
  # add the deployments 
  geom_point(aes(x=time, y=as.numeric(azmean2)), size=2) +
  # add a smooth line and bounds 
  geom_smooth(aes(x=time, y=as.numeric(azmean2)), method=loess, colour="black") +
  # make it pretty
  theme_bw() + ylab("Azimuth - Bearing") + xlab("Time of Day") +
  theme(panel.grid.minor = element_line(size = .2, color = "white")) + theme(panel.grid.major = element_line(size=.5, color = "white")) +
  theme(axis.text.x = element_text(size=12), axis.title.x = element_text(size=12, face="bold",vjust=.01)) +
  theme(axis.text.y = element_text(size=12), axis.title.y = element_text(size=12, face="bold",vjust=.3)) # +
#   facet_grid(.~year) # optional, to split view by year, also uncomment the "+" sign in line above ^
timePlot3





## Sun Compass orientation through ontogeny ---------


dneT <- subset(d, depth..m.=="18" & location=="near" & tide=="ebb")

dneT$azmean <- dneT$azimuth - dneT$mean
dneT = dneT[!is.na(dneT$azmean),] # remove trials with no wind data

# correct for negative angle values produced
for (i in 1:nrow(dneT)){  
  if (dneT[i,c("azmean")]<=0){
    dneT[i,c("azmean")] <- 360 + dneT[i,c("azmean")]
  }
}


clta <- lm.circular(y=dneT$azmean,x=dneT$fish_age..dph.,type="c-l",init=0)


dbTA <- ddply(dneT,.(tide), function(x, bin) {  # 'splitting' by method here b/c theres only 1 method and its easier than changing code to remove ddply
  x$mean <- as.numeric(round_any(x$azmean, bin))
  x$mean[x$mean==360] <- 0
  
  # create a data.frame with count
  counts <- count(x, "mean")
  
  # repeat each point the appropriate number of times
  d <- adply(counts, 1, function(x) {
    data.frame(theta=x$mean, count=1:x$freq)
  }, .expand=F)
  
  # make the scale prettier
  d$count <- 10 + d$count
  
  return(d)
}, bin=bin)

circStatsTA <- circStatsAz(dneT)

for (i in 1:nrow(circStatsTA)){  
  if (circStatsTA[i,c("azBearing")]<=0){
    circStatsTA[i,c("azBearing")] <- 360 + circStatsTA[i,c("azBearing")]
  }
}

pTA <- ggplot(dbTA) + polar() + # set the data and the polar layout
  # draw circles on the plot for the r-axis
  geom_line(data=circle, aes(theta,rho), color="gray") +
  geom_line(data=circle2, aes(theta,rho), color="gray") +
  # add the data
  geom_point(aes(x=theta, y=count), size=1) + # add the bearings as stacked dots
  geom_segment(aes(x=azBearing, y=0, xend=azBearing, yend=r*10, color=signif), size = .6, data=circStatsTA) + # add the mean lines, w/ type determined by significance
  scale_color_manual(values=c("gray", "black")) + # set black line to be significant
  # groom the plot a bit
  scale_x_continuous("", limits=c(0,360), breaks=seq(0,360-1,by=45), labels=c("Towards Azimuth","","","","Away From Azimuth","","","")) + # set the X axis (angle axis)
  theme_bw() + labs(color="Statistically Significant") + # remove the gray background and adjust legend label
  theme(legend.title = element_text(size=8)) + theme(legend.text = element_text(size=6)) + # adjust legend text sizes
  theme(panel.grid.minor = element_line(size = .2, color = "white")) + theme(panel.grid.major = element_line(size=.5, color = "white")) + # remove x-axis grid, make y-axis grid black
  theme(plot.title = element_text(face="bold", size=18)) + theme(axis.title.y = element_text(size=8, face="bold")) + # format the plot title and r-axis label
  theme(axis.text.x = element_text(size=8, vjust=.05)) + theme(axis.text.y = element_blank(), axis.ticks = element_blank()) + # adjust axis text
  theme(strip.text.x = element_text(size=6, face="bold"), strip.background = element_rect(colour="black", fill="#CCCCFF"))  # format facet labels

pTA


# now try with specific age gropus

circStatsAgeAz <- ddply(dneT,.(fish_age..dph.),circStatsAz) # run the statistics
circStatsAgeAz$age <- circStatsAgeAz$fish_age..dph.
for (i in 1:nrow(circStatsAgeAz)){  
  if (circStatsAgeAz[i,c("azBearing")]<=0){
    circStatsAgeAz[i,c("azBearing")] <- 360 + circStatsAgeAz[i,c("azBearing")]
  }
}

# bin angles for plot
t <- dneT
t$theta <- as.numeric(round_any(t$azmean, bin))
t$theta[t$theta==360] <- 0
# create a data.frame with count
counts <- count(t, vars=c("theta","fish_age..dph."))

# repeat each point the appropriate number of times
tbAz <- adply(counts, 1, function(x) {
  data.frame(theta=x$theta,age=x$fish_age..dph., count=1:x$freq)
}, .expand=F)

# make the scale prettier
tbAz$count <- 10 + tbAz$count



poaAz <- ggplot(tbAz) + polar() + # set the data and the polar layout
  # draw circles on the plot for the r-axis
  geom_line(data=circle, aes(theta,rho), color="gray") +
  geom_line(data=circle2, aes(theta,rho), color="gray") +
  # add the data
  geom_point(aes(x=theta, y=count), size=1) + # add the bearings as stacked dots
  geom_segment(aes(x=azBearing, y=0, xend=azBearing, yend=r*10, color=signif), size = .6, data=circStatsAgeAz) + # add the mean lines, w/ type determined by significance
  scale_color_manual(values=c("gray", "black")) + # set black line to be significant
  facet_wrap(~age) +
  # groom the plot a bit
  scale_x_continuous("", limits=c(0,360), breaks=seq(0,360-1,by=45), labels=c("Away from Azimuth","","","","Towards Azimuth","","","")) + # set the X axis (angle axis)
  theme_bw() + labs(color="Statistically Significant") + # remove the gray background and adjust legend label
  theme(legend.title = element_text(size=8)) + theme(legend.text = element_text(size=6)) + # adjust legend text sizes
  theme(panel.grid.minor = element_line(size = .2, color = "white")) + theme(panel.grid.major = element_line(size=.5, color = "white")) + # remove x-axis grid, make y-axis grid black
  theme(plot.title = element_text(face="bold", size=18)) + theme(axis.title.y = element_text(size=8, face="bold")) + # format the plot title and r-axis label
  theme(axis.text.x = element_text(size=8, vjust=.05)) + theme(axis.text.y = element_blank(), axis.ticks = element_blank()) + # adjust axis text
  theme(strip.text.x = element_text(size=6, face="bold"), strip.background = element_rect(colour="black", fill="#CCCCFF"))  # format facet labels

poaAz


# try larger groups
tg <- dneT
tg$ageG <- 0
tg[tg$fish_age..dph. <= 10,c("ageG")] <- "10&Under"
tg[tg$fish_age..dph. > 20,c("ageG")] <- "22&Over"
tg[tg$ageG == 0, c("ageG")] <- "12to20"
tg$ageG <- as.factor(tg$ageG)


circStatsAgeGAz <- ddply(tg,.(ageG),circStatsAz) # run the statistics
circStatsAgeGAz$age <- circStatsAgeGAz$ageG
for (i in 1:nrow(circStatsAgeGAz)){  
  if (circStatsAgeGAz[i,c("azBearing")]<=0){
    circStatsAgeGAz[i,c("azBearing")] <- 360 + circStatsAgeGAz[i,c("azBearing")]
  }
}

# bin angles for plot
tg$theta <- as.numeric(round_any(tg$azmean, bin))
tg$theta[tg$theta==360] <- 0
# create a data.frame with count
counts <- count(tg, vars=c("theta","ageG"))

# repeat each point the appropriate number of times
tbg <- adply(counts, 1, function(x) {
  data.frame(theta=x$theta,age=x$ageG, count=1:x$freq)
}, .expand=F)

# make the scale prettier
tbg$count <- 10 + tbg$count



pgaAz <- ggplot(tbg) + polar() + # set the data and the polar layout
  # draw circles on the plot for the r-axis
  geom_line(data=circle, aes(theta,rho), color="gray") +
  geom_line(data=circle2, aes(theta,rho), color="gray") +
  # add the data
  geom_point(aes(x=theta, y=count), size=1) + # add the bearings as stacked dots
  geom_segment(aes(x=azBearing, y=0, xend=azBearing, yend=r*10, color=signif), size = .6, data=circStatsAgeGAz) + # add the mean lines, w/ type determined by significance
  scale_color_manual(values=c("gray", "black")) + # set black line to be significant
  facet_wrap(~age) +
  # groom the plot a bit
  scale_x_continuous("", limits=c(0,360), breaks=seq(0,360-1,by=45), labels=c("Away from Azimuth","","","","Towards Azimuth","","","")) + # set the X axis (angle axis)
  theme_bw() + labs(color="Statistically Significant") + # remove the gray background and adjust legend label
  theme(legend.title = element_text(size=8)) + theme(legend.text = element_text(size=6)) + # adjust legend text sizes
  theme(panel.grid.minor = element_line(size = .2, color = "white")) + theme(panel.grid.major = element_line(size=.5, color = "white")) + # remove x-axis grid, make y-axis grid black
  theme(plot.title = element_text(face="bold", size=18)) + theme(axis.title.y = element_text(size=8, face="bold")) + # format the plot title and r-axis label
  theme(axis.text.x = element_text(size=5, vjust=.05)) + theme(axis.text.y = element_blank(), axis.ticks = element_blank()) + # adjust axis text
  theme(strip.text.x = element_text(size=6, face="bold"), strip.background = element_rect(colour="black", fill="#CCCCFF"))  # format facet labels

pgaAz


# looking at bimodality
tg2 <- tg
tg2$azmean <- tg$azmean*2
tg2[tg2$azmean >= 360,c("azmean")] <- tg2[tg2$azmean >= 360,c("azmean")] - 360


circStatsAgeGAz2 <- ddply(tg2,.(ageG),circStatsAz) # run the statistics
circStatsAgeGAz2$age <- circStatsAgeGAz$ageG
for (i in 1:nrow(circStatsAgeGAz2)){  
  if (circStatsAgeGAz2[i,c("azBearing")]<=0){
    circStatsAgeGAz2[i,c("azBearing")] <- 360 + circStatsAgeGAz2[i,c("azBearing")]
  }
}

# bin angles for plot
tg2$theta <- as.numeric(round_any(tg2$azmean, bin))
tg2$theta[tg2$theta==360] <- 0
# create a data.frame with count
counts <- count(tg2, vars=c("theta","ageG"))

# repeat each point the appropriate number of times
tbg2 <- adply(counts, 1, function(x) {
  data.frame(theta=x$theta,age=x$ageG, count=1:x$freq)
}, .expand=F)

# make the scale prettier
tbg2$count <- 10 + tbg2$count



pgaAz2 <- ggplot(tbg2) + polar() + # set the data and the polar layout
  # draw circles on the plot for the r-axis
  geom_line(data=circle, aes(theta,rho), color="gray") +
  geom_line(data=circle2, aes(theta,rho), color="gray") +
  # add the data
  geom_point(aes(x=theta, y=count), size=1) + # add the bearings as stacked dots
  geom_segment(aes(x=azBearing, y=0, xend=azBearing, yend=r*10, color=signif), size = .6, data=circStatsAgeGAz2) + # add the mean lines, w/ type determined by significance
  scale_color_manual(values=c("gray", "black")) + # set black line to be significant
  facet_wrap(~age) +
  # groom the plot a bit
  scale_x_continuous("", limits=c(0,360), breaks=seq(0,360-1,by=45), labels=c("Away from Azimuth","","","","Towards Azimuth","","","")) + # set the X axis (angle axis)
  theme_bw() + labs(color="Statistically Significant") + # remove the gray background and adjust legend label
  theme(legend.title = element_text(size=8)) + theme(legend.text = element_text(size=6)) + # adjust legend text sizes
  theme(panel.grid.minor = element_line(size = .2, color = "white")) + theme(panel.grid.major = element_line(size=.5, color = "white")) + # remove x-axis grid, make y-axis grid black
  theme(plot.title = element_text(face="bold", size=18)) + theme(axis.title.y = element_text(size=8, face="bold")) + # format the plot title and r-axis label
  theme(axis.text.x = element_text(size=5, vjust=.05)) + theme(axis.text.y = element_blank(), axis.ticks = element_blank()) + # adjust axis text
  theme(strip.text.x = element_text(size=6, face="bold"), strip.background = element_rect(colour="black", fill="#CCCCFF"))  # format facet labels

pgaAz2




## Orientation and Swimming in Relation to Light/Temp ----------


#get rid of some weird r light values
dl <- subset(d, meanLight <= 2500)


# # Repeat logistic model on proportion of orienters (make sure you've not removed the FALSErs)
# 
# logitFitHobo <- glm(signif ~ tide+depth..m.+location+fish_age..dph.+sky+dateTime, 
#                     data=d, family="binomial")
# summary(logitFitHobo)


# Linear / Circular correlations with mean bearings

clt <- lm.circular(y=d$mean,x=d$meanTemp,type="c-l",init=0)
cll <- lm.circular(y=dl$mean,x=dl$meanLight,type="c-l",init=0)
cll2 <- lm.circular(y=dl$mean,x=dl$maxLight,type="c-l",init=0)
cll3 <- lm.circular(y=dl$mean,x=dl$minLight,type="c-l",init=0)



# Split by environmental context

d2 <- subset(d, fish_age..dph. >17 & fish_age..dph. <25)
d2l <- subset(dl, fish_age..dph. >17 & fish_age..dph. <25)
  
dne <- subset(d2, depth..m.=="18" & location=="near" & tide=="ebb")
dnf <- subset(d2, depth..m.=="18" & location=="near" & tide=="flood")
dfe <- subset(d2, depth..m.=="18" & location=="far" & tide=="ebb")
sne <- subset(d2, depth..m.=="9" & location=="near" & tide=="ebb")

dnel <- subset(d2l, depth..m.=="18" & location=="near" & tide=="ebb")
dnfl <- subset(d2l, depth..m.=="18" & location=="near" & tide=="flood")
dfel <- subset(d2l, depth..m.=="18" & location=="far" & tide=="ebb")
snel <- subset(d2l, depth..m.=="9" & location=="near" & tide=="ebb")

# ---------------

cltdne <- lm.circular(y=dne$mean,x=dne$meanTemp,type="c-l",init=0)
cltdne$p.values
#significant finding - mean temp correlated with mean bearing in deep,near,ebb


circStatsdne <- circStats(dne) # run the statistics
dne$theta <- dne$mean
dne$count <- 11

pT <- ggplot(dne) + polar() + # set the data and the polar layout
  # draw circles on the plot for the r-axis
  geom_line(data=circle, aes(theta,rho), color="gray") +
  geom_line(data=circle2, aes(theta,rho), color="gray") +
  # add the data
  geom_point(aes(x=theta, y=count,fill=meanTemp),pch=21, size=3) + # add the bearings as stacked dots
  scale_fill_gradient(limits=c(84.8, 87), low="blue", high="red", name="Temp (F)") +
  geom_segment(aes(x=bearing, y=0, xend=bearing, yend=r*10, color=signif), size = .6, data=circStatsdne) + # add the mean lines, w/ type determined by significance
  scale_color_manual(values=c("gray", "black"),guide=FALSE) + # set black line to be significant
  # groom the plot a bit
  scale_x_continuous("", limits=c(0,360), breaks=seq(0,360-1,by=45), labels=c("N","","","","S","","","")) + # set the X axis (angle axis)
  theme_bw() + labs(color="Statistically Significant") + # remove the gray background and adjust legend label
  theme(legend.title = element_text(size=8)) + theme(legend.text = element_text(size=6)) + # adjust legend text sizes
  theme(panel.grid.minor = element_line(size = .2, color = "white")) + theme(panel.grid.major = element_line(size=.5, color = "white")) + # remove x-axis grid, make y-axis grid black
  theme(plot.title = element_text(face="bold", size=18)) + theme(axis.title.y = element_text(size=8, face="bold")) + # format the plot title and r-axis label
  theme(axis.text.x = element_text(size=14, vjust=.05,face="bold")) + theme(axis.text.y = element_blank(), axis.ticks = element_blank()) + # adjust axis text
  theme(strip.text.x = element_text(size=6, face="bold"), strip.background = element_rect(colour="black", fill="#CCCCFF")) +# format facet labels
  theme(legend.title=element_text(size=12),legend.text=element_text(size=10))

pT


meanG <- mean(dne$meanTemp)
dne$tempG <- 0
dne[dne$meanTemp <= meanG, c("tempG")] <- "Low Temp"
dne[dne$meanTemp > meanG, c("tempG")] <- "High Temp"
dne$tempG <- as.factor(dne$tempG)

tempStats <- ddply(dne, ~tempG, circStats)


dbt <- ddply(dne,.(tempG), function(x, bin) {  # 'splitting' by method here b/c theres only 1 method and its easier than changing code to remove ddply
  x$mean <- as.numeric(round_any(x$mean, bin))
  x$mean[x$mean==360] <- 0
  
  # create a data.frame with count
  counts <- count(x, "mean")
  
  # repeat each point the appropriate number of times
  d <- adply(counts, 1, function(x) {
    data.frame(theta=x$mean, count=1:x$freq)
  }, .expand=F)
  
  # make the scale prettier
  d$count <- 10 + d$count
  
  return(d)
}, bin=bin)



p2t <- ggplot(dbt) + polar() + # set the data and the polar layout
  # draw circles on the plot for the r-axis
  geom_line(data=circle, aes(theta,rho), color="gray") +
  geom_line(data=circle2, aes(theta,rho), color="gray") +
  # add the data
  geom_point(aes(x=theta, y=count,fill=tempG),pch=21, size=2) + # add the bearings as stacked dots
  scale_fill_manual(values=c("red","blue")) +
  geom_segment(aes(x=bearing, y=0, xend=bearing, yend=r*10, color=signif), size = .6, data=tempStats) + # add the mean lines, w/ type determined by significance
  scale_color_manual(values=c("gray", "black")) + # set black line to be signif
  facet_wrap(~tempG) +
  # groom the plot a bit
  scale_x_continuous("", limits=c(0,360), breaks=seq(0,360-1,by=45), labels=c("N","","","","S","","","")) + # set the X axis (angle axis)
  theme_bw() + labs(color="Statistically Significant") + # remove the gray background and adjust legend label
  theme(legend.title = element_text(size=8)) + theme(legend.text = element_text(size=6)) + # adjust legend text sizes
  theme(panel.grid.minor = element_line(size = .2, color = "white")) + theme(panel.grid.major = element_line(size=.5, color = "white")) + # remove x-axis grid, make y-axis grid black
  theme(plot.title = element_text(face="bold", size=18)) + theme(axis.title.y = element_text(size=8, face="bold")) + # format the plot title and r-axis label
  theme(axis.text.x = element_text(size=8, vjust=.05)) + theme(axis.text.y = element_blank(), axis.ticks = element_blank()) + # adjust axis text
  theme(strip.text.x = element_text(size=6, face="bold"), strip.background = element_rect(colour="black", fill="#CCCCFF"))  # format facet labels

p2t


# todo: how does temp vary through time?
# todo: how does temp vary with current direction/distance?

# is this consistent at all ages?
dneT <- subset(d, depth..m.=="18" & location=="near" & tide=="ebb")
cltdneT <- lm.circular(y=dneT$mean,x=dneT$meanTemp,type="c-l",init=0)
cltdneT$p.values
# NOPE


circStatsdneT <- circStats(dneT) # run the statistics
dneT$theta <- dneT$mean
dneT$count <- 11

pdT <- ggplot(dneT) + polar() + # set the data and the polar layout
  # draw circles on the plot for the r-axis
  geom_line(data=circle, aes(theta,rho), color="gray") +
  geom_line(data=circle2, aes(theta,rho), color="gray") +
  # add the data
  geom_point(aes(x=theta, y=count,fill=meanTemp),pch=21, size=3) + # add the bearings as stacked dots
  scale_fill_gradient(limits=c(84.8, 87), low="blue", high="red", name="Temp (F)") +
  geom_segment(aes(x=bearing, y=0, xend=bearing, yend=r*10, color=signif), size = .6, data=circStatsdneT) + # add the mean lines, w/ type determined by significance
  scale_color_manual(values=c("gray", "black"),guide=FALSE) + # set black line to be significant
  # groom the plot a bit
  scale_x_continuous("", limits=c(0,360), breaks=seq(0,360-1,by=45), labels=c("N","","","","S","","","")) + # set the X axis (angle axis)
  theme_bw() + labs(color="Statistically Significant") + # remove the gray background and adjust legend label
  theme(legend.title = element_text(size=8)) + theme(legend.text = element_text(size=6)) + # adjust legend text sizes
  theme(panel.grid.minor = element_line(size = .2, color = "white")) + theme(panel.grid.major = element_line(size=.5, color = "white")) + # remove x-axis grid, make y-axis grid black
  theme(plot.title = element_text(face="bold", size=18)) + theme(axis.title.y = element_text(size=8, face="bold")) + # format the plot title and r-axis label
  theme(axis.text.x = element_text(size=14, vjust=.05,face="bold")) + theme(axis.text.y = element_blank(), axis.ticks = element_blank()) + # adjust axis text
  theme(strip.text.x = element_text(size=6, face="bold"), strip.background = element_rect(colour="black", fill="#CCCCFF")) +# format facet labels
  theme(legend.title=element_text(size=12),legend.text=element_text(size=10))

pdT


# try groups
tg <- dneT
tg$ageG <- 0
tg[tg$fish_age..dph. <= 10,c("ageG")] <- "10&Under"
tg[tg$fish_age..dph. > 20,c("ageG")] <- "22&Over"
tg[tg$ageG == 0, c("ageG")] <- "12to20"
tg$ageG <- as.factor(tg$ageG)



circStatsAgeG <- ddply(tg,.(ageG),circStats) # run the statistics
circStatsAgeG$age <- circStatsAgeG$ageG
for (i in 1:nrow(circStatsAgeG)){  
  if (circStatsAgeG[i,c("bearing")]<=0){
    circStatsAgeG[i,c("bearing")] <- 360 + circStatsAgeG[i,c("bearing")]
  }
}

cltdneT <- lm.circular(y=dneT$mean,x=dneT$meanTemp,type="c-l",init=0)

tg10 <- lm.circular(y=tg[tg$ageG == "10&Under",c("mean")],x=tg[tg$ageG == "10&Under",c("meanTemp")],type="c-l",init=0)
tg20 <- lm.circular(y=tg[tg$ageG == "12to20",c("mean")],x=tg[tg$ageG == "12to20",c("meanTemp")],type="c-l",init=0)
tg30 <- lm.circular(y=tg[tg$ageG == "22&Over",c("mean")],x=tg[tg$ageG == "22&Over",c("meanTemp")],type="c-l",init=0)


pdoT <- ggplot(tg) + polar() + # set the data and the polar layout
  # draw circles on the plot for the r-axis
  geom_line(data=circle, aes(theta,rho), color="gray") +
  geom_line(data=circle2, aes(theta,rho), color="gray") +
  # add the data
  geom_point(aes(x=theta, y=count,fill=meanTemp),pch=21, size=3) + # add the bearings as stacked dots
  scale_fill_gradient(limits=c(84.8, 87), low="blue", high="red", name="Temp (F)") +
  geom_segment(aes(x=bearing, y=0, xend=bearing, yend=r*10, color=signif), size = .6, data=circStatsAgeG) + # add the mean lines, w/ type determined by significance
  scale_color_manual(values=c("gray", "black"),guide=FALSE) + # set black line to be significant
  # groom the plot a bit
  facet_wrap(~ageG) +
  scale_x_continuous("", limits=c(0,360), breaks=seq(0,360-1,by=45), labels=c("N","","","","S","","","")) + # set the X axis (angle axis)
  theme_bw() + labs(color="Statistically Significant") + # remove the gray background and adjust legend label
  theme(legend.title = element_text(size=8)) + theme(legend.text = element_text(size=6)) + # adjust legend text sizes
  theme(panel.grid.minor = element_line(size = .2, color = "white")) + theme(panel.grid.major = element_line(size=.5, color = "white")) + # remove x-axis grid, make y-axis grid black
  theme(plot.title = element_text(face="bold", size=18)) + theme(axis.title.y = element_text(size=8, face="bold")) + # format the plot title and r-axis label
  theme(axis.text.x = element_text(size=14, vjust=.05,face="bold")) + theme(axis.text.y = element_blank(), axis.ticks = element_blank()) + # adjust axis text
  theme(strip.text.x = element_text(size=6, face="bold"), strip.background = element_rect(colour="black", fill="#CCCCFF")) +# format facet labels
  theme(legend.title=element_text(size=12),legend.text=element_text(size=10))

pdoT



cltdnf <- lm.circular(y=dnf$mean,x=dnf$meanTemp,type="c-l",init=0)
# promising

cltdfe <- lm.circular(y=dfe$mean,x=dfe$meanTemp,type="c-l",init=0)

cltsne <- lm.circular(y=sne$mean,x=sne$meanTemp,type="c-l",init=0)

# ---------------

clldne <- lm.circular(y=dnel$mean,x=dnel$meanLight,type="c-l",init=0)
clldne$p.values
# significant

circStatsdnel <- circStats(dnel) # run the statistics
dnel$theta <- dnel$mean
dnel$count <- 11

pT <- ggplot(dnel) + polar() + # set the data and the polar layout
  # draw circles on the plot for the r-axis
  geom_line(data=circle, aes(theta,rho), color="gray") +
  geom_line(data=circle2, aes(theta,rho), color="gray") +
  # add the data
  geom_point(aes(x=theta, y=count,fill=meanLight),pch=21, size=3) + # add the bearings as stacked dots
  scale_fill_gradient(limits=c(178.5, 1984), low="blue", high="red") +
  geom_segment(aes(x=bearing, y=0, xend=bearing, yend=r*10, color=signif), size = .6, data=circStatsdnel) + # add the mean lines, w/ type determined by significance
  scale_color_manual(values=c("gray", "black")) + # set black line to be significant
  # groom the plot a bit
  scale_x_continuous("", limits=c(0,360), breaks=seq(0,360-1,by=45), labels=c("N","","","","S","","","")) + # set the X axis (angle axis)
  theme_bw() + labs(color="Statistically Significant") + # remove the gray background and adjust legend label
  theme(legend.title = element_text(size=8)) + theme(legend.text = element_text(size=6)) + # adjust legend text sizes
  theme(panel.grid.minor = element_line(size = .2, color = "white")) + theme(panel.grid.major = element_line(size=.5, color = "white")) + # remove x-axis grid, make y-axis grid black
  theme(plot.title = element_text(face="bold", size=18)) + theme(axis.title.y = element_text(size=8, face="bold")) + # format the plot title and r-axis label
  theme(axis.text.x = element_text(size=5, vjust=.05)) + theme(axis.text.y = element_blank(), axis.ticks = element_blank()) + # adjust axis text
  theme(strip.text.x = element_text(size=6, face="bold"), strip.background = element_rect(colour="black", fill="#CCCCFF"))  # format facet labels

pT


meanL <- mean(dnel$meanLight)
dnel$tempL <- 0
dnel[dnel$meanLight <= meanL, c("tempL")] <- "lowLight"
dnel[dnel$meanLight > meanL, c("tempL")] <- "highLight"
dnel$tempG <- as.factor(dnel$tempL)

lightStats <- ddply(dnel, ~tempL, circStats)


dbl <- ddply(dnel,.(tempL), function(x, bin) {  # 'splitting' by method here b/c theres only 1 method and its easier than changing code to remove ddply
  x$mean <- as.numeric(round_any(x$mean, bin))
  x$mean[x$mean==360] <- 0
  
  # create a data.frame with count
  counts <- count(x, "mean")
  
  # repeat each point the appropriate number of times
  d <- adply(counts, 1, function(x) {
    data.frame(theta=x$mean, count=1:x$freq)
  }, .expand=F)
  
  # make the scale prettier
  d$count <- 10 + d$count
  
  return(d)
}, bin=bin)



p2l <- ggplot(dbl) + polar() + # set the data and the polar layout
  # draw circles on the plot for the r-axis
  geom_line(data=circle, aes(theta,rho), color="gray") +
  geom_line(data=circle2, aes(theta,rho), color="gray") +
  # add the data
  geom_point(aes(x=theta, y=count), size=1) + # add the bearings as stacked dots
  geom_segment(aes(x=bearing, y=0, xend=bearing, yend=r*10, color=signif), size = .6, data=lightStats) + # add the mean lines, w/ type determined by significance
  scale_color_manual(values=c("gray", "black")) + # set black line to be signif
  facet_wrap(~tempL) +
  # groom the plot a bit
  scale_x_continuous("", limits=c(0,360), breaks=seq(0,360-1,by=45), labels=c("N","","","","S","","","")) + # set the X axis (angle axis)
  theme_bw() + labs(color="Statistically Significant") + # remove the gray background and adjust legend label
  theme(legend.title = element_text(size=8)) + theme(legend.text = element_text(size=6)) + # adjust legend text sizes
  theme(panel.grid.minor = element_line(size = .2, color = "white")) + theme(panel.grid.major = element_line(size=.5, color = "white")) + # remove x-axis grid, make y-axis grid black
  theme(plot.title = element_text(face="bold", size=18)) + theme(axis.title.y = element_text(size=8, face="bold")) + # format the plot title and r-axis label
  theme(axis.text.x = element_text(size=5, vjust=.05)) + theme(axis.text.y = element_blank(), axis.ticks = element_blank()) + # adjust axis text
  theme(strip.text.x = element_text(size=6, face="bold"), strip.background = element_rect(colour="black", fill="#CCCCFF"))  # format facet labels

p2l



clldnf <- lm.circular(y=dnfl$mean,x=dnfl$meanLight,type="c-l",init=0)
# significant

clldfe <- lm.circular(y=dfel$mean,x=dfel$meanLight,type="c-l",init=0)
# promising

cllsne <- lm.circular(y=snel$mean,x=snel$meanLight,type="c-l",init=0)


# ---------------

cll2dne <- lm.circular(y=dnel$mean,x=dnel$maxLight,type="c-l",init=0)
# significant

cll2dnf <- lm.circular(y=dnfl$mean,x=dnfl$maxLight,type="c-l",init=0)
# significant

cll2dfe <- lm.circular(y=dfel$mean,x=dfel$maxLight,type="c-l",init=0)
# promising

cll2sne <- lm.circular(y=snel$mean,x=snel$maxLight,type="c-l",init=0)

# ---------------

cll3dne <- lm.circular(y=dnel$mean,x=dnel$minLight,type="c-l",init=0)

cll3dnf <- lm.circular(y=dnfl$mean,x=dnfl$minLight,type="c-l",init=0)

cll3dfe <- lm.circular(y=dfel$mean,x=dfel$minLight,type="c-l",init=0)

cll3sne <- lm.circular(y=snel$mean,x=snel$minLight,type="c-l",init=0)






## Precision in relation to light/temp -------------


rmodelt <- lm(r ~ meanTemp, data=d)
summary(rmodelt)


ggplot(d,aes(y=r,x=meanTemp)) + geom_point() + geom_smooth(method="lm") 


# what about only dne?

rmodelt2 <- lm(r ~ meanTemp, data=dne)
summary(rmodelt2)


ggplot(d,aes(y=r,x=meanTemp)) + geom_point() + geom_smooth(method="lm") 



#get rid of some weird r light values
dl <- subset(d, meanLight <= 2500)

rmodell <- lm(r ~ meanLight, data=dl)
summary(rmodell)

ggplot(dl,aes(y=r,x=meanLight)) + geom_point() + geom_smooth(method="lm") + scale_x_continuous(limits=c(0,2500))


rmodell2 <- lm(r ~ meanLight, data=dnel)
summary(rmodell2)


ggplot(dl,aes(y=r,x=meanLight)) + geom_point() + geom_smooth(method="lm") + scale_x_continuous(limits=c(0,2500))





## Swimming speed in relation to light/temp -------------



smodelt <- lm(speed.mean ~ meanTemp, data=d)
summary(smodelt)

ggplot(d,aes(y=speed.mean,x=meanTemp)) + geom_point() + geom_smooth() 



smodell <- lm(speed.mean ~ meanLight, data=dl)
summary(smodell)

ggplot(dl,aes(y=speed.mean,x=meanLight)) + geom_point() + geom_smooth() + scale_x_continuous(limits=c(0,2500))





## Temp & Light in Different Environmental Conditions -----------------


# temp 
temodel <- lm(meanTemp ~ tide+depth..m.+location, data=d2)
summary(temodel)
# significant finding - its warmer at the surface


pte <- ggplot(d2) +
  geom_boxplot(aes(tide,meanTemp,fill=location,linetype=depth..m.))
pte


# light

#get rid of some weird r light values
d2l <- subset(d2, meanLight <= 2500)

lemodel <- lm(meanLight ~ tide+depth..m.+location, data=d2l)
summary(lemodel)
# significant finding - its brighter near the surface

ple <- ggplot(d2l) +
  geom_boxplot(aes(tide,meanLight,fill=location,linetype=depth..m.))
ple




## Multiple Regression for the Precision ------------

rmrmodel <- lm(r ~ fish_age..dph.+meanTemp+meanLight,data=d)
summary(rmrmodel)





## Orientation to the transect ---------------

tran_lat <- 16.8175333
tran_lon <- -88.0766

dneT <- subset(d, depth..m.=="18" & location=="near" & tide=="ebb")

dneT$tran_lat <- tran_lat
dneT$tran_lon <- tran_lon

dneT$dy <- dneT$lat - dneT$tran_lat
dneT$dx <- dneT$lon - dneT$tran_lon
dneT$tran_ang <- atan(dneT$dy/dneT$dx) * 180 / pi

dneT[dneT$dy < 0, c("tran_ang")] <- -dneT[dneT$dy < 0, c("tran_ang")] + 270
dneT[dneT$dy > 0, c("tran_ang")] <- -dneT[dneT$dy > 0, c("tran_ang")] + 270



# correct angles so headings are relative to transect direction
dneT$tran_bearing <- dneT$tran_ang - dneT$mean
dneT = dneT[!is.na(dneT$tran_bearing),] # remove trials with no current data

# correct for negative angle values produced
for (i in 1:nrow(dneT)){  
  if (dneT[i,c("tran_bearing")]<=0){
    dneT[i,c("tran_bearing")] <- 360 + dneT[i,c("tran_bearing")]
  }
}


# repeat for all larvae - SUPMAT PLOT
d$tran_lat <- tran_lat
d$tran_lon <- tran_lon

d$dy <- d$lat - d$tran_lat
d$dx <- d$lon - d$tran_lon
d$tran_ang <- atan(d$dy/d$dx) * 180 / pi

d[d$dy < 0, c("tran_ang")] <- -d[d$dy < 0, c("tran_ang")] + 270
d[d$dy > 0, c("tran_ang")] <- -d[d$dy > 0, c("tran_ang")] + 270



# correct angles so headings are relative to transect direction
d$tran_bearing <- d$tran_ang - d$mean
d = d[!is.na(d$tran_bearing),] # remove trials with no current data

# correct for negative angle values produced
for (i in 1:nrow(d)){  
  if (d[i,c("tran_bearing")]<=0){
    d[i,c("tran_bearing")] <- 360 + d[i,c("tran_bearing")]
  }
}


# define function to look at current bearings
circStatsTran <- function(x) {
  n = nrow(x)
  tranBearing = mean.circular(x$tran_bearing)
  sdbearing = sd.circular(x$tran_bearing)
  meanR = mean(x$r, na.rm=T)
  sdR = sd(x$r, na.rm=T)
  ray = rayleigh.test(na.omit(x$tran_bearing))
  r = ray$statistic
  p = ray$p.value
  signif = p < 0.05
  return(data.frame(n, tranBearing,sdbearing, meanR, sdR, r, p, signif))
}

d2 <- subset(d, fish_age..dph. >17 & fish_age..dph. <25)

dbTran <- ddply(d2,.(tide,depth..m.,location), function(x, bin) {  # 'splitting' by method here b/c theres only 1 method and its easier than changing code to remove ddply
  x$mean <- as.numeric(round_any(x$tran_bearing, bin))
  x$mean[x$mean==360] <- 0
  
  # create a data.frame with count
  counts <- count(x, "mean")
  
  # repeat each point the appropriate number of times
  d <- adply(counts, 1, function(x) {
    data.frame(theta=x$mean, count=1:x$freq)
  }, .expand=F)
  
  # make the scale prettier
  d$count <- 10 + d$count
  
  return(d)
}, bin=bin)



circStatsTra <- ddply(d2,.(tide,depth..m.,location),circStatsTran)
circStatsTra[circStatsTra$tranBearing <= 0, c("tranBearing")] <- circStatsTra[circStatsTra$tranBearing <= 0, c("tranBearing")] + 360


pTran <- ggplot(dbTran) + polar() + # set the data and the polar layout
  # draw circles on the plot for the r-axis
  geom_line(data=circle, aes(theta,rho), color="gray") +
  geom_line(data=circle2, aes(theta,rho), color="gray") +
  # add the data
  geom_point(aes(x=theta, y=count), size=1) + # add the bearings as stacked dots
  geom_segment(aes(x=tranBearing, y=0, xend=tranBearing, yend=r*10, color=signif), size = .6, data=circStatsTra) + # add the mean lines, w/ type determined by significance
  scale_color_manual(values=c("gray", "black")) + # set black line to be significant
  # groom the plot a bit
  facet_wrap(~tide+location+depth..m.) +
  theme_bw() + ylab("Rho-Value") +
  theme(panel.grid.minor = element_line(size = .2, color = "white")) + 
  theme(panel.grid.major = element_line(size=.5, color = "white")) +
  theme(axis.title.y = element_text(size=16, face="bold")) + 
  theme(axis.text.y = element_text(size=16)) +
  scale_y_continuous(breaks=c(0,3,6,10),labels=c(0,0.3,0.6,1)) +
  theme(axis.text.x = element_text(face="bold",size=12, vjust=.05)) + 
  scale_x_continuous("", limits=c(0,360), breaks=seq(0,360-1,by=45), labels=c("Towards Transect","","","","Away from Transect","","","")) + 
  theme(strip.background = element_blank(),strip.text.x = element_blank())

pTran


# repeat analysis with ontogeny groups
tg <- dneT
tg$ageG <- 0
tg[tg$fish_age..dph. <= 10,c("ageG")] <- "10&Under"
tg[tg$fish_age..dph. > 20,c("ageG")] <- "22&Over"
tg[tg$ageG == 0, c("ageG")] <- "12to20"
tg$ageG <- as.factor(tg$ageG)


circStatsAgeT <- ddply(tg,.(ageG),circStatsTran) # run the statistics
circStatsAgeT$age <- circStatsAgeT$ageG
for (i in 1:nrow(circStatsAgeT)){  
  if (circStatsAgeT[i,c("tranBearing")]<=0){
    circStatsAgeT[i,c("tranBearing")] <- 360 + circStatsAgeT[i,c("tranBearing")]
  }
}

# bin angles for plot
tg$theta <- as.numeric(round_any(tg$tran_bearing, bin))
tg$theta[tg$theta==360] <- 0
# create a data.frame with count
counts <- count(tg, vars=c("theta","ageG"))

# repeat each point the appropriate number of times
tbg <- adply(counts, 1, function(x) {
  data.frame(theta=x$theta,age=x$ageG, count=1:x$freq)
}, .expand=F)

# make the scale prettier
tbg$count <- 10 + tbg$count



pgaAz <- ggplot(tbg) + polar() + # set the data and the polar layout
  # draw circles on the plot for the r-axis
  geom_line(data=circle, aes(theta,rho), color="gray") +
  geom_line(data=circle2, aes(theta,rho), color="gray") +
  # add the data
  geom_point(aes(x=theta, y=count), size=0.8) + # add the bearings as stacked dots
  geom_segment(aes(x=tranBearing, y=0, xend=tranBearing, yend=r*10, color=signif), size = .6, data=circStatsAgeT) + # add the mean lines, w/ type determined by significance
  scale_color_manual(values=c("gray", "black")) + # set black line to be significant
  facet_wrap(~age) +
  # groom the plot a bit
  scale_x_continuous("", limits=c(0,360), breaks=seq(0,360-1,by=45), labels=c("Towards Transect","","","","Away from Transect","","","")) + # set the X axis (angle axis)
  theme_bw() + labs(color="Statistically Significant") + # remove the gray background and adjust legend label
  theme(legend.title = element_text(size=8)) + theme(legend.text = element_text(size=6)) + # adjust legend text sizes
  theme(panel.grid.minor = element_line(size = .2, color = "white")) + theme(panel.grid.major = element_line(size=.5, color = "white")) + # remove x-axis grid, make y-axis grid black
  theme(plot.title = element_text(face="bold", size=18)) + theme(axis.title.y = element_text(size=8, face="bold")) + # format the plot title and r-axis label
  theme(axis.text.x = element_text(size=5, vjust=.05)) + theme(axis.text.y = element_blank(), axis.ticks = element_blank()) + # adjust axis text
  theme(strip.text.x = element_text(size=6, face="bold"), strip.background = element_rect(colour="black", fill="#CCCCFF"))  # format facet labels

pgaAz





## Orientation to the wetlab --------------

lab_lat <- 16.8154833
lab_lon <- -88.0816

dneT <- subset(d, depth..m.=="18" & location=="near" & tide=="ebb")

dneT$lab_lat <- lab_lat
dneT$lab_lon <- lab_lon

dneT$dy2 <- dneT$lat - dneT$lab_lat 
dneT$dx2 <- dneT$lon - dneT$lab_lon 
dneT$lab_ang <- atan(dneT$dy2/dneT$dx2) * 180 / pi

dneT[dneT$dy < 0, c("lab_ang")] <- -dneT[dneT$dy < 0, c("lab_ang")] + 270
dneT[dneT$dy > 0, c("lab_ang")] <- -dneT[dneT$dy > 0, c("lab_ang")] + 270


# correct angles so headings are relative to lab direction
dneT$lab_bearing <- dneT$lab_ang - dneT$mean
dneT = dneT[!is.na(dneT$lab_bearing),] # remove trials with no current data

# correct for negative angle values produced
for (i in 1:nrow(dneT)){  
  if (dneT[i,c("lab_bearing")]<=0){
    dneT[i,c("lab_bearing")] <- 360 + dneT[i,c("lab_bearing")]
  }
}

# repeat for all trials - SUPMAT PLOT
d$lab_lat <- lab_lat
d$lab_lon <- lab_lon

d$dy2 <- d$lat - d$lab_lat 
d$dx2 <- d$lon - d$lab_lon 
d$lab_ang <- atan(d$dy2/d$dx2) * 180 / pi

d[d$dy < 0, c("lab_ang")] <- -d[d$dy < 0, c("lab_ang")] + 270
d[d$dy > 0, c("lab_ang")] <- -d[d$dy > 0, c("lab_ang")] + 270


# correct angles so headings are relative to lab direction
d$lab_bearing <- d$lab_ang - d$mean
d = d[!is.na(d$lab_bearing),] # remove trials with no current data

# correct for negative angle values produced
for (i in 1:nrow(d)){  
  if (d[i,c("lab_bearing")]<=0){
    d[i,c("lab_bearing")] <- 360 + d[i,c("lab_bearing")]
  }
}


# define function to look at current bearings
circStatsLab <- function(x) {
  n = nrow(x)
  labBearing = mean.circular(x$lab_bearing)
  sdbearing = sd.circular(x$lab_bearing)
  meanR = mean(x$r, na.rm=T)
  sdR = sd(x$r, na.rm=T)
  ray = rayleigh.test(na.omit(x$lab_bearing))
  r = ray$statistic
  p = ray$p.value
  signif = p < 0.05
  return(data.frame(n, labBearing,sdbearing, meanR, sdR, r, p, signif))
}

d2 <- subset(d, fish_age..dph. >17 & fish_age..dph. <25)

dblab <- ddply(d2,.(tide,depth..m.,location), function(x, bin) {  # 'splitting' by method here b/c theres only 1 method and its easier than changing code to remove ddply
  x$mean <- as.numeric(round_any(x$lab_bearing, bin))
  x$mean[x$mean==360] <- 0
  
  # create a data.frame with count
  counts <- count(x, "mean")
  
  # repeat each point the appropriate number of times
  d <- adply(counts, 1, function(x) {
    data.frame(theta=x$mean, count=1:x$freq)
  }, .expand=F)
  
  # make the scale prettier
  d$count <- 10 + d$count
  
  return(d)
}, bin=bin)



circStatsLa <- ddply(d2,.(tide,depth..m.,location),circStatsLab)
circStatsLa[circStatsLa$labBearing <= 0, c("labBearing")] <- circStatsLa[circStatsLa$labBearing <= 0, c("labBearing")] + 360


plab <- ggplot(dblab) + polar() + # set the data and the polar layout
  # draw circles on the plot for the r-axis
  geom_line(data=circle, aes(theta,rho), color="gray") +
  geom_line(data=circle2, aes(theta,rho), color="gray") +
  # add the data
  geom_point(aes(x=theta, y=count), size=1) + # add the bearings as stacked dots
  geom_segment(aes(x=labBearing, y=0, xend=labBearing, yend=r*10, color=signif), size = .6, data=circStatsLa) + # add the mean lines, w/ type determined by significance
  scale_color_manual(values=c("gray", "black")) + # set black line to be significant
  # groom the plot a bit
  facet_wrap(~tide+location+depth..m.) +
  theme_bw() + ylab("Rho-Value") +
  theme(panel.grid.minor = element_line(size = .2, color = "white")) + 
  theme(panel.grid.major = element_line(size=.5, color = "white")) +
  theme(axis.title.y = element_text(size=16, face="bold")) + 
  theme(axis.text.y = element_text(size=16)) +
  scale_y_continuous(breaks=c(0,3,6,10),labels=c(0,0.3,0.6,1)) +
  theme(axis.text.x = element_text(face="bold",size=12, vjust=.05)) + 
  scale_x_continuous("", limits=c(0,360), breaks=seq(0,360-1,by=45), labels=c("Towards Lab","","","","Away from Lab","","","")) + 
  theme(strip.background = element_blank(),strip.text.x = element_blank())

plab


# repeat analysis with ontogeny groups
tg <- dneT
tg$ageG <- 0
tg[tg$fish_age..dph. <= 10,c("ageG")] <- "10&Under"
tg[tg$fish_age..dph. > 20,c("ageG")] <- "22&Over"
tg[tg$ageG == 0, c("ageG")] <- "12to20"
tg$ageG <- as.factor(tg$ageG)


circStatsAgeL <- ddply(tg,.(ageG),circStatsLab) # run the statistics
circStatsAgeL$age <- circStatsAgeL$ageG
for (i in 1:nrow(circStatsAgeL)){  
  if (circStatsAgeL[i,c("labBearing")]<=0){
    circStatsAgeL[i,c("labBearing")] <- 360 + circStatsAgeL[i,c("labBearing")]
  }
}

# bin angles for plot
tg$theta <- as.numeric(round_any(tg$lab_bearing, bin))
tg$theta[tg$theta==360] <- 0
# create a data.frame with count
counts <- count(tg, vars=c("theta","ageG"))

# repeat each point the appropriate number of times
tbg <- adply(counts, 1, function(x) {
  data.frame(theta=x$theta,age=x$ageG, count=1:x$freq)
}, .expand=F)

# make the scale prettier
tbg$count <- 10 + tbg$count



pgaAz <- ggplot(tbg) + polar() + # set the data and the polar layout
  # draw circles on the plot for the r-axis
  geom_line(data=circle, aes(theta,rho), color="gray") +
  geom_line(data=circle2, aes(theta,rho), color="gray") +
  # add the data
  geom_point(aes(x=theta, y=count), size=0.8) + # add the bearings as stacked dots
  geom_segment(aes(x=labBearing, y=0, xend=labBearing, yend=r*10, color=signif), size = .6, data=circStatsAgeL) + # add the mean lines, w/ type determined by significance
  scale_color_manual(values=c("gray", "black")) + # set black line to be significant
  facet_wrap(~age) +
  # groom the plot a bit
  scale_x_continuous("", limits=c(0,360), breaks=seq(0,360-1,by=45), labels=c("Towards Lab","","","","Away from Lab","","","")) + # set the X axis (angle axis)
  theme_bw() + labs(color="Statistically Significant") + # remove the gray background and adjust legend label
  theme(legend.title = element_text(size=8)) + theme(legend.text = element_text(size=6)) + # adjust legend text sizes
  theme(panel.grid.minor = element_line(size = .2, color = "white")) + theme(panel.grid.major = element_line(size=.5, color = "white")) + # remove x-axis grid, make y-axis grid black
  theme(plot.title = element_text(face="bold", size=18)) + theme(axis.title.y = element_text(size=8, face="bold")) + # format the plot title and r-axis label
  theme(axis.text.x = element_text(size=5, vjust=.05)) + theme(axis.text.y = element_blank(), axis.ticks = element_blank()) + # adjust axis text
  theme(strip.text.x = element_text(size=6, face="bold"), strip.background = element_rect(colour="black", fill="#CCCCFF"))  # format facet labels

pgaAz





## Orientation to the moon azimuth -----------------


# get rid of deployments where the sun is moon is directly overhead
dm <- subset(d, abs(moonAlt) < 85)


# correct angles so headings are relative to current direction
dm$moon_bearing <- dm$mean - dm$moonAz

# correct for negative angle values produced
for (i in 1:nrow(dm)){  
  if (dm[i,c("moon_bearing")]<=0){
    dm[i,c("moon_bearing")] <- 360 + dm[i,c("moon_bearing")]
  }
}

# define function to look at moon bearings
circStatsM <- function(x) {
  n = nrow(x)
  moonBearing = mean.circular(x$moon_bearing)
  sdbearing = sd.circular(x$moon_bearing)
  meanR = mean(x$r, na.rm=T)
  sdR = sd(x$r, na.rm=T)
  ray = rayleigh.test(na.omit(x$moon_bearing))
  r = ray$statistic
  p = ray$p.value
  signif = p < 0.05
  return(data.frame(n, moonBearing,sdbearing, meanR, sdR, r, p, signif))
}

dm$expcode <- 1
dbM <- ddply(dm,.(expcode), function(x, bin) {  # 'splitting' by method here b/c theres only 1 method and its easier than changing code to remove ddply
  x$mean <- as.numeric(round_any(x$moon_bearing, bin))
  x$mean[x$mean==360] <- 0
  
  # create a data.frame with count
  counts <- count(x, "mean")
  
  # repeat each point the appropriate number of times
  d <- adply(counts, 1, function(x) {
    data.frame(theta=x$mean, count=1:x$freq)
  }, .expand=F)
  
  # make the scale prettier
  d$count <- 10 + d$count
  
  return(d)
}, bin=bin)

circStatsMM <- circStatsM(dm) # run the statistics
circStatsMM$moonBearing <- 360 + circStatsMM$moonBearing

pM <- ggplot(dbM) + polar() + # set the data and the polar layout
  # draw circles on the plot for the r-axis
  geom_line(data=circle, aes(theta,rho), color="gray") +
  geom_line(data=circle2, aes(theta,rho), color="gray") +
  # add the data
  geom_point(aes(x=theta, y=count), size=1) + # add the bearings as stacked dots
  geom_segment(aes(x=moonBearing, y=0, xend=moonBearing, yend=r*10, color=signif), size = .6, data=circStatsMM) + # add the mean lines, w/ type determined by significance
  scale_color_manual(values=c("gray", "black")) + # set black line to be signif
  # groom the plot a bit
  scale_x_continuous("", limits=c(0,360), breaks=seq(0,360-1,by=45), labels=c("Towards Moon","","","","Away from Moon","","","")) + # set the X axis (angle axis)
  theme_bw() + labs(color="Statistically Significant") + # remove the gray background and adjust legend label
  theme(legend.title = element_text(size=8)) + theme(legend.text = element_text(size=6)) + # adjust legend text sizes
  theme(panel.grid.minor = element_line(size = .2, color = "white")) + theme(panel.grid.major = element_line(size=.5, color = "white")) + # remove x-axis grid, make y-axis grid black
  theme(plot.title = element_text(face="bold", size=18)) + theme(axis.title.y = element_text(size=8, face="bold")) + # format the plot title and r-axis label
  theme(axis.text.x = element_text(size=10, vjust=.05,face="bold")) + theme(axis.text.y = element_blank(), axis.ticks = element_blank()) + # adjust axis text
  # theme(strip.text.x = element_text(size=6, face="bold"), strip.background = element_rect(colour="black", fill="#CCCCFF"))  # format facet labels
  theme(strip.background = element_blank(), strip.text.x = element_blank()) # or remove them

pM






## Moon orientation by environmental context -----------------


d2m <- subset(dm, fish_age..dph. >17 & fish_age..dph. <25)


db2M <- ddply(d2m,.(depth..m.,tide,location), function(x, bin) {  # 'splitting' by method here b/c theres only 1 method and its easier than changing code to remove ddply
  x$mean <- as.numeric(round_any(x$moon_bearing, bin))
  x$mean[x$mean==360] <- 0
  
  # create a data.frame with count
  counts <- count(x, "mean")
  
  # repeat each point the appropriate number of times
  d <- adply(counts, 1, function(x) {
    data.frame(theta=x$mean, count=1:x$freq)
  }, .expand=F)
  
  # make the scale prettier
  d$count <- 10 + d$count
  
  return(d)
}, bin=bin)

circStats2m <- ddply(d2m,.(depth..m.,tide,location),circStatsM) 
# significant finding! swim against current in the flood tide
for (i in 1:nrow(circStats2m)){  
  if (circStats2m[i,c("moonBearing")]<=0){
    circStats2m[i,c("moonBearing")] <- 360 + circStats2m[i,c("moonBearing")]
  }
}

p2m <- ggplot(db2M) + polar() + # set the data and the polar layout
  # draw circles on the plot for the r-axis
  geom_line(data=circle, aes(theta,rho), color="gray") +
  geom_line(data=circle2, aes(theta,rho), color="gray") +
  # add the data
  geom_point(aes(x=theta, y=count), size=1) + # add the bearings as stacked dots
  geom_segment(aes(x=moonBearing, y=0, xend=moonBearing, yend=r*10, color=signif), size = .6, data=circStats2m) + # add the mean lines, w/ type determined by significance
  scale_color_manual(values=c("gray", "black")) + # set black line to be significant
  facet_wrap(~tide+location+depth..m.) +
  # groom the plot a bit
  scale_x_continuous("", limits=c(0,360), breaks=seq(0,360-1,by=45), labels=c("Towards Moon","","","","Away from Moon","","","")) + # set the X axis (angle axis)
  theme_bw() + labs(color="Statistically Significant") + # remove the gray background and adjust legend label
  theme(legend.title = element_text(size=8)) + theme(legend.text = element_text(size=6)) + # adjust legend text sizes
  theme(panel.grid.minor = element_line(size = .2, color = "white")) + theme(panel.grid.major = element_line(size=.5, color = "white")) + # remove x-axis grid, make y-axis grid black
  theme(plot.title = element_text(face="bold", size=18)) + theme(axis.title.y = element_text(size=8, face="bold")) + # format the plot title and r-axis label
  theme(axis.text.x = element_text(size=5, vjust=.05)) + theme(axis.text.y = element_blank(), axis.ticks = element_blank()) + # adjust axis text
  theme(strip.text.x = element_text(size=6, face="bold"), strip.background = element_rect(colour="black", fill="#CCCCFF"))  # format facet labels

p2m



## Looking only at the best orienters-------------

# whats the avg value?
drmean <- mean(d$r)

# lets try just above that
d <- subset(d, r >= drmean) # leaves 71 observations

# NOW REPEAT EVERYTHING AND NOTHING MATTERS


dneT <- subset(d, depth..m.=="18" & location=="near" & tide=="ebb")


# bin angles for plot
dneT$theta <- as.numeric(round_any(dneT$mean, bin))
dneT$theta[dneT$theta==360] <- 0
# create a data.frame with count
counts <- count(dneT, vars="theta")

# repeat each point the appropriate number of times
db <- adply(counts, 1, function(x) {
  data.frame(theta=x$theta, count=1:x$freq)
}, .expand=F)

# make the scale prettier
db$count <- 10 + db$count


circStatsAll <- circStats(dneT) # run the statistics


pdneT <- ggplot(db) + polar() + # set the data and the polar layout
  # draw circles on the plot for the r-axis
  geom_line(data=circle, aes(theta,rho), color="gray") +
  geom_line(data=circle2, aes(theta,rho), color="gray") +
  # add the data
  geom_point(aes(x=theta, y=count), size=1) + # add the bearings as stacked dots
  geom_segment(aes(x=bearing, y=0, xend=bearing, yend=r*10, color=signif), size = .6, data=circStatsAll) + # add the mean lines, w/ type determined by significance
  scale_color_manual(values=c("gray", "black")) + # set black line to be significant
  # groom the plot a bit
  scale_x_continuous("", limits=c(0,360), breaks=seq(0,360-1,by=45), labels=c("N","","","","S","","","")) + # set the X axis (angle axis)
  theme_bw() + labs(color="Statistically Significant") + # remove the gray background and adjust legend label
  theme(legend.title = element_text(size=8)) + theme(legend.text = element_text(size=6)) + # adjust legend text sizes
  theme(panel.grid.minor = element_line(size = .2, color = "white")) + theme(panel.grid.major = element_line(size=.5, color = "white")) + # remove x-axis grid, make y-axis grid black
  theme(plot.title = element_text(face="bold", size=18)) + theme(axis.title.y = element_text(size=8, face="bold")) + # format the plot title and r-axis label
  theme(axis.text.x = element_text(size=5, vjust=.05)) + theme(axis.text.y = element_blank(), axis.ticks = element_blank()) + # adjust axis text
  theme(strip.text.x = element_text(size=6, face="bold"), strip.background = element_rect(colour="black", fill="#CCCCFF"))  # format facet labels

pdneT





## producing data table to upload for NSF ---------------

dnew <- d

valid_column_names <- make.names(names=names(dnew), unique=TRUE, allow_ = TRUE)
names(dnew) <- valid_column_names

final <- select(dnew,deploy_id,leg,date_start,date_stop,time_start,
                time_stop,disc_id,depth..m.,bin_id,
                batch_id,fish_age..dph.,lat_start,lon_start,
                lat_stop,lon_stop,sky,sea,wind_dir,wind_speed,tide,
                location,start_wp,end_wp,drift_direction,
                drift_distance,In_Behavior,Out_Behavior,n,mean,
                se.mean,kappa,se.kappa,variance,r,p.value,signif,
                turn.n,turn.abs.mean,turn.freq.gt45,speed.n,
                speed.mean,speed.sd,speed.median,meanTemp,maxTemp,
                minTemp,meanLight,maxLight,minLight,rotationRange,
                art,current_strength,current_bearing,wind_bearing,
                zenith,azimuth,azmean,lat,lon,tran_lat,tran_lon,
                tran_ang,tran_bearing,lab_lat,lab_lon,lab_ang,
                lab_bearing,moonAlt,moonAz,moonHor,moon_bearing)

write.csv(final,"NSF_Belize_Data.csv")










## plots for john's publication ------------------------


# try larger groups
tg <- dneT
tg$ageG <- 0
tg[tg$fish_age..dph. <= 10,c("ageG")] <- "10&Under"
tg[tg$fish_age..dph. > 20,c("ageG")] <- "22&Over"
tg[tg$ageG == 0, c("ageG")] <- "12to20"
tg$ageG <- as.factor(tg$ageG)


circStatsAgeG <- ddply(tg,.(ageG),circStats) # run the statistics
circStatsAgeG$age <- circStatsAgeG$ageG
for (i in 1:nrow(circStatsAgeG)){  
  if (circStatsAgeG[i,c("bearing")]<=0){
    circStatsAgeG[i,c("bearing")] <- 360 + circStatsAgeG[i,c("bearing")]
  }
}

# bin angles for plot
tg$theta <- as.numeric(round_any(tg$mean, bin))
tg$theta[tg$theta==360] <- 0
# create a data.frame with count
counts <- count(tg, vars=c("theta","ageG"))

# repeat each point the appropriate number of times
tbg <- adply(counts, 1, function(x) {
  data.frame(theta=x$theta,age=x$ageG, count=1:x$freq)
}, .expand=F)

# make the scale prettier
tbg$count <- 10 + tbg$count



pga <- ggplot(tbg) + polar() + # set the data and the polar layout
  # draw circles on the plot for the r-axis
  geom_line(data=circle, aes(theta,rho), color="gray") +
  geom_line(data=circle2, aes(theta,rho), color="gray") +
  # add the data
  geom_point(aes(x=theta, y=count), size=2) + # add the bearings as stacked dots
  geom_segment(aes(x=bearing, y=0, xend=bearing, yend=r*10), color="black", size = 1, data=circStatsAgeG) + # add the mean lines, w/ type determined by significance
  facet_grid(age~.) +
  # groom the plot a bit
  scale_x_continuous("", limits=c(0,360), breaks=seq(0,360-1,by=45), labels=c("N","","","","S","","","")) + # set the X axis (angle axis)
  theme_bw() +
  theme(panel.grid.minor = element_line(size = .2, color = "white")) + theme(panel.grid.major = element_line(size=.5, color = "white")) + # remove x-axis grid, make y-axis grid black
  theme(plot.title = element_text(face="bold", size=18)) + theme(axis.title.y = element_text(size=16, face="bold")) + # format the plot title and r-axis label
  theme(axis.text.x = element_text(size=14, face="bold", vjust=.05)) + theme(axis.text.y = element_text(size=16)) + # adjust axis text
  scale_y_continuous(breaks=c(0,3,6,10),labels=c(0,0.3,0.6,1)) +
  # theme(strip.text.x = element_text(size=6, face="bold"), strip.background = element_rect(colour="black", fill="#CCCCFF"))  # format facet labels
  theme(strip.background = element_blank(), strip.text.x = element_blank()) # or remove them

pga




