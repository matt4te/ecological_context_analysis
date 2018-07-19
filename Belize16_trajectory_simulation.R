

## model ideas:
# 
# 
# - surface model (random orientation, surface current parameters)




## Required libraries and directory ------------------------
library("ggplot2")
library("truncnorm")
library("circular")
library("plyr")

setwd("C:/Users/mattf/Desktop/RSMAS/ecological_context")


## Important parameters for the simulation -----------------------

# mean Ucrit 4 cm/s (from john majoris), so 2 cm/s routine swimming

# mean bearing to current 187.02 (n=33)
# r value 0.358
# mle.vonmises(dbz$current_bearing)
# mu = 187.2, kappa = 0.7665
# mean r = 0.47, sd r = 0.24

# currents in DNE:
# mle.vonmises(d2[d2$depth..m.==18 & d2$location=="near" & d2$tide=="ebb",c("drift_direction")])
# mu = 227, kappa = 1.487

# currents in DNF:
# mle.vonmises(d2[d2$depth..m.==18 & d2$location=="near" & d2$tide=="flood",c("drift_direction")])
# mu = 332.75, kappa = 1.456

# currents in SNE:
# mle.vonmises

# calculating the current speed:
# dbz$drif_distance is in meters over the course of the deployment
# deployments are 10 minutes, plus 5 minutes acclimation, plus 5 minutes pickup = 25 minutes
# divide by 25 minutes to get in units of m/min
# m/min * (1min/60s) * (100cm/1m) -> cm/s
# speeds <- dbz$drift_distance / 25 * (1/60) * 100
# mean(speeds)
# sd(speeds)
# mean = 3.9, sd = 2.73
# or for DNE only: mean = 4.09, sd = 3.09
# or for DNF only: mean = 4.11, sd = 2.75
# or for SNE only: mean = 5.74, sd = 2.95



## Drawing current / swimming directions from data distribution ----------

# assuming integration time of 1 hour, split into 6-hour chunks (tidal phase)

# simulate swimming angles
swimAngle <- rvonmises(n=1, mu=circular(187.2,type=c("angles"),units=c("degrees")), kappa=0.7665)

# simulate swimming speeds
swimSpeed <- rtruncnorm(n=1,a=0,b=Inf,mean=2,sd=0.5)
# or if using Ucrit
swimFast <- rtruncnorm(n=1,a=0,b=Inf,mean=4,sd=0.5)

# simulate currents in DNE
dneCurrent <- rvonmises(n=1, mu=circular(227,type=c("angles"),units=c("degrees")), kappa=1.487)
dneSpeed <- rtruncnorm(n=1,a=0,b=Inf,mean=4.09,sd=3.09)

# simulate currents in DNF
dnfCurrent <- rvonmises(n=1, mu=circular(332.75,type=c("angles"),units=c("degrees")), kappa=1.456)
dnfSpeed <- rtruncnorm(n=1,a=0,b=Inf,mean=4.11,sd=2.75)

# simulate current speeds, (lower limit of 0)
currentSpeed <- rtruncnorm(n=1,a=0,b=Inf,mean=3.9,sd=2.73)
shallowSpeed <- rtruncnorm(n=1,a=0,b=Inf,mean=5.74,sd=2.95)

# simulate R values - not using this yet
Rscaler <- rtruncnorm(n=1,a=0=b=Inf,mean=0.47,sd=0.24)




## Procedure for combining the current and swimming speeds --------------


# x and y components of current vector
Cy <- as.numeric(currentSpeed * sin(currentAngle*pi/180))
Cx <- as.numeric(currentSpeed * cos(currentAngle*pi/180))

# x and y components of swimming vector
Sy <- as.numeric(swimSpeed * sin(swimAngle*pi/180))
Sx <- as.numeric(swimSpeed * cos(swimAngle*pi/180))

# x and y component of the net vector 
Ty <- Sy + Cy # y component of displacement vector
Tx <- Sx + Cx # x component of displacement vector

# calculate displacement vector 
Tr <- sqrt(Ty^2 + Tx^2) # length of displacement vector
Tt <- atan2(Ty,Tx)*180/pi # angle of displacement vector (in degrees)

# correct for negative angle values produced
if (Tt <0) { Tt <- 360 + Tt}






## Set the parameters for the model run ---------------------------

# set integration step
dt <- 6 # in hours

# set total duration 
dur <- 10 # in days

# initital condition
xi <- 0
yi <- 0

# convert times to useful units
dt2 <- dt * (60/1) * (60/1) # hours to seconds, for dt
dur2 <- dur * (24/1) / dt # days to numbers of integrations, for dur

# create a container for the output
displacements <- data.frame(matrix(nrow=dur2,ncol=11))

# use means or draw from distribution?
simple= "FALSE"
simpleSwim = "FALSE"
noSwim = "FALSE"
fullModel = "FALSE"
randomSwimming = "FALSE"
fastSwim = "FALSE"
surfaceCurrents = "TRUE"





## Run the model ------------------------------------------

xk = xi
yk = yi

if (simple=="TRUE") { # if you're using the simplest model type
  
  # set all the static model parameters
  swimSpeed = 2 # cm/s
  swimAngleRelative = circular(187.2,type=c("angles"),units=c("degrees"))
  currentSpeed = 3.9 # cm/s
  dneCurrent = circular(227,type=c("angles"),units=c("degrees"))
  dnfCurrent = circular(332.75,type=c("angles"),units=c("degrees"))
  
  for (i in 1:dur2){ # then begin the simulation
    
    if (i%%2 == 1) { # for odd integration steps
      currentAngle <- dneCurrent # use the ebb current direction
      tide <- "ebb"
      
    } else if (i%%2 == 0) { # for even integration steps 
      currentAngle <- dnfCurrent # use the flood current direction
      tide <- "flood"
      
    }
    
    # x and y components of current vector
    Cy <- as.numeric(currentSpeed * sin(currentAngle*pi/180))
    Cx <- as.numeric(currentSpeed * cos(currentAngle*pi/180))
    
    # calculating absolute swimming angle based on relative angle
    swimAngle <- currentAngle - swimAngleRelative
    if (swimAngle < 0) { swimANgle <- 360 + swimAngle}
    
    # x and y components of swimming vector 
    # Rscaler <- rtruncnorm(n=1,a=0,b=Inf,mean=0.47,sd=0.24)
    Sy <- as.numeric(swimSpeed * sin(swimAngle*pi/180)) 
    Sx <- as.numeric(swimSpeed * cos(swimAngle*pi/180)) 
    
    # x and y component of the net vector 
    Ty <- Sy + Cy 
    Tx <- Sx + Cx 
    # polar coordinates
    Tr <- sqrt(Ty^2 + Tx^2) # length of displacement vector
    Tt <- atan2(Ty,Tx)*180/pi # angle of displacement vector (in degrees)
    if (Tt <0) { Tt <- 360 + Tt} # correct for negative angle values produced
    
    # move the larva - IMPORTANT - SWITCH Tx / Ty here because of R angle convention
    xkn <- xk + Ty * dt2 # multiply by number of seconds in timestep
    ykn <- yk + Tx * dt2 # multiply by number of seconds in timestep
    
    # collect the output of the integration step
    output <- c(xk,yk,tide,currentSpeed,currentAngle,swimSpeed,swimAngle,Tr,Tt,xkn,ykn)
    # and write it to the output file
    displacements[i,] <- output
    
    # update position for next run
    xk <- xkn
    yk <- ykn
    
  } 
} else if (simpleSwim == "TRUE"){
  
  # set static swimming parameters
  swimSpeed <- 2 # cm/s
  swimAngleRelative = circular(187.2,type=c("angles"),units=c("degrees"))
  
  for (i in 1:dur2){ # then begin the simulation
    
    if (i%%2 == 1) { # for odd integration steps
      
      # get a current angle from the distribution in the ebb tide data
      currentAngle <- rvonmises(n=1, mu=circular(227,type=c("angles"),units=c("degrees")), kappa=1.487)
      tide <- "ebb"
      
    } else if (i%%2 == 0) { # for even integration steps 
      
      # get a current angle from the distribution in the flood tide data
      currentAngle <- rvonmises(n=1, mu=circular(332.75,type=c("angles"),units=c("degrees")), kappa=1.456)
      tide <- "flood"
      
    }
    
    # generate a current speed from the distribution in the data
    currentSpeed <- rtruncnorm(n=1,a=0,b=Inf,mean=3.9,sd=2.73)
    
    # x and y components of current vector
    Cy <- as.numeric(currentSpeed * sin(currentAngle*pi/180))
    Cx <- as.numeric(currentSpeed * cos(currentAngle*pi/180))
    
    # calculating absolute swimming angle based on relative angle
    swimAngle <- currentAngle - swimAngleRelative
    if (swimAngle < 0) { swimANgle <- 360 + swimAngle}
    
    # x and y components of swimming vector
    Sy <- as.numeric(swimSpeed * sin(swimAngle*pi/180))
    Sx <- as.numeric(swimSpeed * cos(swimAngle*pi/180))
    
    # x and y component of the net vector 
    Ty <- Sy + Cy 
    Tx <- Sx + Cx 
    # polar coordinates
    Tr <- sqrt(Ty^2 + Tx^2) # length of displacement vector
    Tt <- atan2(Ty,Tx)*180/pi # angle of displacement vector (in degrees)
    if (Tt <0) { Tt <- 360 + Tt} # correct for negative angle values produced
    
    # move the larva - IMPORTANT - SWITCH Tx / Ty here because of R angle convention
    xkn <- xk + Ty * dt2 # multiply by number of seconds in timestep
    ykn <- yk + Tx * dt2 # multiply by number of seconds in timestep
    
    # collect the output of the integration step
    output <- c(xk,yk,tide,currentSpeed,currentAngle,swimSpeed,swimAngle,Tr,Tt,xkn,ykn)
    # and write it to the output file
    displacements[i,] <- output
    
    # update position for next run
    xk <- xkn
    yk <- ykn
    
  }
} else if (fullModel=="TRUE"){
  
  # we have no static variables
  
  for (i in 1:dur2){ # we just begin the simulation
    
    if (i%%2 == 1) { # for odd integration steps
      
      # get a current angle from the distribution in the ebb tide data
      currentAngle <- rvonmises(n=1, mu=circular(227,type=c("angles"),units=c("degrees")), kappa=1.487)
      tide <- "ebb"
      
    } else if (i%%2 == 0) { # for even integration steps 
      
      # get a current angle from the distribution in the flood tide data
      currentAngle <- rvonmises(n=1, mu=circular(332.75,type=c("angles"),units=c("degrees")), kappa=1.456)
      tide <- "flood"
      
    }
    
    # generate a current speed from the distribution in the data
    currentSpeed <- rtruncnorm(n=1,a=0,b=Inf,mean=3.9,sd=2.73)
    
    # generate a swimming speed from John's data
    swimSpeed <- rtruncnorm(n=1,a=0,b=Inf,mean=2,sd=0.5)
    
    # genereate a swimming angle (relative to current) from DISC data
    swimAngleRelative <- rvonmises(n=1, mu=circular(187.2,type=c("angles"),units=c("degrees")), kappa=0.7665)
    
    # x and y components of current vector
    Cy <- as.numeric(currentSpeed * sin(currentAngle*pi/180))
    Cx <- as.numeric(currentSpeed * cos(currentAngle*pi/180))
    
    # calculating absolute swimming angle based on relative angle
    swimAngle <- currentAngle - swimAngleRelative
    if (swimAngle < 0) { swimANgle <- 360 + swimAngle}
    
    # x and y components of swimming vector
    Sy <- as.numeric(swimSpeed * sin(swimAngle*pi/180))
    Sx <- as.numeric(swimSpeed * cos(swimAngle*pi/180))
    
    # x and y component of the net vector 
    Ty <- Sy + Cy 
    Tx <- Sx + Cx 
    # polar coordinates
    Tr <- sqrt(Ty^2 + Tx^2) # length of displacement vector
    Tt <- atan2(Ty,Tx)*180/pi # angle of displacement vector (in degrees)
    if (Tt <0) { Tt <- 360 + Tt} # correct for negative angle values produced
    
    # move the larva - IMPORTANT - SWITCH Tx / Ty here because of R angle convention
    xkn <- xk + Ty * dt2 # multiply by number of seconds in timestep
    ykn <- yk + Tx * dt2 # multiply by number of seconds in timestep
    
    # collect the output of the integration step
    output <- c(xk,yk,tide,currentSpeed,currentAngle,swimSpeed,swimAngle,Tr,Tt,xkn,ykn)
    # and write it to the output file
    displacements[i,] <- output
    
    # update position for next run
    xk <- xkn
    yk <- ykn
    
  }
} else if (noSwim == "TRUE") {
  
  # we don't set any swimming parameters
  swimSpeed=0
  swimAngle=0
  
  for (i in 1:dur2){ # then begin the simulation
    
    if (i%%2 == 1) { # for odd integration steps
      
      # get a current angle from the distribution in the ebb tide data
      currentAngle <- rvonmises(n=1, mu=circular(227,type=c("angles"),units=c("degrees")), kappa=1.487)
      tide <- "ebb"
      
    } else if (i%%2 == 0) { # for even integration steps 
      
      # get a current angle from the distribution in the flood tide data
      currentAngle <- rvonmises(n=1, mu=circular(332.75,type=c("angles"),units=c("degrees")), kappa=1.456)
      tide <- "flood"
      
    }
    
    # generate a current speed from the distribution in the data
    currentSpeed <- rtruncnorm(n=1,a=0,b=Inf,mean=3.9,sd=2.73)
    
    # x and y components of current vector
    Cy <- as.numeric(currentSpeed * sin(currentAngle*pi/180))
    Cx <- as.numeric(currentSpeed * cos(currentAngle*pi/180))
    
    # x and y component of the net vector 
    Ty <- Cy 
    Tx <- Cx 
    # polar coordinates
    Tr <- sqrt(Ty^2 + Tx^2) # length of displacement vector
    Tt <- atan2(Ty,Tx)*180/pi # angle of displacement vector (in degrees)
    if (Tt <0) { Tt <- 360 + Tt} # correct for negative angle values produced
    
    # move the larva
    xkn <- xk + Tx * dt2 # multiply by number of seconds in timestep
    ykn <- yk + Ty * dt2 # multiply by number of seconds in timestep
    
    # collect the output of the integration step
    output <- c(xk,yk,tide,currentSpeed,currentAngle,swimSpeed,swimAngle,Tr,Tt,xkn,ykn)
    # and write it to the output file
    displacements[i,] <- output
    
    # update position for next run
    xk <- xkn
    yk <- ykn
    
  }
} else if (randomSwimming=="TRUE"){
  
  # we have no static variables
  
  for (i in 1:dur2){ # we just begin the simulation
    
    if (i%%2 == 1) { # for odd integration steps
      
      # get a current angle from the distribution in the ebb tide data
      currentAngle <- rvonmises(n=1, mu=circular(227,type=c("angles"),units=c("degrees")), kappa=1.487)
      tide <- "ebb"
      
    } else if (i%%2 == 0) { # for even integration steps 
      
      # get a current angle from the distribution in the flood tide data
      currentAngle <- rvonmises(n=1, mu=circular(332.75,type=c("angles"),units=c("degrees")), kappa=1.456)
      tide <- "flood"
      
    }
    
    # generate a current speed from the distribution in the data
    currentSpeed <- rtruncnorm(n=1,a=0,b=Inf,mean=3.9,sd=2.73)
    
    # generate a swimming speed from John's data
    swimSpeed <- rtruncnorm(n=1,a=0,b=Inf,mean=2,sd=0.5)
    
    # genereate a random angle from 0-359
    randA <- sample(0:359,1)
    swimAngleRelative <- circular(randA,type=c("angles"),units=c("degrees"))
    
    # x and y components of current vector
    Cy <- as.numeric(currentSpeed * sin(currentAngle*pi/180))
    Cx <- as.numeric(currentSpeed * cos(currentAngle*pi/180))
    
    # calculating absolute swimming angle based on relative angle
    swimAngle <- currentAngle - swimAngleRelative
    if (swimAngle < 0) { swimANgle <- 360 + swimAngle}
    
    # x and y components of swimming vector
    Sy <- as.numeric(swimSpeed * sin(swimAngle*pi/180))
    Sx <- as.numeric(swimSpeed * cos(swimAngle*pi/180))
    
    # x and y component of the net vector 
    Ty <- Sy + Cy 
    Tx <- Sx + Cx 
    # polar coordinates
    Tr <- sqrt(Ty^2 + Tx^2) # length of displacement vector
    Tt <- atan2(Ty,Tx)*180/pi # angle of displacement vector (in degrees)
    if (Tt <0) { Tt <- 360 + Tt} # correct for negative angle values produced
    
    # move the larva - IMPORTANT - SWITCH Tx / Ty here because of R angle convention
    xkn <- xk + Ty * dt2 # multiply by number of seconds in timestep
    ykn <- yk + Tx * dt2 # multiply by number of seconds in timestep
    
    # collect the output of the integration step
    output <- c(xk,yk,tide,currentSpeed,currentAngle,swimSpeed,swimAngle,Tr,Tt,xkn,ykn)
    # and write it to the output file
    displacements[i,] <- output
    
    # update position for next run
    xk <- xkn
    yk <- ykn
    
  }
} else if (fastSwim=="TRUE"){
  
  # we have no static variables
  
  for (i in 1:dur2){ # we just begin the simulation
    
    if (i%%2 == 1) { # for odd integration steps
      
      # get a current angle from the distribution in the ebb tide data
      currentAngle <- rvonmises(n=1, mu=circular(227,type=c("angles"),units=c("degrees")), kappa=1.487)
      tide <- "ebb"
      
    } else if (i%%2 == 0) { # for even integration steps 
      
      # get a current angle from the distribution in the flood tide data
      currentAngle <- rvonmises(n=1, mu=circular(332.75,type=c("angles"),units=c("degrees")), kappa=1.456)
      tide <- "flood"
      
    }
    
    # generate a current speed from the distribution in the data
    currentSpeed <- rtruncnorm(n=1,a=0,b=Inf,mean=3.9,sd=2.73)
    
    # generate a swimming speed from John's data
    swimSpeed <- rtruncnorm(n=1,a=0,b=Inf,mean=4,sd=0.5)
    
    # genereate a swim angle relative to current
    swimAngleRelative <- rvonmises(n=1, mu=circular(187.2,type=c("angles"),units=c("degrees")), kappa=0.7665)
    
    # x and y components of current vector
    Cy <- as.numeric(currentSpeed * sin(currentAngle*pi/180))
    Cx <- as.numeric(currentSpeed * cos(currentAngle*pi/180))
    
    # calculating absolute swimming angle based on relative angle
    swimAngle <- currentAngle - swimAngleRelative
    if (swimAngle < 0) { swimANgle <- 360 + swimAngle}
    
    # x and y components of swimming vector
    Sy <- as.numeric(swimSpeed * sin(swimAngle*pi/180))
    Sx <- as.numeric(swimSpeed * cos(swimAngle*pi/180))
    
    # x and y component of the net vector 
    Ty <- Sy + Cy 
    Tx <- Sx + Cx 
    # polar coordinates
    Tr <- sqrt(Ty^2 + Tx^2) # length of displacement vector
    Tt <- atan2(Ty,Tx)*180/pi # angle of displacement vector (in degrees)
    if (Tt <0) { Tt <- 360 + Tt} # correct for negative angle values produced
    
    # move the larva - IMPORTANT - SWITCH Tx / Ty here because of R angle convention
    xkn <- xk + Ty * dt2 # multiply by number of seconds in timestep
    ykn <- yk + Tx * dt2 # multiply by number of seconds in timestep
    
    # collect the output of the integration step
    output <- c(xk,yk,tide,currentSpeed,currentAngle,swimSpeed,swimAngle,Tr,Tt,xkn,ykn)
    # and write it to the output file
    displacements[i,] <- output
    
    # update position for next run
    xk <- xkn
    yk <- ykn
    
  }
} else if (surfaceCurrents=="TRUE"){
  
  # we don't set any swimming parameters
  swimSpeed=0
  swimAngle=0
  
  for (i in 1:dur2){ # we just begin the simulation
    
    if (i%%2 == 1) { # for odd integration steps
      
      # get a current angle from the distribution in the ebb tide data
      currentAngle <- rvonmises(n=1, mu=circular(227,type=c("angles"),units=c("degrees")), kappa=1.487)
      tide <- "ebb"
      
    } else if (i%%2 == 0) { # for even integration steps
      
      # get a current angle from the distribution in the flood tide data
      currentAngle <- rvonmises(n=1, mu=circular(332.75,type=c("angles"),units=c("degrees")), kappa=1.456)
      tide <- "flood"
      
    }
    
    # generate a current speed from the distribution in the data
    currentSpeed <- rtruncnorm(n=1,a=0,b=Inf,mean=5.74,sd=2.95)
    
    # x and y components of current vector
    Cy <- as.numeric(currentSpeed * sin(currentAngle*pi/180))
    Cx <- as.numeric(currentSpeed * cos(currentAngle*pi/180))
    
    # x and y component of the net vector
    Ty <- Cy
    Tx <- Cx
    # polar coordinates
    Tr <- sqrt(Ty^2 + Tx^2) # length of displacement vector
    Tt <- atan2(Ty,Tx)*180/pi # angle of displacement vector (in degrees)
    if (Tt <0) { Tt <- 360 + Tt} # correct for negative angle values produced
    
    # move the larva - IMPORTANT - SWITCH Tx / Ty here because of R angle convention
    xkn <- xk + Ty * dt2 # multiply by number of seconds in timestep
    ykn <- yk + Tx * dt2 # multiply by number of seconds in timestep
    
    # collect the output of the integration step
    output <- c(xk,yk,tide,currentSpeed,currentAngle,swimSpeed,swimAngle,Tr,Tt,xkn,ykn)
    # and write it to the output file
    displacements[i,] <- output
    
    # update position for next run
    xk <- xkn
    yk <- ykn
    
  }
} else {cat("Check Parameter List - No Model Type Selected")}





## Process the model output ----------------------------

# format data frame
colnames(displacements) <- c("startX","startY","tide","currentSpeed","currentAngle","swimSpeed","swimAngle","displaceDistance","displaceAngle","endX","endY")

# convert final distances to km
distanceX <- as.numeric(displacements[nrow(displacements),c("endX")]) /100000
distanceY <- as.numeric(displacements[nrow(displacements),c("endY")]) /100000

# total distance
distance <- sqrt(distanceX^2 + distanceY^2)




## Repeated iterations of model runs to produce a dispersal kernel -----------

# set number of iterations
iter <- 1000

# create empty vector to hold distances
distances <- data.frame(matrix(nrow=iter,ncol=3))

for (j in 1:iter){
  
  xk = xi
  yk = yi
  
  if (simple=="TRUE") { 
    swimSpeed = 2 # cm/s
    swimAngleRelative = circular(187.2,type=c("angles"),units=c("degrees"))
    currentSpeed = 3.9 # cm/s
    dneCurrent = circular(227,type=c("angles"),units=c("degrees"))
    dnfCurrent = circular(332.75,type=c("angles"),units=c("degrees"))
    for (i in 1:dur2){ 
      if (i%%2 == 1) { 
        currentAngle <- dneCurrent 
        tide <- "ebb"
      } else if (i%%2 == 0) { 
        currentAngle <- dnfCurrent 
        tide <- "flood"
      }
      Cy <- as.numeric(currentSpeed * sin(currentAngle*pi/180))
      Cx <- as.numeric(currentSpeed * cos(currentAngle*pi/180))
      swimAngle <- currentAngle - swimAngleRelative
      if (swimAngle < 0) { swimANgle <- 360 + swimAngle}
      Sy <- as.numeric(swimSpeed * sin(swimAngle*pi/180))
      Sx <- as.numeric(swimSpeed * cos(swimAngle*pi/180))
      Ty <- Sy + Cy 
      Tx <- Sx + Cx 
      Tr <- sqrt(Ty^2 + Tx^2) 
      Tt <- atan2(Ty,Tx)*180/pi
      if (Tt <0) { Tt <- 360 + Tt} 
      xkn <- xk + Ty * dt2 
      ykn <- yk + Tx * dt2 
      output <- c(xk,yk,tide,currentSpeed,currentAngle,swimSpeed,swimAngle,Tr,Tt,xkn,ykn)
      displacements[i,] <- output
      xk <- xkn
      yk <- ykn
    } 
  } else if (simpleSwim == "TRUE"){
    swimSpeed <- 2 
    swimAngleRelative = circular(187.2,type=c("angles"),units=c("degrees"))
    for (i in 1:dur2){
      if (i%%2 == 1) { 
        currentAngle <- rvonmises(n=1, mu=circular(227,type=c("angles"),units=c("degrees")), kappa=1.487)
        tide <- "ebb"
      } else if (i%%2 == 0) { 
        currentAngle <- rvonmises(n=1, mu=circular(332.75,type=c("angles"),units=c("degrees")), kappa=1.456)
        tide <- "flood"
      }
      currentSpeed <- rtruncnorm(n=1,a=0,b=Inf,mean=3.9,sd=2.73)
      Cy <- as.numeric(currentSpeed * sin(currentAngle*pi/180))
      Cx <- as.numeric(currentSpeed * cos(currentAngle*pi/180))
      swimAngle <- currentAngle - swimAngleRelative
      if (swimAngle < 0) { swimANgle <- 360 + swimAngle}
      Sy <- as.numeric(swimSpeed * sin(swimAngle*pi/180))
      Sx <- as.numeric(swimSpeed * cos(swimAngle*pi/180))
      Ty <- Sy + Cy 
      Tx <- Sx + Cx 
      Tr <- sqrt(Ty^2 + Tx^2) 
      Tt <- atan2(Ty,Tx)*180/pi 
      if (Tt <0) { Tt <- 360 + Tt} 
      xkn <- xk + Ty * dt2 
      ykn <- yk + Tx * dt2 
      output <- c(xk,yk,tide,currentSpeed,currentAngle,swimSpeed,swimAngle,Tr,Tt,xkn,ykn)
      displacements[i,] <- output
      xk <- xkn
      yk <- ykn
    }
  } else if (fullModel=="TRUE"){
    for (i in 1:dur2){
      if (i%%2 == 1) { 
        currentAngle <- rvonmises(n=1, mu=circular(227,type=c("angles"),units=c("degrees")), kappa=1.487)
        tide <- "ebb"
      } else if (i%%2 == 0) { 
        currentAngle <- rvonmises(n=1, mu=circular(332.75,type=c("angles"),units=c("degrees")), kappa=1.456)
        tide <- "flood"
      }
      currentSpeed <- rtruncnorm(n=1,a=0,b=Inf,mean=3.9,sd=2.73)
      swimSpeed <- rtruncnorm(n=1,a=0,b=Inf,mean=2,sd=0.5)
      swimAngleRelative <- rvonmises(n=1, mu=circular(187.2,type=c("angles"),units=c("degrees")), kappa=0.7665)
      Cy <- as.numeric(currentSpeed * sin(currentAngle*pi/180))
      Cx <- as.numeric(currentSpeed * cos(currentAngle*pi/180))
      swimAngle <- currentAngle - swimAngleRelative
      if (swimAngle < 0) { swimANgle <- 360 + swimAngle}
      Sy <- as.numeric(swimSpeed * sin(swimAngle*pi/180))
      Sx <- as.numeric(swimSpeed * cos(swimAngle*pi/180))
      Ty <- Sy + Cy 
      Tx <- Sx + Cx 
      Tr <- sqrt(Ty^2 + Tx^2) 
      Tt <- atan2(Ty,Tx)*180/pi 
      if (Tt <0) { Tt <- 360 + Tt} 
      xkn <- xk + Ty * dt2 
      ykn <- yk + Tx * dt2 
      output <- c(xk,yk,tide,currentSpeed,currentAngle,swimSpeed,swimAngle,Tr,Tt,xkn,ykn)
      displacements[i,] <- output
      xk <- xkn
      yk <- ykn
    }
  } else if (noSwim == "TRUE") {
    swimSpeed=0
    swimAngle=0
    for (i in 1:dur2){ 
      if (i%%2 == 1) { 
        currentAngle <- rvonmises(n=1, mu=circular(227,type=c("angles"),units=c("degrees")), kappa=1.487)
        tide <- "ebb"
      } else if (i%%2 == 0) { 
        currentAngle <- rvonmises(n=1, mu=circular(332.75,type=c("angles"),units=c("degrees")), kappa=1.456)
        tide <- "flood"
      }
      currentSpeed <- rtruncnorm(n=1,a=0,b=Inf,mean=3.9,sd=2.73)
      Cy <- as.numeric(currentSpeed * sin(currentAngle*pi/180))
      Cx <- as.numeric(currentSpeed * cos(currentAngle*pi/180)) 
      Ty <- Cy 
      Tx <- Cx 
      Tr <- sqrt(Ty^2 + Tx^2) 
      Tt <- atan2(Ty,Tx)*180/pi
      if (Tt <0) { Tt <- 360 + Tt}
      xkn <- xk + Ty * dt2 
      ykn <- yk + Tx * dt2 
      output <- c(xk,yk,tide,currentSpeed,currentAngle,swimSpeed,swimAngle,Tr,Tt,xkn,ykn)
      displacements[i,] <- output
      xk <- xkn
      yk <- ykn
    }
  } else if (randomSwimming=="TRUE"){   
    for (i in 1:dur2){ 
      if (i%%2 == 1) { 
        currentAngle <- rvonmises(n=1, mu=circular(227,type=c("angles"),units=c("degrees")), kappa=1.487)
        tide <- "ebb"
      } else if (i%%2 == 0) { 
        currentAngle <- rvonmises(n=1, mu=circular(332.75,type=c("angles"),units=c("degrees")), kappa=1.456)
        tide <- "flood"
      }
      currentSpeed <- rtruncnorm(n=1,a=0,b=Inf,mean=3.9,sd=2.73)
      swimSpeed <- rtruncnorm(n=1,a=0,b=Inf,mean=2,sd=0.5)
      randA <- sample(0:359,1)
      swimAngleRelative <- circular(randA,type=c("angles"),units=c("degrees"))
      Cy <- as.numeric(currentSpeed * sin(currentAngle*pi/180))
      Cx <- as.numeric(currentSpeed * cos(currentAngle*pi/180))
      swimAngle <- currentAngle - swimAngleRelative
      if (swimAngle < 0) { swimANgle <- 360 + swimAngle}
      Sy <- as.numeric(swimSpeed * sin(swimAngle*pi/180))
      Sx <- as.numeric(swimSpeed * cos(swimAngle*pi/180))
      Ty <- Sy + Cy 
      Tx <- Sx + Cx 
      Tr <- sqrt(Ty^2 + Tx^2) 
      Tt <- atan2(Ty,Tx)*180/pi 
      if (Tt <0) { Tt <- 360 + Tt} 
      xkn <- xk + Ty * dt2 
      ykn <- yk + Tx * dt2 
      output <- c(xk,yk,tide,currentSpeed,currentAngle,swimSpeed,swimAngle,Tr,Tt,xkn,ykn)
      displacements[i,] <- output
      xk <- xkn
      yk <- ykn
    }
  }else if (fastSwim=="TRUE"){
    for (i in 1:dur2){ 
      if (i%%2 == 1) { 
        currentAngle <- rvonmises(n=1, mu=circular(227,type=c("angles"),units=c("degrees")), kappa=1.487)
        tide <- "ebb"
      } else if (i%%2 == 0) { 
        currentAngle <- rvonmises(n=1, mu=circular(332.75,type=c("angles"),units=c("degrees")), kappa=1.456)
        tide <- "flood"
      }
      currentSpeed <- rtruncnorm(n=1,a=0,b=Inf,mean=3.9,sd=2.73)
      swimSpeed <- rtruncnorm(n=1,a=0,b=Inf,mean=4,sd=0.5)
      swimAngleRelative <- rvonmises(n=1, mu=circular(187.2,type=c("angles"),units=c("degrees")), kappa=0.7665)
      Cy <- as.numeric(currentSpeed * sin(currentAngle*pi/180))
      Cx <- as.numeric(currentSpeed * cos(currentAngle*pi/180))
      swimAngle <- currentAngle - swimAngleRelative
      if (swimAngle < 0) { swimANgle <- 360 + swimAngle}
      Sy <- as.numeric(swimSpeed * sin(swimAngle*pi/180))
      Sx <- as.numeric(swimSpeed * cos(swimAngle*pi/180))
      Ty <- Sy + Cy 
      Tx <- Sx + Cx 
      Tr <- sqrt(Ty^2 + Tx^2) 
      Tt <- atan2(Ty,Tx)*180/pi 
      if (Tt <0) { Tt <- 360 + Tt} 
      xkn <- xk + Ty * dt2 
      ykn <- yk + Tx * dt2 
      output <- c(xk,yk,tide,currentSpeed,currentAngle,swimSpeed,swimAngle,Tr,Tt,xkn,ykn)
      displacements[i,] <- output
      xk <- xkn
      yk <- ykn
    }
  } else if (surfaceCurrents=="TRUE"){
    swimSpeed=0
    swimAngle=0
    for (i in 1:dur2){ 
      if (i%%2 == 1) { 
        currentAngle <- rvonmises(n=1, mu=circular(227,type=c("angles"),units=c("degrees")), kappa=1.487)
        tide <- "ebb"
      } else if (i%%2 == 0) { 
        currentAngle <- rvonmises(n=1, mu=circular(332.75,type=c("angles"),units=c("degrees")), kappa=1.456)
        tide <- "flood"
      }
      currentSpeed <- rtruncnorm(n=1,a=0,b=Inf,mean=5.74,sd=2.95)
      Cy <- as.numeric(currentSpeed * sin(currentAngle*pi/180))
      Cx <- as.numeric(currentSpeed * cos(currentAngle*pi/180))
      Ty <- Cy
      Tx <- Cx
      Tr <- sqrt(Ty^2 + Tx^2) 
      Tt <- atan2(Ty,Tx)*180/pi 
      if (Tt <0) { Tt <- 360 + Tt} 
      xkn <- xk + Ty * dt2 
      ykn <- yk + Tx * dt2 
      output <- c(xk,yk,tide,currentSpeed,currentAngle,swimSpeed,swimAngle,Tr,Tt,xkn,ykn)
      displacements[i,] <- output
      xk <- xkn
      yk <- ykn
    }
  } else {cat("Check Parameter List - No Model Type Selected")}
  
  # extract particle distance
  colnames(displacements) <- c("startX","startY","tide","currentSpeed","currentAngle","swimSpeed","swimAngle","displaceDistance","displaceAngle","endX","endY")
  distanceX <- as.numeric(displacements[nrow(displacements),c("endX")]) /100000
  distanceY <- as.numeric(displacements[nrow(displacements),c("endY")]) /100000
  distance <- sqrt(distanceX^2 + distanceY^2)
  
  # save distance to a file
  distances[j,1] <- distance
  distances[j,2] <- distanceX
  distances[j,3] <- distanceY
  
}



## Some different plots of the dispersal data --------------------------

# format it
colnames(distances) <- c("distance","distanceX","distanceY")

# plot the pdf
ggplot(distances, aes(x=distance)) + 
  geom_density(size=1) + geom_vline(xintercept=mean(distances$distance),color="blue") + 
  geom_vline(xintercept=median(distances$distance),color="red") + theme_bw() +
  theme(axis.text.x = element_text(size=14),axis.title.x = element_text(size=16)) + xlab("Distance (km)")



# histogram of x dispersal
ggplot(distances, aes(x=distanceX)) + 
  geom_density(size=1) + geom_vline(xintercept=mean(distances$distanceX),color="blue") + 
  geom_vline(xintercept=median(distances$distanceX),color="red") + theme_bw() +
  theme(axis.text.x = element_text(size=14),axis.title.x = element_text(size=16)) + xlab("Distance (km)")


# histogram of y dispersal
ggplot(distances, aes(x=distanceY)) + 
  geom_density(size=1) + geom_vline(xintercept=mean(distances$distanceY),color="blue") + 
  geom_vline(xintercept=median(distances$distanceY),color="red") + theme_bw() +
  theme(axis.text.x = element_text(size=14),axis.title.x = element_text(size=16)) + xlab("Distance (km)")




## Plot a comparison of models on a standard grid ---------

# easy version with points
models <- read.csv("model_distances.csv")

ggplot(models) +
  geom_point(aes(x=distanceX,y=distanceY,color=model),size=2) +
  xlim(-18,0) + ylim(0,3) +
  xlab("X Distance (km)") + ylab("Y Distance (km)")


# complex version with error bars
distancesTot <- read.csv("full_distances_file.csv")

distanceSum <- ddply(distancesTot, ~model, summarize,
                     meanX=mean(distanceX),
                     meanY=mean(distanceY),
                     sdX=sd(distanceX),
                     sdY=sd(distanceY),
                     seX=sd(distanceX)/sqrt(length(distanceX)),
                     seY=sd(distanceY)/sqrt(length(distanceY)))

distanceSum <- distanceSum[-4,]

ggplot(distanceSum) +
  geom_point(aes(x=meanX,y=meanY,color=model),size=3) +
  geom_errorbar(aes(x=meanX,ymin=meanY-seY,ymax=meanY+seY,color=model,width=.3),size=1) +
  geom_errorbarh(aes(x=meanX,y=meanY,xmin=meanX-seX,xmax=meanX+seX,color=model,height=.1),size=1) +
  xlab("W-E Displacement (km)") + ylab("N-S Displacement (km)") +
  # xlim(-18,-10) + ylim(1.5,3.25) +
  theme_bw() +
  theme(axis.text = element_text(size=14)) +
  theme(axis.title = element_text(size=16,face="bold",vjust=0.5)) +
  scale_color_manual(labels = c("Deep Currents", "Surface Currents", "Swimming Behavior"),
                     values = c("#F8766D","#C77CFF","#00BFC4")) +
  labs(color="Model Variation") +
  guides(color=guide_legend(title.theme=element_text(face="bold",angle=0)))





## Plot overlapping dispersal kernels----------

# total dispersal
ggplot() + 
  geom_density(data=distancesFull,aes(x=distance),size=1,color="blue") + 
  geom_vline(xintercept=mean(distancesFull$distance),color="blue") +
  geom_density(data=distancesSurf,aes(x=distance),size=1,color="black") + 
  geom_vline(xintercept=mean(distancesSurf$distance),color="black") +
  theme_bw() +
  theme(axis.text.x = element_text(size=14),axis.title.x = element_text(size=16)) + xlab("Distance (km)")


# for x dispersal
ggplot() + 
  geom_density(data=distancesFull,aes(x=distanceX),size=1,color="blue") + 
  geom_vline(xintercept=mean(distancesFull$distanceX),color="blue") +
  geom_density(data=distancesSurf,aes(x=distanceX),size=1,color="black") + 
  geom_vline(xintercept=mean(distancesSurf$distanceX),color="black") +
  theme_bw() +
  theme(axis.text.x = element_text(size=14),axis.title.x = element_text(size=16)) + xlab("Distance (km)")



# for y dispersal
ggplot() + 
  geom_density(data=distancesFull,aes(x=distanceY),size=1,color="blue") + 
  geom_vline(xintercept=mean(distancesFull$distanceY),color="blue") +
  geom_density(data=distancesSurf,aes(x=distanceY),size=1,color="black") + 
  geom_vline(xintercept=mean(distancesSurf$distanceY),color="black") +
  theme_bw() +
  theme(axis.text.x = element_text(size=14),axis.title.x = element_text(size=16)) + xlab("Distance (km)")




## Some Statistics Comparing Models -----------------

anova(lm(distance~model,data=distancesTot))
anova(lm(distanceX~model,data=distancesTot))
anova(lm(distanceY~model,data=distancesTot))


t = pairwise.t.test(distancesTot$distance,distancesTot$model,p.adj="none")
pairwise.t.test(distancesTot$distanceX,distancesTot$model,p.adj="none")
pairwise.t.test(distancesTot$distanceY,distancesTot$model,p.adj="none")

t.test(distancesTot[distancesTot$model=="18m_currents",c("distance")],distancesTot[distancesTot$model=="9m_currents",c("distance")],alternative=c("two.sided"))
t.test(distancesTot[distancesTot$model=="18m_currents",c("distance")],distancesTot[distancesTot$model=="routine_swimming",c("distance")],alternative=c("two.sided"))




## Produce a model based ONLY on real data ----------------

# FIRST - open "Belize16_Ecology_DISC_Analysis.R" and load all data

# split into the two contexts where they orient against the current
de <- subset(d2,tide=="ebb" & depth..m.=="18" & location=="near")
df <- subset(d2,tide=="flood")


# run a modified model

# set number of iterations
iter <- 1000

# create empty vector to hold distances
distances <- data.frame(matrix(nrow=iter,ncol=3))

for (j in 1:iter){
  
  xk = xi
  yk = yi
  
  if (fullModel=="TRUE"){
    
    
    for (i in 1:dur2){
      
      if (i%%2 == 1) { # for odd integration steps
        
        # sample a random row from the ebb table data
        larva <- sample_n(de,1)
        tide <- "ebb"
        
      } else if (i%%2 == 0) { # for even integration steps 
        
        # sample a random row from the flood table data
        larva <- sample_n(df,1)
        tide <- "flood"
        
      }
      
      # grab the current speed from the row of data
      currentSpeed <- larva$drift_distance * 100 / 25 / 60
      
      # grab the current direction from the row of data
      currentAngle <- larva$drift_direction
      
      # x and y components of current vector
      Cy <- as.numeric(currentSpeed * sin(currentAngle*pi/180))
      Cx <- as.numeric(currentSpeed * cos(currentAngle*pi/180))
      
      # generate a swimming speed from John's data
      swimSpeed <- rtruncnorm(n=1,a=0,b=Inf,mean=2,sd=0.5)
      
      # grab the swimming direction from the row of data
      swimAngle <- larva$mean
      
      # x and y components of swimming vector
      Sy <- as.numeric(swimSpeed * sin(swimAngle*pi/180))
      Sx <- as.numeric(swimSpeed * cos(swimAngle*pi/180))
      
      # x and y component of the net vector 
      Ty <- Sy + Cy 
      Tx <- Sx + Cx 
      # polar coordinates
      Tr <- sqrt(Ty^2 + Tx^2) # length of displacement vector
      Tt <- atan2(Ty,Tx)*180/pi # angle of displacement vector (in degrees)
      if (Tt <0) { Tt <- 360 + Tt} # correct for negative angle values produced
      
      # move the larva - IMPORTANT - SWITCH Tx / Ty here because of R angle convention
      xkn <- xk + Ty * dt2 # multiply by number of seconds in timestep
      ykn <- yk + Tx * dt2 # multiply by number of seconds in timestep
      
      # collect the output of the integration step
      output <- c(xk,yk,tide,currentSpeed,currentAngle,swimSpeed,swimAngle,Tr,Tt,xkn,ykn)
      # and write it to the output file
      displacements[i,] <- output
      
      # update position for next run
      xk <- xkn
      yk <- ykn
      
    }
  } else if (noSwim == "TRUE") {
    
    # we don't set any swimming parameters
    swimSpeed=0
    swimAngle=0
    
    for (i in 1:dur2){
      
      if (i%%2 == 1) { # for odd integration steps
        
        # sample a random row from the ebb table data
        larva <- sample_n(de,1)
        tide <- "ebb"
        
      } else if (i%%2 == 0) { # for even integration steps 
        
        # sample a random row from the flood table data
        larva <- sample_n(df,1)
        tide <- "flood"
        
      }
      
      # grab the current speed from the row of data
      currentSpeed <- larva$drift_distance * 100 / 25 / 60
      
      # grab the current direction from the row of data
      currentAngle <- larva$drift_direction
      
      # x and y components of current vector
      Cy <- as.numeric(currentSpeed * sin(currentAngle*pi/180))
      Cx <- as.numeric(currentSpeed * cos(currentAngle*pi/180))
      
      # x and y components of swimming vector
      Sy <- as.numeric(swimSpeed * sin(swimAngle*pi/180))
      Sx <- as.numeric(swimSpeed * cos(swimAngle*pi/180))
      
      # x and y component of the net vector 
      Ty <- Sy + Cy 
      Tx <- Sx + Cx 
      # polar coordinates
      Tr <- sqrt(Ty^2 + Tx^2) # length of displacement vector
      Tt <- atan2(Ty,Tx)*180/pi # angle of displacement vector (in degrees)
      if (Tt <0) { Tt <- 360 + Tt} # correct for negative angle values produced
      
      # move the larva - IMPORTANT - SWITCH Tx / Ty here because of R angle convention
      xkn <- xk + Ty * dt2 # multiply by number of seconds in timestep
      ykn <- yk + Tx * dt2 # multiply by number of seconds in timestep
      
      # collect the output of the integration step
      output <- c(xk,yk,tide,currentSpeed,currentAngle,swimSpeed,swimAngle,Tr,Tt,xkn,ykn)
      # and write it to the output file
      displacements[i,] <- output
      
      # update position for next run
      xk <- xkn
      yk <- ykn
      
    }
  } else {cat("Check Parameter List - No Model Type Selected")}
  
  # extract particle distance
  colnames(displacements) <- c("startX","startY","tide","currentSpeed","currentAngle","swimSpeed","swimAngle","displaceDistance","displaceAngle","endX","endY")
  distanceX <- as.numeric(displacements[nrow(displacements),c("endX")]) /100000
  distanceY <- as.numeric(displacements[nrow(displacements),c("endY")]) /100000
  distance <- sqrt(distanceX^2 + distanceY^2)
  
  # save distance to a file
  distances[j,1] <- distance
  distances[j,2] <- distanceX
  distances[j,3] <- distanceY
  
}


# format it
colnames(distances) <- c("distance","distanceX","distanceY")

# NOW, repeat this for the two variations of the model here!
# and save them separately
modifiedSwim <- distances
modifiedPass <- distances

# then put a marker on each variation and combine them
modifiedSwim$model <- "Swimming Behavior"
modifiedPass$model <- "Passive"
distancesTot <- rbind(modifiedSwim,modifiedPass)
write.csv(distancesTot,"distances_real_data.csv")

# plot it

distancesTot <- read.csv("distances_real_data.csv")

distanceSum <- ddply(distancesTot, ~model, summarize,
                     meanX=mean(distanceX),
                     meanY=mean(distanceY),
                     sdX=sd(distanceX),
                     sdY=sd(distanceY),
                     seX=sd(distanceX)/sqrt(length(distanceX)),
                     seY=sd(distanceY)/sqrt(length(distanceY)))


ggplot(distanceSum) +
  geom_point(aes(x=meanX,y=meanY,color=model),size=3) +
  geom_errorbar(aes(x=meanX,ymin=meanY-seY,ymax=meanY+seY,color=model,width=.3),size=1) +
  geom_errorbarh(aes(x=meanX,y=meanY,xmin=meanX-seX,xmax=meanX+seX,color=model,height=.1),size=1) +
  xlab("X Distance (km)") + ylab("Y Distance (km)") +
  # xlim(-18,-10) + ylim(1.5,3.25) +
  theme_bw() +
  theme(axis.title = element_text(size=16,face="bold",vjust=0.5)) +
  scale_color_manual(labels = c("Passive", "Swimming Behavior"),
                     values = c("#F8766D","#C77CFF")) +
  labs(color="Model Variation") +
  guides(color=guide_legend(title.theme=element_text(face="bold",angle=0)))




t.test(distancesTot[distancesTot$model=="Swimming Behavior",c("distance")],distancesTot[distancesTot$model=="Passive",c("distance")],alternative=c("two.sided"))


