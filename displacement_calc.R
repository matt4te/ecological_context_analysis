

## Calculate displacement vector for each deployment -----------------------------


# convert total drift to cm/s
d$current_strength <- d$drift_distance * 100 / 25 / 60

# x and y components of swimming vector
Sy <- as.numeric(d$speed.mean * sin(d$mean*pi/180))
Sx <- as.numeric(d$speed.mean * cos(d$mean*pi/180))

# x and y components of current vector
Cy <- as.numeric(d$current_strength * sin(d$drift_direction*pi/180))
Cx <- as.numeric(d$current_strength * cos(d$drift_direction*pi/180))

Dy <- Sy + Cy # x component of displacement vector
Dx <- Sx + Cx # y component of displacement vector

Dh <- sqrt(Dy^2 + Dx^2) # length of displacement vector
Dt <- atan(Dy/Dx)*180/pi # angle of displacement vector (in degrees)


# correct for negative angle values produced
Dt[Dt < 0 & !is.na(Dt)] <-360 + Dt[Dt < 0 & !is.na(Dt)]