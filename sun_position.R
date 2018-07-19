#
# Compute the sun position relative to a point in earth
#
# Implemented based on
# 	http://answers.google.com/answers/threadview/id/782886.html
# and
#	http://www.mathworks.com/matlabcentral/fileexchange/4605
#
# (c) Copyright 2009 Jean-Olivier Irisson. GNU General Public License
#
#------------------------------------------------------------

# Input from http://answers.google.com/answers/threadview/id/782886.html
# date = as.POSIXct("2006-11-15 10:35:00")
# lon = 2.28972
# lat = 48.81667

rad <- function(x) { x*pi/180 }
deg <- function(x) { x*180/pi }

sun.position <- function(date, lat, lon)
#
#	Compute the azimuth (bearing) and zenith (angle from vertical) of the sun in degrees at a given time in a given place
# 	date		**GMT** date and time given as a POSIXct object (tz setting does not matter) or a character string that can be coerced to POSIXct (has format %Y-%M-%D %h:%m:%s)
#	lat,lon	position in decimal degrees of latitude and longitude
#
{
	if (! "POSIXct" %in% class(date) ) {
		date = as.POSIXct(date)
	}

	# Check input length
	nDate = length(date)
	nLat = length(lat)
	nLon = length(lon)
	if (nDate != nLat || nLat != nLon) {
		stop("Input date, lat and lon must of be the same size")
	}

	# Remove NAs
	d = data.frame(date, lat, lon)
	OKidx = which(apply(d, 1, function(x) !any(is.na(x))))
	date = date[OKidx]
	lat = lat[OKidx]
	lon = lon[OKidx]

	# Fractional year in degrees (g)
	julDate = julian(date, origin=as.POSIXct(paste(format(date[1],"%Y"),"-01-01 00:00:00",sep="")))
	julDate = as.integer(julDate) + 1
	fracHour = as.numeric(difftime(date, trunc(date, units="days"), units="hours"))
	g = (360/365.25)*(julDate + fracHour/24)
	# convert to radians
	gr = rad(g)

	# Declination of the sun
	D = 0.396372 - 22.91327*cos(gr) + 4.02543*sin(gr) - 0.387205*cos(2*gr) + 0.051967*sin(2*gr) - 0.154527*cos(3*gr) + 0.084798*sin(3*gr)
	Dr = rad(D)

	# Time correction for solar angle
	TC = 0.004297 + 0.107029*cos(gr) - 1.837877*sin(gr) - 0.837378*cos(2*gr) - 2.340475*sin(2*gr)

	# Solar Hour Angle (SHA)
	SHA = (fracHour-12)*15 + lon + TC
	SHAr = rad(SHA)

	# Sun Zenith Angle (SZA):
	latr = rad(lat)
	cosSZA = sin(latr)*sin(Dr)+cos(latr)*cos(Dr)*cos(SHAr)
	SZAr = acos(cosSZA)
	SZA = deg(SZAr)

	# Sun Elevation Angle or Altitude(SEA)
	#SEA = 90 - SZA

	# Sun Azimuth Angle (AZ)
	# # google answer formula
	# cosAZ = (sin(Dr)-sin(latr)*cos(SZAr))/(cos(latr)*sin(SZAr))
	# AZr = acos(cosAZ)
	# AZ = deg(AZr)
	# # but this always gives a result between 0 and 180 i.e. does not distinguish between morning and afternoon
	
	# other formula from the MATLAB code
	AZr = atan2(sin(SHAr), cos(SHAr)*sin(latr)-tan(Dr)*cos(latr)) + pi
	# NB: the +pi conversion is to pass from astronomer notation (westward from south) to navigation notation (eastward from north)
	AZ = deg(AZr)
	# keep between 0 and 360
	AZ = AZ %% 360

	res = data.frame(zenith=rep(NA, nDate), azimuth=NA)
	res$zenith[OKidx] = SZA
	res$azimuth[OKidx] = AZ
	return(res)
}
