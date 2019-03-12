distance.function<-function(start.lat, start.lon, end.lat, end.lon)
{
med.lat <- (start.lat + end.lat)/2
rad.lat <- (pi * med.lat)/180
shrink <- cos(rad.lat)
delta.lat <- end.lat - start.lat
delta.lon <- start.lon - end.lon
mpermile <- 111195
distance <- mpermile * sqrt((delta.lon * shrink)^2 + (delta.lat)^2)
distance
}