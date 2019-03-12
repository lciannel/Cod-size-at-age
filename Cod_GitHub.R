#various GAM models
library(mgcv)
library(maps)
library(mapdata)
library(spacetime)
library(fields)
library(nlme)
library(colorRamps)  ## to create color palletes
library(itsadug)    ##  to create color.legends
library(RColorBrewer)

source('distance.function.R')# Calculate great circle distance given coordinates
source( "vis.gam_COLORS.R")# plot color legend using vis.gam function


#Obtain bathymetry data from NOAA repository: http://maps.ngdc.noaa.gov/viewers/wcs-client/
bathy.dat<-read.table('BeringDepth.txt',sep='')
names(bathy.dat)<-c('lon','lat','depth')
bathy.dat$depth[bathy.dat$depth>0]<-NA#Avoid points above water
head(bathy.dat)
bathy.mat<-matrix(bathy.dat$depth,nrow=length(unique(bathy.dat$lon)),ncol=length(unique(bathy.dat$lat)))[,order(unique(bathy.dat$lat))]
land.mask<-(bathy.mat*0)+1
range(land.mask,na.rm=T)


# IMPORT COD DATA
fishdata<- read.table("Age.txt", header = TRUE) #Age samples
catchdata<- read.table("Catch.txt", header = TRUE) # stations and catch
lengthdata<- read.table("Length.txt", header = TRUE) #length


#Make prediction grid for plotting temperature
lond<-unique(bathy.dat$lon)
latd<-sort(unique(bathy.dat$lat))
predict.grid<-expand.grid(lond,latd)
names(predict.grid)<-c('lon.avg','lat.avg')
for (i in 1: nrow(predict.grid)){
	predict.grid$dist[i]<-min(distance.function(predict.grid$lat[i],predict.grid$lon[i],catchdata$START_LATITUDE,catchdata$START_LONGITUDE))}
head(predict.grid)

# ADDS STATION VARIABLES TO FISHDATA/LENGTHDATA
match <- match(lengthdata$HAULJOIN, catchdata$HAULJOIN)
lengthdata$TEMP <- catchdata$GEAR_TEMPERATURE[match]
lengthdata$LAT <- catchdata$START_LATITUDE[match]
lengthdata$LON <- catchdata$START_LONGITUDE[match]
lengthdata$DAY_MONTH <- catchdata$DAY_MONTH[match]
lengthdata$TIME24HR <- catchdata$TIME24HR[match]
lengthdata$DEPTH <- catchdata$BOTTOM_DEPTH[match]

match <- match(fishdata$HAULJOIN, catchdata$HAULJOIN)
fishdata$TEMP <- catchdata$GEAR_TEMPERATURE[match]
fishdata$LAT <- catchdata$START_LATITUDE[match]
fishdata$LON <- catchdata$START_LONGITUDE[match]
fishdata$COHORT<-fishdata$YEAR-fishdata$AGE
fishdata$DEPTH <- catchdata$BOTTOM_DEPTH[match]
fishdata$STATIONID <- catchdata$STATIONID[match]

catchdata$cpue<-catchdata$WEIGHT/catchdata$DISTANCE_FISHED

#Add average lat and lon by STATIONID to fishdata
lat.avg<-tapply(fishdata$LAT,fishdata$STATIONID,mean)
lon.avg<-tapply(fishdata$LON,fishdata$STATIONID,mean)
index<-match(fishdata$STATIONID,names(lon.avg))
fishdata$lon.avg<-lon.avg[index]
fishdata$lat.avg<-lat.avg[index]

#Find out and add age of individuals in length file by first building an lt-age gam.
#fit a different model each year
years<-sort(unique(fishdata$YEAR))
par(mfrow=c(6,4),mai=c(0,0.1,0,0))
lengthdata$AGE<-NA
#fit annual GAM, make prediction on length file, and show predictions with residuals
for(i in 1:length(years)){
subdata<-fishdata[fishdata$YEAR==years[i],]
age.gam<-gam(AGE~s(LENGTH,by=factor(SEX))+s(LON,LAT,k=10),gamma=1.7,data=subdata)
lengthdata$AGE[lengthdata$YEAR==years[i]]<-predict(age.gam,newdata=lengthdata[lengthdata$YEAR==years[i],])
pred.1<-lengthdata$AGE[lengthdata$YEAR==years[i]&lengthdata$SEX==1]
pred.2<-lengthdata$AGE[lengthdata$YEAR==years[i]&lengthdata$SEX==2]
plot(sort(lengthdata$LENGTH[lengthdata$YEAR==years[i]&lengthdata$SEX==1]),pred.1[order(lengthdata$LENGTH[lengthdata$YEAR==years[i]&lengthdata$SEX==1])],type='l',xlab="",ylab="",col='red',ylim=range(1:10),xlim=c(0,1000))
lines(sort(lengthdata$LENGTH[lengthdata$YEAR==years[i]&lengthdata$SEX==2]),pred.2[order(lengthdata$LENGTH[lengthdata$YEAR==years[i]&lengthdata$SEX==2])],col='blue')
text(300,6,as.character(years[i]),col='red')
points(subdata$LENGTH[subdata$SEX==1], subdata$AGE[subdata$SEX==1],pch=3,cex=0.5,col='red')
points(subdata$LENGTH[subdata$SEX==2], subdata$AGE[subdata$SEX==2]-0.1,pch=3,cex=0.5,col='blue')
}

#Round ages and set age-0 predictions = 1
lengthdata$AGE_ROUND<-round(lengthdata$AGE,0)
lengthdata$AGE_ROUND[lengthdata$AGE_ROUND==0]<-1
table(lengthdata$AGE_ROUND)

#Add FREQUENCY_SUM to catch data: total number of individuals measured
age_groups<-tapply(lengthdata$FREQUENCY,list(lengthdata$HAULJOIN,lengthdata$AGE_ROUND),sum,na.rm=T)
FREQUENCY_SUM<-apply(age_groups,1,sum,na.rm=T)
index<-match(catchdata$HAULJOIN,as.numeric(names(FREQUENCY_SUM)))
catchdata$FREQUENCY_SUM<-NA
catchdata$HAULCHECK<-NA
catchdata$FREQUENCY_SUM<-FREQUENCY_SUM[index]
catchdata$HAULCHECK<-as.numeric(names(FREQUENCY_SUM))[index]

#Exploratory analysis of sample size, for both length and age data
par(mfrow=c(4,1),omi=c(0.1,0.2,0.3,0.1),mai=c(0.2,0.6,0.1,0))
plot(as.numeric(names(tapply(catchdata$cpue,catchdata$YEAR,length))),tapply(catchdata$cpue,catchdata$YEAR,length),xlab="",ylab="Stations sampled",main='',ylim=c(340,400),type='b',cex.axis=1.2,cex.lab=1.2,cex.main=1.5,axes=F)
axis(1,labels=F);axis(2);box()

plot(as.numeric(names(tapply(lengthdata$LENGTH,lengthdata$YEAR,length))),tapply(lengthdata$LENGTH,lengthdata$YEAR,length)/1000,xlab="",ylab="Individuals caught x 1000",main='',ylim=c(6,11),type='b',cex.axis=1.2,cex.lab=1.2,cex.main=1.5,axes=F)
axis(1,labels=F);axis(2);box()

plot(as.numeric(names(tapply(fishdata$LENGTH,fishdata$YEAR,length))),tapply(fishdata$LENGTH,fishdata$YEAR,length),xlab="Year",ylab="Individuals aged",main='',ylim=c(500,2000),type='b',cex.axis=1.2,cex.lab=1.2,cex.main=1.5,axes=F)
axis(1,labels=F);axis(2);box()

plot(as.numeric(names(tapply(fishdata$LENGTH,fishdata$YEAR,length))),round(tapply(fishdata$LENGTH,fishdata$YEAR,length)/tapply(lengthdata$LENGTH,lengthdata$YEAR,length)*100,1),xlab="Year",ylab="% aged",main='',ylim=c(0,25),type='b',cex.axis=1.2,cex.lab=1.2,cex.main=1.5)

#number of stations sampled/year
a<-as.matrix(1*(table(list(fishdata$YEAR,fishdata$STATIONID))>0))
mean(apply(a,1,sum))
range(apply(a,1,sum))

#Fig. 1 of MS: Cod catches on cold and warm years over bottom temperatures and time series of cold pool index
#GAM of bottom temperature for mapping, choose year (e.g., 2003 and 1999)
quartz(width=7,height=7)
par(mfrow=c(2,2))

years<-2003#choose 1999 of 2003

gam.temp<-gam(GEAR_TEMPERATURE~s(START_LONGITUDE,START_LATITUDE),data=catchdata[catchdata$YEAR==years,])
summary(gam.temp)
#2003
#R-sq.(adj) =  0.836   Deviance explained = 84.8%
#GCV score = 0.30825  Scale est. = 0.28439   n = 368
#1999
#R-sq.(adj) =  0.904   Deviance explained = 91.1%
#GCV score = 0.21931  Scale est. = 0.20114   n = 348

#BW version
jet.colors<-colorRampPalette(gray(seq(0.2,1,length=10)))#Defined as jet in myvis.gam
myvis.gam(gam.temp,view=c('START_LONGITUDE','START_LATITUDE'),plot.type="contour",color="jet",too.far=.05,main='',ylab='Latitude',xlab='',zlim=c(-2,7))
fitted.val<-fitted(gam.temp)
gradientLegend(valRange=range(-2,7), color = gray(seq(0.2,1,length=10)), pos = c(-178.2,54.7,-176.7,56.7),coords = TRUE,n.seg = 5, border.col = "black")
#contour(unique(bathy.dat$lon),sort(unique(bathy.dat$lat)),bathy.mat,levels=-200,col='grey',add=T,lwd=2)
map("worldHires",fill=T,col="grey",add=T)
symbols(catchdata$START_LONGITUDE[catchdata$YEAR==years],catchdata$START_LATITUDE[catchdata$YEAR==years],circle=log(catchdata$cpue[catchdata$YEAR==years]+1),add=T,inches=0.05,fg='black',bg=alpha('white',0.7))
#points(fishdata$LON[fishdata$YEAR==years],fishdata$LAT[fishdata$YEAR==years],pch="x")
text(-162,61.5,as.character(years),cex=1.5,col='white')
#STOP: repeat for the second year (2003)

par(new=T)
par(mfrow=c(2,1),mfg=c(2,1),omi=c(0.5,0.1,0.1,0.1),mai=c(0.3,0.9,0.7,0.3))
tmp<-catchdata[catchdata$BOTTOM_DEPTH<100&catchdata$BOTTOM_DEPTH>50,]
cold.pool.index<-tapply(tmp$GEAR_TEMPERATURE,tmp$YEAR,mean,na.rm=T)
plot(as.numeric(names(cold.pool.index)),cold.pool.index,type='b',xlab='Year',ylab=expression(paste("Temperature ("^0,'C)')))

# MODEL 1, Table 1 (reference):
gam_age <- gam(LENGTH ~ s(AGE), data = na.exclude(fishdata[,c('LENGTH','SEX','LAT','LON','TEMP','AGE','COHORT')]))
summary(gam_age)
AIC(gam_age)
# R-sq.(adj) =  0.909   Deviance explained =   91%
#GCV = 3397.8  Scale est. = 3396.6    n = 25213
#276557.3

#Fig. 2
pred<-predict(gam_age,se.fit=T)
resid<-residuals(gam_age)
res.avg<-tapply(resid,subdata$STATIONID,mean)
max.res<-max(abs(res.avg),na.rm=T)
lat.avg<-tapply(subdata$lat.avg,subdata$STATIONID,mean)
lon.avg<-tapply(subdata$lon.avg,subdata$STATIONID,mean)
pred.length<-pred[[1]][order(subdata$AGE)]
pred.se<-pred[[2]][order(subdata$AGE)]

#GAM of bottom temperature for mapping
gam.temp<-gam(GEAR_TEMPERATURE~factor(YEAR)+s(START_LONGITUDE,START_LATITUDE),data=catchdata)
summary(gam.temp)
#R-sq.(adj) =  0.805   Deviance explained = 80.6%
#GCV = 0.73664  Scale est. = 0.73206   n = 8324
dev.new(width=10,height=5)
par(mfrow=c(1,2))
plot(0,-100,ylim=c(0,1200),xlim=c(0,18),lwd=2,col='blue',ylab="Length (mm)",xlab="Age (years)",main="",cex.axis=1.3,cex.lab=1.3,cex.main=1.3)
lines(sort(subdata$AGE),pred.length+1.96*pred.se,lty=2)
lines(sort(subdata$AGE),pred.length-1.96*pred.se,lty=2)
points(subdata$AGE[resid>=0], subdata $LENGTH[resid>=0],pch=16,cex=0.5,fg="black",col='grey')
points(subdata $AGE[resid<0], subdata $LENGTH[resid<0],pch=16,cex=0.5,fg="black",col='black')
lines(sort(subdata $AGE),pred.length,lwd=2)


plot(-1,-1,xlim=range(subdata $lon.avg),ylim=range(subdata $lat.avg),ylab="",xlab="", main="",cex.axis=1.3,cex.lab=1.3,cex.main=1.3)
contour(unique(bathy.dat$lon),sort(unique(bathy.dat$lat)),bathy.mat,levels=-200,col='grey',add=T,lwd=2)
map("worldHires",fill=T,col="grey",add=T)
symbols(lon.avg[res.avg>=0],lat.avg[res.avg>=0],circle=res.avg[res.avg>=0],add=T,inches=0.08*max(res.avg[res.avg>=0],na.rm=T)/max.res,fg='grey',bg='grey')
symbols(lon.avg[res.avg<0],lat.avg[res.avg<0],circle=abs(res.avg[res.avg<0]),add=T,inches=0.08*max(abs(res.avg[res.avg<0]),na.rm=T)/max.res,fg='black',bg='black')
vis.gam(gam.temp,view=c('START_LONGITUDE','START_LATITUDE'),plot.type='contour',add=T,too.far=0.05,color='bw',nCol=0,col='black',levels=c(0,1,2,3),labcex=1.3,lwd=1.5)
symbols(rep(-177.5,4),c(56.5,56,55.5,55),circle=c(85,60,30,10),add=T,inches=0.08*85/max.res,fg='grey10',bg='grey40')
text(rep(-175,4),c(56.5,56,55.5,55),labels=c('85 mm','60 mm','30 mm','10 mm'))


#Fig S5: Map residuals for different age groups separately
dev.new()
par(mfrow=c(3,3))
par(mfrow=c(3,3))
age<-1:9#choose age group here

for(i in 1:length(age)){
	pred<-predict(gam_age,se.fit=T)
	resid<-residuals(gam_age)
	res.avg<-tapply(resid[subdata$AGE==age[i]], subdata$STATIONID[subdata$AGE==age[i]],mean)
	max.res<-max(abs(res.avg),na.rm=T)
	lat.avg<-tapply(subdata$lat.avg[subdata$AGE==age[i]],subdata$STATIONID[subdata$AGE==age[i]],mean)
	lon.avg<-tapply(subdata$lon.avg[subdata$AGE==age[i]],subdata$STATIONID[subdata$AGE==age[i]],mean)
	
	plot(-1,-1,xlim=range(subdata$lon.avg),ylim=range(subdata$lat.avg),ylab="",xlab="", main=paste('Residuals Age-',age[i]))
	contour(unique(bathy.dat$lon),sort(unique(bathy.dat$lat)),bathy.mat,levels=-200,col='grey',add=T,lwd=2)
	map("worldHires",fill=T,col="grey",add=T)
	symbols(lon.avg[res.avg>=0],lat.avg[res.avg>=0],circle=res.avg[res.avg>=0],add=T,inches=0.05*max(res.avg[res.avg>=0],na.rm=T)/max.res,fg='green4',bg='green4')
	symbols(lon.avg[res.avg<0],lat.avg[res.avg<0],circle=abs(res.avg[res.avg<0]),add=T,inches=0.05*max(abs(res.avg[res.avg<0]),na.rm=T)/max.res,fg='red',bg='red')
	vis.gam(gam.temp,view=c('START_LONGITUDE','START_LATITUDE'),plot.type='contour',add=T,too.far=0.05,color='bw',nCol=0,col='black',levels=c(0,1,2,3))}


#Model 2, Table 1: All additive
gam1<-gam(LENGTH~factor(YEAR)+factor(SEX)+s(LON,LAT)+s(TEMP)+s(AGE),data=fishdata)
summary(gam1)
AIC(gam1)
#R-sq.(adj) =  0.929   Deviance explained = 92.9%
#GCV = 2675.5  Scale est. = 2668.3    n = 25213
#270531.4

#Model 3, Table 1: Interaction between age and temperature
gam2<-gam(LENGTH~factor(YEAR)+factor(SEX)+s(LON,LAT)+s(AGE,TEMP),data=fishdata)
summary(gam2)
AIC(gam2)
#R-sq.(adj) =  0.929   Deviance explained =   93%
#GCV = 2655.9  Scale est. = 2647.4    n = 25213
#270346.1


#Model 4, Table 1: Interaction between age and temperature by sex
gam4<- gam(LENGTH ~ factor(YEAR)+ s(LON,LAT) + s(AGE, TEMP, by = factor(SEX)), data = fishdata)
summary(gam4)
AIC(gam4)
#R-sq.(adj) =   0.93   Deviance explained =   93%
#GCV = 2633.6  Scale est. = 2622      n = 25213
#AIC: 270133.1

#Model 5, Table 1: Same as gam4 with COHORT effect
gam4.1<- gam(LENGTH ~ factor(COHORT)+s(LON,LAT) + s(AGE, TEMP, by = factor(SEX)), data = fishdata)
summary(gam4.1)
AIC(gam4.1)
#R-sq.(adj) =  0.927   Deviance explained = 92.8%
#GCV = 2736.3  Scale est. = 2723      n = 25213
#AIC: 271097.8


#Fig. 4: Predict size given age, sex over a range of temperatures
dev.new(width=8,height=4)
newdata<-data.frame(SEX=rep(1,100),AGE=rep(7,100),TEMP=seq(-1,8,length=100),YEAR=rep(2006,100),LAT=rep(58,100),LON=rep(-170,100))
pred<-predict.gam(gam4,newdata=newdata)
plot(newdata$TEMP,pred,type='l',ylim=c(500,900),xlab="Temperature (C)",ylab="Size (mm)",cex.axis=1.2,cex.lab=1.2)

newdata<-data.frame(SEX=rep(2,100),AGE=rep(7,100),TEMP=seq(-1,8,length=100),YEAR=rep(2006,100),LAT=rep(58,100),LON=rep(-170,100))
pred<-predict.gam(gam4,newdata=newdata)
lines(newdata$TEMP,pred,col='grey')

newdata<-data.frame(SEX=rep(1,100),AGE=rep(5,100),TEMP=seq(-1,8,length=100),YEAR=rep(2006,100),LAT=rep(58,100),LON=rep(-170,100))
pred<-predict.gam(gam4,newdata=newdata)
lines(newdata$TEMP,pred,lty=2)

newdata<-data.frame(SEX=rep(2,100),AGE=rep(5,100),TEMP=seq(-1,8,length=100),YEAR=rep(2006,100),LAT=rep(58,100),LON=rep(-170,100))
pred<-predict.gam(gam4,newdata=newdata)
lines(newdata$TEMP,pred,lty=2,col='grey')

newdata<-data.frame(SEX=rep(1,100),AGE=rep(9,100),TEMP=seq(-1,8,length=100),YEAR=rep(2006,100),LAT=rep(58,100),LON=rep(-170,100))
pred<-predict.gam(gam4,newdata=newdata)
lines(newdata$TEMP,pred,lty=3)

newdata<-data.frame(SEX=rep(2,100),AGE=rep(9,100),TEMP=seq(-1,8,length=100),YEAR=rep(2006,100),LAT=rep(58,100),LON=rep(-170,100))
pred<-predict.gam(gam4,newdata=newdata)
lines(newdata$TEMP,pred,lty=3,col='grey')

# FIG. 3 plotting length as a contour plot of position
dev.new(width=10,height=5)
par(mfrow=c(1,2))
vis.gam(gam4, view=c("LON","LAT"), cond=list(SEX="1",AGE=2), plot.type="contour", color='bw', main="", xlab = "longitude", ylab = "latitude", too.far = 0.1,labcex=1.3)
map("worldHires",fill=T,col="grey",add=T)
# symbols(fishdata$LON[male], fishdata$LAT[male], circles = fishdata$LENGTH[male], inches = 0.02, add = T)
vis.gam(gam4, view=c("LON","LAT"), cond=list(SEX="2",AGE=2), plot.type="contour", color="bw", main="", xlab = "longitude", ylab = "latitude", too.far = 0.1,labcex=1.3)
map("worldHires",fill=T,col="grey",add=T)


# Model 6, Table 1: GAMMs with random cohort effect of GAM 4: AIC and StDev of random term indicate that the includion of COHORT is improving model.
gamm1 <- gamm(LENGTH ~ factor(YEAR) + s(AGE, TEMP, by = factor(SEX)) + s(LON,LAT), random = list(COHORT =~ 1),method = "REML", data = fishdata)
summary(gamm1$gam)
summary(gamm1$lme)
#R-sq.(adj) =  0.927   Scale est. = 2570.5    n = 25213
 #AIC      		BIC    		logLik
 # 269860.7 	270161.6 	-134893.3
#         (Intercept) Residual
#StdDev:    13.68149 50.69962


#########BIAS analyses (very long time to run!!!!)################
#Based on Kastelle et al (2018):
# p(right age)  = 0.61
# p(of overestimating age by 1 year) = 0.25
# p(underestimating age by 1 year) = 0.13
# p(overestimating by 2 years) = 0.005
# p(underestimating by 2 years) = 0.005
#Assume that all of this applies to years < 2011, which is when the bias became known within the age and growth lab

#Step 1: randomize age by adding bias
fishdata2<-fishdata
ok<-round(nrow(fishdata[fishdata$YEAR<2012,])*0.61)
ov.1<-round(nrow(fishdata[fishdata$YEAR<2012,])*0.25)
ud.1<-round(nrow(fishdata[fishdata$YEAR<2012,])*0.13)
ov.2<-round(nrow(fishdata[fishdata$YEAR<2012,])*0.005)
ud.2<-round(nrow(fishdata[fishdata$YEAR<2012,])*0.005)
sum(ok,ov.1,ov.2,ud.1,ud.2)==nrow(fishdata[fishdata$YEAR<2012,])#OK if T
bias<-c(rep(0,ok),rep(-1,ov.1),rep(1,ud.1),rep(-2,ov.2),rep(2,ud.2))
#table(bias)

#Step 2: re-fit GAM model and get year effect. Repeat n = 200 times
gam.ref<-gam(LENGTH ~ factor(YEAR)+ s(AGE)-1, data = fishdata)
nd<-data.frame(YEAR=unique(fishdata2$YEAR),AGE=rep(3,length(unique(fishdata2$YEAR))))
n<-200
year.effect<-matrix(NA,nrow=n,ncol=length(unique(fishdata2$YEAR)))
for(i in 1:nrow(year.effect)){
fishdata2$AGE[fishdata2$YEAR<2012]<-fishdata$AGE[fishdata$YEAR<2012]+sample(bias)
gam4.bias<- gam(LENGTH ~ factor(YEAR)+ s(AGE)-1, data = fishdata2)
year.effect[i,]<-predict(gam4.bias,newdata=nd)}

#Step 3: plot results (Fig. S3) 
dev.new(width=8,height=4)
x<-unique(fishdata2$YEAR)
avg<-apply(year.effect,2,mean)
up<-avg+apply(year.effect,2,sd)
lw<-avg-apply(year.effect,2,sd)
ref<-predict(gam.ref,newdata=nd)
plot(x, avg,ylim=range(c(avg, up, lw, ref)),pch=20, xlab="Year", ylab="Length (mm)",main="",cex.lab=1.4,cex.axis=1.4)
arrows(x, lw, x, up, length=0.05, angle=90, code=3)
points(x,ref,pch=3,col='red')
legend(2012,400,legend=c('reference','bootstrap'),col=c('red','black'),pch=c(3,20))
dev.copy(jpeg,'Bias_fig.jpg',width=8,height=4,units='in',res=100)
dev.off()

#Step 4: Calculate mean bias for the perios 1994-2012
mean(predict(gam.ref,newdata=nd)[unique(fishdata2$YEAR)<2012]-avg[unique(fishdata2$YEAR)<2012])

########COHORTS MODELS#######################
#Build new data frame, called age_data, with age-specific cpue

age_fraction<-age_groups/FREQUENCY_SUM
index<-match(as.numeric(names(FREQUENCY_SUM)),catchdata$HAULJOIN)
age_cpue<-age_fraction*catchdata$cpue[index]
age_data<-expand.grid(dimnames(age_cpue)[[1]],dimnames(age_cpue)[[2]])
names(age_data)<-c('HAULJOIN','AGE')
age_data$cpue<-as.vector(age_cpue)
age_data$AGE<-as.numeric(age_data$AGE)
index<-match(age_data$HAULJOIN,catchdata$HAULJOIN)
age_data$cpue[is.na(age_data$cpue)]<-0
age_data$START_LATITUDE<-catchdata$START_LATITUDE[index]
age_data$START_LONGITUDE<-catchdata$START_LONGITUDE[index]
age_data$YEAR<-catchdata$YEAR[index]
age_data$STATIONID<-catchdata$STATIONID[index]
age_data$AGEPLUS<-ifelse(age_data$AGE<=8,age_data$AGE,'9+')
age_data$BIRTHYEAR<-age_data$YEAR-age_data$AGE
age_data$GEAR_TEMPERATURE<-catchdata$GEAR_TEMPERATURE[index]

#Add average lat and lon by STATIONID
lat.avg<-tapply(age_data$START_LATITUDE,age_data$STATIONID,mean)
lon.avg<-tapply(age_data$START_LONGITUDE,age_data$STATIONID,mean)
index<-match(age_data$STATIONID,names(lon.avg))
age_data$lon.avg<-lon.avg[index]
age_data$lat.avg<-lat.avg[index]



#Fit a GAM model and make a image plot of predictions for each age group (Fig S1)
dev.new(width=8,height=7)
par(mfrow=c(3,3),mai=c(0.3,0.3,0.08,0.2))
ages<-unique(age_data$AGEPLUS)
for(i in 1:length(ages)){
	gam.age<-gam(log(cpue+1)~factor(YEAR)+s(START_LONGITUDE,START_LATITUDE),data=age_data[age_data$AGEPLUS==ages[i],])
	vis.gam(gam.age,view=c('START_LONGITUDE','START_LATITUDE'),plot.type="contour",color="topo",too.far=.05,main='',ylab='Latitude',xlab='Longitude')
symbols(age_data$START_LONGITUDE[age_data$AGEPLUS==ages[i]&age_data$cpue>0],age_data$START_LATITUDE[age_data$AGEPLUS==ages[i]&age_data$cpue>0],circles=log((age_data$cpue+1)[age_data$AGEPLUS==ages[i]&age_data$cpue>0]),inches=0.03,bg='grey',fg='black',add=T)
map("worldHires",fill=T,col="grey",add=T)
contour(unique(bathy.dat$lon),sort(unique(bathy.dat$lat)),bathy.mat,levels=-c(200),labcex=0.4,col='black',add=T,)
rsq<-round(summary(gam.age)$r.sq,3)
text(-172.8,55,labels=paste('Age',"-",ages[i]),cex=1)
text(-173.0,55.8,labels=expression(paste('R'^2,' = ','         %')),cex=1)
text(-172.1,55.73,labels=eval(round(rsq*100,1)),cex=1)
}


#For the rest of the cohort analyses, trim data to only include individuals from ages from 1-5, and remove cohorts that are not fully represented within through this age range (i.e., from 1-5)
subdata<-na.exclude(age_data[age_data$AGE<=5,])
cohort_summary<-table(subdata$BIRTHYEAR,subdata$AGE)
index<-apply(1*(cohort_summary==0),1,sum)==0
keepcohort<-sort(unique(subdata$BIRTHYEAR))[index]
index<-subdata$BIRTHYEAR%in%keepcohort
age_data5<-subdata[index,]

#Calculate average temperature, and average position (lat and lon, age-1 only), both weighted by cpue, experienced by each age group within each year (i.e., cohort throughout its life)
tmp1.t<-age_data5$cpue*age_data5$GEAR_TEMPERATURE
tmp1.lat<-age_data5$cpue*age_data5$START_LATITUDE
tmp1.lon<-age_data5$cpue*age_data5$START_LONGITUDE
tmp2.t<-tapply(tmp1.t,list(age_data5$YEAR,age_data5$AGE),sum)/tapply(age_data5$cpue,list(age_data5$YEAR,age_data5$AGE),sum)
tmp2.lat<-tapply(tmp1.lat,list(age_data5$YEAR,age_data5$AGE),sum)/tapply(age_data5$cpue,list(age_data5$YEAR,age_data5$AGE),sum)
tmp2.lon<-tapply(tmp1.lon,list(age_data5$YEAR,age_data5$AGE),sum)/tapply(age_data5$cpue,list(age_data5$YEAR,age_data5$AGE),sum)

cohort<-expand.grid(dimnames(tmp2.t))
cohort$temp<-c(tmp2.t)
cohort$lat<-c(tmp2.lat)
cohort$lon<-c(tmp2.lon)

age_data5$year.age<-paste(age_data5$YEAR,age_data5$AGE,sep='.')
index<-match(age_data5$year.age,paste(cohort[,1],cohort[,2],sep="."))
age_data5$t.age<-cohort$temp[index]
age_data5$lat.age<-cohort$lat[index]
age_data5$lon.age<-cohort$lon[index]


#Fig. S2: map of age-1 distribution for each cohort and center of distribution overlaid on image of bottom temperature 

dev.new(width=12,height=8)
par(mfcol=c(5,4),mai=c(0.23,0.28,0.1,0.23))
for(i in 1:length(keepcohort)){
bt.loess<-loess(GEAR_TEMPERATURE~lon.avg*lat.avg,span=0.08,degree=2,data=age_data5[age_data5$BIRTHYEAR==keepcohort[i]&age_data5$AGE==1,])
bt.pred<-predict(bt.loess,newdata=predict.grid)
bt.pred[predict.grid$dist>30000]<-NA
image(lond,latd,bt.pred,col=tim.colors(100),zlim=c(-2.5,10),ylab='',xlab='',ylim=range(age_data5$START_LATITUDE),xlim=range(age_data5$START_LONGITUDE))
symbols(age_data5$START_LONGITUDE[age_data5$BIRTHYEAR==keepcohort[i]&age_data5$AGE==1&age_data5$cpue>0],age_data5$START_LATITUDE[age_data5$BIRTHYEAR==keepcohort[i]&age_data5$AGE==1&age_data5$cpue>0],circles=log((age_data5$cpue+1)[age_data5$BIRTHYEAR==keepcohort[i]&age_data5$AGE==1&age_data5$cpue>0]),inches=0.025,bg='grey',fg='black',add=T)
points(age_data5$lon.age[age_data5$BIRTHYEAR==keepcohort[i]&age_data5$AGE==1],age_data5$lat.age[age_data5$BIRTHYEAR==keepcohort[i]&age_data5$AGE==1],pch=8,cex=2,col='red')
map("worldHires",fill=T,col="grey",add=T)
contour(unique(bathy.dat$lon),sort(unique(bathy.dat$lat)),bathy.mat,levels=-c(200),labcex=0.4,col='black',add=T,)
text(-176.1,55.7,labels=paste('Age',"-",1),cex=1)
text(-173,55,labels=paste('BY',"=",keepcohort[i],'-','CY','=',I(keepcohort[i]+1)),cex=1)
	}
par(new=T,mfcol=c(5,4),mfg=c(5,4))
plot(1,1,axes=F,pch="")
points(0.62,0.9,pch=8,col='red',cex=2)
text(0.77,0.9,labels='Centroid',cex=1.4)
text(0.78,0.76,labels='BY= Birth Year',cex=1.4)
text(0.83,0.64,labels='CY= Calendar Year',cex=1.4)
par(new=T,mfcol=c(5,4),mfg=c(5,4),omi=c(1.2,0.25,0,0.1))
image.plot(legend.only=T,zlim=c(-2.5,10),col=tim.colors(100),horizontal=T,legend.lab='Bottom temperature (C)',legend.width=3,legend.cex=0.9)


#Make new data file ('name'),  with average temperature from age-0 to age_target
#For age-0 temperature use average middle shelf temp during birth year
#For age-1+ temperature, use length file to find out ambient temperature

tmp3<-tapply(age_data5$t.age,list(age_data5$BIRTHYEAR,age_data5$AGE),mean)
name<-cbind(expand.grid(as.numeric(dimnames(tmp3)[[1]]),as.numeric(dimnames(tmp3)[[2]])),c(tmp3))
names(name)<-c('COHORT','AGE','T.AGE')
name$YEAR<-name$COHORT+name$AGE
name$T.AGE_minus_1[name$AGE>=2]<-name$T.AGE[name$AGE<=4]
name$T.AGE1<-NA;name$T.AGE2<-NA;name$T.AGE3<-NA;name$T.AGE4<-NA

name$T.AGE_minus_1[name$AGE>=2]<-name$T.AGE[name$AGE<=4]
name$T.AGE1[name$AGE==1]<-name$T.AGE[name$AGE==1]
name$T.AGE1[name$AGE==2]<-name$T.AGE[name$AGE==1]
name$T.AGE1[name$AGE==3]<-name$T.AGE[name$AGE==1]
name$T.AGE1[name$AGE==4]<-name$T.AGE[name$AGE==1]
name$T.AGE1[name$AGE==5]<-name$T.AGE[name$AGE==1]
name$T.AGE2[name$AGE==1]<-name$T.AGE[name$AGE==2]
name$T.AGE2[name$AGE==2]<-name$T.AGE[name$AGE==2]
name$T.AGE2[name$AGE==3]<-name$T.AGE[name$AGE==2]
name$T.AGE2[name$AGE==4]<-name$T.AGE[name$AGE==2]
name$T.AGE2[name$AGE==5]<-name$T.AGE[name$AGE==2]
name$T.AGE3[name$AGE==3]<-name$T.AGE[name$AGE==3]
name$T.AGE3[name$AGE==4]<-name$T.AGE[name$AGE==3]
name$T.AGE3[name$AGE==5]<-name$T.AGE[name$AGE==3]
name$T.AGE4[name$AGE==4]<-name$T.AGE[name$AGE==4]
name$T.AGE4[name$AGE==5]<-name$T.AGE[name$AGE==4]
name$T.AGE5[name$AGE==5]<-name$T.AGE[name$AGE==5]

#Import data with average temperature in middle shelf area (50-100m), which also includes pre-1994 years 
st_full<-read.table('Station_Data_1992_2012.csv',sep=',',header=T)

T.avg<-tapply(st_full$BOTTOM_TEMP_CELCIUS[st_full$DEPTH_METERS>50&st_full$DEPTH_METERS<100],st_full$YEAR[st_full$DEPTH_METERS>50&st_full$DEPTH_METERS<100],mean,na.rm=T)
name$T.AGE0<-T.avg[match(name$COHORT,as.numeric(names(T.avg)))]
name$T.AGE_minus_1[name$AGE==1]<-name$T.AGE0[name$AGE==1]


#Average position of age-1 groups from age_data5 (based on length) file, weighted by cpue: proxy for spawning origin of cohort of origin
tmp4<-tapply(age_data5$lat.age,list(age_data5$BIRTHYEAR,age_data5$AGE),mean)
index<-match(name$COHORT,as.numeric(dimnames(tmp4)[[1]]))
name$LAT.AGE1<-tmp4[,1][index]
tmp4<-tapply(age_data5$lon.age,list(age_data5$BIRTHYEAR,age_data5$AGE),mean)
index<-match(name$COHORT,as.numeric(dimnames(tmp4)[[1]]))
name$LON.AGE1<-tmp4[,1][index]

#Calculate the least square line that passes through the center of distribution of age-1 fish
lon1<-tapply(name$LON.AGE1,name$COHORT,mean)
lat1<-tapply(name$LAT.AGE1,name$COHORT,mean)
lm.pos<-lm(lat1~lon1)
summary(lm.pos)
#Residual standard error: 0.2208 on 17 degrees of freedom
#Multiple R-squared:  0.2842,	Adjusted R-squared:  0.2421 
#F-statistic: 6.749 on 1 and 17 DF,  p-value: 0.01877

#Calculate shortest distance of each cohort average position at age-1 from the least square line
lond<-seq(min(lon1),max(lon1),length=100)
latd<-coef(lm.pos)[1]+coef(lm.pos)[2]*lond
#Plot least square line
plot(lon1,lat1,ylim=range(fishdata$LAT),xlim=range(fishdata$LON),ylab="Latitude",xlab="Longitude",main="Position at age-1")
map("worldHires",fill=T,col="grey",add=T)
abline(lm.pos)
points(lond,latd,pch='.',col='red')
#Calculate shortest distance from origin of least square line to projection point of each cohort on the line
cohort.dist<-lon1*NA
for(i in 1:length(lon1)){
	line.index<-order(distance.function(lat1[i],lon1[i],latd,lond))[1]
	cohort.dist[i]<-distance.function(latd[line.index],lond[line.index],max(latd),max(lond))/1000
}
points(lon1[order(cohort.dist)[1]],lat1[order(cohort.dist)[1]],pch='+',col='blue')
points(lon1[order(cohort.dist,decreasing=T)[1]],lat1[order(cohort.dist,decreasing=T)[1]],pch='+',col='green')

#Add new column to the 'name' data frame with the cohort distance from the origin of the least square line
index<-match(name$COHORT,as.numeric(names(cohort.dist)))
name$cohort.dist<-cohort.dist[index]



#Update the fishdata file with the new covariates: Tage-1, position, cohort distance
index<-match(paste(fishdata$COHORT,fishdata$AGE),paste(name$COHORT,name$AGE))
fishdata$T.AGE1<-name$T.AGE1[index]
fishdata$LON.AGE1<-name$LON.AGE1[index]
fishdata$cohort.dist<-name$cohort.dist[index]

#Limit cohort analyses up to age5, only complete cohorts
tmp<-fishdata[fishdata$SEX<3&fishdata$AGE<6,]
cohort_summary<-table(tmp$COHORT,tmp$AGE)
index<-apply(1*(cohort_summary==0),1,sum)==0
keepcohort<-sort(unique(tmp$COHORT))[index]
index<-tmp$COHORT%in%keepcohort
subdata<-tmp[index,]


#COHORTS MODELS: Various shades, shown in table 2
#GAM fit of model with year vs cohort temp effect (Models 1 and 2 in Table 2)
#Year model: Model 1, Table 2
gam1<-gam(LENGTH~factor(YEAR)+factor(SEX)+s(LON,LAT)+s(TEMP,k=3)+s(AGE,k=3),data=subdata)
summary(gam1)
#R-sq.(adj) =    0.9   Deviance explained =   90%
#GCV =   2283  Scale est. = 2275.2    n = 16236
AIC(gam1)
#171634.6

#Cohort model: Model 2, Table 2
gam2<-gam(LENGTH~factor(COHORT)+factor(SEX)+s(LON,LAT)+s(TEMP,k=3)+s(AGE,k=3),data=subdata)
summary(gam2)
#R-sq.(adj) =  0.897   Deviance explained = 89.7%
#GCV = 2352.7  Scale est. = 2345.3    n = 16236
AIC(gam2)
#172122.9

#Plot year and cohort effects (Fig. 5)
dev.new(width=10,height=8)
par(mfrow=c(2,1))
plot.gam(gam1,all.term=T,select=4,scale=0,cex.axis=1.4,cex.lab=1.4,main='Year',ylabs='Length anomalies (mm)',xlabs='Year sampled')
plot.gam(gam2,all.term=T,select=4,scale=0,cex.axis=1.4,cex.lab=1.4,main='Cohort',ylabs='Length anomalies (mm)',xlabs='Year of birth')
dev.copy(jpeg,'Fig5.jpg',width=10,height=8,units='in',res=100)
dev.off()

#Calculate cohort composition for each year
tmp1<-table(list(subdata$COHORT,subdata$YEAR))
tmp.sum<-apply(tmp1,2,sum)
cohort.comp<-tmp1*NA
for(i in 1:ncol(tmp1)){cohort.comp[,i]<-round(tmp1[,i]/tmp.sum[i]*100,3)}
#Check (OK if 100)
apply(cohort.comp,2,sum)


#Assign weights to each year based on the estimated cohort effect of gam2
tmp<-na.exclude(subdata[,c('YEAR','LENGTH','COHORT','SEX','LON','LAT','TEMP','AGE')])
tmp.1<-predict.gam(gam2,type='terms')[,1]
tmp2<-tapply(tmp.1,tmp$COHORT,mean)
cohort.effect<-max(abs(tmp2))+tmp2

year.cohort<-as.numeric(names(tmp.sum))*NA
for(i in 1:ncol(cohort.comp)){
year.cohort[i]<-sum(cohort.comp[,i]*cohort.effect)/100
}

#Correlation of cohort effect with YEAR effect from Model 1 (Figure 6)
tmp<-na.exclude(subdata[,c('YEAR','LENGTH','COHORT','SEX','LON','LAT','TEMP','AGE')])
tmp.1<-predict.gam(gam1,type='terms')[,1]
year.effect<-tapply(tmp.1,tmp$YEAR,mean)

dev.new(width=10,height=5)
par(mfrow=c(1,2),mai=c(1,0.9,0.5,0.2))
image.plot(as.numeric(unlist(dimnames(cohort.comp)[1])),as.numeric(unlist(dimnames(cohort.comp)[2])),cohort.comp,ylab='Year sampled',xlab='Year of birth',col=gray((40:0)/40),smallplot= c(.75,0.8,0.2,.9))

plot(year.cohort,year.effect,xlab='Weighted cohort effect',ylab='Estimated year effect')
lm1<-lm(year.effect~year.cohort)
summary(lm1)
#Residual standard error: 11.26 on 21 degrees of freedom
#Multiple R-squared:  0.751,	Adjusted R-squared:  0.7391 
#F-statistic: 63.33 on 1 and 21 DF,  p-value: 8.964e-08
abline(lm1)

#Correlation between year effect and cold pool
tmp<-catchdata[catchdata$BOTTOM_DEPTH<100&catchdata$BOTTOM_DEPTH>50,]
cold.pool.index<-tapply(tmp$GEAR_TEMPERATURE,tmp$YEAR,mean,na.rm=T)

plot(cold.pool.index,year.effect,xlab='Cold Pool Temp (C)',ylab='Estimated year effect')
lm2<-lm(year.effect~cold.pool.index)
summary(lm2)
#Residual standard error: 22.24 on 21 degrees of freedom
#Multiple R-squared:  0.02763,	Adjusted R-squared:  -0.01867 
#F-statistic: 0.5968 on 1 and 21 DF,  p-value: 0.4484
abline(lm2)
text(0,0,paste("Rsq = ",as.character(round(summary(lm2)$adj.r.squared,2))),cex=1.5)


#Model with T.AGE_1: Model 3, Table 2
gam8<-gam(LENGTH~factor(SEX)+s(LON,LAT)+s(TEMP,k=3)+s(AGE,k=3)+s(T.AGE1,k=3),data=subdata)
summary(gam8)

#R-sq.(adj) =  0.889   Deviance explained = 88.9%
#GCV = 2520.9  Scale est. = 2515.5    n = 16236
AIC(gam8)
#173244.4

#Model with cohort.dist: Model 4, Table 2
gam10.1<-gam(LENGTH~factor(SEX)+s(LON,LAT)+s(TEMP,k=3)+s(AGE,k=3)+s(cohort.dist,k=3),data=subdata)
summary(gam10.1)
plot(gam10.1,pages=1,all.terms=T,scale=0)
#R-sq.(adj) =  0.886   Deviance explained = 88.7%
#GCV = 2586.3  Scale est. = 2580.8    n = 16236
AIC(gam10.1)
#173659.9


#Model with cohort.dist + T.AGE1: Model 5, Table 2
gam12.1<-gam(LENGTH~factor(SEX)+s(LON,LAT)+s(TEMP,k=3)+s(AGE,k=3)+s(cohort.dist,k=3)+s(T.AGE1,k=3),data=subdata)
summary(gam12.1)
#R-sq.(adj) =   0.89   Deviance explained =   89%
#GCV = 2506.5  Scale est. = 2500.8    n = 16236
AIC(gam12.1)
#173151




####Threshold models: Model 6, Table 2
#Additive effects of T.AGE1, and with threshold effect of cohort.dist (Best resonable model)
dist.age1<-tapply(subdata$cohort.dist,subdata$COHORT,mean)
th.values<-seq(quantile(dist.age1,0.2),quantile(dist.age1,0.8),length=20)
aic.dist1<-th.values*NA
for(i in 1:length(th.values)){
	gam.test<-gam(LENGTH~factor(SEX)+s(LON,LAT)+s(TEMP,k=3)+s(AGE,k=3)+s(T.AGE1,k=3)+as.factor(cohort.dist<th.values[i]),data=subdata)
	aic.dist1[i]<-AIC(gam.test)
}

#Fig. S3
dev.new(height=10,width=6)
par(mfrow=c(2,1))
lond<-seq(min(lon1),max(lon1),length=100)
latd<-coef(lm.pos)[1]+coef(lm.pos)[2]*lond
plot(tapply(name$LON.AGE1,name$COHORT,mean),tapply(name$LAT.AGE1,name$COHORT,mean),ylim=range(subdata$LAT),xlim=range(subdata$LON),ylab="Latitude",xlab="Longitude",main="Position at age-1")
map("worldHires",fill=T,col="grey",add=T)
points(name$LON.AGE1[name$cohort.dist<164.4003&name$AGE==1],name$LAT.AGE1[name$cohort.dist<164.4003&name$AGE==1],pch=16, col='blue')
points(name$LON.AGE1[name$cohort.dist>=164.4003&name$AGE==1],name$LAT.AGE1[name$cohort.dist>=164.4003&name$AGE==1],pch=16, col='red')
lines(c(min(lond),max(lond)),c(max(latd),min(latd)),col='black',lwd=2)

plot(th.values,aic.dist1,type='b',xlab='Distance (km)',ylab='AIC')
th.best<-th.values[order(aic.dist1)[1]]
abline(v=th.best,lty=2)

#Fit threshold GAM model with best distance
subdata$cohort.loc<-ifelse(subdata$cohort.dist<th.best,'East','West')
gam23.7<-gam(LENGTH~as.factor(cohort.loc)+factor(SEX)+s(LON,LAT)+s(TEMP,k=3)+s(AGE,k=3)+s(T.AGE1,k=3)-1,data=subdata)
summary(gam23.7)
#R-sq.(adj) =  0.891   Deviance explained = 89.1%
#GCV = 2489.7  Scale est. = 2484.1    n = 16236
AIC(gam23.7)
#173041.6

#Plot position and T.AGE1 effects estimated from gam23.7 (Fig 7 in paper)
dev.new(width=10,height=5)
par(mfrow=c(1,2))
plot(gam23.7,select=5,all.terms=T,scale=0,xlabs='',ylabs='Average length (mm)',cex.axis=1.6,cex.lab=1.6)
plot.gam(gam23.7,select=4,scale=0,ylab='Length anomalies (mm)',xlab='Temperature (C)', cex.axis=1.6,cex.lab=1.6,lwd=2)


