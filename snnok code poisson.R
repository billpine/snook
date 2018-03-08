#code to read FWC SAS files with water quality data

setwd("C:/Users/billpine/Google Drive/GEBF Lone Cabbage Oyster/Fish/Snook")
#setwd("C:/Users/Owner/Google Drive/GEBF Lone Cabbage Oyster/Fish/Snook")

#Note from Caleb
#Zone B is north and Zone C is south

library(sas7bdat)
library(MASS)
library(coefplot)
library(plyr)

install.packages("R2admb")
install.packages("glmmADMB", 
                 repos=c("http://glmmadmb.r-forge.r-project.org/repos",
                         getOption("repos")),
                 type="source")

library(glmmADMB)



dat=read.sas7bdat("ckm_cu_c.sas7bdat")
names(dat)

plot(dat$month, dat$number, pch=as.integer(dat$Zone))

snook=subset(dat,dat$month>3 & dat$month<11 & dat$Zone=="C") #only snook months 3-10 in Zone C
winter=dat[dat$month %in% c(1,2,11,12), ] #only temps from these months


## Observed frequency of occurrence of snook
snook$success=ifelse(snook$number>0,1,0)
mean_occur=aggregate(success~year,data=snook,mean) #note this is the mean of the successes
total_occur=aggregate(number~year,data=snook,sum) #note this is the sum catch
years=as.matrix(1997:2016)
colnames(years)=c("year")
freq_occur=merge(years,mean_occur,all.x=TRUE,all.Y=TRUE)
freq_occur=merge(mean_occur,total_occur,all.x=TRUE,all.Y=TRUE)

#summer water temps same as snook months
summer_temp=aggregate(temperature~year,data=snook,mean)

#winter water temps
winter_temp=aggregate(temperature~year,data=winter,mean)
#winter_temp=merge(freq_occur,winter_temp,all.x=TRUE,all.Y=TRUE)

#all water temps 
all_temp=aggregate(temperature~year,data=dat,mean)

par(mfrow = c(2,2))
plot(all_temp$temperature~all_temp$year)
plot(winter_temp$temperature~winter_temp$year)
plot(summer_temp$temperature~summer_temp$year)



#merge mean snook catch, winter teamp, annual temp
snook_temp=merge(freq_occur,all_temp,by="year",all.x=TRUE,all.Y=TRUE)
snook_temp=merge(snook_temp,winter_temp, by="year",all.x=TRUE,all.Y=TRUE)
snook2=merge(snook_temp,summer_temp, by="year",all.x=TRUE,all.Y=TRUE)


colnames(snook2)=c('year','mean_success','total_catch','all_t','winter_t', 'summer_t')


#read in relative snook abundance from snook assessment 1997-2014

GOM_N<-c(1195493,1197976,1455324,1735345,1486287,1414318,2118084,1950108,
                1587182,1330583,1127938,958978,927989,1031855,1228955,2011836,2772778,3581340)
max(GOM_N) 
year<-c(1997:2014)
rel_N<-GOM_N/max(GOM_N)
plot(year,rel_N)
aa<-cbind(year, rel_N)

snook3=merge(snook2,aa, by="year",all.x=TRUE,all.Y=TRUE)

names(snook3)


colnames(snook3)=c('year','mean_success','total_catch','all_t','winter_t', 'summer_t', 'rel_N')


###########
###########
###########



library(lattice)
library(MASS)

names(snook2)


#make a simple graph snook catch Mar-Oct vs. Summer temp in Zone C
xyplot(snook2$total_catch~snook2$summer_t, xlab= "Summer temperature", ylab="Snook catch Mar-Oct, Zone C")
xyplot(snook3$total_catch~snook3$rel_N, xlab= "Relative GOM Snook N", ylab="Snook catch Mar-Oct, Zone C")

plot((snook3$total_catch/max(snook3$total_catch))~snook3$year, xlab= "Year", ylab="Relative Snook catch Mar-Oct, Zone C")
lines(snook3$rel_N~snook3$year, col='red', type='p')

#Fit simple linear model
M0<-lm(total_catch~summer_t, data=snook2)
summary(M0)
#model isn't significant

#Get residuals and fitted values
E0 <- resid(M0)
F0 <- fitted(M0)

#residuals and fitted regression
par(mfrow = c(1,2), mar = c(5,5,3,2))
plot(x = F0, 
     y = E0,
     xlab = "Fitted values",
     ylab = "Residuals",
     cex.lab = 1.5)
abline(h = 0, lty = 2)

plot(x = snook2$summer_t, 
     y = snook2$total_catch,
     xlab = "Summer Temp",
     ylab = "Snook catch",
     cex.lab = 1.5,
     pch = 16)
abline(M0, lwd = 5)





#Predicted and observed
#Linear model doesn't do a good job. negative fitted values

par(mfrow = c(1,1), mar = c(5,5,3,2))
plot(x = snook2$summer_t, 
     y = snook2$total_catch,
     xlab = "Summer Temperature",
     ylab = "Snook catch",
     cex.lab = 1.5,
     pch = 1,
     ylim = c(-50, 100))
abline(M0, lwd = 5)
abline(h = 0, lty = 2)

range(snook2$summer_t)
md <- seq(24, 30, length = 10)

Beta <- coef(M0)
for (i in 1:10){
  mu <- Beta[1] + Beta[2] * md[i]
  yi <- rnorm(100, mean = mu, sd = summary(M0)$sigma)
  points(jitter(rep(md[i], 100)), jitter(yi), col = grey(0.5), pch = 16, cex = 1)
}

#now switch to Poisson GLM

M1 <- glm(total_catch~summer_t, 
          data = snook2, 
          family = poisson(link = "log"))
summary(M1)


#This is observed snook catch vs summer temp with the fitted Poisson GLM line
#improves over linear regression as it doesn't go to zero
par(mar = c(5,5,2,2))
MyData <- data.frame(summer_t = seq(24, 30, length = 25))
P1 <- predict(M1, newdata = MyData, type = "response")
plot(x = snook2$summer_t,
     y = snook2$total_catch,
     ylim = c(0,100),
     xlab = "Summer temp",
     ylab = "Total snook catch", cex.lab = 1.5)

lines(MyData$summer_t, P1, lwd = 3)



#ok can simulate a bunch of values from a Poisson and then plot those and see how it fits

HL <- seq(24, 30, length = 25)
Beta <- coef(M1)
for (i in 1:25){
  mu <- exp(Beta[1] + Beta[2] * HL[i])
  yi <- rpois(50, lambda= mu)
  points(jitter(rep(HL[i], 50)), 
         jitter(yi), col = grey(0.5), 
         pch = 16, cex = 1)
}

lines(MyData$summer_t, P1, lwd = 3)

#it fits ok, not great, but simulated range is similar to the actual range of data

#
#check for overdispersion
E1<-resid(M1, type="pearson")
N<-nrow(snook2)
p<-length(coef(M1))
Dispersion <- sum(E1^2) / (N - p)
Dispersion
#high overdispersion


##negative binomial

M2 <- glm.nb(total_catch~summer_t, 
          data = snook2)
summary(M2)

#check for overdispersion
E2<-resid(M2, type="pearson")
N<-nrow(snook2)
p2<-length(coef(M2))+1
Dispersion.nb <- sum(E2^2) / (N - p2)
Dispersion.nb
#much lower dispersion


##ok just added snook catch, need to re-run everything as snook3 now. need more covariates
#pscl package

M3 <- glm.nb(total_catch~summer_t+rel_N-1, 
             data = snook3, maxit=300)
summary(M3)

#check for overdispersion
E2<-resid(M2, type="pearson")
N<-nrow(snook2)
p2<-length(coef(M2))+1
Dispersion.nb <- sum(E2^2) / (N - p2)
Dispersion.nb
#much lower dispersion
