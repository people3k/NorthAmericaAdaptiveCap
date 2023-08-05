library(rcarbon)
library(ggplot2)
library(nlme)
library(cowplot)
library(plyr)
library(tidyr)
###Controlled comparison of Central Texas and the Middle Mississippi River Valley

######Box Plots and Violin Plots for Isotopes in Central Texas
d2 <- read.csv("data/CTexIsotope.csv")
dCIII<-subset(d2,CompID1=="Component III")

median(dCIII$DeltaN)

pcar2 <- ggplot(dCIII, aes(factor(CompID3), (DeltaC13Car), fill=BRMIndex))+
  geom_violin()+
  geom_boxplot(notch = FALSE, width=0.2)+
  stat_boxplot(geom ='errorbar')+
  scale_fill_gradient(low ="#F8766D", high = "#619CFF" ) +
  #scale_fill_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  facet_wrap( ~ factor(CompID1))+
  scale_y_continuous(limits=c(-15.5, -4.5))+
  #scale_fill_manual(values=c("#FC4E07", "#00A4CCFF"))+
  theme_bw() +
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28),  plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Time period cal BP", y="Delta 13C apatite", title = "Delta 13C Apatite and Midden Index")
pcar2


pn2 <- ggplot(dCIII, aes(factor(CompID3), (DeltaN), fill=BRMIndex))+
  geom_violin()+
  #geom_boxplot(notch = FALSE, width=0.2)+
  stat_summary(fun=median, geom="point", size=4, color="black")+
  scale_fill_gradient(low ="#F8766D", high = "#619CFF", limits=c(0,.6),  breaks = c(.2, .4, .6)) +
  stat_boxplot(geom ='errorbar')+
  scale_y_continuous(limits=c(6, 13.5))+
  #scale_fill_manual(values=c("#FC4E07", "#00A4CCFF"))+
  theme_bw() +
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28),  plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Time period cal BP", y="Delta 15N collagen", title = "F. Central Texas Delta 15N and Midden Index",  fill = "BRM Index")
pn2

nitro1<-subset(dCIII, CompID3 %in% c('2','3'))
wilcox.test(DeltaN ~ CompID3, data=nitro1)

nitro2<-subset(dCIII, CompID3 %in% c('1','3'))
wilcox.test(DeltaN ~ CompID3, data=nitro2)

nitro3<-subset(dCIII, CompID3 %in% c('1','4'))
wilcox.test(DeltaN ~ CompID3, data=nitro3)

nitro4<-subset(dCIII, CompID3 %in% c('2','4'))
wilcox.test(DeltaN ~ CompID3, data=nitro4)

nitro5<-subset(dCIII, CompID3 %in% c('3','4'))
wilcox.test(DeltaN ~ CompID3, data=nitro5)


##Central Texas Fire-Cracked Rock Midden Analysis===============================================

###########Run cluster analysis on burned rock feature surface area

library(ggplot2)
library(splines)
library(cluster)
library(factoextra)
#BASIC GEOMETRIC DATA ANALYSIS WORKFLOW:
#read in data from csv and make r dataframe
d2 <- read.csv("data/ClusterSize2.csv")
#convert first column to rownames
#rownames(d2) <- d2$bp
#remove the supurflous column of names
#d2$CalBP <- NULL
d2$Site<- NULL
d2$Feature<- NULL
d2$SurfaceArea<- NULL
d2$MedCal<- NULL
d2$CompID<-NULL
d2$PhaseID<-NULL
d2$RockNumber<-NULL
d2$RockWeight<-NULL
###Assess optimal number of clusters
set.seed(12345)
fviz_nbclust(d2, pam, method = "silhouette")+
  theme_classic()
####Run cluster analysis on optimal number
set.seed(12345)
pam.res <- pam(d2, 2)
print(pam.res)
####Add clusters to original data
dd2 <- cbind(d2, cluster = pam.res$cluster)
head(dd2, n = 2)
write.table(dd2, file = "ClustersSurface.csv", sep = ",")

####Visualize cluster results
SurArea2<- read.csv("data/SurfaceMeans2.csv")
SurArea3<-subset(SurArea2, TimeID %in% c('1','2','3', '4'))

sur<-ggplot(SurArea3, aes(x=log(SurfaceArea), fill=factor(SurfaceClust))) +
  geom_histogram(color="white")+
  facet_wrap(~TimeID)+
  scale_y_continuous(limits=c(0,10), breaks=c(0,2,4,6,8,10))+
  theme_bw() +
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28),  plot.title = element_text(size=18, face = "bold"))+
  labs(x = "ln Surface area", y="Frequency", title = "Central Texas Feature Cluster Frequency",  fill = "BRM Cluster ID")
sur

  #geom_hline(data=SurArea3, aes(yintercept=grp.length, color=factor(SurfaceClust)),
            # linetype="dashed")

#Central Texas Radiocarbon Analysis=========================================================================
box<- read.csv("data/FinalRCDTexas3.csv")
box2<- subset(box, Region=="CTx")

###Calibrate the radiocarbon ages
cptcal <- calibrate(x = box2$Age,  errors = box2$Error, calCurves = "intcal20",  normalised = FALSE)
boxbins <- binPrep(sites = box2$Trinomial, ages = box2$Age, h = 100)


####Run analysis for component 3 logistic 3500 to 150
spd.CTx <- spd(cptcal, bins=boxbins, runm=150, timeRange=c(3500,200))
plot(spd.CTx, runm=150, xlim=c(3500,200), type="simple")

PrDens<-spd.CTx$grid$PrDens
calBP<-spd.CTx$grid$calBP

######Four parameter logistic model
#Plot null model logistic
logFit3 <- nls(PrDens~SSfpl(calBP, A, B, xmid, scale), data=spd.CTx$grid,control=nls.control(maxiter=700),start=list(A=0.12, B=0.35, xmid=1711, scale=-100))
# Generate a data frame containing the fitted values
logFitDens3=data.frame(calBP=spd.CTx$grid$calBP,PrDens=SSfpl(input=spd.CTx$grid$calBP,A=coefficients(logFit3)[1],B=coefficients(logFit3)[2],xmid=coefficients(logFit3)[3],scal=coefficients(logFit3)[4]))
# Use the modelTest function (returning the raw simulation output - see below)
LogNull3 <- modelTest(cptcal, errors=box2$Error,nsim=100,
                      timeRange=c(3500, 200), model="custom",predgrid=logFitDens3, bins=boxbins, runm=100, raw=TRUE)

##Plot 4 parameter logistic over data and return logfit3 model summary
plot(LogNull3, xlim = c(3500,250))
lines(Time,predict(logFit3),col="red",lty=2,lwd=3)
summary(logFit3)

##Global significance test of the fit of the 4 parameter logistic model
round(LogNull3$pval,4) #p-value

Pre1<-predict(logFit3)

##Data.Frame of simulation
Con<-cbind(LogNull3$result,predict(logFit3))
write.table(Con, file = "data/CTxPhase3ConfitLogfit.csv", sep = ",",col.names=NA)


Con3<- read.csv("data/CTxPhase3ConfitLogfit.csv")

p3 <- ggplot(Con3,aes(x=(calBP), y=(PrDens))) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70") +
  geom_point(aes(), size=2) +
  #scale_color_gradient(low ="#F8766D", high = "#619CFF" ) +
  geom_line(aes(y=logFit3), color="blue", size=1) +
  theme_bw() +
  scale_x_reverse(breaks=c(3500, 2500, 1500, 500), limits=c(3500,250))+
  # scale_y_continuous(limits=c(0,0.42))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Years Cal BP", y="SPD of radiocarbon ages", title = "B. Central Texas Population")+
  geom_vline(xintercept = 2000)+
  geom_vline(xintercept = 1200)+
  geom_vline(xintercept = 650)
p3

####Sum into 30 year intervals==========================================================

###Sum SPD by 30 year intervals
library(zoo)
# sum and save new csvs.
out50 <- rollapply(Con3,30,(sum),by=30,by.column=TRUE,align='right')
write.table(out50, file = "sums/CTex30.csv", sep = ",", col.names=NA)

CTex2<-read.csv(file="sums/CTex30.csv", header=T)

ctex <- ggplot(CTex2,aes(x=(calBP), y=(PerCap))) +
  #geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70") +
  geom_point(aes(), size=1) +
  #scale_color_gradient(low ="#F8766D", high = "#619CFF") +
  #scale_color_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  #geom_line(aes(y=logFit3), color="blue", size=1) +
  theme_bw() +
  scale_x_reverse(breaks=c(3500, 2500, 1500, 500), limits=c(3500,250))+
  # scale_y_continuous(limits=c(-.3,0.38))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Years Cal BP", y="SPD per capita growth", title = "B. Central Texas SPD Per Capita Growth")+
  geom_line(size=1.25)+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 2000)+
  geom_vline(xintercept = 1200)+
  geom_vline(xintercept = 650)
ctex


##KDE=============================================
US.randates = sampleDates(cptcal, bins=boxbins, nsim=200,verbose=FALSE)
D.ckde = ckde(US.randates,timeRange=c(3500,200),bw=50, normalised = FALSE)
plot(D.ckde,type='multiline')
D.ckde$timeRange

##Write matrix of KDEs as a data frame
Check<-as.data.frame(D.ckde$res.matrix)

##Cbind spd and KDEs and write
dd<-cbind(calBP,PrDens, Check)
write.table(dd, file = "data/CTexKDE50.csv", sep = ",", col.names=NA)

###At this point, I open the spreadsheet, and I add a two new columns: MKDE and PerCap.
#MKDE I calculate as the average of the KDE values at each time step. So, in this case, I add the KDE at cal BP 4000 from V1-V200
#and then divide by 200. This gives the mean KDE value at 4000 cal BP. Once I have done this for all time points, I calculate
#per capita growth as LN(MKDE at t+1) -LN(MKDE at time t). Once this is done, I load the csv file back into R.

#Load KDE data frame with calculated mean KDE and per capita growth
dd2<- read.csv("data/CTexKDE50.csv")
#I then trimm time points with missing values at the end due to KDE bandwidth.
dd2<-na.omit(dd2)

###transform dd2 into long table format to plot all of the simulated KDEs
dd3<-dd2 %>% pivot_longer(cols=c(V1:V200),names_to='Runs', values_to='KDE')

##OK, now we fit a logistic model to the mean KDE. This can be a bit tricky due to the start parameters.
#These need to be adjusted based on the MKDE's values. In this case, the start values work well.
#fit logistic to KDE
logFitKD <- nls(MKDE~SSfpl(calBP, A, B, xmid, scale), data=dd2,control=nls.control(maxiter=700),start=list(A=.00005, B=0.00015, xmid=1711, scale=-100))
summary(logFitKD)
predict<-predict(logFitKD)

##Combine predicted logistic values and KDE table for plotting
dd4<-cbind(dd2, predict)

###Simple plot
plot(dd2$MKDE~dd2$calBP)
lines(dd2$calBP, dd2$MKDE,col="green",lty=2,lwd=3)
lines(dd2$calBP,predict(logFitKD),col="red",lty=2,lwd=3)

###Plot mean KDE and all of the simulated KDE models using ggplot
p1ct<-ggplot(data = dd3) +
  geom_line(aes(x = calBP, y = (KDE*100), group = Runs), color = "purple", alpha = 0.5) +
  geom_point(data = dd4, aes(x = calBP, y = (MKDE*100)), color="black", size=2) +
  geom_line(data = dd2, aes(x = calBP, y = (predict*100)), linetype = "dashed", color="green", size=1.5) +
  scale_x_reverse(breaks=c(3500, 2500, 1500, 500), limits=c(3500,300))+
  theme_bw() +
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Years Cal BP", y="KDE of radiocarbon ages", title = "B. Central Texas KDE and Logistic Model")+
  # annotate("text", x =3500, y = .047, label = "Phase 1", size = 7, angle=90)+
  #annotate("text", x =2000, y = .047, label = "Phase 2", size = 7, angle=90)+
  #annotate("text", x =1020, y = .047, label = "Phase 3", size = 7, angle=90)+
  #annotate("text", x =210, y = .047, label = "Phase 4", size = 7, angle=90)+
  #geom_vline(xintercept = 1746)+
  #geom_vline(xintercept = 1904)+
  # facet_wrap(factor(Region)~.)
  geom_vline(xintercept = 2000)+
  geom_vline(xintercept = 1200)+
  geom_vline(xintercept = 650)
p1ct

dd2c<- read.csv("data/CTexKDE50.csv")
#I then trimm time points with missing values at the end due to KDE bandwidth.
dd2c<-na.omit(dd2c)
dd2c<-subset(dd2c, select=-c(V1:V200))
###Sum SPD by 30 year intervals
library(zoo)
# sum and save new csvs.
out50 <- rollapply(dd2c,30,(sum),by=30,by.column=TRUE,align='right')
write.table(out50, file = "sums/CTexKDE30.csv", sep = ",", col.names=NA)

pc30ct<- read.csv("sums/CTexKDE30.csv")

ctexpc <- ggplot(pc30ct,aes(x=(calBP), y=(PerCap))) +
  #geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70") +
  geom_point(aes(), size=1) +
  #scale_color_gradient(low ="#F8766D", high = "#619CFF") +
  #scale_color_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  #geom_line(aes(y=logFit3), color="blue", size=1) +
  theme_bw() +
  scale_x_reverse(breaks=c(3500, 2500, 1500, 500), limits=c(3500,300))+
  scale_y_continuous(limits=c(-.2,0.25))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Years Cal BP", y="KDE per capita growth", title = "D. Central Texas KDE Per Capita Growth")+
  geom_line(size=1.25)+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 2000)+
  geom_vline(xintercept = 1200)+
  geom_vline(xintercept = 650)
ctexpc



################Middle Mississippi==============================================================================

##Starting with the isotope data.  Load isotope data.
d5 <- read.csv("data/Midwest2.csv")

median(d5$delta15n)

#Box plot of isotopes by cultural historical periods.

pn1 <- ggplot(d5, aes(factor(TimeID), (delta15n), fill=(CultIndex)))+
  geom_violin()+
  stat_summary(fun=median, geom="point", size=4, color="black")+
  #geom_boxplot(notch = FALSE, width=0.2)+
  stat_boxplot(geom ='errorbar')+
  scale_fill_gradient(low ="#F8766D", high = "#619CFF", limits=c(0,.9),  breaks = c(.2, .4, .6, .8)) +
  theme_bw() +
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28),  plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Time period cal BP", y="Delta 15N collagen", title = "E. Miss. R. Valley Delta 15N and Cultivation Index", fill = "Cult. Index")
#facet_wrap(~factor(PeriodBP))
pn1

pcol1 <- ggplot(d5, aes(factor(TimeID), (delta13ccol), fill=(CultIndex)))+
  geom_violin()+
  stat_summary(fun=median, geom="point", size=4, color="black")+
  #geom_boxplot(notch = FALSE, width=0.2)+
  stat_boxplot(geom ='errorbar')+
  scale_fill_gradient(low ="#F8766D", high = "#619CFF", limits=c(0,.9),  breaks = c(.2, .4, .6, .8)) +
  theme_bw() +
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28),  plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Time period cal BP", y="Delta 13C collagen", title = "Miss. R. Valley Delta 13C Collagen and Cultivation Index", fill = "Cult. Index")
#facet_wrap(~factor(PeriodBP))
pcol1

#Subset data for Mann-Whittney U tests and run the tests.

nitro1<-subset(d5, TimeID %in% c('2','3'))
wilcox.test(delta15n ~ TimeID, data=nitro1)

nitro2<-subset(d5, TimeID %in% c('1','3'))
wilcox.test(delta15n ~ TimeID, data=nitro2)

nitro3<-subset(d5, TimeID %in% c('1','4'))
wilcox.test(delta15n ~ TimeID, data=nitro3)

nitro4<-subset(d5, TimeID %in% c('2','4'))
wilcox.test(delta15n ~ TimeID, data=nitro4)

nitro5<-subset(d5, TimeID %in% c('3','4'))
wilcox.test(delta15n ~ TimeID, data=nitro5)

nitro6<-subset(d5, TimeID %in% c('4','5'))
wilcox.test(delta15n ~ TimeID, data=nitro6)

##Cultivation index analysis:

#Load data and subset by inclusion of sites from Illinois and Missouri
dc <- read.csv("data/Milner2023CultIndex.csv")
d3c<-subset(dc, Inclusion2=="1")

#Calculate median cultivation index by cultural historical periods.
meansC1<-aggregate((d3c$CultIndex), list(d3c$PeriodID), FUN=median)
meansC1

###run Mann-Whittney U for cultivation index
dct1<-subset(d3c, PeriodID2 %in% c('1','2'))
dct2<-subset(d3c, PeriodID2 %in% c('2','3'))
dct3<-subset(d3c, PeriodID2 %in% c('3','4'))
dct4<-subset(d3c, PeriodID2 %in% c('4','5'))

wilcox.test(CultIndex ~ PeriodID2, data=dct1)
wilcox.test(CultIndex ~ PeriodID2, data=dct2)
wilcox.test(CultIndex ~ PeriodID2, data=dct3)
wilcox.test(CultIndex ~ PeriodID2, data=dct4)

###Box plots of Cultivation Index by components
pc1 <- ggplot(d3c, aes(factor(PeriodID2), (CultIndex), fill=factor(PeriodID2)))+
  # geom_violin()+
  geom_boxplot(notch = FALSE, width=0.2)+
  stat_boxplot(geom ='errorbar')+
  #stat_summary(fun=median, geom="point", size=4, color="black")+
  #facet_wrap( ~ factor(CompID1))+
  #scale_fill_manual(values=c("#FC4E07", "#00A4CCFF"))+
  theme_bw() +
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28),  plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Time period ID", y="Cultivation index", title = "Cultivation Index in the Middle Mississippi River")
# facet_wrap(~factor(PeriodBP))
pc1

###Box plot comparing sits included and not included in the analysis

pn1 <- ggplot(dc, aes(factor(PeriodID2), (CultIndex), fill=factor(PeriodID2)))+
  geom_violin()+
  geom_boxplot(notch = FALSE, width=0.2)+
  stat_boxplot(geom ='errorbar')+
  facet_wrap( ~ factor(Inclusion2))+
  #scale_fill_manual(values=c("#FC4E07", "#00A4CCFF"))+
  theme_bw() +
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28),  plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Time period ID", y="Cultivation index", title = "Cultivation Index of Excluded and Included Sites")
pn1



###Mississippi Radiocarbon Analysis=======================================================
# Middle Mississippi River Valley
SPD.1<-read.csv(file="data/UpMiss.csv", header=T)

cptcal2 <- calibrate(x = SPD.1$Age,  errors = SPD.1$Error, calCurves = "intcal20",  normalised = FALSE)
bin<- binPrep(sites = SPD.1$SitID, ages = SPD.1$Age, h = 100)


####Run analysis for component 3 logistic 3400 to 150
spd <- spd(cptcal2, bins=NA, runm=150, timeRange=c(3500,200))
plot(spd, runm=100, xlim=c(3500, 200), type="simple")


PrDens<-spd$grid$PrDens
calBP<-spd$grid$calBP

######Four parameter logistic model
#Plot null model logistic
logFit3 <- nls(PrDens~SSfpl(calBP, A, B, xmid, scale), data=spd$grid,control=nls.control(maxiter=1000),start=list(A=2, B=0.25, xmid=1700, scale=-300))
# Generate a data frame containing the fitted values
logFitDens3=data.frame(calBP=spd$grid$calBP,PrDens=SSfpl(input=spd$grid$calBP,A=coefficients(logFit3)[1],B=coefficients(logFit3)[2],xmid=coefficients(logFit3)[3],scal=coefficients(logFit3)[4]))
# Use the modelTest function (returning the raw simulation output - see below)
LogNull3 <- modelTest(cptcal2, errors=SPD.1$Error,nsim=100,
                      timeRange=c(3500, 200), model="custom",predgrid=logFitDens3, bins=bin, runm=150, raw=TRUE)

##Plot model
plot(LogNull3, xlim = c(3500,250))
lines(calBP,predict(logFit3),col="red",lty=2,lwd=3)
summary(logFit3)

plot(PrDens~calBP)
lines(calBP, PrDens,col="green",lty=2,lwd=3)
lines(calBP,predict(logFit3),col="red",lty=2,lwd=3)

##Global significance test of the fit of the exponential model
round(LogNull3$pval,4) #p-value

Pre1<-predict(logFit3)
residuals(logFit3)

plot((Pre1)~Time)

##Data.Frame of simulation
Con<-cbind(LogNull3$result,predict(logFit3))
write.table(Con, file = "data/UpMissPhase3ConfitLogfit150.csv", sep = ",")

miss<-read.csv(file="data/UpMissPhase3ConfitLogfit150.csv", header=T)

pmiss <- ggplot(miss,aes(x=(calBP), y=(PrDens))) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70") +
  geom_point(aes(),size=3) +
  # scale_color_gradient(low ="#F8766D", high = "#619CFF" ) +
  #scale_color_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  geom_line(aes(y=logFit3), color="blue", size=1) +
  theme_bw() +
  scale_x_reverse(breaks=c(3500, 2500, 1500, 500), limits=c(3500,300))+
  # scale_y_continuous(limits=c(0,0.42))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Years Cal BP", y="SPD of radiocarbon ages", title = "A. Miss. R. Valley Population")+
  geom_vline(xintercept = 1900)+
  geom_vline(xintercept = 1250)+
  geom_vline(xintercept = 900)+
  geom_vline(xintercept = 750)
#geom_vline(xintercept = 1300)+
#geom_vline(xintercept = 650)
pmiss

Fig5Rev<-plot_grid(upmiss,ctex, ncol=2, align="hv", axis = "rl")
Fig5Rev

pdf("SI/PCCtexMiss.pdf", width=10, height=12.55)
Fig5Rev
dev.off()

###Sum SPD by 30 year intervals
library(zoo)
# sum and save new csvs.
out50 <- rollapply(miss,30,(sum),by=30,by.column=TRUE,align='right')
write.table(out50, file = "sums/Miss30.csv", sep = ",", col.names=NA)

Miss2<-read.csv(file="sums/Miss30.csv", header=T)

upmiss <- ggplot(Miss2,aes(x=(calBP), y=(PerCap))) +
  #geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70") +
  geom_point(aes(), size=1) +
  #scale_color_gradient(low ="#F8766D", high = "#619CFF") +
  #scale_color_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  #geom_line(aes(y=logFit3), color="blue", size=1) +
  theme_bw() +
  scale_x_reverse(breaks=c(3500, 2500, 1500, 500), limits=c(3500,400))+
  scale_y_continuous(limits=c(-.3,0.2))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Years Cal BP", y="SPD per capita growth", title = "A. Mississippi R. Valley Per Capita Growth")+
  geom_line(size=1.25)+
  geom_vline(xintercept = 1900)+
  geom_vline(xintercept = 1250)+
  geom_vline(xintercept = 900)+
  geom_vline(xintercept = 750)+
  geom_hline(yintercept = 0)
upmiss

misphase <- ggplot(Miss2,aes(x=(PrDens), y=(PerCap))) +
  #geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70") +
  geom_point(aes(color=calBP), size=1) +
  #scale_color_gradient(low ="#F8766D", high = "#619CFF") +
  #scale_color_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  #geom_line(aes(y=logFit3), color="blue", size=1) +
  theme_bw() +
  #scale_x_reverse(breaks=c(3500, 2500, 1500, 500), limits=c(3500,400))+
  scale_y_continuous(limits=c(-.3,0.2))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Years Cal BP", y="SPD per capita growth", title = "A. Miss. R. Valley Per Capita Growth")
# geom_line(size=1.25)
misphase


#### KDE for Miss. Valley US=====================================================================
Miss.randates = sampleDates(cptcal2, bins=bin, nsim=200,verbose=FALSE)
D.ckde = ckde(Miss.randates,timeRange=c(3500,200),bw=50, normalised = FALSE)
plot(D.ckde,type='multiline')

##Write matrix of KDEs as a data frame
Check<-as.data.frame(D.ckde$res.matrix)

##Cbind spd and KDEs and write
dd<-cbind(calBP,PrDens, Check)
write.table(dd, file = "data/MissKDE50.csv", sep = ",", col.names=NA)

###At this point, I open the spreadsheet, and I add a two new columns: MKDE and PerCap.
#MKDE I calculate as the average of the KDE values at each time step. So, in this case, I add the KDE at cal BP 4000 from V1-V200
#and then divide by 200. This gives the mean KDE value at 4000 cal BP. Once I have done this for all time points, I calculate
#per capita growth as LN(MKDE at t+1) -LN(MKDE at time t). Once this is done, I load the csv file back into R.

#Load KDE data frame with calculated mean KDE and per capita growth
dd2m<- read.csv("data/MissKDE50.csv")
#I then trimm time points with missing values at the end due to KDE bandwidth.
dd2m<-na.omit(dd2m)

###transform dd2 into long table format to plot all of the simulated KDEs
dd3m<-dd2m %>% pivot_longer(cols=c(V1:V200),names_to='Runs', values_to='KDE')

##OK, now we fit a logistic model to the mean KDE. This can be a bit tricky due to the start parameters.
#These need to be adjusted based on the MKDE's values. In this case, the start values work well.
#fit logistic to KDE
logFitKDm <- nls(MKDE~SSfpl(calBP, A, B, xmid, scale), data=dd2m,control=nls.control(maxiter=700),start=list(A=.00005, B=0.00015, xmid=1711, scale=-100))
summary(logFitKDm)
predictm<-predict(logFitKDm)

##Combine predicted logistic values and KDE table for plotting
dd4m<-cbind(dd2m, predictm)

###Simple plot
plot(dd2m$MKDE~dd2m$calBP)
lines(dd2m$calBP, dd2m$MKDE,col="green",lty=2,lwd=3)
lines(dd2m$calBP,predict(logFitKDm),col="red",lty=2,lwd=3)

###Plot mean KDE and all of the simulated KDE models using ggplot
p1<-ggplot(data = dd3m) +
  geom_line(aes(x = calBP, y = (KDE*100), group = Runs), color = "purple", alpha = 0.5) +
  geom_point(data = dd4m, aes(x = calBP, y = (MKDE*100)), color="black", size=2) +
  geom_line(data = dd2m, aes(x = calBP, y = (predictm*100)), linetype = "dashed", color="green", size=1.5) +
  scale_x_reverse(breaks=c(3500, 2500, 1500, 500), limits=c(3500,300))+
  theme_bw() +
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Years Cal BP", y="KDE of radiocarbon ages", title = "A. Mississippi R. Valley KDE and Logistic Model")+
  # annotate("text", x =3500, y = .047, label = "Phase 1", size = 7, angle=90)+
  #annotate("text", x =2000, y = .047, label = "Phase 2", size = 7, angle=90)+
  #annotate("text", x =1020, y = .047, label = "Phase 3", size = 7, angle=90)+
  #annotate("text", x =210, y = .047, label = "Phase 4", size = 7, angle=90)+
  # facet_wrap(factor(Region)~.)
  geom_vline(xintercept = 1900)+
  geom_vline(xintercept = 1250)+
  geom_vline(xintercept = 900)+
  geom_vline(xintercept = 750)
p1

dd2c<- read.csv("data/MissKDE50.csv")
#I then trimm time points with missing values at the end due to KDE bandwidth.
dd2c<-na.omit(dd2c)
dd2c<-subset(dd2c, select=-c(V1:V200))
###Sum SPD by 30 year intervals
library(zoo)
# sum and save new csvs.
out50 <- rollapply(dd2c,30,(sum),by=30,by.column=TRUE,align='right')
write.table(out50, file = "sums/MissKDE30.csv", sep = ",", col.names=NA)

pc30<- read.csv("sums/MissKDE30.csv")

Misspc <- ggplot(pc30,aes(x=(calBP), y=(PerCap))) +
  #geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70") +
  geom_point(aes(), size=1) +
  #scale_color_gradient(low ="#F8766D", high = "#619CFF") +
  #scale_color_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  #geom_line(aes(y=logFit3), color="blue", size=1) +
  theme_bw() +
  scale_x_reverse(breaks=c(3500, 2500, 1500, 500), limits=c(3500,300))+
  scale_y_continuous(limits=c(-.2,0.25))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Years Cal BP", y="KDE per capita growth", title = "C. Mississippi R. Valley KDE Per Capita Growth")+
  geom_line(size=1.25)+
  geom_vline(xintercept = 1900)+
  geom_vline(xintercept = 1250)+
  geom_vline(xintercept = 900)+
  geom_vline(xintercept = 750)+
  geom_hline(yintercept = 0)
Misspc


######Plots that compare Central Texas and the Mississippi River Valley
Fig5Rev<-plot_grid(p1, p1ct, Misspc, ctexpc,pn1,pn2, ncol=2, align="hv", axis = "rl")
Fig5Rev

pdf("CtexMiss3.pdf", width=15, height=17.55)
Fig5Rev
dev.off()
#==================================================================================================================================

####ArchaeoGlobe for the lower 48 US states. Run the SPD by region at 100 year smoothing and sum into 30 year intervals
#Also run a KDE with 50 year bandwidth. Sum the mean KDE into 30 year intervals. Need a picture
#of the overall curves and of the per capita growth rates.

##First, map the data==================================================================================================
library(maps)
library(sf)
library(tidyverse)
library(usmap)
#library(mapdata)

###Load raw radiocarbon data for mapping:
SPD<-read.csv(file="RAWNA1Map.csv", header=T)

##Map radiocarbon data by US ArchaeoGlobe Regions
counties<-map_data("state")
#canada <- map_data("worldHires", "Canada")

ArchGlobeMap<-ggplot() +
  geom_polygon(data = counties, mapping = aes(x = long, y = lat, group = group),
               fill = "grey", color = "white") +
  coord_map(projection = "albers", lat0 = 39, lat1 = 45)+
  #geom_polygon(data = canada, aes(x=long, y = lat, group = group),
          #    fill = "white", color="black") +
  theme_bw()+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=20, face = "bold"))+
  labs(x = "Longitude", y="Latitude", title = "US ArchaeoGlobe Regions and Radiocarbon")+
  geom_point(data=SPD, aes(MapLong, Latitude, color=factor(Archaeo_RG)),
             inherit.aes = FALSE, alpha = 0.5, size = 2)
ArchGlobeMap

###Make map of isotope sites using data sets from the case study analysis above.

Cases<-ggplot() +
  geom_polygon(data = counties, mapping = aes(x = long, y = lat, group = group),
               fill = "grey", color = "white") +
  coord_map(projection = "albers", lat0 = 39, lat1 = 45)+
  theme_bw()+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=20, face = "bold"))+
  labs(x = "Longitude", y="Latitude", title = "Case Study Sites With Bone Isotope Data")+
  geom_point(data=d5, aes(Estimated.Longitude, Estimated.Latitude),
             inherit.aes = FALSE, alpha = 0.5, size = 2)+
  geom_point(data=dCIII, aes(Longitude, Latitude),
             inherit.aes = FALSE, alpha = 0.5, size = 2)
Cases

##Load radiocarbon data for Archaeoglobe regional analysis====================================================
SPD<-read.csv(file="RawNA1.csv", header=T)

#SPlit data for SW region
SPD.1<- subset(SPD, Archaeo_ID=="15")

##calibrate
cptcal <- calibrate(x = SPD.1$Age,  errors = SPD.1$Error, calCurves = "intcal20",  normalised = FALSE)
# <- binPrep(sites = SPD.1$SiteID, ages = SPD.1$Age, h = 100)

spd.sw <- spd(cptcal, bins=NA, runm=100, timeRange=c(3500,200))
plot(spd.sw, runm=100, xlim=c(3500,200), type="simple")

PrDens<-spd.sw$grid$PrDens
calBP<-spd.sw$grid$calBP

SW.randates = sampleDates(cptcal,bins=NA,nsim=200,verbose=FALSE)
D.ckde = ckde(SW.randates,timeRange=c(3500,200),bw=50, normalised = FALSE)
plot(D.ckde,type='multiline')

##Write matrix of KDEs as a data frame
Check<-as.data.frame(D.ckde$res.matrix)

##Cbind spd and KDEs and write
dd<-cbind(calBP,PrDens, Check)
write.table(dd, file = "data/AGSWKDE50.csv", sep = ",", col.names=NA)

###At this point, I open the spreadsheet, and I add a two new columns: MKDE and PerCap.
#MKDE I calculate as the average of the KDE values at each time step. So, in this case, I add the KDE at cal BP 3500 from V1-V200
#and then divide by 200. This gives the mean KDE value at 3500 cal BP. Once I have done this for all time points, I calculate
#per capita growth as LN(MKDE at t+1) -LN(MKDE at time t). Once this is done, I load the csv file back into R.

#Load KDE data frame with calculated mean KDE and per capita growth
dd2m<- read.csv("data/AGSWKDE50.csv")
#I then trimm time points with missing values at the end due to KDE bandwidth.
dd2m<-na.omit(dd2m)

###transform dd2 into long table format to plot all of the simulated KDEs
dd3m<-dd2m %>% pivot_longer(cols=c(V1:V200),names_to='Runs', values_to='KDE')

##OK, now we fit a logistic model to the mean KDE. This can be a bit tricky due to the start parameters.
#These need to be adjusted based on the MKDE's values. In this case, the start values work well.
#fit logistic to KDE
logFitKDm <- nls(MKDE~SSfpl(calBP, A, B, xmid, scale), data=dd2m,control=nls.control(maxiter=700),start=list(A=.00005, B=0.00015, xmid=1711, scale=-100))
summary(logFitKDm)
predictm<-predict(logFitKDm)

##Combine predicted logistic values and KDE table for plotting
dd4m<-cbind(dd2m, predictm)
#write.table(dd4, file = "SWKDE50.csv", sep = ",", col.names=NA)

###Simple plot
plot(dd2m$MKDE~dd2m$calBP)
lines(dd2m$calBP, dd2m$MKDE,col="green",lty=2,lwd=3)
lines(dd2m$calBP,predict(logFitKDm),col="red",lty=2,lwd=3)

###Plot mean KDE and all of the simulated KDE models using ggplot
p1<-ggplot(data = dd3m) +
  geom_line(aes(x = calBP, y = (KDE*100), group = Runs), color = "purple", alpha = 0.5) +
  geom_point(data = dd4m, aes(x = calBP, y = (MKDE*100)), color="black", size=2) +
  geom_line(data = dd2m, aes(x = calBP, y = (predictm*100)), linetype = "dashed", color="green", size=1.5) +
  scale_x_reverse(breaks=c(3500, 2500, 1500, 500), limits=c(3500,400))+
  theme_bw() +
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Years Cal BP", y="KDE of radiocarbon ages", title = "A. SW ArchaeoGlobe KDE and Logistic Model")
  # annotate("text", x =3500, y = .047, label = "Phase 1", size = 7, angle=90)+
  #annotate("text", x =2000, y = .047, label = "Phase 2", size = 7, angle=90)+
  #annotate("text", x =1020, y = .047, label = "Phase 3", size = 7, angle=90)+
  #annotate("text", x =210, y = .047, label = "Phase 4", size = 7, angle=90)+
  # facet_wrap(factor(Region)~.)
p1

pdf("SI/SWPop.pdf", width=8, height=6.55)
p1
dev.off()

dd2c<- read.csv("data/AGSWKDE50.csv")
#I then trimm time points with missing values at the end due to KDE bandwidth.
dd2c<-na.omit(dd2c)
dd2c<-subset(dd2c, select=-c(V1:V200))
###Sum SPD by 30 year intervals
library(zoo)
# sum and save new csvs.
out50 <- rollapply(dd2c,30,(sum),by=30,by.column=TRUE,align='right')
write.table(out50, file = "sums/AGSWKDE30.csv", sep = ",", col.names=NA)

pc30sw<- read.csv("sums/AGSWKDE30.csv")

swpc <- ggplot(pc30sw,aes(x=(calBP), y=(PerCap))) +
  #geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70") +
  geom_point(aes(), size=1) +
  #scale_color_gradient(low ="#F8766D", high = "#619CFF") +
  #scale_color_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  #geom_line(aes(y=logFit3), color="blue", size=1) +
  theme_bw() +
  scale_x_reverse(breaks=c(3500, 2500, 1500, 500), limits=c(3500,300))+
  scale_y_continuous(limits=c(-.2,0.25))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Years Cal BP", y="KDE per capita growth", title = "A. SW ArchaeoGlobe KDE Per Capita Growth")+
  geom_line(size=1.25)+
  geom_hline(yintercept = 0)
swpc

swpcs <- ggplot(pc30sw,aes(x=(calBP), y=(spdPerCap))) +
  geom_point(aes(), size=1) +
  #scale_color_gradient(low ="#F8766D", high = "#619CFF") +
  #scale_color_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  #geom_line(aes(y=logFit3), color="blue", size=1) +
  theme_bw() +
  scale_x_reverse(breaks=c(3500, 2500, 1500, 500), limits=c(3500,300))+
  scale_y_continuous(limits=c(-.2,0.25))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Years Cal BP", y="KDE per capita growth", title = "A. SW ArchaeoGlobe KDE Per Capita Growth")+
  geom_line(size=1.25)+
  geom_line(aes(x=calBP, y=PerCap),size=1.25, color="red")+
  geom_hline(yintercept = 0)
swpcs

pdf("SI/SWPerCap.pdf", width=8, height=6.55)
swpcs
dev.off()

####Western US==================================================

###Run logistic for the western US
SPD.3<- subset(SPD, Archaeo_ID=="14")
##calibrate
cptcal <- calibrate(x = SPD.3$Age,  errors = SPD.3$Error, calCurves = "intcal20",  normalised = FALSE)
# <- binPrep(sites = box2$Trinomial, ages = box2$DatesUse, h = 100)

####Run analysis for component 3 logistic 3400 to 150
spd.we<- spd(cptcal, bins=NA, runm=100, timeRange=c(3500,200))
plot(spd.we, runm=100, xlim=c(3500,200), type="simple")


PrDens<-spd.we$grid$PrDens
calBP<-spd.we$grid$calBP

WE.randates = sampleDates(cptcal,bins=NA,nsim=200,verbose=FALSE)
D.ckde = ckde(WE.randates,timeRange=c(3500,200),bw=50, normalised = FALSE)
plot(D.ckde,type='multiline')

##Write matrix of KDEs as a data frame
Check<-as.data.frame(D.ckde$res.matrix)

##Cbind spd and KDEs and write
dd<-cbind(calBP,PrDens, Check)
write.table(dd, file = "data/AGWEKDE50.csv", sep = ",", col.names=NA)

###At this point, I open the spreadsheet, and I add a two new columns: MKDE and PerCap.
#MKDE I calculate as the average of the KDE values at each time step. So, in this case, I add the KDE at cal BP 3500 from V1-V200
#and then divide by 200. This gives the mean KDE value at 3500 cal BP. Once I have done this for all time points, I calculate
#per capita growth as LN(MKDE at t+1) -LN(MKDE at time t). Once this is done, I load the csv file back into R.

#Load KDE data frame with calculated mean KDE and per capita growth
dd2m<- read.csv("data/AGWEKDE50.csv")
#I then trimm time points with missing values at the end due to KDE bandwidth.
dd2m<-na.omit(dd2m)
dd2m<-subset(dd2m,calBP>500)
###transform dd2 into long table format to plot all of the simulated KDEs
dd3m<-dd2m %>% pivot_longer(cols=c(V1:V200),names_to='Runs', values_to='KDE')

##OK, now we fit a logistic model to the mean KDE. This can be a bit tricky due to the start parameters.
#These need to be adjusted based on the MKDE's values. In this case, the start values work well.
#fit logistic to KDE
logFitKDm <- nls(MKDE~SSfpl(calBP, A, B, xmid, scale), data=dd2m,control=nls.control(maxiter=1000),start=list(A=.00005, B=0.00015, xmid=1700, scale=-300))
summary(logFitKDm)
predictm<-predict(logFitKDm)

##Combine predicted logistic values and KDE table for plotting
dd4m<-cbind(dd2m, predictm)
#write.table(dd4, file = "SWKDE50.csv", sep = ",", col.names=NA)

###Simple plot
plot(dd2m$MKDE~dd2m$calBP)
lines(dd2m$calBP, dd2m$MKDE,col="green",lty=2,lwd=3)
lines(dd2m$calBP,predict(logFitKDm),col="red",lty=2,lwd=3)

###Plot mean KDE and all of the simulated KDE models using ggplot
p2<-ggplot(data = dd3m) +
  geom_line(aes(x = calBP, y = (KDE*100), group = Runs), color = "purple", alpha = 0.5) +
  geom_point(data = dd4m, aes(x = calBP, y = (MKDE*100)), color="black", size=2) +
  geom_line(data = dd2m, aes(x = calBP, y = (predictm*100)), linetype = "dashed", color="green", size=1.5) +
  scale_x_reverse(breaks=c(3500, 2500, 1500, 500), limits=c(3500,400))+
  theme_bw() +
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Years Cal BP", y="KDE of radiocarbon ages", title = "B. WE ArchaeoGlobe KDE and Logistic Model")
# annotate("text", x =3500, y = .047, label = "Phase 1", size = 7, angle=90)+
#annotate("text", x =2000, y = .047, label = "Phase 2", size = 7, angle=90)+
#annotate("text", x =1020, y = .047, label = "Phase 3", size = 7, angle=90)+
#annotate("text", x =210, y = .047, label = "Phase 4", size = 7, angle=90)+
# facet_wrap(factor(Region)~.)
p2

pdf("SI/WEPop.pdf", width=8, height=6.55)
p2
dev.off()

dd2c<- read.csv("data/AGWEKDE50.csv")
#I then trimm time points with missing values at the end due to KDE bandwidth.
dd2c<-na.omit(dd2c)
dd2c<-subset(dd2c, select=-c(V1:V200))
###Sum SPD by 30 year intervals
library(zoo)
# sum and save new csvs.
out50 <- rollapply(dd2c,30,(sum),by=30,by.column=TRUE,align='right')
write.table(out50, file = "sums/AGWEKDE30.csv", sep = ",", col.names=NA)

pc30we<- read.csv("sums/AGWEKDE30.csv")
pc30we2<-subset(pc30we,calBP>500)

wepc <- ggplot(pc30we,aes(x=(calBP), y=(PerCap))) +
  #geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70") +
  geom_point(aes(), size=1) +
  #scale_color_gradient(low ="#F8766D", high = "#619CFF") +
  #scale_color_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  #geom_line(aes(y=logFit3), color="blue", size=1) +
  theme_bw() +
  scale_x_reverse(breaks=c(3500, 2500, 1500, 500), limits=c(3500,300))+
  scale_y_continuous(limits=c(-.2,0.25))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Years Cal BP", y="KDE per capita growth", title = "B. WE ArchaeoGlobe KDE Per Capita Growth")+
  geom_line(size=1.25)+
  geom_hline(yintercept = 0)
wepc

wepcs <- ggplot(pc30we,aes(x=(calBP), y=(spdPerCap))) +
  geom_point(aes(), size=1) +
  #scale_color_gradient(low ="#F8766D", high = "#619CFF") +
  #scale_color_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  #geom_line(aes(y=logFit3), color="blue", size=1) +
  theme_bw() +
  scale_x_reverse(breaks=c(3500, 2500, 1500, 500), limits=c(3500,400))+
  scale_y_continuous(limits=c(-.2,0.25))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Years Cal BP", y="KDE per capita growth", title = "B. WE ArchaeoGlobe KDE Per Capita Growth")+
  geom_line(size=1.25)+
  geom_line(aes(x=calBP, y=PerCap),size=1.25, color="red")+
  geom_hline(yintercept = 0)
wepcs

pdf("SI/WEPerCap.pdf", width=8, height=6.55)
wepcs
dev.off()


########Central Plains
SPD.3<- subset(SPD, Archaeo_ID=="17")
#write.table(SPD.3, file = "Northeast.csv", sep = ",")
##calibrate
cptcal <- calibrate(x = SPD.3$Age,  errors = SPD.3$Error, calCurves = "intcal20",  normalised = FALSE)
# <- binPrep(sites = box2$Trinomial, ages = box2$DatesUse, h = 100)

####Run analysis for component 3 logistic 3400 to 150
spd.cp<- spd(cptcal, bins=NA, runm=100, timeRange=c(3500,200))
plot(spd.cp, runm=100, xlim=c(3500,200), type="simple")

PrDens<-spd.cp$grid$PrDens
calBP<-spd.cp$grid$calBP

CP.randates = sampleDates(cptcal,bins=NA,nsim=200,verbose=FALSE)
D.ckde = ckde(CP.randates,timeRange=c(3500,200),bw=50, normalised = FALSE)
plot(D.ckde,type='multiline')

##Write matrix of KDEs as a data frame
Check<-as.data.frame(D.ckde$res.matrix)

##Cbind spd and KDEs and write
dd<-cbind(calBP,PrDens, Check)
write.table(dd, file = "data/AGCPKDE50.csv", sep = ",", col.names=NA)

###At this point, I open the spreadsheet, and I add a two new columns: MKDE and PerCap.
#MKDE I calculate as the average of the KDE values at each time step. So, in this case, I add the KDE at cal BP 3500 from V1-V200
#and then divide by 200. This gives the mean KDE value at 3500 cal BP. Once I have done this for all time points, I calculate
#per capita growth as LN(MKDE at t+1) -LN(MKDE at time t). Once this is done, I load the csv file back into R.

#Load KDE data frame with calculated mean KDE and per capita growth
dd2m<- read.csv("data/AGCPKDE50.csv")
#I then trimm time points with missing values at the end due to KDE bandwidth.
dd2m<-na.omit(dd2m)

###transform dd2 into long table format to plot all of the simulated KDEs
dd3m<-dd2m %>% pivot_longer(cols=c(V1:V200),names_to='Runs', values_to='KDE')

##OK, now we fit a logistic model to the mean KDE. This can be a bit tricky due to the start parameters.
#These need to be adjusted based on the MKDE's values. In this case, the start values work well.
#fit logistic to KDE
logFitKDm <- nls(MKDE~SSfpl(calBP, A, B, xmid, scale), data=dd2m,control=nls.control(maxiter=700),start=list(A=.00005, B=0.00015, xmid=1200, scale=-100))
summary(logFitKDm)
predictm<-predict(logFitKDm)

##Combine predicted logistic values and KDE table for plotting
dd4m<-cbind(dd2m, predictm)
#write.table(dd4, file = "SWKDE50.csv", sep = ",", col.names=NA)

###Simple plot
plot(dd2m$MKDE~dd2m$calBP)
lines(dd2m$calBP, dd2m$MKDE,col="green",lty=2,lwd=3)
lines(dd2m$calBP,predict(logFitKDm),col="red",lty=2,lwd=3)

###Plot mean KDE and all of the simulated KDE models using ggplot
p3<-ggplot(data = dd3m) +
  geom_line(aes(x = calBP, y = (KDE*100), group = Runs), color = "purple", alpha = 0.5) +
  geom_point(data = dd4m, aes(x = calBP, y = (MKDE*100)), color="black", size=2) +
  geom_line(data = dd2m, aes(x = calBP, y = (predictm*100)), linetype = "dashed", color="green", size=1.5) +
  scale_x_reverse(breaks=c(3500, 2500, 1500, 500), limits=c(3500,400))+
  theme_bw() +
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Years Cal BP", y="KDE of radiocarbon ages", title = "C. CP ArchaeoGlobe KDE and Logistic Model")
# annotate("text", x =3500, y = .047, label = "Phase 1", size = 7, angle=90)+
#annotate("text", x =2000, y = .047, label = "Phase 2", size = 7, angle=90)+
#annotate("text", x =1020, y = .047, label = "Phase 3", size = 7, angle=90)+
#annotate("text", x =210, y = .047, label = "Phase 4", size = 7, angle=90)+
# facet_wrap(factor(Region)~.)
p3

pdf("SI/CPPop.pdf", width=8, height=6.55)
p3
dev.off()

dd2c<- read.csv("data/AGCPKDE50.csv")
#I then trimm time points with missing values at the end due to KDE bandwidth.
dd2c<-na.omit(dd2c)
dd2c<-subset(dd2c, select=-c(V1:V200))
###Sum SPD by 30 year intervals
library(zoo)
# sum and save new csvs.
out50 <- rollapply(dd2c,30,(sum),by=30,by.column=TRUE,align='right')
write.table(out50, file = "sums/AGCPKDE30.csv", sep = ",", col.names=NA)

pc30cp<- read.csv("sums/AGCPKDE30.csv")

cppc <- ggplot(pc30cp,aes(x=(calBP), y=(PerCap))) +
  #geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70") +
  geom_point(aes(), size=1) +
  #scale_color_gradient(low ="#F8766D", high = "#619CFF") +
  #scale_color_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  #geom_line(aes(y=logFit3), color="blue", size=1) +
  theme_bw() +
  scale_x_reverse(breaks=c(3500, 2500, 1500, 500), limits=c(3500,300))+
  scale_y_continuous(limits=c(-.2,0.25))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Years Cal BP", y="KDE per capita growth", title = "C. CP ArchaeoGlobe KDE Per Capita Growth")+
  geom_line(size=1.25)+
  geom_hline(yintercept = 0)
cppc

cppcs<- ggplot(pc30cp,aes(x=(calBP), y=(spdPerCap))) +
  geom_point(aes(), size=1) +
  #scale_color_gradient(low ="#F8766D", high = "#619CFF") +
  #scale_color_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  #geom_line(aes(y=logFit3), color="blue", size=1) +
  theme_bw() +
  scale_x_reverse(breaks=c(3500, 2500, 1500, 500), limits=c(3500,400))+
  scale_y_continuous(limits=c(-.2,0.25))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Years Cal BP", y="KDE per capita growth", title = "C. CP ArchaeoGlobe KDE Per Capita Growth")+
  geom_line(size=1.25)+
  geom_line(aes(x=calBP, y=PerCap),size=1.25, color="red")+
  geom_hline(yintercept = 0)
cppcs

pdf("SI/CPPerCap.pdf", width=8, height=6.55)
cppcs
dev.off()


#######Southeastern==========================================================

###Run logistic for the Southeastern US
SPD.3<- subset(SPD, Archaeo_ID=="21")
##calibrate
cptcal <- calibrate(x = SPD.3$Age,  errors = SPD.3$Error, calCurves = "intcal20",  normalised = FALSE)
# <- binPrep(sites = box2$Trinomial, ages = box2$DatesUse, h = 100)


####Run analysis for component 3 logistic 3400 to 150
spd.se <- spd(cptcal, bins=NA, runm=100, timeRange=c(3500,200))
plot(spd.se, runm=100, xlim=c(3500,200), type="simple")

PrDens<-spd.se$grid$PrDens
calBP<-spd.se$grid$calBP

SE.randates = sampleDates(cptcal,bins=NA,nsim=200,verbose=FALSE)
D.ckde = ckde(SE.randates,timeRange=c(3500,200),bw=50, normalised = FALSE)
plot(D.ckde,type='multiline')

##Write matrix of KDEs as a data frame
Check<-as.data.frame(D.ckde$res.matrix)

##Cbind spd and KDEs and write
dd<-cbind(calBP,PrDens, Check)
write.table(dd, file = "data/AGSEKDE50.csv", sep = ",", col.names=NA)

###At this point, I open the spreadsheet, and I add a two new columns: MKDE and PerCap.
#MKDE I calculate as the average of the KDE values at each time step. So, in this case, I add the KDE at cal BP 3500 from V1-V200
#and then divide by 200. This gives the mean KDE value at 3500 cal BP. Once I have done this for all time points, I calculate
#per capita growth as LN(MKDE at t+1) -LN(MKDE at time t). Once this is done, I load the csv file back into R.

#Load KDE data frame with calculated mean KDE and per capita growth
dd2m<- read.csv("data/AGSEKDE50.csv")
#I then trimm time points with missing values at the end due to KDE bandwidth.
dd2m<-na.omit(dd2m)

###transform dd2 into long table format to plot all of the simulated KDEs
dd3m<-dd2m %>% pivot_longer(cols=c(V1:V200),names_to='Runs', values_to='KDE')

##OK, now we fit a logistic model to the mean KDE. This can be a bit tricky due to the start parameters.
#These need to be adjusted based on the MKDE's values. In this case, the start values work well.
#fit logistic to KDE
logFitKDm <- nls(MKDE~SSfpl(calBP, A, B, xmid, scale), data=dd2m,control=nls.control(maxiter=700),start=list(A=.00005, B=0.00015, xmid=1200, scale=-100))
summary(logFitKDm)
predictm<-predict(logFitKDm)

##Combine predicted logistic values and KDE table for plotting
dd4m<-cbind(dd2m, predictm)
#write.table(dd4, file = "SWKDE50.csv", sep = ",", col.names=NA)

###Simple plot
plot(dd2m$MKDE~dd2m$calBP)
lines(dd2m$calBP, dd2m$MKDE,col="green",lty=2,lwd=3)
lines(dd2m$calBP,predict(logFitKDm),col="red",lty=2,lwd=3)

###Plot mean KDE and all of the simulated KDE models using ggplot
p4<-ggplot(data = dd3m) +
  geom_line(aes(x = calBP, y = (KDE*100), group = Runs), color = "purple", alpha = 0.5) +
  geom_point(data = dd4m, aes(x = calBP, y = (MKDE*100)), color="black", size=2) +
  geom_line(data = dd2m, aes(x = calBP, y = (predictm*100)), linetype = "dashed", color="green", size=1.5) +
  scale_x_reverse(breaks=c(3500, 2500, 1500, 500), limits=c(3500,400))+
  theme_bw() +
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Years Cal BP", y="KDE of radiocarbon ages", title = "D. SE ArchaeoGlobe KDE and Logistic Model")
# annotate("text", x =3500, y = .047, label = "Phase 1", size = 7, angle=90)+
#annotate("text", x =2000, y = .047, label = "Phase 2", size = 7, angle=90)+
#annotate("text", x =1020, y = .047, label = "Phase 3", size = 7, angle=90)+
#annotate("text", x =210, y = .047, label = "Phase 4", size = 7, angle=90)+
# facet_wrap(factor(Region)~.)
p4

pdf("SI/SEPop.pdf", width=8, height=6.55)
p4
dev.off()

dd2c<- read.csv("data/AGSEKDE50.csv")
#I then trimm time points with missing values at the end due to KDE bandwidth.
dd2c<-na.omit(dd2c)
dd2c<-subset(dd2c, select=-c(V1:V200))
###Sum SPD by 30 year intervals
library(zoo)
# sum and save new csvs.
out50 <- rollapply(dd2c,30,(sum),by=30,by.column=TRUE,align='right')
write.table(out50, file = "sums/AGSEKDE30.csv", sep = ",", col.names=NA)

pc30se<- read.csv("sums/AGSEKDE30.csv")

sepc <- ggplot(pc30se,aes(x=(calBP), y=(PerCap))) +
  #geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70") +
  geom_point(aes(), size=1) +
  #scale_color_gradient(low ="#F8766D", high = "#619CFF") +
  #scale_color_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  #geom_line(aes(y=logFit3), color="blue", size=1) +
  theme_bw() +
  scale_x_reverse(breaks=c(3500, 2500, 1500, 500), limits=c(3500,300))+
  scale_y_continuous(limits=c(-.2,0.25))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Years Cal BP", y="KDE per capita growth", title = "D. SE ArchaeoGlobe KDE Per Capita Growth")+
  geom_line(size=1.25)+
  geom_hline(yintercept = 0)
sepc

sepcs <- ggplot(pc30se,aes(x=(calBP), y=(spdPerCap))) +
  geom_point(aes(), size=1) +
  #scale_color_gradient(low ="#F8766D", high = "#619CFF") +
  #scale_color_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  #geom_line(aes(y=logFit3), color="blue", size=1) +
  theme_bw() +
  scale_x_reverse(breaks=c(3500, 2500, 1500, 500), limits=c(3500,300))+
  scale_y_continuous(limits=c(-.2,0.25))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Years Cal BP", y="KDE per capita growth", title = "D. SE ArchaeoGlobe KDE Per Capita Growth")+
  geom_line(size=1.25)+
  geom_line(aes(x=calBP, y=PerCap),size=1.25, color="red")+
  geom_hline(yintercept = 0)
sepcs

pdf("SI/SEPerCap.pdf", width=8, height=6.55)
sepcs
dev.off()

#######Midwest==========================================================

SPD.2<- subset(SPD, Archaeo_ID=="18")
##calibrate
cptcal <- calibrate(x = SPD.2$Age,  errors = SPD.2$Error, calCurves = "intcal20",  normalised = FALSE)
# <- binPrep(sites = box2$Trinomial, ages = box2$DatesUse, h = 100)


####Run analysis for component 3 logistic 3400 to 150
spd.mw <- spd(cptcal, bins=NA, runm=100, timeRange=c(3500,200))
plot(spd.mw, runm=100, xlim=c(3500,200), type="simple")

PrDens<-spd.mw$grid$PrDens
calBP<-spd.mw$grid$calBP

MW.randates = sampleDates(cptcal,bins=NA,nsim=200,verbose=FALSE)
D.ckde = ckde(MW.randates,timeRange=c(3500,200),bw=50, normalised = FALSE)
plot(D.ckde,type='multiline')

##Write matrix of KDEs as a data frame
Check<-as.data.frame(D.ckde$res.matrix)

##Cbind spd and KDEs and write
dd<-cbind(calBP,PrDens, Check)
write.table(dd, file = "data/AGMWKDE50.csv", sep = ",", col.names=NA)

###At this point, I open the spreadsheet, and I add a two new columns: MKDE and PerCap.
#MKDE I calculate as the average of the KDE values at each time step. So, in this case, I add the KDE at cal BP 3500 from V1-V200
#and then divide by 200. This gives the mean KDE value at 3500 cal BP. Once I have done this for all time points, I calculate
#per capita growth as LN(MKDE at t+1) -LN(MKDE at time t). Once this is done, I load the csv file back into R.

#Load KDE data frame with calculated mean KDE and per capita growth
dd2m<- read.csv("data/AGMWKDE50.csv")
#I then trimm time points with missing values at the end due to KDE bandwidth.
dd2m<-na.omit(dd2m)

###transform dd2 into long table format to plot all of the simulated KDEs
dd3m<-dd2m %>% pivot_longer(cols=c(V1:V200),names_to='Runs', values_to='KDE')

##OK, now we fit a logistic model to the mean KDE. This can be a bit tricky due to the start parameters.
#These need to be adjusted based on the MKDE's values. In this case, the start values work well.
#fit logistic to KDE
logFitKDm <- nls(MKDE~SSfpl(calBP, A, B, xmid, scale), data=dd2m,control=nls.control(maxiter=700),start=list(A=.00005, B=0.00015, xmid=1200, scale=-100))
summary(logFitKDm)
predictm<-predict(logFitKDm)

##Combine predicted logistic values and KDE table for plotting
dd4m<-cbind(dd2m, predictm)
#write.table(dd4, file = "SWKDE50.csv", sep = ",", col.names=NA)

###Simple plot
plot(dd2m$MKDE~dd2m$calBP)
lines(dd2m$calBP, dd2m$MKDE,col="green",lty=2,lwd=3)
lines(dd2m$calBP,predict(logFitKDm),col="red",lty=2,lwd=3)

###Plot mean KDE and all of the simulated KDE models using ggplot
p5<-ggplot(data = dd3m) +
  geom_line(aes(x = calBP, y = (KDE*100), group = Runs), color = "purple", alpha = 0.5) +
  geom_point(data = dd4m, aes(x = calBP, y = (MKDE*100)), color="black", size=2) +
  geom_line(data = dd2m, aes(x = calBP, y = (predictm*100)), linetype = "dashed", color="green", size=1.5) +
  scale_x_reverse(breaks=c(3500, 2500, 1500, 500), limits=c(3500,400))+
  theme_bw() +
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Years Cal BP", y="KDE of radiocarbon ages", title = "E. MW ArchaeoGlobe KDE and Logistic Model")
# annotate("text", x =3500, y = .047, label = "Phase 1", size = 7, angle=90)+
#annotate("text", x =2000, y = .047, label = "Phase 2", size = 7, angle=90)+
#annotate("text", x =1020, y = .047, label = "Phase 3", size = 7, angle=90)+
#annotate("text", x =210, y = .047, label = "Phase 4", size = 7, angle=90)+
# facet_wrap(factor(Region)~.)
p5

pdf("SI/MWPop.pdf", width=8, height=6.55)
p5
dev.off()

dd2c<- read.csv("data/AGMWKDE50.csv")
#I then trimm time points with missing values at the end due to KDE bandwidth.
dd2c<-na.omit(dd2c)
dd2c<-subset(dd2c, select=-c(V1:V200))
###Sum SPD by 30 year intervals
library(zoo)
# sum and save new csvs.
out50 <- rollapply(dd2c,30,(sum),by=30,by.column=TRUE,align='right')
write.table(out50, file = "sums/AGMWKDE30.csv", sep = ",", col.names=NA)

pc30mw<- read.csv("sums/AGMWKDE30.csv")

mwpc <- ggplot(pc30mw,aes(x=(calBP), y=(PerCap))) +
  #geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70") +
  geom_point(aes(), size=1) +
  #scale_color_gradient(low ="#F8766D", high = "#619CFF") +
  #scale_color_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  #geom_line(aes(y=logFit3), color="blue", size=1) +
  theme_bw() +
  scale_x_reverse(breaks=c(3500, 2500, 1500, 500), limits=c(3500,300))+
  scale_y_continuous(limits=c(-.2,0.25))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Years Cal BP", y="KDE per capita growth", title = "E. MW ArchaeoGlobe KDE Per Capita Growth")+
  geom_line(size=1.25)+
  geom_hline(yintercept = 0)
mwpc

mwpcs <- ggplot(pc30mw,aes(x=(calBP), y=(spdPerCap))) +
  geom_point(aes(), size=1) +
  #scale_color_gradient(low ="#F8766D", high = "#619CFF") +
  #scale_color_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  #geom_line(aes(y=logFit3), color="blue", size=1) +
  theme_bw() +
  scale_x_reverse(breaks=c(3500, 2500, 1500, 500), limits=c(3500,400))+
  scale_y_continuous(limits=c(-.2,0.25))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Years Cal BP", y="KDE per capita growth", title = "E. MW ArchaeoGlobe KDE Per Capita Growth")+
  geom_line(size=1.25)+
  geom_line(aes(x=calBP, y=PerCap),size=1.25, color="red")+
  geom_hline(yintercept = 0)
mwpcs

pdf("SI/MWPerCap.pdf", width=8, height=6.55)
mwpcs
dev.off()

####Northeast Logistic=======================================================================================


SPD.3<- subset(SPD, Archaeo_ID=="19" )
cptcal <- calibrate(x = SPD.3$Age,  errors = SPD.3$Error, calCurves = "intcal20",  normalised = FALSE)
# <- binPrep(sites = box2$Trinomial, ages = box2$DatesUse, h = 100)


#Create SPD
spd.ne<- spd(cptcal, bins=NA, runm=100, timeRange=c(3500,200))
plot(spd.ne, runm=100, xlim=c(3500,200), type="simple")

PrDens<-spd.ne$grid$PrDens
calBP<-spd.ne$grid$calBP

NE.randates = sampleDates(cptcal,bins=NA,nsim=200,verbose=FALSE)
D.ckde = ckde(NE.randates,timeRange=c(3500,200),bw=50, normalised = FALSE)
plot(D.ckde,type='multiline')

##Write matrix of KDEs as a data frame
Check<-as.data.frame(D.ckde$res.matrix)

##Cbind spd and KDEs and write
dd<-cbind(calBP,PrDens, Check)
write.table(dd, file = "data/AGNEKDE50.csv", sep = ",", col.names=NA)

###At this point, I open the spreadsheet, and I add a two new columns: MKDE and PerCap.
#MKDE I calculate as the average of the KDE values at each time step. So, in this case, I add the KDE at cal BP 3500 from V1-V200
#and then divide by 200. This gives the mean KDE value at 3500 cal BP. Once I have done this for all time points, I calculate
#per capita growth as LN(MKDE at t+1) -LN(MKDE at time t). Once this is done, I load the csv file back into R.

#Load KDE data frame with calculated mean KDE and per capita growth
dd2m<- read.csv("data/AGNEKDE50.csv")
#I then trimm time points with missing values at the end due to KDE bandwidth.
dd2m<-na.omit(dd2m)

###transform dd2 into long table format to plot all of the simulated KDEs
dd3m<-dd2m %>% pivot_longer(cols=c(V1:V200),names_to='Runs', values_to='KDE')

##OK, now we fit a logistic model to the mean KDE. This can be a bit tricky due to the start parameters.
#These need to be adjusted based on the MKDE's values. In this case, the start values work well.
#fit logistic to KDE
logFitKDm <- nls(MKDE~SSfpl(calBP, A, B, xmid, scale), data=dd2m,control=nls.control(maxiter=700),start=list(A=.00005, B=0.00015, xmid=1200, scale=-100))
summary(logFitKDm)
predictm<-predict(logFitKDm)

##Combine predicted logistic values and KDE table for plotting
dd4m<-cbind(dd2m, predictm)
#write.table(dd4, file = "SWKDE50.csv", sep = ",", col.names=NA)

###Simple plot
plot(dd2m$MKDE~dd2m$calBP)
lines(dd2m$calBP, dd2m$MKDE,col="green",lty=2,lwd=3)
lines(dd2m$calBP,predict(logFitKDm),col="red",lty=2,lwd=3)

###Plot mean KDE and all of the simulated KDE models using ggplot
p6<-ggplot(data = dd3m) +
  geom_line(aes(x = calBP, y = (KDE*100), group = Runs), color = "purple", alpha = 0.5) +
  geom_point(data = dd4m, aes(x = calBP, y = (MKDE*100)), color="black", size=2) +
  geom_line(data = dd2m, aes(x = calBP, y = (predictm*100)), linetype = "dashed", color="green", size=1.5) +
  scale_x_reverse(breaks=c(3500, 2500, 1500, 500), limits=c(3500,400))+
  theme_bw() +
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Years Cal BP", y="KDE of radiocarbon ages", title = "F. NE ArchaeoGlobe KDE and Logistic Model")
# annotate("text", x =3500, y = .047, label = "Phase 1", size = 7, angle=90)+
#annotate("text", x =2000, y = .047, label = "Phase 2", size = 7, angle=90)+
#annotate("text", x =1020, y = .047, label = "Phase 3", size = 7, angle=90)+
#annotate("text", x =210, y = .047, label = "Phase 4", size = 7, angle=90)+
# facet_wrap(factor(Region)~.)
p6

pdf("SI/NEPop.pdf", width=8, height=6.55)
p6
dev.off()

dd2c<- read.csv("data/AGNEKDE50.csv")
#I then trimm time points with missing values at the end due to KDE bandwidth.
dd2c<-na.omit(dd2c)
dd2c<-subset(dd2c, select=-c(V1:V200))
###Sum SPD by 30 year intervals
library(zoo)
# sum and save new csvs.
out50 <- rollapply(dd2c,30,(sum),by=30,by.column=TRUE,align='right')
write.table(out50, file = "sums/AGNEKDE30.csv", sep = ",", col.names=NA)

pc30ne<- read.csv("sums/AGNEKDE30.csv")

nepc <- ggplot(pc30ne,aes(x=(calBP), y=(PerCap))) +
  #geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70") +
  geom_point(aes(), size=1) +
  #scale_color_gradient(low ="#F8766D", high = "#619CFF") +
  #scale_color_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  #geom_line(aes(y=logFit3), color="blue", size=1) +
  theme_bw() +
  scale_x_reverse(breaks=c(3500, 2500, 1500, 500), limits=c(3500,300))+
  scale_y_continuous(limits=c(-.2,0.25))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Years Cal BP", y="KDE per capita growth", title = "E. MW ArchaeoGlobe KDE Per Capita Growth")+
  geom_line(size=1.25)+
  geom_hline(yintercept = 0)
nepc

nepcs <- ggplot(pc30ne,aes(x=(calBP), y=(spdPerCap))) +
  geom_point(aes(), size=1) +
  #scale_color_gradient(low ="#F8766D", high = "#619CFF") +
  #scale_color_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  #geom_line(aes(y=logFit3), color="blue", size=1) +
  theme_bw() +
  scale_x_reverse(breaks=c(3500, 2500, 1500, 500), limits=c(3500,400))+
  scale_y_continuous(limits=c(-.2,0.25))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Years Cal BP", y="KDE per capita growth", title = "F. NE ArchaeoGlobe KDE Per Capita Growth")+
  geom_line(size=1.25)+
  geom_line(aes(x=calBP, y=PerCap),size=1.25, color="red")+
  geom_hline(yintercept = 0)
nepcs

pdf("SI/NEPerCap.pdf", width=8, height=6.55)
nepcs
dev.off()

####Southern Logistic=========================================================

SPD.3<- subset(SPD, Archaeo_ID=="20" )
##calibrate
cptcal <- calibrate(x = SPD.3$Age,  errors = SPD.3$Error, calCurves = "intcal20",  normalised = FALSE)
# <- binPrep(sites = box2$Trinomial, ages = box2$DatesUse, h = 100)

spd.s<- spd(cptcal, bins=NA, runm=100, timeRange=c(3500,200))
plot(spd.s, runm=100, xlim=c(3500,200), type="simple")

PrDens<-spd.s$grid$PrDens
calBP<-spd.s$grid$calBP

SS.randates = sampleDates(cptcal,bins=NA,nsim=200,verbose=FALSE)
D.ckde = ckde(SS.randates,timeRange=c(3500,200),bw=50, normalised = FALSE)
plot(D.ckde,type='multiline')

##Write matrix of KDEs as a data frame
Check<-as.data.frame(D.ckde$res.matrix)

##Cbind spd and KDEs and write
dd<-cbind(calBP,PrDens, Check)
write.table(dd, file = "data/AGSSKDE50.csv", sep = ",", col.names=NA)

###At this point, I open the spreadsheet, and I add a two new columns: MKDE and PerCap.
#MKDE I calculate as the average of the KDE values at each time step. So, in this case, I add the KDE at cal BP 3500 from V1-V200
#and then divide by 200. This gives the mean KDE value at 3500 cal BP. Once I have done this for all time points, I calculate
#per capita growth as LN(MKDE at t+1) -LN(MKDE at time t). Once this is done, I load the csv file back into R.

#Load KDE data frame with calculated mean KDE and per capita growth
dd2m<- read.csv("data/AGSSKDE50.csv")
#I then trimm time points with missing values at the end due to KDE bandwidth.
dd2m<-na.omit(dd2m)

###transform dd2 into long table format to plot all of the simulated KDEs
dd3m<-dd2m %>% pivot_longer(cols=c(V1:V200),names_to='Runs', values_to='KDE')

##OK, now we fit a logistic model to the mean KDE. This can be a bit tricky due to the start parameters.
#These need to be adjusted based on the MKDE's values. In this case, the start values work well.
#fit logistic to KDE
logFitKDm <- nls(MKDE~SSfpl(calBP, A, B, xmid, scale), data=dd2m,control=nls.control(maxiter=700),start=list(A=.00005, B=0.00015, xmid=1200, scale=-100))
summary(logFitKDm)
predictm<-predict(logFitKDm)

##Combine predicted logistic values and KDE table for plotting
dd4m<-cbind(dd2m, predictm)
#write.table(dd4, file = "SWKDE50.csv", sep = ",", col.names=NA)

###Simple plot
plot(dd2m$MKDE~dd2m$calBP)
lines(dd2m$calBP, dd2m$MKDE,col="green",lty=2,lwd=3)
lines(dd2m$calBP,predict(logFitKDm),col="red",lty=2,lwd=3)

###Plot mean KDE and all of the simulated KDE models using ggplot
p7<-ggplot(data = dd3m) +
  geom_line(aes(x = calBP, y = (KDE*100), group = Runs), color = "purple", alpha = 0.5) +
  geom_point(data = dd4m, aes(x = calBP, y = (MKDE*100)), color="black", size=2) +
  geom_line(data = dd2m, aes(x = calBP, y = (predictm*100)), linetype = "dashed", color="green", size=1.5) +
  scale_x_reverse(breaks=c(3500, 2500, 1500, 500), limits=c(3500,400))+
  theme_bw() +
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Years Cal BP", y="KDE of radiocarbon ages", title = "G. SS ArchaeoGlobe KDE and Logistic Model")
# annotate("text", x =3500, y = .047, label = "Phase 1", size = 7, angle=90)+
#annotate("text", x =2000, y = .047, label = "Phase 2", size = 7, angle=90)+
#annotate("text", x =1020, y = .047, label = "Phase 3", size = 7, angle=90)+
#annotate("text", x =210, y = .047, label = "Phase 4", size = 7, angle=90)+
# facet_wrap(factor(Region)~.)
p7

pdf("SI/SSPop.pdf", width=8, height=6.55)
p7
dev.off()

dd2c<- read.csv("data/AGSSKDE50.csv")
#I then trimm time points with missing values at the end due to KDE bandwidth.
dd2c<-na.omit(dd2c)
dd2c<-subset(dd2c, select=-c(V1:V200))
###Sum SPD by 30 year intervals
library(zoo)
# sum and save new csvs.
out50 <- rollapply(dd2c,30,(sum),by=30,by.column=TRUE,align='right')
write.table(out50, file = "sums/AGSSKDE30.csv", sep = ",", col.names=NA)

pc30ss<- read.csv("sums/AGSSKDE30.csv")

sspc <- ggplot(pc30ss,aes(x=(calBP), y=(PerCap))) +
  #geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70") +
  geom_point(aes(), size=1) +
  #scale_color_gradient(low ="#F8766D", high = "#619CFF") +
  #scale_color_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  #geom_line(aes(y=logFit3), color="blue", size=1) +
  theme_bw() +
  scale_x_reverse(breaks=c(3500, 2500, 1500, 500), limits=c(3500,300))+
  scale_y_continuous(limits=c(-.2,0.25))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Years Cal BP", y="KDE per capita growth", title = "E. SS ArchaeoGlobe KDE Per Capita Growth")+
  geom_line(size=1.25)+
  geom_hline(yintercept = 0)
sspc

sspcs <- ggplot(pc30ss,aes(x=(calBP), y=(spdPerCap))) +
  geom_point(aes(), size=1) +
  #scale_color_gradient(low ="#F8766D", high = "#619CFF") +
  #scale_color_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  #geom_line(aes(y=logFit3), color="blue", size=1) +
  theme_bw() +
  scale_x_reverse(breaks=c(3500, 2500, 1500, 500), limits=c(3500,400))+
  scale_y_continuous(limits=c(-.2,0.25))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Years Cal BP", y="KDE per capita growth", title = "G. SS ArchaeoGlobe KDE Per Capita Growth")+
  geom_line(size=1.25)+
  geom_line(aes(x=calBP, y=PerCap),size=1.25, color="red")+
  geom_hline(yintercept = 0)
sspcs

pdf("SI/SSPerCap.pdf", width=8, height=6.55)
sspcs
dev.off()


###Run logistic for the Rocky Mountain===============================================================================
SPD.1<- subset(SPD, Archaeo_ID=="16")

##calibrate
cptcal <- calibrate(x = SPD.1$Age,  errors = SPD.1$Error, calCurves = "intcal20",  normalised = FALSE)
# <- binPrep(sites = box2$Trinomial, ages = box2$DatesUse, h = 100)

####Run analysis for component 3 logistic 3400 to 250
spd.rm <- spd(cptcal, bins=NA, runm=100, timeRange=c(3500,200))
plot(spd.rm, runm=100, xlim=c(3500,200), type="simple")


PrDens<-spd.rm$grid$PrDens
calBP<-spd.rm$grid$calBP

RM.randates = sampleDates(cptcal,bins=NA,nsim=200,verbose=FALSE)
D.ckde = ckde(RM.randates,timeRange=c(3500,200),bw=50, normalised = FALSE)
plot(D.ckde,type='multiline')

##Write matrix of KDEs as a data frame
Check<-as.data.frame(D.ckde$res.matrix)

##Cbind spd and KDEs and write
dd<-cbind(calBP,PrDens, Check)
write.table(dd, file = "data/AGRMKDE50.csv", sep = ",", col.names=NA)

###At this point, I open the spreadsheet, and I add a two new columns: MKDE and PerCap.
#MKDE I calculate as the average of the KDE values at each time step. So, in this case, I add the KDE at cal BP 3500 from V1-V200
#and then divide by 200. This gives the mean KDE value at 3500 cal BP. Once I have done this for all time points, I calculate
#per capita growth as LN(MKDE at t+1) -LN(MKDE at time t). Once this is done, I load the csv file back into R.

#Load KDE data frame with calculated mean KDE and per capita growth
dd2m<- read.csv("data/AGRMKDE50.csv")
#I then trimm time points with missing values at the end due to KDE bandwidth.
dd2m<-na.omit(dd2m)

###transform dd2 into long table format to plot all of the simulated KDEs
dd3m<-dd2m %>% pivot_longer(cols=c(V1:V200),names_to='Runs', values_to='KDE')

##OK, now we fit a logistic model to the mean KDE. This can be a bit tricky due to the start parameters.
#These need to be adjusted based on the MKDE's values. In this case, the start values work well.
#fit logistic to KDE
logFitKDm <- nls(MKDE~SSfpl(calBP, A, B, xmid, scale), data=dd2m,control=nls.control(maxiter=700),start=list(A=.00005, B=0.00015, xmid=1200, scale=-100))
summary(logFitKDm)
predictm<-predict(logFitKDm)

##Combine predicted logistic values and KDE table for plotting
dd4m<-cbind(dd2m, predictm)
#write.table(dd4, file = "SWKDE50.csv", sep = ",", col.names=NA)

###Simple plot
plot(dd2m$MKDE~dd2m$calBP)
lines(dd2m$calBP, dd2m$MKDE,col="green",lty=2,lwd=3)
lines(dd2m$calBP,predict(logFitKDm),col="red",lty=2,lwd=3)

###Plot mean KDE and all of the simulated KDE models using ggplot
p8<-ggplot(data = dd3m) +
  geom_line(aes(x = calBP, y = (KDE*100), group = Runs), color = "purple", alpha = 0.5) +
  geom_point(data = dd4m, aes(x = calBP, y = (MKDE*100)), color="black", size=2) +
  geom_line(data = dd2m, aes(x = calBP, y = (predictm*100)), linetype = "dashed", color="green", size=1.5) +
  scale_x_reverse(breaks=c(3500, 2500, 1500, 500), limits=c(3500,400))+
  theme_bw() +
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Years Cal BP", y="KDE of radiocarbon ages", title = "H. RM ArchaeoGlobe KDE and Logistic Model")
# annotate("text", x =3500, y = .047, label = "Phase 1", size = 7, angle=90)+
#annotate("text", x =2000, y = .047, label = "Phase 2", size = 7, angle=90)+
#annotate("text", x =1020, y = .047, label = "Phase 3", size = 7, angle=90)+
#annotate("text", x =210, y = .047, label = "Phase 4", size = 7, angle=90)+
# facet_wrap(factor(Region)~.)
p8

pdf("SI/RMPop.pdf", width=8, height=6.55)
p8
dev.off()

dd2c<- read.csv("data/AGRMKDE50.csv")
#I then trimm time points with missing values at the end due to KDE bandwidth.
dd2c<-na.omit(dd2c)
dd2c<-subset(dd2c, select=-c(V1:V200))
###Sum SPD by 30 year intervals
library(zoo)
# sum and save new csvs.
out50 <- rollapply(dd2c,30,(sum),by=30,by.column=TRUE,align='right')
write.table(out50, file = "sums/AGRMKDE30.csv", sep = ",", col.names=NA)

pc30rm<- read.csv("sums/AGRMKDE30.csv")

rmpc <- ggplot(pc30rm,aes(x=(calBP), y=(PerCap))) +
  #geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70") +
  geom_point(aes(), size=1) +
  #scale_color_gradient(low ="#F8766D", high = "#619CFF") +
  #scale_color_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  #geom_line(aes(y=logFit3), color="blue", size=1) +
  theme_bw() +
  scale_x_reverse(breaks=c(3500, 2500, 1500, 500), limits=c(3500,300))+
  scale_y_continuous(limits=c(-.2,0.25))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Years Cal BP", y="KDE per capita growth", title = "E. SS ArchaeoGlobe KDE Per Capita Growth")+
  geom_line(size=1.25)+
  geom_hline(yintercept = 0)
rmpc

rmpcs <- ggplot(pc30rm,aes(x=(calBP), y=(spdPerCap))) +
  geom_point(aes(), size=1) +
  #scale_color_gradient(low ="#F8766D", high = "#619CFF") +
  #scale_color_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  #geom_line(aes(y=logFit3), color="blue", size=1) +
  theme_bw() +
  scale_x_reverse(breaks=c(3500, 2500, 1500, 500), limits=c(3500,400))+
  scale_y_continuous(limits=c(-.2,0.25))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Years Cal BP", y="KDE per capita growth", title = "H. RM ArchaeoGlobe KDE Per Capita Growth")+
  geom_line(size=1.25)+
  geom_line(aes(x=calBP, y=PerCap),size=1.25, color="red")+
  geom_hline(yintercept = 0)
rmpcs

pdf("SI/RMPerCap.pdf", width=8, height=6.55)
rmpcs
dev.off()

###Transitions histograms

Tran<-read.csv(file="data/Transitions.csv", header=T)
Tran2<-subset(Tran, ScaleID==1 & KDEID==1 & Transition.Type=="Positive")
Tran4<-subset(Tran, ScaleID==1 & KDEID==1 & Transition.Type=="Negative")
Tran3<-read.csv(file="data/TransitionsSum.csv", header=T)

pcar1 <- ggplot(Tran4, aes(factor(ID), (Peak), fill=ID))+
  geom_violin()+
  # geom_boxplot(notch = FALSE, width=0.2)+
  stat_summary(fun=median, geom="point", size=4, color="black")+
  stat_boxplot(geom ='errorbar')+
  #facet_wrap( ~ factor(CompID1))+
  #scale_fill_manual(values=c("#FC4E07", "#00A4CCFF"))+
  theme_bw() +
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28),  plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Subsistence ID", y="Valley of net growth periods", title = "B. Valley of Negative Growth and Agriculture")+
  annotate("text", x = 1.5, y = -.10, label = "W=184", size = 8)+
  annotate("text", x =2.55, y = -.10, label = "W=53", size = 8)+
  annotate("text", x =2.1, y = -.15, label = "W=87**", size = 8)
#facet_wrap(~Transition.Type)
pcar1

pcar1g <- ggplot(Tran4, aes(factor(ID), (Generations), fill=ID))+
  geom_violin()+
  # geom_boxplot(notch = FALSE, width=0.2)+
  stat_summary(fun=median, geom="point", size=4, color="black")+
  stat_boxplot(geom ='errorbar')+
  #facet_wrap( ~ factor(CompID1))+
  #scale_fill_manual(values=c("#FC4E07", "#00A4CCFF"))+
  theme_bw() +
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28),  plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Subsistence ID", y="Generations of negative growth", title = "B. Generations of Negative Growth and Agriculture")+
  annotate("text", x = 1.5, y = 12, label = "W=122", size = 8)+
  annotate("text", x =2.55, y = 12, label = "W=18.5", size = 8)+
  annotate("text", x =2.1, y = 15, label = "W=25.2**", size = 8)
#facet_wrap(~Transition.Type)
pcar1g


pcar2 <- ggplot(Tran2, aes(factor(ID), (Generations), fill=ID))+
  geom_violin()+
  stat_summary(fun=median, geom="point", size=4, color="black")+
  stat_boxplot(geom ='errorbar')+
  #facet_wrap( ~ factor(CompID1))+
  #scale_fill_manual(values=c("#FC4E07", "#00A4CCFF"))+
  theme_bw() +
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28),  plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Subsistence ID", y="Generations of growth", title = "A. Generations of Positive Growth and Agriculture")+
  annotate("text", x = 1.5, y = 36, label = "W=98.5", size = 8)+
  annotate("text", x =2.55, y = 36, label = "W=7**", size = 8)+
  annotate("text", x =2.1, y = 45, label = "W=4.5**", size = 8)
#facet_wrap(~Transition.Type)
pcar2

pcar3 <- ggplot(Tran2, aes(factor(ID), (Peak), fill=ID))+
  geom_violin()+
  stat_summary(fun=median, geom="point", size=4, color="black")+
  stat_boxplot(geom ='errorbar')+
  #facet_wrap( ~ factor(CompID1))+
  #scale_fill_manual(values=c("#FC4E07", "#00A4CCFF"))+
  theme_bw() +
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28),  plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Subsistence ID", y="Generations of growth", title = "C. Peak Positive Growth and Agriculture")+
  annotate("text", x = 1.5, y = .12, label = "W=128", size = 8)+
  annotate("text", x =2.55, y = .12, label = "W=10**", size = 8)+
  annotate("text", x =2.1, y = .15, label = "W=22**", size = 8)
#facet_wrap(~Transition.Type)
pcar3

Tran4<-subset(Tran2, ID %in% c('HG','Min Agg'))
wilcox.test(Peak ~ ID, data=Tran4)

Tran5<-subset(Tran2, ID %in% c('Min Agg','Moderate Agg'))
wilcox.test(Peak ~ ID, data=Tran5)

Tran6<-subset(Tran2, ID %in% c('HG','Moderate Agg'))
wilcox.test(Peak ~ ID, data=Tran6)


Tran7<-subset(Tran4, ID %in% c('HG','Min Agg'))
wilcox.test(Generations ~ ID, data=Tran7)

Tran8<-subset(Tran4, ID %in% c('Min Agg','Moderate Agg'))
wilcox.test(Generations ~ ID, data=Tran8)

Tran9<-subset(Tran4, ID %in% c('HG','Moderate Agg'))
wilcox.test(Generations ~ ID, data=Tran9)

Fig5Rev<-plot_grid(pcar2,pcar3,pcar1g,pcar1, ncol=2, align="hv", axis = "rl")
Fig5Rev

pdf("data/DT_USA3.pdf", width=13, height=11.55)
Fig5Rev
dev.off()