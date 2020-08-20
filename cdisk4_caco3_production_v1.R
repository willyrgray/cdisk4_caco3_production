#error propagation to calculate caco3 production from standing stocks and turnover times
#ziveri et al, submitted to nat. comms.

#vioplot must be installed, required for figure
require(vioplot) 

#station data
stations<- c(1,2,3,4,5)
lat<- c(22.754, 27.733, 35.267, 41.750, 49.833)

#set aesthetics
st_cols<- c("firebrick3","darkorange2","seagreen", "dodgerblue3","darkorchid3")
st_pch<- c(1,2,3,4,5)

groups<- c('ptero', 'hetero','foram','cocco')
gr_cols<- c("darkorchid3","mediumpurple3","cyan3", "cadetblue")
gr_pch<- c(17,17,16,16)


#import standing stock data
stdstk<- read.csv('stdstk.csv') #note this is unpublished data

ptero<- stdstk[,2]; ptero_sigma<- stdstk[,6]
hetero<- stdstk[,3]; hetero_sigma<- stdstk[,7]
foram<- stdstk[,4]; foram_sigma<- stdstk[,8]
coco<- stdstk[,5]; coco_sigma<- stdstk[,9]

calcite<- foram + coco
calcite_sigma<- foram_sigma + coco_sigma
aragonite<- ptero + hetero
aragonite_sigma<- ptero_sigma + hetero_sigma
tot<- calcite + aragonite
tot_sigma<- calcite_sigma + aragonite_sigma


###
#stdstk to production 

#turnover time
#set min and max turnover time in days
turnovr<-data.frame(group=c('ptero','hetero','foram','coco'),min=c(5,5,10,1/1.5), max=c(16,16,30,1/0.1))

#seasonal bias satellite PIC
satellite_pic <- read.csv('seawifs_pic_monthly_clim.csv') #seawifs pic monthly clim for each station
satellite_pic2017<- read.csv('seawifs_pic_2017.csv') #seawifs pic aug 2017 for each station
satellite_pic_annual_mean<- as.numeric(colMeans(satellite_pic, na.rm=TRUE))
satellite_pic_aug2017<- as.numeric(satellite_pic2017[8,])
seasonal_bias<- satellite_pic_aug2017/satellite_pic_annual_mean #ratio of pic Aug 2017 compared to annual climatology

#seasonal bias - zooplankton timeseries
seasonal_bias_zoo_HOTS_PAPA<- c(1.2,2) #ratio of aug/annal mean climatological zooplankton abundance from HOTS/PAPA
lat_HOTS_PAPA<-c(22.750, 50.100)
seasonal_bias_zoo<- predict(lm(seasonal_bias_zoo_HOTS_PAPA~lat_HOTS_PAPA), data.frame(lat_HOTS_PAPA=lat)) #extrapolate based on lat


#monte carlo simulation
#normal probability distributions for stdstk, flat probability distributions for turnovr time
mcmc_iterations<- 9999

p<-matrix(, nrow = mcmc_iterations, ncol = length(stations))
h<-matrix(, nrow = mcmc_iterations, ncol = length(stations))
f<-matrix(, nrow = mcmc_iterations, ncol = length(stations))
c<-matrix(, nrow = mcmc_iterations, ncol = length(stations))

ca<-matrix(, nrow = mcmc_iterations, ncol = length(stations))
a<-matrix(, nrow = mcmc_iterations, ncol = length(stations))
s<-matrix(, nrow = mcmc_iterations, ncol = length(stations))


p1<-matrix(, nrow = mcmc_iterations, ncol = length(stations))
h1<-matrix(, nrow = mcmc_iterations, ncol = length(stations))
f1<-matrix(, nrow = mcmc_iterations, ncol = length(stations))
c1<-matrix(, nrow = mcmc_iterations, ncol = length(stations))

ca1<-matrix(, nrow = mcmc_iterations, ncol = length(stations))
a1<-matrix(, nrow = mcmc_iterations, ncol = length(stations))
s1<-matrix(, nrow = mcmc_iterations, ncol = length(stations))

for(i in 1:mcmc_iterations){
	pi<- rnorm(length(ptero),ptero,ptero_sigma)*(365/runif(1,min=turnovr$min[1], max=turnovr$max[1]))
	hi<- rnorm(length(hetero),hetero,hetero_sigma)*(365/runif(1,min=turnovr$min[2], max=turnovr$max[2]))
	fi<- rnorm(length(foram),foram,foram_sigma)*(365/runif(1,min=turnovr$min[3], max=turnovr$max[3]))
	ci<- rnorm(length(coco),coco,coco_sigma)*(365/runif(1,min=turnovr$min[4], max=turnovr$max[4]))
		
	p[i,]<- pi
	h[i,]<- hi
	f[i,]<- fi
	c[i,]<- ci
	
	a[i,]<- pi + hi
	ca[i,]<- ci + fi
	s[i,]<- pi + hi + ci +fi

	#seasonal bias correction
	sbiasi<- seasonal_bias #for cocco/forams
	sbiasi_zoo<- seasonal_bias_zoo #for hetero/pteropods
		
	p1[i,]<- pi/sbiasi_zoo
	h1[i,]<- hi/sbiasi_zoo
	f1[i,]<- fi/sbiasi
	c1[i,]<- ci/sbiasi
		
	ca1[i,]<- (fi + ci)/sbiasi
	a1[i,]<- (pi + hi)/sbiasi_zoo
	s1[i,]<- ((fi + ci)/sbiasi + (pi + hi)/sbiasi_zoo)
}

#annual
ptero_annual<- apply(p, 2, quantile, probs=0.5, na.rm=TRUE); ptero_annual_lwr <- apply(p, 2, quantile, probs=0.025, na.rm=TRUE);ptero_annual_upr <- apply(p, 2, quantile, probs=0.975, na.rm=TRUE)
hetero_annual<- apply(h, 2, quantile, probs=0.5, na.rm=TRUE); hetero_annual_lwr <- apply(h, 2, quantile, probs=0.025, na.rm=TRUE);hetero_annual_upr <- apply(h, 2, quantile, probs=0.975, na.rm=TRUE)
foram_annual<- apply(f, 2, quantile, probs=0.5, na.rm=TRUE); foram_annual_lwr <- apply(f, 2, quantile, probs=0.025, na.rm=TRUE);foram_annual_upr <- apply(f, 2, quantile, probs=0.975, na.rm=TRUE)
coco_annual<- apply(c, 2, quantile, probs=0.5, na.rm=TRUE); coco_annual_lwr <- apply(c, 2, quantile, probs=0.025, na.rm=TRUE);coco_annual_upr <- apply(c, 2, quantile, probs=0.975, na.rm=TRUE) 

calcite_annual<- apply(ca, 2, quantile, probs=0.5, na.rm=TRUE); calcite_annual_lwr <- apply(ca, 2, quantile, probs=0.025, na.rm=TRUE);calcite_annual_upr <- apply(ca, 2, quantile, probs=0.975, na.rm=TRUE)
aragonite_annual<- apply(a, 2, quantile, probs=0.5, na.rm=TRUE); aragonite_annual_lwr <- apply(a, 2, quantile, probs=0.025, na.rm=TRUE);aragonite_annual_upr <- apply(a, 2, quantile, probs=0.975, na.rm=TRUE)

tot_annual<- apply(s, 2, quantile, probs=0.5, na.rm=TRUE); tot_annual_lwr <- apply(s, 2, quantile, probs=0.025, na.rm=TRUE);tot_annual_upr <- apply(s, 2, quantile, probs=0.975, na.rm=TRUE)


annual<- data.frame(stations=stations, 
ptero_annual=ptero_annual,ptero_annual_lwr=ptero_annual_lwr,ptero_annual_upr=ptero_annual_upr,
hetero_annual=hetero_annual,hetero_annual_lwr=hetero_annual_lwr, hetero_annual_upr=hetero_annual_upr,
foram_annual=foram_annual,foram_annual_lwr=foram_annual_lwr,foram_annual_upr=foram_annual_upr,
coco_annual=coco_annual,coco_annual_lwr=coco_annual_lwr,coco_annual_upr=coco_annual_upr,
calcite_annual=calcite_annual, calcite_annual_lwr=calcite_annual_lwr, calcite_annual_upr=calcite_annual_upr,
aragonite_annual=aragonite_annual,aragonite_annual_lwr=aragonite_annual_lwr,aragonite_annual_upr=aragonite_annual_upr,
tot_annual=tot_annual, tot_annual_lwr=tot_annual_lwr,tot_annual_upr=tot_annual_upr)

#export annual values, uncorrectd for seasonality
write.table(annual,file='annual_mg_m2_yr.csv', sep=',',row.names=FALSE)



#annual, seasonally corrected
ptero_annual_cor<- apply(p1, 2, quantile, probs=0.5, na.rm=TRUE); ptero_annual_cor_lwr <- apply(p1, 2, quantile, probs=0.025, na.rm=TRUE);ptero_annual_cor_upr <- apply(p1, 2, quantile, probs=0.975, na.rm=TRUE)
hetero_annual_cor<- apply(h1, 2, quantile, probs=0.5, na.rm=TRUE); hetero_annual_cor_lwr <- apply(h1, 2, quantile, probs=0.025, na.rm=TRUE);hetero_annual_cor_upr <- apply(h1, 2, quantile, probs=0.975, na.rm=TRUE)
foram_annual_cor<- apply(f1, 2, quantile, probs=0.5, na.rm=TRUE); foram_annual_cor_lwr <- apply(f1, 2, quantile, probs=0.025, na.rm=TRUE);foram_annual_cor_upr <- apply(f1, 2, quantile, probs=0.975, na.rm=TRUE)
coco_annual_cor<- apply(c1, 2, quantile, probs=0.5, na.rm=TRUE); coco_annual_cor_lwr <- apply(c1, 2, quantile, probs=0.025, na.rm=TRUE);coco_annual_cor_upr <- apply(c1, 2, quantile, probs=0.975, na.rm=TRUE) 

calcite_annual_cor<- apply(ca1, 2, quantile, probs=0.5, na.rm=TRUE); calcite_annual_cor_lwr <- apply(ca1, 2, quantile, probs=0.025, na.rm=TRUE);calcite_annual_cor_upr <- apply(ca1, 2, quantile, probs=0.975, na.rm=TRUE)
aragonite_annual_cor<- apply(a1, 2, quantile, probs=0.5, na.rm=TRUE); aragonite_annual_cor_lwr <- apply(a1, 2, quantile, probs=0.025, na.rm=TRUE);aragonite_annual_cor_upr <- apply(a1, 2, quantile, probs=0.975, na.rm=TRUE)

tot_annual_cor<- apply(s1, 2, quantile, probs=0.5, na.rm=TRUE); tot_annual_cor_lwr <- apply(s1, 2, quantile, probs=0.025, na.rm=TRUE);tot_annual_cor_upr <- apply(s1, 2, quantile, probs=0.975, na.rm=TRUE)


annual_cor<- data.frame(stations=stations, 
ptero_annual_cor=ptero_annual_cor,ptero_annual_cor_lwr=ptero_annual_cor_lwr,ptero_annual_cor_upr=ptero_annual_cor_upr,
hetero_annual_cor=hetero_annual_cor,hetero_annual_cor_lwr=hetero_annual_cor_lwr, hetero_annual_cor_upr=hetero_annual_cor_upr,
foram_annual_cor=foram_annual_cor,foram_annual_cor_lwr=foram_annual_cor_lwr,foram_annual_cor_upr=foram_annual_cor_upr,
coco_annual_cor=coco_annual_cor,coco_annual_cor_lwr=coco_annual_cor_lwr,coco_annual_cor_upr=coco_annual_cor_upr,
calcite_annual_cor=calcite_annual_cor, calcite_annual_cor_lwr=calcite_annual_cor_lwr, calcite_annual_cor_upr=calcite_annual_cor_upr,
aragonite_annual_cor=aragonite_annual_cor,aragonite_annual_cor_lwr=aragonite_annual_cor_lwr,aragonite_annual_cor_upr=aragonite_annual_cor_upr,
tot_annual_cor=tot_annual_cor, tot_annual_cor_lwr=tot_annual_cor_lwr,tot_annual_cor_upr=tot_annual_cor_upr)

#covert to mol/m2/yr
annual_cor_mol<- annual_cor/100/1000
#export annual values, correctd for seasonality in mol/m2/yr
write.table(annual_cor_mol,file='annual_cor_mol_m2_yr.csv', sep=',',row.names=FALSE)




###
#plot!
#dimensions
h=5.7
w=5.7
p<- dev.new(width=w, height=h); 
par(mfrow=c(2,2))
par(ps = 10, cex = 1, cex.main = 1); par(mgp=c(2.25,0.5,0)); par(las=1);par(tck=-0.02)
par(mar=c(3.5,3.5,1,1),xpd=TRUE)


#stndstk
plot(stations, tot, xlab='station', ylab=expression(CaCO[3]~(mg/m^2)), col=adjustcolor('grey99',alpha=0.01), pch=st_pch, xlim=c(0.5,5.5), ylim=c(0,5250))

for(i in 1:length(stations)){
	for(j in 1:length(groups)){
		lines(c(stations[i]+0.375-0.15*j, stations[i]+0.375-0.15*j), c(stdstk[i,j+1]-stdstk[i,j+5], stdstk[i,j+1]+stdstk[i,j+5]), col=adjustcolor(gr_cols[j],alpha=0.8), pch=gr_pch[j])
		points(stations[i]+0.375-0.15*j, stdstk[i,j+1], col=adjustcolor(gr_cols[j],alpha=0.8), pch=gr_pch[j],cex=1.1)	
	}
}

for(i in 1:length(stations)){
		lines(c(stations[i], stations[i]), c(tot[i]-tot_sigma[i], tot[i]+tot_sigma[i]), col=adjustcolor('grey37',alpha=0.9)) 
}

points(stations, tot, col=adjustcolor('grey37',alpha=0.9), pch=15,cex=1.3)

legend('topleft', pch=c(gr_pch,15),col=c(gr_cols,'grey43'),c(groups,'total'), cex=0.9,bty='n') #add legend


#turnover time
par(mar=c(3.5,3.5,1,1),xpd=TRUE)
plot(-1, -1, xlab='', ylab='turnover time (days)', col=adjustcolor('grey99',alpha=0.01), type='p', ylim=c(0,30),xlim=c(0,5), axes=FALSE);
axis(2, at=seq(0,30,by=5),labels=TRUE)

for(i in 1:length(groups)){
rect(i-0.1,turnovr$min[i],i+0.1,turnovr$max[i], col=adjustcolor(gr_cols[i],alpha=0.8),lty=1, lwd=1.5, border=NA)
}
text(c(1,2,3,4),-5,c('ptero', 'hetero','foram','cocco'),srt=45,cex=0.9)
 

#daily production, uncorrected for seasonality
par(mar=c(3.5,3.5,1,1),xpd=FALSE)
plot(stations, tot_annual/365, xlab='station', ylab=expression(CaCO[3]~(mg/m^2/day)), col=adjustcolor('grey99',alpha=0.01), pch=st_pch, xlim=c(0.5,5.5), ylim=c(0,3500))

for(i in 1:length(stations)){
vioplot(s[,i]/365, col=adjustcolor('grey17', alpha=0.3), horizontal=FALSE, at=i, add=TRUE,lty=1,pchMed=15, colMed=adjustcolor('grey17', alpha=0.7),cexMed=1.5,lineCol=adjustcolor('grey17', alpha=0.5),border=NA)
}

for(i in 1:length(stations)){
	lines(c(stations[i]-0.225,stations[i]-0.225), c(ptero_annual_lwr[i]/365,ptero_annual_upr[i])/365,col=adjustcolor(gr_cols[1],alpha=0.9), lwd=0.9,lty=1)
		lines(c(stations[i]-0.225/2,stations[i]-0.225/2), c(hetero_annual_lwr[i]/365,hetero_annual_upr[i])/365,col=adjustcolor(gr_cols[2],alpha=0.9), lwd=0.9,lty=1)
				lines(c(stations[i]+0.225/2,stations[i]+0.225/2), c(foram_annual_lwr[i]/365,foram_annual_upr[i])/365,col=adjustcolor(gr_cols[3],alpha=0.9), lwd=0.9,lty=1)
								lines(c(stations[i]+0.225,stations[i]+0.225), c(coco_annual_lwr[i]/365,coco_annual_upr[i])/365,col=adjustcolor(gr_cols[4],alpha=0.9), lwd=0.9,lty=1)	
	}
points(stations-0.225, ptero_annual/365,col=adjustcolor(gr_cols[1],alpha=0.9), pch=17, cex=1.1)
points(stations-0.225/2, hetero_annual/365,col=adjustcolor(gr_cols[2],alpha=0.9), pch=17, cex=1.1)
points(stations+0.225/2, foram_annual/365,col=adjustcolor(gr_cols[3],alpha=0.9), pch=16, cex=1.1)
points(stations+0.225, coco_annual/365,col=adjustcolor(gr_cols[4],alpha=0.9), pch=16, cex=1.1)


#annual production, corected for seasonality
par(mar=c(3.5,3.5,1,1),xpd=FALSE)
plot(stations, tot_annual_cor/100/1000, xlab='station', ylab=expression(CaCO[3]~(mol/m^2/yr)), col=adjustcolor('grey99',alpha=0.01), pch=st_pch, xlim=c(0.5,5.5), ylim=c(0,5))

for(i in 1:length(stations)){
vioplot(s1[,i]/100/1000, col=adjustcolor('grey17', alpha=0.3), horizontal=FALSE, at=i, add=TRUE,lty=1,pchMed=15, colMed=adjustcolor('grey17', alpha=0.7),lineCol=adjustcolor('grey17', alpha=0.5),border=NA)
}

for(i in 1:length(stations)){
	lines(c(stations[i]-0.225,stations[i]-0.225), c(ptero_annual_cor_lwr[i]/100/1000,ptero_annual_cor_upr[i])/100/1000,col=adjustcolor(gr_cols[1],alpha=0.9), lwd=0.9,lty=1)
		lines(c(stations[i]-0.225/2,stations[i]-0.225/2), c(hetero_annual_cor_lwr[i]/100/1000,hetero_annual_cor_upr[i])/100/1000,col=adjustcolor(gr_cols[2],alpha=0.9), lwd=0.9,lty=1)
				lines(c(stations[i]+0.225/2,stations[i]+0.225/2), c(foram_annual_cor_lwr[i]/100/1000,foram_annual_cor_upr[i])/100/1000,col=adjustcolor(gr_cols[3],alpha=0.9), lwd=0.9,lty=1)
								lines(c(stations[i]+0.225,stations[i]+0.225), c(coco_annual_cor_lwr[i]/100/1000,coco_annual_cor_upr[i])/100/1000,col=adjustcolor(gr_cols[4],alpha=0.9), lwd=0.9,lty=1)
		
	}
points(stations-0.225, ptero_annual_cor/100/1000,col=adjustcolor(gr_cols[1],alpha=0.9), pch=17, cex=1.1)
points(stations-0.225/2, hetero_annual_cor/100/1000,col=adjustcolor(gr_cols[2],alpha=0.9), pch=17, cex=1.1)
points(stations+0.225/2, foram_annual_cor/100/1000,col=adjustcolor(gr_cols[3],alpha=0.9), pch=16, cex=1.1)
points(stations+0.225, coco_annual_cor/100/1000,col=adjustcolor(gr_cols[4],alpha=0.9), pch=16, cex=1.1)

text(c(1.5,1.5), c(4.9,4.5),c('seasonally','adjusted'), col='grey67')