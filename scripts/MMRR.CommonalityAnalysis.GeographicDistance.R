#=============================================================================================================#
# Script created by Glenn Seeholzer, contact at seeholzer.glenn@gmail.com
# Script created in version R 3.3.2 
# This script: runs univariate and multivariate matrix regressions testing for isolation-by-adaptation, 
#			isolation-by-environment, and isolation-by-distance in Cranioleuca antisiensis. Where 
#			geographic distance is used in place of dispersal distance.  
#			Runs a commonality analysis assesing the relative important of these processes in shaping 
#			genetic structure of C. antisiensis
#			Produces Figure S8 & Table S2
# Usage notes: run line by line
#=============================================================================================================#
#set working directory
setwd('~/Dropbox/CRAN/Chapter.2/Seeholzer&Brumfield_ms_revision1/DataArchival/')

#load packages and supporting functions
library(plyr)
library(adegenet)
library(StAMPP)
library(gdistance)
library(ecodist)
library(maptools)
source('scripts/supporting.functions/FUN.add.alpha.R', chdir = TRUE)
source('scripts/supporting.functions/MMRR.R', chdir = F)


####	import phenotypic and meta-data
data = read.delim('Data.Table.txt',stringsAsFactors=F)

####	import genetic data
gen = get(load('SNPdata/SNPdata.genind.Rdata'))

####	specify the populations
tmp = data[data$ID %in% indNames(gen), ]
tmp = tmp[match(indNames(gen),tmp$ID), ]
strata(gen) = data.frame(pop = tmp[,'Population'])
setPop(gen) = ~pop


####	Calculate pairwise population global Fst matrix
####	because of minor variations in each run of stamppFst, it is advisable to upload the pairwise Fst matrix used in the original publication to ensure exact duplication of results presented in paper. 
fst = read.delim('Fst.all.05outliersexcluded.txt')

####	however an fst matrix almost identical to that above can be generated with the following commands 
# outliers = readLines('BayEnv_outlier_loci.txt') #import list of outlier loci to exclude from analysis
# gen.new = gen[loc=!(locNames(gen) %in% outliers)] #remove outlier from genind object 
# source('scripts/supporting.functions/genind2stampp.R', chdir = TRUE)
# gdata = genind2stampp(gen.new)
# set.seed(1234)
# fst = stamppFst(gdata,nboots=1) #this step will take a long time to run

####	Clean Data
####	remove outgroups Cranioleuca curatata and Cranioleuca erythrops from all datasets
data = data[!(data$Population %in% c('erythrops','curtata')), ]
gen = gen[i = which(!grepl('cra.cur',indNames(gen)))]
fst = fst[rownames(fst) != 'curtata',colnames(fst) != 'curtata']

####	Create final phenotypic data set
data = data[data$Population %in% rownames(fst) & data$ID %in% indNames(gen), ]

#--------------------------------------------------------------------#
#--------------------------------------------------------------------#
#--------------------------------------------------------------------#
#------------ Create raw population means for predictors ------------#
#--------------------------------------------------------------------#
#--------------------------------------------------------------------#
#--------------------------------------------------------------------#


####	ENV 
####	import bioclim and enviro rasters
files <- list.files(path="bioclim.veg.layers", pattern='grd', full.names=TRUE)
BIOCLIM = stack(files[grep('bio\\d',files)]) #create raster stack of BIOCLIM rasters
envdata = data.frame(pop=data[,'Population'],as.data.frame(extract(BIOCLIM, data[,c('Longitude','Latitude')]))) #extract BIOCLIM data for each individual
tmp = aggregate(. ~ pop, data= envdata,mean,na.rm=T) #take mean for each locality
my.prc = prcomp(tmp[,-1], center=T, scale = T, na.omit=T) #do PCA
pc.values = predict(my.prc)[,1:ncol(tmp[,-1])]
pop.env.raw = data.frame(pc.values)
rownames(pop.env.raw) = tmp[,1]

####	ALT 
tmp = aggregate(. ~ Population, data= data[,c('Population','Bioclim.Elevation')],mean,na.rm=T)
pop.alt = as.matrix(tmp[,'Bioclim.Elevation'])
rownames(pop.alt) = tmp[,'Population']
pop.alt.raw = pop.alt[rownames(fst),] #put in same order as fst matrix

####	MASS 
####	remove juveniles and individuals without mass data
tmp = aggregate(. ~ Population, data= data[!is.na(data$Mass) & !(data$Age %in% c('j')),c('Population','Mass')],mean,na.rm=T) 
pop.mass = as.matrix(tmp[,'Mass'])
rownames(pop.mass) = tmp[,'Population']
pop.mass.raw = pop.mass[rownames(fst),] #put in same order as dgen.fst.all

####	Wing.length
####	remove juveniles and individuals without size data
tmp = aggregate(. ~ Population, data= data[!is.na(data$Wing.length) & !(data$Age %in% c('j')) , c('Population','Wing.length')],mean,na.rm=T)
pop.wing = as.matrix(tmp[,'Wing.length'])
rownames(pop.wing) = tmp[,'Population']
order = rownames(fst)[rownames(fst) %in% rownames(pop.wing)] #missing populations without morphological measurment 
pop.wing.raw = pop.wing[order,] #put in same order as fst

####	PLUM (PC1.plumage)
tmp = aggregate(. ~ Population, data= data[,c('Population','PC1.plumage')],mean,na.rm=T)
pop.plumage = as.matrix(tmp[,'PC1.plumage'])
rownames(pop.plumage) = tmp[,'Population']
order = rownames(fst)[rownames(fst) %in% rownames(pop.plumage)]
pop.plumage.raw = pop.plumage[order,] #put in same order as dgen.fst.all


####	GREAT CIRCLE DISTANCE MATRIX
tmp = aggregate(. ~ Population, data= data[,c('Population','Longitude','Latitude')],mean,na.rm=T)
pop.loc = as.matrix(tmp[,c('Longitude','Latitude')])
rownames(pop.loc) = tmp[,'Population']
pop.loc['CorBlanca',] = c(-77.50920,-9.345465) #Population CorBlanca is a population aggregated from two distinct collecting localities separated by 35 kilomters. Averaging their coordinates centers this pouplation in a glacier, so we've used the coordinates of the collecting locality with the most samples (n=4).
geo.raw = spDists(pop.loc,longlat=TRUE) ; colnames(geo.raw) = rownames(geo.raw) = rownames(pop.loc)


####	DISPERSAL DISTANCE MATRICES
####	This section of the script takes a long to run. I recommended loading the distances matrices from the text files in the working directory using the code below. 
cost.raw = read.delim('cost.raw.txt',stringsAsFactors=F)
cost.masked = read.delim('cost.masked.txt',stringsAsFactors=F)

# If you would like to run the code to produce these distance matrices, un-comment each line below and run it.
####	Shortcut for R Console on Mac: Highlight the entire section and press Command+Option+' to remove comments from a large block of text

# ####RAW MAXENT output
# r = raster('SDM/MAXENT_SDM_Cranioleuca_antisiensis.asc')
# proj4string(r) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")

# rnew = r
# tr = transition(rnew,mean,directions=16)
# trC = geoCorrection(tr,'c',scl=T)
# save(trC,file='outfiles/cost.trC.raw.Rdata')

# dcost = costDistance(trC,pop.loc)
# cost.raw = as.matrix(dcost)
# write.table(cost.raw,file='outfiles/cost.raw.txt',row.names=T,col.names=T,quote=F,sep='\t')


# ####ELEVATIONAL FILTER + QGIS EDITS
# puna = readShapeSpatial('SDM/TNC_CentralAndeanWetPuna.shp')
# rnew2 = mask(r,puna,inverse=T, updatevalue=0.01)

# tr = transition(rnew2,mean,directions=16)
# trC = geoCorrection(tr,'c',scl=T)
# save(trC,file='outfiles/cost.trC.masked.Rdata')

# dcost = costDistance(trC,pop.loc)
# cost.masked = as.matrix(dcost)
# write.table(cost.masked,file='outfiles/cost.masked.txt',row.names=T,col.names=T,quote=F,sep='\t')

#------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------#
#------------ Create distance matrices of population data for predictors ------------#
#------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------#

# Libraries
library(yhat) # Commonality analysis

#functions for formatting predictor matrices and standardizing the values so coefficients comparable as beta-weights
standardize = function(x){tmp=(x-mean(x,na.rm=T))/sd(x,na.rm=T);diag(tmp)=0;return(tmp)}

prep.standardize <- function(x){x=data.matrix(x)
					x=x[lower.tri(x, diag = FALSE)]
					x=(x-mean(x))/sqrt(var(x)) 
					x}

prep <- function(x){x=data.matrix(x)
					x=x[lower.tri(x, diag = FALSE)]
					x}

################################################################################################
####	Raw Distance Matrices
################################################################################################
mFst = as.matrix(fst)
mFst[upper.tri(mFst)] = mFst[lower.tri(mFst)]

mGeo = as.matrix(geo.raw)
mCost.raw = as.matrix(cost.raw)
mCost.masked = as.matrix(cost.masked)
mAlt = as.matrix(dist(pop.alt.raw))
mMass = as.matrix(dist(pop.mass.raw))
mEnv = as.matrix(dist(pop.env.raw))
mPlumage = as.matrix(dist(pop.plumage.raw))
mWing = as.matrix(dist(pop.wing.raw))

################################################################################################
####	Lists of distance matrices formatted for Multiple Matrix Regression on Distance Matrices (MMRR, Wang 2013)
################################################################################################
 
#for univariate and regressions
mDist = mGeo
predictors = list(GEO=mDist,ELE=mAlt,ENV=mEnv,MASS=mMass,PLUM=mPlumage)
predictors = lapply(predictors,standardize)


################################################################################################
####	Distance matrices formatted for Commonality Analysis
################################################################################################
####	for commonality analysis
FST = prep.standardize(mFst)
GEO = prep.standardize(mDist) 
ELE = prep.standardize(mAlt)
ENV = prep.standardize(mEnv)
MASS = prep.standardize(mMass)
WING = prep.standardize(mWing)
PLUM = prep.standardize(mPlumage)



#-------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------#
#------------ Run Univariate, Multivariate and Commonality Analyses ------------#
#--------------------------------ALL POPULATIONS--------------------------------#
#-------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------#



############################
############################
#	Table 4 - Multivariate Regression and Commonality analysis For Body Mass
############################
############################
multi.table = c()
resBootstrap.list = list()

#MMRR to get regression coefficients and significance values for multivariate model
#data for MMRR
x = predictors[!(names(predictors) %in% c('SIZE','PLUM'))] #can't include plumage, besides, explains almost none of the variance
fit = MMRR(Y = mFst, X = x, nperm=1000)

#Do Commonality Analysis to get Unique, Common, and Total contribution of each predictor variable
data=data.frame(FST = FST, GEO = GEO, ELE = ELE, ENV = ENV, MASS = MASS)
ca = regr(lm(FST ~ GEO + MASS + ELE + ENV,data=data))

#add in the common
colnames = c('model','predictor','coefficient','tstatistic','pvalue','Unique','Common','Total')
table = data.frame(matrix(ncol=length(colnames),nrow=length(x)))
colnames(table) = colnames

model = paste0('Fst ~ ',paste(names(x),collapse=' + '))
table[,'model'] = c(model,paste0('r2 = ',round(fit$r.squared,2)),'','')
table[,'predictor'] = names(x)
table[,'coefficient'] = round(fit$coefficients[-1],2)
table[,'tstatistic'] = round(fit$tstatistic[-1],2)
table[,'pvalue'] = round(fit$tpvalue[-1],2)
table[,'Unique'] = round(ca$Commonality_Data$CCTotalbyVar[,'Unique'],2)
table[,'Common'] = round(ca$Commonality_Data$CCTotalbyVar[,'Common'],2)
table[,'Total']  = round(ca$Commonality_Data$CCTotalbyVar[,'Total'],2)
table = apply(table,2,as.character)

multi.table = rbind(multi.table,table,rep('',ncol(table)))

#this corresponds to the upper table for all populations in Table 4
table.path.multi = 'outfiles/Table S2.txt'
write.table(multi.table,file = table.path.multi,quote=F,col.names=T,row.names=F,sep='\t')



############################
############################
#Commonality analysis For Body Mass with confidence intervals
############################
############################
#bootstrap procedure modified from code in supplementary of Prunier et al. (2015), which was based on methods in Peterman et al. (2014)

#	Prunier JG, Colyn M, Legendre X, Nimon KF, Flamand MC (2015) Multicollinearity in spatial genetics: Separating the wheat from the chaff using commonality analyses. Molecular Ecology, 24, 263–283.
#	Peterman WE, Connette GM, Semlitsch RD, Eggert LS (2014) Ecological resistance surfaces predict fine-scale genetic differentiation in a terrestrial woodland salamander. Molecular Ecology, 23, 2402–2413.

nperm=1000
n.predictors = 4
ncombos = ((2^n.predictors)-1)
boot=matrix(data = 0, nrow = nperm, ncol = ncombos)
resBootstrap=data.frame(n=rep(0,ncombos),o=rep(0,ncombos),l=rep(0,ncombos),u=rep(0,ncombos),p=rep(0,ncombos))
n=ncol(mFst)
sn=0.9*n

for (i in 1:nperm){
	rarray= sort(sample(n,sn,replace=F))
	mmFst = mFst[rarray,rarray][lower.tri(mFst[rarray,rarray],diag=F)]
	mmDIST= prep(mDist[rarray,rarray])
	mmELE = prep(mAlt[rarray,rarray])
	mmENV = prep(mEnv[rarray,rarray])
	mmMASS =prep(mMass[rarray,rarray])
	
	comm=regr(lm(mmFst ~ mmDIST + mmMASS + mmELE + mmENV))
	boot[i,]=comm$Commonality_Data$CC[c(1:ncombos),1]
	print(i)
}
for (i in 1:ncombos){
	q=quantile(boot[,i], c(.025,.975))
	resBootstrap[i,1]=i
	resBootstrap[i,3]=q[1]
	resBootstrap[i,4]=q[2]
}

resBootstrap[,2]=ca$Commonality_Data$CC[c(1:ncombos),1]
resBootstrap[,5]=ca$Commonality_Data$CC[c(1:ncombos),2]
rownames(resBootstrap) = rownames(ca$Commonality_Data$CC)[1:ncombos]
resBootstrap[,c('o','l','u')] = resBootstrap[,c('o','l','u')]

new.rownames = gsub('\\s|[?!Unique$|Common$|to$|and$]','',rownames(resBootstrap))
rownames(resBootstrap) = new.rownames

pdf('outfiles/Figure S8.pdf',width=8,height=7, bg='transparent')

tmp = resBootstrap
#dev.new(width=8,height=7)
tmp = round(tmp,3)
tmp$n = rev(tmp$n)

par(mar=c(4,2,4,20))
xlim = c(min(tmp[,3]), max(tmp[,4]))
plot(tmp[,2],tmp[,1],xlim=xlim,font=5,lab=c(ncombos, 7, 1),xaxt="n",yaxt="n",cex.lab=1,xlab='',ylab='',cex=1.5) 
arrows(tmp[,3],tmp[,1],tmp[,4],tmp[,1],code=3,length=0.05,angle=90,lwd=1)
abline(v=0,lty=2)
axis(1,cex.axis=1)
mtext('Correlation Coefficient',1,line=2.5,cex=1.2)
#yaxis
x=cbind(1:15,rev(rownames(tmp)))
axis(4,at=x[,1],labels=x[,2],cex.axis=1,las=2)
mtext('Predictor Sets',2,cex=1.2,line=.5)
#upper xaxis
seq(0,round(max(tmp[,4]),2),by=round(max(tmp[,4]),2)/5)
at = seq(0,round(max(tmp[,4]),2),by=round(max(tmp[,4]),2)/5)
labels= round(at/sum(tmp[,'o']),2)
axis(3,at=at,labels=labels,cex.axis=1)
mtext('% Total',3,line=2.5,cex=1.2)
#coefficients
mtext('Coefficient',4,at = 16,las=2,line=12,adj=.5)
axis(4,line=12,at=tmp[,1],labels=tmp[,'o'],tick=F,cex.axis=1,las=2,hadj=1)
mtext('% Total',4,at = 16,las=2,line=17,adj=.5)
axis(4,line=17,at=tmp[,1],labels=tmp[,'p'],tick=F,cex.axis=1,las=2,hadj=1)

mtext('Total',4,at = -.5,las=2,line=8,adj=.5)
mtext(sum(tmp[,'o']),4,at = -.5,las=2,line=12,adj=.5)
mtext(100,4,at = -.5,las=2,line=17,adj=.5)

dev.off()










