#=============================================================================================================#
# Script created by Glenn Seeholzer, contact at seeholzer.glenn@gmail.com
# Script created in version R 3.3.2  
# This script: Generates Fst data for Pst-Fst
# Usage notes: 	Step 3 of 3 in Pst-Fst analysis pipeline
#				Entire script is automated and can be run all at once. In R console or Rstudio, just 
#				highlight the entire document and send to console. 
#				Produces Figure 4 (upper and lower panels) and Table 2.
#=============================================================================================================#
#set working directory
setwd('~/Dropbox/CRAN/Chapter.2/Seeholzer&Brumfield_ms_revision1/DataArchival/')

#load supporting functions 
source('~/Dropbox/FUN.add.alpha.R', chdir = TRUE)


#Load datasets
######	All Populations
Fst.sum.all = read.delim('outfiles/Fst.table.all.txt',stringsAsFactors=F)
Pst.all = read.delim('outfiles/Pst.table.all.txt',stringsAsFactors=F)

######	SW Slope Excluded
Fst.sum.SWslopeExcluded = read.delim('outfiles/Fst.table.SWslopeExcluded.txt',stringsAsFactors=F)
Pst.SWslopeExcluded = read.delim('outfiles/Pst.table.SWslopeExcluded.txt',stringsAsFactors=F)


####################
####################
####################
####################
####################
#	Table 2
####################
####################		
####################		
####################		
####################		

#all populations
sum = Pst.all[,c('trait','sample.size','varB','lower.varB','upper.varB','varW','lower.varW','upper.varW','Pst','lower.CI','upper.CI','CV.CHratio')]
sum[,-c(1,2)] = round(sum[,-c(1,2)],2)
sum = apply(sum,2,as.character)
sum = apply(sum,2,trimws)

traits = c('body mass','wing length','plumage')
samplesize = sum[,'sample.size']
varB = paste0(sum[,'varB'],' (',sum[,'lower.varB'],', ',sum[,'upper.varB'],')')
varW = paste0(sum[,'varW'],' (',sum[,'lower.varW'],', ',sum[,'upper.varW'],')')
pst = paste0(sum[,'Pst'],' (',sum[,'lower.CI'],', ',sum[,'upper.CI'],')')
ch2 = sum[,'CV.CHratio']

tmp = Fst.sum.all
colnames(tmp) = c('dataset','Fst.global','Fst.lower.CI','Fst.upper.CI')
tmp = round(tmp[,-1],3)
tmp = apply(tmp,1:2, function(x) sprintf("%.2f", round(x,2)))
fst = paste0(tmp[,'Fst.global'],' (',tmp[,'Fst.lower.CI'],', ',tmp[,'Fst.upper.CI'],')')

sum.all = cbind(traits,samplesize,varB,varW,pst,fst,ch2)


#SW slope excluded
sum = Pst.SWslopeExcluded[,c('trait','sample.size','varB','lower.varB','upper.varB','varW','lower.varW','upper.varW','Pst','lower.CI','upper.CI','CV.CHratio')]
sum[,-c(1,2)] = round(sum[,-c(1,2)],2)
sum = apply(sum,2,as.character)
sum = apply(sum,2,trimws)

traits = c('body mass','wing length','plumage')
samplesize = sum[,'sample.size']
varB = paste0(sum[,'varB'],' (',sum[,'lower.varB'],', ',sum[,'upper.varB'],')')
varW = paste0(sum[,'varW'],' (',sum[,'lower.varW'],', ',sum[,'upper.varW'],')')
pst = paste0(sum[,'Pst'],' (',sum[,'lower.CI'],', ',sum[,'upper.CI'],')')
ch2 = sum[,'CV.CHratio']

tmp = Fst.sum.SWslopeExcluded
colnames(tmp) = c('dataset','Fst.global','Fst.lower.CI','Fst.upper.CI')
tmp = round(tmp[,-1],3)
tmp = apply(tmp,1:2, function(x) sprintf("%.2f", round(x,2)))
fst = paste0(tmp[,'Fst.global'],' (',tmp[,'Fst.lower.CI'],', ',tmp[,'Fst.upper.CI'],')')

sum.SWslopeExcluded = cbind(traits,samplesize,varB,varW,pst,fst,ch2)

sum = rbind(sum.all,'',sum.SWslopeExcluded)
sum = cbind('',sum)
sum[,1] = c('All Populations','','','','SW Slope Excluded','','')


write.table(sum,'outfiles/Table 2.txt',quote=F,col.names=T,row.names=F,sep='\t')




####################
####################
####################
####################
####################
#	Figure 4, Pst-Fst Plot
#	The first section makes the upper panel (all populations) of Figure 4,
#	and the second section make the lower panel (SW slope populations excluded) of Figure 4
#	These panels were combined in Adobe Illustrator for the final version of Figure 4.
####################
####################		
####################		
####################		
####################

#----------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------
#	Figure 4 Upper Panel : All populations
#----------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------

Fst.sum = Fst.sum.all
Pst = Pst.all
ALL = TRUE; SWslopeExcluded= FALSE

pdf('outfiles/Figure 4 Upper Panel - All Populations.pdf',width=12,height=4.5,bg='transparent')
	
	key = c('Mass','Wing Length','Plumage')
	names(key) = c('Mass','Wing.length','PC1.plumage')
	
	#dev.new(width=12,height=4.5)
	mat = matrix(c(1:3), nrow=1, ncol=3, byrow=T)
	mat = cbind(5,rbind(mat,4))
	mat[2,1] = 6
	nf = layout(mat=mat,heights=c(4,.5),widths=c(1,4,4,4))
	#layout.show(nf)
	
	i = 1
	for(i in 1:3){
		
		cat(Pst[i,'trait'],'\n')
		
		if(Pst[i,'trait'] == 'Mass'){
			
			if(ALL){
				tmp = get(load('outfiles/Fst.Mass.all.Rdata'))	
				}
			if(SWslopeExcluded){
				tmp = get(load('outfiles/Fst.Mass.SWslopeExcluded.Rdata'))	
			}
			
			Fst.global = tmp$wc.global
			Fst.all = tmp$wc.per.loc
			Fst.global.CI = Fst.sum[Fst.sum$trait == 'Mass',]
	
		}else if(Pst[i,'trait'] == 'Wing.length'){
		
			if(ALL){
				tmp = get(load('outfiles/Fst.Wing.length.all.Rdata'))	
				}
			if(SWslopeExcluded){
				tmp = get(load('outfiles/Fst.Wing.length.SWslopeExcluded.Rdata'))	
			}
						
			Fst.global = tmp$wc.global
			Fst.all = tmp$wc.per.loc
			Fst.global.CI = Fst.sum[Fst.sum$trait == 'Wing.length',]	
				
		}else if(Pst[i,'trait'] == 'PC1.plumage'){
		
			if(ALL){
				tmp = get(load('outfiles/Fst.PC1.plumage.all.Rdata'))	
				}
			if(SWslopeExcluded){
				tmp = get(load('outfiles/Fst.PC1.plumage.SWslopeExcluded.Rdata'))	
			}
						
			Fst.global = tmp$wc.global
			Fst.all = tmp$wc.per.loc
			Fst.global.CI = Fst.sum[Fst.sum$trait == 'PC1.plumage',]	
				
		}
		
		
		######################################################
		##################	plotting
		######################################################
		#curve functions
		Pst.mode = function(CHratio) (CHratio*Pst[i,'varB']) / ((CHratio*Pst[i,'varB'])+(2*Pst[i,'varW']))
		Pst.upper = function(CHratio) (CHratio*Pst[i,'upper.varB']) / ((CHratio*Pst[i,'upper.varB'])+(2*Pst[i,'lower.varW']))
		Pst.lower = function(CHratio) (CHratio*Pst[i,'lower.varB']) / ((CHratio*Pst[i,'lower.varB'])+(2*Pst[i,'upper.varW']))
			
	
		#margins for the combined plot
		par(mar=c(4,1,3,2),xpd=FALSE)
		
		#specify upper limit of x-axis
		to = 2 
		#set up blank plot
		curve(Pst.mode,from = 0, to=to,lwd=2,add=F,xlim=c(0,to),ylim=c(0,1),type='n',axes=F,xlab='',ylab='')
		
		
		#add axes
		axis(1,lwd=2) 
		axis(2,lwd=2)
		text = key[Pst[i,'trait']]
		
		mtext(paste0(letters[i],'. ',text),side=3,line=.5,cex=1.25,adj=0)
		
		#add curves
		upper.curve = curve(Pst.upper,from = 0, to=to,lwd=1,lty=2,add=T,col='blue')
		lower.curve = curve(Pst.lower,from = 0, to=to,lwd=1,lty=2,add=T,col='blue')
		polygon(c(upper.curve$x,rev(lower.curve$x)),c(upper.curve$y,rev(lower.curve$y)),border=F,col=add.alpha('skyblue',.5))
		mode.curve = curve(Pst.mode,from = 0, to=to,lwd=2,add=T,col=add.alpha('blue',.95))
		
		#add vertical histogram of distribution of Fst across genome
		dens = hist(Fst.all,plot=F)$density
		dens = dens/(4*max(dens))
		barplot.col = add.alpha('red',.75)
		barplot(dens, width=.05, space=0, border=F, col=barplot.col, horiz=T,axes=F, add=T,offset=Pst[i,'CV.CHratio']) # barplot
		barplot(-dens, width=.05, space=0, border=F, col=barplot.col, horiz=T,axes=F, add=T,offset=Pst[i,'CV.CHratio']) # barplot
		
		abline(h= Fst.global.CI[,'upper'],lty=5,col=add.alpha('red',.95))
		abline(h= Fst.global.CI[,'global'],lty=1,lwd=2,col=add.alpha('red',.95))
		abline(h= Fst.global.CI[,'lower'],lty=5,col=add.alpha('red',.95))
	
		abline(v=Pst[i,'CV.CHratio'],lty=6,lwd=2,col='black')
	
		
	}	
		par(mar=c(0,0,0,0))
		plot(1,1,type='n',axes=FALSE, ann=FALSE)
		mtext(bquote(italic(c)/italic(h)^2),side=1,line=-3,cex=2)
		mtext('less heritable',side=1,line=-4,cex=1.2,adj = 0)
		mtext('more heritable',side=1,line=-4,cex=1.2,adj = 1)
	
		par(mar=c(0,0,0,0))
		plot(1,type='n',axes=FALSE, ann=FALSE)
		mtext(bquote(italic(P)[ST]),side=2,line=-5,cex=2,col=add.alpha('blue',.75),adj=.35)
		mtext('or',side=2,line=-5,cex=1.5,col='black')
		mtext(bquote(italic(F)[ST]),side=2,line=-5,cex=2,col=add.alpha('red',.75),adj=.65)		
	
dev.off()


#----------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------
#	Figure 4 Lower Panel : SW Slope Excluded
#----------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------

Fst.sum = Fst.sum.SWslopeExcluded
Pst = Pst.SWslopeExcluded
ALL = F; SWslopeExcluded= T

pdf('outfiles/Figure 4 Lower Panel - SW slope excluded.pdf',width=12,height=4.5,bg='transparent')
	
	key = c('Mass','Wing Length','Plumage')
	names(key) = c('Mass','Wing.length','PC1.plumage')
	
	#dev.new(width=12,height=4.5)
	mat = matrix(c(1:3), nrow=1, ncol=3, byrow=T)
	mat = cbind(5,rbind(mat,4))
	mat[2,1] = 6
	nf = layout(mat=mat,heights=c(4,.5),widths=c(1,4,4,4))
	#layout.show(nf)
	
	i = 1
	for(i in 1:3){
		
		cat(Pst[i,'trait'],'\n')
		
		if(Pst[i,'trait'] == 'Mass'){
			
			if(ALL){
				tmp = get(load('outfiles/Fst.Mass.all.Rdata'))	
				}
			if(SWslopeExcluded){
				tmp = get(load('outfiles/Fst.Mass.SWslopeExcluded.Rdata'))	
			}
			
			Fst.global = tmp$wc.global
			Fst.all = tmp$wc.per.loc
			Fst.global.CI = Fst.sum[Fst.sum$trait == 'Mass',]
	
		}else if(Pst[i,'trait'] == 'Wing.length'){
		
			if(ALL){
				tmp = get(load('outfiles/Fst.Wing.length.all.Rdata'))	
				}
			if(SWslopeExcluded){
				tmp = get(load('outfiles/Fst.Wing.length.SWslopeExcluded.Rdata'))	
			}
						
			Fst.global = tmp$wc.global
			Fst.all = tmp$wc.per.loc
			Fst.global.CI = Fst.sum[Fst.sum$trait == 'Wing.length',]	
				
		}else if(Pst[i,'trait'] == 'PC1.plumage'){
		
			if(ALL){
				tmp = get(load('outfiles/Fst.PC1.plumage.all.Rdata'))	
				}
			if(SWslopeExcluded){
				tmp = get(load('outfiles/Fst.PC1.plumage.SWslopeExcluded.Rdata'))	
			}
						
			Fst.global = tmp$wc.global
			Fst.all = tmp$wc.per.loc
			Fst.global.CI = Fst.sum[Fst.sum$trait == 'PC1.plumage',]	
				
		}
		
		
		######################################################
		##################	plotting
		######################################################
		#curve functions
		Pst.mode = function(CHratio) (CHratio*Pst[i,'varB']) / ((CHratio*Pst[i,'varB'])+(2*Pst[i,'varW']))
		Pst.upper = function(CHratio) (CHratio*Pst[i,'upper.varB']) / ((CHratio*Pst[i,'upper.varB'])+(2*Pst[i,'lower.varW']))
		Pst.lower = function(CHratio) (CHratio*Pst[i,'lower.varB']) / ((CHratio*Pst[i,'lower.varB'])+(2*Pst[i,'upper.varW']))
			
	
		#margins for the combined plot
		par(mar=c(4,1,3,2),xpd=FALSE)
		
		#specify upper limit of x-axis
		to = 2 
		#set up blank plot
		curve(Pst.mode,from = 0, to=to,lwd=2,add=F,xlim=c(0,to),ylim=c(0,1),type='n',axes=F,xlab='',ylab='')
		
		
		#add axes
		axis(1,lwd=2) 
		axis(2,lwd=2)
		text = key[Pst[i,'trait']]
		
		mtext(paste0(letters[i],'. ',text),side=3,line=.5,cex=1.25,adj=0)
		
		#add curves
		upper.curve = curve(Pst.upper,from = 0, to=to,lwd=1,lty=2,add=T,col='blue')
		lower.curve = curve(Pst.lower,from = 0, to=to,lwd=1,lty=2,add=T,col='blue')
		polygon(c(upper.curve$x,rev(lower.curve$x)),c(upper.curve$y,rev(lower.curve$y)),border=F,col=add.alpha('skyblue',.5))
		mode.curve = curve(Pst.mode,from = 0, to=to,lwd=2,add=T,col=add.alpha('blue',.95))
		
		#add vertical histogram of distribution of Fst across genome
		dens = hist(Fst.all,plot=F)$density
		dens = dens/(4*max(dens))
		barplot.col = add.alpha('red',.75)
		barplot(dens, width=.05, space=0, border=F, col=barplot.col, horiz=T,axes=F, add=T,offset=Pst[i,'CV.CHratio']) # barplot
		barplot(-dens, width=.05, space=0, border=F, col=barplot.col, horiz=T,axes=F, add=T,offset=Pst[i,'CV.CHratio']) # barplot
		
		abline(h= Fst.global.CI[,'upper'],lty=5,col=add.alpha('red',.95))
		abline(h= Fst.global.CI[,'global'],lty=1,lwd=2,col=add.alpha('red',.95))
		abline(h= Fst.global.CI[,'lower'],lty=5,col=add.alpha('red',.95))
	
		abline(v=Pst[i,'CV.CHratio'],lty=6,lwd=2,col='black')
	
		
	}	
		par(mar=c(0,0,0,0))
		plot(1,1,type='n',axes=FALSE, ann=FALSE)
		mtext(bquote(italic(c)/italic(h)^2),side=1,line=-3,cex=2)
		mtext('less heritable',side=1,line=-4,cex=1.2,adj = 0)
		mtext('more heritable',side=1,line=-4,cex=1.2,adj = 1)
	
		par(mar=c(0,0,0,0))
		plot(1,type='n',axes=FALSE, ann=FALSE)
		mtext(bquote(italic(P)[ST]),side=2,line=-5,cex=2,col=add.alpha('blue',.75),adj=.35)
		mtext('or',side=2,line=-5,cex=1.5,col='black')
		mtext(bquote(italic(F)[ST]),side=2,line=-5,cex=2,col=add.alpha('red',.75),adj=.65)		
	
dev.off()
