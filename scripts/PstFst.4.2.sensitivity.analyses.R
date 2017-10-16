#=============================================================================================================#
# Script created by Glenn Seeholzer, contact at seeholzer.glenn@gmail.com
# Script created in version R 3.3.2  
# This script: Generates Figure S6
# Usage notes: 	Step 2 of 2. 
# 				Run line-by-line
#=============================================================================================================#
#set working directory
setwd('~/Dropbox/CRAN/Chapter.2/Seeholzer&Brumfield_ms_revision1/DataArchival/')

~/Dropbox/CRAN/Chapter.2/Seeholzer&Brumfield_ms_revision1/DataArchival/outfiles/

files = list.files('outfiles/sensitivity.analysis.outfiles',full.names=T)
results.list = list()
for(i in 1:length(files)){
	cat(i,'\n')
	results.list[[i]] = get(load(files[grep(paste0('sensitivity.analysis.',i,'.Rdata'),files)]))
}

mass.ratios 			= sapply(results.list,function(x) x[[3]][x[[3]] == 'Mass','CV.CHratio'] )
Wing.length.ratios 	= sapply(results.list,function(x) x[[3]][x[[3]] == 'Wing.length','CV.CHratio'] )
PC1.plumage.ratios	 	= sapply(results.list,function(x) x[[3]][x[[3]] == 'PC1.plumage','CV.CHratio'] )



pdf('outfiles/Figure S6.pdf',width=12,height=4.5,bg='transparent')

#dev.new(width=12,height=4.5)
mat = matrix(c(1:3), nrow=1, ncol=3, byrow=T)
mat = cbind(5,rbind(mat,4))
mat[2,1] = 6
nf = layout(mat=mat,heights=c(4,.25),widths=c(.5,4,4,4))
#layout.show(nf)
	
par(mar=c(4,1,3,2),xpd=FALSE)
	
# Mass
ratios = mass.ratios
breaks = seq(0,max(ratios)+0.2,0.05)
hist(ratios[c(-1,-2)],breaks=breaks,xlab='',main='',axes=F)
axis(1)
axis(2)
abline(v= ratios[1],lty=1,lwd=3)
abline(v= ratios[2],lty=2,lwd=3)
mtext('a) Body Mass',side=3,line=1,cex=1.25,adj=0)


# PC1.morphology
ratios = Wing.length.ratios
breaks = seq(0,max(ratios)+0.2,0.05)
hist(ratios[c(-1,-2)],breaks=breaks,xlab='',main='',axes=F)
axis(1)
abline(v= ratios[1],lty=1,lwd=3)
abline(v= ratios[2],lty=2,lwd=3)
mtext('b) Wing Length',side=3,line=1,cex=1.25,adj=0)


#  PC1.plumage
ratios = PC1.plumage.ratios
breaks = seq(0,max(ratios)+0.2,0.05)
hist(ratios[c(-1,-2)],breaks=breaks,xlab='',main='',axes=F)
axis(1)
abline(v= ratios[1],lty=1,lwd=3)
abline(v= ratios[2],lty=2,lwd=3)
mtext('c) Plumage',side=3,line=1,cex=1.25,adj=0)

par(mar=c(0,0,0,0))
plot(1,1,type='n',axes=FALSE, ann=FALSE)
mtext(bquote(italic(c)/italic(h)^2),side=1,line=-1,cex=2)

par(mar=c(0,0,0,0))
plot(1,1,type='n',axes=FALSE, ann=FALSE)
mtext("Frequency",side=2,line=-2.5,cex=2)

dev.off()



