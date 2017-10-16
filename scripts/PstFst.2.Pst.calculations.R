#=============================================================================================================#
# Script created by Glenn Seeholzer, contact at seeholzer.glenn@gmail.com
# Script created in version R 3.3.2 
# This script: Generates Fst data for Pst-Fst
# Usage notes: 	Step 2 of 3 in Pst-Fst analysis pipeline
#				Entire script is automated in a for loop and can be run all at once. BUT, must specify which
#				dataset to run it on: All Populations or SW Slope Excluded
#				This is accomplished by commenting out the two lines pertaining to each dataset in the
#				"Alternative Dataset" section and leaving the other uncommented.
#				Then, in R console or Rstudio, just highlight the entire document and send to console. 
#=============================================================================================================#
#set working directory
setwd('~/Dropbox/CRAN/Chapter.2/Seeholzer&Brumfield_ms_revision1/DataArchival/')

#load packages and supporting functions
library(MCMCglmm)

# Key to some terms 
#random effect = Population
#fixed effect = Sex
#response = Mass
#G is for the random effects, must put Population as random effect so you can estimate it's variance component
#R is for the residual variance

####	import phenotypic and meta-data
data = read.delim('Data.Table.txt',stringsAsFactors=F)
####	subset to get rid of un-sexed individuals, individuals without geSWslopeExcludedtic data, non-adults, and non-antisiensis
data.all = data[data$Sex %in% c('m','f') & !(data$Age %in% 'j') & !(data$Population %in% c('erythrops','curtata')), ]


#----------------------------------------------
#----------------------------------------------
#------------ Alternative Datasets ------------
#----------------------------------------------
#----------------------------------------------
#----------------------------------------------
#----------------------------------------------
#----------------------------------------------
#----------------------------------------------
#----------------------------------------------

######	all populations

data = data.all
ALL = TRUE; SWslopeExcluded= FALSE

######	SW Slope Excluded

# data = data.all[!(data.all$Population %in% c('Macate','Pacar','Quichas','Huamatanga')),]
# ALL = FALSE; SWslopeExcluded= TRUE


#----------------------------------------------
#----------------------------------------------
#----------------------------------------------
#----------------------------------------------
#----------------------------------------------
#----------------------------------------------
#----------------------------------------------
#----------------------------------------------
#----------------------------------------------
#----------------------------------------------
#----------------------------------------------
#----------------------------------------------
#----------------------------------------------
#----------------------------------------------
#----------------------------------------------
#----------------------------------------------
#----------------------------------------------
#----------------------------------------------




####	Make empty Fst table
colnames = c('trait','global','lower','upper')
traits = c('Mass','Wing.length','PC1.plumage')
Fst = data.frame(matrix(nrow=length(traits),ncol=length(colnames)))
colnames(Fst) = colnames
Fst$trait = traits

####	Make empty Pst table
traits = c('Mass','Wing.length','PC1.plumage')
colnames = c('trait','sample.size','varB','lower.varB','upper.varB','varW','lower.varW','upper.varW','Pst','lower CI','upper CI','CV.CHratio')
Pst.table = data.frame(matrix(nrow=length(traits),ncol=length(colnames)))
colnames(Pst.table) = colnames

i = 1
for(i in 1:length(traits)){
	
	#get the Fst values and the individuals used with genetic data
	
	if(traits[i] == 'Mass'){	# body mass
		if(ALL){
			tmp = get(load('outfiles/Fst.Mass.all.Rdata'))	
			}
		if(SWslopeExcluded){
			tmp = get(load('outfiles/Fst.Mass.SWslopeExcluded.Rdata'))	
		}
	
	}else if(traits[i] %in% c('Wing.length')){ #morph variable
		if(ALL){
			tmp = get(load('outfiles/Fst.Wing.length.all.Rdata'))	
			}
		if(SWslopeExcluded){
			tmp = get(load('outfiles/Fst.Wing.length.SWslopeExcluded.Rdata'))	
		}
	
	}else if(traits[i] == 'PC1.plumage'){ #plumage variable
		if(ALL){
			tmp = get(load('outfiles/Fst.PC1.plumage.all.Rdata'))	
			}
		if(SWslopeExcluded){
			tmp = get(load('outfiles/Fst.PC1.plumage.SWslopeExcluded.Rdata'))	
		}
	}
	
	if(traits[i] %in% 'Mass'){
		tmp.data = data[,c('ID','Population','Sex','Month',traits[i])]
		tmp.data = tmp.data[complete.cases(tmp.data),]
	}else if(traits[i] %in% c('Wing.length','PC1.plumage')){
		tmp.data = data[,c('ID','Population','Sex',traits[i])]
		tmp.data = tmp.data[complete.cases(tmp.data),]
	}
	
####	subset to just the individuals with genetic data
	tmp.data = tmp.data[tmp.data$ID %in% tmp$IDs,]	
	
####	get sample size table for summary table
	sampsize = table(tmp.data$Population)
					
	cat('\n\n',traits[i],'\n\n')
####	Create formula and prior vector for MCMCglmm	
	set.seed(545)
	if(traits[i] %in% 'Mass'){
			formula = formula(paste0(traits[i],'~Sex+Month'))
			prior = list(R = list(R1 = list(V = 1, n = 0),R2 = list(V = 1, n = 0), R3 = list(V = 1, n = 0)),G = list(G1 = list(V = 1, n = 1)))	
		}else if(traits[i] %in% c('Wing.length','PC1.plumage')){
			formula = formula(paste0(traits[i],'~Sex'))
			prior = list(R = list(V = 1, n = 0), G = list(G1 = list(V = 1, n = 1)))
		}

####	run Bayesian generalized linear mixed model		
	fit = MCMCglmm(formula, random = ~Population, data = tmp.data, family = "gaussian", prior = prior, nitt=100000, burnin=20000, verbose = F)
	sum = summary(fit)
####	get vector of parameter values			
	param = cbind(mode=posterior.mode(fit$VCV),HPDinterval(fit$VCV, 0.95))
	param = rbind(param,1)
	rownames(param) = c('varB','varW','CHratio')
####	Key for object names	
	#varB = variance between populations
	#varW = variance within populations
	#CHratio = ratio of the proportion of variation due to additive genetic effects across populations (c) to the proportion within populations (h2, i.e. heretability)  
	#CV.CHratio = critcal value of CHratio at which Pst no longer exceeds the upper confidence interval of Fst
	Pst = function(varB,varW,CHratio) (CHratio*varB)/((CHratio*varB)+(2*varW))

####	the lower estimate of varB must be paired with the upper 
####	estimate of varW to get the lower confidence limit for Pst
####	vice versa for the upper limit
	
	Pst.mode = function(CHratio) (CHratio*param['varB','mode']) / ((CHratio*param['varB','mode'])+(2*param['varW','mode']))
	Pst.upper = function(CHratio) (CHratio*param['varB','upper'])/((CHratio*param['varB','upper'])+(2*param['varW','lower']))
	Pst.lower = function(CHratio) (CHratio*param['varB','lower'])/((CHratio*param['varB','lower'])+(2*param['varW','upper']))
	
	######################################################
	##################	Fst
	######################################################

	##########################################		
	if(traits[i] == 'Mass'){
	##########################################		
		
		Fst.global = tmp$wc.global
		Fst.all = tmp$wc.per.loc
		
		Fst[Fst$trait %in% 'Mass',-1] = c(Fst.global,quantile(Fst.all, c(.05, .95),na.rm=T))
		#make formula
		critical = function(Fst,varB,varW) (-2*Fst*varW)/(varB*(Fst-1))
		
		upper.Fst = Fst[Fst$trait == 'Mass','upper']
		CV.CHratio = critical(upper.Fst,param['varB','lower'],param['varW','upper'])
			
	##########################################		
			}else if(traits[i] == 'Wing.length'){ #morph variables
	##########################################		

		Fst.global = tmp$wc.global
		Fst.all = tmp$wc.per.loc
		
		Fst[Fst$trait %in% 'Wing.length',-1] = c(Fst.global,quantile(Fst.all, c(.05, .95),na.rm=T))
		#make formula
		critical = function(Fst,varB,varW) (-2*Fst*varW)/(varB*(Fst-1))
		
		upper.Fst = Fst[Fst$trait %in% 'Wing.length','upper']
		CV.CHratio = critical(upper.Fst,param['varB','lower'],param['varW','upper'])

	##########################################		
			}else if(traits[i] == 'PC1.plumage'){
	##########################################	
				
		Fst.global = tmp$wc.global
		Fst.all = tmp$wc.per.loc
		
		Fst[Fst$trait %in% 'PC1.plumage',-1] = c(Fst.global,quantile(Fst.all, c(.05, .95),na.rm=T))
		#make formula
		critical = function(Fst,varB,varW) (-2*Fst*varW)/(varB*(Fst-1))
		
		upper.Fst = Fst[Fst$trait == 'PC1.plumage','upper']
		CV.CHratio = critical(upper.Fst,param['varB','lower'],param['varW','upper'])

			}
		
	######################################################
	##################	add Pst data to table
	######################################################
	Pst.table[i,'trait'] = traits[i]
	Pst.table[i,'sample.size'] = paste0(length(sampsize),' (',round(mean(sampsize),1),', ',min(sampsize),'-',max(sampsize),')')
	Pst.table[i,'varB'] = param['varB','mode']
	Pst.table[i,'lower.varB'] = param['varB','lower']
	Pst.table[i,'upper.varB'] = param['varB','upper']
	Pst.table[i,'varW'] = param['varW','mode']
	Pst.table[i,'lower.varW'] = param['varW','lower']
	Pst.table[i,'upper.varW'] = param['varW','upper']
	Pst.table[i,'Pst'] = Pst.mode(1)
	Pst.table[i,'lower CI'] = Pst.lower(1)
	Pst.table[i,'upper CI'] = Pst.upper(1)
	Pst.table[i,'CV.CHratio'] = CV.CHratio	
	print(Pst.table)
	cat('\n\n')
}


		if(ALL){
			write.table(Fst,'outfiles/Fst.table.all.txt',col.names=T,row.names=F,quote=F,sep='\t')
			write.table(Pst.table,'outfiles/Pst.table.all.txt',col.names=T,row.names=F,quote=F,sep='\t')
			}
		if(SWslopeExcluded){
			write.table(Fst,'outfiles/Fst.table.SWslopeExcluded.txt',col.names=T,row.names=F,quote=F,sep='\t')
			write.table(Pst.table,'outfiles/Pst.table.SWslopeExcluded.txt',col.names=T,row.names=F,quote=F,sep='\t')
		}
		
		

		
		