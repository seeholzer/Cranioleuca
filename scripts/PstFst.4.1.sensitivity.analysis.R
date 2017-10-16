#=============================================================================================================#
# Script created by Glenn Seeholzer, contact at seeholzer.glenn@gmail.com
# Script created in version R 3.3.2  
# This script: Generates Fst and Pst tables for Pst-Fst sensitivity analysis described in Discussion section "Pst-Fst"
# Usage notes: 	Step 1 of 2. 
# 				Best to run this terminal as the calculations take a while. I've written the script to run
#				in parallel on Mac machines using a foreach loop from the package foreach.
#				Simply enter the following command to make the script run
#				Rscript {insert path to DataArchival directory here}/scripts/PstFst.4.1.sensitivity.analysis.R
#				Once script has finished run the second R script: PstFst.4.2.sensitivity.analyses.R to produce figures
#=============================================================================================================#
#set working directory
setwd('~/Dropbox/CRAN/Chapter.2/Seeholzer&Brumfield_ms_revision1/DataArchival/')

#create directory for sensitivity analysis data files
dir.create('outfiles/sensitivity.analysis.outfiles/', showWarnings = FALSE)

#load packages and supporting functions
library(adegenet)
library(hierfstat)
library(MCMCglmm)
library(caTools)

####	import phenotypic and meta-data
data = read.delim('Data.Table.txt',stringsAsFactors=F)
#subset to get rid of un-sexed individuals, individuals without genetic data, non-adults, and non-antisiensis
data = data[data$Sex %in% c('m','f') & !(data$Age %in% 'j') & !(data$Population %in% c('erythrops','curtata')), ]


gen = get(load('SNPdata/SNPdata.genind.Rdata'))
#specify the populations
tmp = data[data$ID %in% indNames(gen),c('ID','Population')]
tmp = tmp[match(indNames(gen),tmp$ID), ]
strata(gen) = data.frame(pop = tmp[,'Population'])
setPop(gen) = ~pop

#gen = gen[loc=sample(1:5154,50)] #subsample so runs quicker

outliers = readLines('BayEnv_outlier_loci.txt')


set.seed(1234)
results.list = list()
traits = c('Mass','Wing.length','PC1.plumage')


##################
#	Setup parallel processing
library(foreach)
library(doMC)
registerDoMC(8)  #change to the number of parallel processes you want running  
##################

i = 1
results.list = foreach(i=1:102) %dopar% {
#for(i in 1:16){
	if(i == 1){
		pops.excluded = c()		
	}else if(i == 2){
		pops.excluded = c('Macate','Pacar','Quichas','Huamatanga')
	}else{
		nonSWslope.populations = unique(data$Population)[!(unique(data$Population) %in% c('Macate','Pacar','Quichas','Huamatanga'))]
		all.combos = combs(nonSWslope.populations,4)
		pops.excluded = all.combos[sample(1:nrow(all.combos),1), ]
	}

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
	t = 1
	for(t in 1:length(traits)){
	
		cat(i,' ',traits[t],'\n')

	################################################################################################
	####	Fst
		IDs = data[!is.na(data[,traits[t]]) & !(data$Population %in% pops.excluded),'ID'] 
		gen.new = gen[indNames(gen)[which(indNames(gen) %in% IDs)]] #subset to get same individuals as in data
		#subset to get only pops with 4 or more individuals
		tab = table(gen.new@pop)
		pops.to.remove = names(tab[tab <= 3 ])
		gen.new = gen.new[pop=popNames(gen.new)[!(popNames(gen.new) %in% pops.to.remove)]]
		#remove outlier loci
		gen.new = gen.new[loc=locNames(gen.new)[!(locNames(gen.new) %in% outliers)]]
		#get Fst stats
		gen2 = genind2hierfstat(gen.new)
		wc.all = wc(gen2)
		wc.per.loc = wc.all$per.loc$FST ; names(wc.per.loc) = locNames(gen.new)
		wc.global = wc.all$FST
		IDs = indNames(gen.new)
		fst.results = list(wc.per.loc=wc.per.loc,wc.global=wc.global,IDs=IDs)

	################################################################################################
	####	Pst
		
		####	get data
	
		if(traits[t] %in% 'Mass'){
			tmp.data = data[,c('ID','Population','Sex','Month',traits[t])]
			tmp.data = tmp.data[complete.cases(tmp.data),]
		}else if(traits[t] %in% c('Wing.length','PC1.plumage')){
			tmp.data = data[,c('ID','Population','Sex',traits[t])]
			tmp.data = tmp.data[complete.cases(tmp.data),]
		}
		
		####	subset to just the individuals with genetic data
			tmp.data = tmp.data[tmp.data$ID %in% fst.results$IDs,]	
			
		####	get sample size table for summary table
			sampsize = table(tmp.data$Population)
							
		####	Create formula and prior vector for MCMCglmm	
			set.seed(545)
			if(traits[t] %in% 'Mass'){
					formula = formula(paste0(traits[t],'~Sex+Month'))
					prior = list(R = list(R1 = list(V = 1, n = 0),R2 = list(V = 1, n = 0), R3 = list(V = 1, n = 0)),G = list(G1 = list(V = 1, n = 1)))	
				}else if(traits[t] %in% c('Wing.length','PC1.plumage')){
					formula = formula(paste0(traits[t],'~Sex'))
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
			if(traits[t] == 'Mass'){
			##########################################		
				
				Fst.global = fst.results$wc.global
				Fst.all = fst.results$wc.per.loc
				
				Fst[Fst$trait %in% 'Mass',-1] = c(Fst.global,quantile(Fst.all, c(.05, .95),na.rm=T))
				#make formula
				critical = function(Fst,varB,varW) (-2*Fst*varW)/(varB*(Fst-1))
				
				upper.Fst = Fst[Fst$trait == 'Mass','upper']
				CV.CHratio = critical(upper.Fst,param['varB','lower'],param['varW','upper'])
					
			##########################################		
					}else if(traits[t] %in% c('Wing.length')){ #morph variables
			##########################################		
		
				Fst.global = fst.results$wc.global
				Fst.all = fst.results$wc.per.loc
				
				Fst[Fst$trait %in% 'Wing.length',-1] = c(Fst.global,quantile(Fst.all, c(.05, .95),na.rm=T))
				#make formula
				critical = function(Fst,varB,varW) (-2*Fst*varW)/(varB*(Fst-1))
				
				upper.Fst = Fst[Fst$trait == 'Wing.length','upper']
				CV.CHratio = critical(upper.Fst,param['varB','lower'],param['varW','upper'])
				
			##########################################		
					}else if(traits[t] == 'PC1.plumage'){
			##########################################	
						
				Fst.global = fst.results$wc.global
				Fst.all = fst.results$wc.per.loc
				
				Fst[Fst$trait %in% 'PC1.plumage',-1] = c(Fst.global,quantile(Fst.all, c(.05, .95),na.rm=T))
				#make formula
				critical = function(Fst,varB,varW) (-2*Fst*varW)/(varB*(Fst-1))
				
				upper.Fst = Fst[Fst$trait == 'PC1.plumage','upper']
				CV.CHratio = critical(upper.Fst,param['varB','lower'],param['varW','upper'])
		
					}
				
			######################################################
			##################	add Pst data to table
			######################################################
			Pst.table[t,'trait'] = traits[t]
			Pst.table[t,'sample.size'] = paste0(length(sampsize),' (',round(mean(sampsize),1),', ',min(sampsize),'-',max(sampsize),')')
			Pst.table[t,'varB'] = param['varB','mode']
			Pst.table[t,'lower.varB'] = param['varB','lower']
			Pst.table[t,'upper.varB'] = param['varB','upper']
			Pst.table[t,'varW'] = param['varW','mode']
			Pst.table[t,'lower.varW'] = param['varW','lower']
			Pst.table[t,'upper.varW'] = param['varW','upper']
			Pst.table[t,'Pst'] = Pst.mode(1)
			Pst.table[t,'lower CI'] = Pst.lower(1)
			Pst.table[t,'upper CI'] = Pst.upper(1)
			Pst.table[t,'CV.CHratio'] = CV.CHratio	
		
	} # t (trait) loop

	results = list(pops.excluded = pops.excluded, Fst.table = Fst, Pst.table = Pst.table)
	results
	save(results,file=paste0('outfiles/sensitivity.analysis.outfiles/sensitivity.analysis.',i,'.Rdata'))

}	# random sample i loop










