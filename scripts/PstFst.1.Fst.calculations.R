#=============================================================================================================#
# Script created by Glenn Seeholzer, contact at seeholzer.glenn@gmail.com
# Script created in version R 3.3.2 
# This script: Generates Fst data for Pst-Fst
# Usage notes: 	Step 1 of 3 in Pst-Fst analysis pipeline 
# 				Best to run in terminal as the calculations take a while 
#=============================================================================================================#
#set working directory
setwd('~/Dropbox/CRAN/Chapter.2/Seeholzer&Brumfield_ms_revision1/DataArchival/')

#load packages and supporting functions
library(adegenet)
library(hierfstat)

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

################################################################################################
################################################################################################
##	Body mass
################################################################################################
################################################################################################
cat('\n\n\n\n Fst calculations: Body mass ...\n')

################################################################################################
####	Global Fst for all populations
	IDs = data[!is.na(data$Mass),'ID'] 
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
	results = list(wc.per.loc=wc.per.loc,wc.global=wc.global,IDs=IDs)
	save(results,file='outfiles/Fst.mass.all.Rdata')

################################################################################################
####	Global Fst excluding SW slope
	IDs = data[!is.na(data$Mass) & !(data$Population %in% c('Macate','Pacar','Quichas','Huamatanga')),'ID'] 
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
	results = list(wc.per.loc=wc.per.loc,wc.global=wc.global,IDs=IDs)
	save(results,file='outfiles/Fst.mass.SWslopeExcluded.Rdata')

################################################################################################
################################################################################################
##	Wing.length
################################################################################################
################################################################################################
cat('\n\n\n\n Fst calculations: Wing.length (size)...\n')

################################################################################################
####	Global Fst for all populations
	IDs = data[!is.na(data$Wing.length) ,'ID'] 
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
	results = list(wc.per.loc=wc.per.loc,wc.global=wc.global,IDs=IDs)
	save(results,file='outfiles/Fst.Wing.length.all.Rdata')

################################################################################################
####	Global Fst excluding SW slope
	IDs = data[!is.na(data$Wing.length) & !(data$Population %in% c('Macate','Pacar','Quichas','Huamatanga')),'ID'] 
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
	results = list(wc.per.loc=wc.per.loc,wc.global=wc.global,IDs=IDs)
	save(results,file='outfiles/Fst.Wing.length.SWslopeExcluded.Rdata')
	
################################################################################################
################################################################################################
##	PC1.plumage
################################################################################################
################################################################################################
cat('\n\n\n\n Fst calculations: PC1.plumage...\n')

################################################################################################
####	Global Fst for all populations
	IDs = data[!is.na(data$PC1.plumage) ,'ID'] 
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
	results = list(wc.per.loc=wc.per.loc,wc.global=wc.global,IDs=IDs)
	save(results,file='outfiles/Fst.PC1.plumage.all.Rdata')

################################################################################################
####	Global Fst excluding SW slope
	IDs = data[!is.na(data$PC1.plumage) & !(data$Population %in% c('Macate','Pacar','Quichas','Huamatanga')),'ID'] 
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
	results = list(wc.per.loc=wc.per.loc,wc.global=wc.global,IDs=IDs)
	save(results,file='outfiles/Fst.PC1.plumage.SWslopeExcluded.Rdata')


	