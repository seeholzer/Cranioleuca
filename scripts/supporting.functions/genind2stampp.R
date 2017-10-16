#=============================================================================================================#
# Script created by Glenn Seeholzer, contact at seeholzer.glenn@gmail.com
# Script created in version R 3.2.3 
# This script: wrapper script for to convert genind object to STAMPP input file format
# Usage notes: requires adegenet, StAMPP, and hierfstat
#=============================================================================================================#

genind2stampp = function(gen){
	require(adegenet)
	require(StAMPP)
	require(hierfstat)
	df = genind2hierfstat(gen)
	df[is.na(df)] = '-9'
	df = apply(df,1:2,function(x)gsub('2','B',x))
	df = apply(df,1:2,function(x)gsub('1','A',x))
	df = data.frame(ID=indNames(gen),df)
	head(df)
	#convert to STAMPP object
	input = data.frame(df[,1:2],Ploidy=2,Format='BiA',df[,-c(1:2)])
	input = input[order(input[,'pop']), ]
	gdata = stamppConvert(input, "r")
	return(gdata)
}
