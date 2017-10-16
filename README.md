This repository contains all the data and scripts necessary to duplicate the Pst-Fst and 
commonality analyses in Seeholzer and Brumfield (2017). The scripts in the directory 
'scripts' will duplicate Figures 4-6, S8 and the results in Tables 2-4, and write them to 
the outfile directory. Below are descriptions of each sub-directory or file in 
hierarchical and alphabetical order.




### SDM - Directory containing MAXENT models and presence localities.
	MAXENT_SDM_Cranioleuca_antisiensis.asc
			Raw logistic output from the MAXENT model. See Supporting Methods for more 
			information on how this model was created.

	TNC_CentralAndeanWetPuna.shp
			Shapefile of the Central Andean Wet Puna ecoregion used to mask 
			MAXENT_SDM_Cranioleuca_antisiensis.asc

	presence.localities.txt
			Tab-delimited text file with coordinates in decimal degrees for presence 
			localities used in MAXENT model. See Supplementary methods for more 
			information on the data sources and filtering schemes used to obtain this 
			final list. 

### SNPdata
	Directory containing SNP data in various formats
		
		SNPdata.txt
				Contains unphased SNPs for 5154 loci across 175 individuals include 
				three Cranioleuca curtata as outgroups (cra.cur.B43876, cra.cur.B6032, 
				cra.cur.B8175). Columns correspond to locus ID (1 through 5154), the 
				key for the identification of the nucleotides of the SNPs, and the IDs 
				of the 175 individuals included in the dataset. The IDs correspond to 
				the IDs in Data.Table.txt.

		STRUCTURE.input.file.txt
				This the SNPdata.txt file converted into the input format for 
				STRUCTURE. Although we did not use STRUCTURE in our analyses, it is a 
				common input file format that can be converted easily to the input 
				formats for the R package Adegenet and for ADMIXTURE.	

		SNPdata.genind.Rdata
				R data object in genind format (package adgenenet). This is the 
				primary file used in the R scripts for Pst-Fst and commonality 
				analysis.

### GBS_data
	This directory contains the original data files sent by the Cornell Institute 
	for Genomic Diversity to the authors. This data is the output of the UNEAK pipeline,
	not the raw GBS reads. The authors conducted additional processing as explained
	in the methods of Seeholzer & Brumfield 2017.


### bioclim.veg.layers
	Directory containing grid files containing rasters of altitude, 19 bioclim 
	variables, and two vegetation indices (TREE and NDVI) used in various analyses. 
	These rasters are cropped to the study area for Cranioleuca antisiensis.

### outfiles
	Directory to where outfiles from the scripts will be written.

### scripts
	Directory containing all scripts necessary to dupicate Pst-Fst and commonality 
	analysis. Each script should be run line-by-line.

		supporting.functions
				Directory containing functions called in the scripts above
		
		MMRR.CommonalityAnalysis.DispersalDistance.R
				R script to duplicate results in Figures 5, 6, S2, S3, S7, 
				and Tables 3 & 4
		
		MMRR.CommonalityAnalysis.Geographic.Distance.R
				R script to duplicate results in Figure S8 & Table S2 where geographic 
				distance is substituted for dispersal distance.
				
		PstFst.1.Fst.calculations.R
				R script to create Fst matrices used in following PstFst scripts.
		
		PstFst.2.Pst.calculations.R
				R script to calculate Pst values, confidence intervals, and critical
				values of c/h2.
		
		PstFst.3.PstFst.plotting.R
				Produces Figure 4 (upper and lower panels) and Table 2.
			
		PstFst.4.1.sensitivity.analysis.R
				Generates Fst and Pst data tables for Pst-Fst sensitivity analysis 
				described in Discussion section "Pst-Fst"

		PstFst.4.2.sensitivity.analysis.R
				Generates Figure S6 from the output of PstFst.4.1.sensitivity.analysis.R

### BayEnv_outlier_loci.txt
		Tab-delimited text file with locus names of outliers identified in BayEnv run

### Data.Table.txt
		Tab-delimited text file with geographic, morphologic, and plumage data for all 
		individuals used in study. Meaning of selected column names below
			ID - unique individual ID
			Population Code - Numeric code used for each populations (1-19) 
			Transect Position - Position in km of each population along an orthogonal
				regression of each population's coordinates with origin in north and 
				terminus in south
			Verbatim Elevation - elevation recorded on specimen tags (generally from GPS)
			Bioclim Elevation - elevation for each lat/long from BioClim 30-arcsec raster  
			Age - a = adult, j = juvenile, u = unknown (but likely adult)
			Sex - m = male, f = female
			PC1 morphology - PC1 for all the linear morphological measurements
			PC1 plumage - PC1 of plumage data

### Fst.all.05outliersexcluded.txt
		Tab-delimited text file with pairwise global Fst matrix between all populations 
		based after excluding 'outlier' loci

### README
		This document. 

### cost.raw.txt
		Tab-delimited text file of least cost path distances between populations based on 
		MAXENT_SDM_Cranioleuca_antisiensis.asc

### cost.refined.txt
		Tab-delimited text file of least cost path distances between populations based on 
		MAXENT_SDM_Cranioleuca_antisiensis.asc masked with TNC_CentralAndeanWetPuna.shp 
		(see lines 121-148 of script MMRR.CommonalityAnalysis.R)
		