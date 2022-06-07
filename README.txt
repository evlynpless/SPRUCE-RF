Predicting environmental and ecological drivers of human population structure

Summary:
Have you ever encountered a migration surface and wondered why some areas were barriers or corridors to gene flow? SPRUCE is a flexible machine learning approach for determining what spatial variables are most predictive of a migration surface, particularly those created by the MAPS software (https://github.com/halasadi/MAPS).
Here we provide a guide and associated scripts for implementing SPRUCE. This is a flexible pipeline that can be adjusted according to the needs of the user.

Software and prerequisites:

Programs used in dealing with genomic data include Plink, SHAPEIT2, hap-ibd, and MAPS. Helpful tools for editing environmental data include GDAL and PKTOOLs. Various R packages are required (especially "raster", "plotmaps", and "randomForestSRC"); the specific packages required are listed in each script. These scripts are meant to provide full transparency about scripting procedure, and they will need updating to be run with new genetic and environmental data. 

Scripts and materials included:

SPRUCE_Guide.txt 

1B_Phase_and_Call_IBD.sh

1C_Interpolate_Position_to_cM.R

1E_IBDfile_to_sims.R

3_SPRUCE_RandomForest.R

GenerateSpatialVariable (folder)