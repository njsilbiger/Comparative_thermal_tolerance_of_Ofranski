[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3333002.svg)](https://doi.org/10.5281/zenodo.3333002)
# Comparative thermal performance of the reef-building coral *Orbicella franksi* at its latitudinal range limits

### Authors: Nyssa J. Silbiger*, Gretchen Goodbody-Gringley, John F. Bruno, Hollie M. Putnam

### Journal: Marine Biology 

### Link: (https://link.springer.com/article/10.1007/s00227-019-3573-6)[https://link.springer.com/article/10.1007/s00227-019-3573-6]

### Funding:
This research was funded in part by the National Science Foundation (grant OCE #1737071 to JFB), California State University, and the Pembroke Foundation International (to GGG and HMP).

### Folders:

**Data**\
BermudaRates.csv (calculated photosynthesis and resipration rates for Bermuda)
PanamaRates.csv (calculated photosynthesis and resipration rates for Bermuda)

- **Mapping**\
Shapefiles for Panama and Bermuda for the map in Figure 2

- **TemperatureData**\
Temperature time series data for Panama and Bermuda from long-term stations

- **Bermuda**
  - **MetaData**\
All metadata files used in this analysis\
Nubbin_Sample_Info_T0_Bermuda_QC (Combined sample info for each individual nubbin)\
SA_Master_QC.csv (Surface area data for the corals)\
WaterChem_TPC_Bermuda_QC.csv (raw pH, salinity, and water volume for each sample)\

  - **Respirometry**\
All the raw respirometry files from the PreSens output

  - **PI_curve**\
All the raw respirometry data and metadata for creating the PI curves in Bermuda

  - **TAData**\
NECData.csv (raw TA data and calculated net calcification rates for the Bermuda corals)

- **Panama**
  - **MetaData**\
Nubbin_Sample_Info_T0_Panama_QC.csv (Combined sample info for each individual nubbin)\
SurfaceArea_QC.csv (Coral Surface area data for Panama)\

  - **Respirometry**\
All the raw respirometry files from the PreSens output

  - **PI_curve**\
All the raw respirometry data and metadata for creating the PI curves in Panama

**Output**
Diagnostics from all of the models including Bayesian P-values and quantiles.
- **Bermuda**
  - **Photo_Resp_Output**\
Files for the raw output of the respirometry rates from the LoLinR package for Bermuda

  - **PIC_Output**\
Files for the raw output of the respirometry rates from the LoLinR package for Bermuda PI curves

- **Panama**
  - **Photo_Resp_Output**\
Files for the raw output of the respirometry rates from the LoLinR package for Panama

  - **PIC_Output**\
Files for the raw output of the respirometry rates from the LoLinR package for Panama PI curves

- **MSPlots**\
Plots used in the manuscript

- **traceplots**\
Trace plots and autocorrelation plots to check Bayesian output

**Scripts**
- All_PI_Curves.R (calculate photosynthesis-irradiance curves for both Panama and Bermuda)
- BayesInParallel.R (R script for running models for photosynthesis and respiration between the two locations using parallel computing)
- JAGSnestedLocation (JAGS code for comparing rates between locations)
- JAGSnetstedRates (JAGS code for comparing rates between different organismal functions)
- modelbyrates.R (R script for running models comparing among the three organismal functions)
- Segmented_Respirometry_Bermuda.R (R script to take raw respirometry data and calculate rates for Bermuda)
- SegmentedResp_Panama.R (R script to take raw respirometry data and calculate rates for Panama)
- TemperatureScript.R (R script to make the temperature plots and the maps)
