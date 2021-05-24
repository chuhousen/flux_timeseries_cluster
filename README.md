# Flux harmonic analysis across AmeriFlux network
## Summary
This project aims to use the key time-series features (e.g., diurnal, seasonal dynamics) of flux & met variables, to explore the similarity among sites across AmeriFlux.

### File sources
All available AmeriFlux BASE files (i.e., latest version 20210331), and down-select based on the data availability. 

### Target variables
NEE, LE, H, NETRAD, SW_IN, TA, SWC, VPD, USTAR. If multiple levels/locations are present at a site, use only the top-level & aggregated one. Use only the non-filled variables from the data files.

### Processing summary
- Select sites that have minimal 3 year data record.   
- Filter data by expected physical ranges. 
- Most variables are used as provided, except for NEE:
  - Use NEE directly if NEE is provided.
  - If NEE is not provided, assume NEE = FC + SC from following sources:
    - FC, filtered using site-specific ustar thresholds. 
    - For non-forest sites, use SC if provided or assume SC = 0.
    - For forest sites, use SC if provided. If SC is not provided, then calculate SC based on top-level CO2.
      - Drop a forest site if neither SC nor CO2 (& its measurement height) is provided. 
- Aggregate data into the desired composite diurnal-seasonal time series
  - Downscale to hourly (24 points per day), 
  - Group into 24 non-overlapped windows in a calendar year (~15 days per window, e.g., Jan 1-15, 16-30...).
  - The median of the diurnal time series for each window are calculated. 
  - All available years are used together.
  - The final composite time series is a multiyear composite diurnal-seasonal time series with 24*24 time steps (points).  
- Fill short remaining gaps by interpolation within each window and then across the windows. 
  - Skip a site-variable if more than 20% of gaps
- Down-select sites meeting the following data availability
  - At least 3 out of NEE, LE, H, USTAR available
  - At least 4 out of SW_IN, NETRAD, TA, VPD, SWC available
  - The latest version (20210331) contain 248 sites.
- Cluster time series across all available sites
  - Cluster based on each variable
  - Cluster based on multiple variables together
    - All flux variables (NEE, LE, H, USTAR)
    - All met variables (SW_IN, TA, VPD, SWC)
    - All flux and met variables
  - Dynamic Time Warping (DTW) is adopted to calculate the similarity (or distance)
  - Hierarchical clustering is done based on the DTW distance (average linkage) and to generate the polygenesis trees.
  - The Cluster Validity Indices are used to trim the tree.  
  
----

## File Directory
- R\ all R functions and workflows
- diurnal-seasonality\<version>\
  - ALL_BASE_site_list.csv: A full list of sites, site general information, data availability (before/after gap-filling)
  - ALL_BASE_site_short_list.csv: A short list of sites, aftering down-selecting based on data availability.
  - <SITE_ID>_<ORIGINAL_TIME_RESOLTION>_<BASE_VERSION>_MEDIAN.csv: Composite median of diurnal-seasonal time series of a site
- cluster-output\<version>\
  - ALL_BASE_site_list2.csv: A similar full list of sites, with additional info of ecoregions and clustering groups.
    - Ecoregions are pulled from the following sources, based on a site’s geo-location. 
      - North America Terrestrial Ecological Regions (http://www.cec.org/tools-and-resources/map-files/terrestrial-ecoregions-level-iii)
      - Global Terrestrial Ecological Regions (https://databasin.org/datasets/68635d7c77f1475f9b6c1d1dbe0a4c4c)
    - Clustering groups are the corresponding groups after trimming the dendrograms. Note: it’s per-variable based.
  - AMF-diurnal-seasonal-cluster-<TARGET_VARIABLE>.png: Diurnal-seasonal time series by clustering groups, a mean composite diurnal-seasonal time series is provided in colors for each group. The site-specific diurnal-seasonal time series is provided in a grey color. 
  - AMF-diurnal-seasonal-cluster-<TARGET_VARIABLE>_subgroup_<GROUP_NUMBER>.png 
    - Break-down figures for clustering groups that have more than 5 sites. The group number is specified in the filename. Diurnal-seasonal time series from all sites in the same cluster group.
    - Group-specific mean-composite time series is provided in black. 
  - AMF-diurnal-seasonal-cluster-<TARGET_VARIABLE>-tree.png: Clustering trees color-coded by groups
  - AMF-diurnal-seasonal-cluster-<TARGET_VARIABLE>-compiled_ts.csv: Diurnal-seasonal time series used in clustering
  - AMF-diurnal-seasonal-cluster-<TARGET_VARIABLE>-distMatrix.csv: Distance matrix used in clustering

----

## Collaborator
- Housen Chu: Lawrence Berkeley National Lab
- David E Reed: University of Science and Arts of Oklahoma
- Brad G Peter: The University of Alabama

