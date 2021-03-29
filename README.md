# Flux harmonic analysis across AmeriFlux network
## Summary
This project aims to use the key time-series features (e.g., diurnal, seasonal dynamics) of flux & met variables, to explore the similarity among sites across AmeriFlux.

### File sources
All available AmeriFlux BASE files (i.e., 302 sites, version 20210319), and down-select based on the data availability (see Table).  

### Target variables
NEE, LE, H, NETRAD, SW_IN, TA, SWC, VPD, USTAR, WS. If multi-levels/-locations are present at a site, use only the top-level/aggregated one. 

### Processing summary
- All data are first filtered by expected physical ranges. 
- Most variables are used as provided, except for NEE:
  - Use NEE directly if NEE is provided.
  - Assume NEE = FC + SC if NEE is not provided, and use FC and SC from following sources:
    - FC, filtered using site-specific ustar thresholds. 
    - For non-forest sites, use SC if it’s provided. Or, assume SC = 0.
    - For forest sites, use SC if it’s provided. If SC is not provided, then calculate SC based on top-level CO2.
      - Drop a forest site if neither SC nor CO2 (& its measurement height) is provided. 
  - The data are aggregated into the desired composite diurnal-seasonal time series, i.e., to reduce the data points, fill the short gaps, & keep the key features of the time series (i.e., diurnal & seasonal dynamics)
  - For this version, the original time series is first downscaled to hourly (24 points per day), and then grouped into 24 non-overlapped windows in a calendar year (~15 days per window, e.g., Jan 1-15, 16-30...). The median of the diurnal time series for each window are calculated. All available years are used together.
  - The final composite time series is a multi-year composite diurnal-seasonal time series with 24*24 time steps (points), starting with 0h-1h, 1h-2h, 2h-3h… of the first window (Jan 1st-15th), then 0h-1h… of the second window (Jan 16th-30th), and all the way to 22h-23h, 23h-0h of the last window (DOY 346-365).  
  - The composite time series, if there’s any remaining gap, is then filled by interpolation within each window and then across the windows. If there’s more than 20% of gaps within each window and 20% of gaps across windows, then don’t do gap-filling.
  - The time series clustering is done for each of the target variables on the composite time series from all available sites (excluding those still have remaining gaps in the composite time series)
  - Dynamic Time Warping (DTW) is adopted to calculate the similarity (or distance) between each pair of the composite time series. Then the hierarchical clustering is done based on the DTW distance (average linkage) and to generate the polygenesis trees. Finally, all available sites are grouped into 30 groups after trimming the dendrograms.  

## File Directory
- R\ all R functions and workflows
- diurnal-seasonality\<version>\
  - **ALL_BASE_site_list.csv**: A full list of sites, site general information, data availability (before/after gap-filling)
  - **<SITE_ID>_<ORIGINAL_TIME_RESOLTION>_<BASE_VERSION>_diurnalseasonal-#.png**: Composite diurnal-seasonal time series of target variables, showing the median, 25th, and 75th percentile 
  - **<SITE_ID>_<ORIGINAL_TIME_RESOLTION>_<BASE_VERSION>_MEDIAN.csv**: Composite median of diurnal-seasonal time series of target variables
- cluster-output\<version>\
  - **ALL_BASE_site_list2.csv**: A full list of sites, site general information, data availability (before/after gap-filling), Ecoregions, and clustering groups.
    - Ecoregions are pulled from the following sources, based on a site’s geo-location. 
      - North America Terrestrial Ecological Regions (http://www.cec.org/tools-and-resources/map-files/terrestrial-ecoregions-level-iii)
      - Global Terrestrial Ecological Regions (https://databasin.org/datasets/68635d7c77f1475f9b6c1d1dbe0a4c4c)
    - Clustering groups are the corresponding 30 groups after trimming the dendrograms. Note: it’s per-variable based.
  - **AMF-diurnal-seasonal-cluster-<TARGET_VARIABLE>.png**: Diurnal-seasonal time series by clustering groups, a mean composite diurnal-seasonal time series is provided in colors for each group. The site-specific diurnal-seasonal time series is provided in a grey color. 
  - **AMF-diurnal-seasonal-cluster-<TARGET_VARIABLE>_subgroup_<GROUP_NUMBER>.png** 
    - Break-down figures for clustering groups that have more than 5 sites. The group number is specified in the filename. Diurnal-seasonal time series from all sites in the same cluster group.
    - Group-specific mean-composite time series is provided in black. 
  - **AMF-diurnal-seasonal-cluster-<TARGET_VARIABLE>-tree.png**: Clustering trees color-coded by groups
  - **AMF-diurnal-seasonal-cluster-<TARGET_VARIABLE>-compiled_ts.csv**: Diurnal-seasonal time series used in clustering
  - **AMF-diurnal-seasonal-cluster-<TARGET_VARIABLE>-distMatrix.csv**: Distance matrix used in clustering
----

## Collaborator
- David E Reed: University of Science and Arts of Oklahoma
- Housen Chu: Lawrence Berkeley National Lab
- Brad G Peter: The University of Alabama

