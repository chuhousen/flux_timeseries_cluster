# Flux harmonic analysis across AmeriFlux network
## Summary
This project aims to use the key time-series features (e.g., diurnal, seasonal dynamics) of flux & met variables, to explore the similarity among sites across AmeriFlux.

### File sources
All available AmeriFlux BASE files (i.e., latest version 20230922), and down-select based on the data availability. 

### Target variables
NEE, LE, H, NETRAD, TA, SWC, VPD, USTAR. If multiple levels/locations are present at a site, use only the top-level & aggregated one. Use only the non-filled variables from the data files.

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
  - Downscale to desired resolutions (e.g., 24, 12, 8 points per day), 
  - Group into desired non-overlapped windows in a calendar year (e.g., 15, 7, 5 days per window).
  - The median of the diurnal time series for each window are calculated. 
  - All available years are used together.
  - The final composite time series is a multiyear composite diurnal-seasonal time series (e.g., 24*24 time steps for hourly & 15 day aggregation).  
- Fill short remaining gaps by interpolation within each window and then across the windows. 
  - Skip a site-variable if more than 20% of gaps
- Down-select sites meeting the following data availability
  - At least 1 out of all 8 target variables
- Cluster time series across all available sites
  - Cluster based on each variable
  - Cluster based on multiple variables together
    - All flux variables (NEE, LE, H, USTAR)
    - All met variables (NETRAD, TA, VPD, SWC)
    - All flux and met variables
  - Dynamic Time Warping (DTW) & Euclidean distances are adopted to calculate the similarity (or distance)
  - Hierarchical clustering is done based on the DTW distance (average linkage) and to generate the polygenesis trees.
  - The Cluster Validity Indices are used to trim the tree.  
  
----

## File Directory
- R\ all R functions and workflows
- diurnal-seasonality\<version>\
- cluster-output\<version>\
  
----

## Collaborator
- Housen Chu: Lawrence Berkeley National Lab
- David E Reed: Yale University
- Brad G Peter: University of Arkansas

