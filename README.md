# Network of Networks: Time-Series Clustering of AmeriFlux Sites
## Summary
This project aims to use the key time-series features (e.g., diurnal, seasonal dynamics) of flux & met variables, to explore the similarity among sites across AmeriFlux. This repository supplements the manuscript *Network of Networks: Time-Series Clustering of AmeriFlux Sites* by Reed, Chu, Peter, Chen, et al. [in prep]. Journal TBD.

### File sources
All available AmeriFlux BASE files (i.e., accessed on 2023/09/22), and down-select based on the data availability. 

### Target variables
- Flux variables: net ecosystem exchange (NEE) of carbon dioxide, sensible heat flux (H), latent heat flux (LE), friction velocity (USTAR, i.e., a measure of momentum flux).
- Environmental variables: net radiation (NETRAD), air temperature (TA), vapor pressure deficit (VPD), soil water content (SWC). 
- If multiple levels/locations are present at a site, use only the top-level & aggregated one. Use only the non-filled variables from the data files.

### Processing summary
- Select sites that have minimal 3 year data record (at least 3 growing seasons).   
- Filter data by expected physical ranges. 
- Most variables are used as provided, except for NEE:
  - Use NEE directly if NEE is provided.
  - If NEE is not provided, assume NEE = FC + SC from following sources:
    - FC, filtered using site-specific USTAR thresholds. 
    - For non-forest sites, use SC if provided or assume SC = 0.
    - For forest sites, use SC if provided or calculate SC from top-level CO2.
      - Drop a forest site if neither SC nor CO2 (& its measurement height) is provided. 
- Aggregate data into the desired composite diurnal-seasonal time series
  - Downscale to desired resolutions (e.g., 24 (hourly), 12 (bi-hourly), 8 (3-hourly) time steps per day), 
  - Group data into desired non-overlapped windows in a calendar year (e.g., 15, 7, 5 days per window).
  - The median of the diurnal time series for each window are calculated. 
  - All available years are used together.
  - The final composite time series is a multi-year composite diurnal-seasonal time series (e.g., 12 * 52 = 624 time steps for bi-hourly & 7-day aggregation).  
- Fill short remaining gaps by interpolation within each window and then across the windows. 
  - Skip a site-variable if more than 20% of gaps
- Down-select sites meeting the following data availability
  - At least 1 out of all 4 flux variables
- Cluster time series across all available sites
  - Calculate Dynamic Time Warping (DTW) distance as a measure of similarities among each site's time series
  - Adopt hierarchical clustering to construct the hierarchy of site clusters based on the sites’ DTW distances
  - Cluster based on each variable (univariate)
  - Cluster based on all flux variables (multivariate, NEE, LE, H, USTAR)
  - The Cluster Validity Indices (CVI) are used to trim the tree.  
    - Six CVIs include Sil, CH, DB, DBstar, D, COP
    - Search optimal cluster number between 30 to 75
- Examine dependency of sites' DTW distances on IGBP classifications, ecoregions, and spatial proximity.
- Calculate Harmonic Uniqueness Parameter (HUP) 
  - Mean of normalized DTW distances to all other sites
  - Use all target flux variables, i.e., NEE, LE, H, USTAR
  
----

## File Directory
- /R all R functions and workflows
- /ecoregion source files of ecoregions

----

## Extended Resources
- An interactive web map https://cartoscience.users.earthengine.app/view/flux-networks

----

## Collaborator
- Housen Chu: Lawrence Berkeley National Lab
- David E Reed: Yale University
- Brad G Peter: University of Arkansas

