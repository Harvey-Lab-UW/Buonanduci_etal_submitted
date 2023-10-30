# Buonanduci_etal_submitted
Data and code accompanying the manuscript 'Scaling severe fire patterns across fire sizes yields insights for data-sparse and infrequent-fire regimes' by Buonanduci, Donato, Halofsky, Kennedy, and Harvey; submitted for publication in Ecosphere.


[![CC BY 4.0][cc-by-shield]][cc-by]

This information is licensed under a
[Creative Commons Attribution 4.0 International License][cc-by]. Any user of these data ("User" hereafter) is required to cite it appropriately in any publication that results from its use. These data may be actively used by others for ongoing research, so coordination may be necessary to prevent duplicate publication. The User is urged to contact the authors of these data for questions about methodology or results.  The User is encouraged to consider collaboration or co-authorship with authors where appropriate. Misinterpretation of data may occur if used out of context of the original study. Substantial efforts are made to ensure accuracy of the data and documentation, however complete accuracy of data sets cannot be guaranteed. All data are made available as is. Data may be updated periodically and it is the responsibility of the User to check for new versions of the data. The authors and the repository where these data were obtained shall not be liable for damages resulting from any use or misinterpretation of the data.

[![CC BY 4.0][cc-by-image]][cc-by]

[cc-by]: http://creativecommons.org/licenses/by/4.0/
[cc-by-image]: https://i.creativecommons.org/l/by/4.0/88x31.png
[cc-by-shield]: https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg



For reproducibility, the following files are made available:

## Data files

#### fire_metrics.csv
This file contains high-severity patch size and structure metrics for all fire events included in the analysis. The following columns are included:

- **Fire_ID**: unique fire identifier, as assigned by MTBS
- **Fire_Name**: fire name, as assigned by MTBS
- **Year**: year in which fire occurred
- **NW_Cascadia**: binary indicator for whether fire occurred within Northwestern Cascadia (1) or elsewhere within the Northwest US (0)
- **Fire_Regime**: historical fire regime; either "Low" (frequent and low severity), "Mixed" (moderately frequent and mixed severity), or "High" (infrequent and high severity), as classified by LANDFIRE Fire Regime Group
- **prp_forest**: proportion of burned area that is forested, as classified by LANDFIRE Environmental Site Potential
- **fire_area**: total area burned (ha)
- **log_fire_area**: log10-transformed total area burned
- **HS_area**: total area burned at high severity (ha)
- **HS_prop**: proportion of area burned at high severity
- **HS_forest_area**: total forested area burned at high severity (ha)
- **HS_forest_prop**: proportion of area that was both forested and burned at high severity
- **patch_area_mean**: arithmetic mean of high-severity patch sizes (ha)
- **patch_area_AW_mean**: area-weighted mean of high-severity patch sizes (ha)
- **patch_area_gt1_mean**: arithmetic mean of high-severity patch sizes exceeding 1 ha in size
- **patch_area_gt1_AW_mean**: area-weighted mean of high-severity patch sizes exceeding 1 ha in size
- **log_patch_area_AW_mean**: log10-transformed area-weighted mean of high-severity patch sizes
- **log_patch_area_gt1_AW_mean**: log10-transformed area-weighted mean of high-severity patch sizes exceeding 1 ha in size
- **beta**: beta parameter for high-severity patch size distribution
- **psi**: psi parameter for high-severity patch size distribution
- **n_patches**: number of high-severity patches
- **n_patches_gte1**: number of high-severity patches equal to or exceeding 1 ha in size
- **total_core**: total core area (ha)
- **log_total_core**: log10-transformed total core area
- **SDC**: stand-replacing decay coefficient; parameter for distance-to-seed distribution
- **log_SDC**: log10-transformed SDC parameter


## Spatial data files

#### fire_perims.gpkg
This file contains polygon perimeters for all fire events. Coordinate reference system is NAD83 / Conus Albers (EPSG:5070). The following columns are included:

- **Fire_ID**: unique fire identifier, as assigned by MTBS
- **Fire_Name**: fire name, as assigned by MTBS
- **Year**: year in which fire occurred
- **NW_Cascadia**: binary indicator for whether fire occurred within Northwestern Cascadia (1) or elsewhere within the Northwest US (0)
- **Fire_Regime**: historical fire regime; either "Low" (frequent and low severity), "Mixed" (moderately frequent and mixed severity), or "High" (infrequent and high severity), as classified by LANDFIRE Fire Regime Group


#### region_NW.gpkg
This file contains a polygon perimeter for the Northwest US study region. Coordinate reference system is NAD83 / Conus Albers (EPSG:5070). 


#### region_NW_Cascadia.gpkg
This file contains a polygon perimeter for the Northwestern Cascadia study region. Coordinate reference system is NAD83 / Conus Albers (EPSG:5070). 


#### states.gpkg
This file contains US state outlines, for plotting purposes. Coordinate reference system is NAD83 / Conus Albers (EPSG:5070). 


#### WA_OR_line.gpkg
This file contains the border between Washington and Oregon as a polyline, for plotting purposes. Coordinate reference system is NAD83 / Conus Albers (EPSG:5070). 


#### CCrest.gpkg
This file contains the Cascade Crest as a polyline, for plotting purposes. Coordinate reference system is NAD83 / Conus Albers (EPSG:5070). 


## R script files

#### install_quantregGrowth.R
Code for installing the appropriate version of the `quantregGrowth` package. Version `1.4-0` of the `quantregGrowth` package must be installed by the user to ensure all R scripts within this repository run as intended.

#### figure_study_region.R
Code for producing Figure 1 in the manuscript.

#### figure_simulation_workflow.R
Code demonstrating the simulation workflow and producing components of Figure 2 in the manuscript.

#### figure_scaling.R
Code for fitting quantile regression models and producing Figure 3 in the manuscript.

#### figure_simulation_core_DTS_dist.R
Code for simulating core areas and distance-to-seed distributions. Also produces Figure 4 and Table 2 in the manuscript.

#### figure_simulation_patch_dist.R
Code for simulating patch size distributions. Also produces Figure 5 and Table 3 in the manuscript.

#### figure_NWC_quantile_loss.R
Code for calculating quantile loss (prediction error) for the Northwestern Cascadia data using each fire regime-specific model. Produces Appendix S2, Figure S1.

#### truncated_lnorm_functions.R
Functions for fitting, plotting, and simulating from truncated lognormal patch size distributions.
