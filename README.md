# Holte_MixedLayerDepth
Adaptation of the Holte et al. (2009) ARGO blended mixed layer depth code to discrete profiles. MATLAB CODE.

Holte, James and Lynne Talley, 2009: A new algorithm for finding mixed layer depths with applications to Argo data and Subantarctic Mode Water formation, Journal of Atmospheric and Oceanic Technology 26(9), 1920-1939, doi:10.1175/2009JTECHO543.1.

Holte et al. (2009) provided a suite of Matlab routines specifically to support users of the ARGO profiling float datasests. Their blended algorithm (their file findmld.m) is adapted here to work as a standalone function (findmld_2.m), with inputs salinity, temperature and potential density. The code uses threshold and gradient approaches to estimate mixed layer depth and, in this application, returns the deepest MLD estimates from the temperature-, salinity- and density-derived values. Refer to Holte et al. (2009) for details of the calculations and choice of constants.

You are free to use and adapt this code. 
Please make sure you acknowledge Holte et al. (2009)!
