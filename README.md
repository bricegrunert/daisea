# daisea2
A package for decomposing total non-water absorption spectra into CDOM + non-algal particulates (adg) and phytoplankton absorption (aph).

DAISEA2 is an update to the original DAISEA package. It utilizes genetic algorithms combined with Gaussian decomposition to separate total non-water absorption spectra into an estimate of adg and aph. DAISEA2 allows the user to select either an exponential or hyperbolic relationship for CDM. It is also able to handle spectra starting at or near 400 nm (e.g. AC-S data), rather than requiring total absorption down to 350 nm. 

The primary function is daisea2 (Derivative Analysis and Iterative Spectral Evaluation of Absorption). 

Please see the manuscript: A hyperspectral approach for retrieving inherent optical properties, phytoplankton pigments, and associated uncertainties from non-water absorption (Grunert, Ciochetto and Mouw, https://doi.org/10.3389/fmars.2025.1549312). 

## **Primary function**
**daisea2**

*Sub-functions*

*daisea2_get_started*
Use this function first! It will get you started using provided example data. 

*daisea_adg_ga*
Genetic algorithm for initial fit of adg. 

*daisea_apply_cdom_model*
Calculate predicted value (yhat) after fitting adg. 

*daisea_combined_fit_ga*
Genetic algorithm for combined fit of adg and Gaussian decomposition of aph. 

*daisea_combined_model_ga*
Models of adg and Gaussian curves for use by the genetic algorithms. 

*daisea_find_wavelengths*
Checks for critical wavelengths required by daisea2. For example, data must exist between 400 nm and 700 nm. 

*daisea_fit_cdom_nap*
Least squares fitting of adg. 

*daisea_gauss_peak_locations*
Initial Gaussian decomposition of aph after initial fit of adg. 

*daisea_gaussian_fit_ga*
Genetic algorithm for Gaussian decomposition of aph alone. 

*daisea_gaussian_fit*
Optimizes a least squares fit of Gaussian decomposition of aph after daisea_gauss_peak_locations.

*daisea_gaussian_model_ga*
DAISEA models for Gaussian fits with genetic algorithm. 

*daisea_initial_adg*
Identify portion of the curve most associated with adg and perform initial least squares fit. 

*daisea_wavelength_interval_check*
DAISEA2 requires a consistent wavelength interval. It works best with an interval of 1 nm. This function checks to make sure the wavelength interval is consistent and deals with NaNs as per user input. 

*daisea_models*
Models for adg and Gaussian spectra. 

## **How to cite this package**
This algorithm/package was released with an accompanying publication: A hyperspectral approach for retrieving inherent optical properties, phytoplankton pigments, and associated uncertainties from non-water absorption (Grunert, Ciochetto and Mouw, https://doi.org/10.3389/fmars.2025.1549312). Please cite this paper when using this package for presentations, publications, etc.
