# daisea
A package for decomposing total non-water absorption spectra into CDM and phytoplankton absorption.

The daisea package utilizes derivative analysis, iterative spectral evaluation and Gaussian decomposition of total non-water absorption spectra to provide an estimate of colored detrital material (CDM) absorption and phytoplankton absorption.

The primary function is daisea (Derivative Analysis and Iterative Spectral Evaluation of Absorption) as described in the manuscript Deriving inherent optical properties from decomposition of hyperspectral non-water absorption (Grunert et al. 2019, https://doi.org/10.1016/j.rse.2019.03.004). We encourage users to refer to this document for a detailed description of steps along with a flow chart; additional information on function output can be found in this document as well.  As the main executable function, daisea utilizes a variety of custom sub-functions as well as the following Matlab Toolboxes: Curve Fitting Toolbox and Signal Processing Toolbox. Input spectra should be hyperspectral, ideally at a wavelength range of ≤ 5 nm. Spectral range should extend to approximately 350 nm and up to 690 nm or greater. All function files include a detailed description of required input variables.

## **Primary function**
**daisea**

*Sub-functions*

*build_daisea_model*

Function utilizes an estimated single exponential model and pre-defined number of estimated Gaussian curves to fit total non-water absorption using a least squares approach.

*build_phy_model*

Function optimizes a pre-defined number of estimated Gaussian curves for a given signal using a least squares approach.

*cdom_model_lam0_noK*

Function models a colored dissolved organic matter (CDOM), non-algal particulate (NAP), or CDM absorption spectra with a single exponential following a defined initial wavelength, absorption at that wavelength, and no offset value (K value).

*spectral_slope_nogauss_noK*

Function estimates spectral slope (S) for a CDOM, NAP or CDM absorption spectra assuming the spectra can be fit with a single exponential model, no Gaussian components and K value.

*cdom_model_0gaussian_noK*

Function models an absorption spectra best approximated with a single exponential and no Gaussian components with no K value.

*cdom_model_1gaussian_noK*

Function models an absorption spectra best approximated with a single exponential and one Gaussian component with no K value.

*cdom_model_2gaussian_noK*

Function models an absorption spectra best approximated with a single exponential and two Gaussian components with no K value.

*cdom_model_3gaussian_noK*

Function models an absorption spectra best approximated with a single exponential and three Gaussian components with no K value.

*cdom_model_4gaussian_noK*

Function models an absorption spectra best approximated with a single exponential and four Gaussian components with no K value.

*cdom_model_5gaussian_noK*

Function models an absorption spectra best approximated with a single exponential and five Gaussian components with no K value.

*cdom_model_6gaussian_noK*

Function models an absorption spectra best approximated with a single exponential and six Gaussian components with no K value.

*cdom_model_7gaussian_noK*

Function models an absorption spectra best approximated with a single exponential and seven Gaussian components with no K value.

*cdom_model_8gaussian_noK*

Function models an absorption spectra best approximated with a single exponential and eight Gaussian components with no K value.

*cdom_model_9gaussian_noK*

Function models an absorption spectra best approximated with a single exponential and nine Gaussian components with no K value.

*cdom_model_10gaussian_noK*

Function models an absorption spectra best approximated with a single exponential and ten Gaussian components with no K value.

*cdom_model_11gaussian_noK*

Function models an absorption spectra best approximated with a single exponential and eleven Gaussian components with no K value.

*cdom_model_12gaussian_noK*

Function models an absorption spectra best approximated with a single exponential and twelve Gaussian components with no K value.

*cdom_model_13gaussian_noK*

Function models an absorption spectra best approximated with a single exponential and thirteen Gaussian components with no K value.

*cdom_model_14gaussian_noK*

Function models an absorption spectra best approximated with a single exponential and fourteen Gaussian components with no K value.

*cdom_model_15gaussian_noK*

Function models an absorption spectra best approximated with a single exponential and fifteen Gaussian components with no K value.

*cdom_model_16gaussian_noK*

Function models an absorption spectra best approximated with a single exponential and sixteen Gaussian components with no K value.

*gauss*

Function models a single Gaussian curve.

*gauss1*

Function models a single Gaussian curve.

*gauss2*

Function models two Gaussian curves.

*gauss3*

Function models three Gaussian curves.

*gauss4*

Function models four Gaussian curves.

*gauss5*

Function models five Gaussian curves.

*gauss6*

Function models six Gaussian curves.

*gauss7*

Function models seven Gaussian curves.

*gauss8*

Function models eight Gaussian curves.

*gauss9*

Function models nine Gaussian curves.

*gauss10*

Function models ten Gaussian curves.

*gauss11*

Function models eleven Gaussian curves.

*gauss12*

Function models twelve Gaussian curves.

*gauss13*

Function models thirteen Gaussian curves.

*gauss14*

Function models fourteen Gaussian curves.

*gauss15*

Function models fifteen Gaussian curves.

*gauss16*

Function models sixteen Gaussian curves.


## **How to use example spectra**

for ii=1:length(example_spectra)
	output(ii)=daisea(example_spectra(ii).at_nw,example_spectra(ii).wavelength,350,700);
end

## **How to plot results**

for ii=1:length(output)
figure;
plot(output(ii).lam,output(ii).at_nw,’--k’)
hold on;
plot(output(ii).lam,output(ii).agd_final_estimate,’-‘,’color’,[140 81 10]/255)
plot(output(ii).lam, output(ii).agd_final_estimate+output(ii).aphy_final_estimate,’-‘,’color’,[35 139 69]/255)
end

## **How to cite this package**
This algorithm/package was released with an accompanying publication, Deriving inherent optical properties from decomposition of hyperspectral non-water absorption (Grunert et al. 2019, https://doi.org/10.1016/j.rse.2019.03.004). Please cite this paper when using this package for presentations, publications, etc.
