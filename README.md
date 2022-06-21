# FAIW-PCA

% Fiber architecture informed synchrony mapping
% Syntax:  calculate_SYN(f_fmri_total, f_dwi, f_nodif_mask, par);
% Author: Yu zhao
% email: yuzhao923@gmail.com
% Jun. 2022; 
% This toolbox was created for the following paper: 
% Detection of functional activity in brain white matter using fiber architecture informed synchrony mapping. DOI: 10.1101/2022.02.23.481698 

% Fiber architecture informed gragh is generated using the fuction,generate_odf_informed_gragh,
% which adapted from the previous report (Abramian, David, et al. "Diffusion-informed
% spatial smoothing of fMRI data in white matter using spectral graph filters." Neuroimage 237 (2021): 118095.)

### Software prerequisits

- [Matlab](https://se.mathworks.com/products/matlab.html)
- [SPM toolbox](https://www.fil.ion.ucl.ac.uk/spm/): used for loading and saving NIfTI volumes.
- [DSI Studio](http://dsi-studio.labsolver.org/): used for generating diffusion ODFs from DWI data.

### Setup

1. Install SPM and add it to your Matlab path.
2. Install DSI Studio.
3. Download and extract DSS. Add to your Matlab path the code folder and all subfolders.
4. DSI Studio is used from Matlab through `system` calls. The path of the DSI Studio executable has to be provided in the `DSI_STUDIO_PATH` environmental variable. This can be set manually outside Matlab, but the easiest way is to do it from Matlab. Modify the provided `setup_dsi.m` function to have the path point to your DSI Studio executable, and make sure to run this function before using DSS every time Matlab is started. The process can be automated by calling this function from the `startup.m` script.

### Usage

The main function of the toolbox is `mian_FAIW_PCA.m`:
