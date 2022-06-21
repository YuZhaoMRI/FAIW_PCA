% Fiber architecture informed synchrony mapping
% Syntax:  calculate_SYN(f_fmri_total, f_dwi, f_nodif_mask, par);
% Author: Yu zhao
% email: yuzhao923@gmail.com
% Jun. 2022; 
% This toolbox was created for the following paper: 
% Detection of functional activity in brain white matter using fiber architecture informed synchrony mapping. DOI: 10.1101/2022.02.23.481698 

% Fiber architecture informed gragh is generated using the fuction, generate_odf_informed_gragh,
% which adapted from the previous report (Abramian, David, et al. "Diffusion-informed
% spatial smoothing of fMRI data in white matter using spectral graph filters." Neuroimage 237 (2021): 118095.)

 

root='T:\HCP_3T_64S\preprocessed_data_mni2acpc\117930';
% f_dwi - Path of NIfTI volume series containing DWI data (can be
%  compressed).
  f_dwi  =fullfile(root, 'Diffusion\data.nii.gz');


% f_mask - Path of uncompressed NIfTI volume containing binary white
%      matter mask.
f_nodif_mask =fullfile(root, 'Diffusion\nodif_brain_mask.nii');

% if you change the parameters of generate_odf_informed_gragh, please do this
f_odf=fullfile(root, 'Diffusion\odfs.mat');
delete(f_odf)

 

%    fMRI volumes need to be registered to the DWI space!!!
%    f_out - Path of  fMRI volume series. Data in which the nuisances of head
%    motion and CSF have been regressed out.
%    
%    Note that smooting and bandpass filtering will be done in our codes.
%    the two HCP runs with opposite phase encoding directions will be
%    concatenated.
% motor task
f_fmri_motor1 = fullfile(root, 'tfMRI_MOTOR\tfMRI_MOTOR_LR_ACPC1.25mm.nii.gz');
f_fmri_motor2 = fullfile(root, 'tfMRI_MOTOR\tfMRI_MOTOR_RL_ACPC1.25mm.nii.gz');
% rest state
f_fmri_rest1 = fullfile(root, 'rfMRI_REST1\rfMRI_REST1_LR_ACPC1.25mm.nii.gz');
f_fmri_rest2 = fullfile(root, 'rfMRI_REST1\rfMRI_REST1_RL_ACPC1.25mm.nii.gz');

f_fmri_total={f_fmri_motor1,f_fmri_motor2,f_fmri_rest1,f_fmri_rest2};





%% estimate of SYN map


par.tau=0.30;% The duration of graphicl diffusion
par.num_neigh=64; %% For white matter, ~95% of the signals in graphicl diffusion are inclused in 64 voxels   
par.FWHM=3;%smoothing [mm]
par.band=[0.01 0.10];%% Bandpass filtering [Hz]
par.num_tc=284;% Time frames of the motor-task fMRI run  
par.alpha=0.905;%% Soft threshold for graph weights. Higher values result in more
%                  anisotropic filters. Takes values in [0,1] (default: 0.9).
par.neigh_graph=5;%  Size of neighborhood shell defining the set of neighbors for
%                 each vertex. Takes values in {3,5} (default: 5).



generate_odf_informed_gragh(f_nodif_mask, f_dwi, par.alpha, par.neigh_graph);
calculate_SYN(root,f_fmri_total, f_dwi, f_nodif_mask, par);



 
 
 