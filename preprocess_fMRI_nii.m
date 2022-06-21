function [S_voi_1,S_voi_2]=preprocess_fMRI_nii(f_fmri_total,mask_odf,I_mask,FWHM,bf,N4)

fprintf('Read fMRI data.\n');
%% read motor-task fMRI(nii format) data
f_fmri_1=f_fmri_total{1};
f_fmri_2=f_fmri_total{2};
% fprintf('Read fMRI data.\n');
S_4D=niftiread(f_fmri_1);
info=niftiinfo(f_fmri_1);
TR=info.PixelDimensions(4);
voxelsize=info.PixelDimensions(1);
S_4D=S_4D(:,:,:,1:N4); % To keep the same time frames as that of resting-state data
[N1,N2,N3,N4]=size(S_4D); 

%% Important!!! The coordinates is adjusted for the graph dreived from DSI-studio 
  S_4D=flip(S_4D,1);
  
  
%  S_3D= squeeze(S_4D(:,:,:,1));

for ii=1:N4
    S_4D(:,:,:,ii)=imgaussfilt3(S_4D(:,:,:,ii).*mask_odf,FWHM/voxelsize/2.355);
end
S_ind=reshape(S_4D,[N1*N2*N3,N4]);
S_voi=S_ind(I_mask,:);
clear S_4D ;
S_voi=S_voi';


%% bandpass filtering and normalzaition
S_voi = bandpass(S_voi,bf,1/TR);        
mean_voi=mean(S_voi);
S_voi=S_voi-repmat(mean_voi,N4,1);
std_voi=std(S_voi);
std_voi(std_voi==0)=1;
S_voi=S_voi./repmat(std_voi,N4,1); 
S_voi_1=S_voi;


S_4D=niftiread(f_fmri_2);
info=niftiinfo(f_fmri_2);
TR=info.PixelDimensions(4);
voxelsize=info.PixelDimensions(1);
S_4D=S_4D(:,:,:,1:N4);

[N1,N2,N3,N4]=size(S_4D); 
 S_4D=flip(S_4D,1);


for ii=1:N4
    S_4D(:,:,:,ii)=imgaussfilt3(S_4D(:,:,:,ii).*mask_odf,FWHM/voxelsize/2.355);
end

S_ind=reshape(S_4D,[N1*N2*N3,N4]);
S_voi=S_ind(I_mask,:);
clear S_4D;
S_voi=S_voi';



S_voi = bandpass(S_voi,bf,1/TR);        
mean_voi=mean(S_voi);
S_voi=S_voi-repmat(mean_voi,N4,1);
std_voi=std(S_voi);
std_voi(std_voi==0)=1;
S_voi=S_voi./repmat(std_voi,N4,1);


S_voi_1=[S_voi_1;S_voi];%% concatenate two runs (LR AND RL)
std_voi_1=std(S_voi_1);
std_voi_1(std_voi_1==0)=1;
S_voi_1=S_voi_1./repmat(std_voi_1,2*N4,1); 







%% read motor-task fMRI(nii format) data

f_fmri_1=f_fmri_total{3};
f_fmri_2=f_fmri_total{4};

S_4D=niftiread(f_fmri_1);
info=niftiinfo(f_fmri_1);
TR=info.PixelDimensions(4);
voxelsize=info.PixelDimensions(1);
S_4D=S_4D(:,:,:,1:N4);

[N1,N2,N3,N4]=size(S_4D); 
 S_4D=flip(S_4D,1);
%  S_3D= squeeze(S_4D(:,:,:,1));

for ii=1:N4
    S_4D(:,:,:,ii)=imgaussfilt3(S_4D(:,:,:,ii).*mask_odf,FWHM/voxelsize/2.355);
end

S_ind=reshape(S_4D,[N1*N2*N3,N4]);
S_voi=S_ind(I_mask,:);
clear S_4D;
S_voi=S_voi';
S_voi = bandpass(S_voi,bf,1/TR);        
mean_voi=mean(S_voi);
S_voi=S_voi-repmat(mean_voi,N4,1);
std_voi=std(S_voi);
std_voi(std_voi==0)=1;
S_voi=S_voi./repmat(std_voi,N4,1); 
S_voi_2=S_voi;


S_4D=niftiread(f_fmri_2);
info=niftiinfo(f_fmri_2);
TR=info.PixelDimensions(4);
voxelsize=info.PixelDimensions(1);
S_4D=S_4D(:,:,:,1:N4);

[N1,N2,N3,N4]=size(S_4D); 
 S_4D=flip(S_4D,1);
 S_3D= squeeze(S_4D(:,:,:,1));

for ii=1:N4
    S_4D(:,:,:,ii)=imgaussfilt3(S_4D(:,:,:,ii).*mask_odf,FWHM/voxelsize/2.355);
end

S_ind=reshape(S_4D,[N1*N2*N3,N4]);
S_voi=S_ind(I_mask,:);
clear S_4D;
S_voi=S_voi';


S_voi = bandpass(S_voi,bf,1/TR);        
mean_voi=mean(S_voi);
S_voi=S_voi-repmat(mean_voi,N4,1);
std_voi=std(S_voi);
std_voi(std_voi==0)=1;
S_voi=S_voi./repmat(std_voi,N4,1);


S_voi_2=[S_voi_2;S_voi];
std_voi_2=std(S_voi_2);
std_voi_2(std_voi_2==0)=1;
S_voi_2=S_voi_2./repmat(std_voi_2,2*N4,1); 


