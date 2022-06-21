function dss_filter_fmri(f_out, f_fmri_total, mask, G, tau, Seg_p)
%DSS_FILTER_FMRI Perform graph filtering on fMRI volume series.
%Smooth fMRI volume series using graph heat kernel filtering.
%
% Syntax:  dss_filter_fmri(f_fmri,f_mask,G,tau)
%
% Inputs:
%    f_out - Path of output filtered fMRI volume series. Data will be
%      stored in an uncompressed NIfTI file.
%    f_fmri - Path of uncompressed NIfTI volume series containing fMRI
%      data to filter.
%    f_mask - Path of uncompressed NIfTI volume containing binary white
%      matter mask.
%    G - White matter graph structure.
%    tau - Size parameter for filter kernel. Larger values result in larger
%      filters. Takes real non-negative values.
%
% Examples:
%    graph_filter_fmri('data/filtered.nii', data/fmri.nii', ...
%      'data/mask.nii', G, 4)
%
% See also: DSS_CREATE_GRAPH, GRAPH_FILTERING

% Author: David Abramian, adapted from Martin Larsson, 2017
% Department of Biomedical Engineering, Linköping University, Sweden
% email: david.abramian@liu.se
% May 2021; Last revision: 13-May-2021


% % check input
% assert(endsWith(f_fmri, '.nii', 'IgnoreCase', true), ...
%     'Expected .nii file\n%s', f_fmri)
% % assert(endsWith(f_mask, '.nii', 'IgnoreCase', true), ...
% %     'Expected .nii file\n%s', f_mask)
% assert(exist(f_fmri,'file')==2, 'File does not exist\n%s', f_fmri)
% assert(exist(f_mask,'file')==2, 'File does not exist\n%s', f_mask)
assert(all(isfield(G, ["L","l_max"])), 'Missing fields in graph structure')
assert(tau >= 0, 'tau must be a non-negative scalar')

% load white matter mask

dim = size(mask);
I_mask = find(mask);

% construct heat kernel filter
% kernel = cell(length(tau), 1);
% for i = 1:length(tau)
%     kernel{i} = @(x) exp(-x*tau(i));
% end
kernel = @(x) exp(-x*tau);
cheb_ord = 15;

% load volume headers
% headers = ml_load_nifti(f_fmri);
% n_vols = length(headers);

% % load fMRI volumes
% signals = zeros(length(I_mask), n_vols);
% for i = 1:n_vols
%     progresss(i, n_vols, 'Loading fMRI volumes... ')
%     
%     [h, vol] = ml_load_nifti(f_fmri, 1);
%     signals(:,i) = vol(I_mask);
% end
% progresss(i+1, n_vols, 'Loading fMRI volumes... ')
FWHM=3.0;
f_fmri=f_fmri_total{1};
tic
S_4D=niftiread(f_fmri);
toc
info=niftiinfo(f_fmri);
TR=info.PixelDimensions(4);
voxelsize=info.PixelDimensions(1);

fprintf('Read and smooth data with FWHM=%dmm\n',FWHM);
S_4D_t=S_4D(:,:,:,1:284);
S_4D=S_4D_t;
[N1,N2,N3,N4]=size(S_4D); 
for ii=1:N4
    S_4D(:,:,:,ii)=imgaussfilt3(S_4D(:,:,:,ii).*mask,FWHM/voxelsize/2.355);
end


[N1,N2,N3,N4]=size(S_4D);
S_ind=reshape(S_4D,[N1*N2*N3,N4]);
S_voi=S_ind(I_mask,:);
clear S_4D S_ind S_4D_t;
S_voi=S_voi';
S_voi = bandpass(S_voi,[0.01 0.12],1/TR);        
mean_voi=mean(S_voi);
S_voi=S_voi-repmat(mean_voi,N4,1);
std_voi=std(S_voi);
std_voi(std_voi==0)=1;
S_voi_1=S_voi./repmat(std_voi,N4,1); 
clear S_voi;



f_fmri=f_fmri_total{2};
tic
S_4D=niftiread(f_fmri);
toc
info=niftiinfo(f_fmri);
TR=info.PixelDimensions(4);
voxelsize=info.PixelDimensions(1);

fprintf('Read and smooth data with FWHM=%dmm\n',FWHM);
S_4D_t=S_4D(:,:,:,1:284);
S_4D=S_4D_t;
[N1,N2,N3,N4]=size(S_4D); 
for ii=1:N4
    S_4D(:,:,:,ii)=imgaussfilt3(S_4D(:,:,:,ii).*mask,FWHM/voxelsize/2.355);
end


[N1,N2,N3,N4]=size(S_4D);
S_ind=reshape(S_4D,[N1*N2*N3,N4]);
S_voi=S_ind(I_mask,:);
clear S_4D S_ind S_4D_t;
S_voi=S_voi';
S_voi = bandpass(S_voi,[0.01 0.12],1/TR);        
mean_voi=mean(S_voi);
S_voi=S_voi-repmat(mean_voi,N4,1);
std_voi=std(S_voi);
std_voi(std_voi==0)=1;
S_voi_2=S_voi./repmat(std_voi,N4,1); 
clear S_voi;

% f_fin=[f_fmri(1:end-4) '_fin'];
% load(f_fin);
% 
% fin=tc_fin*S_voi;
% 
% 
% mean_voi=mean(fin);
% fin=fin-repmat(mean_voi,size(fin,1),1);
% std_voi=std(fin);
% std_voi(std_voi==0)=1;
% fin=fin./repmat(std_voi,size(fin,1),1); 




% perform graph heat kernel filtering
fprintf('Performing calculation of HOMO:\n')







% vn=size(S_n,2);

% for ii=1:N1
%     slice_index=ii
%     for jj=1:N2
%         for kk=1:N3
%             key=mask(ii,jj,kk);
%             if key==1 
% 
%                 tem=zeros(size(mask));
%     %              tem(66,122,61)=1;
%                 tem(ii,jj,kk)=1;
%                 sig=tem(I_mask);
%                 coeff = graph_filtering(G.L, G.l_max, kernel, sig, cheb_ord);
%                 tem(I_mask)=coeff;
%                 [w_sort,ind_sort]=sort(coeff,'descend');
%                 vox_n=64;
%                 nei_ind=ind_sort(1:vox_n);
%                 nei_w=w_sort(1:vox_n);
%                 sum_w=sum(nei_w(2:end));
%                 if sum_w<0.0001
%                    continue 
%                 end
% 
%                 nei_w=nei_w/sum_w;
%                 nei_w(1)=0.1;
%                 nei_w=sqrt(nei_w);
%                 mat_w=repmat(nei_w',size(S_voi,1),1); 
%     %            mat_w=repmat(nei_w',size(fin,1),1); 
%                 S_n=S_voi(:,nei_ind);
% 
%     %           S_n=fin(:,nei_ind);
%                 S_n=S_n.*mat_w;
%     %             S_n(:,1)=S_n(:,1)*4;
% 
% 
%                 cor=S_n'*S_n;
%                 A=svd(cor);
%     %             A=svd(S_n);        
%                 HOMO(ii,jj,kk)=A(1);
% 
%             end 
%         end
%     end
% end


% mask(:,:,[1:68 70:end])=0;
% mask(1:90,:,:)=0;
% mask(:,1:90,:)=0;

tic
% delete(gcp('nocreate'));
% MyPar=parpool(8);
% delete(MyPar);
ind_mask=find(mask);


if Seg_p(2)<Seg_p(1)
   R=((Seg_p(2)-1)*ceil(length(ind_mask)/Seg_p(1))+1):1:Seg_p(2)*ceil(length(ind_mask)/Seg_p(1));
else
   R=((Seg_p(2)-1)*ceil(length(ind_mask)/Seg_p(1))+1):1:length(ind_mask);   
end
HOMO_1=zeros(N1,N2,N3);
HOMO_2=zeros(N1,N2,N3);
HOMO_id_1=zeros(length(ind_mask),1);
HOMO_id_2=zeros(length(ind_mask),1);
L=G.L;
l_max=G.l_max;


assert(size(L,2) == size(L,1), 'A is not square.');
arange = [0 l_max];
% approximate kernel with Chebyshev polynomial
c = sgwt_cheby_coeff(kernel, cheb_ord, cheb_ord+1, arange);



for ind=R
    [ii,jj,kk]=ind2sub([N1,N2,N3],ind_mask(ind));
                tem=zeros(N1,N2,N3);
%                tem(50,78,41)=1;
                tem(ii,jj,kk)=1;
                sig=tem(I_mask);
                % coeff = graph_filtering(L, l_max, kernel, sig, cheb_ord);
                
                % decompose signal and return low-pass component
                tic
                coeff = sgwt_cheby_op(sig, L, c, arange);
                toc
%                  tem(I_mask)=coeff;
                
                tic
                vox_n=64;
                x=single(coeff);
                nei_ind=zeros(vox_n,1);
                nei_w=zeros(vox_n,1);
                for ii=1:vox_n
                    [nei_w(ii), nei_ind(ii)] = max(x);
                    x(nei_ind(ii)) = 0;
                end


                sum_w=sum(nei_w(2:end));
                if sum_w<0.0001
                   continue 
                end
                fra_c=0.2;
                nei_w=(1-fra_c)*nei_w/sum_w;
                nei_w(1)=fra_c;
                nei_w=sqrt(nei_w);
                mat_w=repmat(nei_w',size(S_voi_1,1),1); 

                S_n=S_voi_1(:,nei_ind);
                S_n=S_n.*mat_w;
                cor=S_n'*S_n/size(S_voi_1,1);
                A=svd(cor);        
                HOMO_id_1(ind)=A(1);
                
                S_n=S_voi_2(:,nei_ind);
                S_n=S_n.*mat_w;
                cor=S_n'*S_n/size(S_voi_2,1);
                A=svd(cor);        
                HOMO_id_2(ind)=A(1);
                toc
                h=1;
                if mod(ind,200)==0
                progresss(ind, R(end), '... ')
                end
end

HOMO_1(ind_mask)=HOMO_id_1;
HOMO_2(ind_mask)=HOMO_id_2;





f_fmri=f_fmri_total{1};
f_homo=[f_fmri(1:end-4) '_homo2_' num2str(Seg_p(2)) '.mat'];
save(f_homo,'HOMO_1');
f_fmri=f_fmri_total{2};
f_homo=[f_fmri(1:end-4) '_homo2_' num2str(Seg_p(2)) '.mat'];
save(f_homo,'HOMO_2');

toc

% figure;
% imshow(squeeze(log(1+abs(((HOMO1(:,:,69)-200)/200)))));
% imcontrast
% 
% 
% 
% 
% coeff = graph_filtering(G.L, G.l_max, kernel, signals, cheb_ord);
% mask_odf(75,52,61)=1000;
% figure;imshow(squeeze(mask_odf(:,:,61)))
% imcontrast
% 
% % save filtered volumes
% for i = 1:n_vols
%     progresss(i, n_vols, 'Saving filtered volumes... ')
%     
%     vol_out = nan(dim);
%     vol_out(I_mask) = coeff(:,i);
%     
%     h.fname = f_out;
%     h.n(1) = i;
%     spm_write_vol(h, vol_out);
% end
% progresss(i+1, n_vols, 'Saving filtered volumes... ')


end

