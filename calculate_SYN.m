function calculate_SYN(root,f_fmri_total, f_dwi, f_mask,par)
% Fiber architecture informed synchrony mapping
% Syntax:  calculate_SYN(root,f_fmri_total, f_dwi, f_nodif_mask, par);

%% parameters 
tau=par.tau;
num_neigh=par.num_neigh;
FWHM=par.FWHM;
bf=par.band;
N4=par.num_tc;
neigh_graph=par.neigh_graph;
alpha=par.alpha;

%% check input
dwi_path = fileparts(f_dwi);
f_graph = fullfile(dwi_path, ['graph_a', num2str(100*alpha, '%.0f'), '_n', num2str(neigh_graph), '_WB.mat']);
assert(exist(f_graph,'file')==2, 'File does not exist\n%s', f_graph)
load(f_graph, 'G')
mask_path=fullfile(dwi_path, 'mask_odf');
load(mask_path,'mask_odf');
dim = G.dim;
I_mask = find(mask_odf);
mask_odf=single(mask_odf);




%% the 2rd step of the preprocess fMRI data:smoothing,spatial filtering, and normalizing

f_fmri_1=f_fmri_total{1};
f_path = fileparts(f_fmri_1);
f_data = [f_path,'\FWHM',num2str(FWHM), '_s_voi_WB.mat'];

if ~exist(f_data, 'file')
   [S_voi_1,S_voi_2]=preprocess_fMRI_nii(f_fmri_total,mask_odf,I_mask,FWHM,bf,N4);
    save(f_data,'S_voi_1','S_voi_2','-v7.3');   
else
    load(f_data);

end



fprintf('Performing syn mapping')


SYN_1=zeros(dim);
SYN_2=zeros(dim);
SYN_id_1=zeros(length(I_mask),1);
SYN_id_2=zeros(length(I_mask),1);


W_t=G.W;
n_vox = length(I_mask);
I_vox = zeros(dim);
I_vox(I_mask) = 1:n_vox;
I_vox=single(I_vox);
tem=single(zeros(dim));
vox_n=num_neigh;
TN=2*N4;
N_w=5;%% the size of volume larly determins the speed of calculation 




for ind=1:1:length(I_mask)
    [ii,jj,kk]=ind2sub(dim,I_mask(ind));
    if (kk-N_w>0)&&kk==65

        tem=tem*0;
        tem(ii-N_w:ii+N_w,jj-N_w:jj+N_w,kk-N_w:kk+N_w)=1;
        I_V=find(tem.*mask_odf>0);
        I_V2=I_vox(I_V);
        W=single(full(W_t(I_V2,I_V2)));

        D=diag(sum(W,1));
        L = W-D; %


%% simulating diffusion

        tem=tem*0;
        tem(ii,jj,kk)=1;
        s=tem(I_V);
%         [~, m_ind] = max(s);   

        deta=0.01;  % step of iteration    
        dL=L.*deta;
        N_iter=ceil(tau/deta);
         
        
        for iter=1:1:N_iter
            s=s+dL*s;
        end
        


%% The simulation of the virtual diffusion can be implemented in the following codes
%% But slow!!
%     [u, v]=eig(L,'balance','vector');
%     tem=tem*0;
%     tem(ii,jj,kk)=1;
%     s=tem(I_V);
%     flt =exp(-tau*(v));
%     sf=u*(flt.*(u'*s));
%     tem(I_V)=sf;
%     


    




%% prepare spatial windows
    nei_ind=zeros(vox_n,1);
    nei_w=zeros(vox_n,1);
    for i_n=1:vox_n
        [nei_w(i_n), tem_ind] = max(s);
        s(tem_ind) = 0;
        nei_ind(i_n)=I_V(tem_ind);
    end
    nei_ind=I_vox(nei_ind);
    sum_w=sum(nei_w);
    nei_w=nei_w/sum_w;
    nei_w=sqrt(nei_w);
    mat_w=repmat(nei_w',TN,1); 


%% calcualte SYN1
    S_n=S_voi_1(:,nei_ind);
    S_n=S_n.*mat_w;
    cor=S_n'*S_n;
    A=svd(cor);        
    SYN_id_1(ind)=max(A);
    

    

 % calcualte SYN2
    S_n=S_voi_2(:,nei_ind);
    S_n=S_n.*mat_w;
    cor=S_n'*S_n;
    A=svd(cor);        
    SYN_id_2(ind)=max(A);


    end
end





SYN_1(I_mask)=SYN_id_1/(TN-1);
SYN_2(I_mask)=SYN_id_2/(TN-1);



%% svae SYN 1
data_nii=load_nii(f_mask);
data_nii.img=flip(SYN_1,1);
f_fmri=f_fmri_total{1};
f_SYN_1=[root, '\syn_tau' num2str(tau*10)  'ACPC_' num2str(N4),'.nii'];
save_nii(data_nii,f_SYN_1);



%% save SYN 2
data_nii=load_nii(f_mask);
data_nii.img=flip(SYN_2,1);
f_fmri=f_fmri_total{3};
f_SYN_2=[root,'\syn_tau' num2str(tau*10)  'ACPC_' num2str(N4),'.nii'];
save_nii(data_nii,f_SYN_2);





end

