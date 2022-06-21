function G = dss_create_graph(f_odfs, f_mask, alpha, neigh)
%DSS_CREATE_GRAPH Create graph strucutre from WM mask and ODFs.
%Creates a structure containing the white matter graph adjacency and
%Laplacian matrices from diffusion ODFs. The graph strucutre is saved to
%a MAT file for reuse.
%
% Syntax:  G = dss_create_graph(f_odfs,f_mask,alpha,neigh)
%
% Inputs:
%    f_odfs - Path of MAT file containing diffusion ODFs for every position
%      in the mask.
%    f_mask - Path of uncompressed NIfTI volume containing binary white
%      matter mask.
%    alpha - Soft threshold for graph weights. Higher values result in more
%      anisotropic filters. Takes values in [0,1] (default: 0.9).
%    neigh - Size of neighborhood shell defining the set of neighbors for
%      each vertex. Takes values in {3,5} (default: 5).
%
% Outputs:
%    G - Structure containing white matter graph. Includes the adjacency
%      and Laplacian matrices, as well as the voxel indices of the graph
%      vertices.
%
% Examples:
%    f_odfs = dss_create_graph('data/odfs.mat', 'data/mask.nii')
%    f_odfs = dss_create_graph('data/odfs.mat', 'data/mask.nii', 0.8, 3)
%
% See also: DSS_ADJACENCY_MATRIX,  SGWT_LAPLACIAN

% Author: David Abramian, adapted from Martin Larsson, 2017
% Department of Biomedical Engineering, Linköping University, Sweden
% email: david.abramian@liu.se
% May 2021; Last revision: 13-May-2021


% check input
assert(endsWith(f_odfs, '.mat', 'IgnoreCase', true), ...
    'Expected .mat file\n%s', f_odfs)
assert(endsWith(f_mask, '.nii', 'IgnoreCase', true), ...
    'Expected .nii file\n%s', f_mask)
assert(exist(f_odfs,'file')==2, 'File does not exist\n%s', f_odfs)
assert(exist(f_mask,'file')==2, 'File does not exist\n%s', f_mask)
if exist('alpha', 'var') && ~isempty(alpha)
    assert(0 <= alpha && alpha <= 1, 'alpha takes values in [0,1]')
else
    alpha = 0.9;  % default
end
if exist('neigh', 'var') && ~isempty(neigh)
    assert(any(neigh == [3, 5]), 'neigh takes values in {3,5}')
else
    neigh = 5;  % default
end

f_path = fileparts(f_odfs);
f_graph = fullfile(f_path, ['graph_a', num2str(100*alpha, '%.0f'), '_n', num2str(neigh), '_WB.mat']);

if exist(f_graph, 'file')
    % load existing white matter graph
    fprintf('Loading existing WM graph\n')
    load(f_graph, 'G');
else
    % load white matter mask
%     [~, mask] = ml_load_nifti(f_mask);
%     nii_mask=load_nii(f_mask);% zhaoyu
%     mask= nii_mask.img;         %zhaoyu
f_odfs_dir=dir(f_odfs);
mask_dir=[f_odfs_dir.folder '\mask_odf'];
load(mask_dir);
mask=mask_odf;


%     mask = logical(mask);
    I_mask = find(mask);
    
    fprintf('Creating WM graph and Laplacian\n')
    
    % create graph structure
    G = struct();
    G.alpha = alpha;
    G.neigh = neigh;
    G.dim = size(mask);
    G.N = length(I_mask);
    G.indices = I_mask;
    
    % create graph adjacency matrix
    G.A = dss_adjacency_matrix(f_odfs, mask, alpha, neigh);
    
    % create graph Laplacian matrix
    [G.L, G.W]= sgwt_laplacian(G.A,'opt','normalized');
    G.l_max = min(2,eigs(G.L,1,'lm',struct('tol',5e-3,'p',10,'disp',0))*1.02);
    
    % save graph structure
    save(f_graph, 'G','-v7.3');
end

