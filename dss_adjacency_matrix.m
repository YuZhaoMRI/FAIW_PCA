function A = dss_adjacency_matrix(f_odfs, mask, alpha, neigh)
%DSS_ADJACENCY_MATRIX Create white matter graph adjacency matrix.
%Creates white matter grpah adjacency matrix from white matter mask and
%diffusion ODFs. The adjacency mask provides a full description of the
%graph and encodes the white matter microstructure.
%
% Syntax:  A = dss_adjacency_matrix(f_odfs,mask,alpha,neigh)
%
% Inputs:
%    f_odfs - Path of MAT file containing diffusion ODFs for every position
%      in the mask.
%    mask - 3D volume of binary white matter mask.
%    alpha - Soft threshold for graph weights. Higher values result in more
%      anisotropic filters. Takes values in [0,1] (default: 0.9).
%    neigh - Size of neighborhood shell defining the set of neighbors for
%      each vertex. Takes values in {3,5} (default: 5).
%
% Outputs:
%    A - Sparse adjacency matrix of white matter graph.
%
% Examples:
%    A = dss_adjacency_matrix('data/odfs.mat', mask)
%    A = dss_adjacency_matrix('data/odfs.mat', mask, 0.8, 3)
%
% See also: DSS_CREATE_GRAPH


% Author: David Abramian, adapted from Martin Larsson, 2017
% Department of Biomedical Engineering, Linköping University, Sweden
% email: david.abramian@liu.se
% May 2021; Last revision: 13-May-2021


% check input
assert(endsWith(f_odfs, '.mat', 'IgnoreCase', true), ...
    'Expected .mat file\n%s', f_odfs)
assert(exist(f_odfs,'file')==2, 'File does not exist\n%s', f_odfs)
if exist('alpha', 'var') && ~isempty(alpha)
    assert(0 <= alpha && alpha <= 1, 'alpha takes values in [0,1]')
else
    alpha = 0.90;  % default:0.9
end
if exist('neigh', 'var') && ~isempty(neigh)
    assert(any(neigh == [3, 5]), 'neigh takes values in {3,5}')
else
    neigh = 5;  % default:5
end

dim = size(mask);
I_mask = find(mask);

% load ODFs
[odfs, I_odfs, odf_vertices] = load_odfs(f_odfs, true);
odfs = [odfs; odfs];

% I_mask=I_odfs;

% check correspondence of mask and ODFs
assert(length(I_mask) == size(odfs, 2), ...
    'Number of voxels in mask does not match number of ODFs')
assert(isequal(I_mask, I_odfs), ...
    'Mask and ODFs are misaligned')

n_vox = length(I_mask);
I_vox = zeros(dim);
I_vox(I_mask) = 1:n_vox;

% define neighborhood and direction of edges
[dirs, norm_dirs, ~, n_par] = ml_create_neighborhood(neigh, true);
neighborhood = size(dirs, 2);

% find all neighboring voxels
ci = repelem(I_mask, neighborhood, 1);
cs = zeros(size(ci, 1), 3);
[cs(:,1), cs(:,2), cs(:,3)] = ind2sub(dim, ci);
ns = cs + repmat(dirs', n_vox, 1);
ni = sub2ind(dim, ns(:,1), ns(:,2), ns(:,3));
ni = reshape(ni, neighborhood, n_vox);
ci = I_vox(ci);
ni = I_vox(ni);
lni = logical(ni);

% calculate Pdiff from ODF spherical harmonics
omega = 4*pi/n_par;  % we still want 98 for 5x5x5 neighborhood
[odfs_sum, Nt] = hb_soliangodf(odfs, odf_vertices, norm_dirs, omega);
p_diff = (odfs_sum/Nt) .* lni;
p_diff = 0.5*p_diff ./ repmat(max(p_diff),neighborhood,1);

% create adjacency matrix from Pdiff
A = sparse(ci(lni), ni(lni), p_diff(lni), n_vox, n_vox);
A = A+A';

% apply soft threshold on edge weights
I_adj = full(find(A));
edge_weights = full(A(I_adj));

mu = @(x,a,b) (x*(1-a)).^b ./ ((x*(1-a)).^b + (a*(1-x)).^b);  % tunable sigmoid function

weights = mu(edge_weights, alpha, 50);
edge_weights = edge_weights .* weights;

% rebuild adjacency matrix
[i,j] = ind2sub([n_vox, n_vox], I_adj);
A = sparse(i, j, edge_weights);


