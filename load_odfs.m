function [odfs, I_odfs, odf_vertices, odf_faces] = load_odfs(f_odfs, flip_y)
%LOAD_ODFS Load ODF data generated with DSI Studio.
%
% Syntax:  [odfs,I_odfs,odf_vertices,odf_faces] = load_odfs(f_odfs,flip_y)
%
% Inputs:
%    f_odfs - Path of MAT file containing diffusion ODFs for every position
%      in the mask.
%    flip_y - Logical specifying whether to flip the loaded data along the
%      Y-axis. This is required to match the voxel space of DSI Studio
%      (-X,-Y,Z) with that of the HCP dataset (-X,Y,Z).
%
% Outputs:
%    odfs - Array containing the sampled ODF values for every voxel.
%    I_odfs - Indices of voxels with ODFs within the brain volume.
%    odf_vertices - Array containing 3D coordinates of ODF sample points.
%    odf_faces - Array containing the set of ODF vertex triplets that
%      would allow plotting the ODF surfaces as triangle meshes.
%
% Examples:
%    [odfs, I_odfs, odf_vertices] = load_odfs('data/odfs.mat', true)
%      First two outputs are enough for creating WM graph.
%    [odfs, I_odfs, odf_vertices, odf_faces] = ...
%      load_odfs('data/odfs.mat', true)
%      ODF faces are required for plotting the ODFs.

% Author: David Abramian, adapted from Martin Larsson, 2017
% Department of Biomedical Engineering, Linköping University, Sweden
% email: david.abramian@liu.se
% May 2021; Last revision: 13-May-2021


% check input
assert(endsWith(f_odfs, '.mat', 'IgnoreCase', true), ...
    'Expected .mat file\n%s', f_odfs)
assert(exist(f_odfs,'file')==2, 'File does not exist\n%s', f_odfs)

% load ODF data
load(f_odfs, 'odfs', 'I_odfs', 'odf_vertices', 'odf_faces', 'dimension');

if nargin > 1 && flip_y
    % find voxel indices of flipped volume
    I_vox = zeros(dimension);
    I_vox(I_odfs) = 1:length(I_odfs);
    I_vox_inv = flip(I_vox, 2);
    I_odfs = find(I_vox_inv);
    
    % rearrange ODFs according to permutation
    perm = I_vox_inv(I_odfs);
    odfs = odfs(:, perm);
    
    % flip Y dimension for ODF vertices
    odf_vertices(2,:) = -odf_vertices(2, :);
end

