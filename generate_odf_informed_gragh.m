function generate_odf_informed_gragh(f_mask, f_dwi, alpha, neigh)


% check input
% assert(nargin == 7, 'Expected 8 arguments.')
assert(endsWith(f_mask, '.nii', 'IgnoreCase', true), ...
    'Expected .nii file\n%s', f_mask)
assert(endsWith(f_dwi, [".nii", ".nii.gz"], 'IgnoreCase', true), ...
    'Expected .nii or .nii.gz file\n%s', f_dwi)
assert(exist(f_mask,'file')==2, 'File does not exist\n%s', f_mask)
assert(exist(f_dwi,'file')==2, 'File does not exist\n%s', f_dwi)



f_path = fileparts(f_dwi);
f_graph = fullfile(f_path, ['graph_a', num2str(100*alpha, '%.0f'), '_n', num2str(neigh), '_WB.mat']);

if ~exist(f_graph, 'file')
    f_odfs = fullfile(f_path, 'odfs.mat');
    if ~exist(f_odfs, 'file')
    % generate ODFs
        f_odfs = generate_odfs(f_dwi, f_mask);
    end
    
    % create white matter graph
    G = dss_create_graph(f_odfs, f_mask, alpha, neigh);
else
    % load existing white matter graph
    load(f_graph, 'G')
end



