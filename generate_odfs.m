function f_odfs = generate_odfs(f_dwi, f_mask)
%GENERATE_ODFS Generate ODFs from diffusion data using DSI Studio.
%Generates diffusion ODFs from diffusion weighted imaging (DWI) data by
%using system calls to the DSI Studio software. The ODFs are saved to a MAT
%file for reuse. The path to the DSI Studio software needs to be set in the
%DSI_STUDIO_PATH environmental variable.
%
% Syntax:  f_odfs = generate_odfs(f_mask,f_diff,dsi_path)
%
% Inputs:
%    f_dwi - Path of NIfTI volume series containing DWI data (can be
%      compressed).
%    f_mask - Path of uncompressed NIfTI volume containing binary white
%      matter mask.
%
% Outputs:
%    f_odfs - Path of MAT file containing generated ODFs. The ODFs are
%      saved in file 'odfs.mat' in the directory of the DWI data.
%
% Example:
%    f_odfs = generate_odfs('data/dwi.nii.gz', 'data/mask.nii')

% Author: David Abramian
% Department of Biomedical Engineering, Linköping University, Sweden
% email: david.abramian@liu.se
% May 2021; Last revision: 13-May-2021


% check input
dsi_path = getenv('DSI_STUDIO_PATH');
assert(~isempty(dsi_path), 'DSI_STUDIO_PATH environmental variable is not set')
assert(endsWith(f_dwi, [".nii", ".nii.gz"], 'IgnoreCase', true), ...
    'Expected .nii or .nii.gz file\n%s', f_dwi)
assert(endsWith(f_mask, '.nii', 'IgnoreCase', true), ...
    'Expected .nii file\n%s', f_mask)
assert(exist(f_dwi, 'file')==2, 'File does not exist\n%s', f_dwi)
assert(exist(f_mask, 'file')==2, 'File does not exist\n%s', f_mask)
assert(exist(dsi_path, 'file')==2, 'File does not exist\n%s', dsi_path)

f_dir=dir(f_dwi);
f_bval=[f_dir.folder '\bvals'];
% f_bvec=[f_dir.folder '\bvecs_MNI152'];
f_bvec=[f_dir.folder '\bvecs'];
assert(exist(f_bval, 'file')==2, 'File does not exist\n%s', f_bval)
assert(exist(f_bvec, 'file')==2, 'File does not exist\n%s', f_bvec)



f_path = fileparts(f_dwi);
% f_src = fullfile(f_path, 'odfs.src.gz');
f_src = [f_dwi '_BW.src.gz'];
% f_fib = fullfile(f_path, 'odfs.src.gz.odf8.f5rec.bal.012fy.rdi.gqi.1.25.fib.gz');
f_fib = [f_dwi '_BW.src.gz.odf.gqi.1.25.fib.gz'];
f_odfs = fullfile(f_path, 'odfs.mat');

if ~exist(f_odfs,'file')
    % generate .src file
    if ~exist(f_src, 'file')

        fprintf('DSI Studio: generating source file\n')
        command = [dsi_path, ' --action=src --source=', f_dwi, ' --bval=',f_bval ' --bvec=',f_bvec,' --output=', f_src];
        [status, results] = system(command);
        assert(status == 0, 'Failed to create .src file:\n%s', results)
    end
    
    % generate .fib file
    if ~exist(f_fib, 'file')
        fprintf('DSI Studio: constructing ODFs\n')
        command = [dsi_path, ' --action=rec --source=', f_src, ' --thread=2 --mask=', f_mask, ' --method=4 --param0=1.25 --odf_order=8 --num_fiber=5 --check_btable=1 --record_odf=1'];
        [status, results] = system(command);
        assert(status == 0, 'Failed to create .fib file:\n%s', results)
    end
    
    % extract .fib file and rename to .mat file
    gunzip(f_fib);
   
    movefile(f_fib(1:end-3), f_odfs)
    
    % process ODF file
    fprintf('Processing ODFs\n')
    mask_size=size(niftiread(f_mask));
    process_odfs(f_odfs,mask_size)
    
    % delete intermediate files
%     delete(f_src)  %zhaoyu
%     delete(f_fib)  %zhaoyu
else
    fprintf('ODF file already exists\n')
end

