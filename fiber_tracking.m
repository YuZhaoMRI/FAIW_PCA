

dsi_path = getenv('DSI_STUDIO_PATH');
f_src='T:\HCP_3T_64S\preprocessed_data_mni2acpc\100206\Diffusion\data.nii.gz.src.gz.odf.gqi.1.25.fib.gz';
f_reg='T:\HCP_3T_64S\preprocessed_data_mni2acpc\100206\Structural\Motor_mask.nii';
command = [dsi_path, ' --action=trk --source=', f_src, ' --parameter_id=c9A99193FF304353Fb803FcbF041b1643A08601dcaCDCC4C3Ec ', '--seed=',f_reg,' --export=tdi'];

[status, results] = system(command);
assert(status == 0, 'Failed!\n%s', results)



nii=load_untouch_nii('dattdi.nii.gz');
img_dti=nii.img;
MatrixUser

