%% Run this script after manual reorientation of structural from current subject folder
% nifti image needed. 
% Change log:
% 04/27/2018: Adapted for spm12

addpath('C:\Program Files\spm12')

cd('C:\data')
disp('Select subject folder on C:\data');
d = uigetdir;
load([d filesep 'config.mat']); % needed to load allrois

cd(dicom_dir)

spm_path = fileparts(which('spm'));

%% copy masks to dicom_dir for later normalization
for i=1:2 % copy masks to dicom_dir for later normalization
    copyfile(allrois{i}, dicom_dir);
    [~,fn, fne] = fileparts(allrois{i});
    myrois{i} = [fn fne];
end

matlabbatch{1}.spm.spatial.preproc.channel.vols(1) = cellstr(spm_select('FPList','.','^s.*\.nii$')); 
matlabbatch{1}.spm.spatial.preproc.channel.biasreg = 0.001; 
matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = 60; 
matlabbatch{1}.spm.spatial.preproc.channel.write = [0 1]; 
matlabbatch{1}.spm.spatial.preproc.tissue(1).tpm = {[spm_path filesep 'tpm' filesep 'TPM.nii,1']}; 
matlabbatch{1}.spm.spatial.preproc.tissue(1).ngaus = 1; 
matlabbatch{1}.spm.spatial.preproc.tissue(1).native = [1 0]; 
matlabbatch{1}.spm.spatial.preproc.tissue(1).warped = [0 0]; 
matlabbatch{1}.spm.spatial.preproc.tissue(2).tpm = {[spm_path filesep 'tpm' filesep 'TPM.nii,2']}; 
matlabbatch{1}.spm.spatial.preproc.tissue(2).ngaus = 1; 
matlabbatch{1}.spm.spatial.preproc.tissue(2).native = [1 0]; 
matlabbatch{1}.spm.spatial.preproc.tissue(2).warped = [0 0]; 
matlabbatch{1}.spm.spatial.preproc.tissue(3).tpm = {[spm_path filesep 'tpm' filesep 'TPM.nii,3']}; 
matlabbatch{1}.spm.spatial.preproc.tissue(3).ngaus = 2; 
matlabbatch{1}.spm.spatial.preproc.tissue(3).native = [1 0]; 
matlabbatch{1}.spm.spatial.preproc.tissue(3).warped = [0 0]; 
matlabbatch{1}.spm.spatial.preproc.tissue(4).tpm = {[spm_path filesep 'tpm' filesep 'TPM.nii,4']}; 
matlabbatch{1}.spm.spatial.preproc.tissue(4).ngaus = 3; 
matlabbatch{1}.spm.spatial.preproc.tissue(4).native = [1 0]; 
matlabbatch{1}.spm.spatial.preproc.tissue(4).warped = [0 0]; 
matlabbatch{1}.spm.spatial.preproc.tissue(5).tpm = {[spm_path filesep 'tpm' filesep 'TPM.nii,5']}; 
matlabbatch{1}.spm.spatial.preproc.tissue(5).ngaus = 4; 
matlabbatch{1}.spm.spatial.preproc.tissue(5).native = [1 0]; 
matlabbatch{1}.spm.spatial.preproc.tissue(5).warped = [0 0]; 
matlabbatch{1}.spm.spatial.preproc.tissue(6).tpm = {[spm_path filesep 'tpm' filesep 'TPM.nii,6']}; 
matlabbatch{1}.spm.spatial.preproc.tissue(6).ngaus = 2; 
matlabbatch{1}.spm.spatial.preproc.tissue(6).native = [0 0]; 
matlabbatch{1}.spm.spatial.preproc.tissue(6).warped = [0 0]; 
matlabbatch{1}.spm.spatial.preproc.warp.mrf = 1; 
matlabbatch{1}.spm.spatial.preproc.warp.cleanup = 1; 
matlabbatch{1}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2]; 
matlabbatch{1}.spm.spatial.preproc.warp.affreg = 'mni'; 
matlabbatch{1}.spm.spatial.preproc.warp.fwhm = 0; 
matlabbatch{1}.spm.spatial.preproc.warp.samp = 3; 
matlabbatch{1}.spm.spatial.preproc.warp.write = [1 0]; % write deformation fields: [1 0] write inverse, [0 1] write forward (not needed)

segbatch = matlabbatch;

clear matlabbatch;
spm_figure('GetWin', 'Graphics');
spm_image('init', segbatch{1}.spm.spatial.preproc.channel.vols{1});

spm_jobman('run', segbatch);

iy_file = spm_select('FPList','.','^iy_s.*.nii$');

matlabbatch{1}.spm.spatial.normalise.write.subj.def = cellstr(iy_file); 
matlabbatch{1}.spm.spatial.normalise.write.subj.resample = cellstr(myrois)';
matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70 
                                                          78 76 85]; 
matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = [3 3 3]; 
matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 5; 
matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = 'w';

normmaskbatch = matlabbatch;
clear matlabbatch;

spm_jobman('run', normmaskbatch);

%% illustrate the result
spm_figure('GetWin', 'Graphics');
spm_clf;
H = spm_orthviews('Image', segbatch{1}.spm.spatial.preproc.channel.vols{1});
for ix=1:length(myrois)
    ims(ix) = cellstr(spm_select('FPList','.',['^w' myrois{ix}]));
end
for i=1:numel(ims)
    col = zeros(1,3);
    col(mod(i,3)+1) = 1;
    spm_orthviews('AddColouredImage', H, ims{i}, col);
end
spm_orthviews('Redraw');

%% copy the normalized ROIs into the "data" (outputdir) folder
for i=1:numel(ims)
    copyfile(ims{i}, outputdir);
end
csf_image = spm_select('FPList','.','^c3s.*\.nii$');
copyfile(csf_image,[outputdir filesep 'c3s_csf.nii'])
copyfile(segbatch{1}.spm.spatial.preproc.channel.vols{1},[outputdir filesep 's',subjid,'_structural.nii']) % here we just put the structural into the "data"-folder!

%% tell the user that now the functional experiment can start
Status_text.String = 'Starting structural procedure.. done';
Status_text.String = char(Status_text.String, 'Everything should be ready for functional run.');