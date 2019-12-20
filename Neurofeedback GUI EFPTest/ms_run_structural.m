function ms_run_structural(series, out_dir, allrois, dicom_dir)
% this function runs all the structural stuff

% change log:
% 2017/12/04: Saving of output files (GM, WM, bias corrected) during segmentation step disabled by setting last parameter in vector array to 0
% 2018/03/27: Changed for compatibility with spm12; consider to delete spm8 code
% 2018/03/27: Omit normalization of structural; consider to delete this code
% 2018/03/27: CSF-image

% this is to update the gui
Status_text = findobj('Tag','Status_text');
Status_text.String = 'Starting structural procedure..';
%% load the settings
load([out_dir filesep 'config.mat']);

%% copy masks to dicom_dir for later normalization
for i=1:2 % copy masks to dicom_dir for later normalization
    copyfile(allrois{i}, dicom_dir);
    [~,fn, fne] = fileparts(allrois{i});
    myrois{i} = [fn fne];
end

%% all structural stuff: dicom import, segmentation, nomralization etc.
cd(dicom_dir) % not sure if all the stuff should be saved here or in "data"
filter_regex = ['^001_', sprintf('%06d',series), '_.*'];
matlabbatch{1}.spm.util.dicom.data = cellstr(spm_select('FPList',dicom_dir,filter_regex));
matlabbatch{1}.spm.util.dicom.root = 'flat';
matlabbatch{1}.spm.util.dicom.outdir = {'.'};
matlabbatch{1}.spm.util.dicom.convopts.format = 'nii';
matlabbatch{1}.spm.util.dicom.convopts.icedims = 0;

dcmbatch = matlabbatch; 
clear matlabbatch;
spm_jobman('initcfg'); 
spm_jobman('run', dcmbatch);

% segment
% matlabbatch{1}.spm.spatial.preproc.data = cellstr(spm_select('FPList','.','^s.*\.nii$')); 
% matlabbatch{1}.spm.spatial.preproc.output.GM = [0 0 0];
% matlabbatch{1}.spm.spatial.preproc.output.WM = [0 0 0];
% matlabbatch{1}.spm.spatial.preproc.output.CSF = [0 0 0];
% matlabbatch{1}.spm.spatial.preproc.output.biascor = 0;
% matlabbatch{1}.spm.spatial.preproc.output.cleanup = 0;
% spmdir = spm('Dir');
% 
% matlabbatch{1}.spm.spatial.preproc.opts.tpm = {
%                                                [spmdir '/tpm/grey.nii']
%                                                [spmdir '/tpm/white.nii']
%                                                [spmdir '/tpm/csf.nii']
%                                                };
% matlabbatch{1}.spm.spatial.preproc.opts.ngaus = [2
%                                                  2
%                                                  2
%                                                  4];
% matlabbatch{1}.spm.spatial.preproc.opts.regtype = 'mni';
% matlabbatch{1}.spm.spatial.preproc.opts.warpreg = 1;
% matlabbatch{1}.spm.spatial.preproc.opts.warpco = 25;
% matlabbatch{1}.spm.spatial.preproc.opts.biasreg = 0.0001;
% matlabbatch{1}.spm.spatial.preproc.opts.biasfwhm = 60;
% matlabbatch{1}.spm.spatial.preproc.opts.samp = 3;
% matlabbatch{1}.spm.spatial.preproc.opts.msk = {''};

% spm12:
spm_path = fileparts(which('spm'));

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
% spm_image('init', segbatch{1}.spm.spatial.preproc.data{1});
spm_image('init', segbatch{1}.spm.spatial.preproc.channel.vols{1});

spm_jobman('run', segbatch);

% normalize structural to mni % NOT NEEDED !?
% snfile = spm_select('FPList','.','^s.*_seg_sn\.mat$');
% y_file = spm_select('FPList','.','^y_s.*.nii$');
% 
% % matlabbatch{1}.spm.spatial.normalise.write.subj.matname = cellstr(snfile);
% % matlabbatch{1}.spm.spatial.normalise.write.subj.resample = segbatch{1}.spm.spatial.preproc.data;
% % 
% % matlabbatch{1}.spm.spatial.normalise.write.roptions.preserve = 0;
% % matlabbatch{1}.spm.spatial.normalise.write.roptions.bb = [-78 -112 -50
% %                                                           78 76 85];
% % matlabbatch{1}.spm.spatial.normalise.write.roptions.vox = [1 1 1];
% % matlabbatch{1}.spm.spatial.normalise.write.roptions.interp = 0;
% % matlabbatch{1}.spm.spatial.normalise.write.roptions.wrap = [0 0 0];
% % matlabbatch{1}.spm.spatial.normalise.write.roptions.prefix = 'w';
% 
% matlabbatch{1}.spm.spatial.normalise.write.subj.def = cellstr(y_file);
% matlabbatch{1}.spm.spatial.normalise.write.subj.resample = segbatch{1}.spm.spatial.preproc.channel.vols;
% matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70 
%                                                           78 76 85]; 
% matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = [1 1 1]; 
% matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 0; 
% matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = 'w';
% 
% normmaskbatch = matlabbatch;
% clear matlabbatch;
% spm_jobman('run', normmaskbatch);

% next move masks to native space
% snfile = spm_select('FPList','.','^s.*_seg_inv_sn\.mat$');
 
% matlabbatch{1}.spm.spatial.normalise.write.subj.matname = cellstr(snfile);
% matlabbatch{1}.spm.spatial.normalise.write.subj.resample = myrois; 
% matlabbatch{1}.spm.spatial.normalise.write.roptions.preserve = 0;
% matlabbatch{1}.spm.spatial.normalise.write.roptions.bb = [-78 -112 -50
%                                                           78 76 85];
% matlabbatch{1}.spm.spatial.normalise.write.roptions.vox = [3 3 3];
% matlabbatch{1}.spm.spatial.normalise.write.roptions.interp = 0;
% matlabbatch{1}.spm.spatial.normalise.write.roptions.wrap = [0 0 0];
% matlabbatch{1}.spm.spatial.normalise.write.roptions.prefix = 'w';

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
    copyfile(ims{i}, out_dir);
end
csf_image = spm_select('FPList','.','^c3s.*\.nii$');
copyfile(csf_image,[out_dir filesep 'c3s_csf.nii'])
copyfile(segbatch{1}.spm.spatial.preproc.channel.vols{1},[out_dir filesep 's',subjid,'_structural.nii']) % here we just put the structural into the "data"-folder!

%% tell the user that now the functional experiment can start
Status_text.String = 'Starting structural procedure.. done';
Status_text.String = char(Status_text.String, 'Everything should be ready for functional run.');