function [outputdir, allrois, configpath] = ms_prerun_config(dicom_dir, SetupFile)
%% all we need in the beginning is the "dicom_dir"; which is the rawdata dir
% idea: there is a ms_SetExpProps_*.m file for each study, for example
% ms_SetExpProps_AmyUp.m; To switch a study or test other settings simply
% create another m-file and replace it here!

% [outputdir, allrois] = ms_SetExpProps_EFPpre(dicom_dir);
% to define a setup file we use an eval
eval(['[outputdir, allrois] = ' SetupFile '(dicom_dir);'])

configpath = [outputdir filesep 'config.mat']; % I think I don't want that!

%% old stuff! 
% % homedir = 'C:\AmyUp\';
% homedir = 'D:\AmyUp\';
% addpath(homedir); % why?
% % cd([homedir,'\rawdata'])
% 
% % dicom_dir = uigetdir;
% timeout = 600;
% nrROIs = 2;
% TR = 2;
% wait_detrend = 24; % wait X scans before detrending starts, needed for filter stabilizatioon
% % window_size = 22; % sliding window for regression, needs to be saved
% 
% mask_dir = 'C:\BrainVoyager\Masks'; 
% 
% allrois = {[mask_dir '\Right Amygdala 25.nii'],...
%     [mask_dir '\rect.nii']}; 
% 
% roispec{1}.srcimg = 'wRight Amygdala 25.nii'; % used in ms_run_functional_kalman_detrend
% roispec{2}.srcimg = 'wrect.nii'; % used in ms_run_functional_kalman_detrend
% 
% % realign options (not checked) % used in ms_run_functional_kalman_detrend
% rflags.quality = 0.5; % decrease if needed, default is 0.9
% rflags.fwhm = 5;
% rflags.rtm = 0;
% rflags.interp =  1; %trilinear
% rflags.sep = 8; % increase if possible, default is 4
% rflags.graphics = 0;
% 
% % reslice options (not checked) % used in ms_run_functional_kalman_detrend
% rsflags = struct('mask', false,...
%                 'mean', false,...
%                 'interp', 2,... % spline order (0 = nearest neighbour, 1=trilinear, 2 = 2nd order bspline 
%                 'which', 1,...
%                 'wrap', [0; 0; 0],...
%                 'prefix', 'r');
%             
% % Filter settings % used in ms_run_functional_kalman_detrend
% % Kalman covariance matrices and threshold definition from Yuri Koushs Script "test_feedbackSP_subsequentROIs.m"
% S.Q = 0; 
% S.P = S.Q; 
% S.x = 0;
% S(1:nrROIs) = S;
% fposDer = zeros(nrROIs,1);
% fnegDer = zeros(nrROIs,1);
% 
% % prepare output directory
% % subjid=dicom_dir(end-5:end); % retrieve subject id from data folder
% tmp_dirname = regexp(dicom_dir,filesep,'split'); tmp_dirname = tmp_dirname{end}; % MS: that is a bit more robust!
% loc=strfind(tmp_dirname,'.'); subjid = tmp_dirname(loc(end)+1:end); % retrieve subject id from data folder
% outputdir=[homedir filesep 'data' filesep subjid];
% if ~exist(outputdir,'dir')
%     mkdir(outputdir);
% end
% cd(outputdir) % we switch to here because we want to save something in the config.mat file..
% 
% % figure set-up
% screen = [1600 900]; % resolution, for graphics
% 
% panel_brain_ratio = 566.67/816.94;
% height_panel_brain = screen(2)*2/3;
% width_panel_brain = height_panel_brain*panel_brain_ratio;
% pos_panel_brain = [screen(1)-width_panel_brain 0 width_panel_brain height_panel_brain];
% pos_panel_timecourse = [screen(1)/3 screen(2)*2/3 (screen(1))*2/3 (screen(2))/3];
% 
% maxvolumes = 467;
% 
% % save experiment settings
% rtconfig = struct('dicom_dir', dicom_dir,...
%                     'timeout', timeout,...
%                     'maxvolumes',maxvolumes,...
%                     'port', 8082); %8082
%                
% save config % problem here is: everything is stored!! what do we really need?
% configpath = [outputdir filesep 'config.mat'];