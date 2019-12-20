function [outputdir, allrois] = ms_SetExpProps_EFPpre(dicom_dir)
% this function is used to set and save the general setup of a
% study/experiment

%% since the file-structure follows a simple scheme we only need the "dicom-dir" (rawData dir)
% e.g. dicom_dir = C:\AmyUp\rawdata\20160621.211011.211011
% then: homedir = C:\AmyUp\; outputdir = D:\AmyUp\data\211011 (last folder is subjectID!)
[p,n]=fileparts(dicom_dir);
homedir = [fileparts(p) filesep]; % define homedir; last filesep needed?
loc=strfind(n,'.'); subjid = n(loc(end)+1:end); % retrieve subject id from data folder
outputdir=[homedir 'data' filesep subjid filesep]; % output directory
if ~exist(outputdir,'dir'); mkdir(outputdir); end
clear p n loc

%% expected MPRAGE settings
expMPRAGEseries = 2; expMPRAGEslices = 192;
%% expected Functional series
expFUNCTseries=6;

%% set paradigm parameters
% define duration of experimental blocks/trials (unit: TR)
paradigm.WaitBegin = 25; % number of TRs to wait until start of experiment (i.e. first condition)
paradigm.FT = 15; % condition: Finger Tapping (FT)
paradigm.Rest = 30; % condition: Rest
paradigm.Regulate = 35; % condition: Regulate
paradigm.NumBlocks = 2; % A block is a sequence of trials defined above
paradigm.WaitEnd = 3; % additional TRs after last TR of last trial
paradigm.Label = {'Rest', 'Regulate', 'FT'}; % that means that Rest has to get a 0, 'Regulate' a 1, 'FT' gets a 2 and so on!

paradigm.Block = [zeros(paradigm.Rest,1); ones(paradigm.Regulate,1); ones(paradigm.FT,1)*2];
paradigm.Scheme = [zeros(paradigm.WaitBegin,1); repmat(paradigm.Block,[paradigm.NumBlocks 1]); zeros(paradigm.WaitEnd,1)];

maxvolumes = length(paradigm.Scheme); % needed for later run_functional

%% figure set-up
screen = [1920 1080]; % resolution, for graphics; adapt to screen dimensions
min_factor = 1.25; % factor chosen in the Anzeigeeinstellungen to improve readibility on display
screen = screen/1.25;

panel_brain_ratio = 566.67/816.94;
height_panel_brain = screen(2)*2/3;
width_panel_brain = height_panel_brain*panel_brain_ratio;
pos_panel_brain = [screen(1)-width_panel_brain*2 0 width_panel_brain height_panel_brain];
pos_panel_timecourse = [screen(1)/3 screen(2)*2/3 (screen(1))*2/3 (screen(2))/3];

%% the following are the settings for the functional run
mask_dir = 'C:\BrainVoyager\Masks'; 
allrois = {[mask_dir '\Right Amygdala 25.nii'],...
            [mask_dir '\rect.nii']};

% used for the server/functional run
timeout = 600;
nrROIs = 2;
TR = 2;

wait_detrend = 18; % wait at least 18 scans before detrending starts, needed for filter stabilization

roispec{1}.srcimg = 'wRight Amygdala 25.nii'; % used in ms_run_functional_kalman_detrend

roispec{2}.srcimg = 'wrect.nii'; % used in ms_run_functional_kalman_detrend

% realign options (not checked) % used in ms_run_functional_kalman_detrend
rflags.quality = 0.5; % decrease if needed, default is 0.9
rflags.fwhm = 5;
rflags.rtm = 0;
rflags.interp =  1; %trilinear
rflags.sep = 8; % increase if possible, default is 4
rflags.graphics = 0;

% reslice options (not checked) % used in ms_run_functional_kalman_detrend
rsflags = struct('mask', false,...
                 'mean', false,...
                 'interp', 2,... % spline order (0 = nearest neighbour, 1=trilinear, 2 = 2nd order bspline 
                 'which', 1,...
                 'wrap', [0; 0; 0],...
                 'prefix', 'r');
            
% coregister options
corflags = struct('graphics', 0,...
                  'sep', 3); % 5
            
% Filter settings % used in ms_run_functional_kalman_detrend
% Kalman covariance matrices and threshold definition from Yuri Koushs Script "test_feedbackSP_subsequentROIs.m"
S.Q = 0; 
S.P = S.Q; 
S.x = 0;
S(1:nrROIs) = S;
fposDer = zeros(nrROIs,1);
fnegDer = zeros(nrROIs,1);

%initialize movement regressors
rp = [0 0 0 0 0 0]; 

% save experiment settings
rtconfig = struct('dicom_dir', dicom_dir,...
                    'timeout', timeout,...
                    'maxvolumes',maxvolumes,...
                    'port', 8082); %8082
                
%% jump into the outputdir and save everything as config.mat
cd(outputdir)
save config
fprintf('PreRun: Saved settings in config.mat in %s\n', outputdir)