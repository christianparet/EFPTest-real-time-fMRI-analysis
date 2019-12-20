function ms_run_functional_kalman_detrend(series_number, out_dir, ax)

% change log
% 2017/12/04: Wrong syntax in realignment and reslicing procedure was corrected to realign to first image
% 2018/01/05: Added coregistration routine, developed by Tomás Slavícek.  First functional image is coregistered to structural image. Header of functional image is updated with coregistration parameters.
% 2018/03/27: Preload dicom routine to accelerate dicom import, developed by Tomás Slavícek.

% This function "works" completely in the output directory!
cd(out_dir);
%% IF RUNNING SUBJECT change to dryrun = 0; and fakeServer = 0; !!!
dryrun = 0; fakeServer = 0;
if dryrun || fakeServer 
    answer = questdlg('You run a "dry run" and/or a fake server! Is that okay?', 'Check Settings', 'It''s okay', 'Change that!', 'Change that!');
        if ~isempty(strfind(answer,'Change')); edit(mfilename); return; end
end
%% load the settings
load([out_dir filesep 'config.mat']);

%% Preload dicom
% Routine suggested by Tomás Slavícek to speed up dicom import. To use it, change code of spm_dicom_header.m as shown here:

% Original:
% 
% function dict = readdict(P)
% if nargin<1, P = 'spm_dicom_dict.mat'; end
% try
%      dict = load(P);
% catch
%     fprintf('\nUnable to load the file "%s".\n', P);
%     rethrow(lasterror);
% end
% 
% 
% Modified:
% 
% function dict = readdict(P)
% global spm_dicom_dict_var
% if nargin<1, P = 'spm_dicom_dict.mat'; end
% try
%     if ~isempty(spm_dicom_dict_var)
%         dict = spm_dicom_dict_var;
%     else
%         dict = load(P);
%     end
% catch
%     fprintf('\nUnable to load the file "%s".\n', P);
%     rethrow(lasterror);
% end

global spm_dicom_dict_var
spm_dicom_dict_var = load([spm('dir') filesep 'spm_dicom_dict.mat']);

%% start the functional stuff
spm('defaults','fmri');
spm_jobman('initcfg');

% Variables for timing monitoring
delta_img=[];
delta_prepro=[];
delta_ROIcompute=[];

try
    dicom_dir = rtconfig.dicom_dir;
    timeout = rtconfig.timeout;
    maxvolumes = rtconfig.maxvolumes;
catch
    disp('error in configuration structure');
    return
end

current_image_count = 1;

disp(dicom_dir);

V = []; % MS: good practice to initialize variables first.. but doesn't make much sense in the current form, imho
P = {};

%% Figure
update_orthview_new(['s',subjid,'_structural.nii'], roispec);
s=gcf;
set(s,'OuterPosition',pos_panel_brain)

%% draw it in the paradigm preview
hold(ax,'on')

sigtarget=plot(0,0,'r','LineWidth',2, 'Parent',ax);
sigcontrolroi=plot(0,0,'g','LineWidth',2, 'Parent',ax);
sigdiff=plot(0,0,'b-','LineWidth',2, 'Parent',ax);
% try
%     for ix=length(ax.Children):-1:4
%         try
%             le{ix}=char(ax.Children(ix).DisplayName);
%         catch
%             le{ix}=[];
%         end
%     end
%     le(1:3)=[]; ole=char(flip(le));
%     legend(ax, char(ole, 'target ROI','control ROI','delta'),'Location','northwestoutside')
%     catch
%     legend(ax,'rest','feedback', 'target ROI','control ROI','delta','Location','northwestoutside')
% end
% in Matlab2017b we don't need the stuff above! new legend entries are
% automatically added! since we added the 3 sigs, name them!

ax.Legend.String(end-2:end)={'target ROI','control ROI','cleaned signal'}; % I don't like that.. but I don't know if the sigdiff is maybe important!

hold(ax,'off')
pause(.025);
%% configure server

if (rtconfig.port > 0)
    tc = tcpip('0.0.0.0', rtconfig.port, 'NetworkRole', 'server');
    set(tc,'InputBufferSize', 4096);
    set(tc,'OutputBufferSize', 4096);
    set(tc,'Timeout', 30);
    fprintf(1, 'waiting for network connection\n');
    if ~fakeServer; fopen(tc); end
    fprintf(1, 'network open\n');
else
    fprintf('no network, running\n');
end

%% Real-time fMRI
run_start_time = clock;
volfile_present_time = run_start_time;
time_img = clock;

while (1)
    expected_fn = [dicom_dir '/001_' sprintf('%06d',series_number), '_', sprintf('%06d',current_image_count), '.dcm'];
    if dryrun; end% dry run
    lastwarn('');
    if (etime(clock, volfile_present_time) > timeout)
        disp('volume timeout, finishing');
        if rtconfig.port > 0
            fclose(tc);
        end
        save([subjid '_matout'], 'V');
        return;
    end
    
    if (exist(expected_fn, 'file'))
        
        fprintf('\nimage count: %d\n', current_image_count);
        
        delta_img=[delta_img etime(clock,time_img)];
        disp(['delta img = ', num2str(delta_img(end))]);
        time_img = clock;
        tic;
        volfile_present_time = clock;
        pause(0.025); % wait 25 ms
        try
            hdr = spm_dicom_headers(expected_fn);
            out = spm_dicom_convert(hdr, 'all', 'flat', 'nii');
        catch
            disp('dicom read error');
            continue;
        end
        msg =  lastwarn;
        if (~strcmp(msg, ''))
            disp('dicom read warning');
            continue;
        end
        if current_image_count == 1
            first_image = out.files{1};
            % coregistration
            coreg_params_disp = spm_coreg(['s',subjid,'_structural.nii'],first_image,corflags);
            M = inv(spm_matrix(coreg_params_disp));
            MM = spm_get_space(first_image);
            spm_get_space(deblank(first_image), M*MM); % wo output it "creates" first_image; applies coreg parameter
            fprintf(2,'\ncoregistered by x,y,z (mm): ');
            fprintf(2,'%4.1f ', coreg_params_disp(1:3));
            fprintf('\n');
            
            update_orthview_new(first_image, roispec);
            
            P{1} = first_image;
        else
            % realign
            spm_realign(strvcat(first_image, out.files{1}), rflags); 
            spm_reslice(strvcat(first_image, out.files{1}), rsflags);
            
            % read estimated movement parameters                    
            fid = fopen(spm_file(first_image, 'prefix','rp_', 'ext','.txt'),'r');
            movepar = textscan(fid, '%f %f %f %f %f %f', 'Delimiter','\t','HeaderLines',1);
            fclose(fid);
            rp = [rp; cell2mat(movepar)];
            
            save([out_dir,'rp'],'rp');
        end
        
        delta_prepro = [delta_prepro toc];
        disp(['delta prepro = ',num2str(delta_prepro(end))]);
        tic;
        
        try
            if current_image_count == 1
                fname = out.files{1};
            else
                [p,fn,e] = fileparts(out.files{1});
                fname = [p '/r' fn '.nii'];
                P{current_image_count} = fname;
            end
            [outr, res] = SNiP_tbxvol_extract_fast(fname , roispec, 'none');
        catch
            disp('extract error');
            continue;
        end
        
        % Kalman filter
        for roi=1:2

            V(current_image_count, roi) = mean(outr{roi});

%             if current_image_count == 0 % MS: does that make sense at all? current_image_count is initialized as "1" and only can grow 
%                 tmp_std = 0;
%             elseif current_image_count < 3 % MS: this elseif is the same as the else!!
%                 tmp_std = std(V(1:current_image_count,roi))';
%             else
%                 tmp_std = std(V(1:current_image_count,roi))';
%             end

            tmp_std = std(V(1:current_image_count,roi))'; % CP: changed; code can be removed

            if TR == 1 % MS: bad naming of TR in a fMRI thing ;) CP: think not; this means the repetition time
                % lambda = R/Q = 4 : approx. equivalent to the Butterworth Fc = 0.155 Hz, TR = 1s, Fs = 1Hz
                S(roi).Q = tmp_std.^2;
                S(roi).R = 4*S(roi).Q;
                Th_K(roi) = .9*tmp_std;
            elseif TR == 2
                % lambda = R/Q = 1.95 : approx. equivalent to the Butterworth Fc = 0.106 Hz, TR = 2s, Fs = 0.5Hz
                S(roi).Q = tmp_std.^2;
                S(roi).R = 1.95*S(roi).Q; % .001 for despiking only
                Th_K(roi) = .9*tmp_std;   %   2   for despiking only
            end

            [out_kalm_Sample(current_image_count,roi), S(roi), fposDer(roi), fnegDer(roi)] = kalman_spike(Th_K(roi), V(current_image_count,roi), S(roi), fposDer(roi), fnegDer(roi));

        end
            
%         if (current_image_count == 1) % MS: why is that here? place it in line 175 or so in the "if current_image_count == 1" condition. CP: done
%             update_orthview_new(first_image, roispec);
%         end
       
        if current_image_count > wait_detrend          
            %% Method: Subtract filtered control ROI signal from filtered amygdala signal, after detrending and percent signal change calculation 
            % This method returns values that can directly be used for
            % feedback. No baseline normalization or standardization
            % necessary.
            detrend_out_kalm_Sample = detrend(out_kalm_Sample(wait_detrend+1:end,:));

            R(1) = detrend_out_kalm_Sample(end,1)/mean(V(:,1))*100; % percent signal change
            R(2) = detrend_out_kalm_Sample(end,2)/mean(V(:,2))*100; % percent signal change
            R(3) = R(1)-R(2);
                        
        else
            R = [0.0 0.0 0.0];
        end
        
        Rfinal(1,current_image_count)=R(1);
        Rfinal(2,current_image_count)=R(2);
        Rfinal(3,current_image_count)=R(3);
        
        sigtarget.XData = 1:length(Rfinal(1,1:end));
        sigtarget.YData = Rfinal(1,1:end);
        
        sigcontrolroi.XData = 1:length(Rfinal(2,1:end));
        sigcontrolroi.YData = Rfinal(2,1:end);
        
        sigdiff.XData = 1:length(Rfinal(3,1:end));
        sigdiff.YData = Rfinal(3,1:end);
        
        if rtconfig.port > 0
            try
                fprintf(tc, '%4.2f ', R(3));
            catch
                fprintf(2, 'error writing to tcp socket\n');
            end
        end
                     
        delta_ROIcompute = [delta_ROIcompute toc];
        disp(['delta ROI compute = ',num2str(delta_ROIcompute(end))]);
        disp(['extracted ROI value = ',num2str(R(3))]);
        
        current_image_count = current_image_count + 1;
    end
    
    if current_image_count > maxvolumes % MS: I really don't like this "maxvolumes".. there is no other way?
        disp('maxvolumes exceed, finish');
        if rtconfig.port > 0
            fclose(tc);
        end
        break
    end
end