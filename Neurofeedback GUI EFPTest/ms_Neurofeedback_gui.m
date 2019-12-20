function varargout = ms_Neurofeedback_gui(varargin)
% MS_NEUROFEEDBACK_GUI MATLAB code for ms_Neurofeedback_gui.fig
%      MS_NEUROFEEDBACK_GUI, by itself, creates a new MS_NEUROFEEDBACK_GUI or raises the existing
%      singleton*.
%
%      H = MS_NEUROFEEDBACK_GUI returns the handle to a new MS_NEUROFEEDBACK_GUI or the handle to
%      the existing singleton*.
%
%      MS_NEUROFEEDBACK_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MS_NEUROFEEDBACK_GUI.M with the given input arguments.
%
%      MS_NEUROFEEDBACK_GUI('Property','Value',...) creates a new MS_NEUROFEEDBACK_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ms_Neurofeedback_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ms_Neurofeedback_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ms_Neurofeedback_gui

% Last Modified by GUIDE v2.5 01-Feb-2018 15:05:44

addpath('C:\Program Files\spm12') % added by Chris Paret on May 2nd 2018S

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ms_Neurofeedback_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @ms_Neurofeedback_gui_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before ms_Neurofeedback_gui is made visible.
function ms_Neurofeedback_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ms_Neurofeedback_gui (see VARARGIN)

% Choose default command line output for ms_Neurofeedback_gui
handles.output = hObject;

% set a timer to check for MPRAGE slices; start the timer after the preRun
handles.timer_checkMPRAGEslices = timer;
handles.timer_checkMPRAGEslices.Name = 'checkMPRAGEslices';
handles.timer_checkMPRAGEslices.Period = 2; % check every 2 seconds
handles.timer_checkMPRAGEslices.ExecutionMode = 'fixedSpacing';
handles.timer_checkMPRAGEslices.TimerFcn = @check_MPRAGEslices_byTimer;
% handles.timer_checkMPRAGEslices.TasksToExecute = 1; start(handles.timer_checkMPRAGEslices) % for testing

% make sure that the GUI is set in path
addpath(fileparts(which('ms_Neurofeedback_gui')))

% maybe spm is also a good idea
if isempty(which('spm')) 
    h=warndlg('SPM was not found! Please set path.');
    waitfor(h);
    d=uigetdir; addpath(d);
end

% to add the possibility of setting an Experiment Setup File first add an
% empty field
handles.ExpSetUpFile=[];

% for checking if everything is ready
handles.SubjectFolderSet=0;
handles.ExpSetupSet=0;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ms_Neurofeedback_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);

function check_MPRAGEslices_byTimer(mTimer,~)
% this function searches for the MPRAGE slices periodically; by a timer
% if all slices are present the user is asked to go on:
% if the answer is Yes: stop timer and run structural
% if the answer is Cancel: only stop timer 
% if you say No the timer will go on! in case the set expected slices is
% too low
%% for realtime: first we need to get the real actual set values 
% to make this work you have to set the figure visibility to on!
handles.DicomOutDir_text=findobj('Tag','DicomOutDir_text');
handles.SubjectFolder_text=findobj('Tag','SubjectFolder_text');
handles.MPRAGE_ScanNum_edit=findobj('Tag','MPRAGE_ScanNum_edit');
handles.ROIs_listbox = findobj('Tag','ROIs_listbox');
handles.expectNumSlices_edit = findobj('Tag','expectNumSlices_edit');
handles.Status_text = findobj('Tag','Status_text');

% if isempty(handles.Status_text); disp('MPRAGE Timer: no figure found'); return; end
%% set all parameters 
series = str2double(handles.MPRAGE_ScanNum_edit.String); % structural 3D series number
out_dir = handles.DicomOutDir_text.String; % Dicom output dir
dicom_dir = handles.SubjectFolder_text.String; % subject/rawdcm dir
allrois = cellstr(handles.ROIs_listbox.String)'; % ROIs
expectedSlices = str2double(handles.expectNumSlices_edit.String);
handles.Status_text.String = 'Checking for MPRAGE slices..';
%% check IMAs
DCMs = ms_ScanDCMfolder(dicom_dir);

if ~isempty(DCMs) % at this point that shouldn't happen anyway.. just to be sure
    if any([DCMs.seriesNumber] == series)
        slicesMPR = length(DCMs([DCMs.seriesNumber]==series).DCMs);
        handles.Status_text.String = sprintf('Checking for MPRAGE slices of series %d.. found %d/%d\n', series, slicesMPR, expectedSlices) ;
        if slicesMPR == expectedSlices
            answer = questdlg(sprintf('All MPRAGE slices are present. %d / %d\n Proceed with structural run?', slicesMPR, str2double(handles.expectNumSlices_edit.String)));
            if strcmp(answer,'Yes'); stop(mTimer); ms_run_structural(series, out_dir, allrois, dicom_dir); end % if the answer is Yes: stop timer and run structural
            if strcmp(answer,'Cancel'); stop(mTimer); end % if the answer is Cancel: only stop timer 
            % that means: if you say NO the timer will go on!
        end
    else
        handles.Status_text.String = sprintf('Checking for MPRAGE slices of series %d.. no such series yet\n', series) ;
    end
end

% --- Outputs from this function are returned to the command line.
function varargout = ms_Neurofeedback_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in RunFunctional_pushbutton.
function RunFunctional_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to RunFunctional_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%% set all parameters and run "run_functional_kalman_detrend"
series_number = str2double(handles.fMRI_ScanNum_edit.String);
out_dir = handles.DicomOutDir_text.String; % output dir
plotParadigm(handles) % redraw the paradigm preview
ms_run_functional_kalman_detrend(series_number, out_dir, handles.ParadigmPreview_axes)

% --- Executes on button press in RunStructural_pushbutton.
function RunStructural_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to RunStructural_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% set all parameters 
series = str2double(handles.MPRAGE_ScanNum_edit.String); % structural 3D series number
out_dir = handles.DicomOutDir_text.String; % Dicom output dir
dicom_dir = handles.SubjectFolder_text.String;
allrois = handles.allrois; % ROIs
%% check if all slices are there
DCMs = ms_ScanDCMfolder(dicom_dir);
if ~isempty(DCMs) % at this point that shouldn't happen anyway.. just to be sure
    slicesMPR = length(DCMs([DCMs.seriesNumber]==series).DCMs);
    if slicesMPR ~= str2double(handles.expectNumSlices_edit.String)
        answer = questdlg(sprintf('Number of MPRAGE slices is not correct. %d / %d\n Proceed anyway?', slicesMPR, str2double(handles.expectNumSlices_edit.String)));
        if ~strcmp(answer,'Yes'); return; end % if the answer is anything than Yes.. just end here
    end
end
%% run "run_structural"
ms_run_structural(series, out_dir, allrois, dicom_dir)


% --- Executes on button press in PreRun_pushbutton.
function PreRun_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to PreRun_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%% first check if the folder exists!
if exist(handles.SubjectFolder_text.String,'dir')
    [handles.DicomOutDir_text.String, handles.allrois, handles.configMat] = ms_prerun_config(handles.SubjectFolder_text.String, handles.ExpSetUpFile);
    handles.ROIs_listbox.String = handles.allrois;
    load([handles.DicomOutDir_text.String filesep 'config.mat']); 
    if exist('expMPRAGEseries','var'); handles.MPRAGE_ScanNum_edit.String = num2str(expMPRAGEseries); end
    if exist('expMPRAGEslices','var'); handles.expectNumSlices_edit.String = num2str(expMPRAGEslices); end
    if exist('expFUNCTseries','var'); handles.fMRI_ScanNum_edit.String = num2str(expFUNCTseries); end
    stop(timerfind('Name','checkMPRAGEslices')); start(handles.timer_checkMPRAGEslices); % start the timer for checking the MPRAGE slices
    plotParadigm(handles);
    guidata(hObject, handles); % needed for allrois.. maybe we can put it somewhere directly into the GUI
else
    warndlg('Cannot find the subject folder! Please check your settings')
end


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over SubjectFolder_text.
function SubjectFolder_text_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to SubjectFolder_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% disp('hi')
% d = uigetdir('C:\AmyUp\rawdata'); % set this in extProp.m file!!
d = uigetdir('C:\BrainVoyager');
if d ~= 0 % 0 means it was canceled
    handles.SubjectFolder_text.String = d;
    handles.SubjectFolderSet=1;
    guidata(hObject,handles)
    checkExpPreSettings(handles); % simple function to check if everything is ready to go 
end

function checkExpPreSettings(handles)
if handles.SubjectFolderSet && handles.ExpSetupSet
    % turn on the buttons
    handles.PreRun_pushbutton.Enable = 'on';
    handles.RunStructural_pushbutton.Enable = 'on';
    handles.RunFunctional_pushbutton.Enable = 'on';
end

function MPRAGE_ScanNum_edit_Callback(hObject, eventdata, handles)
% hObject    handle to MPRAGE_ScanNum_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MPRAGE_ScanNum_edit as text
%        str2double(get(hObject,'String')) returns contents of MPRAGE_ScanNum_edit as a double


% --- Executes during object creation, after setting all properties.
function MPRAGE_ScanNum_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MPRAGE_ScanNum_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function expectNumSlices_edit_Callback(hObject, eventdata, handles)
% hObject    handle to expectNumSlices_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of expectNumSlices_edit as text
%        str2double(get(hObject,'String')) returns contents of expectNumSlices_edit as a double


% --- Executes during object creation, after setting all properties.
function expectNumSlices_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to expectNumSlices_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function fMRI_ScanNum_edit_Callback(hObject, eventdata, handles)
% hObject    handle to fMRI_ScanNum_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fMRI_ScanNum_edit as text
%        str2double(get(hObject,'String')) returns contents of fMRI_ScanNum_edit as a double


% --- Executes during object creation, after setting all properties.
function fMRI_ScanNum_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fMRI_ScanNum_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in ROIs_listbox.
function ROIs_listbox_Callback(hObject, eventdata, handles)
% hObject    handle to ROIs_listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ROIs_listbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ROIs_listbox


% --- Executes during object creation, after setting all properties.
function ROIs_listbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ROIs_listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
stop(timerfind('Name','checkMPRAGEslices'))
delete(timerfind('Name','checkMPRAGEslices'))
% Hint: delete(hObject) closes the figure
delete(hObject);

function plotParadigm(handles)
handles.ParadigmPreview_axes.Visible = 'on';
load([handles.DicomOutDir_text.String 'config.mat'])
% plotting
u=unique(paradigm.Scheme); c=[1 1 1];
for ix=1:length(u)
    area((paradigm.Scheme==u(ix))*4-2,'BaseValue',-2,'FaceColor',c-0.2*(ix-1),'LineStyle','none','Parent',handles.ParadigmPreview_axes); % zeros
    if ix==1; hold on; end;
end
hold off
% area((paradigm.Scheme+1==1)*4-2,'BaseValue',-2,'FaceColor',[1 1 1],'LineStyle','none','Parent',handles.ParadigmPreview_axes); % zeros
% hold on;
% area((paradigm.Scheme==1)*4-2,'BaseValue',-2,'FaceColor',[0.8 0.8 0.8],'LineStyle','none','Parent',handles.ParadigmPreview_axes); % ones
% if any(paradigm.Scheme==2)
%     area((paradigm.Scheme==2)*4-2,'BaseValue',-2,'FaceColor',[.5 .5 .5],'LineStyle','none','Parent',handles.ParadigmPreview_axes); % twos
% end
% hold off
if isfield(paradigm,'Label')
    legend(paradigm.Label,'Location','northwestoutside')
else
    legend('view','regulate','Location','northwestoutside')
end
handles.ParadigmPreview_axes.XLim = [0 length(paradigm.Scheme)];
handles.ParadigmPreview_axes.YTick =[];
[~,locs]=findpeaks(double(paradigm.Scheme~=0));
handles.ParadigmPreview_axes.XTick =[locs; length(paradigm.Scheme)];
handles.ParadigmPreview_axes.XColor= [0.5 0 0];
% handles.ParadigmPreview_axes.XAxis.TickLength = [0.0100 0.0250]*2; does
% not work in 2014b


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over ExpSetUp_text.
function ExpSetUp_text_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to ExpSetUp_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%% here we choose the "SetExpProps" file 
filt =[fileparts(which('ms_Neurofeedback_gui.m')) filesep '*etExpProps*.m'];
[f,d]=uigetfile(filt, 'Select the Experiment Setup File');
if d~=0
    addpath(d); % add path of this "function"
    handles.ExpSetUp_text.String = f;
    handles.ExpSetUpFile = f(1:end-2); % make a function out of it; remove .m
    handles.ExpSetupSet=1;
    guidata(hObject,handles)
    checkExpPreSettings(handles)
end


% --- Executes on selection change in ROIs_listbox.
function listbox2_Callback(hObject, eventdata, handles)
% hObject    handle to ROIs_listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ROIs_listbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ROIs_listbox


% --- Executes during object creation, after setting all properties.
function listbox2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ROIs_listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
