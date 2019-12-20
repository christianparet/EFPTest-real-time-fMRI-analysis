function varargout = ms_gui_defineParadigm(varargin)
% MS_GUI_DEFINEPARADIGM MATLAB code for ms_gui_defineParadigm.fig
%      MS_GUI_DEFINEPARADIGM, by itself, creates a new MS_GUI_DEFINEPARADIGM or raises the existing
%      singleton*.
%
%      H = MS_GUI_DEFINEPARADIGM returns the handle to a new MS_GUI_DEFINEPARADIGM or the handle to
%      the existing singleton*.
%
%      MS_GUI_DEFINEPARADIGM('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MS_GUI_DEFINEPARADIGM.M with the given input arguments.
%
%      MS_GUI_DEFINEPARADIGM('Property','Value',...) creates a new MS_GUI_DEFINEPARADIGM or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ms_gui_defineParadigm_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ms_gui_defineParadigm_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ms_gui_defineParadigm

% Last Modified by GUIDE v2.5 02-Oct-2017 10:29:14

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ms_gui_defineParadigm_OpeningFcn, ...
                   'gui_OutputFcn',  @ms_gui_defineParadigm_OutputFcn, ...
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


% --- Executes just before ms_gui_defineParadigm is made visible.
function ms_gui_defineParadigm_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ms_gui_defineParadigm (see VARARGIN)

% Choose default command line output for ms_gui_defineParadigm
handles.output = hObject;

% table stuff
handles.paradigm_uitable.ColumnName = {'Type','Reps'};
% define a block which consists of view, rest, regulate,...
handles.paradigm_uitable.Data{1,1}='Wait Begin';
handles.paradigm_uitable.Data{2,1}='View';
handles.paradigm_uitable.Data{3,1}='Rest';
handles.paradigm_uitable.Data{4,1}='Regulate';
handles.paradigm_uitable.Data{5,1}='#Blocks';
handles.paradigm_uitable.Data{6,1}='Wait End';
handles.paradigm_uitable.Data(:,2)={5, 18, 21, 18, 24, 2};
handles.paradigm_uitable.ColumnEditable = [false, true];
% plot something
doPlot(handles)
% how to deal with the data?
% simply saving that stuff somewhere? there is a config.mat file.. maybe
% there?
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ms_gui_defineParadigm wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ms_gui_defineParadigm_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes when entered data in editable cell(s) in paradigm_uitable.
function paradigm_uitable_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to paradigm_uitable (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)

doPlot(handles)

function doPlot(handles)
data = handles.paradigm_uitable.Data;
waitBegin = data{contains(data(:,1), 'Wait Begin'),2};
waitEnd = data{contains(data(:,1), 'Wait End'),2};
viewRep = data{contains(data(:,1), 'View'),2};
restRep = data{contains(data(:,1), 'Rest'),2};
regRep = data{contains(data(:,1), 'Regulate'),2};
numBlocks = data{contains(data(:,1), '#Blocks'),2};

% define a block; not really clear how a block is defined.. especially the
% time between the blocks
% block: view - rest - regulate - rest
block = [ones(viewRep,1); zeros(restRep,1); ones(regRep,1)*2; zeros(restRep,1);];

% put all together
paradigm = [zeros(waitBegin,1); repmat(block,[numBlocks 1]); zeros(waitEnd,1)];
area((paradigm==1)-2,'BaseValue',-2,'FaceColor',[0 0 0],'LineStyle','none','Parent',handles.axes1); % view
hold on
area((paradigm==2)-2,'BaseValue',-2,'FaceColor',[.5 .5 .5],'LineStyle','none','Parent',handles.axes1); % regulate
hold off
legend('view','regulate','Location','northwestoutside')
handles.axes1.XLim = [0 length(paradigm)];
handles.axes1.YTick =[];
[~,locs]=findpeaks(double(paradigm~=0));
handles.axes1.XTick =[locs; length(paradigm)];
handles.axes1.XColor= [0.5 0 0];
handles.axes1.XAxis.TickLength = [0.0100 0.0250]*2;
