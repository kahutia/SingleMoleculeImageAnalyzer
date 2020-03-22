%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%             Trotsky V8               by SHK        %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

function varargout = trotsky(varargin)
% TROTSKY MATLAB code for trotsky.fig
%      TROTSKY, by itself, creates a new TROTSKY or raises the existing
%      singleton*.
%
%      H = TROTSKY returns the handle to a new TROTSKY or the handle to the
%      existing singleton*.
%
%      TROTSKY('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TROTSKY.M with the given input arguments.
%
%      TROTSKY('Property','Value',...) creates a new TROTSKY or raises the
%      existing singleton*.  Starting from the left, property value pairs
%      are applied to the GUI before trotsky_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property
%      application stop.  All inputs are passed to trotsky_OpeningFcn via
%      varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help trotsky

% Last Modified by GUIDE v2.5 10-Feb-2020 17:00:58

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @trotsky_OpeningFcn, ...
    'gui_OutputFcn',  @trotsky_OutputFcn, ...
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


% --- Executes just before trotsky is made visible.
function trotsky_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn. hObject    handle to
% figure eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA) varargin
% command line arguments to trotsky (see VARARGIN)

% Choose default command line output for trotsky
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes trotsky wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = trotsky_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT); hObject
% handle to figure eventdata  reserved - to be defined in a future version
% of MATLAB handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



% --- Executes on button press in CreatMap.
function CreatMap_Callback(hObject, eventdata, handles)
% hObject    handle to CreatMap (see GCBO) eventdata  reserved - to be
% defined in a future version of MATLAB handles    structure with handles
% and user data (see GUIDATA)
% localarea_size  = str2double(get(handles.local_area_size,'String'));
% threshold_index = str2double(get(handles.threshold_index,'String'));
mapcreator(hObject, handles);


% --- Executes on button press in clearmap.
function clearmap_Callback(hObject, eventdata, handles)
% hObject    handle to clearmap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isfield(handles,'mapcoef')
    handles=rmfield(handles,'mapcoef');
    add_log(handles,'Map coefficient were cleared');
    guidata(hObject,handles);
end


% --- Executes on button press in LoadMap.
function LoadMap_Callback(hObject, eventdata, handles)
% hObject    handle to LoadMap (see GCBO) eventdata  reserved - to be
% defined in a future version of MATLAB handles    structure with handles
% and user data (see GUIDATA) disp('load mapping coefficients');

if ~isfield(handles,'pma_root_path')
    handles.pma_root_path='';
end

[mapfilename, mappathname]=uigetfile(fullfile(handles.pma_root_path, '*.mapcoef'), 'mapping file');
handles.pma_root_path=mappathname;

mapfid=fopen(fullfile(mappathname, mapfilename),'r');
% version 1
% map=zeros(2,7);
% for j=1:2
%     for i=1:7
%         fscanf(mapfid,'%s',1);
%         map(j,i)=fscanf(mapfid,'%g',1);
%     end
% end
% version 2
map=zeros(2,9);
for j=1:2
    for i=1:10
        fscanf(mapfid,'%s',1);
        map(j,i)=fscanf(mapfid,'%g',1);
    end
end

handles.mapcoef=map;
guidata(hObject,handles);
fclose(mapfid);


%% update log
% logtext = get(handles.Logbox,'String');
% [num_entry, ~]=size(logtext);
% if num_entry>9,    logtext=logtext(2:end);    num_entry=num_entry-1;    end
% logtext{num_entry+1}=['mapping file loaded. (' mappathname mapfilename ')'];
% set(handles.Logbox,'String',logtext);
add_log(handles,['mapping file loaded. (' mappathname mapfilename ')']);    

% --- Executes on button press in Loadpmabutton.
function Loadpmabutton_Callback(hObject, eventdata, handles)
% hObject    handle to Loadpmabutton (see GCBO) eventdata  reserved - to be
% defined in a future version of MATLAB handles    structure with handles
% and user data (see GUIDATA)

if isfield('pmaloaded',handles)
    rmfield('pmaloaded',handles);
end

if isfield(handles,'pma_root_path')
%     disp('pg data available');
%     disp(handles.pma_root_path);
else
    handles.pma_root_path='';
end
[pmafilename, pmapathname]=uigetfile(fullfile(handles.pma_root_path, '*.*'), 'image files','MultiSelect','on');
handles.pma_root_path=pmapathname;
if iscell(pmafilename)
    handles.pmaloaded.filename=pmafilename;
    handles.pmaloaded.pathname=pmapathname;
    [~, num_pmafiles]=size(pmafilename);
else
    pmafilename=cellstr(pmafilename);
    handles.pmaloaded.filename=pmafilename;
    handles.pmaloaded.pathname=pmapathname;
    num_pmafiles=1;
end
handles.num_pmafiles_all=num_pmafiles;
handles.file_mode='individual files';
guidata(hObject,handles);

%% update log
add_log(handles,[num2str(num_pmafiles) ' image files loaded.']);    


% --- Executes on button press in Loadpmafolder.
function Loadpmafolder_Callback(hObject, eventdata, handles)
% hObject    handle to Loadpmafolder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isfield(handles,'pmaloaded')
    handles=rmfield(handles,'pmaloaded');
end

if isfield(handles,'pma_root_path')
%     disp('pg data available');
    disp(handles.pma_root_path);
else
    handles.pma_root_path='';
end
[pmapathname]=uigetdir(handles.pma_root_path, 'Select folder containing image files');

add_log(handles,'Searching for image files...');    guidata(hObject,handles);

subdirs_tmp=get_pma_all_subdir(pmapathname);

% remove bead and LeonResult folder
vid=true(length(subdirs_tmp),1);
for tfi=1:length(subdirs_tmp)
    if contains(subdirs_tmp{tfi},'bead') || contains(subdirs_tmp{tfi},'bead')
        vid(tfi)=false;
    end
end
subdirs=subdirs_tmp(vid);
    

handles.pma_root_path=pmapathname;
handles.pma_sub_dir_info=subdirs;
num_pmafiles=[];
for dri=1:length(subdirs)
    %         f_info=[dir([pmapathname '\' subdirs(dri).name '\*.tif']); dir([pmapathname '\' subdirs(dri).name '\*.pma'])];
    f_info_tmp=[dir(fullfile(subdirs{dri}, '*.tif')); dir(fullfile(subdirs{dri}, '*.pma')); dir(fullfile(subdirs{dri}, '*.shf'))] ;
    
    % remove images generated by Tiff_array_viewer (keyword='_TAV')
    vid1=true(length(f_info_tmp),1);
    for tfi=1:length(f_info_tmp)
        if contains(f_info_tmp(tfi).name,'_TAV')
            vid1(tfi)=false;
        end
    end
%     f_info=f_info_tmp(vid);
    
    % remove Tiff appendix file
    vid2=true(length(f_info_tmp),1);
    for tfi=1:length(f_info_tmp)
        if contains(f_info_tmp(tfi).name,'_X')
            vid2(tfi)=false;
        end
    end
    f_info=f_info_tmp(vid1 & vid2);
    
    handles.pmafileinfo{dri}=f_info;
    handles.pmapathinfo{dri}=subdirs{dri};
    num_pmafiles(dri)=length(f_info);
end
handles.num_pmafiles_all=sum(num_pmafiles(:));
handles.file_mode='file folders';
add_log(handles,'Searching for image files...');
add_log(handles,[num2str(sum(num_pmafiles(:))) ' image files in ' num2str(length(subdirs)) ' sub-folders found.']);
guidata(hObject,handles);



% --- Executes on button press in Backup.
function Backup_Callback(hObject, eventdata, handles)
% hObject    handle to Backup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


if ~isfield(handles,'pma_root_path')
    handles.pma_root_path='';
end

backup_leon_result(handles.pma_root_path);



% --- Executes on button press in show_stat.
function show_stat_Callback(hObject, eventdata, handles)
% hObject    handle to show_stat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isfield(handles,'pma_root_path')
%     disp('pg data available');
    disp(handles.pma_root_path);
else
    handles.pma_root_path='';
end
% ReadTrotsky(handles.pma_root_path)
rootpath=handles.pma_root_path;
run('ReadTrotsky.m')

% --- Executes on button press in processallbutton.
function processallbutton_Callback(hObject, eventdata, handles)
% hObject    handle to processallbutton (see GCBO) eventdata  reserved - to
% be defined in a future version of MATLAB handles    structure with
% handles and user data (see GUIDATA)
warning('off')
process_pmafiles(handles,hObject);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function local_area_size_Callback(hObject, eventdata, handles)
% hObject    handle to local_area_size (see GCBO) eventdata  reserved - to
% be defined in a future version of MATLAB handles    structure with
% handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function local_area_size_CreateFcn(hObject, eventdata, handles)
% hObject    handle to local_area_size (see GCBO) eventdata  reserved - to
% be defined in a future version of MATLAB handles    empty - handles not
% created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function Maskwidth_Callback(hObject, eventdata, handles)
% hObject    handle to Maskwidth (see GCBO) eventdata  reserved - to be defined
% in a future version of MATLAB handles    structure with handles and user
% data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function Maskwidth_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Maskwidth (see GCBO) eventdata  reserved - to be defined
% in a future version of MATLAB handles    empty - handles not created
% until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function first_frame_peakfinder_Callback(hObject, eventdata, handles)
% hObject    handle to first_frame_peakfinder (see GCBO) eventdata  reserved -
% to be defined in a future version of MATLAB handles    structure with
% handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function first_frame_peakfinder_CreateFcn(hObject, eventdata, handles)
% hObject    handle to first_frame_peakfinder (see GCBO) eventdata  reserved -
% to be defined in a future version of MATLAB handles    empty - handles
% not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function time_bin_Callback(hObject, eventdata, handles)
% hObject    handle to time_bin (see GCBO) eventdata  reserved - to be
% defined in a future version of MATLAB handles    structure with handles
% and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of time_bin as text
%        str2double(get(hObject,'String')) returns contents of time_bin as
%        a double


% --- Executes during object creation, after setting all properties.
function time_bin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to time_bin (see GCBO) eventdata  reserved - to be
% defined in a future version of MATLAB handles    empty - handles not
% created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function threshold_index_Callback(hObject, eventdata, handles)
% hObject    handle to threshold_index (see GCBO) eventdata  reserved - to
% be defined in a future version of MATLAB handles    structure with
% handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of threshold_index as text
%        str2double(get(hObject,'String')) returns contents of
%        threshold_index as a double


% --- Executes during object creation, after setting all properties.
function threshold_index_CreateFcn(hObject, eventdata, handles)
% hObject    handle to threshold_index (see GCBO) eventdata  reserved - to
% be defined in a future version of MATLAB handles    empty - handles not
% created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function fixedBG_d_Callback(hObject, eventdata, handles)
% hObject    handle to fixedBG_d (see GCBO) eventdata  reserved - to be defined
% in a future version of MATLAB handles    structure with handles and user
% data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of fixedBG_d as text
%        str2double(get(hObject,'String')) returns contents of fixedBG_d as a
%        double


% --- Executes during object creation, after setting all properties.
function fixedBG_d_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fixedBG_d (see GCBO) eventdata  reserved - to be defined
% in a future version of MATLAB handles    empty - handles not created
% until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in DriftCorrection.
function DriftCorrection_Callback(hObject, eventdata, handles)
% hObject    handle to DriftCorrection (see GCBO) eventdata  reserved - to
% be defined in a future version of MATLAB handles    structure with
% handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of DriftCorrection



% --- Executes on button press in subfolder.
function subfolder_Callback(hObject, eventdata, handles)
% hObject    handle to subfolder (see GCBO) eventdata  reserved - to be
% defined in a future version of MATLAB handles    structure with handles
% and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of subfolder


% --- Executes on button press in MergedImg.
function MergedImg_Callback(hObject, eventdata, handles)
% hObject    handle to MergedImg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in AcceptorImg.
function AcceptorImg_Callback(hObject, eventdata, handles)
% hObject    handle to AcceptorImg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in DonorImg.
function DonorImg_Callback(hObject, eventdata, handles)
% hObject    handle to DonorImg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes on button press in GaussFilter.
function GaussFilter_Callback(hObject, eventdata, handles)
% hObject    handle to GaussFilter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


function Npixel2sum_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function Npixel2sum_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Npixel2sum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Num_frame_peakfind_Callback(hObject, eventdata, handles)
% hObject    handle to Num_frame_peakfind (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --- Executes during object creation, after setting all properties.
function Num_frame_peakfind_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Num_frame_peakfind (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function co_local_tolerance_Callback(hObject, eventdata, handles)
% hObject    handle to co_local_tolerance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of co_local_tolerance as text
%        str2double(get(hObject,'String')) returns contents of co_local_tolerance as a double


% --- Executes during object creation, after setting all properties.
function co_local_tolerance_CreateFcn(hObject, eventdata, handles)
% hObject    handle to co_local_tolerance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in show_images.
function show_images_Callback(hObject, eventdata, handles)
% hObject    handle to show_images (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of show_images



function ColorEnhancerGreen_Callback(hObject, eventdata, handles)
% hObject    handle to ColorEnhancerGreen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ColorEnhancerGreen as text
%        str2double(get(hObject,'String')) returns contents of ColorEnhancerGreen as a double


% --- Executes during object creation, after setting all properties.
function ColorEnhancerGreen_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ColorEnhancerGreen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ColorEnhancerRed_Callback(hObject, eventdata, handles)
% hObject    handle to ColorEnhancerRed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ColorEnhancerRed as text
%        str2double(get(hObject,'String')) returns contents of ColorEnhancerRed as a double


% --- Executes during object creation, after setting all properties.
function ColorEnhancerRed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ColorEnhancerRed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function Logbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Logbox (see GCBO) eventdata  reserved - to be
% defined in a future version of MATLAB handles    empty - handles not
% created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes during object creation, after setting all properties.
function fixedBG_a_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fixedBG_a (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes during object creation, after setting all properties.
function img_bg_green_CreateFcn(hObject, eventdata, handles)
% hObject    handle to img_bg_green (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes during object creation, after setting all properties.
function img_bg_red_CreateFcn(hObject, eventdata, handles)
% hObject    handle to img_bg_red (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function img_bg_green_Callback(hObject, eventdata, handles)
% hObject    handle to img_bg_green (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of img_bg_green as text
%        str2double(get(hObject,'String')) returns contents of img_bg_green as a double


% --- Executes on button press in img_autoscale.
function img_autoscale_Callback(hObject, eventdata, handles)
% hObject    handle to img_autoscale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of img_autoscale


% --- Executes on button press in peak_all_fr.
function peak_all_fr_Callback(hObject, eventdata, handles)
% hObject    handle to peak_all_fr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of peak_all_fr


function img_bg_red_Callback(hObject, eventdata, handles)
% hObject    handle to img_bg_red (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of img_bg_red as text
%        str2double(get(hObject,'String')) returns contents of img_bg_red as a double


% --- Executes on button press in ext_trace.
function ext_trace_Callback(hObject, eventdata, handles)
% hObject    handle to ext_trace (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ext_trace

function Logbox_Callback(hObject, eventdata, handles)
% hObject    handle to Logbox (see GCBO) eventdata  reserved - to be
% defined in a future version of MATLAB handles    structure with handles
% and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of Logbox as text
%        str2double(get(hObject,'String')) returns contents of Logbox as a
%        double


% --- Executes on button press in img_pivot.
function img_pivot_Callback(hObject, eventdata, handles)
% hObject    handle to img_pivot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of img_pivot

% 
% % --- Executes on button press in radiobutton34.
% function radiobutton34_Callback(hObject, eventdata, handles)
% % hObject    handle to radiobutton34 (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% 
% % Hint: get(hObject,'Value') returns toggle state of radiobutton34
% 
% 
% 
% function fixedBG_a_Callback(hObject, eventdata, handles)
% % hObject    handle to fixedBG_a (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% 
% % Hints: get(hObject,'String') returns contents of fixedBG_a as text
% %        str2double(get(hObject,'String')) returns contents of fixedBG_a as a double


% --- Executes on button press in Use_saved_img.
function Use_saved_img_Callback(hObject, eventdata, handles)
% hObject    handle to Use_saved_img (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Use_saved_img
