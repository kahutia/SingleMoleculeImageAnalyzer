function process_pmafiles(handles,hObject)

%% Check if image file is available
if ~isfield(handles,'file_mode')
    add_log(handles,'ERROR: no image file selected');
    return
end

%% Get map coeff.

if ~isfield(handles,'mapcoef')
    add_log(handles,'No mapping file available. Running in COUNT mode');
    handles.ch_mode='count';
    map=[];
else
    add_log(handles,'Running in FRET mode');
    handles.ch_mode='FRET';    
    map=handles.mapcoef;
end
% try
%     map=handles.mapcoef;
% catch err
%     if (strcmp(err.identifier,'MATLAB:nonExistentField'))
%         add_log(handles,'No mapping file available. Running in COUNT mode');        
%         ch_mode='count';
%         map=[];
%     end
% end
guidata(hObject,handles);

%% Run over the files/dir
handles.num_pmafiles_processed=0;
guidata(hObject,handles);
if isfield('file_mode',handles)
    add_log(handles,'ERROR: Lode image file first!!');
else
    if strcmp(handles.file_mode,'individual files')
        try
            pmafilename=handles.pmaloaded.filename;
            pmapathname=handles.pmaloaded.pathname;
        catch err
            if (strcmp(err.identifier,'MATLAB:nonExistentField'))
                add_log(handles,'ERROR: load image file first!!!');
            else
                % Display any other errors as usual.
                add_log(handles,'ERROR: Something is wrong...');rethrow(err);
            end
        end
        process_selected_pmafile(handles,pmafilename,pmapathname,map);
        
    elseif strcmp(handles.file_mode,'file folders')
        for dri=1:length(handles.pma_sub_dir_info)
            tmpname={};
            for j=1:length(handles.pmafileinfo{dri})
                tmpname{j}=handles.pmafileinfo{dri}(j).name;
            end
            add_log(handles,['working on: ' handles.pmapathinfo{dri}]);
            handles=process_selected_pmafile(handles,tmpname,handles.pmapathinfo{dri},map,hObject);
        end
    else
        add_log(handles,'ERROR: Something is wrong...');rethrow(err);
    end
end
add_log(handles,'Job finished!');


function handles=process_selected_pmafile(handles,pmafilename,pmapathname,map,hObject)
%     %% Get info. of image file loaded
[~, numfiles]=size(pmafilename);
if numfiles~=0
    %% make subfolder
    subfolder= get(handles.subfolder,'Value');
    if subfolder
        timetag=clock;
        timetag=[num2str(timetag(1)) num2str(timetag(2),'%02d') num2str(timetag(3),'%02d') ...
            num2str(timetag(4),'%02d') num2str(timetag(5),'%02d') num2str(floor(timetag(6)),'%02d') ];
        subfolder_name=['LeonResult.' timetag];
        mkdir(pmapathname,subfolder_name);
    else
        subfolder_name='';
    end
    
    %% Work on each pma files

    for cur_pma=1:numfiles
        handles.num_pmafiles_processed=handles.num_pmafiles_processed+1;
        handles=Leonembeded_v8(handles, pmafilename{cur_pma}, pmapathname, subfolder_name, map);
        add_log(handles,[num2str(handles.num_pmafiles_processed) '/' num2str(handles.num_pmafiles_all) ' (' num2str(round(handles.num_pmafiles_processed/handles.num_pmafiles_all*100)) ' %) done.']);
    end

end
add_log(handles,[num2str(numfiles) ' image files analyzed.']);
% guidata(hObject,handles);