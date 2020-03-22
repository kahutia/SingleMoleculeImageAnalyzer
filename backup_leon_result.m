
function backup_leon_result(root_path)

if nargin<1
    root_path='';
end

Leon_pathname=uigetdir(root_path,'Select folder');

%% get all the sub directory names
subdirs={};
subdirs=get_subdirs(Leon_pathname,subdirs);

%% select Leon Result folders
subdirs_leon={};
N_leon=0;
for i=1:length(subdirs)
    if ~isempty(strfind(subdirs{i},'LeonResult')) && isempty(strfind(subdirs{i},'Leon_backup'))
        N_leon=N_leon+1;
        subdirs_leon{N_leon,1}=subdirs{i};
    end
end

%% Work on the Leon Results
% user_ans=menu([num2str(length(subdirs_leon)) ' Leon result found.'],'Remove','Back up','Do nothing');
user_ans=get_user_choice([num2str(length(subdirs_leon)) ' Leon result found.'],Leon_pathname,subdirs_leon);

% if user_ans==1
%     disp('Removing files....');
%     user_confirm=menu('Are you sure to remove?', 'Yes','No');
%     if user_confirm==1
%         remove_dirs(subdirs_leon);   
%     end
% elseif user_ans==2
%     disp('Moving files to ');
%     move_dirs(Leon_pathname,subdirs_leon) 
% else
%     disp('halted!');
% end

function user_ans=get_user_choice(q_str,Leon_pathname,subdirs_leon)
hdl_trotsky=gcf;
hdl_trotsky.Units='pixels';
d=dialog('Position',[hdl_trotsky.Position(1)+300 hdl_trotsky.Position(2)+250 190 100]);

txt = uicontrol('Parent',d,...
    'Style','text',...
    'Position',[10 40 190 40],...
    'String',q_str);

% btn = uicontrol('Parent',d,...
%     'Position',[10 20 60 25],...
%     'String','Remove',...
%     'Enable','inactive',...
%     'Callback',{@remove_dirs,gcf,subdirs_leon});

btn = uicontrol('Parent',d,...
    'Position',[25 20 65 25],...
    'String','Back up',...
    'Callback',{@move_dirs,gcf,Leon_pathname,subdirs_leon});

btn = uicontrol('Parent',d,...
    'Position',[100 20 65 25],...
    'String','Do nothing',...
    'Callback','delete(gcf)');
user_ans=3;
           

function subdirs=get_subdirs(c_dirname,subdirs)
% Note that this is a recursive function

subdirs=[subdirs; c_dirname];

% get directories only
tmp_dirs=dir(c_dirname);
tmp_dirs=tmp_dirs(3:end);

for dri=1:length(tmp_dirs)
    if tmp_dirs(dri).isdir   % work only for the directories
        subdirs=get_subdirs([c_dirname '\' tmp_dirs(dri).name],subdirs);
    end
end


function move_dirs(src,eventdata,dia_hdl,Leon_pathname,subdirs)
close(dia_hdl);
disp('Moving files to ');
N_dir=length(subdirs);


%% make subfolder
if N_dir>0
    timetag=clock;
    timetag=[num2str(timetag(1)) num2str(timetag(2),'%02d') num2str(timetag(3),'%02d') ...
        num2str(timetag(4),'%02d') num2str(timetag(5),'%02d') num2str(floor(timetag(6)),'%02d') ];
    subfolder_name=['Leon_backup' timetag '\'];
    mkdir(Leon_pathname,['Leon_backup' timetag]);
    
    disp(['Backup root: ' Leon_pathname '\' subfolder_name]);
    N_file_rm=0;
    for i=1:N_dir
%         if strfind(subdirs{i},'LeonResult')
            reserved_name=strrep(subdirs{i},Leon_pathname,'');
            disp(['moving: ' reserved_name]);
            movefile(subdirs{i},[Leon_pathname '\' subfolder_name reserved_name],'f');
            N_file_rm=N_file_rm+1;
%         end
    end
else
    N_file_rm=0;
end
disp([num2str(N_file_rm) ' directories moved!!!']);



function remove_dirs(src,eventdata,dia_hdl,subdirs)
close(dia_hdl);
disp('Removing files....');

N_dir=length(subdirs);

N_file_rm=0;
for i=1:N_dir
%     if strfind(subdirs{i},'LeonResult')
        disp(['removing: ' subdirs{i}]);
        rmdir(subdirs{i},'s');
        N_file_rm=N_file_rm+1;
%     end
end

disp([num2str(N_file_rm) ' directories removed!!!']);


