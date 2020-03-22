
function subdirs=get_pma_all_subdir(Leon_pathname)
if nargin<1
    Leon_pathname=uigetdir();
end

subdirs={};
subdirs=get_subdirs(Leon_pathname,subdirs);
subdirs=subdirs(2:end);


function subdirs=get_subdirs(c_dirname,subdirs)
% disp('start_new_call:');
% disp(c_dirname);
if isempty(strfind(c_dirname,'Leon'))
    subdirs=[subdirs; c_dirname];

    tmp_dirs=dir(c_dirname);

    % get directories only
    tmp_dirs=tmp_dirs(3:end);
    for dri=1:length(tmp_dirs)
        if tmp_dirs(dri).isdir   % work only for the directories
            subdirs=get_subdirs([c_dirname '\' tmp_dirs(dri).name],subdirs);
        end
    end
end
