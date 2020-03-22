function add_log(handles,new_entry,add_mode)
max_N_line=200;
if nargin < 3
    add_mode=1; %% shift downward by one line
end

logtext = get(handles.Logbox,'String');

[num_entry, ~]=size(logtext);
if add_mode==1
    if num_entry<max_N_line
    	logtext(2:num_entry+1)=logtext; 
    else
        logtext(2:max_N_line)=logtext(1:max_N_line-1); 
    end
end
logtext{1}=new_entry;

set(handles.Logbox,'String',logtext);

drawnow;

