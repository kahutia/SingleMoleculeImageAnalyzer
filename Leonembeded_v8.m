function handles=Leonembeded_v8(handles, pmafilename, pmapathname, subfolder_name, map)

%% input arguments preprocessing
f.pmafilename=pmafilename;
f.pmapathname=pmapathname;
f.subfolder_name=subfolder_name;

%% get parameters from GUI main and data files
[f, img, uip]=get_parameters(handles,f);
add_log(handles,[num2str(img.N_frame_all) ' frames in ' f.pmapathname f.pmafilename],-1);

%% check image selection and mapping method is correct
if isempty(map) && (strcmp(uip.imgselected,'Merged') || strcmp(uip.imgselected,'Co-local'))
    add_log(handles,'No mapping file available. Execution halted!!');
    errordlg('No mapping file available. Execution halted!!');
    return;
end

%% make Full map matrix (NOTE: this can be moved to mapcreator and map coefficient loader for better performance)
[D2Amap,A2Dmap]=calc_map_matrix(img.XpixelSize, img.YpixelSize/2,map);

%% find molecules
[handles,f,img, uip, mol]=get_mol_pos(handles,f,img, uip, map, A2Dmap, D2Amap);
mol=find_unique_pos(handles,uip,mol,map);

%% Get time traces
if uip.ext_trace
    [handles, mol]=get_traces(handles,f,img, uip, mol, map);
    save_binary_tr(handles,f,mol,img,map);
end

%% Save position of molecules
tmp_out=mol.pairsfound;
save(fullfile(f.pmapathname, f.subfolder_name, [strrep(strrep(f.pmafilename,'.pma',''),'.tif','') '.pos']),'tmp_out','-ascii','-append');
if isfield(mol,'pairs_firstfr')
    tmp_out=[mol.pairs_firstfr mol.pairs_lastfr];
    save(fullfile(f.pmapathname, f.subfolder_name, [strrep(strrep(f.pmafilename,'.pma',''),'.tif','') '.frames']),'tmp_out','-ascii','-append');
end

%% Save parameters
logfid=fopen(fullfile(f.pmapathname, f.subfolder_name, [strrep(f.pmafilename,'.pma','') '_leonlog.txt']),'wt');
fprintf(logfid,'Image autoscale : %d\n', uip.autoscale);
fprintf(logfid,'Image enhancer green ch : %d\n', uip.image_enhancer_green);
fprintf(logfid,'Image enhancer red ch : %d\n', uip.image_enhancer_red);
fprintf(logfid,'Image background green ch : %d\n', uip.image_BG_green);
fprintf(logfid,'Image background red ch : %d\n', uip.image_BG_red);

fprintf(logfid,'image size : (%d x %d) %d pixels\n', img.XpixelSize, img.YpixelSize, img.framesize);
fprintf(logfid,'first frame for the initial image analysis : %d\n', uip.first_frame_peakfinder);
fprintf(logfid,'# frame averaged for the initial image analysis : %d\n', uip.Num_frame_peakfinder);
fprintf(logfid,'Trace length : %d\n', img.N_frame_all);
fprintf(logfid,'channel selected : %s\n', uip.imgselected);

fprintf(logfid,'Local area size (pixel): %d\n', uip.localarea_size);
fprintf(logfid,'gaussian filter : %d\n', uip.use_gauss_filter);
fprintf(logfid,'co-local tolerance (pixel): %d\n', uip.co_local_tolerance);

if strcmp(uip.BGname,'Fixed value')
    fprintf(logfid,'Background correction method : %s\n',...
        [uip.BGname '(d:' num2str(uip.fixedBG_d) ', a:' num2str(uip.fixedBG_d) ')'] );
else
    fprintf(logfid,'Background correction method : %s\n', uip.BGname);
end

fprintf(logfid,'threshold method : %s\n', uip.threshold_method);
fprintf(logfid,'threshold index : %d\n', uip.threshold_index);

fprintf(logfid,'Mask type : %s\n', uip.MaskName);
fprintf(logfid,'Mask width: %g\n', uip.MaskWidth);

fprintf(logfid,'Drift correction : %s\n', uip.Driftname);
fprintf(logfid,'Peaks in all fr? %g\n', uip.peaks_all_fr);
fprintf(logfid,'peaks in multi frames : %s\n', uip.peaks_in_multi_fr);
fprintf(logfid,'Type of Peak : %s\n',uip.peak_type);
fclose(logfid);

pause(1);

%% save result summary
N_fid=fopen(fullfile(f.pmapathname, f.subfolder_name, [strrep(pmafilename,'.pma','') '_leon_Nmol.txt']),'wt');
fprintf(N_fid,'Left ch intensity : \t%d\n', mol.ch_int_avr_D(1));
fprintf(N_fid,'Right ch intensity : \t%d\n', mol.ch_int_avr_A(1));

if strcmp(uip.imgselected,'Co-local')
    fprintf(N_fid,'Co-local count : \t%d\n', mol.num_paired);
% 	mol.N_unpaired_eachfr(new_t_tag,1:2)=[num_peak num_peak_extra];
    fprintf(N_fid,'Left ch count : \t%d\n', mol.N_unpaired_eachfr(1,1));
    fprintf(N_fid,'Right ch count : \t%d\n', mol.N_unpaired_eachfr(1,2));
else
    fprintf(N_fid,'Count : \t%d\n', mol.num_paired);
end

fclose(N_fid);

N_fid=fopen(fullfile(f.pmapathname, f.subfolder_name, [strrep(pmafilename,'.pma','') '_leon_Nmol_all_frames.txt']),'wt');
if strcmp(uip.imgselected,'Co-local')
    output=[mol.N_pairs_eachfr mol.N_unpaired_eachfr(:,1) mol.N_unpaired_eachfr(:,2)];
    fprintf(N_fid,'%d%t%d%t%d%t\n', output);
else
%     fprintf(N_fid,'%d\n', mol.N_pairs_eachfr);
    output=[mol.N_pairs_eachfr; mol.ch_int_avr_D; mol.ch_int_avr_A];
%     fprintf(N_fid,'N_mol%tDint%tAint');
    fprintf(N_fid,'%d %f %f\n', output);
end

fclose(N_fid);


% end of main function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% start of subfuctions

%%
function mol=find_unique_pos(handles,uip,mol,map)
if isempty(map),    isFRET=0;
else, isFRET=1;end

tic
add_log(handles,'finding unique molecules...',1);

% methods=;

switch uip.peaks_in_multi_fr
    case 'Single image'
        mol.pairsfound=mol.pairs_eachfr{1};
        mol.num_paired=mol.N_pairs_eachfr(1);
    case 'Flatten image'
        %% Find unique molecules
        N_fr=length(mol.N_pairs_eachfr);

        % put the molecules found in the first frame into the unique list
        mol.uni_pairs_pos=mol.pairs_eachfr{1};
        mol.uni_pairs_firstfr=mol.uni_pairs_pos(:,1)*0+uip.first_frame_peakfinder;
        mol.uni_pairs_bg=mol.pairs_eachfr_bg{1};
        
        % run over all the frames
        for fri=2:N_fr
            % run over each molecule in the frame
            for mid=1:mol.N_pairs_eachfr(fri)
                c_x=mol.pairs_eachfr{fri}(mid,1);
                c_y=mol.pairs_eachfr{fri}(mid,2);
                % check if c_mol is already exist in unique pair list
                flg=1;
                for uid=1:length(mol.uni_pairs_pos)
                    if ((c_x - mol.uni_pairs_pos(uid,1))^2 + (c_y - mol.uni_pairs_pos(uid,2))^2 ) <= uip.co_local_tolerance^2
                        flg=0;break;
                    end
                end
                if flg     % add new molecule
                    if ~isFRET
                        mol.uni_pairs_pos=[mol.uni_pairs_pos; [c_x,c_y]];
                        mol.uni_pairs_bg=[mol.uni_pairs_bg; mol.pairs_eachfr_bg{fri}(mid,:)];
                        mol.uni_pairs_firstfr(end+1) = (fri-1)*uip.Num_frame_peakfinder + uip.first_frame_peakfinder;
                    else
                        mol.uni_pairs_pos=[mol.uni_pairs_pos; [c_x,c_y, mol.pairs_eachfr{fri}(mid,3), mol.pairs_eachfr{fri}(mid,4)]];
                        mol.uni_pairs_bg=[mol.uni_pairs_bg; mol.pairs_eachfr_bg{fri}(mid,:)];
                        mol.uni_pairs_firstfr(end+1) = (fri-1)*uip.Num_frame_peakfinder + uip.first_frame_peakfinder;
                    end
                end
            end
        end

        mol.pairsfound=mol.uni_pairs_pos;
        mol.pairs_eachfr_bg{1}=mol.uni_pairs_bg;
        mol.num_paired=length(mol.uni_pairs_firstfr);
        
    case 'Connected molecules'
        %% Find and connect molecules aross the frames
        N_fr=length(mol.N_pairs_eachfr);

        % put the molecules found in the first frame into the live molecule list
        live_pairs_pos=mol.pairs_eachfr{1};
        live_pairs_firstfr=live_pairs_pos(:,1)*0+uip.first_frame_peakfinder;
        dead_pairs_pos=[];
        dead_pairs_firstfr=[];
        dead_pairs_lastfr=[];
        
        % run over all the frames
        for fri=2:N_fr
            % run over each molecule in the frame
            new_mol=[];
            live_confirmed=false(1,length(live_pairs_pos(:,1)));
            
            for mid=1:mol.N_pairs_eachfr(fri)
                c_x=mol.pairs_eachfr{fri}(mid,1);
                c_y=mol.pairs_eachfr{fri}(mid,2);
                
                % check if c_mol is already exist in the previous frame
                flg=1;
                for lid=1:length(live_pairs_pos(:,1))
                    if ((c_x - live_pairs_pos(lid,1))^2 + (c_y - live_pairs_pos(lid,2))^2 ) <= uip.co_local_tolerance^2
                        flg=0;break;
                    end
                end
                if flg     % add new molecule to the live molecule list
                    if ~isFRET
                        new_mol=[new_mol; [c_x,c_y]];
                    else
                        new_mol=[new_mol; [c_x,c_y, mol.pairs_eachfr{fri}(mid,3), mol.pairs_eachfr{fri}(mid,4)]];
                    end
                else 
                    live_confirmed(lid)=true;
                end
            end
            
            % move dead mol to the dead list
            dead_pairs_pos=[dead_pairs_pos' live_pairs_pos(~live_confirmed,:)']';
            dead_pairs_firstfr=[dead_pairs_firstfr' live_pairs_firstfr(~live_confirmed)']';
            dead_pairs_lastfr=[dead_pairs_lastfr' live_pairs_firstfr(~live_confirmed)'*0+(fri-2)*uip.Num_frame_peakfinder + uip.first_frame_peakfinder]';
            
            % add new mol to the live list
            live_pairs_pos=[live_pairs_pos(live_confirmed,:)' new_mol']';
            [N_new_mol,~]=size(new_mol);
            live_pairs_firstfr=[live_pairs_firstfr(live_confirmed)' (1:N_new_mol)*0+(fri-1)*uip.Num_frame_peakfinder + uip.first_frame_peakfinder]';
        end
        
        % Add still alive pairs to the dead list
        dead_pairs_pos=[dead_pairs_pos' live_pairs_pos']';
        dead_pairs_firstfr=[dead_pairs_firstfr' live_pairs_firstfr']';
        dead_pairs_lastfr=[dead_pairs_lastfr' (1:length(live_pairs_firstfr))*0+(fri-1)*uip.Num_frame_peakfinder + uip.first_frame_peakfinder]';
            
        mol.pairsfound=dead_pairs_pos;
        mol.num_paired=length(dead_pairs_pos(:,1));
        mol.pairs_firstfr=dead_pairs_firstfr;
        mol.pairs_lastfr=dead_pairs_lastfr;
end


add_log(handles,[num2str(mol.num_paired) ' molecules found. (served in ' num2str(toc) 's)'],-1)


%%
function save_binary_tr(handles,f,mol,img,map)
tracefilename=[strrep(strrep(f.pmafilename,'.pma',''),'.tif','') '.traces2'];

% If a trace file with the same name is exist, back up the exist trace file
A=dir(fullfile(f.pmapathname, f.subfolder_name));
[num_file, ~]=size(A);
for i=1:num_file
    if strcmp(A(i).name,tracefilename)
        timetag=clock;
        timetag=[num2str(timetag(1)) '.' num2str(timetag(2)) '.' num2str(timetag(3)) ' '  ...
            num2str(timetag(4)) 'h' num2str(timetag(5)) 'm' num2str(floor(timetag(6))) 's'];
        copyfile([pmapathname A(i).name],[pmapathname A(i).name '.back' timetag]);
    end
end

% save the trace info
tracefid = fopen (fullfile(f.pmapathname, f.subfolder_name, tracefilename),'w');
fwrite(tracefid,img.N_frame_all,'int32');
fwrite(tracefid,mol.num_paired*2,'int32');

% save the trace
flg=0;
trace_output=zeros(1,img.N_frame_all*mol.num_paired*2);
if isempty(map)
    mol.acceptor=mol.donor*0;
end
for cur_frame=1:img.N_frame_all
    for cur_mol=1:mol.num_paired
        flg=flg+1;
        trace_output(flg)=mol.donor(cur_frame,cur_mol);
        flg=flg+1;
        trace_output(flg)=mol.acceptor(cur_frame,cur_mol);
    end
end

fwrite(tracefid,trace_output,'int16');
fclose(tracefid);

add_log(handles,['trace files are saved (' tracefilename ')']);

%%
function [handles, mol]=get_traces(handles,f,img, uip, mol, map)
% check if map is available
if isempty(map)
    isFRET=0;
else
    isFRET=1;
end

%% get background and make Mask
% add_log(handles,'checking background...',1);
% mol.background=get_bg(mol,uip,img, isFRET);
add_log(handles,'Making masks...',0);
mol.mask=get_mask(mol, uip , img, isFRET);

%% get time traces
tic
Donor=zeros(img.N_frame_all,mol.num_paired);
Acceptor=zeros(img.N_frame_all,mol.num_paired);

switch f.imgtype
    case {'pma','shf'}
    fseek(f.pmafid, 4, 'bof');    % skip the head;
end

start_flg=1;
for api=1:f.N_appended+1
    end_flg=start_flg+img.N_frame_per_file(api)-1;
    gfr2fid(start_flg:end_flg)=api;
    gfr2lfr(start_flg:end_flg)=1:img.N_frame_per_file(api);
    start_flg=end_flg+1;
end

% read all frames
skip_notificatoin=20;
for curr_frame=1:img.N_frame_all
    % update log at every 10th frame
    if rem(curr_frame,20)==1
        add_log(handles,['Extracting intensity from frames ' num2str(curr_frame) ' to ' num2str(curr_frame+skip_notificatoin-1)...
            ' (' num2str(floor((curr_frame+skip_notificatoin-1)/img.N_frame_all*100)) ' %)'],-1);
    end
    
    % read a frame
    switch f.imgtype
        case 'pma'
            pmatmp=fread(f.pmafid,img.framesize,'uint8=>uint8');
            pmatmp=reshape(pmatmp,img.XpixelSize,img.YpixelSize,1);
            pmatmp=flipud(rot90(pmatmp));
            pmatmp=mean(pmatmp,3);

        case 'shf'
            pmatmp=fread(f.pmafid,img.framesize,'uint16=>uint16','ieee-be');
            pmatmp=reshape(pmatmp,img.XpixelSize,img.YpixelSize,1);
            pmatmp=double(rot90(pmatmp,3));
            pmatmp=flipud(pmatmp);
        case {'tif_solis','tif_imageJ'}
            f.TifLink(gfr2fid(curr_frame)).setDirectory(gfr2lfr(curr_frame));
            pmatmp=double(f.TifLink(gfr2fid(curr_frame)).read());
            pmatmp=fliplr(pmatmp);
    end
    
    % calc intensity of each molecules
    for i=1:mol.num_paired
        %%%% donor %%%%
        subframe=pmatmp(mol.pairsfound(i,1)-uip.localarea_size:mol.pairsfound(i,1)+uip.localarea_size,...
            mol.pairsfound(i,2)-uip.localarea_size:mol.pairsfound(i,2)+uip.localarea_size);
        Donor(curr_frame,i)=sum(sum((double(subframe)-mol.pairs_eachfr_bg{1}(i,1)).*mol.mask.D(:,:,i)));
    end
    if isFRET
        for i=1:mol.num_paired
            %%%% acceptor %%%%
            subframe=pmatmp(mol.pairsfound(i,3)-uip.localarea_size:mol.pairsfound(i,3)+uip.localarea_size,...
                mol.pairsfound(i,4)-uip.localarea_size:mol.pairsfound(i,4)+uip.localarea_size);
            Acceptor(curr_frame,i)=sum(sum((double(subframe)-mol.pairs_eachfr_bg{1}(i,2)).*mol.mask.A(:,:,i)));
        end
    end
    
    
    
end

mol.donor=Donor;
if isFRET,    mol.acceptor=Acceptor;    end

t=toc;
add_log(handles,[num2str(curr_frame) ' frames are processed. (served in ' num2str(t) 'sec)'],-1);

switch f.imgtype
    case {'pma','shf'}
        fclose(f.pmafid);
    case {'tif_solis','tif_imageJ'}
        f.TifLink(1).close();
end

%%
function Mask=get_mask(mol, uip , img, isFRET)
Mask.D=zeros(uip.box_size,uip.box_size,mol.num_paired);
switch uip.MaskName
    case 'Gaussian (given)'
        xx=meshgrid(-uip.localarea_size:uip.localarea_size);
        yy=xx';
        tmpMask=zeros(uip.box_size,uip.box_size,9);
        flg=0;
        for x0=-0.5:0.5:+0.5
            for y0=-0.5:0.5:+0.5
                flg=flg+1;
                tmpMask(:,:,flg)=2*exp( -((xx-x0).^2+(yy-y0).^2) ./ uip.MaskWidth );
            end
        end
        
        for i=1:mol.num_paired
            %%%% donor %%%%
            subframe=img.ini_frame(mol.pairsfound(i,1)-uip.localarea_size:mol.pairsfound(i,1)+uip.localarea_size,...
                mol.pairsfound(i,2)-uip.localarea_size:mol.pairsfound(i,2)+uip.localarea_size);

            for flg=1:9
                cur_score=sum(sum((double(subframe)-mol.pairs_eachfr_bg{1}(i,1)).*tmpMask(:,:,flg)));
                if flg==1
                    score=cur_score;
                    Mask.D(:,:,i)=tmpMask(:,:,flg);
                else
                    if cur_score > score
                        score=cur_score;
                        Mask.D(:,:,i)=tmpMask(:,:,flg);
                    end
                end
            end
            
            %%%% acceptor %%%%
            if isFRET
                subframe=img.ini_frame(mol.pairsfound(i,3)-uip.localarea_size:mol.pairsfound(i,3)+uip.localarea_size,...
                    mol.pairsfound(i,4)-uip.localarea_size:mol.pairsfound(i,4)+uip.localarea_size);
                
                for flg=1:9
                    cur_score=sum(sum((double(subframe)-mol.pairs_eachfr_bg{1}(i,2)).*tmpMask(:,:,flg)));
                    if flg==1
                        score=cur_score;
                        Mask.A(:,:,i)=tmpMask(:,:,flg);
                    else
                        if cur_score > score
                            score=cur_score;
                            Mask.A(:,:,i)=tmpMask(:,:,flg);
                        end
                    end
                end
            end
        end
        
    case 'Circular'
        Mask.D=zeros(uip.box_size,uip.box_size,mol.num_paired);
        for mid=1:mol.num_paired
            for i=1:uip.box_size
                for j=1:uip.box_size
                    if (i-uip.localarea_size-1)^2 + (j-uip.localarea_size-1)^2 <= uip.MaskWidth^2
                        Mask.D(i,j,mid)=1;
                    end
                end
            end
        end
        if isFRET
            Mask.A=Mask.D;
        end
end

%%
function [handles,f,img, uip, mol]=get_mol_pos(handles,f,img, uip, map, A2Dmap, D2Amap)
tic
add_log(handles,'peak finding...');

%% prepare figures
% if get(handles.show_images,'Value')
% 
% end


%% peak finding from all frames
if uip.peaks_all_fr
    last_frame_peakfinder=img.N_frame_all-uip.Num_frame_peakfinder+1;
else
    last_frame_peakfinder=uip.first_frame_peakfinder;
end

% make global frame number for each tif file
start_flg=1;
for api=1:f.N_appended+1
    end_flg=start_flg+img.N_frame_per_file(api)-1;
    gfr2fid(start_flg:end_flg)=api;
    gfr2lfr(start_flg:end_flg)=1:img.N_frame_per_file(api);
    
    start_flg=img.N_frame_per_file(api)+1;
end

%% Read frames for peakfinder
new_t_tag=0;
add_log(handles,'reading pma...',-1);
for fri=uip.first_frame_peakfinder:uip.Num_frame_peakfinder:last_frame_peakfinder
    new_t_tag=new_t_tag+1;

    if uip.use_saved_img
%         frame=imread([f.pmapathname f.pmafilename '_TAV.tif']);
%         frame=double(flipud(frame));
        tmp_data=load([f.pmapathname f.pmafilename '_TAV.raw'],'-mat');
        frame=tmp_data.raw_image;
%         frame=double(flipud(frame));
    else
        switch f.imgtype
            case 'pma'
                pmatmp=fread(f.pmafid,img.framesize*uip.Num_frame_peakfinder,'uint8=>uint8');
                pmatmp=double(reshape(pmatmp,img.XpixelSize,img.YpixelSize,uip.Num_frame_peakfinder));
                pmatmp=flipud(rot90(pmatmp));
                frame_org=mean(pmatmp,3);

                switch uip.peak_type
                    case 'mean'
                        frame=frame_org;
                    case 'Allan'
                        frame=.5*sqrt((pmatmp(:,:,2:end)-pmatmp(:,:,1:end-1)).^2);
                    case 'std'
                        frame=std(pmatmp,0,3);
                end

            %%
            case 'shf'
                pmatmp=fread(f.pmafid,img.framesize*uip.Num_frame_peakfinder,'uint16=>uint16', 'ieee-be');
                pmatmp=reshape(pmatmp,img.XpixelSize,img.YpixelSize,uip.Num_frame_peakfinder);
                frame_org=mean(double(pmatmp),3);
                switch uip.peak_type
                    case 'mean'
                        frame=frame_org;
                    case 'Allan'
                        frame=double(pmatmp);
                        frame=.5*sqrt((frame(:,:,2:end)-frame(:,:,1:end-1)).^2);
                    case 'std'
                        frame=std(double(pmatmp),0,3);
                end
                frame_org=flipud(rot90(frame_org,3));
                frame=flipud(rot90(frame,3));
            case {'tif_solis','tif_imageJ'}
                if strcmp(uip.peak_type,'std') || strcmp(uip.peak_type,'Allan')
                    errordlg('Peak finding from std/Allan images is not supported for tif format.');
                    return;
                end
                for gfr=fri:fri+uip.Num_frame_peakfinder-1
                    f.TifLink(gfr2fid(gfr)).setDirectory(gfr2lfr(gfr));
                    tmp=double(f.TifLink(gfr2fid(gfr)).read());
                    if gfr==fri
                        frame=tmp;
                    else
                        frame=frame+tmp;
                    end
                end
                frame=fliplr(frame)/uip.Num_frame_peakfinder;
        end
    end
    
    if uip.image_pivot
        frame=frame';
    end
	D_image=frame(:,1:img.YpixelSize/2);
    A_image=frame(:,img.YpixelSize/2+1:img.YpixelSize);
    
%     add_log(handles,'peak finding',-1);    
    %% bg corraction and rescaling (for image display only)
    if uip.autoscale
        uip.image_BG_green=round(median(D_image(:)));
        uip.image_BG_red=round(median(A_image(:)));
        uip.image_enhancer_green=(round(max(D_image(:)))-uip.image_BG_green );
        uip.image_enhancer_red=(round(max(A_image(:)))-uip.image_BG_red );
        
        handles.ColorEnhancerGreen.String=num2str(uip.image_enhancer_green);
        handles.ColorEnhancerRed.String=num2str(uip.image_enhancer_red);
        handles.img_bg_green.String=num2str(uip.image_BG_green);
        handles.img_bg_red.String=num2str(uip.image_BG_red);
    end
    D_image_bg=(D_image-uip.image_BG_green)/uip.image_enhancer_green*255;
    A_image_bg=(A_image-uip.image_BG_red)/uip.image_enhancer_red*255;
    
    frame_bg=zeros(img.XpixelSize,img.YpixelSize);
    frame_bg(:,1:img.YpixelSize/2)=D_image_bg;
    frame_bg(:,img.YpixelSize/2+1:img.YpixelSize)=A_image_bg;
    
    %% Create Merged image;
    if ~isempty(map)
        [Intensity_merged, ~]=mergeimage(D_image,A_image,map);
        [Intensity_merged_bg, RGB_merged,~]=mergeimage(D_image_bg,A_image_bg,map);
    end
    
    %% get the intensity level of each chanel
    tmp_margin=10;
    tmp=frame(tmp_margin:end-tmp_margin,tmp_margin:img.YpixelSize/2-tmp_margin);
    mol.ch_int_avr_D(new_t_tag)=mean(tmp(:));
    mol.ch_int_median_D(new_t_tag)=median(tmp(:));
    tmp=frame(tmp_margin:end-tmp_margin,img.YpixelSize/2+tmp_margin:end-tmp_margin);
    mol.ch_int_avr_A(new_t_tag)=mean(tmp(:));
    mol.ch_int_median_A(new_t_tag)=median(tmp(:));
    
    %% peak finding
    switch uip.imgselected
        case 'Merged'
            tmp_img=Intensity_merged;
            tmp_bg=(uip.fixedBG_d+uip.fixedBG_a)/2;
        case {'Donor','Co-local'}
            tmp_img=D_image;
            tmp_bg=uip.fixedBG_d;
        case 'Acceptor'
            tmp_img=A_image;
            tmp_bg=uip.fixedBG_a;
        case 'Whole'
            tmp_img=frame;
            tmp_bg=(uip.fixedBG_d+uip.fixedBG_a)/2;
    end
    
    % cal temp bg_matrix and find peaks
    [bg_matrix,bg_matrix_std]=get_bg(tmp_img,uip.BGname,tmp_bg,uip.localarea_size);
    [mol_pos, num_peak, ~]=peakfinderembeded(tmp_img,uip.localarea_size,uip.threshold_method,uip.threshold_index,...
        uip.use_gauss_filter,uip.MaskWidth,bg_matrix,bg_matrix_std);

    % complete bg_matrix
    switch uip.imgselected
        case 'Merged'
            bg_matrix=[];
            [tmp_bg_matrix,~]=get_bg(D_image,uip.BGname,uip.fixedBG_d,uip.localarea_size);
            bg_matrix(:,1:img.YpixelSize/2)=tmp_bg_matrix;
            [tmp_bg_matrix,~]=get_bg(A_image,uip.BGname,uip.fixedBG_a,uip.localarea_size);
            bg_matrix(:,img.YpixelSize/2+1:img.YpixelSize)=tmp_bg_matrix;
        case {'Donor','Co-local'}
            [tmp_bg_matrix,tmp_bg_matrix_std]=get_bg(A_image,uip.BGname,uip.fixedBG_a,uip.localarea_size);
            bg_matrix(:,img.YpixelSize/2+1:img.YpixelSize)=tmp_bg_matrix;
        case 'Acceptor'
            [tmp_bg_matrix,~]=get_bg(D_image,uip.BGname,uip.fixedBG_a,uip.localarea_size);
            bg_matrix(:,img.YpixelSize/2+1:img.YpixelSize)=bg_matrix;
            bg_matrix(:,1:img.YpixelSize/2)=tmp_bg_matrix;
        case 'whole'
            % nothing to do
    end
    
    % re-define backgrounds 
%     bg_matrix=[];
%     [tmp_bg_matrix,~]=get_bg(D_image,uip.BGname,uip.fixedBG_d,uip.localarea_size);
% 	bg_matrix(:,1:img.YpixelSize/2)=tmp_bg_matrix;
%     [tmp_bg_matrix,tmp_bg_matrix_std]=get_bg(A_image,uip.BGname,uip.fixedBG_a,uip.localarea_size);
% 	bg_matrix(:,img.YpixelSize/2+1:img.YpixelSize)=tmp_bg_matrix;
    
    if strcmp(uip.imgselected,'Co-local')
        [mol_pos_extra, num_peak_extra, ~]=peakfinderembeded(A_image,uip.localarea_size,uip.threshold_method,uip.threshold_index,...
                                            uip.use_gauss_filter,uip.MaskWidth,bg_matrix(:,img.YpixelSize/2+1:img.YpixelSize),tmp_bg_matrix_std);
        mol_pos_extra(:,2)=mol_pos_extra(:,2)+img.YpixelSize/2;
    end
    
    %% re-calculate backgrounds if peak is choosen from std/Allan image
    switch uip.peak_type
        case {'std','Allan'}
            bg_matrix=[];
            [tmp_bg_matrix,~]=get_bg(frame_org(:,1:img.YpixelSize/2),uip.BGname,uip.fixedBG_d,uip.localarea_size);
            bg_matrix(:,1:img.YpixelSize/2)=tmp_bg_matrix;
            [tmp_bg_matrix,~]=get_bg(frame_org(:,img.YpixelSize/2+1:img.YpixelSize),uip.BGname,uip.fixedBG_a,uip.localarea_size);
            bg_matrix(:,img.YpixelSize/2+1:img.YpixelSize)=tmp_bg_matrix;
    end
    
    %% update backgrounds
    background=1:num_peak;
    for mli=1:num_peak
        background(mli)=bg_matrix(mol_pos(mli,1),mol_pos(mli,2));
    end
    if strcmp(uip.imgselected,'Co-local')    
        background_extra=1:num_peak_extra;
        for mli=1:num_peak_extra
            if mol_pos_extra(mli,1) > 0 && mol_pos_extra(mli,1) < img.XpixelSize && ...
                    mol_pos_extra(mli,2) > img.YpixelSize/2 && mol_pos_extra(mli,2) < img.YpixelSize 
                background_extra(mli)=bg_matrix(mol_pos_extra(mli,1),mol_pos_extra(mli,2));
            else
                background_extra(mli)=0;
            end
        end    
    end
    
    
    %% peak pairing
    if isempty(map)
        pairsfound=mol_pos;
        num_paired=num_peak;
        pairsfound_bg=background';
        if num_paired~=0
            if strcmp(uip.imgselected,'Acceptor')
                pairsfound(:,2)=pairsfound(:,2)+img.YpixelSize/2;
            end
        end
    else
        pairsfound=[];
        pairsfound_bg=[];
        if strcmp(uip.imgselected,'Acceptor')
            num_paired=0;
            for i=1:num_peak
                mappedx= A2Dmap(mol_pos(i,1),mol_pos(i,2),1);
                mappedy= A2Dmap(mol_pos(i,1),mol_pos(i,2),2);
                % exclude molecules that are too close to the edge
                if (mappedx>uip.box_size && mappedx<img.XpixelSize-uip.box_size &&...
                        mappedy>uip.box_size && mappedy<img.YpixelSize/2-uip.box_size)
                    num_paired=num_paired+1;
                    pairsfound(num_paired,1)=mappedx;
                    pairsfound(num_paired,2)=mappedy;
                    pairsfound(num_paired,3)=mol_pos(i,1);
                    pairsfound(num_paired,4)=mol_pos(i,2)+img.YpixelSize/2;
                    pairsfound_bg(num_paired,1)=bg_matrix(mappedx,mappedy);
                    pairsfound_bg(num_paired,2)=background(i);
                end
            end
        else
            num_paired=0;
            for i=1:num_peak
                mappedx=D2Amap(mol_pos(i,1),mol_pos(i,2),1);
                mappedy=D2Amap(mol_pos(i,1),mol_pos(i,2),2)+img.YpixelSize/2;
                % exclude molecules that are too close to the edge
                if (mappedx>uip.box_size && mappedx<img.XpixelSize-uip.box_size &&...
                        mappedy>img.YpixelSize/2+uip.box_size && mappedy<img.YpixelSize-uip.box_size)
                    num_paired=num_paired+1;
                    pairsfound(num_paired,1)=mol_pos(i,1);
                    pairsfound(num_paired,2)=mol_pos(i,2);
                    pairsfound(num_paired,3)=mappedx;
                    pairsfound(num_paired,4)=mappedy;
                    pairsfound_bg(num_paired,1)=background(i);
                    pairsfound_bg(num_paired,2)=bg_matrix(mappedx,mappedy);
                end
            end
        end
        
    end
    
    if strcmp(uip.imgselected,'Co-local')
        % remove molecules that are not detected in the acceptor channel
        pairsfound_unpaired=pairsfound;
        pairsfound_bg_unpaired=pairsfound_bg;
        pairsfound=[];
        pairsfound_bg=[];
        num_unpaired=num_paired;
        num_paired=0;
        for i=1:num_unpaired
            for gfr=1:num_peak_extra
                if (pairsfound_unpaired(i,3)-mol_pos_extra(gfr,1))^2 < uip.co_local_tolerance^2 && ...
                   (pairsfound_unpaired(i,4)-mol_pos_extra(gfr,2))^2 < uip.co_local_tolerance^2
                    num_paired=num_paired+1;
                    pairsfound=[pairsfound; [squeeze(pairsfound_unpaired(i,1:2)) squeeze(mol_pos_extra(gfr,1:2))]];
                    pairsfound_bg(num_paired,[1 2])=[background(i) background_extra(gfr)];
                    break
                end
            end
        end
    end
    
    mol.pairs_eachfr{new_t_tag}=pairsfound;
    mol.N_pairs_eachfr(new_t_tag)=num_paired;
    mol.pairs_eachfr_bg{new_t_tag}=pairsfound_bg;
    
    if strcmp(uip.imgselected,'Co-local')
        mol.N_unpaired_eachfr(new_t_tag,1:2)=[num_peak num_peak_extra];
    end
    
    %% Mark the images
    circled=frame_bg;
    if num_paired>0
        if ~isempty(map)
            if strcmp(uip.imgselected,'Co-local')
                circled=add_circle(circled,pairsfound_unpaired(:,1),pairsfound_unpaired(:,2),100);
                circled=add_circle(circled,mol_pos_extra(:,1),mol_pos_extra(:,2),100);
            end
            circled=add_circle(circled,pairsfound(:,3),pairsfound(:,4),170);
            merged_circled=add_circle(Intensity_merged_bg,pairsfound(:,1),pairsfound(:,2),170);
        end
        circled=add_circle(circled,pairsfound(:,1),pairsfound(:,2),170);
    else
        if ~isempty(map)
            merged_circled=Intensity_merged_bg;
        end
    end
    
    %% Show images
    if get(handles.show_images,'Value')
        %%% prepare figure
        fhd_Leon_img_org=findobj('Tag','Leon_img_org');
        if isempty(fhd_Leon_img_org)
            fhd_Leon_img_org=figure('Name','Original image','Tag','Leon_img_org','Numbertitle','off','Units','pixels','Position',[0 50 512 512]);
            axes('position',[0 0 1 1]); colormap(img.chrcol)
        end
        
        fhd_Leon_img_org_circled=findobj('Tag','Leon_img_org_circled');
        if isempty(fhd_Leon_img_org_circled)
            fhd_Leon_img_org_circled=figure('Name','Original image','Tag','Leon_img_org_circled','Numbertitle','off','Units','pixels','Position',[0 50 512 512]);
            axes('position',[0 0 1 1]); colormap(img.chrcol)
        end
        
        if ~isempty(map)
            fhd_Leon_img_mer=findobj('Tag','Leon_img_mer');
            if isempty(fhd_Leon_img_mer)
                fhd_Leon_img_mer=figure('Name','Merged image','Tag','Leon_img_mer','Numbertitle','off','Units','pixels','Position',[512 50 256 512]);
                axes('position',[0 0 1 1]); colormap(img.chrcol)
            end
            fhd_Leon_img_mer_circled=findobj('Tag','Leon_img_mer_circled');
            if isempty(fhd_Leon_img_mer_circled)
                fhd_Leon_img_mer_circled=figure('Name','Merged image','Tag','Leon_img_mer_circled','Numbertitle','off','Units','pixels','Position',[512 50 256 512]);
                axes('position',[0 0 1 1]); colormap(img.chrcol)
            end
        end
        
        %%% update figure
        set(groot,'CurrentFigure',fhd_Leon_img_org);
        image(frame_bg);    axis image;
        text(5,img.YpixelSize/20,['fr=' num2str(fri) ', N=' num2str(num_paired)],'Color','y');

        set(groot,'CurrentFigure',fhd_Leon_img_org_circled);
        image(circled);    axis image;
        text(5,img.YpixelSize/20,['fr=' num2str(fri) ', N=' num2str(num_paired)],'Color','y');
        if ~isempty(map)
            set(groot,'CurrentFigure',fhd_Leon_img_mer);
            image(Intensity_merged_bg);axis image;
%             image(RGB_merged*15);axis image;    % multiplying 15 is to make image brighter
            set(groot,'CurrentFigure',fhd_Leon_img_mer_circled);
            image(merged_circled);axis image;
        end
        drawnow;
    end
	add_log(handles,['working on frame #' num2str(fri) ' of ' num2str(last_frame_peakfinder) '. (' num2str(num_paired) ' molecules found)'],-1);
    
    
    %% Save images
    if fri==uip.first_frame_peakfinder
        % archive the first image for later use (ex. BG correction)
        img.ini_frame=frame;
        
        % save image
        save_image(frame_bg,'_org',f,img.chrcol);
        save_image(circled,'_orgO',f,img.chrcol);
        if ~isempty(map)
            save_image(RGB_merged*15,'_merRGB',f,img.chrcol);
            save_image(Intensity_merged_bg,'_mer',f,img.chrcol);
            save_image(merged_circled,'_merO',f,img.chrcol);
        end
    end
end

add_log(handles,[num2str(new_t_tag*uip.Num_frame_peakfinder) ' frames analyzed in ' num2str(toc) 's.'],-1);

%%
pause(1)
if get(handles.show_images,'Value')
    close(fhd_Leon_img_org);
    close(fhd_Leon_img_org_circled);
    if ~isempty(map)
        close(fhd_Leon_img_mer);
        close(fhd_Leon_img_mer_circled);
    end
end

%%
function save_image(frame,suffix,f,chrcol)
%% shared header
tagstruct.Photometric = Tiff.Photometric.RGB;
tagstruct.BitsPerSample = 8;
tagstruct.SamplesPerPixel = 3;
tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
tagstruct.Software = 'L.Trotsky V6';

[~,~,color_depth]=size(frame);
if color_depth==1
    imdata = uint8(ind2rgb(uint8(frame),chrcol*255));
else
    imdata = uint8(frame*256);
end

if strcmp(f.imgtype,'tif_solis') || strcmp(f.imgtype,'tif_imageJ')
    fileext='tif';
else
    fileext=f.imgtype;
end

% tifname=[f.pmapathname f.subfolder_name strrep(strrep(f.pmafilename,'.pma',''),'.tif','') suffix '.tif'];
tifname=fullfile(f.pmapathname, f.subfolder_name, [ strrep(f.pmafilename,['.' fileext],'') suffix '.tif']);
t = Tiff(tifname,'w');
tagstruct.ImageLength = size(imdata,1);
tagstruct.ImageWidth = size(imdata,2);
t.setTag(tagstruct);
t.write(imdata);
t.close();

%%
function [f, img, uip]=get_parameters(handles,f)
add_log(handles,'Getting file info...');

%% Load Color map
img.chrcol=load('singlecolormap.dat');

%%  mask type and BG method
uip.MaskTypeHdl=get(handles.collectionmethod,'SelectedObject');
uip.MaskName=get(uip.MaskTypeHdl,'String');
uip.BGmethodHdl=get(handles.BGmethod,'SelectedObject');
uip.BGname=get(uip.BGmethodHdl,'String');
uip.fixedBG_d = str2double(get(handles.fixedBG_d,'String'));
uip.fixedBG_a = str2double(get(handles.fixedBG_a,'String'));

%% thresholds
uip.threshold_method=get(get(handles.th_method,'SelectedObject'),'String');
uip.threshold_index = str2double(get(handles.threshold_index,'String'));

%% peak type
uip.peak_type=get(get(handles.peak_type,'SelectedObject'),'String');

%% Get drift correction method
uip.Driftname='No correction';

%% Get the image choosen for peak finding
PeakFindimgHdl=get(handles.PeakImg,'SelectedObject');
uip.imgselected=get(PeakFindimgHdl,'String');

%% Get other parameters
uip.autoscale=get(handles.img_autoscale,'Value');
uip.co_local_tolerance = str2double(get(handles.co_local_tolerance,'String'));
uip.Num_frame_peakfinder = str2double(get(handles.Num_frame_peakfind,'String'));
uip.localarea_size = str2double(get(handles.local_area_size,'String'));
uip.box_size=uip.localarea_size*2+1;

uip.MaskWidth= str2double(get(handles.Maskwidth,'String'));
uip.use_gauss_filter= get(handles.GaussFilter,'Value'); % 0:off, 1: on
uip.use_saved_img= get(handles.Use_saved_img,'Value'); % 0:off, 1: on
uip.ext_trace= get(handles.ext_trace,'Value'); % 0:off, 1: on

peaks_in_multi_fr_Hdl=get(handles.peaks_in_multi_fr,'SelectedObject');
uip.peaks_in_multi_fr=get(peaks_in_multi_fr_Hdl,'String');
uip.peaks_all_fr= get(handles.peak_all_fr,'Value'); % 0:off, 1: on
if uip.peaks_all_fr
    uip.first_frame_peakfinder = 1;
else
    uip.peaks_in_multi_fr='Single image';
    uip.first_frame_peakfinder = str2double(get(handles.first_frame_peakfinder,'String'));
end

uip.image_enhancer_green=str2double(get(handles.ColorEnhancerGreen,'String'));
uip.image_enhancer_red=str2double(get(handles.ColorEnhancerRed,'String'));
uip.image_BG_green=str2double(get(handles.img_bg_green,'String'));
uip.image_BG_red=str2double(get(handles.img_bg_red,'String'));
uip.image_pivot=get(handles.img_pivot,'Value');




%% GET file info
f.imgtype=f.pmafilename(end-2:end);
if strcmp(f.imgtype,'tif')
    %% make a dialog box to let user to choose the type
    if 0
        f.imgtype='tif_solis';
    else
        f.imgtype='tif_imageJ';
    end
end
f.N_appended=0;

% if f.pmapathname(end)~='\',    f.pmapathname=[f.pmapathname '\']; end

switch f.imgtype
    case {'shf','pma'}
    	f.pmafid=fopen(fullfile(f.pmapathname, f.pmafilename),'r');
    case 'tif_solis'
        
        f.TifLink(1) = Tiff(fullfile(f.pmapathname, f.pmafilename), 'r');
        % check if tif file is splitted
        tmp_file_list=dir(fullfile(f.pmapathname, [f.pmafilename(1:end-4) '_X*' '.tif']));
        if ~isempty(tmp_file_list)
            f.N_appended=length(tmp_file_list);
        end
        for api=1:f.N_appended
            f.TifLink(api+1) = Tiff(fullfile(f.pmapathname, [f.pmafilename(1:end-4) '_X' num2str(api) '.tif']), 'r');
        end

    case 'tif_imageJ'
        f.TifLink(1) = Tiff(fullfile(f.pmapathname, f.pmafilename), 'r');
        % check if tif file is splitted
        tmp_file_list=dir(fullfile(f.pmapathname, [f.pmafilename(1:end-8) '_*' '.ome.tif']));
        if ~isempty(tmp_file_list)
            f.N_appended=length(tmp_file_list);
        end
        for api=1:f.N_appended
            f.TifLink(api+1) = Tiff(fullfile(f.pmapathname, tmp_file_list(api).name), 'r');
        end
end


%% Other Parameter setting
switch f.imgtype
    case 'pma'
        img.XpixelSize=fread(f.pmafid,1,'uint16');
        img.YpixelSize=fread(f.pmafid,1,'uint16');

        fseek(f.pmafid, 0, 'eof');
        pmafilesize=ftell(f.pmafid);
        img.N_frame_per_file=(pmafilesize-4)/img.XpixelSize/img.YpixelSize;   % 4 for header. (2byte * 2)
        fseek(f.pmafid, 4, 'bof');
        img.N_frame_all=img.N_frame_per_file;
    case 'shf'
        img.XpixelSize=fread(f.pmafid,1,'uint16','ieee-be');
        img.YpixelSize=fread(f.pmafid,1,'uint16','ieee-be');

        fseek(f.pmafid, 0, 'eof');
        pmafilesize=ftell(f.pmafid);
        img.N_frame_per_file=(pmafilesize-4)/(2*img.XpixelSize*img.YpixelSize);   % 4 for header. (2byte * 2)
        img.N_frame_all=img.N_frame_per_file;
        fseek(f.pmafid, 4, 'bof');
    case {'tif_solis','tif_imageJ'}
        tif_info = imfinfo([f.pmapathname f.pmafilename]);
        img.XpixelSize=tif_info(1).Width;
        img.YpixelSize=tif_info(1).Height;

        img.N_frame_per_file(1)=numel(tif_info);
        for api=1:f.N_appended
%             tif_info_app = imfinfo([f.pmapathname f.pmafilename(1:end-4) '_X' num2str(api) '.tif']);
            tif_info_app = imfinfo([f.pmapathname tmp_file_list(api).name]);
            img.N_frame_per_file(api+1)=numel(tif_info_app);
        end
        img.N_frame_all=sum(img.N_frame_per_file);    
end

img.framesize=img.XpixelSize*img.YpixelSize;