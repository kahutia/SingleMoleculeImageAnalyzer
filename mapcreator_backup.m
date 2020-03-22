% close all; clear all;
function mapcreator(hObject,handles,arg1,arg2)
%% argment inputcheck
if nargin<4 % this is only for backward compatibility
    arg1=1;
    arg2=2;
end

%%
warning off;

%% parameters
% uip.localarea_size = str2double(get(handles.local_area_size,'String'));
uip.localarea_size = 2;
uip.Maskwidth = str2double(get(handles.Maskwidth,'String'));
uip.co_local_tolerance= str2double(get(handles.co_local_tolerance,'String'));
uip.use_gauss_filter= get(handles.GaussFilter,'Value'); % 0:off, 1: on; 
NumFrameProc=10;

uip.threshold_method=get(get(handles.th_method,'SelectedObject'),'String');
uip.threshold_index = str2double(get(handles.threshold_index,'String'));
uip.threshold_index = 3;

% uip.BGmethodHdl=get(handles.BGmethod,'SelectedObject');
% uip.BGname=get(uip.BGmethodHdl,'String');
uip.BGname='Channel median';
uip.fixedBG_d = str2double(get(handles.fixedBG_d,'String'));
uip.fixedBG_a = str2double(get(handles.fixedBG_a,'String'));

uip.image_enhancer_green= str2num(get(handles.ColorEnhancerGreen,'String')); 
uip.image_enhancer_red= str2num(get(handles.ColorEnhancerRed,'String')); 
uip.extra_acceptor_enhancer=1;

uip.art_bg=25;  % this is to make the cross hair line visible

uip.imgselected='mapping';



%% File open

add_log(handles,'Reading bead image...');

if ~isfield(handles,'pma_root_path')
    handles.pma_root_path='D:\DNA pulldown\data\MS\180315 bead 540 560 - Copy';
end


[filename, pathname]=uigetfile([handles.pma_root_path '\*.*'], 'CCD images');
handles.pma_root_path=pathname;
uip.filename=filename;
ispma=strcmp(filename(end-3:end),'pma');


if ispma
    fid=fopen(strcat(pathname,filename),'r');

    uip.XpixelSize=fread(fid,1,'uint16');
    uip.YpixelSize=fread(fid,1,'uint16');

    %% PMA read
    frame=zeros(uip.XpixelSize,uip.YpixelSize);
    for i=1:NumFrameProc
        tmp=fread(fid,uip.XpixelSize*uip.uip.YpixelSize,'uint8');
        tmp=reshape(tmp,uip.XpixelSize,uip.YpixelSize);
        frame=frame+tmp;
    end
    frame=frame'/i;


else % tif
    TifLink = Tiff([pathname filename], 'r');
%     frame=zeros(TifLink.ImageLength,TifLink.ImageWidth,'uint16');
    for i=1:NumFrameProc
        TifLink.setDirectory(i);
        tmp=TifLink.read();
        if i==1
            frame=tmp;
        else
            frame=frame+tmp;
        end
    end
    frame=fliplr(double(frame))/i;
    TifLink.close();
    [uip.YpixelSize,uip.XpixelSize]=size(frame);
end

uip.frame=frame;clear('frame');

%% prepare figure
hdl_mcfig=figure('Units','pixels','Position',[0 50 512 512]);clf;
ax1=axes('position',[0.025 0 0.9 .95]);

acceptor_sld = uicontrol('Style', 'slider',...
        'Min',0,'Max',100,'Value',uip.extra_acceptor_enhancer,...
        'Units','normalized',...
        'Position', [0.94 0.5 0.05 0.4],...
        'Callback', {@callback_acceptor_sld, uip});

flg_togo = uicontrol('Style', 'pushbutton',...
        'Units','normalized',...
        'String', 'Choose Pairs',...
        'Position', [0.7 0.95 0.25 0.04],...
        'Callback', {@callback_flg_togo,uip,hdl_mcfig,hObject,handles});
    
    
%% Draw image
draw_fig(hdl_mcfig,uip); 
guidata(hObject, handles);
end

function callback_flg_togo(sObject, eventdata, uip,hdl_mcfig,hObject,handles)
%% update figure
uip=draw_fig(hdl_mcfig,uip);    % this is to update uip after GUI action
mapped_frame=uip.circled_frame;
circle_color=150;

%% get position of 3 molecules (6pt) for mapping
flg_pt=zeros(4,2);  % odd for donor, even for acceptor
for i=1:8
    [flg_pt(i,2), flg_pt(i,1)]=ginput(1);
    flg_pt=floor(flg_pt);
    
    % find local maxima
    subframe=uip.frame(flg_pt(i,1)-uip.localarea_size:flg_pt(i,1)+uip.localarea_size,...
        flg_pt(i,2)-uip.localarea_size:flg_pt(i,2)+uip.localarea_size);
    [tmp, tmpindex1]=max(subframe);
    [~, tmpindex2]=max(tmp);
    flg_pt(i,1)=tmpindex1(tmpindex2)-uip.localarea_size-1+flg_pt(i,1);
    flg_pt(i,2)=tmpindex2-uip.localarea_size-1+flg_pt(i,2);
    
    % draw circle to the selected molecule
    uip.circled_frame=add_circle(uip.circled_frame,flg_pt(i,1),flg_pt(i,2),circle_color);
    image(uip.circled_frame+uip.art_bg);
%     imagesc(circled_frame);
    axis image;drawnow;
end


msgbox('Calculating map coefficients..','notice');


%% Get temporary map coef. with the 3 molecules
x=[flg_pt(1,:) flg_pt(3,:) flg_pt(5,:) flg_pt(7,:)]';
y=[flg_pt(2,:) flg_pt(4,:) flg_pt(6,:) flg_pt(8,:)]';
ft = fittype( 'mapfn3mol( x, a0, b0, a1, b1, a2, b2)' );
f = fit( x, y, ft, 'StartPoint', [0, 0, 1, 1, 1, 1] );


%% find all molecules
% find in donor channel
bg_matrix1=get_bg(uip.frame(:,1:256),uip.BGname,uip.fixedBG_d,uip.localarea_size);
[mol_pos1, num_mol1]=peakfinderembeded(uip.frame(:,1:256),uip.localarea_size,uip.threshold_method,uip.threshold_index,...
    uip.use_gauss_filter,uip.Maskwidth,bg_matrix1);
% find in acceptor channel
bg_matrix2=get_bg(uip.frame(:,256:end),uip.BGname,uip.fixedBG_a,uip.localarea_size);
[mol_pos2, num_mol2]=peakfinderembeded(uip.frame(:,256:end),uip.localarea_size,uip.threshold_method,uip.threshold_index,...
    uip.use_gauss_filter,uip.Maskwidth,bg_matrix2);

mol_pos2(:,2)=mol_pos2(:,2)+uip.XpixelSize/2;

mol_pos=[mol_pos1; mol_pos2];
num_mol=num_mol1+num_mol2;


%% Find pairs from the molecules with temporary map coef.
num_paired=0;
pairsfound=zeros(num_mol,4);
for i=1:num_mol
    if mol_pos(i,2) < 250
        mappedx= round(f.a0 + f.a1*mol_pos(i,1) + f.a2*mol_pos(i,2));
        mappedy= round(f.b0 + f.b1*mol_pos(i,1) + f.b2*mol_pos(i,2));
        
        % check if there is acceptor at the mapped position
        for k=1:num_mol
            if 0 % perfect match
                if (mol_pos(k,1)==mappedx && mol_pos(k,2)==mappedy)
                    num_paired=num_paired+1;
                    pairsfound(num_paired,1)=mol_pos(i,1);
                    pairsfound(num_paired,2)=mol_pos(i,2);
                    pairsfound(num_paired,3)=mappedx;
                    pairsfound(num_paired,4)=mappedy;
                end
            else % matching with tolerance
                tol=uip.co_local_tolerance+1;  % px
                if ( (mol_pos(k,1)>=mappedx-tol && mol_pos(k,1)<=mappedx+tol )&& (mol_pos(k,2)>=mappedy-tol && mol_pos(k,2)<=mappedy+tol))
                    num_paired=num_paired+1;
                    pairsfound(num_paired,1)=mol_pos(i,1);
                    pairsfound(num_paired,2)=mol_pos(i,2);
                    pairsfound(num_paired,3)=mappedx;
                    pairsfound(num_paired,4)=mappedy;
                end
            end
        end
    end
end
pairsfound=pairsfound(1:num_paired,:);

%% Get map coef with the pairs found in 3-moleule-fit
for i=1:num_paired
    x(i*2-1,:)=pairsfound(i,1);
    x(i*2,:)=pairsfound(i,2);
    y(i*2-1,:)=pairsfound(i,3);
    y(i*2,:)=pairsfound(i,4);
end

ft = fittype( 'mapfn( x, a0, a1, a2, a3, a4, a5, a6, b0, b1, b2, b3, b4, b5, b6)' );
f = fit( x, y, ft, 'StartPoint', [0, 1, 1, 0.1, 0.1 0.01, 0.01, 0, 1, 1, 0.1, 0.1 0.01, 0.01] );

%% Find the pair of all the molecules with complete map coef.
for i=1:num_mol
    if mol_pos(i,2) < 250
        mappedx= round(f.a0 + f.a1*mol_pos(i,1) + f.a2*mol_pos(i,2) + f.a3*mol_pos(i,1)^2 + f.a4*mol_pos(i,2)^2 + f.a5*mol_pos(i,1)^3 + f.a6*mol_pos(i,2)^3);
        mappedy= round(f.b0 + f.b1*mol_pos(i,1) + f.b2*mol_pos(i,2) + f.b3*mol_pos(i,1)^2 + f.b4*mol_pos(i,2)^2 + f.b5*mol_pos(i,1)^3 + f.b6*mol_pos(i,2)^3);
        
        % check if there is acceptor at the mapped position
        tol=uip.co_local_tolerance;
        for k=1:num_mol
            %             if (mol_pos(k,1)==mappedx && mol_pos(k,2)==mappedy),
            if (mol_pos(k,1)>mappedx-tol && mol_pos(k,1)<mappedx+tol && mol_pos(k,2)>mappedy-tol && mol_pos(k,2)<mappedy+tol)
                mapped_frame=add_circle(mapped_frame,mol_pos(i,1),mol_pos(i,2),circle_color);
                mapped_frame=add_circle(mapped_frame,mappedx,mappedy,circle_color);
            end
        end
    end
end

%% save mapping coefficient
writefilename=[handles.pma_root_path strrep(uip.filename,'.pma','') '.mapcoef'];
mapfid=fopen(writefilename,'w');
A=coeffnames(f);
B=coeffvalues(f);
for i=1:size(A)
    fprintf(mapfid,'%s %g\n',A{i},B(i));
end
map(1,:)=B(1:size(A)/2);
map(2,:)=B(size(A)/2+1:size(A));
fclose(mapfid);

%% update map info
handles.mapcoef=map;
guidata(hObject,handles);

%% Create Merged image;
% D_image=frame(:,1:YpixelSize/2);
% A_image=frame(:,YpixelSize/2+1:YpixelSize);

[Intensity_merged, RGB_merged, ~]=map_image(uip.frame_left,uip.frame_right,map);


%% plot and save image file

% Load Color map
chrcol=load('singlecolormap.dat');

figure('Units','pixels','Position',[512 50 256 512]);  hdl12=gcf;
axes('position',[0 0 1 1]); colormap(chrcol)
image(Intensity_merged);axis image; drawnow;

figure('Units','pixels','Position',[768 50 256 512]);  hdl13=gcf;
axes('position',[0 0 1 1]); colormap(chrcol)
image(RGB_merged);axis image; drawnow;



figure(hdl_mcfig);clf;
axes('position',[0 0 1 1]); colormap(chrcol)
image(mapped_frame);
axis image;
colormap(chrcol);
% colormap(gray);
drawnow;

% writefilename=[pathname filename '.tif'];

tagstruct.Photometric = Tiff.Photometric.RGB;
tagstruct.BitsPerSample = 8;
tagstruct.SamplesPerPixel = 3;
tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
tagstruct.Software = 'L.Trotsky';
% original image
imdata = uint8(ind2rgb(uint8(uip.frame/255),chrcol*255));
tifname=[handles.pma_root_path strrep(uip.filename,'.pma','') '_1.tif'];
t = Tiff(tifname,'w');
tagstruct.ImageLength = size(imdata,1);
tagstruct.ImageWidth = size(imdata,2);
t.setTag(tagstruct);
t.write(imdata);
t.close();
% circled image
% tmp_image(:,1:floor(XpixelSize)/2)=mapped_frame(:,1:floor(XpixelSize)/2)*image_enhancer_red;
% tmp_image(:,ceil(XpixelSize)/2:XpixelSize)=mapped_frame(:,ceil(XpixelSize)/2:XpixelSize)*image_enhancer_red;
imdata = uint8(ind2rgb(uint8(mapped_frame),chrcol*255));
tifname=[handles.pma_root_path strrep(uip.filename,'.pma','') '_2.tif'];
t = Tiff(tifname,'w');
tagstruct.ImageLength = size(imdata,1);
tagstruct.ImageWidth = size(imdata,2);
t.setTag(tagstruct);
t.write(imdata);
t.close();
% merged image
% imdata = uint8(ind2rgb(uint8(Intensity_merged),chrcol*255));
imdata = uint8(RGB_merged*256);
tifname=[handles.pma_root_path strrep(uip.filename,'.pma','') '_3.tif'];
t = Tiff(tifname,'w');
tagstruct.ImageLength = size(imdata,1);
tagstruct.ImageWidth = size(imdata,2);
t.setTag(tagstruct);
t.write(imdata);
t.close();



%% get user confirmation
add_log(handles,'mapping coefficients are loaded');
tmp=msgbox('Mapping file generated!','notice','replace');
waitfor(tmp); 

close(hdl_mcfig);
close(hdl12);
close(hdl13);
guidata(hObject, handles);

end

function uip=draw_fig(fhd,uip)
figure(fhd);

uip.extra_acceptor_enhancer=fhd.Children(2).Value;
frame=uip.frame-min(uip.frame(:));
frame=frame/max(frame(:))*255;
frame(frame<0)=0;
uip.circled_frame=frame;    % this is to allocate memory

uip.frame_left=frame(:,1:floor(uip.XpixelSize/2))*uip.image_enhancer_green/255;
uip.circled_frame(:,1:floor(uip.XpixelSize/2))=uip.frame_left;
uip.frame_right=frame(:,ceil(uip.XpixelSize/2):uip.XpixelSize)*uip.image_enhancer_red*uip.extra_acceptor_enhancer/255;
uip.circled_frame(:,ceil(uip.XpixelSize/2):uip.XpixelSize)=uip.frame_right;

% frame_left=frame_left-min(frame_left(:));
% frame_left=frame_left/max(frame_left(:));
% frame_right=frame_right-min(frame_right(:));
% frame_right=frame_right/max(frame_right(:));

uip.circled_frame(:,255:256)=255;   % add a line bet. donor and acceptor channel
image(uip.circled_frame+uip.art_bg);
axis image;axis off;
colormap('gray');% colormap(chrcol);

end

function callback_acceptor_sld(sObject, eventdata, uip)
% uip.extra_acceptor_enhancer=sObject.Value;
draw_fig(sObject.Parent,uip);
end


function [Intensity_merged, RGB_merged, FRET_merged]=map_image(D_image, A_image, map)

[X_size, Y_size]=size(D_image);
Intensity_merged=zeros(X_size,Y_size);
FRET_merged=zeros(X_size,Y_size);
for x=1:X_size
    for y=1:Y_size
        mappedx= round(map(1,1) +...
            map(1,2)*x + map(1,3)*y +...
            map(1,4)*x^2 + map(1,5)*y^2 +...
            map(1,6)*x^3 + map(1,7)*y^3);
        mappedy= round(map(2,1) +...
            map(2,2)*x + map(2,3)*y +...
            map(2,4)*x^2 + map(2,5)*y^2 +...
            map(2,6)*x^3 + map(2,7)*y^3);
        mappedy=mappedy-Y_size;
        if mappedx < 1, mappedx = 1; end
        if mappedy < 1, mappedy = 1; end
        if mappedx > X_size, mappedx = X_size; end
        if mappedy > Y_size, mappedy = Y_size; end
        Intensity_merged(x,y)=(D_image(x,y)+A_image(mappedx,mappedy))/2;
        FRET_merged(x,y)=A_image(mappedx,mappedy)/(D_image(x,y)+A_image(mappedx,mappedy));
        mapped_red(x,y)=A_image(mappedx,mappedy);
    end
end

RGB_merged(:,:,1)=mapped_red/255*1.5;
RGB_merged(:,:,2)=D_image/255*1.5;
RGB_merged(:,:,3)=zeros(X_size, Y_size);
% 
% 
% bg=median(mapped_red(:));
% RGB_merged(:,:,1)=(mapped_red-bg)/(max(mapped_red(:))-bg)*1.5;
% bg=median(D_image(:));
% RGB_merged(:,:,2)=(D_image-bg)/(max(D_image(:))-bg)*1.5;
% RGB_merged(:,:,3)=zeros(X_size, Y_size);

end





