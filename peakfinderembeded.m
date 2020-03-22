%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Peak finder for TROTSKY V7 by SHK
% NOT compatible with earier version of Trotsky
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1) is local maximum in the local area? ex. 9x9 box?
% 2) is the local maximum lager than 'threshold_index'*'standard deviation
%        of local area'?

function [mol_pos, num_mol, background]=peakfinderembeded(frame,localarea_size,threshold_method,threshold_index,use_gauss_filter,sigma,bg_matrix,bg_matrix_std)
% function [mol_pos, num_mol, background]=peakfinderembeded(frame,uip,~)

%% parameters
border=localarea_size+1;
border_wide=localarea_size*2+1;
box_size=localarea_size*2+1;

[Xsize, Ysize]=size(frame);
num_mol=0;

%% 2D-Gaussian Filter filtering of the image

if use_gauss_filter
    % Make Gaussian Mask
    FixedGmask=zeros(box_size,box_size);
    for i=1:box_size
        for j=1:box_size
            FixedGmask(j,i)=exp( -((i-localarea_size-1).^2+(j-localarea_size-1).^2) ./ sigma );
        end
    end
    FixedGmask=FixedGmask / sum(sum(FixedGmask));
    
    % Convolute frame with 2D-Guassian Mask
    frameConvo=zeros(Xsize,Ysize);
    for i=border:Xsize-border
        for j=border:Ysize-border
            subframe=frame(i-localarea_size:i+localarea_size,...
                j-localarea_size:j+localarea_size);
            frameConvo(i,j)=sum(sum(subframe.*FixedGmask));
        end
    end
    frame=frameConvo;
end


%% Define threshold
switch threshold_method
    case 'constant'
        th(1:Xsize,1:Ysize) = threshold_index;
    case 'x sigma'
        th=zeros(Xsize,Ysize);
        for i=border_wide:Xsize-border_wide
            for j=border_wide:Ysize-border_wide
                th(i,j)=bg_matrix_std(i,j)*threshold_index;
            end
        end
end

%% Find Local Maxima
mol_pos=[];
background=[];


for i=border_wide:Xsize-border_wide
    for j=border_wide:Ysize-border_wide
        subframe=frame(i-localarea_size:i+localarea_size,j-localarea_size:j+localarea_size);
        subframe_wide=frame(i-localarea_size*2:i+localarea_size*2,j-localarea_size*2:j+localarea_size*2);
        
        [max_valX, max_index]=max(subframe);
        [~, max_indexY]=max(max_valX);
        
        % is local maxima?
        if max_index(max_indexY)==localarea_size+1 && max_indexY==localarea_size+1

            % above threshold?
            if frame(i,j) > th + bg_matrix(i,j)
                % accept as a molecule if the intensity of the centor box
                % is larger than background*threshold_index
                num_mol=num_mol+1;
                mol_pos(num_mol,1)=i;
                mol_pos(num_mol,2)=j;
                background(num_mol)=bg_matrix(i,j);
            end
        end
    end
end
