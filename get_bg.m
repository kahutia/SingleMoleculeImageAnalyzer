function [background,background_std]=get_bg(frame,bg_method,userdef_bg,localarea_size)
[Xsize,Ysize]=size(frame);
box_size_wide=localarea_size*4+1;

switch bg_method
    case 'Fixed value'
        background(1:Xsize,1:Ysize) = userdef_bg;
        background_std(1:Xsize,1:Ysize) = std(frame(:));
    case 'Channel median'
        background(1:Xsize,1:Ysize) = median(frame(:));
        background_std(1:Xsize,1:Ysize) = std(frame(:));
    case 'Local median'
        background=zeros(Xsize,Ysize);
        background_std=zeros(Xsize,Ysize);
        border_wide=localarea_size*2+1;
        borderbox_wide=zeros(box_size_wide);
        borderbox_wide([1:localarea_size box_size_wide-localarea_size+1:box_size_wide],:)=1;
        borderbox_wide(:,[1:localarea_size box_size_wide-localarea_size+1:box_size_wide])=1;
        borderbox_wide_ind=borderbox_wide(:)==1;

        for i=border_wide:Xsize-border_wide
            for j=border_wide:Ysize-border_wide
                subframe_wide=frame(i-localarea_size*2:i+localarea_size*2,j-localarea_size*2:j+localarea_size*2);
                cur_borderbox=double(subframe_wide(borderbox_wide_ind));
                background(i,j)=median(cur_borderbox);
                background_std(i,j)=std(cur_borderbox);
            end
        end
end