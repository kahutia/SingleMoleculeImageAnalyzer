function frame=add_circle(frame, i, j,circle_color)
% circle_color=170;
[Xsize, Ysize]=size(frame);
if isscalar(i)
    if i>3 && i < Xsize-3 && j>3 && j<Ysize-3
        frame(i-3,j-1)=circle_color;
        frame(i-3,j+1)=circle_color;
        frame(i-1,j-3)=circle_color;
        frame(i-1,j+3)=circle_color;
        frame(i+1,j-3)=circle_color;
        frame(i+1,j+3)=circle_color;
        frame(i+3,j-1)=circle_color;
        frame(i+3,j+1)=circle_color;
    end
else
    num_mol=size(i);
    iarray=i;
    jarray=j;
    for curr_mol=1:num_mol
        i=iarray(curr_mol);
        j=jarray(curr_mol);
        if i>3 && i < Xsize-3 && j>3 && j<Ysize-3
            frame(i-3,j-1)=circle_color;
            frame(i-3,j+1)=circle_color;
            frame(i-1,j-3)=circle_color;
            frame(i-1,j+3)=circle_color;
            frame(i+1,j-3)=circle_color;
            frame(i+1,j+3)=circle_color;
            frame(i+3,j-1)=circle_color;
            frame(i+3,j+1)=circle_color;
        end
    end
end
